from collections import defaultdict, Counter
import csv
import json
import gzip
import logging
import pandas
import sys

hits_fields = [
        'qseqid','qlen','qstart','qend',
        'gene_id','slen','sstart','send',
        'gene_desc','rna_acc','protein_acc','gene_symbol',
        'staxids','score' # smaller score is better hit
]

def cast_fields(r) :
    to_int = ('qlen','qstart','qend','slen','sstart','send')
    for f in to_int :
        r[f] = int(r[f])

    r['score'] = float(r['score'])
    return r

def annotate_transcripts(
        transcripts, # names of all transcripts
        in_hits, # list of dicts
        taxa, # dict of k:v as <taxclass name>:<set of taxids>
        taxclass_priority, # list of (int, taxclass name), smaller int is higher priority
        cluster_map, # dict of k:v as <seqname>:<parent seqname>
        source # string specifying source of hit, e.g. blast or minimap2
    ):

    stats = Counter()

    stats['transcripts'] = len(transcripts)
    logging.info('read transcriptome, {} sequences'.format(len(transcripts)))

    stats['num cluster sequences'] = len(cluster_map)
    stats['num clusters'] = len(set([_['parent_tid'] for _ in cluster_map.values()]))
    stats['num hits processed'] = len(in_hits)

    # annot is a dict of k:v as <seqname>:<list of hits>
    annot = defaultdict(list)
    for i,r in enumerate(in_hits) :
        annot[r['qseqid']].append(cast_fields(r))
    stats['num transcripts with hits'] = len(annot)

    # read hits output
    # expect the following fields in hits
    def pick_best_hit(hits) :

        min_score = min([_['score'] for _ in hits])
        best_hits = [_ for _ in hits if _['score'] == min_score]

        def filter_by_title_text(l,text) :
            kept = [_ for _ in best_hits if text not in _['gene_desc']]
            if len(kept) > 0 :
                return kept
            return l

        # prefer hits that don't have PREDICTED in the title
        if len(best_hits) > 1 :
            best_hits = filter_by_title_text(hits,'PREDICTED')

        # prefer hits that don't have LOW QUALITY PROTEIN in the title
        if len(best_hits) > 1 :
            best_hits = filter_by_title_text(hits,'LOW QUALITY PROTEIN')

        # prefer hits that don't have hypothetical protein in the title
        if len(best_hits) > 1 :
            best_hits = filter_by_title_text(hits,'hypothetical protein')

        # otherwise just return first (i.e. random) hit
        return best_hits[0]

    # then process each transcript hits to map to a single sseqid
    # also compute the taxonomic class, and make sure there is only one
    # per transcript
    hits = defaultdict(list)
    for tid, tid_hits in annot.items() :
        best_hit = pick_best_hit(tid_hits)
        best_hit['source'] = source

        # set the annot key for this tid to the best hit for use in
        # no hit resolution next step
        annot[tid] = best_hit

        # some proteins have multiple taxon ids, not sure why
        classes = []
        for taxonid in best_hit['staxids'].split(';') :
            classes.extend([k for k,v in taxa.items() if taxonid in v])
        classes = list(set(classes))

        if len(classes) == 0 : # no taxid, despite hit? map to other
            stats['no hit classes found'] += 1
            classes.append('org_other')

        if len(classes) > 1 : # multiple classes

            priority_class = min((p,c) for p,c in taxclass_priority if c in classes)

            classes = [priority_class[1]]

            stats['multiple taxclasses'] += 1

        assert len(classes) == 1

        best_hit['taxclass'] = classes[0]
        hits[best_hit['gene_id']].append(best_hit)

    # now attempt to map tids with no hits to their parent tids
    nohit_tids = set(transcripts).difference(set(annot))

    for tid in nohit_tids :
        tid_cluster = cluster_map[tid]
        parent_tid = tid_cluster['parent_tid']
        if parent_tid in annot :

            # the parent sequence of this tid's cluster has a hit, assign it to tid
            parent_hit = annot[parent_tid]

            tid_cluster.update(parent_hit)

            annot[tid] = parent_hit

            hits[parent_hit['gene_id']].append(tid_cluster)

            stats['nohit mapped to parent hit'] += 1

        else :
            stats['unmapped nohits'] += 1

    # find the new nohit tids after assigning to parent hit
    nohit_tids = set(transcripts).difference(set(annot))

    return {
            'hits': hits,
            'nohits': nohit_tids,
            'annot': annot,
            'stats': stats
           }

def clustered_transcript_hits(
        transcripts, # names of transcripts to group by cluster
        cluster_map # dict of k:v as <seqname>:<parent seqname>
    ) :

    hits = defaultdict(list)
    for tid in transcripts:
        tid_cluster = cluster_map[tid]

        hits[tid_cluster['parent_tid']].append(tid_cluster)

    return hits

def hits_to_gtf(hits,fn) :

    # construct the gtf file from hits
    gtf_fields = ['qseqid','sstart','send','source']

    with open(fn,'wt') as g :
        g_out = csv.writer(g,delimiter='\t', quoting=csv.QUOTE_NONE,quotechar='')
        for pid,hits in hits.items() :
            df = pandas.DataFrame(hits)
            df.index = df['qseqid']

            attrs = {_:df.iloc[0].get(_) for _ in df.columns if _ not in gtf_fields}
            attrs.setdefault('taxclass','nohit')
            attr_str = ''.join(['{} "{}";'.format(k,v) for k,v in attrs.items()])
            g_out.writerow([
                pid,df.iloc[0].get('source','.'),'gene',
                int(df['sstart'].min()),int(df['send'].max()),
                '.','.','.',
                attr_str
            ])
            for tid,r in df.sort_values('sstart').iterrows()  :

                attrs = {
                    'gene_id': pid,
                    'transcript_id': r['qseqid'],
                    'taxa_ids': r.get('staxids','NA'),
                    'taxclass': r.get('taxclass','nohit'),
                }

                g_out.writerow([
                    pid,r.get('source','.'),'transcript',
                    int(r['sstart']),int(r['send']),
                    '.','.','.',
                    ''.join(['{} "{}";'.format(k,v) for k,v in attrs.items()])
                ])

from itertools import chain
def merge_hits(*hits) :

    all_hit_keys = set(chain(*[_.keys() for _ in hits]))

    all_hits = {}
    for key in all_hit_keys :
        all_hits[key] = list(chain(*[_.get(key,[]) for _ in hits]))

    return all_hits
