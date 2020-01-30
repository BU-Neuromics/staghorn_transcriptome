import csv
from pathlib import Path
import re
import sys

if len(sys.argv) != 3 :
    print('Usage: python sort_salmon_by_taxclass.py <gtf_fn> <salmon_fn>')
    sys.exit(1)

gtf = Path(sys.argv[1])
cnts = Path(sys.argv[2])

patt = re.compile('taxclass "([^"]*)"')

taxclass_map = {}
taxclasses = set()
with open(gtf,'rt') as f :
    print('reading GTF to build a transcript id -> taxclass map')
    for r in csv.reader(f,delimiter='\t') :
        tid = r[0]
        feature = r[2]
        if feature == 'gene' :
            taxclass = patt.search(r[-1]).group(1)
            taxclasses.add(taxclass)
            taxclass_map[tid] = taxclass

with open(cnts,'rt') as f :

    print('sorting counts into respective taxclass files')

    f = csv.reader(f)
    header = next(f)

    # open all the file pointers
    fps = {}
    for tc in taxclasses :
        tc_fn = '{}__{}{}'.format(cnts.stem,tc,cnts.suffix)
        fps[tc] = open(tc_fn,'wt')
        fps[tc].write(','.join(header)+'\n')

    # write each count into the appropriate taxclass file
    for r in f :
        taxclass = taxclass_map[r[0]]
        fps[taxclass].write(','.join(r)+'\n')

[_.close() for _ in fps.values()]

print('done')
