#!/bin/bash

diamond makedb --in /projectnb/staghorn/nr.faa \
    --db nr \
    --taxonmap prot.accession2taxid.gz \
    --taxonnodes nodes.dmp \
    -v -p 28 -t /scratch
