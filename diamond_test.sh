#!/bin/bash

diamond blastx -v --compress 1 -f 6 \
    -o blast_test.tab.gz --sensitive \
    --query Calabash_test.fa \
    --db /projectnb/staghorn/cab/labadorf/nr \
    -p 16 \
    --taxonlist 9606 \
    -t $TMPDIR/
