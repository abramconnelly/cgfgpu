#!/bin/bash

bandconc="../../../../helper/band-concordance"

afn="../hu826751-GS03052-DNA_B01-035e.band"
bfn="../hu34D5B9-035e.band"

for s in {0..34} ; do
  for S in `seq $s 34` ; do
    echo $s $S
    $bandconc -s $s -S $S $afn $bfn
  done
done
