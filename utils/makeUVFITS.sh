#!/bin/bash

## change variables
projID='AGBT17B_360'
session='03'
obj='NGC4258_Field'

for i in {0..6}
do
	bm=$i
	idlToSdfits -j -l -o $projID"_"$session"_"$obj"_Beam"$bm"_baseliend.UVFITS" $projID"_"$session"_"$obj"_Beam"$bm"_edge_Baselined_ss.fits"
	idlToSdfits -j -l -o $projID"_"$session"_"$obj"_Beam"$bm"_noBaseline.UVFITS" $projID"_"$session"_"$obj"_Beam"$bm"_edge_noBaseline_ss.fits"
done
