#!/bin/tcsh
mkdir results
foreach x (`seq 1 1 $2` )
    ./runSimForConservationAndMaxItersAndTrainingSeqsAndProfileLengths_DNA.csh $1 $x $3 $4 &
end
#EOF

