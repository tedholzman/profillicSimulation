#!/bin/tcsh
foreach x (`seq 1 1 $2` )
 ./runSimForConservationAndMaxItersAndTrainingSeqs_AA.csh $1 $x $3 &
end
#EOF

