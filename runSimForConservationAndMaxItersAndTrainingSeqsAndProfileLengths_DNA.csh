#!/bin/tcsh
time dist/profillicSimulation_DNA ProfillicSimulation_runSimForConservationAndMaxItersAndTrainingSeqs.cfg --seed=98103 --saveResultsParentDirectory=results/DNA_runSimForConservation_$1_andMaxIters_$2_AndTrainingSeqs_$3_AndProfileLengths_$4 --conservationRates=$1 --maxIterations=$2 --euclideanDistanceMinimum_iteration=1E-10 --numTrainingSequencesPerProfiles=$3 --profileLengths=$4

