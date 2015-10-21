#!/bin/tcsh
time dist/profillicSimulation_AA ProfillicSimulation_runSimForConservationAndMaxItersAndTrainingSeqs_revised.cfg --seed=98103 --saveResultsParentDirectory=results_revised/AA_runSimForConservation_$1_andMaxIters_$2_AndTrainingSeqs_$3 --conservationRates=$1 --maxIterations=$2 --euclideanDistanceMinimum_iteration=1E-10 --numTrainingSequencesPerProfiles=$3

