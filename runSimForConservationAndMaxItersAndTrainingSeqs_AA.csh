#!/usr/bin/tcsh
time bin/darwin-4.2.1/release/profillicSimulation_AA --seed=98103 --saveResultsParentDirectory=AA_runSimForConservation_$1_andMaxIters_$2_AndTrainingSeqs_$3 --conservationRates=$1 --maxIterations=$2 --euclideanDistanceMinimum_iteration=1E-10 --numTrainingSequencesPerProfiles=$3

