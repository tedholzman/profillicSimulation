#include "ProfuseTest.hpp"

using namespace galosh;

#include <seqan/basic.h>

// For Saturn, a CHUD profiler that I have on my mac.
#ifdef __HAVE_SATURN
#include <Saturn.h>
#endif /* __HAVE_SATURN */

/////////////
// Required for linking with muscle (even though we don't access it from ProfileTrainer.cpp).
#ifdef __HAVE_MUSCLE
int g_argc;
char **g_argv;
#endif // __HAVE_MUSCLE

int
main ( int argc, char **argv )
{
#ifdef __PROFUSE_USE_AMINOS
  typedef seqan::AminoAcid20 ResidueType;
  typedef seqan::AminoAcid SequenceResidueType;
#else // __PROFUSE_USE_AMINOS .. else
  typedef seqan::Dna ResidueType;
  typedef seqan::Iupac SequenceResidueType;
#endif // __PROFUSE_USE_AMINOS .. else ..

  //typedef bfloat ProbabilityType;
  //typedef logspace ProbabilityType;
  //typedef floatrealspace ProbabilityType; // For speed (vs doublerealspace)
  typedef doublerealspace ProbabilityType;
  
  typedef bfloat ScoreType; // Preferred
  //typedef logspace ScoreType; // SLOWer than bfloat
  //typedef realspace ScoreType; // Only for very few & small sequences
  
  typedef bfloat MatrixValueType;
  //typedef logspace MatrixValueType;
  //typedef doublerealspace MatrixValueType; // For speed (vs bfloat)
  //typedef floatrealspace MatrixValueType;


  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType> profuse_test;

  profuse_test.m_parameters.numTrueProfiles = 1;
  profuse_test.m_parameters.numStartingProfiles = 9;// SEE BELOW at startWithPositionsDrawnFromPrior.

  profuse_test.m_parameters.saveResultsToFile = true;
  profuse_test.m_parameters.saveResultsParentDirectory = "./test_results_hiatus_DNA_edc8edlaplf.0125";

  // TODO: ERE I AM. Working on getting 2 profiles to work...
  profuse_test.m_parameters.numProfiles = 1;//2;
  profuse_test.m_parameters.sharedPositionRate = .5; // For profile trees: the probability for each position that it differs between parent and child.

  // Ammendment: *this is the key: train the globals too!  why?!*
  //   profuse_test.m_parameters.trainProfileGlobals = true;//false;
  //  It seems (I see again/recall) that one key is starting the global params sufficiently low.  For some reason the training works best when indels seem extremely unlikely.

  vector<uint32_t> pl( 1 );
  //pl[ 0 ] = 10;
  //pl[ 1 ] = 50;
  //pl[ 2 ] = 100;
  //pl[ 0 ] = 50;
  pl[ 0 ] = 100;
  //pl[ 0 ] = 25;
  //pl[ 0 ] = 200;
  //pl[ 2 ] = 100;
  //pl[ 2 ] = 200;
  profuse_test.m_parameters.profileLengths = &pl;

  vector<uint32_t> tspp( 1 );
  //tspp[ 0 ] = 100;
  //tspp[ 0 ] = 10;
  //tspp[ 1 ] = 25;
  //tspp[ 2 ] = 50;
  //tspp[ 0 ] = 50;
  tspp[ 0 ] = 100;
  //tspp[ 1 ] = 200;
  profuse_test.m_parameters.numTrainingSequencesPerProfiles = &tspp;

  profuse_test.m_parameters.numTestingSequencesPerProfile = 100;//profuse_test.m_parameters.numTrainingSequencesPerProfile;

  vector<double> cr( 7 );
  //cr[ 0 ] = .5;

  //cr[ 0 ] = .3;
  //cr[ 1 ] = .5;
  //cr[ 2 ] = .7;
  //cr[ 3 ] = .9;

  cr[ 0 ] = .3;
  cr[ 1 ] = .4;
  cr[ 2 ] = .5;
  cr[ 3 ] = .6;
  cr[ 4 ] = .7;
  cr[ 5 ] = .8;
  cr[ 6 ] = .9;

  //cr[ 0 ] = .1;
  //cr[ 1 ] = .2;
  //cr[ 2 ] = .3;
  //cr[ 3 ] = .4;
  //cr[ 4 ] = .5;
  //cr[ 5 ] = .6;
  //cr[ 6 ] = .7;
  //cr[ 7 ] = .8;
  //cr[ 8 ] = .9;

  ////cr[ 3 ] = .6;
  profuse_test.m_parameters.conservationRates = &cr;

  // Make the expected number of deletions be .5 or 1.0 per sequence.
  // Note that this will apply to the insertions, too, unless
  // m_parameters.useDeletionsForInsertionsParameters is set to false.
  vector<double> edc( 1 );
  //edc[ 0 ] = 1.0;
  //edc[ 0 ] = 4.0;
  edc[ 0 ] = 8.0;
  //edc[ 1 ] = 1.5;
  //edc[ 1 ] = 1.0;
  //edc[ 2 ] = 1.5;
  //edc[ 0 ] = .5;
  //edc[ 1 ] = .75;
  //edc[ 2 ] = 1.0;
  //edc[ 0 ] = .1;
  ////edc[ 1 ] = .5;
  //edc[ 1 ] = .3;
  //edc[ 2 ] = .5;
  //edc[ 3 ] = .7;
  //edc[ 4 ] = .9;
  profuse_test.m_parameters.expectedDeletionsCounts = &edc;

  // Make the expected number of insertions be .5 or 1.0 per sequence.
  // Note that this is not used unless
  // m_parameters.useDeletionsForInsertionsParameters is set to false.
//  vector<double> eic( 1 );
//  eic[ 0 ] = .5;
//  //eic[ 1 ] = 1.0;
//  profuse_test.m_parameters.expectedInsertionsCounts = &eic;

  // Make the expected length of each deletion be ( profile_length / 20 ) or (
  // profile_length / 10 )...
  // Note that this will apply to the insertions, too, unless
  // m_parameters.useDeletionsForInsertionsParameters is set to false.
  vector<double> edlaplf( 1 );
  //edlaplf[ 0 ] = 0.1;
  //edlaplf[ 0 ] = .025;
  edlaplf[ 0 ] = .0125;
  //edlaplf[ 1 ] = .15;
  //edlaplf[ 1 ] = .10;
  //edlaplf[ 2 ] = .15;
  profuse_test.m_parameters.expectedDeletionLengthAsProfileLengthFractions = &edlaplf;

  // Make the expected length of each insertion be ( profile_length / 20 ) or (
  // profile_length / 10 )...
  // Note that this is not used unless
  // m_parameters.useDeletionsForInsertionsParameters is set to false.
//  vector<double> eilaplf( 1 );
//  eilaplf[ 0 ] = .05;
//  //eilaplf[ 1 ] = .1;
//  profuse_test.m_parameters.expectedInsertionLengthAsProfileLengthFractions = &eilaplf;

  // ..(or 1.25, whichever is larger).
  profuse_test.m_parameters.minExpectedDeletionLength = 1.25;
  profuse_test.m_parameters.minExpectedInsertionLength = 1.25;

  profuse_test.m_parameters.preAlignInsertion = .05;
  profuse_test.m_parameters.postAlignInsertion = .05;

  // trainer Parameters
  profuse_test.m_parameters.trainProfilePositions = true;
  profuse_test.m_parameters.trainProfileGlobals = true;
  profuse_test.m_parameters.startWithUniformGlobals = // see startWithGlobalsDrawnFromPrior below too
    true;
  profuse_test.m_parameters.startWithUniformGlobals_scalar = 10.0;
  profuse_test.m_parameters.startWithUniformGlobals_maxNtoN = .2;
  profuse_test.m_parameters.startWithUniformGlobals_maxBtoD = .2;
  profuse_test.m_parameters.startWithUniformGlobals_maxMtoI = .05;//.2;
  profuse_test.m_parameters.startWithUniformGlobals_maxMtoD = .05;//.2;
  profuse_test.m_parameters.startWithUniformGlobals_maxItoI = .5;
  profuse_test.m_parameters.startWithUniformGlobals_maxDtoD = .5;
  profuse_test.m_parameters.startWithUniformGlobals_maxCtoC = .2;
  profuse_test.m_parameters.startWithUniformPositions = // see startWithPositionsDrawnFromPrior below too
    true;
  profuse_test.m_parameters.maxIterations = 1000;
  profuse_test.m_parameters.maxPositionCycles = 1; // For conditional bw: when not using globals, or just to ensure 1 pos cycle at a time
  profuse_test.m_parameters.maxPositionCycles_globals = 1; // For conditional bw: when not using globals, or just to ensure 1 pos cycle at a time
  profuse_test.m_parameters.scorePercentChangeMinimum_position_cycle = 1;//.1;
  profuse_test.m_parameters.scorePercentChangeMinimum_iteration = .01;
  // TODO: REMOVE
  //profuse_test.m_parameters.euclideanDistanceMinimum_iteration = 5E-7; // 5E-7 is great, though slow...

  // Note: Whenever usePriors is false, alwaysAccept should be irrelevant
  // (except for numerical issues, which of course exist): it will always accept
  // because the score is mathematically guaranteed to not go down.
  // Note: When using usePriors and are doing lengthadjust, I recommend turning
  // alwaysAccept *off*.
  profuse_test.m_parameters.alwaysAccept = false;//true;

  // TODO: REMOVE?
  profuse_test.m_parameters.maxBaumWelchInverseScalar = 0; // Straight-up bw.
  profuse_test.m_parameters.maxBaumWelchInverseScalar_globals = 0; // Straight-up bw.
  //profuse_test.m_parameters.minBaumWelchInverseScalar = 4.0;

  // Lengthadjust trainer parameters
  profuse_test.m_parameters.testLengthadjust = true;
  profuse_test.m_parameters.useAlignmentProfiles = true; // For now I think this is necessary.
  profuse_test.m_parameters.numIterationsBetweenLengthChanges = 0;
  profuse_test.m_parameters.proposeDeletingThreshold = .5;//.01;//.0125;//.025;
  profuse_test.m_parameters.proposeDeletingThreshold_increment = .0005;//.00005;//.005; //.0125; //5E-5;//.00625;
  profuse_test.m_parameters.proposeInsertingThreshold =
    profuse_test.m_parameters.proposeDeletingThreshold;// / seqan::ValueSize<ResidueType>::VALUE; // TODO: ERE I AM. PUT BACK?
  profuse_test.m_parameters.proposeInsertingThreshold_increment =
    profuse_test.m_parameters.proposeDeletingThreshold_increment; // TODO: ERE I AM.  PUT BACK? / seqan::ValueSize<ResidueType>::VALUE;
  profuse_test.m_parameters.proposeInsertingPreAlignThreshold = //.35;
    profuse_test.m_parameters.proposeInsertingThreshold;
  profuse_test.m_parameters.proposeInsertingPostAlignThreshold = //.35;
    profuse_test.m_parameters.proposeInsertingThreshold;

  /// sampler parameters
  profuse_test.m_parameters.sampleProfilePositions =
    profuse_test.m_parameters.trainProfilePositions;//true;
  profuse_test.m_parameters.sampleProfileGlobals =
    profuse_test.m_parameters.trainProfileGlobals;//false;
  profuse_test.m_parameters.numChains = 4;
  profuse_test.m_parameters.minGibbsIterations = 8000;
  profuse_test.m_parameters.gibbsIterationsIncrement =
    profuse_test.m_parameters.minGibbsIterations;
  profuse_test.m_parameters.maxGibbsIterations = 8000;
  profuse_test.m_parameters.burnInFraction = .975;
  profuse_test.m_parameters.reportGibbsMean = true;
  profuse_test.m_parameters.reportGibbsMode = false;

  profuse_test.m_parameters.usePriors = false;//true;
  // For (only) an even distribution, set startWithPositionsDrawnFromPrior = false here and also use startWithUniformPositions == false (see above).  If you do that, remember also to set numStartingProfiles = 1 (above).  But you could use alsoStartWithEvenPositions to make the 0th startingProfile even().
  // Note that we now allow both startWithUniformPositions and startWithPositionsDrawnFromPrior, in which case the first half of the random starting profiles will have positions drawn from the prior and the second half of the random starting profiles will have uniform() positions.  Note that if alsoStartWithEvenPositions, then the very first starting profile (index 0) will have even() positions.
  profuse_test.m_parameters.startWithPositionsDrawnFromPrior = true;//false;
  profuse_test.m_parameters.alsoStartWithEvenPositions = true;
  profuse_test.m_parameters.startWithGlobalsDrawnFromPrior = false;
  profuse_test.m_parameters.priorStrength = 100;
  profuse_test.m_parameters.priorStrength_internal_transitions = 10;

  // Tests to run
  profuse_test.m_parameters.coutTruepath = true;
  profuse_test.m_parameters.testTrueProfile = true;
  profuse_test.m_parameters.coutTrueProfile = true;
  profuse_test.m_parameters.testStartingProfile = true;
  profuse_test.m_parameters.coutStartingProfile = true;
  profuse_test.m_parameters.testUnconditionalProfile = true;
  profuse_test.m_parameters.coutUnconditionalProfile = true;
  profuse_test.m_parameters.testUnconditionalWithFixedStartingGlobalsProfile = false;
  profuse_test.m_parameters.coutUnconditionalWithFixedStartingGlobalsProfile = false;
  profuse_test.m_parameters.testUnconditionalWithFixedTrueGlobalsProfile = false;
  profuse_test.m_parameters.coutUnconditionalWithFixedTrueGlobalsProfile = false;
  profuse_test.m_parameters.testConditionalThenUnconditionalProfile = false;
  profuse_test.m_parameters.coutConditionalThenUnconditionalProfile = false;
  profuse_test.m_parameters.testUnconditionalThenConditionalProfile = false;
  profuse_test.m_parameters.coutUnconditionalThenConditionalProfile = false;
  profuse_test.m_parameters.testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile = false;
  profuse_test.m_parameters.coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile = false;
  profuse_test.m_parameters.testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile = false;
  profuse_test.m_parameters.coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile = false;
  profuse_test.m_parameters.testConditionalGibbsProfile = false;//true;
  profuse_test.m_parameters.coutConditionalGibbsProfile = false;//true;
  profuse_test.m_parameters.testUnconditionalGibbsProfile = false;//true;
  profuse_test.m_parameters.coutUnconditionalGibbsProfile = false;//true;

#ifdef ALLOW_BOLTZMANN_GIBBS
  profuse_test.m_parameters.testBaldi = false;
  profuse_test.m_parameters.testBaldiSiegel = true;
  profuse_test.m_parameters.siegelEpsilonScaleFactor = 1.5;
  profuse_test.m_parameters.siegelMaxRefiningThePeakSteps_positions = 1000;
  profuse_test.m_parameters.siegelRefiningThePeakStepsConvergenceThreshold = 1E-5;
#else // !ALLOW_BOLTZMANN_GIBBS
  profuse_test.m_parameters.testBaldi = false;
  profuse_test.m_parameters.testBaldiSiegel = false;
#endif // ALLOW_BOLTZMANN_GIBBS

  // TODO: REMOVE
  //profuse_test.m_parameters.debug = DEBUG_All;
  //profuse_test.m_parameters.verbosity = VERBOSITY_High;
  //profuse_test.m_parameters.verbosity = VERBOSITY_Low;
  profuse_test.m_parameters.verbosity = VERBOSITY_Meta;
        
  // TODO: REMOVE
  profuse_test.m_parameters.useRabinerScaling = false;

  // TODO: REMOVE
  //profuse_test.m_random.setSeed( 1304905685 );

  // Do it.
  profuse_test.start();
 
// For Saturn, a CHUD profiler that I have on my mac.
#ifdef __HAVE_SATURN
  /* Stop the back-end. */ 
  stopSaturn (); 
#endif /* __HAVE_SATURN */

} // main( int, char** )
