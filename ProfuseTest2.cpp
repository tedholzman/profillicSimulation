/**
 * \file ProfuseTest2.cpp
 * \author D'Oleris Paul Thatcher Edlefsen   paul@galosh.org with some additions
 * by Ted Holzman
 * \par Library:
 * Galosh ProfillicSimulation
 * \brief Simulation tests of profillic system.
 * \copyright &copy; 2008, 2011, 2012 by Paul T. Edlefsen, Fred Hutchinson Cancer
 *    Research Center.
 *  All rights reserved.
 ***/
#include "ProfuseTest2.hpp"

using namespace galosh;

#include <seqan/basic.h>

// For Saturn, a CHUD profiler that I have on my mac.
#ifdef __HAVE_SATURN
#include <Saturn.h>
#endif /* __HAVE_SATURN */

/* ////////////
// Required for linking with muscle (even though we don't access it from ProfileTrainer.cpp).
*/
#ifdef __HAVE_MUSCLE
int g_argc;
char **g_argv;
#endif // __HAVE_MUSCLE

// TAH 6/12 macros for defining and accessing commandline parameters.
#include "CommandlineParameters.hpp"

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


  /**
   * arithmetic systems
   */
  typedef doublerealspace ProbabilityType;  /// choices: bfloat, logspace, floatrealspace
  typedef bfloat ScoreType; /// Preferred; otherwise logspace, realspace
  typedef bfloat MatrixValueType;

  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType> profuse_test(argc,argv);

    // TODO: ERE I AM. Working on getting 2 profiles to work...
    // Amendment: *this is the key: train the globals too!  why?!*
    //   profuse_test.m_parameters.trainProfileGlobals = true;//false;
    //  It seems (I see again/recall) that one key is starting the global params sufficiently low.  For some reason the training works best when indels seem extremely unlikely.


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
  
#ifdef ALLOW_BOLTZMANN_GIBBS
//fixyfix
  profuse_test.m_parameters.testBaldi = false;
  profuse_test.m_parameters.testBaldiSiegel = true;
  profuse_test.m_parameters.siegelEpsilonScaleFactor = 1.5;
  profuse_test.m_parameters.siegelMaxRefiningThePeakSteps_positions = 1000;
  profuse_test.m_parameters.siegelRefiningThePeakStepsConvergenceThreshold = 1E-5;
#endif // ALLOW_BOLTZMANN_GIBBS

/**
 * Some options are not directly native to the m_parameters member of profuse_test.
 * ProfuseTest::Parameters inherits from ProfileGibbs::Parameters <-
 *    ProfileTreeTrainer::Parameters <- ProfileTrainer::Parameters <-
 *       ProlificParameters::Parameters <-
 *          DynamicProgramming::Parameters <-
 *             Galosh::Parameters
 *
 * At this point we are only replacing native ProfuseTest parameters with
 * entries in ProfuseTest.m_profusetest_options_map.  So there are several
 * non-native m_parameters entries
 */
  profuse_test.m_parameters.maxIterations = 1000;
  profuse_test.m_parameters.maxPositionCycles = 1; // For conditional bw: when not using globals, or just to ensure 1 pos cycle at a time
  profuse_test.m_parameters.maxPositionCycles_globals = 1; // For conditional bw: when not using globals, or just to ensure 1 pos cycle at a time
  profuse_test.m_parameters.scorePercentChangeMinimum_position_cycle = 1;//.1;
  profuse_test.m_parameters.scorePercentChangeMinimum_iteration = .01;
  // Note: Whenever usePriors is false, alwaysAccept should be irrelevant
  // (except for numerical issues, which of course exist): it will always accept
  // because the score is mathematically guaranteed to not go down.
  // Note: When using usePriors and are doing lengthadjust, I recommend turning
  // alwaysAccept *off*.
  profuse_test.m_parameters.alwaysAccept = false;//true;
  profuse_test.m_parameters.maxBaumWelchInverseScalar = 0;                // Straight-up bw.
  profuse_test.m_parameters.maxBaumWelchInverseScalar_globals = 0;        // Straight-up bw.
  profuse_test.m_parameters.numIterationsBetweenLengthChanges = 0;
  profuse_test.m_parameters.usePriors = false;
  // TODO: REMOVE
  //profuse_test.m_parameters.debug = DEBUG_All;
  //profuse_test.m_parameters.verbosity = VERBOSITY_High;
  //profuse_test.m_parameters.verbosity = VERBOSITY_Low;
  profuse_test.m_parameters.verbosity = VERBOSITY_Meta;

  profuse_test.m_parameters.proposeDeletingThreshold = .5;                //.01#.0125#.025
  profuse_test.m_parameters.proposeDeletingThreshold_increment = .0005;                             //.00005#.005 #.0125 #5E-5#.00625
  profuse_test.m_parameters.proposeInsertingThreshold = .5;                   // = proposeDeletingThreshold
  profuse_test.m_parameters.proposeInsertingThreshold_increment = .0005;      // = proposeDeletingThreshold_increment
  profuse_test.m_parameters.proposeInsertingPreAlignThreshold = 0.5;          // = proposeInsertingThreshold    #.35
  profuse_test.m_parameters.proposeInsertingPostAlignThreshold = 0.5;         // = proposeInsertingThreshold    #.35

  // sampler parameters
  profuse_test.m_parameters.sampleProfilePositions = true;                    // = trainProfilePositions
  profuse_test.m_parameters.sampleProfileGlobals = true;                      // = trainProfileGlobals          #false
  profuse_test.m_parameters.numChains = 4;
  profuse_test.m_parameters.minGibbsIterations = 8000;
  profuse_test.m_parameters.gibbsIterationsIncrement = 8000;                  // = minGibbsIterations
  profuse_test.m_parameters.maxGibbsIterations = 8000;
  profuse_test.m_parameters.burnInFraction = .975;

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
