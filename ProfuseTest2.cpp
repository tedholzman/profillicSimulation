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

/**
 * Some options are not directly native to the m_parameters member of profuse_test.
 * ProfuseTest::Parameters inherits from ProfileGibbs::Parameters <-
 *    ProfileTreeTrainer::Parameters <- ProfileTrainer::Parameters <-
 *       ProlificParameters::Parameters <-
 *          DynamicProgramming::Parameters <-
 *             Galosh::Parameters
 *
 * At this point we are only replacing native ProfuseTest parameters with
 * entries in Parameters.m_options_map.  So there are several
 * non-native m_parameters entries
 */
//see:ProfuseTest.cfg  profuse_test.m_parameters.maxIterations = 1000;
//see:ProfuseTest.cfg  profuse_test.m_parameters.maxPositionCycles = 1; // For conditional bw: when not using globals, or just to ensure 1 pos cycle at a time
//see:ProfuseTest.cfg  profuse_test.m_parameters.maxPositionCycles_globals = 1; // For conditional bw: when not using globals, or just to ensure 1 pos cycle at a time

//see:ProfuseTest.cfg  profuse_test.m_parameters.scorePercentChangeMinimum_position_cycle = 1;//.1;
//see:ProfuseTest.cfg  profuse_test.m_parameters.scorePercentChangeMinimum_iteration = .01;

  // Note: Whenever usePriors is false, alwaysAccept should be irrelevant
  // (except for numerical issues, which of course exist): it will always accept
  // because the score is mathematically guaranteed to not go down.
  // Note: When using usePriors and are doing lengthadjust, I recommend turning
  // alwaysAccept *off*.
//see:ProfuseTest.cfg  profuse_test.m_parameters.alwaysAccept = false;//true;
//see:ProfuseTest.cfg  profuse_test.m_parameters.maxBaumWelchInverseScalar = 0;                // Straight-up bw.
//see:ProfuseTest.cfg  profuse_test.m_parameters.maxBaumWelchInverseScalar_globals = 0;        // Straight-up bw.
//see:ProfuseTest.cfg  profuse_test.m_parameters.numIterationsBetweenLengthChanges = 0;
//see:ProfuseTest.cfg  profuse_test.m_parameters.usePriors = false;

  // TODO: REMOVE
  //profuse_test.m_parameters.debug = DEBUG_All;
  //profuse_test.m_parameters.verbosity = VERBOSITY_High;
  //profuse_test.m_parameters.verbosity = VERBOSITY_Low;
  profuse_test.m_parameters.verbosity = VERBOSITY_Meta;

  // TODO: Expose these in ProfileTrainerParameters.hpp, and use ProfuseTest.cfg for them.
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
