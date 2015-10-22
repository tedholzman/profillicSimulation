/**
 * \file ProfillicSimulation.cpp
 * \author D'Oleris Paul Thatcher Edlefsen   paul@galosh.org with some additions
 * by Ted Holzman
 * \par Library:
 * Galosh ProfillicSimulation
 * \brief Simulation tests of profillic system.
 * \copyright &copy; 2008, 2011, 2012 by Paul T. Edlefsen, Fred Hutchinson Cancer
 *    Research Center.
 *  All rights reserved.
 ***/
#include "ProfillicSimulation.hpp"

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
  typedef seqan::AminoAcid20 SequenceResidueType;
#else // __PROFUSE_USE_AMINOS .. else
  typedef seqan::Dna ResidueType;
  typedef seqan::Dna SequenceResidueType;
#endif // __PROFUSE_USE_AMINOS .. else ..


  /**
   * arithmetic systems
   */
  typedef doublerealspace ProbabilityType;  /// choices: bfloat, logspace, floatrealspace
  typedef bfloat ScoreType; /// Preferred; otherwise logspace, realspace
  typedef bfloat MatrixValueType;

    ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType> profuse_test(argc,argv);

    // TODO: ERE I AM. Working on getting 2 profiles to work...
    // Amendment: *this is the key: train the globals too!  why?!*
    //   profuse_test.m_parameters.trainProfileGlobals = true;//false;
    //  It seems (I see again/recall) that one key is starting the global params sufficiently low.  For some reason the training works best when indels seem extremely unlikely.

  // TODO: REMOVE
  //profuse_test.m_parameters.debug = DEBUG_All;
  //profuse_test.m_parameters.verbosity = VERBOSITY_High;
  //profuse_test.m_parameters.verbosity = VERBOSITY_Low;
  profuse_test.m_parameters.verbosity = VERBOSITY_Meta;

  // Do it.
  profuse_test.start();
 
// For Saturn, a CHUD profiler that I have on my mac.
#ifdef __HAVE_SATURN
  /* Stop the back-end. */ 
  stopSaturn (); 
#endif /* __HAVE_SATURN */

} // main( int, char** )
