/**
 * \file ProfuseTest2.hpp
 * \author D'Oleris Paul Thatcher Edlefsen   paul@galosh.org with some serious modifications
 * by Ted Holzman
 * \par Library:
 * Galosh ProfillicSimulation
 * \brief Simulation tests of profuse system.
 * \copyright &copy; 2008, 2011, 2012 by Paul T. Edlefsen, Fred Hutchinson Cancer
 *    Research Center.
 *  All rights reserved.
 *****************************************************************************/

#if     _MSC_VER > 1000
#pragma once
#endif

#ifndef __GALOSH_PROFUSETEST_HPP__
#define __GALOSH_PROFUSETEST_HPP__

#include "Parameters.hpp"
using galosh::Parameters;
using galosh::DebugLevel;
using galosh::VerbosityLevel;

#include "Profile.hpp"
using galosh::ProfileTreeRoot;
using galosh::ProfileTreeInternalNode;

#include "ProfileTree.hpp"
using galosh::ProfileTree;

#include "ProfileTreeTrainer.hpp"
using galosh::ProfileTreeTrainer;

#include "ProfileGibbs.hpp"
using galosh::ProfileGibbs;

#include "ProfileTrainer.hpp"
using galosh::ProfileTrainer;

#include "Sequence.hpp"
using galosh::Sequence;

#include "Random.hpp"
using galosh::Random;

#include "DynamicProgramming.hpp"
using galosh::DynamicProgramming;

#include <string>
using std::string;
#include <iostream>
using std::cout;
using std::endl;
#include <ctime>
using std::clock_t;
using std::clock;

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include "boost/filesystem.hpp"

#include <boost/bind.hpp>
//#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/join.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;
#include "CommandlineParameters.hpp"
/// TAH 6/12 predefining DEFAULT_OPTIONS_DESCRIPION and DEFAULT_VARIABLES_MAP will
/// replace those values in the CommandlineParameters macros
#define DEFAULT_OPTIONS_DESCRIPTION m_profusetest_options
#define DEFAULT_VARIABLES_MAP       m_parameters.m_options_map

namespace galosh {

template <class ResidueType,
          class ProbabilityType,
          class ScoreType,
          class MatrixValueType,
          class SequenceResidueType>
  class ProfuseTest {
  public:
    typedef Sequence<SequenceResidueType> SequenceType;

    class Parameters :
       public ProfileGibbs<ProfileTreeRoot<ResidueType, ProbabilityType>,ScoreType,MatrixValueType,SequenceResidueType>::Parameters
    {
    public:
       po::options_description m_profusetest_options;

    private:
      typedef typename ProfileGibbs<ProfileTreeRoot<ResidueType, ProbabilityType>,ScoreType,MatrixValueType,SequenceResidueType>::Parameters profile_gibbs_parameters_t;
      // Boost serialization
      friend class boost::serialization::access;
      template<class Archive>
      void serialize ( Archive & ar, const unsigned int /* file_version */ )
      {
        // save/load base class information
//        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( profile_gibbs_parameters_t );
//        ar & m_options_map;
/**
 * \note This is a hack based on Paul's suggestion that the real use for the serialization
 * is for making a permanent file copy of parameter values, not for deserialization.
 * The "correct" way to do this would be to serialize (and potentially deserialize) the
 * the variables_map object (m_options_map) within the Parameters object.
 * Unfortunately, the variables_map object ultimately stores objects of type boost:any,
 * which are not serializable in the general case. However, since we tend to store only
 * doubles, unsigned ints and strings, it may be serializable for these particular cases.
 * That counts as a \todo for the current time.
 */

//#define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) \
//        ar & boost::serialization::make_nvp(#NAME, ((galosh::Parameters)(*this)).m_options_map[#NAME].as<TYPE>())
//        #include "ProfuseTestOptions.hpp"
//#undef GALOSH_DEF_OPT
      } // serialize Parameters

    public:
  
      /// PARAMETERS

      // TODO: REMOVE. TESTING.
      string configFile;
      uint32_t seed;
      bool useRabinerScaling;

      bool saveResultsToFile;
  #define DEFAULT_saveResultsToFile true

      string saveResultsParentDirectory;
  #define DEFAULT_saveResultsParentDirectory "."

      string resultsFilePrefix;
  #define DEFAULT_resultsFilePrefix "ProfuseTest."

      string tabFileSuffix;
  #define DEFAULT_tabFileSuffix ".tab"

      string parametersFileSuffix;
  #define DEFAULT_parametersFileSuffix ".Parameters.xml"

      bool saveTrueProfileTrees;
  #define DEFAULT_saveTrueProfileTrees true

      string trueProfileTreeFileSuffix;
  #define DEFAULT_trueProfileTreeFileSuffix ".true.ProfileTree.xml"

      bool saveStartingProfiles;
  #define DEFAULT_saveStartingProfiles true

      string startingProfileTreeFileSuffix;
  #define DEFAULT_startingProfileTreeFileSuffix ".starting.ProfileTree.xml"

      bool saveTestProfiles;
  #define DEFAULT_saveTestProfiles true

      string testProfileTreeFileSuffix;
  #define DEFAULT_testProfileTreeFileSuffix ".ProfileTree.xml"

      bool savePatternSequences;
  #define DEFAULT_savePatternSequences true

      string patternSequencesFileSuffix;
  #define DEFAULT_patternSequencesFileSuffix ".pattern_sequences.fasta"

      bool saveTests;
  #define DEFAULT_saveTests false

      string testsFileSuffix;
  #define DEFAULT_testsFileSuffix ".tests"

      bool saveTrainingSequences;
  #define DEFAULT_saveTrainingSequences true

      string trainingSequencesFileSuffix;
  #define DEFAULT_trainingSequencesFileSuffix ".training_sequences.fasta"

      bool saveTestingSequences;
  #define DEFAULT_saveTestingSequences true

      string testingSequencesFileSuffix;
  #define DEFAULT_testingSequencesFileSuffix ".testing_sequences.fasta"

      bool saveTrueTrainingAlignments;
  #define DEFAULT_saveTrueTrainingAlignments true

      string trainingTrueAlignmentsFileSuffix;
  #define DEFAULT_trainingTrueAlignmentsFileSuffix ".training_alignments.fasta"

      bool saveTrueTestingAlignments;
  #define DEFAULT_saveTrueTestingAlignments true

      string trueTestingAlignmentsFileSuffix;
  #define DEFAULT_trueTestingAlignmentsFileSuffix ".testing_alignments.fasta"


      // 1 was the beginning
      // 2 is after I modified the conditional_then_unconditional_root stuff to use the globals from the conditional, but the position-specific values from the starting profile.
      // 3 is after I added unconditional_with_fixed_starting_globals
      // 4 is after I added unconditional_with_fixed_starting_globals_then_with_fixed_positions
      // 5 is after I added startWithUniformGlobals
      // 6 is after I added startWithUniformGlobals_scalar
      // 7 is after I added cout[Test] params
      // 8 is after I added startWithUniformPositions, startWithPositionsDrawnFromPrior, etc.
      // 9 is after I added convert_tab_output_to_log_double
      // 10 is after I added CPU time
      uint32_t saveFileVersion;
  #define DEFAULT_saveFileVersion 10

      uint32_t numProfiles;
  #define DEFAULT_numProfiles 1

      /**
       * When making the pattern sequence (and the true root profile from it),
       * use this length.
       *
       */
      vector<uint32_t> profileLengths;
#define DEFAULT_profileLengths vector<uint32_t>()

      /**
       * When numProfiles > 1, what fraction of the positions of each child
       * should be shared with its parent?
       */
      double sharedPositionRate;
  #define DEFAULT_sharedPositionRate 1.0

      /**
       * Use this number of sequences when training.
       *
       */
      vector<uint32_t> numTrainingSequencesPerProfiles;
#define DEFAULT_numTrainingSequencesPerProfiles vector<uint32_t>()

      uint32_t numTestingSequencesPerProfile;
  #define DEFAULT_numTestingSequencesPerProfile 20

      /**
       * When making the true root profile from the pattern sequence, use this
       * probability for the pattern sequence base at each position, and divide
       * the remaining probability evenly among the remaining bases.
       *
       */
      vector<double> conservationRates;
#define DEFAULT_conservationRates vector<double>()

      /**
       * Lock the indel parameters of the true profile to be the same for
       * insertions as for deletions?  This makes the expectedInsertionsCounts,
       * expectedInsertionLengthAsProfileLengthFractions, and
       * minExpectedInsertionLength unused, since the corresponding deletion
       * values will be used instead.  It also reduces the number of tests by
       * reducing the number of possible combinations (since deletions and
       * insertions will go in lock step).
       */
      bool useDeletionsForInsertionsParameters;
  #define DEFAULT_useDeletionsForInsertionsParameters true

      /**
       * The deletionOpen value of the true profile will be set to (
       * expectedDeletionsCount / profileLength ).  If
       * useDeletionsForInsertionsParameters is true, the insertionOpen value
       * of the true profile will also be set to ( expectedDeletionsCount /
       * profileLength ).
       *
       * @see useDeletionsForInsertionsParameters
       */
      vector<double> expectedDeletionsCounts;
#define DEFAULT_expectedDeletionsCounts vector<double>()

      /**
       * If useDeletionsForInsertionsParameters is false, the insertionOpen
       * value of the true profile will be set to ( expectedInsertionsCount /
       * profileLength ).
       *
       * @see useDeletionsForInsertionsParameters
       */
      vector<double> expectedInsertionsCounts;
  #define DEFAULT_expectedInsertionsCounts vector<double>()

      /**
       * The deletionExtension value of the true profile will be the minimum of
       * ( 1.0 / ( expectedDeletionLengthAsProfileLengthFraction *
       * profileLength ) ) and ( 1.0 / minExpectedDeletionLength ).  If
       * useDeletionsForInsertionsParameters is true, the insertionExtension
       * value of the true profile will also be set to be the minimum of ( 1.0
       * / ( expectedDeletionLengthAsProfileLengthFraction * profileLength ) )
       * and ( 1.0 / minExpectedDeletionLength ).
       *
       * @see useDeletionsForInsertionsParameters
       */
      vector<double> expectedDeletionLengthAsProfileLengthFractions;
  #define DEFAULT_expectedDeletionLengthAsProfileLengthFractions vector<double>()

      /**
       * If useDeletionsForInsertionsParameters is false, the
       * insertionExtension value of the true profile will be the minimum of (
       * 1.0 / ( expectedInsertionLengthAsProfileLengthFraction * profileLength
       * ) ) and ( 1.0 / minExpectedInsertionLength ).
       *
       * @see useDeletionsForInsertionsParameters
       */
      vector<double> expectedInsertionLengthAsProfileLengthFractions;
  #define DEFAULT_expectedInsertionLengthAsProfileLengthFractions vector<double>()

      /**
       * The deletionExtension value of the true profile will be the minimum of
       * ( 1.0 / ( expectedDeletionLengthAsProfileLengthFraction *
       * profileLength ) ) and ( 1.0 / minExpectedDeletionLength ).  If
       * useDeletionsForInsertionsParameters is true, the insertionExtension
       * value of the true profile will also be the minimum of ( 1.0 / (
       * expectedDeletionLengthAsProfileLengthFraction * profileLength ) ) and
       * ( 1.0 / minExpectedDeletionLength ).
       *
       * @see useDeletionsForInsertionsParameters
       */
      double minExpectedDeletionLength;
  #define DEFAULT_minExpectedDeletionLength 1.25

      /**
       * If useDeletionsForInsertionsParameters is false, the
       * insertionExtension value of the true profile will be the minimum of (
       * 1.0 / ( expectedInsertionLengthAsProfileLengthFraction * profileLength
       * ) ) and ( 1.0 / minExpectedInsertionLength ).
       *
       * @see useDeletionsForInsertionsParameters
       */
      double minExpectedInsertionLength;
  #define DEFAULT_minExpectedInsertionLength 1.25

      /**
       * The preAlignInsertion value of the true profile.
       */
      double preAlignInsertion;
  #define DEFAULT_preAlignInsertion .01

      /**
       * The postAlignInsertion value of the true profile.
       */
      double postAlignInsertion;
  #define DEFAULT_postAlignInsertion .01

      /**
       * The effective number of sequences "observed" a priori.  Note that we
       * use a different prior strength for main-model transitions: see
       * priorStrength_internal_transitions.
       */
      float priorStrength;
  #define DEFAULT_priorStrength 1.0f

      /**
       * The effective number of sequences "observed" a priori, for main-model
       * transitions.
       */
      float priorStrength_internal_transitions;
  #define DEFAULT_priorStrength_internal_transitions 10.0f

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * M->M transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      float priorMtoM;
  #define DEFAULT_priorMtoM .95f

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * M->I transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      float priorMtoI;
  #define DEFAULT_priorMtoI .025f

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * M->D transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      float priorMtoD;
  #define DEFAULT_priorMtoD .025f

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * I->M transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      float priorItoM;
  #define DEFAULT_priorItoM .05f

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * I->I transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      float priorItoI;
  #define DEFAULT_priorItoI .95f

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * D->M transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      float priorDtoM;
  #define DEFAULT_priorDtoM .95f

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * D->D transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      float priorDtoD;
  #define DEFAULT_priorDtoD .05f

      /**
       * Additionally report the overall mean of all chains found while
       * performing Gibbs sampling?  The best profile is always reported, which
       * may be the overall mean, the mode, or the mean of one of the chains.
       */
      bool reportGibbsMean;
  #define DEFAULT_reportGibbsMean false

      /**
       * Additionally report the mode found while
       * performing Gibbs sampling?  The best profile is always reported, which
       * may be the overall mean, the mode, or the mean of one of the chains.
       *
       * Note that it takes some extra time to store the mode (this turns on
       * saveGibbsMode in the ProfileGibbs class).
       */
      bool reportGibbsMode;
  #define DEFAULT_reportGibbsMode false

      /**
       * We do the whole thing a number of different times, starting over with
       * a new pattern sequence.
       */
      uint32_t numTrueProfiles;
  #define DEFAULT_numTrueProfiles 4

      /**
       * For each true root profile, we run the trainers from a number of
       * different starting profiles.  This is that number.
       */
      uint32_t numStartingProfiles;
  #define DEFAULT_numStartingProfiles 4

      /**
       * If startWithGlobalsDrawnFromPrior is not true, and
       * if startWithUniformGlobals is true, then we set the global values of
       * the startingProfile to random values between 0 and
       * min(startWithUniformGlobals_scalar times the true
       * values,startWithUniformGlobals_maxXtoY).  If it is false, we start
       * with the known, true globals.
       *
       * @see startWithUniformGlobals_scalar
       */
      bool startWithUniformGlobals;
  #define DEFAULT_startWithUniformGlobals false

      /**
       * @see startWithUniformGlobals
       */
      double startWithUniformGlobals_scalar;
  #define DEFAULT_startWithUniformGlobals_scalar 2.0

      /**
       * @see startWithUniformGlobals
       */
      double startWithUniformGlobals_maxNtoN;
  #define DEFAULT_startWithUniformGlobals_maxNtoN .2

      /**
       * @see startWithUniformGlobals
       */
      double startWithUniformGlobals_maxBtoD;
  #define DEFAULT_startWithUniformGlobals_maxBtoD .2

      /**
       * @see startWithUniformGlobals
       */
      double startWithUniformGlobals_maxMtoI;
  #define DEFAULT_startWithUniformGlobals_maxMtoI .2

      /**
       * @see startWithUniformGlobals
       */
      double startWithUniformGlobals_maxMtoD;
  #define DEFAULT_startWithUniformGlobals_maxMtoD .2

      /**
       * @see startWithUniformGlobals
       */
      double startWithUniformGlobals_maxItoI;
  #define DEFAULT_startWithUniformGlobals_maxItoI .5

      /**
       * @see startWithUniformGlobals
       */
      double startWithUniformGlobals_maxDtoD;
  #define DEFAULT_startWithUniformGlobals_maxDtoD .5

      /**
       * @see startWithUniformGlobals
       */
      double startWithUniformGlobals_maxCtoC;
  #define DEFAULT_startWithUniformGlobals_maxCtoC .2

      /**
       * If startWithUniformPositions is true, then we set the
       * position-specific values of the startingProfile to random values
       * between 0 and 1.  If it is false, we start with the known, true
       * parameter values.  Note that if startWithPositionsDrawnFromPrior is
       * also true, then the first half of the starting profiles will start
       * with positions drawn from the prior and the second half will start
       * with uniform() positions (possibly excluding the index-0 starting
       * profile, if alsoStartWithEvenPositions is true).
       *
       * @see startWithPositionsDrawnFromPrior
       * @see alsoStartWithEvenPositions
       */
      bool startWithUniformPositions;
  #define DEFAULT_startWithUniformPositions false

      /**
       * If startWithGlobalsDrawnFromPrior is true, the
       * global values of the starting profile will be drawn from the prior.
       *
       * @see startWithUniformGlobals
       */
      bool startWithGlobalsDrawnFromPrior;
  #define DEFAULT_startWithGlobalsDrawnFromPrior false

      /**
       * If startWithPositionsDrawnFromPrior is true, the
       * position-specific values of the starting profile will be drawn from
       * the prior... but see the notes in startWithUniformPositions.
       *
       * @see startWithUniformPositions
       * @see alsoStartWithEvenPositions
       */
      bool startWithPositionsDrawnFromPrior;
  #define DEFAULT_startWithPositionsDrawnFromPrior false

      /**
       * Calculate the viterbi scores after training each profile?  Note that
       * this is in addition to the forward scores, which we always calculate.
       */
      bool testViterbi;
  #define DEFAULT_testViterbi true

      /**
       * Write the viterbi scores to STDOUT?
       */
      bool coutViterbi;
  #define DEFAULT_coutViterbi false

      /**
       * Calculate the truepath scores after training each profile?  Note that
       * this is in addition to the forward scores, which we always calculate.
       */
      bool testTruepath;
  #define DEFAULT_testTruepath true

      /**
       * Write the truepath scores to STDOUT?
       */
      bool coutTruepath;
  #define DEFAULT_coutTruepath false

      /**
       * Calculate the SKL distance between the training profiles and the
       * true profile?
       */
      bool calculateSymmeterizedKullbackLeiblerDistancesToTrue;
  #define DEFAULT_calculateSymmeterizedKullbackLeiblerDistancesToTrue true

      /**
       * Calculate the SKL distance between the training profiles and the
       * starting profile?
       */
      bool calculateSymmeterizedKullbackLeiblerDistancesToStarting;
  #define DEFAULT_calculateSymmeterizedKullbackLeiblerDistancesToStarting false

      /**
       * Calculate the SKL distance between the training profiles and the
       * conditional profile?
       */
      bool coutDistances;
  #define DEFAULT_coutDistances true

      /**
       * Calculate SKL profile-profile alignments between each trained profile
       * and the true profile?
       */
      bool calculateProfileProfileAlignments;
  #define DEFAULT_calculateProfileProfileAlignments true

      /**
       * The cost of a gap open when performing SKL profile-profile alignements.
       */
      double profileProfileIndelOpenCost;
  #define DEFAULT_profileProfileIndelOpenCost .25

      /**
       * The cost of a gap extension when performing SKL profile-profile alignements.
       */
      double profileProfileIndelExtensionCost;
  #define DEFAULT_profileProfileIndelExtensionCost .25

      /**
       * Calculate forward (etc) scores for the true profile, too?
       */
      bool testTrueProfile;
  #define DEFAULT_testTrueProfile true

      /**
       * Write forward (etc) scores for the true profile to STDOUT
       * during training?
       */
      bool coutTrueProfile;
  #define DEFAULT_coutTrueProfile true

      /**
       * Calculate forward (etc) scores for the pre-training profile, too?
       */
      bool testStartingProfile;
  #define DEFAULT_testStartingProfile true

      /**
       * Write forward (etc) scores for the pre-training profile to STDOUT
       * during training?
       */
      bool coutStartingProfile;
  #define DEFAULT_coutStartingProfile true

      /**
       * Also train the unconditional profile, and calculate forward (etc)
       * scores using it?
       */
      bool testUnconditionalProfile;
  #define DEFAULT_testUnconditionalProfile true

      /**
       * Write the forward (etc) scores for the unconditional profile to STDOUT
       * during training?
       */
      bool coutUnconditionalProfile;
  #define DEFAULT_coutUnconditionalProfile true

      /**
       * Also train the unconditional (with fixed starting globals) profile,
       * and calculate forward (etc) scores using it?
       */
      bool testUnconditionalWithFixedStartingGlobalsProfile;
  #define DEFAULT_testUnconditionalWithFixedStartingGlobalsProfile false

      /**
       * Write the forward (etc) scores for the unconditional (with fixed
       * starting globals) profile to STDOUT during training?
       */
      bool coutUnconditionalWithFixedStartingGlobalsProfile;
  #define DEFAULT_coutUnconditionalWithFixedStartingGlobalsProfile false

      /**
       * Also train the unconditional (with fixed true globals) profile,
       * and calculate forward (etc) scores using it?
       */
      bool testUnconditionalWithFixedTrueGlobalsProfile;
  #define DEFAULT_testUnconditionalWithFixedTrueGlobalsProfile false

      /**
       * Write the forward (etc) scores for the unconditional (with fixed
       * true globals) profile to STDOUT during training?
       */
      bool coutUnconditionalWithFixedTrueGlobalsProfile;
  #define DEFAULT_coutUnconditionalWithFixedTrueGlobalsProfile false

      /**
       * Also train the "conditional, then unconditional" profile, and
       * calculate forward (etc) scores using it?
       */
      bool testConditionalThenUnconditionalProfile;
  #define DEFAULT_testConditionalThenUnconditionalProfile false

      /**
       * Write the forward (etc) scores for the "conditional, then
       * unconditional" profile to STDOUT during training?
       */
      bool coutConditionalThenUnconditionalProfile;
  #define DEFAULT_coutConditionalThenUnconditionalProfile false

      /**
       * Also train the "unconditional, then conditional" profile, and
       * calculate forward (etc) scores using it?
       */
      bool testUnconditionalThenConditionalProfile;
  #define DEFAULT_testUnconditionalThenConditionalProfile false

      /**
       * Write the forward (etc) scores for the "unconditional, then
       * conditional" profile to STDOUT during training?
       */
      bool coutUnconditionalThenConditionalProfile;
  #define DEFAULT_coutUnconditionalThenConditionalProfile false

      /**
       * Also train the "unconditional (first with fixed starting globals, then
       * with fixed positions)" profile, and calculate forward (etc) scores
       * using it?
       */
      bool testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile;
  #define DEFAULT_testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile false

      /**
       * Write the forward (etc) scores for the unconditional (first with fixed
       * starting globals, then with fixed positions) profile to STDOUT during
       * training?
       */
      bool coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile;
  #define DEFAULT_coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile false

      /**
       * Also train the "unconditional (first with fixed true globals, then
       * with fixed positions)" profile, and calculate forward (etc) scores
       * using it?
       */
      bool testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile;
  #define DEFAULT_testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile false

      /**
       * Write the forward (etc) scores for the unconditional (first with fixed
       * true globals, then with fixed positions) profile to STDOUT during
       * training?
       */
      bool coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile;
  #define DEFAULT_coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile false

      /**
       * Also use the best profile found after Gibbs sampling using Conditional
       * ("Per Position") Gibbs and calculate forward (etc) scores using it?
       */
      bool testConditionalGibbsProfile;
  #define DEFAULT_testConditionalGibbsProfile false

      /**
       * Write the forward (etc) scores for the best profile found (for
       * Conditional Gibbs) to STDOUT during training?
       */
      bool coutConditionalGibbsProfile;
  #define DEFAULT_coutConditionalGibbsProfile true

      /**
       * Also use the best profile found after Gibbs sampling using
       * Unconditional ("Simple") Gibbs and calculate forward (etc) scores
       * using it?
       */
      bool testUnconditionalGibbsProfile;
  #define DEFAULT_testUnconditionalGibbsProfile false

      /**
       * Write the forward (etc) scores for the best profile found profile (for
       * Unconditional Gibbs) to STDOUT during training?
       */
      bool coutUnconditionalGibbsProfile;
  #define DEFAULT_coutUnconditionalGibbsProfile true

      /**
       * Also use lengthadjust to train the profiles, and calculate forward (etc) scores using it?
       */
      bool testLengthadjust;
#define DEFAULT_testLengthadjust false

      /**
       * Also use Baldi to train the profiles, and calculate forward (etc) scores using it?
       */
      bool testBaldi;
#define DEFAULT_testBaldi false

      /**
       * Also use Baldi / Siegel to train the profiles, and calculate forward (etc) scores using it?
       */
      bool testBaldiSiegel;
#define DEFAULT_testBaldiSiegel false

      /**
       * If startWithUniformPositions or startWithPositionsDrawnFromPrior are
       * true, then under normal circumstances numStartingProfiles profiles
       * would be randomly generated; if alsoStartWithEvenPositions is *also*
       * true, then the first starting profile (index 0) will have even
       * positions rather than randomly-generated position emission values.
       */
      bool alsoStartWithEvenPositions;
#define DEFAULT_alsoStartWithEvenPositions false

      Parameters ();
      virtual ~Parameters () {};
    
      // Copy constructor
      template <class AnyParameters>
      Parameters ( const AnyParameters & copy_from );
    
      // Copy constructor/operator
      template <class AnyParameters>
      Parameters &
      operator= (
        const AnyParameters & copy_from
      );
    
      template <class AnyParameters>
      void
      copyFromNonVirtual (
        AnyParameters const & copy_from
      );

      template <class AnyParameters>
      void
      copyFromNonVirtualDontDelegate (
        AnyParameters const & copy_from
      );

      virtual
      void
      copyFrom ( const Parameters & copy_from );
    
      virtual
      void
      resetToDefaults ();

      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        Parameters const& parameters
      )
      {
        parameters.writeParameters( os );
        return os;
      } // friend operator<< ( basic_ostream &, Parameters const& )

      template<class CharT, class Traits>
      void
      writeParameters (
        std::basic_ostream<CharT,Traits>& os
      ) const;

      template <typename T> inline void COUT_IT ( T val ) { cout << val; }
      template <typename T, typename A> inline void COUT_IT ( vector<T,A> val ) { cout << "<vector>"; }

#undef PROFUSETEST_DEFAULT_TMP_ARRAY_TO_VECTOR
#undef GALOSH_DEF_OPT
#define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) inline void SET_##NAME( const TYPE new_value ) { /* cout << "SETTING " << #NAME; */ this->NAME = new_value; /* cout << " to "; COUT_IT( this->NAME ); cout << endl; */ }
       #include "ProfuseTestOptions.hpp"
#undef GALOSH_DEF_OPT

    }; // End inner class Parameters

    template <class ParametersType>
    class ParametersModifierTemplate :
      public ProfileGibbs<ProfileTreeRoot<ResidueType, ProbabilityType>,ScoreType,MatrixValueType,SequenceResidueType>::template ParametersModifierTemplate<ParametersType>
    {
      typedef typename ProfileGibbs<ProfileTreeRoot<ResidueType, ProbabilityType>,ScoreType,MatrixValueType,SequenceResidueType>::template ParametersModifierTemplate<ParametersType> base_parameters_modifier_t; 

      // Boost serialization
    private:
      friend class boost::serialization::access;
      template<class Archive>
      void serialize ( Archive & ar, const unsigned int /* file_version */ )
      {
        // save/load base class information.  This will serialize the
        // parameters too.
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( base_parameters_modifier_t );

        // Serialize the new isModified_ stuff
        ar & BOOST_SERIALIZATION_NVP( isModified_saveResultsToFile );
        ar & BOOST_SERIALIZATION_NVP( isModified_saveResultsParentDirectory );
        ar & BOOST_SERIALIZATION_NVP( isModified_resultsFilePrefix );
        ar & BOOST_SERIALIZATION_NVP( isModified_tabFileSuffix );
        ar & BOOST_SERIALIZATION_NVP( isModified_parametersFileSuffix );
        ar & BOOST_SERIALIZATION_NVP( isModified_saveTrueProfileTrees );
        ar & BOOST_SERIALIZATION_NVP( isModified_trueProfileTreeFileSuffix );
        ar & BOOST_SERIALIZATION_NVP( isModified_saveStartingProfiles );
        ar & BOOST_SERIALIZATION_NVP( isModified_startingProfileTreeFileSuffix );
        ar & BOOST_SERIALIZATION_NVP( isModified_saveTestProfiles );
        ar & BOOST_SERIALIZATION_NVP( isModified_testProfileTreeFileSuffix );
        ar & BOOST_SERIALIZATION_NVP( isModified_savePatternSequences );
        ar & BOOST_SERIALIZATION_NVP( isModified_patternSequencesFileSuffix );
        ar & BOOST_SERIALIZATION_NVP( isModified_saveTests );
        ar & BOOST_SERIALIZATION_NVP( isModified_patternSequencesFileSuffix );
        ar & BOOST_SERIALIZATION_NVP( isModified_saveTrainingSequences );
        ar & BOOST_SERIALIZATION_NVP( isModified_trainingSequencesFileSuffix );
        ar & BOOST_SERIALIZATION_NVP( isModified_saveTestingSequences );
        ar & BOOST_SERIALIZATION_NVP( isModified_testingSequencesFileSuffix );
        ar & BOOST_SERIALIZATION_NVP( isModified_saveTrueTrainingAlignments );
        ar & BOOST_SERIALIZATION_NVP( isModified_trainingTrueAlignmentsFileSuffix );
        ar & BOOST_SERIALIZATION_NVP( isModified_saveTrueTestingAlignments );
        ar & BOOST_SERIALIZATION_NVP( isModified_trueTestingAlignmentsFileSuffix );
        ar & BOOST_SERIALIZATION_NVP( isModified_saveFileVersion );
        ar & BOOST_SERIALIZATION_NVP( isModified_numProfiles );
        ar & BOOST_SERIALIZATION_NVP( isModified_profileLengths );
        ar & BOOST_SERIALIZATION_NVP( isModified_sharedPositionRate );
        ar & BOOST_SERIALIZATION_NVP( isModified_numTrainingSequencesPerProfiles );
        ar & BOOST_SERIALIZATION_NVP( isModified_numTestingSequencesPerProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_conservationRates );
        ar & BOOST_SERIALIZATION_NVP( isModified_useDeletionsForInsertionsParameters );
        ar & BOOST_SERIALIZATION_NVP( isModified_expectedDeletionsCounts );
        ar & BOOST_SERIALIZATION_NVP( isModified_expectedInsertionsCounts );
        ar & BOOST_SERIALIZATION_NVP( isModified_expectedDeletionLengthAsProfileLengthFractions );
        ar & BOOST_SERIALIZATION_NVP( isModified_expectedInsertionLengthAsProfileLengthFractions );
        ar & BOOST_SERIALIZATION_NVP( isModified_minExpectedDeletionLength );
        ar & BOOST_SERIALIZATION_NVP( isModified_minExpectedInsertionLength );
        ar & BOOST_SERIALIZATION_NVP( isModified_preAlignInsertion );
        ar & BOOST_SERIALIZATION_NVP( isModified_postAlignInsertion );
        ar & BOOST_SERIALIZATION_NVP( isModified_priorStrength );
        ar & BOOST_SERIALIZATION_NVP( isModified_priorStrength_internal_transitions );
        ar & BOOST_SERIALIZATION_NVP( isModified_priorMtoM );
        ar & BOOST_SERIALIZATION_NVP( isModified_priorMtoI );
        ar & BOOST_SERIALIZATION_NVP( isModified_priorMtoD );
        ar & BOOST_SERIALIZATION_NVP( isModified_priorItoM );
        ar & BOOST_SERIALIZATION_NVP( isModified_priorItoI );
        ar & BOOST_SERIALIZATION_NVP( isModified_priorDtoM );
        ar & BOOST_SERIALIZATION_NVP( isModified_priorDtoD );
        ar & BOOST_SERIALIZATION_NVP( isModified_reportGibbsMean );
        ar & BOOST_SERIALIZATION_NVP( isModified_reportGibbsMode );
        ar & BOOST_SERIALIZATION_NVP( isModified_numTrueProfiles );
        ar & BOOST_SERIALIZATION_NVP( isModified_numStartingProfiles );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformGlobals );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformGlobals_scalar );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformGlobals_maxNtoN );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformGlobals_maxBtoD );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformGlobals_maxMtoI );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformGlobals_maxMtoD );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformGlobals_maxItoI );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformGlobals_maxDtoD );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformGlobals_maxCtoC );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformPositions );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithGlobalsDrawnFromPrior );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithPositionsDrawnFromPrior );
        ar & BOOST_SERIALIZATION_NVP( isModified_testViterbi );
        ar & BOOST_SERIALIZATION_NVP( isModified_coutViterbi );
        ar & BOOST_SERIALIZATION_NVP( isModified_testTruepath );
        ar & BOOST_SERIALIZATION_NVP( isModified_coutTruepath );
        ar & BOOST_SERIALIZATION_NVP( isModified_calculateSymmeterizedKullbackLeiblerDistancesToTrue );
        ar & BOOST_SERIALIZATION_NVP( isModified_calculateSymmeterizedKullbackLeiblerDistancesToStarting );
        ar & BOOST_SERIALIZATION_NVP( isModified_coutDistances );
        ar & BOOST_SERIALIZATION_NVP( isModified_calculateProfileProfileAlignments );
        ar & BOOST_SERIALIZATION_NVP( isModified_profileProfileIndelOpenCost );
        ar & BOOST_SERIALIZATION_NVP( isModified_profileProfileIndelExtensionCost );
        ar & BOOST_SERIALIZATION_NVP( isModified_testTrueProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_coutTrueProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_testStartingProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_coutStartingProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_testUnconditionalProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_coutUnconditionalProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_testUnconditionalWithFixedStartingGlobalsProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_coutUnconditionalWithFixedStartingGlobalsProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_testUnconditionalWithFixedTrueGlobalsProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_coutUnconditionalWithFixedTrueGlobalsProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_testConditionalThenUnconditionalProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_coutConditionalThenUnconditionalProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_testUnconditionalThenConditionalProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_coutUnconditionalThenConditionalProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_testConditionalGibbsProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_coutConditionalGibbsProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_testUnconditionalGibbsProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_coutUnconditionalGibbsProfile );
        ar & BOOST_SERIALIZATION_NVP( isModified_testLengthadjust );
        ar & BOOST_SERIALIZATION_NVP( isModified_testBaldi );
        ar & BOOST_SERIALIZATION_NVP( isModified_testBaldiSiegel );
        ar & BOOST_SERIALIZATION_NVP( isModified_alsoStartWithEvenPositions );
      } // serialize( Archive &, const unsigned int )

    public:
  
      /// isModified flags for Parameters
      bool isModified_saveResultsToFile;

      bool isModified_saveResultsParentDirectory;

      bool isModified_resultsFilePrefix;

      bool isModified_tabFileSuffix;

      bool isModified_parametersFileSuffix;

      bool isModified_saveTrueProfileTrees;

      bool isModified_trueProfileTreeFileSuffix;

      bool isModified_saveStartingProfiles;

      bool isModified_startingProfileTreeFileSuffix;

      bool isModified_saveTestProfiles;

      bool isModified_testProfileTreeFileSuffix;

      bool isModified_savePatternSequences;

      bool isModified_patternSequencesFileSuffix;

      bool isModified_saveTests;

      bool isModified_testsFileSuffix;

      bool isModified_saveTrainingSequences;

      bool isModified_trainingSequencesFileSuffix;

      bool isModified_saveTestingSequences;

      bool isModified_testingSequencesFileSuffix;

      bool isModified_saveTrueTrainingAlignments;

      bool isModified_trainingTrueAlignmentsFileSuffix;

      bool isModified_saveTrueTestingAlignments;

      bool isModified_trueTestingAlignmentsFileSuffix;

      bool isModified_saveFileVersion;

      bool isModified_numProfiles;

      bool isModified_profileLengths;

      /**
       * When numProfiles > 1, what fraction of the positions of each child
       * should be shared with its parent?
       */
      bool isModified_sharedPositionRate;

      bool isModified_numTrainingSequencesPerProfiles;

      bool isModified_numTestingSequencesPerProfile;

      /**
       * When making the true root profile from the pattern sequence, use this
       * probability for the pattern sequence base at each position, and divide
       * the remaining probability evenly among the remaining bases.
       */
      bool isModified_conservationRates;

      /**
       * Lock the indel parameters of the true profile to be the same for
       * insertions as for deletions?  This makes the expectedInsertionsCounts,
       * expectedInsertionLengthAsProfileLengthFractions, and
       * minExpectedInsertionLength unused, since the corresponding deletion
       * values will be used instead.  It also reduces the number of tests by
       * reducing the number of possible combinations (since deletions and
       * insertions will go in lock step).
       */
      bool isModified_useDeletionsForInsertionsParameters;

      /**
       * The deletionOpen value of the true profile will be set to (
       * expectedDeletionsCounts / profileLength ).
       */
      bool isModified_expectedDeletionsCounts;

      /**
       * The insertionOpen value of the true profile will be set to (
       * expectedInsertionsCounts / profileLength ).
       */
      bool isModified_expectedInsertionsCounts;

      /**
       * The deletionExtension value of the true profile will be the minimum of
       * ( 1.0 / ( expectedDeletionLengthAsProfileLengthFractions *
       * profileLength ) ) and ( 1.0 / minExpectedDeletionLength ).
       */
      bool isModified_expectedDeletionLengthAsProfileLengthFractions;

      /**
       * The insertionExtension value of the true profile will be the minimum of
       * ( 1.0 / ( expectedInsertionLengthAsProfileLengthFractions *
       * profileLength ) ) and ( 1.0 / minExpectedInsertionLength ).
       */
      bool isModified_expectedInsertionLengthAsProfileLengthFractions;

      /**
       * The deletionExtension value of the true profile will be the minimum of
       * ( 1.0 / ( expectedDeletionLengthAsProfileLengthFractions *
       * profileLength ) ) and ( 1.0 / minExpectedDeletionLength ).
       */
      bool isModified_minExpectedDeletionLength;

      /**
       * The insertionExtension value of the true profile will be the minimum of
       * ( 1.0 / ( expectedInsertionLengthAsProfileLengthFractions *
       * profileLength ) ) and ( 1.0 / minExpectedInsertionLength ).
       */
      bool isModified_minExpectedInsertionLength;

      /**
       * The preAlignInsertion value of the true profile.
       */
      bool isModified_preAlignInsertion;

      /**
       * The postAlignInsertion value of the true profile.
       */
      bool isModified_postAlignInsertion;

      /**
       * The effective number of sequences "observed" a priori.
       */
      bool isModified_priorStrength;

      /**
       * The effective number of sequences "observed" a priori.
       */
      bool isModified_priorStrength_internal_transitions;

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * M->M transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      bool isModified_priorMtoM;

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * M->I transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      bool isModified_priorMtoI;

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * M->D transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      bool isModified_priorMtoD;

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * I->M transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      bool isModified_priorItoM;

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * I->I transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      bool isModified_priorItoI;

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * D->M transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      bool isModified_priorDtoM;

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * D->D transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      bool isModified_priorDtoD;

      /**
       * Additionally report the overall mean of all chains found while
       * performing Gibbs sampling?  The best profile is always reported, which
       * may be the overall mean, the mode, or the mean of one of the chains.
       */
      bool isModified_reportGibbsMean;

      /**
       * Additionally report the mode found while
       * performing Gibbs sampling?  The best profile is always reported, which
       * may be the overall mean, the mode, or the mean of one of the chains.
       *
       * Note that it takes some extra time to store the mode (this turns on
       * saveGibbsMode in the ProfileGibbs class).
       */
      bool isModified_reportGibbsMode;

      /**
       * We do the whole thing a number of different times, starting over with
       * a new pattern sequence.
       */
      bool isModified_numTrueProfiles;

      /**
       * For each true root profile, we run the trainers from a number of
       * different starting profiles.  This is that number.
       */
      bool isModified_numStartingProfiles;

      bool isModified_startWithUniformGlobals;

      bool isModified_startWithUniformGlobals_scalar;

      bool isModified_startWithUniformGlobals_maxNtoN;

      bool isModified_startWithUniformGlobals_maxBtoD;

      bool isModified_startWithUniformGlobals_maxMtoI;

      bool isModified_startWithUniformGlobals_maxMtoD;

      bool isModified_startWithUniformGlobals_maxItoI;

      bool isModified_startWithUniformGlobals_maxDtoD;

      bool isModified_startWithUniformGlobals_maxCtoC;

      bool isModified_startWithUniformPositions;

      bool isModified_startWithGlobalsDrawnFromPrior;

      bool isModified_startWithPositionsDrawnFromPrior;

      /**
       * Calculate the viterbi scores after training each profile?  Note that
       * this is in addition to the forward scores, which we always calculate.
       */
      bool isModified_testViterbi;

      /**
       * Write the viterbi scores to STDOUT?
       */
      bool isModified_coutViterbi;

      /**
       * Calculate the truepath scores after training each profile?  Note that
       * this is in addition to the forward scores, which we always calculate.
       */
      bool isModified_testTruepath;

      /**
       * Write the truepath scores to STDOUT?
       */
      bool isModified_coutTruepath;

      /**
       * Calculate the SKL distance between the training profiles and the
       * true profile?
       */
      bool isModified_calculateSymmeterizedKullbackLeiblerDistancesToTrue;

      /**
       * Calculate the SKL distance between the training profiles and the
       * starting profile?
       */
      bool isModified_calculateSymmeterizedKullbackLeiblerDistancesToStarting;

      /**
       * Write the SKL distances to STDOUT?
       */
      bool isModified_coutDistances;

      bool isModified_calculateProfileProfileAlignments;

      bool isModified_profileProfileIndelOpenCost;

      bool isModified_profileProfileIndelExtensionCost;

      /**
       * Calculate forward (etc) scores for the true profile, too?
       */
      bool isModified_testTrueProfile;

      /**
       * Write forward (etc) scores for the true profile to STDOUT
       * during training?
       */
      bool isModified_coutTrueProfile;

      /**
       * Calculate forward (etc) scores for the pre-training profile, too?
       */
      bool isModified_testStartingProfile;

      /**
       * Write forward (etc) scores for the pre-training profile to STDOUT
       * during training?
       */
      bool isModified_coutStartingProfile;

      /**
       * Also train the unconditional profile, and calculate forward (etc)
       * scores using it?
       */
      bool isModified_testUnconditionalProfile;

      /**
       * Write the forward (etc) scores for the unconditional profile to STDOUT
       * during training?
       */
      bool isModified_coutUnconditionalProfile;

      /**
       * Also train the unconditional (with fixed starting globals) profile,
       * and calculate forward (etc) scores using it?
       */
      bool isModified_testUnconditionalWithFixedStartingGlobalsProfile;

      /**
       * Write the forward (etc) scores for the unconditional (with fixed
       * starting globals) profile to STDOUT during training?
       */
      bool isModified_coutUnconditionalWithFixedStartingGlobalsProfile;

      /**
       * Also train the unconditional (with fixed true globals) profile,
       * and calculate forward (etc) scores using it?
       */
      bool isModified_testUnconditionalWithFixedTrueGlobalsProfile;

      /**
       * Write the forward (etc) scores for the unconditional (with fixed
       * true globals) profile to STDOUT during training?
       */
      bool isModified_coutUnconditionalWithFixedTrueGlobalsProfile;

      /**
       * Also train the "conditional, then unconditional" profile, and
       * calculate forward (etc) scores using it?
       */
      bool isModified_testConditionalThenUnconditionalProfile;

      /**
       * Write the forward (etc) scores for the "conditional, then
       * unconditional" profile to STDOUT during training?
       */
      bool isModified_coutConditionalThenUnconditionalProfile;

      /**
       * Also train the "unconditional, then conditional" profile, and
       * calculate forward (etc) scores using it?
       */
      bool isModified_testUnconditionalThenConditionalProfile;

      /**
       * Write the forward (etc) scores for the "unconditional, then
       * conditional" profile to STDOUT during training?
       */
      bool isModified_coutUnconditionalThenConditionalProfile;

      /**
       * Also train the "unconditional (first with fixed starting globals, then
       * with fixed positions)" profile, and calculate forward (etc) scores
       * using it?
       */
      bool isModified_testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile;

      /**
       * Write the forward (etc) scores for the unconditional (first with fixed
       * starting globals, then with fixed positions) profile to STDOUT during
       * training?
       */
      bool isModified_coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile;

      /**
       * Also train the "unconditional (first with fixed true globals, then
       * with fixed positions)" profile, and calculate forward (etc) scores
       * using it?
       */
      bool isModified_testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile;

      /**
       * Write the forward (etc) scores for the unconditional (first with fixed
       * true globals, then with fixed positions) profile to STDOUT during
       * training?
       */
      bool isModified_coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile;

      bool isModified_testConditionalGibbsProfile;
      bool isModified_coutConditionalGibbsProfile;
      bool isModified_testUnconditionalGibbsProfile;
      bool isModified_coutUnconditionalGibbsProfile;

      bool isModified_testLengthadjust;

      bool isModified_testBaldi;

      bool isModified_testBaldiSiegel;

      bool isModified_alsoStartWithEvenPositions;

      ParametersModifierTemplate ();
    
      // Copy constructor
      template <class AnyParametersModifierTemplate>
      ParametersModifierTemplate ( const AnyParametersModifierTemplate & copy_from );
    
      // Copy constructor/operator
      template <class AnyParametersModifierTemplate>
      ParametersModifierTemplate & operator= (
        const AnyParametersModifierTemplate & copy_from
      );
    
      template <class AnyParametersModifierTemplate>
      void
      copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      );

      template <class AnyParametersModifierTemplate>
      void
      isModified_copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      );

      void
      reset ();

      void
      isModified_reset ();

      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        ParametersModifierTemplate const& parameters_modifier
      )
      {
        parameters_modifier.writeParametersModifier( os );

        return os;
      } // friend operator<< ( basic_ostream &, ParametersModifierTemplate const& )

      template<class CharT, class Traits>
      void
      writeParametersModifier (
        std::basic_ostream<CharT,Traits>& os
      ) const;

      template <class AnyParameters>
      void
      applyModifications ( AnyParameters & target_parameters );

    }; // End inner class ParametersModifierTemplate

    typedef ParametersModifierTemplate<typename ProfuseTest::Parameters> ParametersModifier;

    typedef ProfileTreeRoot<ResidueType, ProbabilityType> RootType;
    typedef ProfileTreeRoot<ResidueType, ProbabilityType> InternalNodeType;
    typedef ProfileTree<ResidueType, ProbabilityType, InternalNodeType > ProfileTreeType;
    typedef typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template MultipleAlignment<RootType, SequenceResidueType> RootTypeMultipleAlignment;
    typedef typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template MultipleAlignment<InternalNodeType, SequenceResidueType> InternalNodeTypeMultipleAlignment;
    typedef typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::RowVector RowVector;
    typedef typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer SequentialAccessContainer;

    class Test
    {
      // Boost serialization
    private:
      friend class boost::serialization::access;
      template<class Archive>
      void serialize ( Archive & ar, const unsigned int /* file_version */ )
      {
        ar & BOOST_SERIALIZATION_NVP( name );
        ar & BOOST_SERIALIZATION_NVP( isRun );
        ar & BOOST_SERIALIZATION_NVP( isCout );
        ar & BOOST_SERIALIZATION_NVP( isGibbs );
        ar & BOOST_SERIALIZATION_NVP( coutLeftBrace );
        ar & BOOST_SERIALIZATION_NVP( coutRightBrace );
        ar & BOOST_SERIALIZATION_NVP( parametersModifier );
        ar & BOOST_SERIALIZATION_NVP( startingGlobalsTest );
        ar & BOOST_SERIALIZATION_NVP( startingPositionsTest );
      } // serialize( Archive &, const unsigned int )

    public:

      /**
       * The name of this test.  This will be the base of the names used in the
       * output file and in the STDOUT display.
       */
      string name;

      /**
       * Should this test be run?
       */
      bool isRun;

      /**
       * Should this test be displayed to STDOUT while the tests are running?
       */
      bool isCout;

      /**
       * Does this test use ProfileGibbs instead of ProfileTrainer?
       */
      bool isGibbs;

      /**
       * When displaying the score to cout, precede the score(s) of this test
       * by this string.
       *
       * @see coutRightBrace
       */
      string coutLeftBrace;

      /**
       * When displaying the score to cout, succeed the score(s) of this test
       * by this string.
       *
       * @see coutLeftBrace
       */
      string coutRightBrace;

      /**
       * The parameters modifier for this test.
       */
      typename ProfileGibbs<RootType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifier parametersModifier;

      /**
       * A pointer to the test whose (trained) profile is to be used for the
       * starting values of the global parameters for this test.  This pointer
       * can be NULL to indicate that the usual "starting profile" should be
       * used, or it can point to this test to indicate that the "true profile"
       * should be used.
       */
      Test * startingGlobalsTest;

      /**
       * A pointer to the test whose (trained) profile is to be used for the
       * starting values of the position-specific parameters for this test.
       * This pointer can be NULL to indicate that the usual "starting profile"
       * should be used, or it can point to this test to indicate that the
       * "true profile" should be used.
       */
      Test * startingPositionsTest;

      Test () :
        name( "undefined" ),
        isRun( false ),
        isCout( false ),
        isGibbs( false ),
        coutLeftBrace( "" ),
        coutRightBrace( "" ),
        startingGlobalsTest( NULL ),
        startingPositionsTest( NULL )
      {
      } //<init>()

    }; // End inner class Test

    template <class Serializable>
    static void
    writeXML (
      Serializable const & profuse_object,
      const char * filename
    )
    {
      // TODO: PUT BACK.  Why is this causing a problem?  Multiple threads?
      // make an archive
      std::ofstream ofs( filename );
      assert( ofs.good() );
      boost::archive::xml_oarchive oa( ofs );
      oa << BOOST_SERIALIZATION_NVP( profuse_object );
      
    } // writeXML( Serializable const & )
    
    template <class Serializable>
    static void
    readXML (
      Serializable & profuse_object,
      const char * filename
    )
    {
      // open the archive
      std::ifstream ifs( filename );
      assert( ifs.good() );
      boost::archive::xml_iarchive ia( ifs );
    
      // restore the profile from the archive
      ia >> BOOST_SERIALIZATION_NVP( profuse_object );
    } // readXML( Serializable &, const char * )

    /**
     * The parameters currently being used.
     */
    Parameters m_parameters;

  // If true, the .tab output file will contain values converted using toLogDouble(..).
  // TODO: Make this into a parameter
    static const bool convert_tab_output_to_log_double = false;

    /**
     * For convenience we number the tests from 0 to LAST_TEST_ID.
     */
     static const uint32_t LAST_TEST_ID = 21;

    /**
     * For convenience we number the tests from 0 to LAST_TEST_ID.
     *
     * This is an enum of the test ids, from 0 to LAST_TEST_ID.
     */
    static const enum
    {
      TEST_ID_true = 0,
      TEST_ID_starting,
      TEST_ID_conditional,
      TEST_ID_unconditional,
      TEST_ID_unconditional_with_fixed_starting_globals,
      TEST_ID_unconditional_with_fixed_true_globals,
      TEST_ID_conditional_then_unconditional,
      TEST_ID_unconditional_then_conditional,
      TEST_ID_unconditional_with_fixed_starting_globals_then_with_fixed_positions,
      TEST_ID_unconditional_with_fixed_true_globals_then_with_fixed_positions,
      TEST_ID_gibbs_conditional,
      TEST_ID_gibbs_unconditional,
      TEST_ID_lengthadjust_conditional,
      TEST_ID_lengthadjust_unconditional,
      TEST_ID_baldi_conditional,
      TEST_ID_baldi_unconditional,
      TEST_ID_baldi_lengthadjust_conditional,
      TEST_ID_baldi_lengthadjust_unconditional,
      TEST_ID_baldi_siegel_conditional,
      TEST_ID_baldi_siegel_unconditional,
      TEST_ID_baldi_siegel_lengthadjust_conditional,
      TEST_ID_baldi_siegel_lengthadjust_unconditional
    } TEST_ID;

    /**
     * The Random object we're using.
     */
    Random m_random;

    /**
     * Construct a profuse test object, using the
     * time to get a random seed.
     */  
    ProfuseTest ();

    ///TAH 6/12 constructor from commandline options
   ProfuseTest ( int argc, char **argv );

    /**
     * Construct a profuse test object, using the provided seed.
     */  
    ProfuseTest (
      uint32_t const seed
    );

    void
    start ();

  }; // End class ProfuseTest

  //======//// potentially non-inline implementations ////========//


  ////// Class galosh::ProfuseTest::Parameters ////
  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_INIT
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      Parameters ()
      {
        if( DEFAULT_debug >= DEBUG_All ) {
          cout << "[debug] ProfuseTest::Parameters::<init>()" << endl;
        } // End if DEBUG_All
        resetToDefaults();
#undef GALOSH_DEF_OPT
#define PROFUSETEST_DEFAULT_TMP_ARRAY_TO_VECTOR
#define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP)          \
        DEFAULT_OPTIONS_DESCRIPTION.add_options()(#NAME,po::value<TYPE>()->default_value(DEFAULTVAL)->notifier( boost::bind(&( Parameters::SET_##NAME ), this, _1) ) TMP_EXTRA_STUFF,HELP)
        //DEFAULT_OPTIONS_DESCRIPTION.add_options()(#NAME,po::value<TYPE>()->default_value(DEFAULTVAL) TMP_EXTRA_STUFF,HELP)
        #include "ProfuseTestOptions.hpp"  /// define all the commandline options for this module
#undef GALOSH_DEF_OPT
      } // <init>()

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class AnyParameters>
  GALOSH_INLINE_INIT
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      // Copy constructor
      Parameters ( const AnyParameters & copy_from )
      {
        //if( static_cast<galosh::Parameters>( copy_from ).debug >= DEBUG_All ) {
        //  cout << "[debug] ProfuseTest::Parameters::<init>( copy_from )" << endl;
        //} // End if DEBUG_All
        copyFromNonVirtual( copy_from );
      } // <init>( AnyParameters const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class AnyParameters>
  GALOSH_INLINE_TRIVIAL
  typename ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters &
      // Copy constructor/operator
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      operator= (
        const AnyParameters & copy_from
      )
      {
        if( copy_from.debug >= DEBUG_All ) {
          cout << "[debug] ProfuseTest::Parameters::operator=( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
        return *this;
      } // operator=( AnyParameters const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class AnyParameters>
  GALOSH_INLINE_TRIVIAL
  void
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      copyFromNonVirtual (
        AnyParameters const & copy_from
      )
      {
        ProfileGibbs<ProfileTreeRoot<ResidueType, ProbabilityType>,ScoreType,MatrixValueType,SequenceResidueType>::Parameters::copyFromNonVirtual( copy_from );
        //if( copy_from.debug >= DEBUG_All ) {
        //  cout << "[debug] ProfuseTest::Parameters::copyFromNonVirtual( copy_from )" << endl;
        //} // End if DEBUG_All
        copyFromNonVirtualDontDelegate( copy_from );
      } // copyFromNonVirtual( AnyParameters const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class AnyParameters>
  GALOSH_INLINE_COPY
  void
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      copyFromNonVirtualDontDelegate (
        AnyParameters const & copy_from
      )
      {
        saveResultsToFile =                            copy_from.saveResultsToFile;
        saveResultsParentDirectory =                           copy_from.saveResultsParentDirectory;
        resultsFilePrefix =                            copy_from.resultsFilePrefix;
        tabFileSuffix =                            copy_from.tabFileSuffix;
        parametersFileSuffix =                            copy_from.parametersFileSuffix;
        saveTrueProfileTrees =                            copy_from.saveTrueProfileTrees;
        trueProfileTreeFileSuffix =                            copy_from.trueProfileTreeFileSuffix;
        saveStartingProfiles =                            copy_from.saveStartingProfiles;
        startingProfileTreeFileSuffix =                            copy_from.startingProfileTreeFileSuffix;
        saveTestProfiles =                            copy_from.saveTestProfiles;
        testProfileTreeFileSuffix =                            copy_from.testProfileTreeFileSuffix;
        savePatternSequences =                            copy_from.savePatternSequences;
        patternSequencesFileSuffix =                            copy_from.patternSequencesFileSuffix;
        saveTests =                            copy_from.saveTests;
        testsFileSuffix =                            copy_from.testsFileSuffix;
        saveTrainingSequences =                            copy_from.saveTrainingSequences;
        trainingSequencesFileSuffix =                            copy_from.trainingSequencesFileSuffix;
        saveTestingSequences =                            copy_from.saveTestingSequences;
        testingSequencesFileSuffix =                            copy_from.testingSequencesFileSuffix;
        saveTrueTrainingAlignments =                            copy_from.saveTrueTrainingAlignments;
        trainingTrueAlignmentsFileSuffix =                            copy_from.trainingTrueAlignmentsFileSuffix;
        saveTrueTestingAlignments =                            copy_from.saveTrueTestingAlignments;
        trueTestingAlignmentsFileSuffix =                            copy_from.trueTestingAlignmentsFileSuffix;
        saveFileVersion =                              copy_from.saveFileVersion;
        numProfiles =                                  copy_from.numProfiles;
        profileLengths =                                copy_from.profileLengths;
        sharedPositionRate =                           copy_from.sharedPositionRate;
        numTrainingSequencesPerProfiles =               copy_from.numTrainingSequencesPerProfiles;
        numTestingSequencesPerProfile =                   copy_from.numTestingSequencesPerProfile;
        conservationRates =                             copy_from.conservationRates;
        useDeletionsForInsertionsParameters =                             copy_from.useDeletionsForInsertionsParameters;
        expectedDeletionsCounts =                          copy_from.expectedDeletionsCounts;
        expectedInsertionsCounts =                          copy_from.expectedInsertionsCounts;
        expectedDeletionLengthAsProfileLengthFractions =                          copy_from.expectedDeletionLengthAsProfileLengthFractions;
        expectedInsertionLengthAsProfileLengthFractions =                          copy_from.expectedInsertionLengthAsProfileLengthFractions;
        minExpectedDeletionLength =                          copy_from.minExpectedDeletionLength;
        minExpectedInsertionLength =                          copy_from.minExpectedInsertionLength;
        preAlignInsertion =                          copy_from.preAlignInsertion;
        postAlignInsertion =                          copy_from.postAlignInsertion;
        priorStrength =                          copy_from.priorStrength;
        priorStrength_internal_transitions =                          copy_from.priorStrength_internal_transitions;
        priorMtoM =                          copy_from.priorMtoM;
        priorMtoI =                          copy_from.priorMtoI;
        priorMtoD =                          copy_from.priorMtoD;
        priorItoM =                          copy_from.priorItoM;
        priorItoI =                          copy_from.priorItoI;
        priorDtoM =                          copy_from.priorDtoM;
        priorDtoD =                          copy_from.priorDtoD;
        reportGibbsMean =                          copy_from.reportGibbsMean;
        reportGibbsMode =                          copy_from.reportGibbsMode;
        numTrueProfiles =                          copy_from.numTrueProfiles;
        numStartingProfiles =                          copy_from.numStartingProfiles;
        startWithUniformGlobals =                          copy_from.startWithUniformGlobals;
        startWithUniformGlobals_scalar =                          copy_from.startWithUniformGlobals_scalar;
        startWithUniformGlobals_maxNtoN =                          copy_from.startWithUniformGlobals_maxNtoN;
        startWithUniformGlobals_maxBtoD =                          copy_from.startWithUniformGlobals_maxBtoD;
        startWithUniformGlobals_maxMtoI =                          copy_from.startWithUniformGlobals_maxMtoI;
        startWithUniformGlobals_maxMtoD =                          copy_from.startWithUniformGlobals_maxMtoD;
        startWithUniformGlobals_maxItoI =                          copy_from.startWithUniformGlobals_maxItoI;
        startWithUniformGlobals_maxDtoD =                          copy_from.startWithUniformGlobals_maxDtoD;
        startWithUniformGlobals_maxCtoC =                          copy_from.startWithUniformGlobals_maxCtoC;
        startWithUniformPositions =                          copy_from.startWithUniformPositions;
        startWithGlobalsDrawnFromPrior =                          copy_from.startWithGlobalsDrawnFromPrior;
        startWithPositionsDrawnFromPrior =                          copy_from.startWithPositionsDrawnFromPrior;
        testViterbi =                                   copy_from.testViterbi;
        coutViterbi =                                   copy_from.coutViterbi;
        testTruepath =                                   copy_from.testTruepath;
        coutTruepath =                                   copy_from.coutTruepath;
        calculateSymmeterizedKullbackLeiblerDistancesToTrue =        copy_from.calculateSymmeterizedKullbackLeiblerDistancesToTrue;
        calculateSymmeterizedKullbackLeiblerDistancesToStarting =        copy_from.calculateSymmeterizedKullbackLeiblerDistancesToStarting;
        coutDistances =     copy_from.coutDistances;
        calculateProfileProfileAlignments =   copy_from.calculateProfileProfileAlignments;
        profileProfileIndelOpenCost =   copy_from.profileProfileIndelOpenCost;
        profileProfileIndelExtensionCost =   copy_from.profileProfileIndelExtensionCost;
        testTrueProfile =                          copy_from.testTrueProfile;
        coutTrueProfile =                          copy_from.coutTrueProfile;
        testStartingProfile =                          copy_from.testStartingProfile;
        coutStartingProfile =                          copy_from.coutStartingProfile;
        testUnconditionalProfile =                     copy_from.testUnconditionalProfile;
        coutUnconditionalProfile =                     copy_from.coutUnconditionalProfile;
        testUnconditionalWithFixedStartingGlobalsProfile =                     copy_from.testUnconditionalWithFixedStartingGlobalsProfile;
        coutUnconditionalWithFixedStartingGlobalsProfile =                     copy_from.coutUnconditionalWithFixedStartingGlobalsProfile;
        testUnconditionalWithFixedTrueGlobalsProfile =                     copy_from.testUnconditionalWithFixedTrueGlobalsProfile;
        coutUnconditionalWithFixedTrueGlobalsProfile =                     copy_from.coutUnconditionalWithFixedTrueGlobalsProfile;
        testConditionalThenUnconditionalProfile =      copy_from.testConditionalThenUnconditionalProfile;
        coutConditionalThenUnconditionalProfile =      copy_from.coutConditionalThenUnconditionalProfile;
        testUnconditionalThenConditionalProfile =      copy_from.testUnconditionalThenConditionalProfile;
        coutUnconditionalThenConditionalProfile =      copy_from.coutUnconditionalThenConditionalProfile;
        testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile =  copy_from.testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile;
        coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile =  copy_from.coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile;
        testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile =  copy_from.testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile;
        coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile =  copy_from.coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile;
        testConditionalGibbsProfile =                     copy_from.testConditionalGibbsProfile;
        coutConditionalGibbsProfile =                     copy_from.coutConditionalGibbsProfile;
        testUnconditionalGibbsProfile =                     copy_from.testUnconditionalGibbsProfile;
        coutUnconditionalGibbsProfile =                     copy_from.coutUnconditionalGibbsProfile;
        testLengthadjust = copy_from.testLengthadjust;
        testBaldi = copy_from.testBaldi;
        testBaldiSiegel = copy_from.testBaldiSiegel;
        alsoStartWithEvenPositions = copy_from.alsoStartWithEvenPositions;
      } // copyFromNonVirtualDontDelegate( AnyParameters const & )


  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_TRIVIAL
      void
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      copyFrom ( const Parameters & copy_from )
      {
        copyFromNonVirtual( copy_from );
      } // copyFrom( Parameters const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_REINITIALIZE
  void
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      resetToDefaults ()
      {
        ProfileGibbs<ProfileTreeRoot<ResidueType, ProbabilityType>,ScoreType,MatrixValueType,SequenceResidueType>::Parameters::resetToDefaults();
        // TODO: Why isn't the compiler finding "debug" in galosh::Parameters?
        //if( debug >= DEBUG_All ) {
        //  cout << "[debug] ProfuseTest::Parameters::resetToDefaults()" << endl;
        //} // End if DEBUG_All

        saveResultsToFile =                            DEFAULT_saveResultsToFile;
        saveResultsParentDirectory =                           DEFAULT_saveResultsParentDirectory;
        resultsFilePrefix =                            DEFAULT_resultsFilePrefix;
        tabFileSuffix =                            DEFAULT_tabFileSuffix;
        parametersFileSuffix =                            DEFAULT_parametersFileSuffix;
        saveTrueProfileTrees =                            DEFAULT_saveTrueProfileTrees;
        trueProfileTreeFileSuffix =                            DEFAULT_trueProfileTreeFileSuffix;
        saveStartingProfiles =                            DEFAULT_saveStartingProfiles;
        startingProfileTreeFileSuffix =                            DEFAULT_startingProfileTreeFileSuffix;
        saveTestProfiles =                            DEFAULT_saveTestProfiles;
        testProfileTreeFileSuffix =                            DEFAULT_testProfileTreeFileSuffix;
        savePatternSequences =                            DEFAULT_savePatternSequences;
        patternSequencesFileSuffix =                            DEFAULT_patternSequencesFileSuffix;
        saveTests =                            DEFAULT_saveTests;
        testsFileSuffix =                            DEFAULT_testsFileSuffix;
        saveTrainingSequences =                            DEFAULT_saveTrainingSequences;
        trainingSequencesFileSuffix =                            DEFAULT_trainingSequencesFileSuffix;
        saveTestingSequences =                            DEFAULT_saveTestingSequences;
        testingSequencesFileSuffix =                            DEFAULT_testingSequencesFileSuffix;
        saveTrueTrainingAlignments =                            DEFAULT_saveTrueTrainingAlignments;
        trainingTrueAlignmentsFileSuffix =                            DEFAULT_trainingTrueAlignmentsFileSuffix;
        saveTrueTestingAlignments =                            DEFAULT_saveTrueTestingAlignments;
        trueTestingAlignmentsFileSuffix =                            DEFAULT_trueTestingAlignmentsFileSuffix;
        saveFileVersion =                              DEFAULT_saveFileVersion;
        numProfiles =                                  DEFAULT_numProfiles;
        profileLengths =                                DEFAULT_profileLengths;
        sharedPositionRate =                           DEFAULT_sharedPositionRate;
        numTrainingSequencesPerProfiles =               DEFAULT_numTrainingSequencesPerProfiles;
        numTestingSequencesPerProfile =                   DEFAULT_numTestingSequencesPerProfile;
        conservationRates =                             DEFAULT_conservationRates;
        useDeletionsForInsertionsParameters =                             DEFAULT_useDeletionsForInsertionsParameters;
        expectedDeletionsCounts =                          DEFAULT_expectedDeletionsCounts;
        expectedInsertionsCounts =                          DEFAULT_expectedInsertionsCounts;
        expectedDeletionLengthAsProfileLengthFractions =                          DEFAULT_expectedDeletionLengthAsProfileLengthFractions;
        expectedInsertionLengthAsProfileLengthFractions =                          DEFAULT_expectedInsertionLengthAsProfileLengthFractions;
        minExpectedDeletionLength =                          DEFAULT_minExpectedDeletionLength;
        minExpectedInsertionLength =                          DEFAULT_minExpectedInsertionLength;
        preAlignInsertion =                          DEFAULT_preAlignInsertion;
        postAlignInsertion =                          DEFAULT_postAlignInsertion;
        priorStrength =                          DEFAULT_priorStrength;
        priorStrength_internal_transitions =                          DEFAULT_priorStrength_internal_transitions;
        priorMtoM =                          DEFAULT_priorMtoM;
        priorMtoI =                          DEFAULT_priorMtoI;
        priorMtoD =                          DEFAULT_priorMtoD;
        priorItoM =                          DEFAULT_priorItoM;
        priorItoI =                          DEFAULT_priorItoI;
        priorDtoM =                          DEFAULT_priorDtoM;
        priorDtoD =                          DEFAULT_priorDtoD;
        reportGibbsMean =                          DEFAULT_reportGibbsMean;
        reportGibbsMode =                          DEFAULT_reportGibbsMode;
        numTrueProfiles =                          DEFAULT_numTrueProfiles;
        numStartingProfiles =                          DEFAULT_numStartingProfiles;
        startWithUniformGlobals =                          DEFAULT_startWithUniformGlobals;
        startWithUniformGlobals_scalar =                          DEFAULT_startWithUniformGlobals_scalar;
        startWithUniformGlobals_maxNtoN =                          DEFAULT_startWithUniformGlobals_maxNtoN;
        startWithUniformGlobals_maxBtoD =                          DEFAULT_startWithUniformGlobals_maxBtoD;
        startWithUniformGlobals_maxMtoI =                          DEFAULT_startWithUniformGlobals_maxMtoI;
        startWithUniformGlobals_maxMtoD =                          DEFAULT_startWithUniformGlobals_maxMtoD;
        startWithUniformGlobals_maxItoI =                          DEFAULT_startWithUniformGlobals_maxItoI;
        startWithUniformGlobals_maxDtoD =                          DEFAULT_startWithUniformGlobals_maxDtoD;
        startWithUniformGlobals_maxCtoC =                          DEFAULT_startWithUniformGlobals_maxCtoC;
        startWithUniformPositions =                          DEFAULT_startWithUniformPositions;
        startWithGlobalsDrawnFromPrior =                          DEFAULT_startWithGlobalsDrawnFromPrior;
        startWithPositionsDrawnFromPrior =                          DEFAULT_startWithPositionsDrawnFromPrior;
        testViterbi =                                   DEFAULT_testViterbi;
        coutViterbi =                                   DEFAULT_coutViterbi;
        testTruepath =                                   DEFAULT_testTruepath;
        coutTruepath =                                   DEFAULT_coutTruepath;
        calculateSymmeterizedKullbackLeiblerDistancesToTrue =        DEFAULT_calculateSymmeterizedKullbackLeiblerDistancesToTrue;
        calculateSymmeterizedKullbackLeiblerDistancesToStarting =        DEFAULT_calculateSymmeterizedKullbackLeiblerDistancesToStarting;
        coutDistances =     DEFAULT_coutDistances;
        calculateProfileProfileAlignments =   DEFAULT_calculateProfileProfileAlignments;
        profileProfileIndelOpenCost =   DEFAULT_profileProfileIndelOpenCost;
        profileProfileIndelExtensionCost =   DEFAULT_profileProfileIndelExtensionCost;
        testTrueProfile =                          DEFAULT_testTrueProfile;
        coutTrueProfile =                          DEFAULT_coutTrueProfile;
        testStartingProfile =                          DEFAULT_testStartingProfile;
        coutStartingProfile =                          DEFAULT_coutStartingProfile;
        testUnconditionalProfile =                     DEFAULT_testUnconditionalProfile;
        coutUnconditionalProfile =                     DEFAULT_coutUnconditionalProfile;
        testUnconditionalWithFixedStartingGlobalsProfile =                     DEFAULT_testUnconditionalWithFixedStartingGlobalsProfile;
        coutUnconditionalWithFixedStartingGlobalsProfile =                     DEFAULT_coutUnconditionalWithFixedStartingGlobalsProfile;
        testUnconditionalWithFixedTrueGlobalsProfile =                     DEFAULT_testUnconditionalWithFixedTrueGlobalsProfile;
        coutUnconditionalWithFixedTrueGlobalsProfile =                     DEFAULT_coutUnconditionalWithFixedTrueGlobalsProfile;
        testConditionalThenUnconditionalProfile =      DEFAULT_testConditionalThenUnconditionalProfile;
        coutConditionalThenUnconditionalProfile =      DEFAULT_coutConditionalThenUnconditionalProfile;
        testUnconditionalThenConditionalProfile =      DEFAULT_testUnconditionalThenConditionalProfile;
        coutUnconditionalThenConditionalProfile =      DEFAULT_coutUnconditionalThenConditionalProfile;
        testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile =  DEFAULT_testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile;
        coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile =  DEFAULT_coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile;
        testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile =  DEFAULT_testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile;
        coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile =  DEFAULT_coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile;
        testConditionalGibbsProfile =                     DEFAULT_testConditionalGibbsProfile;
        coutConditionalGibbsProfile =                     DEFAULT_coutConditionalGibbsProfile;
        testUnconditionalGibbsProfile =                     DEFAULT_testUnconditionalGibbsProfile;
        coutUnconditionalGibbsProfile =                     DEFAULT_coutUnconditionalGibbsProfile;
        testLengthadjust = DEFAULT_testLengthadjust;
        testBaldi = DEFAULT_testBaldi;
        testBaldiSiegel = DEFAULT_testBaldiSiegel;
        alsoStartWithEvenPositions = DEFAULT_alsoStartWithEvenPositions;
      } // resetToDefaults()

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
  void
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      writeParameters (
        std::basic_ostream<CharT,Traits>& os
      ) const
      {
        ///  \todo fixyfix ITERATE through m_...whatever

        ProfileGibbs<ProfileTreeRoot<ResidueType, ProbabilityType>,ScoreType,MatrixValueType,SequenceResidueType>::Parameters::writeParameters( os );
        os << endl;

        os << "[ProfuseTest]" << endl;

        os << "saveResultsToFile = " <<                           saveResultsToFile << endl;
        os << "saveResultsParentDirectory = " <<                  saveResultsParentDirectory << endl;
        os << "resultsFilePrefix = " <<                          resultsFilePrefix << endl;
        os << "tabFileSuffix = " <<                          tabFileSuffix << endl;
        os << "parametersFileSuffix = " <<                          parametersFileSuffix << endl;
        os << "saveTrueProfileTrees = " <<                          saveTrueProfileTrees << endl;
        os << "trueProfileTreeFileSuffix = " <<                          trueProfileTreeFileSuffix << endl;
        os << "saveStartingProfiles = " <<                          saveStartingProfiles << endl;
        os << "startingProfileTreeFileSuffix = " <<                          startingProfileTreeFileSuffix << endl;
        os << "saveTestProfiles = " <<                          saveTestProfiles << endl;
        os << "testProfileTreeFileSuffix = " <<                          testProfileTreeFileSuffix << endl;
        os << "savePatternSequences = " <<                          savePatternSequences << endl;
        os << "patternSequencesFileSuffix = " <<                          patternSequencesFileSuffix << endl;
        os << "saveTests = " <<                          saveTests << endl;
        os << "testsFileSuffix = " <<                          testsFileSuffix << endl;
        os << "saveTrainingSequences = " <<                          saveTrainingSequences << endl;
        os << "trainingSequencesFileSuffix = " <<                          trainingSequencesFileSuffix << endl;
        os << "saveTestingSequences = " <<                          saveTestingSequences << endl;
        os << "testingSequencesFileSuffix = " <<                          testingSequencesFileSuffix << endl;
        os << "saveTrueTrainingAlignments = " <<                          saveTrueTrainingAlignments << endl;
        os << "trainingTrueAlignmentsFileSuffix = " <<                          trainingTrueAlignmentsFileSuffix << endl;
        os << "saveTrueTestingAlignments = " <<                          saveTrueTestingAlignments << endl;
        os << "trueTestingAlignmentsFileSuffix = " <<                          trueTestingAlignmentsFileSuffix << endl;
        os << "saveFileVersion = " <<                             saveFileVersion << endl;
        os << "numProfiles = " <<                                 numProfiles << endl;

        if( profileLengths.size() == 0 ) {
          os << "profileLengths = NULL" << endl;
        } else {
          os << "profileLengths = { ";
          for( uint32_t pl_i = 0; pl_i < profileLengths.size(); pl_i++ ) {
            if( pl_i > 0 ) {
              os << ", ";
            }
            os << profileLengths[ pl_i ];
          } // End foreach conservation rate..
          os << "}" << endl;
              } // End if profileLengths == NULL .. else ..
        os << "sharedPositionRate = " <<                          sharedPositionRate << endl;
        if( numTrainingSequencesPerProfiles.size() == 0 ) {
          os << "numTrainingSequencesPerProfiles = NULL" << endl;
              } else {
          os << "numTrainingSequencesPerProfiles = { ";
          for( uint32_t pl_i = 0; pl_i < numTrainingSequencesPerProfiles.size(); pl_i++ ) {
            if( pl_i > 0 ) {
              os << ", ";
            }
            os << numTrainingSequencesPerProfiles[ pl_i ];
          } // End foreach conservation rate..
          os << "}" << endl;
              } // End if numTrainingSequencesPerProfiles == NULL .. else ..
        os << "numTestingSequencesPerProfile = " <<                  numTestingSequencesPerProfile << endl;
        if( profileLengths.size() > 0 ) {
          os << "profileLengths = NULL" << endl;
              } else {
          os << "profileLengths = { ";
          for( uint32_t pl_i = 0; pl_i < profileLengths.size(); pl_i++ ) {
            if( pl_i > 0 ) {
              os << ", ";
            }
            os << profileLengths[ pl_i ];
          } // End foreach conservation rate..
          os << "}" << endl;
              } // End if profileLengths == NULL .. else ..

        if( expectedDeletionsCounts.size() == 0 ) {
          os << "expectedDeletionsCounts = NULL" << endl;
              } else {
          os << "expectedDeletionsCounts = { ";
          for( uint32_t cr_i = 0; cr_i < expectedDeletionsCounts.size(); cr_i++ ) {
            if( cr_i > 0 ) {
              os << ", ";
            }
            os << expectedDeletionsCounts[ cr_i ];
          } // End foreach conservation rate..
          os << "}" << endl;
              } // End if expectedDeletionsCounts == NULL .. else ..

        if( expectedInsertionsCounts.size() == 0 ) {
          os << "expectedInsertionsCounts = NULL" << endl;
              } else {
          os << "expectedInsertionsCounts = { ";
          for( uint32_t cr_i = 0; cr_i < expectedInsertionsCounts.size(); cr_i++ ) {
            if( cr_i > 0 ) {
              os << ", ";
            }
            os << expectedInsertionsCounts[ cr_i ];
          } // End foreach conservation rate..
          os << "}" << endl;
              } // End if expectedInsertionsCounts == NULL .. else ..

        if( expectedDeletionLengthAsProfileLengthFractions.size() == 0 ) {
          os << "expectedDeletionLengthAsProfileLengthFractions = NULL" << endl;
              } else {
          os << "expectedDeletionLengthAsProfileLengthFractions = { ";
          for( uint32_t cr_i = 0; cr_i < expectedDeletionLengthAsProfileLengthFractions.size(); cr_i++ ) {
            if( cr_i > 0 ) {
              os << ", ";
            }
            os << expectedDeletionLengthAsProfileLengthFractions[ cr_i ];
          } // End foreach conservation rate..
          os << "}" << endl;
              } // End if expectedDeletionLengthAsProfileLengthFractions == NULL .. else ..

        if( expectedInsertionLengthAsProfileLengthFractions.size() == 0 ) {
          os << "expectedInsertionLengthAsProfileLengthFractions = NULL" << endl;
              } else {
          os << "expectedInsertionLengthAsProfileLengthFractions = { ";
          for( uint32_t cr_i = 0; cr_i < expectedInsertionLengthAsProfileLengthFractions.size(); cr_i++ ) {
            if( cr_i > 0 ) {
              os << ", ";
            }
            os << expectedInsertionLengthAsProfileLengthFractions[ cr_i ];
          } // End foreach conservation rate..
          os << "}" << endl;
        } // End if expectedInsertionLengthAsProfileLengthFractions == NULL .. else ..

        os << "minExpectedDeletionLength = " <<                         minExpectedDeletionLength << endl;
        os << "minExpectedInsertionLength = " <<                         minExpectedInsertionLength << endl;
        os << "preAlignInsertion = " <<                         preAlignInsertion << endl;
        os << "postAlignInsertion = " <<                         postAlignInsertion << endl;
        os << "priorStrength = " <<                         priorStrength << endl;
        os << "priorStrength_internal_transitions = " <<                         priorStrength_internal_transitions << endl;
        os << "priorMtoM = " <<                         priorMtoM << endl;
        os << "priorMtoI = " <<                         priorMtoI << endl;
        os << "priorMtoD = " <<                         priorMtoD << endl;
        os << "priorItoM = " <<                         priorItoM << endl;
        os << "priorItoI = " <<                         priorItoI << endl;
        os << "priorDtoM = " <<                         priorDtoM << endl;
        os << "priorDtoD = " <<                         priorDtoD << endl;
        os << "reportGibbsMean = " <<                         reportGibbsMean << endl;
        os << "reportGibbsMode = " <<                         reportGibbsMode << endl;
        os << "numTrueProfiles = " <<                         numTrueProfiles << endl;
        os << "numStartingProfiles = " <<                         numStartingProfiles << endl;
        os << "startWithUniformGlobals = " <<                         startWithUniformGlobals << endl;
        os << "startWithUniformGlobals_scalar = " <<                         startWithUniformGlobals_scalar << endl;
        os << "startWithUniformGlobals_maxNtoN = " <<                         startWithUniformGlobals_maxNtoN << endl;
        os << "startWithUniformGlobals_maxBtoD = " <<                         startWithUniformGlobals_maxBtoD << endl;
        os << "startWithUniformGlobals_maxMtoI = " <<                         startWithUniformGlobals_maxMtoI << endl;
        os << "startWithUniformGlobals_maxMtoD = " <<                         startWithUniformGlobals_maxMtoD << endl;
        os << "startWithUniformGlobals_maxItoI = " <<                         startWithUniformGlobals_maxItoI << endl;
        os << "startWithUniformGlobals_maxDtoD = " <<                         startWithUniformGlobals_maxDtoD << endl;
        os << "startWithUniformGlobals_maxCtoC = " <<                         startWithUniformGlobals_maxCtoC << endl;
        os << "startWithUniformPositions = " <<                         startWithUniformPositions << endl;
        os << "startWithGlobalsDrawnFromPrior = " <<                         startWithGlobalsDrawnFromPrior << endl;
        os << "startWithPositionsDrawnFromPrior = " <<                         startWithPositionsDrawnFromPrior << endl;
        os << "testViterbi = " <<                                  testViterbi << endl;
        os << "coutViterbi = " <<                                  coutViterbi << endl;
        os << "testTruepath = " <<                                  testTruepath << endl;
        os << "coutTruepath = " <<                                  coutTruepath << endl;
        os << "calculateSymmeterizedKullbackLeiblerDistancesToTrue = " <<       calculateSymmeterizedKullbackLeiblerDistancesToTrue << endl;
        os << "calculateSymmeterizedKullbackLeiblerDistancesToStarting = " <<       calculateSymmeterizedKullbackLeiblerDistancesToStarting << endl;
        os << "coutDistances = " <<    coutDistances << endl;
        os << "calculateProfileProfileAlignments = " <<  calculateProfileProfileAlignments << endl;
        os << "profileProfileIndelOpenCost = " <<  profileProfileIndelOpenCost << endl;
        os << "profileProfileIndelExtensionCost = " <<  profileProfileIndelExtensionCost << endl;
        os << "testTrueProfile = " <<                         testTrueProfile << endl;
        os << "coutTrueProfile = " <<                         coutTrueProfile << endl;
        os << "testStartingProfile = " <<                         testStartingProfile << endl;
        os << "coutStartingProfile = " <<                         coutStartingProfile << endl;
        os << "testUnconditionalProfile = " <<                    testUnconditionalProfile << endl;
        os << "coutUnconditionalProfile = " <<                    coutUnconditionalProfile << endl;
        os << "testUnconditionalWithFixedStartingGlobalsProfile = " <<                    testUnconditionalWithFixedStartingGlobalsProfile << endl;
        os << "coutUnconditionalWithFixedStartingGlobalsProfile = " <<                    coutUnconditionalWithFixedStartingGlobalsProfile << endl;
        os << "testUnconditionalWithFixedTrueGlobalsProfile = " <<                    testUnconditionalWithFixedTrueGlobalsProfile << endl;
        os << "coutUnconditionalWithFixedTrueGlobalsProfile = " <<                    coutUnconditionalWithFixedTrueGlobalsProfile << endl;
        os << "testConditionalThenUnconditionalProfile = " <<     testConditionalThenUnconditionalProfile << endl;
        os << "coutConditionalThenUnconditionalProfile = " <<     coutConditionalThenUnconditionalProfile << endl;
        os << "testUnconditionalThenConditionalProfile = " <<     testUnconditionalThenConditionalProfile << endl;
        os << "coutUnconditionalThenConditionalProfile = " <<     coutUnconditionalThenConditionalProfile << endl;
        os << "testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile = " << testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile << endl;
        os << "coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile = " << coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile << endl;
        os << "testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile = " << testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile << endl;
        os << "coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile = " << coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile << endl;
        os << "testConditionalGibbsProfile = " <<                    testConditionalGibbsProfile << endl;
        os << "coutConditionalGibbsProfile = " <<                    coutConditionalGibbsProfile << endl;
        os << "testUnconditionalGibbsProfile = " <<                    testUnconditionalGibbsProfile << endl;
        os << "coutUnconditionalGibbsProfile = " <<                    coutUnconditionalGibbsProfile << endl;
        os << "testLengthadjust = " << testLengthadjust << endl;
        os << "testBaldi = " << testBaldi << endl;
        os << "testBaldiSiegel = " << testBaldiSiegel << endl;
        os << "alsoStartWithEvenPositions = " << alsoStartWithEvenPositions << endl;
      } // writeParameters ( basic_ostream & )

  ////// Class galosh::ProfuseTest::ParametersModifierTemplate ////
  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class ParametersType>
  GALOSH_INLINE_INIT
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      ParametersModifierTemplate ()
      {
        if( base_parameters_modifier_t::parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProfuseTest::ParametersModifierTemplate::<init>()" << endl;
        } // End if DEBUG_All
        isModified_reset();
      } // <init>()

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_INIT
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      // Copy constructor
      ParametersModifierTemplate ( const AnyParametersModifierTemplate & copy_from )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProfuseTest::ParametersModifierTemplate::<init>( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
      } // <init>( AnyParametersModifierTemplate const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_TRIVIAL
  typename ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::template ParametersModifierTemplate<ParametersType> &
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      // Copy constructor/operator
  operator= (
        const AnyParametersModifierTemplate & copy_from
      )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProfuseTest::ParametersModifierTemplate::operator=( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
        return *this;
      } // operator=( AnyParametersModifierTemplate const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_COPY
  void
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProfuseTest::ParametersModifierTemplate::copyFromNonVirtual( copy_from )" << endl;
        } // End if DEBUG_All

        isModified_copyFromNonVirtual( copy_from );

        base_parameters_modifier_t::parameters.copyFromNonVirtual( copy_from.parameters );
      } // copyFromNonVirtual( AnyParametersModifierTemplate const & )


  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_COPY
  void
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      isModified_copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      )
      {
        base_parameters_modifier_t::isModified_copyFromNonVirtual( copy_from );

        isModified_saveResultsToFile =                            copy_from.isModified_saveResultsToFile;
        isModified_saveResultsParentDirectory =                           copy_from.isModified_saveResultsParentDirectory;
        isModified_resultsFilePrefix =                           copy_from.isModified_resultsFilePrefix;
        isModified_tabFileSuffix =                           copy_from.isModified_tabFileSuffix;
        isModified_parametersFileSuffix =                           copy_from.isModified_parametersFileSuffix;
        isModified_saveTrueProfileTrees =                           copy_from.isModified_saveTrueProfileTrees;
        isModified_trueProfileTreeFileSuffix =                           copy_from.isModified_trueProfileTreeFileSuffix;
        isModified_saveStartingProfiles =                           copy_from.isModified_saveStartingProfiles;
        isModified_startingProfileTreeFileSuffix =                           copy_from.isModified_startingProfileTreeFileSuffix;
        isModified_saveTestProfiles =                           copy_from.isModified_saveTestProfiles;
        isModified_testProfileTreeFileSuffix =                           copy_from.isModified_testProfileTreeFileSuffix;
        isModified_savePatternSequences =                           copy_from.isModified_savePatternSequences;
        isModified_patternSequencesFileSuffix =                           copy_from.isModified_patternSequencesFileSuffix;
        isModified_saveTests =                           copy_from.isModified_saveTests;
        isModified_testsFileSuffix =                           copy_from.isModified_testsFileSuffix;
        isModified_saveTrainingSequences =                           copy_from.isModified_saveTrainingSequences;
        isModified_trainingSequencesFileSuffix =                           copy_from.isModified_trainingSequencesFileSuffix;
        isModified_saveTestingSequences =                           copy_from.isModified_saveTestingSequences;
        isModified_testingSequencesFileSuffix =                           copy_from.isModified_testingSequencesFileSuffix;
        isModified_saveTrueTrainingAlignments =                           copy_from.isModified_saveTrueTrainingAlignments;
        isModified_trainingTrueAlignmentsFileSuffix =                           copy_from.isModified_trainingTrueAlignmentsFileSuffix;
        isModified_saveTrueTestingAlignments =                           copy_from.isModified_saveTrueTestingAlignments;
        isModified_trueTestingAlignmentsFileSuffix =                           copy_from.isModified_trueTestingAlignmentsFileSuffix;
        isModified_saveFileVersion =                              copy_from.isModified_saveFileVersion;
        isModified_numProfiles =                                  copy_from.isModified_numProfiles;
        isModified_profileLengths =                                copy_from.isModified_profileLengths;
        isModified_sharedPositionRate =                           copy_from.isModified_sharedPositionRate;
        isModified_numTrainingSequencesPerProfiles =               copy_from.isModified_numTrainingSequencesPerProfiles;
        isModified_numTestingSequencesPerProfile =                   copy_from.isModified_numTestingSequencesPerProfile;
        isModified_conservationRates =                             copy_from.isModified_conservationRates;
        isModified_useDeletionsForInsertionsParameters =                             copy_from.isModified_useDeletionsForInsertionsParameters;
        isModified_expectedDeletionsCounts =                          copy_from.isModified_expectedDeletionsCounts;
        isModified_expectedInsertionsCounts =                          copy_from.isModified_expectedInsertionsCounts;
        isModified_expectedDeletionLengthAsProfileLengthFractions =                          copy_from.isModified_expectedDeletionLengthAsProfileLengthFractions;
        isModified_expectedInsertionLengthAsProfileLengthFractions =                          copy_from.isModified_expectedInsertionLengthAsProfileLengthFractions;
        isModified_minExpectedDeletionLength =                          copy_from.isModified_minExpectedDeletionLength;
        isModified_minExpectedInsertionLength =                          copy_from.isModified_minExpectedInsertionLength;
        isModified_preAlignInsertion =                          copy_from.isModified_preAlignInsertion;
        isModified_postAlignInsertion =                          copy_from.isModified_postAlignInsertion;
        isModified_priorStrength =                          copy_from.isModified_priorStrength;
        isModified_priorStrength_internal_transitions =                          copy_from.isModified_priorStrength_internal_transitions;
        isModified_priorMtoM =                          copy_from.isModified_priorMtoM;
        isModified_priorMtoI =                          copy_from.isModified_priorMtoI;
        isModified_priorMtoD =                          copy_from.isModified_priorMtoD;
        isModified_priorItoM =                          copy_from.isModified_priorItoM;
        isModified_priorItoI =                          copy_from.isModified_priorItoI;
        isModified_priorDtoM =                          copy_from.isModified_priorDtoM;
        isModified_priorDtoD =                          copy_from.isModified_priorDtoD;
        isModified_reportGibbsMean =                          copy_from.isModified_reportGibbsMean;
        isModified_reportGibbsMode =                          copy_from.isModified_reportGibbsMode;
        isModified_numTrueProfiles =                          copy_from.isModified_numTrueProfiles;
        isModified_numStartingProfiles =                          copy_from.isModified_numStartingProfiles;
        isModified_startWithUniformGlobals =                          copy_from.isModified_startWithUniformGlobals;
        isModified_startWithUniformGlobals_scalar =                          copy_from.isModified_startWithUniformGlobals_scalar;
        isModified_startWithUniformGlobals_maxNtoN =                          copy_from.isModified_startWithUniformGlobals_maxNtoN;
        isModified_startWithUniformGlobals_maxBtoD =                          copy_from.isModified_startWithUniformGlobals_maxBtoD;
        isModified_startWithUniformGlobals_maxMtoI =                          copy_from.isModified_startWithUniformGlobals_maxMtoI;
        isModified_startWithUniformGlobals_maxMtoD =                          copy_from.isModified_startWithUniformGlobals_maxMtoD;
        isModified_startWithUniformGlobals_maxItoI =                          copy_from.isModified_startWithUniformGlobals_maxItoI;
        isModified_startWithUniformGlobals_maxDtoD =                          copy_from.isModified_startWithUniformGlobals_maxDtoD;
        isModified_startWithUniformGlobals_maxCtoC =                          copy_from.isModified_startWithUniformGlobals_maxCtoC;
        isModified_startWithUniformPositions =                          copy_from.isModified_startWithUniformPositions;
        isModified_startWithGlobalsDrawnFromPrior =                          copy_from.isModified_startWithGlobalsDrawnFromPrior;
        isModified_startWithPositionsDrawnFromPrior =                          copy_from.isModified_startWithPositionsDrawnFromPrior;
        isModified_testViterbi =                                   copy_from.isModified_testViterbi;
        isModified_coutViterbi =                                   copy_from.isModified_coutViterbi;
        isModified_testTruepath =                                   copy_from.isModified_testTruepath;
        isModified_coutTruepath =                                   copy_from.isModified_coutTruepath;
        isModified_calculateSymmeterizedKullbackLeiblerDistancesToTrue =        copy_from.isModified_calculateSymmeterizedKullbackLeiblerDistancesToTrue;
        isModified_calculateSymmeterizedKullbackLeiblerDistancesToStarting =        copy_from.isModified_calculateSymmeterizedKullbackLeiblerDistancesToStarting;
        isModified_coutDistances =     copy_from.isModified_coutDistances;
        isModified_calculateProfileProfileAlignments =   copy_from.isModified_calculateProfileProfileAlignments;
        isModified_profileProfileIndelOpenCost =   copy_from.isModified_profileProfileIndelOpenCost;
        isModified_profileProfileIndelExtensionCost =   copy_from.isModified_profileProfileIndelExtensionCost;
        isModified_testTrueProfile =                          copy_from.isModified_testTrueProfile;
        isModified_coutTrueProfile =                          copy_from.isModified_coutTrueProfile;
        isModified_testStartingProfile =                          copy_from.isModified_testStartingProfile;
        isModified_coutStartingProfile =                          copy_from.isModified_coutStartingProfile;
        isModified_testUnconditionalProfile =                     copy_from.isModified_testUnconditionalProfile;
        isModified_coutUnconditionalProfile =                     copy_from.isModified_coutUnconditionalProfile;
        isModified_testUnconditionalWithFixedStartingGlobalsProfile =                     copy_from.isModified_testUnconditionalWithFixedStartingGlobalsProfile;
        isModified_coutUnconditionalWithFixedStartingGlobalsProfile =                     copy_from.isModified_coutUnconditionalWithFixedStartingGlobalsProfile;
        isModified_testUnconditionalWithFixedTrueGlobalsProfile =                     copy_from.isModified_testUnconditionalWithFixedTrueGlobalsProfile;
        isModified_coutUnconditionalWithFixedTrueGlobalsProfile =                     copy_from.isModified_coutUnconditionalWithFixedTrueGlobalsProfile;
        isModified_testConditionalThenUnconditionalProfile =      copy_from.isModified_testConditionalThenUnconditionalProfile;
        isModified_coutConditionalThenUnconditionalProfile =      copy_from.isModified_coutConditionalThenUnconditionalProfile;
        isModified_testUnconditionalThenConditionalProfile =      copy_from.isModified_testUnconditionalThenConditionalProfile;
        isModified_coutUnconditionalThenConditionalProfile =      copy_from.isModified_coutUnconditionalThenConditionalProfile;
        isModified_testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile =  copy_from.isModified_testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile;
        isModified_coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile =  copy_from.isModified_coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile;
        isModified_testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile =  copy_from.isModified_testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile;
        isModified_coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile =  copy_from.isModified_coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile;
        isModified_testConditionalGibbsProfile =                     copy_from.isModified_testConditionalGibbsProfile;
        isModified_coutConditionalGibbsProfile =                     copy_from.isModified_coutConditionalGibbsProfile;
        isModified_testUnconditionalGibbsProfile =                     copy_from.isModified_testUnconditionalGibbsProfile;
        isModified_coutUnconditionalGibbsProfile =                     copy_from.isModified_coutUnconditionalGibbsProfile;
        isModified_testLengthadjust = copy_from.isModified_testLengthadjust;
        isModified_testBaldi = copy_from.isModified_testBaldi;
        isModified_testBaldiSiegel = copy_from.isModified_testBaldiSiegel;
        isModified_alsoStartWithEvenPositions = copy_from.isModified_alsoStartWithEvenPositions;
      } // isModified_copyFromNonVirtual( AnyParametersModifierTemplate const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class ParametersType>
  GALOSH_INLINE_TRIVIAL
  void
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      reset ()
      {
        isModified_reset();
        base_parameters_modifier_t::parameters.resetToDefaults();
      } // reset()

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class ParametersType>
  GALOSH_INLINE_REINITIALIZE
  void
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      isModified_reset ()
      {
        base_parameters_modifier_t::isModified_reset();

        isModified_saveResultsToFile = false;
        isModified_saveResultsParentDirectory = false;
        isModified_resultsFilePrefix = false;
        isModified_tabFileSuffix = false;
        isModified_parametersFileSuffix = false;
        isModified_saveTrueProfileTrees = false;
        isModified_trueProfileTreeFileSuffix = false;
        isModified_saveStartingProfiles = false;
        isModified_startingProfileTreeFileSuffix = false;
        isModified_saveTestProfiles = false;
        isModified_testProfileTreeFileSuffix = false;
        isModified_savePatternSequences = false;
        isModified_patternSequencesFileSuffix = false;
        isModified_saveTests = false;
        isModified_patternSequencesFileSuffix = false;
        isModified_saveTrainingSequences = false;
        isModified_trainingSequencesFileSuffix = false;
        isModified_saveTestingSequences = false;
        isModified_testingSequencesFileSuffix = false;
        isModified_saveTrueTrainingAlignments = false;
        isModified_trainingTrueAlignmentsFileSuffix = false;
        isModified_saveTrueTestingAlignments = false;
        isModified_trueTestingAlignmentsFileSuffix = false;
        isModified_saveFileVersion = false;
        isModified_numProfiles = false;
        isModified_profileLengths = false;
        isModified_sharedPositionRate = false;
        isModified_numTrainingSequencesPerProfiles = false;
        isModified_numTestingSequencesPerProfile = false;
        isModified_conservationRates = false;
        isModified_useDeletionsForInsertionsParameters = false;
        isModified_expectedDeletionsCounts = false;
        isModified_expectedInsertionsCounts = false;
        isModified_expectedDeletionLengthAsProfileLengthFractions = false;
        isModified_expectedInsertionLengthAsProfileLengthFractions = false;
        isModified_minExpectedDeletionLength = false;
        isModified_minExpectedInsertionLength = false;
        isModified_preAlignInsertion = false;
        isModified_postAlignInsertion = false;
        isModified_priorStrength = false;
        isModified_priorStrength_internal_transitions = false;
        isModified_priorMtoM = false;
        isModified_priorMtoI = false;
        isModified_priorMtoD = false;
        isModified_priorItoM = false;
        isModified_priorItoI = false;
        isModified_priorDtoM = false;
        isModified_priorDtoD = false;
        isModified_reportGibbsMean = false;
        isModified_reportGibbsMode = false;
        isModified_numTrueProfiles = false;
        isModified_numStartingProfiles = false;
        isModified_startWithUniformGlobals = false;
        isModified_startWithUniformGlobals_scalar = false;
        isModified_startWithUniformGlobals_maxNtoN =  false;
        isModified_startWithUniformGlobals_maxBtoD =  false;
        isModified_startWithUniformGlobals_maxMtoI =  false;
        isModified_startWithUniformGlobals_maxMtoD =  false;
        isModified_startWithUniformGlobals_maxItoI =  false;
        isModified_startWithUniformGlobals_maxDtoD =  false;
        isModified_startWithUniformGlobals_maxCtoC =  false;
        isModified_startWithUniformPositions =        false;
        isModified_startWithGlobalsDrawnFromPrior =   false;
        isModified_startWithPositionsDrawnFromPrior = false;
        isModified_testViterbi = false;
        isModified_coutViterbi = false;
        isModified_testTruepath = false;
        isModified_coutTruepath = false;
        isModified_calculateSymmeterizedKullbackLeiblerDistancesToTrue = false;
        isModified_calculateSymmeterizedKullbackLeiblerDistancesToStarting = false;
        isModified_coutDistances = false;
        isModified_calculateProfileProfileAlignments = false;
        isModified_profileProfileIndelOpenCost = false;
        isModified_profileProfileIndelExtensionCost = false;
        isModified_testTrueProfile = false;
        isModified_coutTrueProfile = false;
        isModified_testStartingProfile = false;
        isModified_coutStartingProfile = false;
        isModified_testUnconditionalProfile = false;
        isModified_coutUnconditionalProfile = false;
        isModified_testUnconditionalWithFixedStartingGlobalsProfile = false;
        isModified_coutUnconditionalWithFixedStartingGlobalsProfile = false;
        isModified_testUnconditionalWithFixedTrueGlobalsProfile = false;
        isModified_coutUnconditionalWithFixedTrueGlobalsProfile = false;
        isModified_testConditionalThenUnconditionalProfile = false;
        isModified_coutConditionalThenUnconditionalProfile = false;
        isModified_testUnconditionalThenConditionalProfile = false;
        isModified_coutUnconditionalThenConditionalProfile = false;
        isModified_testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile = false;
        isModified_coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile = false;
        isModified_testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile = false;
        isModified_coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile = false;
        isModified_testConditionalGibbsProfile = false;
        isModified_coutConditionalGibbsProfile = false;
        isModified_testUnconditionalGibbsProfile = false;
        isModified_coutUnconditionalGibbsProfile = false;
        isModified_testLengthadjust = false;
        isModified_testBaldi = false;
        isModified_testBaldiSiegel = false;
        isModified_alsoStartWithEvenPositions = false;
      } // isModified_reset()

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType,
            typename SequenceResidueType>
  template <class ParametersType>
  template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
  void
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      writeParametersModifier (
        std::basic_ostream<CharT,Traits>& os
      ) const
      {
        //base_parameters_modifier_t::operator<<( os, parameters_modifier );
        base_parameters_modifier_t::writeParametersModifier( os );
        os << endl;

        os << "[ProfuseTest]" << endl;
        if( isModified_saveResultsToFile ) {
          os << "saveResultsToFile = " <<                           base_parameters_modifier_t::parameters.saveResultsToFile << endl;
        }
        if( isModified_saveResultsParentDirectory ) {
          os << "saveResultsParentDirectory = " <<                          base_parameters_modifier_t::parameters.saveResultsParentDirectory << endl;
        }
        if( isModified_resultsFilePrefix ) {
          os << "resultsFilePrefix = " <<                          base_parameters_modifier_t::parameters.resultsFilePrefix << endl;
        }
        if( isModified_tabFileSuffix ) {
          os << "tabFileSuffix = " <<                          base_parameters_modifier_t::parameters.tabFileSuffix << endl;
        }
        if( isModified_parametersFileSuffix ) {
          os << "parametersFileSuffix = " <<                          base_parameters_modifier_t::parameters.parametersFileSuffix << endl;
        }
        if( isModified_saveTrueProfileTrees ) {
          os << "saveTrueProfileTrees = " <<                          base_parameters_modifier_t::parameters.saveTrueProfileTrees << endl;
        }
        if( isModified_trueProfileTreeFileSuffix ) {
          os << "trueProfileTreeFileSuffix = " <<                          base_parameters_modifier_t::parameters.trueProfileTreeFileSuffix << endl;
        }
        if( isModified_saveTrueProfileTrees ) {
          os << "saveStartingProfiles = " <<                          base_parameters_modifier_t::parameters.saveStartingProfiles << endl;
        }
        if( isModified_startingProfileTreeFileSuffix ) {
          os << "startingProfileTreeFileSuffix = " <<                          base_parameters_modifier_t::parameters.startingProfileTreeFileSuffix << endl;
        }
        if( isModified_saveTestProfiles ) {
          os << "saveTestProfiles = " <<                          base_parameters_modifier_t::parameters.saveTestProfiles << endl;
        }
        if( isModified_testProfileTreeFileSuffix ) {
          os << "testProfileTreeFileSuffix = " <<                          base_parameters_modifier_t::parameters.testProfileTreeFileSuffix << endl;
        }
        if( isModified_savePatternSequences ) {
          os << "savePatternSequences = " <<                          base_parameters_modifier_t::parameters.savePatternSequences << endl;
        }
        if( isModified_testsFileSuffix ) {
          os << "testsFileSuffix = " <<                          base_parameters_modifier_t::parameters.patternSequencesFileSuffix << endl;
        }
        if( isModified_saveTests ) {
          os << "saveTests = " <<                          base_parameters_modifier_t::parameters.saveTests << endl;
        }
        if( isModified_patternSequencesFileSuffix ) {
          os << "testsFileSuffix = " <<                          base_parameters_modifier_t::parameters.testsFileSuffix << endl;
        }
        if( isModified_saveTrainingSequences ) {
          os << "saveTrainingSequences = " <<                          base_parameters_modifier_t::parameters.saveTrainingSequences << endl;
        }
        if( isModified_trainingSequencesFileSuffix ) {
          os << "trainingSequencesFileSuffix = " <<                          base_parameters_modifier_t::parameters.trainingSequencesFileSuffix << endl;
        }
        if( isModified_saveTestingSequences ) {
          os << "saveTestingSequences = " <<                          base_parameters_modifier_t::parameters.saveTestingSequences << endl;
        }
        if( isModified_testingSequencesFileSuffix ) {
          os << "testingSequencesFileSuffix = " <<                          base_parameters_modifier_t::parameters.testingSequencesFileSuffix << endl;
        }
        if( isModified_saveTrueTrainingAlignments ) {
          os << "saveTrueTrainingAlignments = " <<                          base_parameters_modifier_t::parameters.saveTrueTrainingAlignments << endl;
        }
        if( isModified_trainingTrueAlignmentsFileSuffix ) {
          os << "trainingTrueAlignmentsFileSuffix = " <<                          base_parameters_modifier_t::parameters.trainingTrueAlignmentsFileSuffix << endl;
        }
        if( isModified_saveTrueTestingAlignments ) {
          os << "saveTrueTestingAlignments = " <<                          base_parameters_modifier_t::parameters.saveTrueTestingAlignments << endl;
        }
        if( isModified_trueTestingAlignmentsFileSuffix ) {
          os << "trueTestingAlignmentsFileSuffix = " <<                          base_parameters_modifier_t::parameters.trueTestingAlignmentsFileSuffix << endl;
        }
        if( isModified_saveFileVersion ) {
          os << "saveFileVersion = " <<                             base_parameters_modifier_t::parameters.saveFileVersion << endl;
        }
        if( isModified_numProfiles ) {
          os << "numProfiles = " <<                                 base_parameters_modifier_t::parameters.numProfiles << endl;
        }
        if( isModified_profileLengths ) {
          if( base_parameters_modifier_t::parameters.profileLengths == NULL ) {
            os << "profileLengths = NULL" << endl;
          } else {
            os << "profileLengths = { ";
            for( uint32_t pl_i = 0; pl_i < base_parameters_modifier_t::parameters.profileLengths->size(); pl_i++ ) {
              if( pl_i > 0 ) {
                os << ", ";
              }
              os << ( *base_parameters_modifier_t::parameters.profileLengths )[ pl_i ];
            } // End foreach conservation rate..
            os << "}" << endl;
          } // End if profileLengths == NULL .. else ..
        }
        if( isModified_sharedPositionRate ) {
          os << "sharedPositionRate = " <<                          base_parameters_modifier_t::parameters.sharedPositionRate << endl;
        }
        if( isModified_numTrainingSequencesPerProfiles ) {
          if( base_parameters_modifier_t::parameters.numTrainingSequencesPerProfiles == NULL ) {
            os << "numTrainingSequencesPerProfiles = NULL" << endl;
          } else {
            os << "numTrainingSequencesPerProfiles = { ";
            for( uint32_t pl_i = 0; pl_i < base_parameters_modifier_t::parameters.numTrainingSequencesPerProfiles->size(); pl_i++ ) {
              if( pl_i > 0 ) {
                os << ", ";
              }
              os << ( *base_parameters_modifier_t::parameters.numTrainingSequencesPerProfiles )[ pl_i ];
            } // End foreach conservation rate..
            os << "}" << endl;
          } // End if numTrainingSequencesPerProfiles == NULL .. else ..
        }
        if( isModified_numTestingSequencesPerProfile ) {
          os << "numTestingSequencesPerProfile = " <<                  base_parameters_modifier_t::parameters.numTestingSequencesPerProfile << endl;
        }
        if( isModified_conservationRates ) {
          if( base_parameters_modifier_t::parameters.conservationRates == NULL ) {
            os << "conservationRates = NULL" << endl;
          } else {
            os << "conservationRates = { ";
            for( uint32_t cr_i = 0; cr_i < base_parameters_modifier_t::parameters.conservationRates->size(); cr_i++ ) {
              if( cr_i > 0 ) {
                os << ", ";
              }
              os << ( *base_parameters_modifier_t::parameters.conservationRates )[ cr_i ];
            } // End foreach conservation rate..
            os << "}" << endl;
          } // End if conservationRates == NULL .. else ..
        }
        if( isModified_expectedDeletionsCounts ) {
          if( base_parameters_modifier_t::parameters.expectedDeletionsCounts == NULL ) {
            os << "expectedDeletionsCounts = NULL" << endl;
          } else {
            os << "expectedDeletionsCounts = { ";
            for( uint32_t cr_i = 0; cr_i < base_parameters_modifier_t::parameters.expectedDeletionsCounts->size(); cr_i++ ) {
              if( cr_i > 0 ) {
                os << ", ";
              }
              os << ( *base_parameters_modifier_t::parameters.expectedDeletionsCounts )[ cr_i ];
            } // End foreach conservation rate..
            os << "}" << endl;
          } // End if expectedDeletionsCounts == NULL .. else ..
        }
        if( isModified_expectedInsertionsCounts ) {
          if( base_parameters_modifier_t::parameters.expectedInsertionsCounts == NULL ) {
            os << "expectedInsertionsCounts = NULL" << endl;
          } else {
            os << "expectedInsertionsCounts = { ";
            for( uint32_t cr_i = 0; cr_i < base_parameters_modifier_t::parameters.expectedInsertionsCounts->size(); cr_i++ ) {
              if( cr_i > 0 ) {
                os << ", ";
              }
              os << ( *base_parameters_modifier_t::parameters.expectedInsertionsCounts )[ cr_i ];
            } // End foreach conservation rate..
            os << "}" << endl;
          } // End if expectedInsertionsCounts == NULL .. else ..
        }
        if( isModified_expectedDeletionLengthAsProfileLengthFractions ) {
          if( base_parameters_modifier_t::parameters.expectedDeletionLengthAsProfileLengthFractions == NULL ) {
            os << "expectedDeletionLengthAsProfileLengthFractions = NULL" << endl;
          } else {
            os << "expectedDeletionLengthAsProfileLengthFractions = { ";
            for( uint32_t cr_i = 0; cr_i < base_parameters_modifier_t::parameters.expectedDeletionLengthAsProfileLengthFractions->size(); cr_i++ ) {
              if( cr_i > 0 ) {
                os << ", ";
              }
              os << ( *base_parameters_modifier_t::parameters.expectedDeletionLengthAsProfileLengthFractions )[ cr_i ];
            } // End foreach conservation rate..
            os << "}" << endl;
          } // End if expectedDeletionLengthAsProfileLengthFractions == NULL .. else ..
        }
        if( isModified_expectedInsertionLengthAsProfileLengthFractions ) {
          if( base_parameters_modifier_t::parameters.expectedInsertionLengthAsProfileLengthFractions == NULL ) {
            os << "expectedInsertionLengthAsProfileLengthFractions = NULL" << endl;
          } else {
            os << "expectedInsertionLengthAsProfileLengthFractions = { ";
            for( uint32_t cr_i = 0; cr_i < base_parameters_modifier_t::parameters.expectedInsertionLengthAsProfileLengthFractions->size(); cr_i++ ) {
              if( cr_i > 0 ) {
                os << ", ";
              }
              os << ( *base_parameters_modifier_t::parameters.expectedInsertionLengthAsProfileLengthFractions )[ cr_i ];
            } // End foreach conservation rate..
            os << "}" << endl;
          } // End if expectedInsertionLengthAsProfileLengthFractions == NULL .. else ..
        }
        if( isModified_minExpectedDeletionLength ) {
          os << "minExpectedDeletionLength = " <<                         base_parameters_modifier_t::parameters.minExpectedDeletionLength << endl;
        }
        if( isModified_minExpectedInsertionLength ) {
          os << "minExpectedInsertionLength = " <<                         base_parameters_modifier_t::parameters.minExpectedInsertionLength << endl;
        }
        if( isModified_preAlignInsertion ) {
          os << "preAlignInsertion = " <<                         base_parameters_modifier_t::parameters.preAlignInsertion << endl;
        }
        if( isModified_postAlignInsertion ) {
          os << "postAlignInsertion = " <<                         base_parameters_modifier_t::parameters.postAlignInsertion << endl;
        }
        if( isModified_priorStrength ) {
          os << "priorStrength = " <<                         base_parameters_modifier_t::parameters.priorStrength << endl;
        }
        if( isModified_priorStrength_internal_transitions ) {
          os << "priorStrength_internal_transitions = " <<                         base_parameters_modifier_t::parameters.priorStrength_internal_transitions << endl;
        }
        if( isModified_priorMtoM ) {
          os << "priorMtoM = " <<                         base_parameters_modifier_t::parameters.priorMtoM << endl;
        }
        if( isModified_priorMtoI ) {
          os << "priorMtoI = " <<                         base_parameters_modifier_t::parameters.priorMtoI << endl;
        }
        if( isModified_priorMtoD ) {
          os << "priorMtoD = " <<                         base_parameters_modifier_t::parameters.priorMtoD << endl;
        }
        if( isModified_priorItoM ) {
          os << "priorItoM = " <<                         base_parameters_modifier_t::parameters.priorItoM << endl;
        }
        if( isModified_priorItoI ) {
          os << "priorItoI = " <<                         base_parameters_modifier_t::parameters.priorItoI << endl;
        }
        if( isModified_priorDtoM ) {
          os << "priorDtoM = " <<                         base_parameters_modifier_t::parameters.priorDtoM << endl;
        }
        if( isModified_priorDtoD ) {
          os << "priorDtoD = " <<                         base_parameters_modifier_t::parameters.priorDtoD << endl;
        }
        if( isModified_reportGibbsMean ) {
          os << "reportGibbsMean = " <<                         base_parameters_modifier_t::parameters.reportGibbsMean << endl;
        }
        if( isModified_reportGibbsMode ) {
          os << "reportGibbsMode = " <<                         base_parameters_modifier_t::parameters.reportGibbsMode << endl;
        }
        if( isModified_numTrueProfiles ) {
          os << "numTrueProfiles = " <<                         base_parameters_modifier_t::parameters.numTrueProfiles << endl;
        }
        if( isModified_numStartingProfiles ) {
          os << "numStartingProfiles = " <<                         base_parameters_modifier_t::parameters.numStartingProfiles << endl;
        }
        if( isModified_startWithUniformGlobals ) {
          os << "startWithUniformGlobals = " <<                         base_parameters_modifier_t::parameters.startWithUniformGlobals << endl;
        }
        if( isModified_startWithUniformGlobals_scalar ) {
          os << "startWithUniformGlobals_scalar = " <<                         base_parameters_modifier_t::parameters.startWithUniformGlobals_scalar << endl;
        }
        if( isModified_startWithUniformGlobals_maxNtoN ) {
          os << "startWithUniformGlobals_maxNtoN = " <<                         base_parameters_modifier_t::parameters.startWithUniformGlobals_maxNtoN << endl;
        }
        if( isModified_startWithUniformGlobals_maxBtoD ) {
          os << "startWithUniformGlobals_maxBtoD = " <<                         base_parameters_modifier_t::parameters.startWithUniformGlobals_maxBtoD << endl;
        }
        if( isModified_startWithUniformGlobals_maxMtoI ) {
          os << "startWithUniformGlobals_maxMtoI = " <<                         base_parameters_modifier_t::parameters.startWithUniformGlobals_maxMtoI << endl;
        }
        if( isModified_startWithUniformGlobals_maxMtoD ) {
          os << "startWithUniformGlobals_maxMtoD = " <<                         base_parameters_modifier_t::parameters.startWithUniformGlobals_maxMtoD << endl;
        }
        if( isModified_startWithUniformGlobals_maxItoI ) {
          os << "startWithUniformGlobals_maxItoI = " <<                         base_parameters_modifier_t::parameters.startWithUniformGlobals_maxItoI << endl;
        }
        if( isModified_startWithUniformGlobals_maxDtoD ) {
          os << "startWithUniformGlobals_maxDtoD = " <<                         base_parameters_modifier_t::parameters.startWithUniformGlobals_maxDtoD << endl;
        }
        if( isModified_startWithUniformGlobals_maxCtoC ) {
          os << "startWithUniformGlobals_maxCtoC = " <<                         base_parameters_modifier_t::parameters.startWithUniformGlobals_maxCtoC << endl;
        }
        if( isModified_startWithUniformPositions ) {
          os << "startWithUniformPositions = " <<                         base_parameters_modifier_t::parameters.startWithUniformPositions << endl;
        }
        if( isModified_startWithGlobalsDrawnFromPrior ) {
          os << "startWithGlobalsDrawnFromPrior = " <<                         base_parameters_modifier_t::parameters.startWithGlobalsDrawnFromPrior << endl;
        }
        if( isModified_startWithPositionsDrawnFromPrior ) {
          os << "startWithPositionsDrawnFromPrior = " <<                         base_parameters_modifier_t::parameters.startWithPositionsDrawnFromPrior << endl;
        }
        if( isModified_testViterbi ) {
          os << "testViterbi = " <<                                  base_parameters_modifier_t::parameters.testViterbi << endl;
        }
        if( isModified_coutViterbi ) {
          os << "coutViterbi = " <<                                  base_parameters_modifier_t::parameters.coutViterbi << endl;
        }
        if( isModified_testTruepath ) {
          os << "testTruepath = " <<                                  base_parameters_modifier_t::parameters.testTruepath << endl;
        }
        if( isModified_coutTruepath ) {
          os << "coutTruepath = " <<                                  base_parameters_modifier_t::parameters.coutTruepath << endl;
        }
        if( isModified_calculateSymmeterizedKullbackLeiblerDistancesToTrue ) {
          os << "calculateSymmeterizedKullbackLeiblerDistancesToTrue = " <<       base_parameters_modifier_t::parameters.calculateSymmeterizedKullbackLeiblerDistancesToTrue << endl;
        }
        if( isModified_calculateSymmeterizedKullbackLeiblerDistancesToStarting ) {
          os << "calculateSymmeterizedKullbackLeiblerDistancesToStarting = " <<       base_parameters_modifier_t::parameters.calculateSymmeterizedKullbackLeiblerDistancesToStarting << endl;
        }
        if( isModified_coutDistances ) {
          os << "coutDistances = " <<    base_parameters_modifier_t::parameters.coutDistances << endl;
        }
        if( isModified_calculateProfileProfileAlignments ) {
          os << "calculateProfileProfileAlignments = " <<  base_parameters_modifier_t::parameters.calculateProfileProfileAlignments << endl;
        }
        if( isModified_profileProfileIndelOpenCost ) {
          os << "profileProfileIndelOpenCost = " <<  base_parameters_modifier_t::parameters.profileProfileIndelOpenCost << endl;
        }
        if( isModified_profileProfileIndelExtensionCost ) {
          os << "profileProfileIndelExtensionCost = " <<  base_parameters_modifier_t::parameters.profileProfileIndelExtensionCost << endl;
        }
        if( isModified_testTrueProfile ) {
          os << "testTrueProfile = " <<                         base_parameters_modifier_t::parameters.testTrueProfile << endl;
        }
        if( isModified_coutTrueProfile ) {
          os << "coutTrueProfile = " <<                         base_parameters_modifier_t::parameters.coutTrueProfile << endl;
        }
        if( isModified_testStartingProfile ) {
          os << "testStartingProfile = " <<                         base_parameters_modifier_t::parameters.testStartingProfile << endl;
        }
        if( isModified_coutStartingProfile ) {
          os << "coutStartingProfile = " <<                         base_parameters_modifier_t::parameters.coutStartingProfile << endl;
        }
        if( isModified_testUnconditionalProfile ) {
          os << "testUnconditionalProfile = " <<                    base_parameters_modifier_t::parameters.testUnconditionalProfile << endl;
        }
        if( isModified_coutUnconditionalProfile ) {
          os << "coutUnconditionalProfile = " <<                    base_parameters_modifier_t::parameters.coutUnconditionalProfile << endl;
        }
        if( isModified_testUnconditionalWithFixedStartingGlobalsProfile ) {
          os << "testUnconditionalWithFixedStartingGlobalsProfile = " <<                    base_parameters_modifier_t::parameters.testUnconditionalWithFixedStartingGlobalsProfile << endl;
        }
        if( isModified_coutUnconditionalWithFixedStartingGlobalsProfile ) {
          os << "coutUnconditionalWithFixedStartingGlobalsProfile = " <<                    base_parameters_modifier_t::parameters.coutUnconditionalWithFixedStartingGlobalsProfile << endl;
        }
        if( isModified_testUnconditionalWithFixedTrueGlobalsProfile ) {
          os << "testUnconditionalWithFixedTrueGlobalsProfile = " <<                    base_parameters_modifier_t::parameters.testUnconditionalWithFixedTrueGlobalsProfile << endl;
        }
        if( isModified_coutUnconditionalWithFixedTrueGlobalsProfile ) {
          os << "coutUnconditionalWithFixedTrueGlobalsProfile = " <<                    base_parameters_modifier_t::parameters.coutUnconditionalWithFixedTrueGlobalsProfile << endl;
        }
        if( isModified_testConditionalThenUnconditionalProfile ) {
          os << "testConditionalThenUnconditionalProfile = " <<     base_parameters_modifier_t::parameters.testConditionalThenUnconditionalProfile << endl;
        }
        if( isModified_coutConditionalThenUnconditionalProfile ) {
          os << "coutConditionalThenUnconditionalProfile = " <<     base_parameters_modifier_t::parameters.coutConditionalThenUnconditionalProfile << endl;
        }
        if( isModified_testUnconditionalThenConditionalProfile ) {
          os << "testUnconditionalThenConditionalProfile = " <<     base_parameters_modifier_t::parameters.testUnconditionalThenConditionalProfile << endl;
        }
        if( isModified_coutUnconditionalThenConditionalProfile ) {
          os << "coutUnconditionalThenConditionalProfile = " <<     base_parameters_modifier_t::parameters.coutUnconditionalThenConditionalProfile << endl;
        }
        if( isModified_testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile ) {
          os << "testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile = " << base_parameters_modifier_t::parameters.testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile << endl;
        }
        if( isModified_coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile ) {
          os << "coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile = " << base_parameters_modifier_t::parameters.coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile << endl;
        }
        if( isModified_testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile ) {
          os << "testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile = " << base_parameters_modifier_t::parameters.testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile << endl;
        }
        if( isModified_coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile ) {
          os << "coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile = " << base_parameters_modifier_t::parameters.coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile << endl;
        }
        if( isModified_testConditionalGibbsProfile ) {
          os << "testConditionalGibbsProfile = " <<                    base_parameters_modifier_t::parameters.testConditionalGibbsProfile << endl;
        }
        if( isModified_coutConditionalGibbsProfile ) {
          os << "coutConditionalGibbsProfile = " <<                    base_parameters_modifier_t::parameters.coutConditionalGibbsProfile << endl;
        }
        if( isModified_testUnconditionalGibbsProfile ) {
          os << "testUnconditionalGibbsProfile = " <<                    base_parameters_modifier_t::parameters.testUnconditionalGibbsProfile << endl;
        }
        if( isModified_coutUnconditionalGibbsProfile ) {
          os << "coutUnconditionalGibbsProfile = " <<                    base_parameters_modifier_t::parameters.coutUnconditionalGibbsProfile << endl;
        }
        if( isModified_testLengthadjust ) {
          os << "testLengthadjust = " << base_parameters_modifier_t::parameters.testLengthadjust << endl;
        }
        if( isModified_testBaldi ) {
          os << "testBaldi = " << base_parameters_modifier_t::parameters.testBaldi << endl;
        }
        if( isModified_testBaldiSiegel ) {
          os << "testBaldiSiegel = " << base_parameters_modifier_t::parameters.testBaldiSiegel << endl;
        }
        if( isModified_alsoStartWithEvenPositions ) {
          os << "alsoStartWithEvenPositions = " << base_parameters_modifier_t::parameters.alsoStartWithEvenPositions << endl;
        }
      } // writeParametersModifier ( basic_ostream & ) const


  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class ParametersType>
  template <class AnyParameters>
  GALOSH_INLINE_PARAMETERSMODIFIER_APPLY_MODIFICATIONS
  void
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      applyModifications ( AnyParameters & target_parameters )
      {
        base_parameters_modifier_t::applyModifications( target_parameters );

        if( isModified_saveResultsToFile ) {
          target_parameters.saveResultsToFile =
            base_parameters_modifier_t::parameters.saveResultsToFile;
        }
        if( isModified_saveResultsParentDirectory ) {
          target_parameters.saveResultsParentDirectory =
            base_parameters_modifier_t::parameters.saveResultsParentDirectory;
        }
        if( isModified_resultsFilePrefix ) {
          target_parameters.resultsFilePrefix =
            base_parameters_modifier_t::parameters.resultsFilePrefix;
        }
        if( isModified_tabFileSuffix ) {
          target_parameters.tabFileSuffix =
            base_parameters_modifier_t::parameters.tabFileSuffix;
        }
        if( isModified_parametersFileSuffix ) {
          target_parameters.parametersFileSuffix =
            base_parameters_modifier_t::parameters.parametersFileSuffix;
        }
        if( isModified_saveTrueProfileTrees ) {
          target_parameters.saveTrueProfileTrees =
            base_parameters_modifier_t::parameters.saveTrueProfileTrees;
        }
        if( isModified_trueProfileTreeFileSuffix ) {
          target_parameters.trueProfileTreeFileSuffix =
            base_parameters_modifier_t::parameters.trueProfileTreeFileSuffix;
        }
        if( isModified_saveStartingProfiles ) {
          target_parameters.saveStartingProfiles =
            base_parameters_modifier_t::parameters.saveStartingProfiles;
        }
        if( isModified_startingProfileTreeFileSuffix ) {
          target_parameters.startingProfileTreeFileSuffix =
            base_parameters_modifier_t::parameters.startingProfileTreeFileSuffix;
        }
        if( isModified_saveTestProfiles ) {
          target_parameters.saveTestProfiles =
            base_parameters_modifier_t::parameters.saveTestProfiles;
        }
        if( isModified_testProfileTreeFileSuffix ) {
          target_parameters.testProfileTreeFileSuffix =
            base_parameters_modifier_t::parameters.testProfileTreeFileSuffix;
        }
        if( isModified_savePatternSequences ) {
          target_parameters.savePatternSequences =
            base_parameters_modifier_t::parameters.savePatternSequences;
        }
        if( isModified_patternSequencesFileSuffix ) {
          target_parameters.patternSequencesFileSuffix =
            base_parameters_modifier_t::parameters.testsFileSuffix;
        }
        if( isModified_saveTests ) {
          target_parameters.saveTests =
            base_parameters_modifier_t::parameters.saveTests;
        }
        if( isModified_testsFileSuffix ) {
          target_parameters.testsFileSuffix =
            base_parameters_modifier_t::parameters.testsFileSuffix;
        }
        if( isModified_saveTrainingSequences ) {
          target_parameters.saveTrainingSequences =
            base_parameters_modifier_t::parameters.saveTrainingSequences;
        }
        if( isModified_trainingSequencesFileSuffix ) {
          target_parameters.trainingSequencesFileSuffix =
            base_parameters_modifier_t::parameters.trainingSequencesFileSuffix;
        }
        if( isModified_saveTestingSequences ) {
          target_parameters.saveTestingSequences =
            base_parameters_modifier_t::parameters.saveTestingSequences;
        }
        if( isModified_testingSequencesFileSuffix ) {
          target_parameters.testingSequencesFileSuffix =
            base_parameters_modifier_t::parameters.testingSequencesFileSuffix;
        }
        if( isModified_saveTrueTrainingAlignments ) {
          target_parameters.saveTrueTrainingAlignments =
            base_parameters_modifier_t::parameters.saveTrueTrainingAlignments;
        }
        if( isModified_trainingTrueAlignmentsFileSuffix ) {
          target_parameters.trainingTrueAlignmentsFileSuffix =
            base_parameters_modifier_t::parameters.trainingTrueAlignmentsFileSuffix;
        }
        if( isModified_saveTrueTestingAlignments ) {
          target_parameters.saveTrueTestingAlignments =
            base_parameters_modifier_t::parameters.saveTrueTestingAlignments;
        }
        if( isModified_trueTestingAlignmentsFileSuffix ) {
          target_parameters.trueTestingAlignmentsFileSuffix =
            base_parameters_modifier_t::parameters.trueTestingAlignmentsFileSuffix;
        }
        if( isModified_saveFileVersion ) {
          target_parameters.saveFileVersion =
            base_parameters_modifier_t::parameters.saveFileVersion;
        }
        if( isModified_numProfiles ) {
          target_parameters.numProfiles =
            base_parameters_modifier_t::parameters.numProfiles;
        }
        if( isModified_profileLengths ) {
          target_parameters.profileLengths =
            base_parameters_modifier_t::parameters.profileLengths;
        }
        if( isModified_sharedPositionRate ) {
          target_parameters.sharedPositionRate =
            base_parameters_modifier_t::parameters.sharedPositionRate;
        }
        if( isModified_numTrainingSequencesPerProfiles ) {
          target_parameters.numTrainingSequencesPerProfiles =
            base_parameters_modifier_t::parameters.numTrainingSequencesPerProfiles;
        }
        if( isModified_numTestingSequencesPerProfile ) {
          target_parameters.numTestingSequencesPerProfile =
            base_parameters_modifier_t::parameters.numTestingSequencesPerProfile;
        }
        if( isModified_conservationRates ) {
          target_parameters.conservationRates =
            base_parameters_modifier_t::parameters.conservationRates;
        }
        if( isModified_expectedDeletionsCounts ) {
          target_parameters.expectedDeletionsCounts =
            base_parameters_modifier_t::parameters.expectedDeletionsCounts;
        }
        if( isModified_expectedInsertionsCounts ) {
          target_parameters.expectedInsertionsCounts =
            base_parameters_modifier_t::parameters.expectedInsertionsCounts;
        }
        if( isModified_expectedDeletionLengthAsProfileLengthFractions ) {
          target_parameters.expectedDeletionLengthAsProfileLengthFractions =
            base_parameters_modifier_t::parameters.expectedDeletionLengthAsProfileLengthFractions;
        }
        if( isModified_expectedInsertionLengthAsProfileLengthFractions ) {
          target_parameters.expectedInsertionLengthAsProfileLengthFractions =
            base_parameters_modifier_t::parameters.expectedInsertionLengthAsProfileLengthFractions;
        }
        if( isModified_minExpectedDeletionLength ) {
          target_parameters.minExpectedDeletionLength =
            base_parameters_modifier_t::parameters.minExpectedDeletionLength;
        }
        if( isModified_minExpectedInsertionLength ) {
          target_parameters.minExpectedInsertionLength =
            base_parameters_modifier_t::parameters.minExpectedInsertionLength;
        }
        if( isModified_preAlignInsertion ) {
          target_parameters.preAlignInsertion =
            base_parameters_modifier_t::parameters.preAlignInsertion;
        }
        if( isModified_postAlignInsertion ) {
          target_parameters.postAlignInsertion =
            base_parameters_modifier_t::parameters.postAlignInsertion;
        }
        if( isModified_priorStrength ) {
          target_parameters.priorStrength =
            base_parameters_modifier_t::parameters.priorStrength;
        }
        if( isModified_priorStrength_internal_transitions ) {
          target_parameters.priorStrength_internal_transitions =
            base_parameters_modifier_t::parameters.priorStrength_internal_transitions;
        }
        if( isModified_priorMtoM ) {
          target_parameters.priorMtoM =
            base_parameters_modifier_t::parameters.priorMtoM;
        }
        if( isModified_priorMtoI ) {
          target_parameters.priorMtoI =
            base_parameters_modifier_t::parameters.priorMtoI;
        }
        if( isModified_priorMtoD ) {
          target_parameters.priorMtoD =
            base_parameters_modifier_t::parameters.priorMtoD;
        }
        if( isModified_priorItoM ) {
          target_parameters.priorItoM =
            base_parameters_modifier_t::parameters.priorItoM;
        }
        if( isModified_priorItoI ) {
          target_parameters.priorItoI =
            base_parameters_modifier_t::parameters.priorItoI;
        }
        if( isModified_priorDtoM ) {
          target_parameters.priorDtoM =
            base_parameters_modifier_t::parameters.priorDtoM;
        }
        if( isModified_priorDtoD ) {
          target_parameters.priorDtoD =
            base_parameters_modifier_t::parameters.priorDtoD;
        }
        if( isModified_reportGibbsMean ) {
          target_parameters.reportGibbsMean =
            base_parameters_modifier_t::parameters.reportGibbsMean;
        }
        if( isModified_reportGibbsMode ) {
          target_parameters.reportGibbsMode =
            base_parameters_modifier_t::parameters.reportGibbsMode;
        }
        if( isModified_numTrueProfiles ) {
          target_parameters.numTrueProfiles =
            base_parameters_modifier_t::parameters.numTrueProfiles;
        }
        if( isModified_numStartingProfiles ) {
          target_parameters.numStartingProfiles =
            base_parameters_modifier_t::parameters.numStartingProfiles;
        }
        if( isModified_startWithUniformGlobals ) {
          target_parameters.startWithUniformGlobals =
            base_parameters_modifier_t::parameters.startWithUniformGlobals;
        }
        if( isModified_startWithUniformGlobals_scalar ) {
          target_parameters.startWithUniformGlobals_scalar =
            base_parameters_modifier_t::parameters.startWithUniformGlobals_scalar;
        }
        if( isModified_startWithUniformGlobals_maxNtoN ) {
          target_parameters.startWithUniformGlobals_maxNtoN =
            base_parameters_modifier_t::parameters.startWithUniformGlobals_maxNtoN;
        }
        if( isModified_startWithUniformGlobals_maxBtoD ) {
          target_parameters.startWithUniformGlobals_maxBtoD =
            base_parameters_modifier_t::parameters.startWithUniformGlobals_maxBtoD;
        }
        if( isModified_startWithUniformGlobals_maxMtoI ) {
          target_parameters.startWithUniformGlobals_maxMtoI =
            base_parameters_modifier_t::parameters.startWithUniformGlobals_maxMtoI;
        }
        if( isModified_startWithUniformGlobals_maxMtoD ) {
          target_parameters.startWithUniformGlobals_maxMtoD =
            base_parameters_modifier_t::parameters.startWithUniformGlobals_maxMtoD;
        }
        if( isModified_startWithUniformGlobals_maxItoI ) {
          target_parameters.startWithUniformGlobals_maxItoI =
            base_parameters_modifier_t::parameters.startWithUniformGlobals_maxItoI;
        }
        if( isModified_startWithUniformGlobals_maxDtoD ) {
          target_parameters.startWithUniformGlobals_maxDtoD =
            base_parameters_modifier_t::parameters.startWithUniformGlobals_maxDtoD;
        }
        if( isModified_startWithUniformGlobals_maxCtoC ) {
          target_parameters.startWithUniformGlobals_maxCtoC =
            base_parameters_modifier_t::parameters.startWithUniformGlobals_maxCtoC;
        }
        if( isModified_startWithUniformPositions ) {
          target_parameters.startWithUniformPositions =
            base_parameters_modifier_t::parameters.startWithUniformPositions;
        }
        if( isModified_startWithGlobalsDrawnFromPrior ) {
          target_parameters.startWithGlobalsDrawnFromPrior =
            base_parameters_modifier_t::parameters.startWithGlobalsDrawnFromPrior;
        }
        if( isModified_startWithPositionsDrawnFromPrior ) {
          target_parameters.startWithPositionsDrawnFromPrior =
            base_parameters_modifier_t::parameters.startWithPositionsDrawnFromPrior;
        }
        if( isModified_testViterbi ) {
          target_parameters.testViterbi =
            base_parameters_modifier_t::parameters.testViterbi;
        }
        if( isModified_coutViterbi ) {
          target_parameters.coutViterbi =
            base_parameters_modifier_t::parameters.coutViterbi;
        }
        if( isModified_testTruepath ) {
          target_parameters.testTruepath =
            base_parameters_modifier_t::parameters.testTruepath;
        }
        if( isModified_coutTruepath ) {
          target_parameters.coutTruepath =
            base_parameters_modifier_t::parameters.coutTruepath;
        }
        if( isModified_calculateSymmeterizedKullbackLeiblerDistancesToTrue ) {
          target_parameters.calculateSymmeterizedKullbackLeiblerDistancesToTrue =
            base_parameters_modifier_t::parameters.calculateSymmeterizedKullbackLeiblerDistancesToTrue;
        }
        if( isModified_calculateSymmeterizedKullbackLeiblerDistancesToStarting ) {
          target_parameters.calculateSymmeterizedKullbackLeiblerDistancesToStarting =
            base_parameters_modifier_t::parameters.calculateSymmeterizedKullbackLeiblerDistancesToStarting;
        }
        if( isModified_coutDistances ) {
          target_parameters.coutDistances =
            base_parameters_modifier_t::parameters.coutDistances;
        }
        if( isModified_calculateProfileProfileAlignments ) {
          target_parameters.calculateProfileProfileAlignments =
            base_parameters_modifier_t::parameters.calculateProfileProfileAlignments;
        }
        if( isModified_profileProfileIndelOpenCost ) {
          target_parameters.profileProfileIndelOpenCost =
            base_parameters_modifier_t::parameters.profileProfileIndelOpenCost;
        }
        if( isModified_profileProfileIndelExtensionCost ) {
          target_parameters.profileProfileIndelExtensionCost =
            base_parameters_modifier_t::parameters.profileProfileIndelExtensionCost;
        }
        if( isModified_testTrueProfile ) {
          target_parameters.testTrueProfile =
            base_parameters_modifier_t::parameters.testTrueProfile;
        }
        if( isModified_coutTrueProfile ) {
          target_parameters.coutTrueProfile =
            base_parameters_modifier_t::parameters.coutTrueProfile;
        }
        if( isModified_testStartingProfile ) {
          target_parameters.testStartingProfile =
            base_parameters_modifier_t::parameters.testStartingProfile;
        }
        if( isModified_coutStartingProfile ) {
          target_parameters.coutStartingProfile =
            base_parameters_modifier_t::parameters.coutStartingProfile;
        }
        if( isModified_testUnconditionalProfile ) {
          target_parameters.testUnconditionalProfile =
            base_parameters_modifier_t::parameters.testUnconditionalProfile;
        }
        if( isModified_coutUnconditionalProfile ) {
          target_parameters.coutUnconditionalProfile =
            base_parameters_modifier_t::parameters.coutUnconditionalProfile;
        }
        if( isModified_testUnconditionalWithFixedStartingGlobalsProfile ) {
          target_parameters.testUnconditionalWithFixedStartingGlobalsProfile =
            base_parameters_modifier_t::parameters.testUnconditionalWithFixedStartingGlobalsProfile;
        }
        if( isModified_coutUnconditionalWithFixedStartingGlobalsProfile ) {
          target_parameters.coutUnconditionalWithFixedStartingGlobalsProfile =
            base_parameters_modifier_t::parameters.coutUnconditionalWithFixedStartingGlobalsProfile;
        }
        if( isModified_testUnconditionalWithFixedTrueGlobalsProfile ) {
          target_parameters.testUnconditionalWithFixedTrueGlobalsProfile =
            base_parameters_modifier_t::parameters.testUnconditionalWithFixedTrueGlobalsProfile;
        }
        if( isModified_coutUnconditionalWithFixedTrueGlobalsProfile ) {
          target_parameters.coutUnconditionalWithFixedStartingGlobalsProfile =
            base_parameters_modifier_t::parameters.coutUnconditionalWithFixedStartingGlobalsProfile;
        }
        if( isModified_testConditionalThenUnconditionalProfile ) {
          target_parameters.testConditionalThenUnconditionalProfile =
            base_parameters_modifier_t::parameters.testConditionalThenUnconditionalProfile;
        }
        if( isModified_coutConditionalThenUnconditionalProfile ) {
          target_parameters.coutConditionalThenUnconditionalProfile =
            base_parameters_modifier_t::parameters.coutConditionalThenUnconditionalProfile;
        }
        if( isModified_testUnconditionalThenConditionalProfile ) {
          target_parameters.testUnconditionalThenConditionalProfile =
            base_parameters_modifier_t::parameters.testUnconditionalThenConditionalProfile;
        }
        if( isModified_coutUnconditionalThenConditionalProfile ) {
          target_parameters.coutUnconditionalThenConditionalProfile =
            base_parameters_modifier_t::parameters.coutUnconditionalThenConditionalProfile;
        }
        if( isModified_testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile ) {
          target_parameters.testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile =
            base_parameters_modifier_t::parameters.testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile;
        }
        if( isModified_coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile ) {
          target_parameters.coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile =
            base_parameters_modifier_t::parameters.coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile;
        }
        if( isModified_testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile ) {
          target_parameters.testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile =
            base_parameters_modifier_t::parameters.testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile;
        }
        if( isModified_coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile ) {
          target_parameters.coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile =
            base_parameters_modifier_t::parameters.coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile;
        }
        if( isModified_testConditionalGibbsProfile ) {
          target_parameters.testConditionalGibbsProfile =
            base_parameters_modifier_t::parameters.testConditionalGibbsProfile;
        }
        if( isModified_coutConditionalGibbsProfile ) {
          target_parameters.coutConditionalGibbsProfile =
            base_parameters_modifier_t::parameters.coutConditionalGibbsProfile;
        }
        if( isModified_testUnconditionalGibbsProfile ) {
          target_parameters.testUnconditionalGibbsProfile =
            base_parameters_modifier_t::parameters.testUnconditionalGibbsProfile;
        }
        if( isModified_coutUnconditionalGibbsProfile ) {
          target_parameters.coutUnconditionalGibbsProfile =
            base_parameters_modifier_t::parameters.coutUnconditionalGibbsProfile;
        }

        if( isModified_testLengthadjust ) {
          target_parameters.testLengthadjust =
            base_parameters_modifier_t::parameters.testLengthadjust;
        }
        if( isModified_testBaldi ) {
          target_parameters.testBaldi =
            base_parameters_modifier_t::parameters.testBaldi;
        }
        if( isModified_testBaldiSiegel ) {
          target_parameters.testBaldiSiegel =
            base_parameters_modifier_t::parameters.testBaldiSiegel;
        }
        if( isModified_alsoStartWithEvenPositions ) {
          target_parameters.alsoStartWithEvenPositions =
            base_parameters_modifier_t::parameters.alsoStartWithEvenPositions;
        }
      } // applyModifications( Parameters & )

  ////// Class galosh::ProfuseTest ////
  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_INIT
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::
    ProfuseTest (
    ) :
      m_parameters(),
      m_random( static_cast<uint32_t>( std::time( NULL ) ) )
    {
      if( m_parameters.debug >= DEBUG_All ) {
        cout << "[debug] ProfuseTest::<init>()" << endl;
      } // End if DEBUG_All
      // Do nothing else
    } // <init>()

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_INIT
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::
  /**
   * Construct a profuse test object, using the given command-line arguments.
   */  
  ProfuseTest ( const int argc, char ** const argv ) :
      m_parameters(),
      m_random( static_cast<uint32_t>( std::time( NULL ) ) )
  {
    m_parameters.m_profusetest_options.add_options()( "help", "This help message." );
    po::store( po::parse_command_line( argc, argv, m_parameters.m_profusetest_options ), m_parameters.m_options_map );
    // TODO: REMOVE
    cout << "Config file is: " << m_parameters.m_options_map[ "configFile" ].as<string>() << endl;
    ifstream configFile( m_parameters.m_options_map[ "configFile" ].as<string>().c_str() );
    if( configFile ) {
      try {
        po::store( parse_config_file( configFile, m_parameters.m_profusetest_options ), m_parameters.m_options_map );
      } catch( const std::exception& e ) {
        std::cerr << std::endl << "ERROR PARSING ProfuseTest CONFIG FILE: " << e.what() << std::endl;
        exit( 1 );
      }
      configFile.close();
    }
    po::notify( m_parameters.m_options_map );
    if( m_parameters.m_options_map.count("help") ) {
      cout << m_parameters.m_profusetest_options << endl;
      exit( 0 );
    }
    if( m_parameters.m_options_map.count("seed") ) {
      if( m_parameters.seed != 0 ) {
        m_random.setSeed( m_parameters.seed );
      }
    }
    // TODO: REMOVE
    //cout << "Hi, ok, done constructing" << endl;
  } // <init>( const int argc, char ** const argv )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_INIT
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::
    /**
     * Construct a profuse test object, using the provided seed.
     */  
    ProfuseTest (
      uint32_t const seed
    ) :
      m_parameters(),
      m_random( seed )
    {
      if( m_parameters.debug >= DEBUG_All ) {
        cout << "[debug] ProfuseTest::<init>( uint32_t )" << endl;
      } // End if DEBUG_All
      // Do nothing else
    } // <init>( uint32_t )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_PROFUSETEST_START
  void
  ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::
    start ()
    {
      bool be_verbose = false;//( m_parameters.verbosity >= VERBOSITY_Meta );

      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DirichletMixtureMatchEmissionPrior<float> matchEmissionPrior;

      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DirichletMixtureGlobalPrior<float> globalPrior;

      // Set up the priors
      if( m_parameters.usePriors || m_parameters.startWithPositionsDrawnFromPrior ) {
        matchEmissionPrior.reinitializeToEven( m_parameters.priorStrength );
        //matchEmissionPrior.reinitializeToEven( ( .5f * m_parameters.priorStrength ) );
      } // End if m_parameters.usePriors || m_parameters.startWithPositionsDrawnFromPrior
      if( m_parameters.usePriors || m_parameters.startWithGlobalsDrawnFromPrior ) {
        globalPrior.reinitializeToEven( m_parameters.priorStrength );
        // NOTE: We will do additional set-up of the global prior for each
        // change in the profile length, since we need to adjust for profile
        // length.

        // TODO: Make these parameters
        double priorStrength_flanking_self_transitions = 100;
        double priorStrength_flanking_other_transitions = 100;
        globalPrior[ 0 ][ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] = ( priorStrength_flanking_self_transitions * .01 );
        globalPrior[ 0 ][ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ] = ( priorStrength_flanking_other_transitions * .99 );
#ifdef USE_DEL_IN_DEL_OUT
      // TODO: Create training_parameters_template.priorWtoW
      globalPrior[ 0 ][ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toDeletionOut ] = ( m_parameters.priorStrength_flanking_other_transitions * .5 );
      // TODO: Create training_parameters_template.priorWtoE
      globalPrior[ 0 ][ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toEnd ] = ( priorStrength_flanking_other_transitions * .5 );
      globalPrior[ 0 ][ Transition::fromBegin ][ TransitionFromBegin::toDeletion ] = ( priorStrength_flanking_other_transitions * .01 );
      // TODO: Create training_parameters_template.priorBtoZ
      globalPrior[ 0 ][ Transition::fromBegin ][ TransitionFromBegin::toDeletionIn ] = ( priorStrength_flanking_other_transitions * .50 );
      // TODO: Create training_parameters_template.priorZtoZ
      globalPrior[ 0 ][ Transition::fromDeletionIn ][ TransitionFromDeletionIn::toDeletionIn ] = ( priorStrength_flanking_other_transitions * .5 );
      // TODO: Create training_parameters_template.priorZtoM
      globalPrior[ 0 ][ Transition::fromDeletionIn ][ TransitionFromDeletionIn::toMatch ] = ( priorStrength_flanking_other_transitions * .5 );
#else
        globalPrior[ 0 ][ Transition::fromBegin ][ TransitionFromBegin::toDeletion ] = ( priorStrength_flanking_other_transitions * .01 );
        globalPrior[ 0 ][ Transition::fromBegin ][ TransitionFromBegin::toMatch ] = ( priorStrength_flanking_other_transitions * .99 );
#endif // USE_DEL_IN_DEL_OUT .. else ..
        globalPrior[ 0 ][ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] = ( priorStrength_flanking_self_transitions * .01 );
        globalPrior[ 0 ][ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] = ( priorStrength_flanking_other_transitions * .99 );
      } // End if usePriors

      vector<Test> tests( LAST_TEST_ID + 1 );
      //for( test_id = 0; test_id <= LAST_TEST_ID; test_id++ ) {
      //  tests[ test_id ].name = testNames[ test_id ];
      //  tests[ test_id ].isRun = runTest[ test_id ];
      //  tests[ test_id ].isCout = coutTest[ test_id ];
      //  tests[ test_id ].parametersModifier = testParametersModifier[ test_id ];
      //  if( ( test_id == TEST_ID_true ) ||
      //      ( test_id == TEST_ID_starting ) ) {
      //    tests[ test_id ].startingGlobalsTest = NULL; // ignored
      //    tests[ test_id ].startingPositionsTest = NULL; // ignored
      //  } else { // test_id is not true or starting
      //    if( testStartingGlobals[ test_id ] == TEST_ID_starting ) {
      //      tests[ test_id ].startingGlobalsTest = NULL; // use starting profile
      //    } else if( testStartingGlobals[ test_id ] == TEST_ID_true ) {
      //      tests[ test_id ].startingGlobalsTest =
      //        &tests[ test_id ]; // use true profile
      //    } else {
      //      tests[ test_id ].startingGlobalsTest =
      //        &tests[ testStartingGlobals[ test_id ] ];
      //    }
      //
      //    if( testStartingPositions[ test_id ] == TEST_ID_starting ) {
      //      tests[ test_id ].startingPositionsTest = NULL; // use starting profile
      //    } else if( testStartingPositions[ test_id ] == TEST_ID_true ) {
      //      tests[ test_id ].startingPositionsTest =
      //        &tests[ test_id ]; // use true profile
      //    } else {
      //      tests[ test_id ].startingPositionsTest =
      //        &tests[ testStartingPositions[ test_id ] ];
      //    }
      //  } // End if test_id is true or starting .. else ..
      //  if( test_id == TEST_ID_conditional ) {
      //    tests[ test_id ].coutLeftBrace = "(";
      //    tests[ test_id ].coutRightBrace = ")";
      //  } else if( test_id == TEST_ID_unconditional ) {
      //    tests[ test_id ].coutLeftBrace = "<";
      //    tests[ test_id ].coutRightBrace = ">";
      //  }
      //} // End foreach test_id, set up tests[ test_id ]


      // Set up the convenient arrays...
      /**
       * For convenience we number the tests from 0 to LAST_TEST_ID.
       *
       * If runTest[ test_num ] is true, we are running the corresponding test;
       * the user controls this via the appropriate parameter (in the
       * Parameters class).  Note that some tests will be run that the user did
       * not specify in parameters: any test that is the starting value for a
       * test that is to be run will also be run.  So for example, if
       * parameters.testUnconditionalThenConditionalProfile is true, then
       * runTest[ TEST_ID_unconditional ] will be set to true even if
       * parameters.testUnconditionalProfile is false.
       *
       * @see testStartingGlobals
       * @see testStartingPositions
       */
      //bool runTest[ LAST_TEST_ID + 1 ];

      /**
       * For convenience we number the tests from 0 to LAST_TEST_ID.
       *
       * Test test_num will use the test profile with id testStartingGlobals[
       * test_num ] as the starting value for the globals, if runTest[ test_num
       * ] is true.  For the true profile (TEST_ID_true) and starting profile
       * (TEST_ID_starting), this value is ignored.
       *
       * NOTE: The order of the tests matters.  No test should have a test_id
       * lower than its testStartingGlobals value, and only the true profile
       * and the starting profile should have test_ids equal to their
       * testStartingGlobals values.
       * 
       * @see testStartingPositions
       */
      //uint32_t testStartingGlobals[ LAST_TEST_ID + 1 ];

      /**
       * For convenience we number the tests from 0 to LAST_TEST_ID.
       *
       * Test test_num will use the test profile with id testStartingPositions[
       * test_num ] as the starting value for the positions, if runTest[
       * test_num ] is true.  For the true profile (TEST_ID_true) and starting
       * profile (TEST_ID_starting), this value is ignored.
       *
       * NOTE: The order of the tests matters.  No test should have a test_id
       * lower than its testStartingPositions value, and only the true profile
       * and the starting profile should have test_ids equal to their
       * testStartingPositions values.
       * 
       * @see testStartingGlobals.
       */
      //uint32_t tests[ LAST_TEST_ID + 1 ].startingPositionsTest;

      /**
       * For convenience we number the tests from 0 to LAST_TEST_ID.
       *
       * We write a subset of results to STDOUT (cout) during training.  These
       * bools indicate which tests we write to cout.  Note that all tests are
       * written to the file, regardless of the value of coutTest[ test_num ].
       *
       * If runTest[ test_num ] and coutTest[ test_num ] are both true, 
       * then put the test results of test test_num on the STDOUT stream too.
       */
      //bool coutTest[ LAST_TEST_ID + 1 ];

      /**
       * For convenience we number the tests from 0 to LAST_TEST_ID.
       *
       * The parameters modifier for each test.
       */
      //typename ProfileTrainer<RootType, ScoreType, MatrixValueType>::ParametersModifier tests[ LAST_TEST_ID + 1 ].parametersModifier;

      //string testNames[ LAST_TEST_ID + 1 ] =
      //{
      //  "true",
      //  "starting",
      //  "conditional",
      //  "unconditional",
      //  "unconditional_with_fixed_starting_globals",
      //  "unconditional_with_fixed_true_globals",
      //  "conditional_then_unconditional",
      //  "unconditional_then_conditional",
      //  "unconditional_with_fixed_starting_globals_then_with_fixed_positions",
      //  "unconditional_with_fixed_true_globals_then_with_fixed_positions"
      //}; // test_names

      // Set up the tests.
      if( m_parameters.testTrueProfile ) {
        tests[ TEST_ID_true ].name = "true";
        tests[ TEST_ID_true ].isRun = true;
        if( m_parameters.coutTrueProfile ) {
          tests[ TEST_ID_true ].isCout = true;
        } else {
          tests[ TEST_ID_true ].isCout = false;
        }
      } else {
        tests[ TEST_ID_true ].isRun = false;
      }
      tests[ TEST_ID_true ].isGibbs = false;
      tests[ TEST_ID_true ].startingGlobalsTest = NULL; // ignored
      tests[ TEST_ID_true ].startingPositionsTest = NULL; // ignored
      // No modifications, since we don't train the true profile.
      // tests[ TEST_ID_true ].parametersModifier

      if( m_parameters.testStartingProfile ) {
        tests[ TEST_ID_starting ].name = "starting";
        tests[ TEST_ID_starting ].isRun = true;
        if( m_parameters.coutStartingProfile ) {
          tests[ TEST_ID_starting ].isCout = true;
        } else {
          tests[ TEST_ID_starting ].isCout = false;
        }
      } else {
        tests[ TEST_ID_starting ].isRun = false;
      }
      tests[ TEST_ID_starting ].isGibbs = false;
      tests[ TEST_ID_starting ].startingGlobalsTest =
        NULL; // ignored, but see below for startWithUniformGlobals
      tests[ TEST_ID_starting ].startingPositionsTest = NULL; // ignored
      // No modifications, since we don't train the starting profile.
      // tests[ TEST_ID_starting ].parametersModifier

      // Always test conditional BW
      tests[ TEST_ID_conditional ].name = "conditional";
      tests[ TEST_ID_conditional ].isRun = true;
      // Always write conditional BW results to STDOUT
      tests[ TEST_ID_conditional ].isCout = true;
      tests[ TEST_ID_conditional ].isGibbs = false;
      tests[ TEST_ID_conditional ].coutLeftBrace = "(";
      tests[ TEST_ID_conditional ].coutRightBrace = ")";
      // Start with the starting profile globals and positions
      tests[ TEST_ID_conditional ].startingGlobalsTest = NULL;
      tests[ TEST_ID_conditional ].startingPositionsTest = NULL;
      // Modifiers for conditional test
      // inherit these
      //tests[ TEST_ID_conditional ].parametersModifier.parameters.trainProfilePositions = true;
      //tests[ TEST_ID_conditional ].parametersModifier.isModified_trainProfilePositions = true;
      //if( m_parameters.trainProfileGlobals ) {
      //  tests[ TEST_ID_conditional ].parametersModifier.parameters.trainProfileGlobals = true;
      //  tests[ TEST_ID_conditional ].parametersModifier.isModified_trainProfileGlobals = true;
      //}
      tests[ TEST_ID_conditional ].parametersModifier.parameters.useUnconditionalBaumWelch = false;
      tests[ TEST_ID_conditional ].parametersModifier.isModified_useUnconditionalBaumWelch = true;

      // TODO: REMOVE
      //tests[ TEST_ID_conditional ].parametersModifier.parameters.verbosity = VERBOSITY_All;
      //tests[ TEST_ID_conditional ].parametersModifier.isModified_verbosity = true;

      if( m_parameters.testUnconditionalProfile ) {
        tests[ TEST_ID_unconditional ].name = "unconditional";
        tests[ TEST_ID_unconditional ].isRun = true;
      } else {
        tests[ TEST_ID_unconditional ].isRun = false;
      }
      if( m_parameters.coutUnconditionalProfile ) {
        tests[ TEST_ID_unconditional ].isCout = true;
      } else {
        tests[ TEST_ID_unconditional ].isCout = false;
      }
      tests[ TEST_ID_unconditional ].isGibbs = false;
      tests[ TEST_ID_unconditional ].coutLeftBrace = "<";
      tests[ TEST_ID_unconditional ].coutRightBrace = ">";
      // Start with the starting profile globals and positions
      tests[ TEST_ID_unconditional ].startingGlobalsTest = NULL;
      tests[ TEST_ID_unconditional ].startingPositionsTest = NULL;
      // Modifiers for unconditional test
      // inherit these
      //tests[ TEST_ID_unconditional ].parametersModifier.parameters.trainProfilePositions = true;
      //tests[ TEST_ID_unconditional ].parametersModifier.isModified_trainProfilePositions = true;
      //if( m_parameters.trainProfileGlobals ) {
      //  tests[ TEST_ID_unconditional ].parametersModifier.parameters.trainProfileGlobals = true;
      //  tests[ TEST_ID_unconditional ].parametersModifier.isModified_trainProfileGlobals = true;
      //}
      tests[ TEST_ID_unconditional ].parametersModifier.parameters.useUnconditionalBaumWelch = true;
      tests[ TEST_ID_unconditional ].parametersModifier.isModified_useUnconditionalBaumWelch = true;

      // TODO: REMOVE
      //tests[ TEST_ID_unconditional ].parametersModifier.parameters.verbosity = VERBOSITY_All;
      //tests[ TEST_ID_unconditional ].parametersModifier.isModified_verbosity = true;
      //tests[ TEST_ID_unconditional ].parametersModifier.parameters.trainProfileGlobals = false;
      //tests[ TEST_ID_unconditional ].parametersModifier.isModified_trainProfileGlobals = true;

      if( m_parameters.testUnconditionalWithFixedStartingGlobalsProfile ) {
        tests[ TEST_ID_unconditional_with_fixed_starting_globals ].name = "unconditional_with_fixed_starting_globals";
        tests[ TEST_ID_unconditional_with_fixed_starting_globals ].isRun = true;
      } else {
        tests[ TEST_ID_unconditional_with_fixed_starting_globals ].isRun = false;
      }
      if( m_parameters.coutUnconditionalWithFixedStartingGlobalsProfile ) {
        tests[ TEST_ID_unconditional_with_fixed_starting_globals ].isCout = true;
      } else {
        tests[ TEST_ID_unconditional_with_fixed_starting_globals ].isCout = false;
      }
      tests[ TEST_ID_unconditional_with_fixed_starting_globals ].isGibbs = false;
      // Start with the starting profile globals and positions
      tests[ TEST_ID_unconditional_with_fixed_starting_globals ].startingGlobalsTest =
        NULL;
      tests[ TEST_ID_unconditional_with_fixed_starting_globals ].startingPositionsTest =
        NULL;
      // Modifiers for unconditional_with_fixed_starting_globals test
      tests[ TEST_ID_unconditional_with_fixed_starting_globals ].parametersModifier.parameters.trainProfilePositions = true;
      tests[ TEST_ID_unconditional_with_fixed_starting_globals ].parametersModifier.isModified_trainProfilePositions = true;
      tests[ TEST_ID_unconditional_with_fixed_starting_globals ].parametersModifier.parameters.trainProfileGlobals = false;
      tests[ TEST_ID_unconditional_with_fixed_starting_globals ].parametersModifier.isModified_trainProfileGlobals = true;
      tests[ TEST_ID_unconditional_with_fixed_starting_globals ].parametersModifier.parameters.useUnconditionalBaumWelch = true;
      tests[ TEST_ID_unconditional_with_fixed_starting_globals ].parametersModifier.isModified_useUnconditionalBaumWelch = true;

      if( m_parameters.testUnconditionalWithFixedTrueGlobalsProfile ) {
        tests[ TEST_ID_unconditional_with_fixed_true_globals ].name = "unconditional_with_fixed_true_globals";
        tests[ TEST_ID_unconditional_with_fixed_true_globals ].isRun = true;
      } else {
        tests[ TEST_ID_unconditional_with_fixed_true_globals ].isRun = false;
      }
      if( m_parameters.coutUnconditionalWithFixedTrueGlobalsProfile ) {
        tests[ TEST_ID_unconditional_with_fixed_true_globals ].isCout = true;
      } else {
        tests[ TEST_ID_unconditional_with_fixed_true_globals ].isCout = false;
      }
      tests[ TEST_ID_unconditional_with_fixed_true_globals ].isGibbs = false;
      // Start with the true profile globals
      tests[ TEST_ID_unconditional_with_fixed_true_globals ].startingGlobalsTest =
        &tests[ TEST_ID_unconditional_with_fixed_true_globals ]; // true
      // Start with the starting profile positions
      tests[ TEST_ID_unconditional_with_fixed_true_globals ].startingPositionsTest =
        NULL;
      // Modifiers for unconditional_with_fixed_true_globals test
      tests[ TEST_ID_unconditional_with_fixed_true_globals ].parametersModifier.parameters.trainProfilePositions = true;
      tests[ TEST_ID_unconditional_with_fixed_true_globals ].parametersModifier.isModified_trainProfilePositions = true;
      tests[ TEST_ID_unconditional_with_fixed_true_globals ].parametersModifier.parameters.trainProfileGlobals = false;
      tests[ TEST_ID_unconditional_with_fixed_true_globals ].parametersModifier.isModified_trainProfileGlobals = true;
      tests[ TEST_ID_unconditional_with_fixed_true_globals ].parametersModifier.parameters.useUnconditionalBaumWelch = true;
      tests[ TEST_ID_unconditional_with_fixed_true_globals ].parametersModifier.isModified_useUnconditionalBaumWelch = true;

      if( m_parameters.testConditionalThenUnconditionalProfile ) {
        tests[ TEST_ID_conditional_then_unconditional ].name = "conditional_then_unconditional";
        tests[ TEST_ID_conditional_then_unconditional ].isRun = true;
      } else {
        tests[ TEST_ID_conditional_then_unconditional ].isRun = false;
      }
      if( m_parameters.coutConditionalThenUnconditionalProfile ) {
        tests[ TEST_ID_conditional_then_unconditional ].isCout = true;
      } else {
        tests[ TEST_ID_conditional_then_unconditional ].isCout = false;
      }
      tests[ TEST_ID_conditional_then_unconditional ].isGibbs = false;
      // Start with the conditional profile globals and positions
      tests[ TEST_ID_conditional_then_unconditional ].startingGlobalsTest =
        &tests[ TEST_ID_conditional ];
      tests[ TEST_ID_conditional_then_unconditional ].startingPositionsTest =
        &tests[ TEST_ID_conditional ];
      // Modifiers for conditional_then_unconditional test
      tests[ TEST_ID_conditional_then_unconditional ].parametersModifier.parameters.trainProfilePositions = true;
      tests[ TEST_ID_conditional_then_unconditional ].parametersModifier.isModified_trainProfilePositions = true;
      if( m_parameters.trainProfileGlobals ) {
        tests[ TEST_ID_conditional_then_unconditional ].parametersModifier.parameters.trainProfileGlobals = true;
        tests[ TEST_ID_conditional_then_unconditional ].parametersModifier.isModified_trainProfileGlobals = true;
      }
      tests[ TEST_ID_conditional_then_unconditional ].parametersModifier.parameters.useUnconditionalBaumWelch = true;
      tests[ TEST_ID_conditional_then_unconditional ].parametersModifier.isModified_useUnconditionalBaumWelch = true;

      if( m_parameters.testUnconditionalThenConditionalProfile ) {
        tests[ TEST_ID_unconditional_then_conditional ].name = "unconditional_then_conditional";
        tests[ TEST_ID_unconditional_then_conditional ].isRun = true;
      } else {
        tests[ TEST_ID_unconditional_then_conditional ].isRun = false;
      }
      if( m_parameters.coutUnconditionalThenConditionalProfile ) {
        tests[ TEST_ID_unconditional_then_conditional ].isCout = true;
      } else {
        tests[ TEST_ID_unconditional_then_conditional ].isCout = false;
      }
      tests[ TEST_ID_unconditional_then_conditional ].isGibbs = false;
      // Start with the unconditional profile globals and positions
      tests[ TEST_ID_unconditional_then_conditional ].startingGlobalsTest =
        &tests[ TEST_ID_unconditional ];
      tests[ TEST_ID_unconditional_then_conditional ].startingPositionsTest =
        &tests[ TEST_ID_unconditional ];
      // Modifiers for unconditional_then_conditional test
      tests[ TEST_ID_unconditional_then_conditional ].parametersModifier.parameters.trainProfilePositions = true;
      tests[ TEST_ID_unconditional_then_conditional ].parametersModifier.isModified_trainProfilePositions = true;
      if( m_parameters.trainProfileGlobals ) {
        tests[ TEST_ID_unconditional_then_conditional ].parametersModifier.parameters.trainProfileGlobals = true;
        tests[ TEST_ID_unconditional_then_conditional ].parametersModifier.isModified_trainProfileGlobals = true;
      }
      tests[ TEST_ID_unconditional_then_conditional ].parametersModifier.parameters.useUnconditionalBaumWelch = false;
      tests[ TEST_ID_unconditional_then_conditional ].parametersModifier.isModified_useUnconditionalBaumWelch = true;

      if( m_parameters.testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile ) {
        tests[ TEST_ID_unconditional_with_fixed_starting_globals_then_with_fixed_positions ].name = "unconditional_with_fixed_starting_globals_then_with_fixed_positions";
        tests[ TEST_ID_unconditional_with_fixed_starting_globals_then_with_fixed_positions ].isRun = true;
      } else {
        tests[ TEST_ID_unconditional_with_fixed_starting_globals_then_with_fixed_positions ].isRun = false;
      }
      if( m_parameters.coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile ) {
        tests[ TEST_ID_unconditional_with_fixed_starting_globals_then_with_fixed_positions ].isCout = true;
      } else {
        tests[ TEST_ID_unconditional_with_fixed_starting_globals_then_with_fixed_positions ].isCout = false;
      }
      tests[ TEST_ID_unconditional_with_fixed_starting_globals_then_with_fixed_positions ].isGibbs = false;
      // Start with the unconditional_with_fixed_starting_globals profile
      // globals and positions
      tests[ TEST_ID_unconditional_with_fixed_starting_globals_then_with_fixed_positions ].startingGlobalsTest =
        &tests[ TEST_ID_unconditional_with_fixed_starting_globals ];
      tests[ TEST_ID_unconditional_with_fixed_starting_globals_then_with_fixed_positions ].startingPositionsTest =
        &tests[ TEST_ID_unconditional_with_fixed_starting_globals ];
      // Modifiers for unconditional_with_fixed_starting_globals_then_with_fixed_positions test
      tests[ TEST_ID_unconditional_with_fixed_starting_globals_then_with_fixed_positions ].parametersModifier.parameters.trainProfilePositions = false;
      tests[ TEST_ID_unconditional_with_fixed_starting_globals_then_with_fixed_positions ].parametersModifier.isModified_trainProfilePositions = true;
      if( m_parameters.trainProfileGlobals ) {
        tests[ TEST_ID_unconditional_with_fixed_starting_globals_then_with_fixed_positions ].parametersModifier.parameters.trainProfileGlobals = true;
        tests[ TEST_ID_unconditional_with_fixed_starting_globals_then_with_fixed_positions ].parametersModifier.isModified_trainProfileGlobals = true;
      }
      tests[ TEST_ID_unconditional_with_fixed_starting_globals_then_with_fixed_positions ].parametersModifier.parameters.useUnconditionalBaumWelch = true;
      tests[ TEST_ID_unconditional_with_fixed_starting_globals_then_with_fixed_positions ].parametersModifier.isModified_useUnconditionalBaumWelch = true;

      if( m_parameters.testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile ) {
        tests[ TEST_ID_unconditional_with_fixed_true_globals_then_with_fixed_positions ].name = "unconditional_with_fixed_true_globals_then_with_fixed_positions";
        tests[ TEST_ID_unconditional_with_fixed_true_globals_then_with_fixed_positions ].isRun = true;
      } else {
        tests[ TEST_ID_unconditional_with_fixed_true_globals_then_with_fixed_positions ].isRun = false;
      }
      if( m_parameters.coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile ) {
        tests[ TEST_ID_unconditional_with_fixed_true_globals_then_with_fixed_positions ].isCout = true;
      } else {
        tests[ TEST_ID_unconditional_with_fixed_true_globals_then_with_fixed_positions ].isCout = false;
      }
      tests[ TEST_ID_unconditional_with_fixed_true_globals_then_with_fixed_positions ].isGibbs = false;
      // Start with the unconditional_with_fixed_true_globals profile
      // globals and positions
      tests[ TEST_ID_unconditional_with_fixed_true_globals_then_with_fixed_positions ].startingGlobalsTest =
        &tests[ TEST_ID_unconditional_with_fixed_true_globals ];
      tests[ TEST_ID_unconditional_with_fixed_true_globals_then_with_fixed_positions ].startingPositionsTest =
        &tests[ TEST_ID_unconditional_with_fixed_true_globals ];
      // Modifiers for unconditional_with_fixed_true_globals_then_with_fixed_positions test
      tests[ TEST_ID_unconditional_with_fixed_true_globals_then_with_fixed_positions ].parametersModifier.parameters.trainProfilePositions = false;
      tests[ TEST_ID_unconditional_with_fixed_true_globals_then_with_fixed_positions ].parametersModifier.isModified_trainProfilePositions = true;
      if( m_parameters.trainProfileGlobals ) {
        tests[ TEST_ID_unconditional_with_fixed_true_globals_then_with_fixed_positions ].parametersModifier.parameters.trainProfileGlobals = true;
        tests[ TEST_ID_unconditional_with_fixed_true_globals_then_with_fixed_positions ].parametersModifier.isModified_trainProfileGlobals = true;
      }
      tests[ TEST_ID_unconditional_with_fixed_true_globals_then_with_fixed_positions ].parametersModifier.parameters.useUnconditionalBaumWelch = true;
      tests[ TEST_ID_unconditional_with_fixed_true_globals_then_with_fixed_positions ].parametersModifier.isModified_useUnconditionalBaumWelch = true;

      if( m_parameters.testConditionalGibbsProfile ) {
        tests[ TEST_ID_gibbs_conditional ].name = "cGibbs";
        tests[ TEST_ID_gibbs_conditional ].isRun = true;
      } else {
        tests[ TEST_ID_gibbs_conditional ].isRun = false;
      }
      if( m_parameters.coutConditionalGibbsProfile ) {
        tests[ TEST_ID_gibbs_conditional ].isCout = true;
      } else {
        tests[ TEST_ID_gibbs_conditional ].isCout = false;
      }
      tests[ TEST_ID_gibbs_conditional ].isGibbs = true;
      tests[ TEST_ID_gibbs_conditional ].coutLeftBrace = "`";
      tests[ TEST_ID_gibbs_conditional ].coutRightBrace = "'";
      // Start with the starting profile globals and positions
      tests[ TEST_ID_gibbs_conditional ].startingGlobalsTest =
        NULL;
      tests[ TEST_ID_gibbs_conditional ].startingPositionsTest =
        NULL;
      // Modifiers for conditional gibbs test
      // inherit these
      //tests[ TEST_ID_gibbs_conditional ].parametersModifier.parameters.sampleProfilePositions = true;
      //tests[ TEST_ID_gibbs_conditional ].parametersModifier.isModified_sampleProfilePositions = true;
      //if( m_parameters.trainProfileGlobals ) {
      //  tests[ TEST_ID_gibbs_conditional ].parametersModifier.parameters.sampleProfileGlobals = true;
      //  tests[ TEST_ID_gibbs_conditional ].parametersModifier.isModified_sampleProfileGlobals = true;
      //} else {
      //  tests[ TEST_ID_gibbs_conditional ].parametersModifier.parameters.sampleProfileGlobals = false;
      //  tests[ TEST_ID_gibbs_conditional ].parametersModifier.isModified_sampleProfileGlobals = false;
      //}
      tests[ TEST_ID_gibbs_conditional ].parametersModifier.parameters.useUnconditionalGibbs = false;
      tests[ TEST_ID_gibbs_conditional ].parametersModifier.isModified_useUnconditionalGibbs = true;
      if( m_parameters.reportGibbsMode ) {
        tests[ TEST_ID_gibbs_conditional ].parametersModifier.parameters.saveGibbsMode = true;
        tests[ TEST_ID_gibbs_conditional ].parametersModifier.isModified_saveGibbsMode = true;
      } // End if reportGibbsMode

      if( m_parameters.testUnconditionalGibbsProfile ) {
        tests[ TEST_ID_gibbs_unconditional ].name = "uGibbs";
        tests[ TEST_ID_gibbs_unconditional ].isRun = true;
      } else {
        tests[ TEST_ID_gibbs_unconditional ].isRun = false;
      }
      if( m_parameters.coutUnconditionalGibbsProfile ) {
        tests[ TEST_ID_gibbs_unconditional ].isCout = true;
      } else {
        tests[ TEST_ID_gibbs_unconditional ].isCout = false;
      }
      tests[ TEST_ID_gibbs_unconditional ].isGibbs = true;
      tests[ TEST_ID_gibbs_unconditional ].coutLeftBrace = "*";
      tests[ TEST_ID_gibbs_unconditional ].coutRightBrace = "*";
      // Start with the starting profile globals and positions
      tests[ TEST_ID_gibbs_unconditional ].startingGlobalsTest =
        // TODO: Put back: NULL;
        &tests[ TEST_ID_gibbs_unconditional ]; // true
      tests[ TEST_ID_gibbs_unconditional ].startingPositionsTest = NULL;
      // Modifiers for unconditional gibbs test
      // inherit these, don't modify them.
      //tests[ TEST_ID_gibbs_unconditional ].parametersModifier.parameters.sampleProfilePositions = true;
      //tests[ TEST_ID_gibbs_unconditional ].parametersModifier.isModified_sampleProfilePositions = true;
      //if( m_parameters.trainProfileGlobals ) {
      //  tests[ TEST_ID_gibbs_unconditional ].parametersModifier.parameters.sampleProfileGlobals = true;
      //  tests[ TEST_ID_gibbs_unconditional ].parametersModifier.isModified_sampleProfileGlobals = true;
      //} else {
      //  tests[ TEST_ID_gibbs_unconditional ].parametersModifier.parameters.sampleProfileGlobals = false;
      //  tests[ TEST_ID_gibbs_unconditional ].parametersModifier.isModified_sampleProfileGlobals = false;
      //}
      tests[ TEST_ID_gibbs_unconditional ].parametersModifier.parameters.useUnconditionalGibbs = true;
      tests[ TEST_ID_gibbs_unconditional ].parametersModifier.isModified_useUnconditionalGibbs = true;
      if( m_parameters.reportGibbsMode ) {
        tests[ TEST_ID_gibbs_unconditional ].parametersModifier.parameters.saveGibbsMode = true;
        tests[ TEST_ID_gibbs_unconditional ].parametersModifier.isModified_saveGibbsMode = true;
      } // End if reportGibbsMode

      if( m_parameters.testLengthadjust ) {
        tests[ TEST_ID_lengthadjust_conditional ].name = "conditionalLA";
        tests[ TEST_ID_lengthadjust_conditional ].isRun = true;
        // Always write conditional BW results to STDOUT
        tests[ TEST_ID_lengthadjust_conditional ].isCout = true;
        tests[ TEST_ID_lengthadjust_conditional ].isGibbs = false;
        tests[ TEST_ID_lengthadjust_conditional ].coutLeftBrace = "(*";
        tests[ TEST_ID_lengthadjust_conditional ].coutRightBrace = "*)";
        // Start with the starting profile globals and positions
        tests[ TEST_ID_lengthadjust_conditional ].startingGlobalsTest = NULL;
        tests[ TEST_ID_lengthadjust_conditional ].startingPositionsTest = NULL;
        // Modifiers for lengthadjust_conditional test
        tests[ TEST_ID_lengthadjust_conditional ].parametersModifier.parameters.proposeProfileLengthChanges = true;
        tests[ TEST_ID_lengthadjust_conditional ].parametersModifier.isModified_proposeProfileLengthChanges = true;
        tests[ TEST_ID_lengthadjust_conditional ].parametersModifier.parameters.useUnconditionalBaumWelch = false;
        tests[ TEST_ID_lengthadjust_conditional ].parametersModifier.isModified_useUnconditionalBaumWelch = true;
        
        // TODO: REMOVE
        //tests[ TEST_ID_lengthadjust_conditional ].parametersModifier.parameters.verbosity = VERBOSITY_All;
        //tests[ TEST_ID_lengthadjust_conditional ].parametersModifier.isModified_verbosity = true;
       } // End if testLengthadjust
        
      if( m_parameters.testLengthadjust ) {
        if( m_parameters.testUnconditionalProfile ) {
          tests[ TEST_ID_lengthadjust_unconditional ].name = "unconditionalLA";
          tests[ TEST_ID_lengthadjust_unconditional ].isRun = true;
        } else {
          tests[ TEST_ID_lengthadjust_unconditional ].isRun = false;
        }
        if( m_parameters.coutUnconditionalProfile ) {
          tests[ TEST_ID_lengthadjust_unconditional ].isCout = true;
        } else {
          tests[ TEST_ID_lengthadjust_unconditional ].isCout = false;
        }
        tests[ TEST_ID_lengthadjust_unconditional ].isGibbs = false;
        tests[ TEST_ID_lengthadjust_unconditional ].coutLeftBrace = "<*";
        tests[ TEST_ID_lengthadjust_unconditional ].coutRightBrace = "*>";
        // Start with the starting profile globals and positions
        tests[ TEST_ID_lengthadjust_unconditional ].startingGlobalsTest = NULL;
        tests[ TEST_ID_lengthadjust_unconditional ].startingPositionsTest = NULL;
        // Modifiers for lengthadjust_unconditional test
        tests[ TEST_ID_lengthadjust_unconditional ].parametersModifier.parameters.proposeProfileLengthChanges = true;
        tests[ TEST_ID_lengthadjust_unconditional ].parametersModifier.isModified_proposeProfileLengthChanges = true;
        tests[ TEST_ID_lengthadjust_unconditional ].parametersModifier.parameters.useUnconditionalBaumWelch = true;
        tests[ TEST_ID_lengthadjust_unconditional ].parametersModifier.isModified_useUnconditionalBaumWelch = true;
      } // End if testLengthadjust

      // mark
      if( m_parameters.testBaldi ) {
        // Always test conditional Baldi if testBaldi is true.
        tests[ TEST_ID_baldi_conditional ].name = "conditionalBaldi";
        tests[ TEST_ID_baldi_conditional ].isRun = true;
        // Always write conditional BW results to STDOUT
        tests[ TEST_ID_baldi_conditional ].isCout = true;
        tests[ TEST_ID_baldi_conditional ].isGibbs = false;
        tests[ TEST_ID_baldi_conditional ].coutLeftBrace = "(~";
        tests[ TEST_ID_baldi_conditional ].coutRightBrace = "~)";
        // Start with the starting profile globals and positions
        tests[ TEST_ID_baldi_conditional ].startingGlobalsTest = NULL;
        tests[ TEST_ID_baldi_conditional ].startingPositionsTest = NULL;
        // Modifiers for conditional Baldi test
        tests[ TEST_ID_baldi_conditional ].parametersModifier.parameters.useUnconditionalBaumWelch = false;
        tests[ TEST_ID_baldi_conditional ].parametersModifier.isModified_useUnconditionalBaumWelch = true;
        tests[ TEST_ID_baldi_conditional ].parametersModifier.parameters.baldiLearningRate = 1;
        tests[ TEST_ID_baldi_conditional ].parametersModifier.isModified_baldiLearningRate = true;
        tests[ TEST_ID_baldi_conditional ].parametersModifier.parameters.baldiTemperature = 1;
        tests[ TEST_ID_baldi_conditional ].parametersModifier.isModified_baldiTemperature = true;
        tests[ TEST_ID_baldi_conditional ].parametersModifier.parameters.baldiHybrid = false;
        tests[ TEST_ID_baldi_conditional ].parametersModifier.isModified_baldiHybrid = true;
        tests[ TEST_ID_baldi_conditional ].parametersModifier.parameters.siegelMaxFindingThePeakAttempts_positions = 0;
        tests[ TEST_ID_baldi_conditional ].parametersModifier.isModified_siegelMaxFindingThePeakAttempts_positions = true;
        tests[ TEST_ID_baldi_conditional ].parametersModifier.parameters.minBaumWelchInverseScalar = 0;
        tests[ TEST_ID_baldi_conditional ].parametersModifier.isModified_minBaumWelchInverseScalar = true;
        tests[ TEST_ID_baldi_conditional ].parametersModifier.parameters.maxBaumWelchInverseScalar = 0;
        tests[ TEST_ID_baldi_conditional ].parametersModifier.isModified_maxBaumWelchInverseScalar = true;
  
        // TODO: REMOVE
        //tests[ TEST_ID_baldi_conditional ].parametersModifier.parameters.verbosity = VERBOSITY_All;
        //tests[ TEST_ID_baldi_conditional ].parametersModifier.isModified_verbosity = true;
  
        if( m_parameters.testUnconditionalProfile ) {
          tests[ TEST_ID_baldi_unconditional ].name = "unconditionalBaldi";
          tests[ TEST_ID_baldi_unconditional ].isRun = true;
        } else {
          tests[ TEST_ID_baldi_unconditional ].isRun = false;
        }
        if( m_parameters.coutUnconditionalProfile ) {
          tests[ TEST_ID_baldi_unconditional ].isCout = true;
        } else {
          tests[ TEST_ID_baldi_unconditional ].isCout = false;
        }
        tests[ TEST_ID_baldi_unconditional ].isGibbs = false;
        tests[ TEST_ID_baldi_unconditional ].coutLeftBrace = "<~";
        tests[ TEST_ID_baldi_unconditional ].coutRightBrace = "~>";
        // Start with the starting profile globals and positions
        tests[ TEST_ID_baldi_unconditional ].startingGlobalsTest = NULL;
        tests[ TEST_ID_baldi_unconditional ].startingPositionsTest = NULL;
        // Modifiers for unconditional Baldi test
        tests[ TEST_ID_baldi_unconditional ].parametersModifier.parameters.useUnconditionalBaumWelch = true;
        tests[ TEST_ID_baldi_unconditional ].parametersModifier.isModified_useUnconditionalBaumWelch = true;
        tests[ TEST_ID_baldi_unconditional ].parametersModifier.parameters.baldiLearningRate = 1;
        tests[ TEST_ID_baldi_unconditional ].parametersModifier.isModified_baldiLearningRate = true;
        tests[ TEST_ID_baldi_unconditional ].parametersModifier.parameters.baldiTemperature = 1;
        tests[ TEST_ID_baldi_unconditional ].parametersModifier.isModified_baldiTemperature = true;
        tests[ TEST_ID_baldi_unconditional ].parametersModifier.parameters.baldiHybrid = false;
        tests[ TEST_ID_baldi_unconditional ].parametersModifier.isModified_baldiHybrid = true;
        tests[ TEST_ID_baldi_unconditional ].parametersModifier.parameters.siegelMaxFindingThePeakAttempts_positions = 0;
        tests[ TEST_ID_baldi_unconditional ].parametersModifier.isModified_siegelMaxFindingThePeakAttempts_positions = true;
        tests[ TEST_ID_baldi_unconditional ].parametersModifier.parameters.minBaumWelchInverseScalar = 0;
        tests[ TEST_ID_baldi_unconditional ].parametersModifier.isModified_minBaumWelchInverseScalar = true;
        tests[ TEST_ID_baldi_unconditional ].parametersModifier.parameters.maxBaumWelchInverseScalar = 0;
        tests[ TEST_ID_baldi_unconditional ].parametersModifier.isModified_maxBaumWelchInverseScalar = true;
  
        // TODO: REMOVE
        //tests[ TEST_ID_baldi_unconditional ].parametersModifier.parameters.verbosity = VERBOSITY_All;
        //tests[ TEST_ID_baldi_unconditional ].parametersModifier.isModified_verbosity = true;
        if( m_parameters.testLengthadjust ) {
          tests[ TEST_ID_baldi_lengthadjust_conditional ].name = "conditionalBaldiLA";
          tests[ TEST_ID_baldi_lengthadjust_conditional ].isRun = true;
          // Always write conditional BW results to STDOUT
          tests[ TEST_ID_baldi_lengthadjust_conditional ].isCout = true;
          tests[ TEST_ID_baldi_lengthadjust_conditional ].isGibbs = false;
          tests[ TEST_ID_baldi_lengthadjust_conditional ].coutLeftBrace = "(*~";
          tests[ TEST_ID_baldi_lengthadjust_conditional ].coutRightBrace = "~*)";
          // Start with the starting profile globals and positions
          tests[ TEST_ID_baldi_lengthadjust_conditional ].startingGlobalsTest = NULL;
          tests[ TEST_ID_baldi_lengthadjust_conditional ].startingPositionsTest = NULL;
          // Modifiers for lengthadjust_conditional Baldi test
          tests[ TEST_ID_baldi_lengthadjust_conditional ].parametersModifier.parameters.proposeProfileLengthChanges = true;
          tests[ TEST_ID_baldi_lengthadjust_conditional ].parametersModifier.isModified_proposeProfileLengthChanges = true;
          tests[ TEST_ID_baldi_lengthadjust_conditional ].parametersModifier.parameters.useUnconditionalBaumWelch = false;
          tests[ TEST_ID_baldi_lengthadjust_conditional ].parametersModifier.isModified_useUnconditionalBaumWelch = true;
          tests[ TEST_ID_baldi_lengthadjust_conditional ].parametersModifier.parameters.baldiLearningRate = 1;
          tests[ TEST_ID_baldi_lengthadjust_conditional ].parametersModifier.isModified_baldiLearningRate = true;
          tests[ TEST_ID_baldi_lengthadjust_conditional ].parametersModifier.parameters.baldiTemperature = 1;
          tests[ TEST_ID_baldi_lengthadjust_conditional ].parametersModifier.isModified_baldiTemperature = true;
          tests[ TEST_ID_baldi_lengthadjust_conditional ].parametersModifier.parameters.baldiHybrid = false;
          tests[ TEST_ID_baldi_lengthadjust_conditional ].parametersModifier.isModified_baldiHybrid = true;
          tests[ TEST_ID_baldi_lengthadjust_conditional ].parametersModifier.parameters.siegelMaxFindingThePeakAttempts_positions = 0;
          tests[ TEST_ID_baldi_lengthadjust_conditional ].parametersModifier.isModified_siegelMaxFindingThePeakAttempts_positions = true;
          tests[ TEST_ID_baldi_lengthadjust_conditional ].parametersModifier.parameters.minBaumWelchInverseScalar = 0;
          tests[ TEST_ID_baldi_lengthadjust_conditional ].parametersModifier.isModified_minBaumWelchInverseScalar = true;
          tests[ TEST_ID_baldi_lengthadjust_conditional ].parametersModifier.parameters.maxBaumWelchInverseScalar = 0;
          tests[ TEST_ID_baldi_lengthadjust_conditional ].parametersModifier.isModified_maxBaumWelchInverseScalar = true;
          
          // TODO: REMOVE
          //tests[ TEST_ID_baldi_lengthadjust_conditional ].parametersModifier.parameters.verbosity = VERBOSITY_All;
          //tests[ TEST_ID_baldi_lengthadjust_conditional ].parametersModifier.isModified_verbosity = true;
        } // End if testLengthadjust
          
        if( m_parameters.testLengthadjust ) {
          if( m_parameters.testUnconditionalProfile ) {
            tests[ TEST_ID_baldi_lengthadjust_unconditional ].name = "unconditionalBaldiLA";
            tests[ TEST_ID_baldi_lengthadjust_unconditional ].isRun = true;
          } else {
            tests[ TEST_ID_baldi_lengthadjust_unconditional ].isRun = false;
          }
          if( m_parameters.coutUnconditionalProfile ) {
            tests[ TEST_ID_baldi_lengthadjust_unconditional ].isCout = true;
          } else {
            tests[ TEST_ID_baldi_lengthadjust_unconditional ].isCout = false;
          }
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].isGibbs = false;
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].coutLeftBrace = "<*~";
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].coutRightBrace = "~*>";
          // Start with the starting profile globals and positions
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].startingGlobalsTest = NULL;
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].startingPositionsTest = NULL;
          // Modifiers for lengthadjust_unconditional Baldi test
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].parametersModifier.parameters.proposeProfileLengthChanges = true;
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].parametersModifier.isModified_proposeProfileLengthChanges = true;
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].parametersModifier.parameters.useUnconditionalBaumWelch = true;
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].parametersModifier.isModified_useUnconditionalBaumWelch = true;
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].parametersModifier.parameters.baldiLearningRate = 1;
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].parametersModifier.isModified_baldiLearningRate = true;
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].parametersModifier.parameters.baldiTemperature = 1;
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].parametersModifier.isModified_baldiTemperature = true;
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].parametersModifier.parameters.baldiHybrid = false;
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].parametersModifier.isModified_baldiHybrid = true;
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].parametersModifier.parameters.siegelMaxFindingThePeakAttempts_positions = 0;
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].parametersModifier.isModified_siegelMaxFindingThePeakAttempts_positions = true;
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].parametersModifier.parameters.minBaumWelchInverseScalar = 0;
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].parametersModifier.isModified_minBaumWelchInverseScalar = true;
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].parametersModifier.parameters.maxBaumWelchInverseScalar = 0;
          tests[ TEST_ID_baldi_lengthadjust_unconditional ].parametersModifier.isModified_maxBaumWelchInverseScalar = true;
        } // End if testLengthadjust
      } // End if testBaldi

      if( m_parameters.testBaldiSiegel ) {
        // Always test conditional Baldi / Siegel if testBaldiSiegel is true.
        tests[ TEST_ID_baldi_siegel_conditional ].name = "conditionalBaldiSiegel";
        tests[ TEST_ID_baldi_siegel_conditional ].isRun = true;
        // Always write conditional BW results to STDOUT
        tests[ TEST_ID_baldi_siegel_conditional ].isCout = true;
        tests[ TEST_ID_baldi_siegel_conditional ].isGibbs = false;
        tests[ TEST_ID_baldi_siegel_conditional ].coutLeftBrace = "(~~";
        tests[ TEST_ID_baldi_siegel_conditional ].coutRightBrace = "~~)";
        // Start with the starting profile globals and positions
        tests[ TEST_ID_baldi_siegel_conditional ].startingGlobalsTest = NULL;
        tests[ TEST_ID_baldi_siegel_conditional ].startingPositionsTest = NULL;
        // Modifiers for conditional Baldi / Siegel test
        tests[ TEST_ID_baldi_siegel_conditional ].parametersModifier.parameters.useUnconditionalBaumWelch = false;
        tests[ TEST_ID_baldi_siegel_conditional ].parametersModifier.isModified_useUnconditionalBaumWelch = true;
        tests[ TEST_ID_baldi_siegel_conditional ].parametersModifier.parameters.baldiLearningRate = 1;
        tests[ TEST_ID_baldi_siegel_conditional ].parametersModifier.isModified_baldiLearningRate = true;
        tests[ TEST_ID_baldi_siegel_conditional ].parametersModifier.parameters.baldiTemperature = 1;
        tests[ TEST_ID_baldi_siegel_conditional ].parametersModifier.isModified_baldiTemperature = true;
        tests[ TEST_ID_baldi_siegel_conditional ].parametersModifier.parameters.baldiHybrid = false;
        tests[ TEST_ID_baldi_siegel_conditional ].parametersModifier.isModified_baldiHybrid = true;
        //tests[ TEST_ID_baldi_siegel_conditional ].parametersModifier.parameters.siegelMaxFindingThePeakAttempts_positions = 10000;
        //tests[ TEST_ID_baldi_siegel_conditional ].parametersModifier.isModified_siegelMaxFindingThePeakAttempts_positions = true;
        //tests[ TEST_ID_baldi_siegel_conditional ].parametersModifier.parameters.siegelEpsilonScaleFactor = 2;
        //tests[ TEST_ID_baldi_siegel_conditional ].parametersModifier.isModified_siegelEpsilonScaleFactor = true;
        tests[ TEST_ID_baldi_siegel_conditional ].parametersModifier.parameters.minBaumWelchInverseScalar = 0;
        tests[ TEST_ID_baldi_siegel_conditional ].parametersModifier.isModified_minBaumWelchInverseScalar = true;
        tests[ TEST_ID_baldi_siegel_conditional ].parametersModifier.parameters.maxBaumWelchInverseScalar = 0;
        tests[ TEST_ID_baldi_siegel_conditional ].parametersModifier.isModified_maxBaumWelchInverseScalar = true;
  
        // TODO: REMOVE
        //tests[ TEST_ID_baldi_siegel_conditional ].parametersModifier.parameters.verbosity = VERBOSITY_All;
        //tests[ TEST_ID_baldi_siegel_conditional ].parametersModifier.isModified_verbosity = true;
  
        if( m_parameters.testUnconditionalProfile ) {
          tests[ TEST_ID_baldi_siegel_unconditional ].name = "unconditionalBaldiSiegel";
          tests[ TEST_ID_baldi_siegel_unconditional ].isRun = true;
        } else {
          tests[ TEST_ID_baldi_siegel_unconditional ].isRun = false;
        }
        if( m_parameters.coutUnconditionalProfile ) {
          tests[ TEST_ID_baldi_siegel_unconditional ].isCout = true;
        } else {
          tests[ TEST_ID_baldi_siegel_unconditional ].isCout = false;
        }
        tests[ TEST_ID_baldi_siegel_unconditional ].isGibbs = false;
        tests[ TEST_ID_baldi_siegel_unconditional ].coutLeftBrace = "<~~";
        tests[ TEST_ID_baldi_siegel_unconditional ].coutRightBrace = "~~>";
        // Start with the starting profile globals and positions
        tests[ TEST_ID_baldi_siegel_unconditional ].startingGlobalsTest = NULL;
        tests[ TEST_ID_baldi_siegel_unconditional ].startingPositionsTest = NULL;
        // Modifiers for unconditional Baldi / Siegel test
        tests[ TEST_ID_baldi_siegel_unconditional ].parametersModifier.parameters.useUnconditionalBaumWelch = true;
        tests[ TEST_ID_baldi_siegel_unconditional ].parametersModifier.isModified_useUnconditionalBaumWelch = true;
        tests[ TEST_ID_baldi_siegel_unconditional ].parametersModifier.parameters.baldiLearningRate = 1;
        tests[ TEST_ID_baldi_siegel_unconditional ].parametersModifier.isModified_baldiLearningRate = true;
        tests[ TEST_ID_baldi_siegel_unconditional ].parametersModifier.parameters.baldiTemperature = 1;
        tests[ TEST_ID_baldi_siegel_unconditional ].parametersModifier.isModified_baldiTemperature = true;
        tests[ TEST_ID_baldi_siegel_unconditional ].parametersModifier.parameters.baldiHybrid = false;
        tests[ TEST_ID_baldi_siegel_unconditional ].parametersModifier.isModified_baldiHybrid = true;
        //tests[ TEST_ID_baldi_siegel_unconditional ].parametersModifier.parameters.siegelMaxFindingThePeakAttempts_positions = 10000;
        //tests[ TEST_ID_baldi_siegel_unconditional ].parametersModifier.isModified_siegelMaxFindingThePeakAttempts_positions = true;
        //tests[ TEST_ID_baldi_siegel_unconditional ].parametersModifier.parameters.siegelEpsilonScaleFactor = 2;
        //tests[ TEST_ID_baldi_siegel_unconditional ].parametersModifier.isModified_siegelEpsilonScaleFactor = true;
        tests[ TEST_ID_baldi_siegel_unconditional ].parametersModifier.parameters.minBaumWelchInverseScalar = 0;
        tests[ TEST_ID_baldi_siegel_unconditional ].parametersModifier.isModified_minBaumWelchInverseScalar = true;
        tests[ TEST_ID_baldi_siegel_unconditional ].parametersModifier.parameters.maxBaumWelchInverseScalar = 0;
        tests[ TEST_ID_baldi_siegel_unconditional ].parametersModifier.isModified_maxBaumWelchInverseScalar = true;
  
        // TODO: REMOVE
        //tests[ TEST_ID_baldi_siegel_unconditional ].parametersModifier.parameters.verbosity = VERBOSITY_All;
        //tests[ TEST_ID_baldi_siegel_unconditional ].parametersModifier.isModified_verbosity = true;

        if( m_parameters.testLengthadjust ) {
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].name = "conditionalBaldiSiegelLA";
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].isRun = true;
          // Always write conditional BW results to STDOUT
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].isCout = true;
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].isGibbs = false;
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].coutLeftBrace = "(*~~";
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].coutRightBrace = "~~*)";
          // Start with the starting profile globals and positions
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].startingGlobalsTest = NULL;
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].startingPositionsTest = NULL;
          // Modifiers for lengthadjust_conditional Baldi / Siegel test
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].parametersModifier.parameters.proposeProfileLengthChanges = true;
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].parametersModifier.isModified_proposeProfileLengthChanges = true;
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].parametersModifier.parameters.useUnconditionalBaumWelch = false;
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].parametersModifier.isModified_useUnconditionalBaumWelch = true;
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].parametersModifier.parameters.baldiLearningRate = 1;
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].parametersModifier.isModified_baldiLearningRate = true;
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].parametersModifier.parameters.baldiTemperature = 1;
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].parametersModifier.isModified_baldiTemperature = true;
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].parametersModifier.parameters.baldiHybrid = false;
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].parametersModifier.isModified_baldiHybrid = true;
          //tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].parametersModifier.parameters.siegelMaxFindingThePeakAttempts_positions = 10000;
          //tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].parametersModifier.isModified_siegelMaxFindingThePeakAttempts_positions = true;
          //tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].parametersModifier.parameters.siegelEpsilonScaleFactor = 2;
          //tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].parametersModifier.isModified_siegelEpsilonScaleFactor = true;
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].parametersModifier.isModified_minBaumWelchInverseScalar = true;
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].parametersModifier.parameters.maxBaumWelchInverseScalar = 0;
          tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].parametersModifier.isModified_maxBaumWelchInverseScalar = true;
          
          // TODO: REMOVE
          //tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].parametersModifier.parameters.verbosity = VERBOSITY_All;
          //tests[ TEST_ID_baldi_siegel_lengthadjust_conditional ].parametersModifier.isModified_verbosity = true;
        } // End if testLengthadjust
        
        if( m_parameters.testLengthadjust ) {
          if( m_parameters.testUnconditionalProfile ) {
            tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].name = "unconditionalBaldiSiegelLA";
            tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].isRun = true;
          } else {
            tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].isRun = false;
          }
          if( m_parameters.coutUnconditionalProfile ) {
            tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].isCout = true;
          } else {
            tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].isCout = false;
          }
          tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].isGibbs = false;
          tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].coutLeftBrace = "<*~~";
          tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].coutRightBrace = "~~*>";
          // Start with the starting profile globals and positions
          tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].startingGlobalsTest = NULL;
          tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].startingPositionsTest = NULL;
          // Modifiers for lengthadjust_unconditional Baldi / Siegel test
          tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].parametersModifier.parameters.proposeProfileLengthChanges = true;
          tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].parametersModifier.isModified_proposeProfileLengthChanges = true;
          tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].parametersModifier.parameters.useUnconditionalBaumWelch = true;
          tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].parametersModifier.isModified_useUnconditionalBaumWelch = true;
          tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].parametersModifier.parameters.baldiLearningRate = 1;
          tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].parametersModifier.isModified_baldiLearningRate = true;
          tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].parametersModifier.parameters.baldiTemperature = 1;
          tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].parametersModifier.isModified_baldiTemperature = true;
          tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].parametersModifier.parameters.baldiHybrid = false;
          tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].parametersModifier.isModified_baldiHybrid = true;
          //tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].parametersModifier.parameters.siegelMaxFindingThePeakAttempts_positions = 10000;
          //tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].parametersModifier.isModified_siegelMaxFindingThePeakAttempts_positions = true;
          //tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].parametersModifier.parameters.siegelEpsilonScaleFactor = 2;
          //tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].parametersModifier.isModified_siegelEpsilonScaleFactor = true;
          tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].parametersModifier.parameters.minBaumWelchInverseScalar = 0;
          tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].parametersModifier.isModified_minBaumWelchInverseScalar = true;
          tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].parametersModifier.parameters.maxBaumWelchInverseScalar = 0;
          tests[ TEST_ID_baldi_siegel_lengthadjust_unconditional ].parametersModifier.isModified_maxBaumWelchInverseScalar = true;
        } // End if testLengthadjust
      } // End if testBaldiSiegel
      // endmark


      uint32_t test_id;
    
      // Now go back through the tests and turn on tests[ test_num ].isRun for
      // any test that is required for another test with tests[ .. ].isRun ==
      // true.
      bool cycle_run_tests_again = true;
      while ( cycle_run_tests_again ) {
        cycle_run_tests_again = false;
        for( test_id = 0; test_id <= LAST_TEST_ID; test_id++ ) {
    // TODO: REMOVE
          //cout << "Here I am about to crash: test_id is " << test_id << endl;
          //cout << "Test name: " << tests[ test_id ].name << endl;

          if( !tests[ test_id ].isRun ) {
            continue;
          }
            // TODO: REMOVE
          //cout << "Test is run!" << endl;
          if( tests[ test_id ].startingGlobalsTest &&
              !tests[ test_id ].startingGlobalsTest->isRun ) {
            // TODO: REMOVE
            //cout << "Test is run with startingGlobalsTest that heretofore wouldn't have been run..." << endl;
            tests[ test_id ].startingGlobalsTest->isRun = true;
            cycle_run_tests_again = true;
          }
          if( tests[ test_id ].startingPositionsTest &&
              !tests[ test_id ].startingPositionsTest->isRun ) {
            // TODO: REMOVE
            //cout << "Test is run with startingPositionsTest that heretofore wouldn't have been run..." << endl;
            tests[ test_id ].startingPositionsTest->isRun = true;
            cycle_run_tests_again = true;
          }
        } // End foreach test_id
      } // End while cycle_run_tests_again

      uint32_t last_test_id = tests.size() - 1;

      if( be_verbose ) {
        cout << "Training with seed " << m_random.getSeed() << ":" << endl;
      } // End if be_verbose

      // Results
      string run_unique_id =
        ( "v" + boost::lexical_cast<string>( m_parameters.saveFileVersion ) + "_seed" + boost::lexical_cast<string>( m_random.getSeed() ) + "_ProfuseTest2" );
      fs::path dirname =
        ( static_cast<fs::path>( m_parameters.saveResultsParentDirectory ) /
          run_unique_id );
      std::ofstream tab_stream;
      if(m_parameters.saveResultsToFile ) {
        if( !fs::exists(m_parameters.saveResultsParentDirectory ) ) {
          cout << "Creating directory " << m_parameters.saveResultsParentDirectory << "..";
          cout.flush();
          fs::create_directory(m_parameters.saveResultsParentDirectory);
          cout << ".done." << endl;
//        } else {
//          cout << "Directory " << m_parameters.saveResultsParentDirectory << " exists." << endl;
        }
        if( !fs::exists( dirname ) ) {
          cout << "Creating directory " << dirname << "..";
          cout.flush();
          fs::create_directory( dirname );
          cout << ".done." << endl;
        } else {
          cout << "Directory " << dirname << " exists." << endl;
          cout << "Please use a different seed or erase that directory." << endl;
          return;
        }

        // TODO: REPLACE by saving a reusable ProfuseTest.cfg file with all of the current options.
        // // Parameters file
        // fs::path parameters_filename =
        //   ( m_parameters.resultsFilePrefix +
        //     run_unique_id +
        //     m_parameters.parametersFileSuffix );
        // ofstream parameters_stream( ( dirname / parameters_filename ).string().c_str() );
        // assert( parameters_stream.good() );
        // boost::archive::xml_oarchive parameters_oa( parameters_stream );
        // //parameters_oa << BOOST_SERIALIZATION_NVP( m_parameters );  //fixyfix
        // //parameters_oa << BOOST_SERIALIZATION_NVP( m_parameters );
        // //parameters_oa << 
        // //  boost::serialization::make_nvp( "m_parameters", m_parameters );
        // //parameters_oa << 
        // //  boost::serialization::make_nvp( "Parameters", m_parameters );
        // //parameters_oa.close();
        // //parameters_stream.close();

        if( m_parameters.saveTests ) {
          fs::path tests_filename =
            ( m_parameters.resultsFilePrefix +
              run_unique_id +
              m_parameters.testsFileSuffix );
              writeXML( tests, ( dirname / tests_filename ).string().c_str() );
        } // End if saveTests

        fs::path tab_filename =
          ( m_parameters.resultsFilePrefix +
            run_unique_id +
            m_parameters.tabFileSuffix );
        cout << endl << "For filename: " << tab_filename << endl;
        tab_stream.open( ( dirname / tab_filename ).string().c_str() );

        // print the tab header
        tab_stream << "conservation_rate\t";
        tab_stream << "profile_length\t";
        tab_stream << "num_training_sequences_per_profile\t";
        tab_stream << "expected_deletions_count\t";
        tab_stream << "expected_insertions_count\t";
        tab_stream << "expected_deletion_length_as_profile_length_fraction\t";
        tab_stream << "expected_insertion_length_as_profile_length_fraction\t";
        if( m_parameters.numTrueProfiles > 1 ) {
          tab_stream << "true_profile_id\t";
          cout << "true_profile_id ";
        }

        if( m_parameters.calculateSymmeterizedKullbackLeiblerDistancesToTrue ) {
          for( test_id = 0; test_id <= last_test_id; test_id++ ) {
            if( tests[ test_id ].isRun ) {
              tab_stream << "SKLPositions_" << tests[ test_id ].name << "_to_true\t";
              if( ( test_id != TEST_ID_true ) && m_parameters.coutDistances && tests[ test_id ].isCout ) {
                cout << "SKLPositions_" << tests[ test_id ].name << "_to_true ";
              }
              tab_stream << "SKLExceptPositions_" << tests[ test_id ].name << "_to_true\t";
              if( ( test_id != TEST_ID_true ) && m_parameters.coutDistances && tests[ test_id ].isCout ) {
                cout << "SKLExceptPositions_" << tests[ test_id ].name << "_to_true ";
              }
            }
          } // End foreach test_id, print "SKL*_[name]_to_true\t";
        } // End if calculateSymmeterizedKullbackLeiblerDistancesToTrue
        if( m_parameters.calculateSymmeterizedKullbackLeiblerDistancesToStarting ) {
          for( test_id = 0; test_id <= last_test_id; test_id++ ) {
            if( tests[ test_id ].isRun ) {
              tab_stream << "SKLPositions_" << tests[ test_id ].name << "_to_starting\t";
              if( ( !m_parameters.calculateSymmeterizedKullbackLeiblerDistancesToTrue || ( test_id != TEST_ID_true ) ) && ( test_id != TEST_ID_starting ) && m_parameters.coutDistances && tests[ test_id ].isCout ) {
                cout << "SKLPositions_" << tests[ test_id ].name << "_to_starting ";
              }
              tab_stream << "SKLExceptPositions_" << tests[ test_id ].name << "_to_starting\t";
              if( ( !m_parameters.calculateSymmeterizedKullbackLeiblerDistancesToTrue || ( test_id != TEST_ID_true ) ) && ( test_id != TEST_ID_starting ) && m_parameters.coutDistances && tests[ test_id ].isCout ) {
                cout << "SKLExceptPositions_" << tests[ test_id ].name << "_to_starting ";
              }
            } // End if isRun
          } // End foreach test_id, print "SKL*_[name]_to_starting\t";
        } // End if calculateSymmeterizedKullbackLeiblerDistancesToStarting

        if( m_parameters.calculateProfileProfileAlignments ) {
          for( test_id = 0; test_id <= last_test_id; test_id++ ) {
            if(
              ( test_id == TEST_ID_true ) ||
              ( test_id == TEST_ID_starting ) ||
              !tests[ test_id ].isRun
            ) {
              continue;
            }
            tab_stream << "PPAlign_SKL_" << tests[ test_id ].name << "_to_true\t";
            if( tests[ test_id ].isCout ) {
              cout << "PPAlign_SKL_" << tests[ test_id ].name << "_to_true ";
            }
            tab_stream << "PPAlign_numPositions_" << tests[ test_id ].name << "_to_true\t";
            if( tests[ test_id ].isCout ) {
              cout << "PPAlign_numPositions_" << tests[ test_id ].name << "_to_true ";
            }
            if( true ) {
              tab_stream << "length_" << tests[ test_id ].name << "\t";
              if( tests[ test_id ].isCout ) {
                cout << "length_" << tests[ test_id ].name << " ";
              }
            } // End if true
          } // End foreach test_id, ...
        } // End if calculateProfileProfileAlignments

        for( test_id = 0; test_id <= last_test_id; test_id++ ) {
          if( tests[ test_id ].isRun ) {
            if( !( ( test_id == TEST_ID_true ) || ( test_id == TEST_ID_starting ) ) ) {
              tab_stream << "training_" << tests[ test_id ].name << "_iters\t";
              if( tests[ test_id ].isCout ) {
                cout << "training_" << tests[ test_id ].name << "_iters";
                cout << " ";
              }
              tab_stream << "training_" << tests[ test_id ].name << "_CPUtime\t";
              if( tests[ test_id ].isCout ) {
                cout << "training_" << tests[ test_id ].name << "_CPUtime";
                cout << " ";
              }
            } // End if !( TEST_ID_true || TEST_ID_starting )
            tab_stream << "training_" << tests[ test_id ].name << "_forward\t";
            if( tests[ test_id ].isCout ) {
              cout << tests[ test_id ].coutLeftBrace;
              cout << "training_" << tests[ test_id ].name << "_forward";
              cout << tests[ test_id ].coutRightBrace;
              cout << " ";
            }
            if( tests[ test_id ].isGibbs ) {
              if( m_parameters.reportGibbsMean ) {
                tab_stream << "training_" << tests[ test_id ].name << "_mean_forward\t";
              }
              if( m_parameters.reportGibbsMode ) {
                tab_stream << "training_" << tests[ test_id ].name << "_mode_forward\t";
              }
              if( tests[ test_id ].isCout ) {
                if( m_parameters.reportGibbsMean ) {
                  cout << tests[ test_id ].coutLeftBrace;
                  cout << "training_" << tests[ test_id ].name << "_mean_forward";
                  cout << tests[ test_id ].coutRightBrace;
                  cout << " ";
                } // End if reportGibbsMean
                if( m_parameters.reportGibbsMode ) {
                  cout << tests[ test_id ].coutLeftBrace;
                  cout << "training_" << tests[ test_id ].name << "_mode_forward";
                  cout << tests[ test_id ].coutRightBrace;
                  cout << " ";
                } // End if reportGibbsMode
              }
            } // End if isGibbs
          } // End if isRun
        } // End foreach test_id, print "training_[name]_forward\t";

        if( m_parameters.testViterbi ) {
          for( test_id = 0; test_id <= last_test_id; test_id++ ) {
            if( tests[ test_id ].isRun ) {
              tab_stream << "training_" << tests[ test_id ].name << "_viterbi\t";
              if( m_parameters.coutViterbi && tests[ test_id ].isCout ) {
                cout << "training_" << tests[ test_id ].name << "_viterbi ";
              }
            }
          } // End foreach test_id, print "training_[name]_viterbi\t";

        } // End if testViterbi

        if( m_parameters.testTruepath ) {
          for( test_id = 0; test_id <= last_test_id; test_id++ ) {
            if( tests[ test_id ].isRun ) {
              tab_stream << "training_" << tests[ test_id ].name << "_truepath\t";
              if( m_parameters.coutTruepath && tests[ test_id ].isCout ) {
                cout << "training_" << tests[ test_id ].name << "_truepath ";
              }
            }
          } // End foreach test_id, print "training_[name]_truepath\t";

        } // End if testTruepath

        cout << "<==--#\t" << endl << "#--==>\t";
        for( test_id = 0; test_id <= last_test_id; test_id++ ) {
          if( tests[ test_id ].isRun ) {
            tab_stream << "test_" << tests[ test_id ].name << "_forward\t";
            if( tests[ test_id ].isCout ) {
              cout << tests[ test_id ].coutLeftBrace;
              cout << tests[ test_id ].coutLeftBrace;
              cout << "test_" << tests[ test_id ].name << "_forward";
              cout << tests[ test_id ].coutRightBrace;
              cout << tests[ test_id ].coutRightBrace;
              cout << " ";
            }
            if( tests[ test_id ].isGibbs ) {
              if( m_parameters.reportGibbsMean ) {
                tab_stream << "test_" << tests[ test_id ].name << "_mean_forward\t";
              }
              if( m_parameters.reportGibbsMode ) {
                tab_stream << "test_" << tests[ test_id ].name << "_mode_forward\t";
              }
              if( tests[ test_id ].isCout ) {
                if( m_parameters.reportGibbsMean ) {
                  cout << tests[ test_id ].coutLeftBrace;
                  cout << tests[ test_id ].coutLeftBrace;
                  cout << "test_" << tests[ test_id ].name << "_mean_forward";
                  cout << tests[ test_id ].coutRightBrace;
                  cout << tests[ test_id ].coutRightBrace;
                  cout << " ";
                } // End if reportGibbsMean
                if( m_parameters.reportGibbsMode ) {
                  cout << tests[ test_id ].coutLeftBrace;
                  cout << tests[ test_id ].coutLeftBrace;
                  cout << "test_" << tests[ test_id ].name << "_mode_forward";
                  cout << tests[ test_id ].coutRightBrace;
                  cout << tests[ test_id ].coutRightBrace;
                  cout << " ";
                } // End if reportGibbsMode
              } // End if isCout
            } // End if isGibbs
          } // End if isRun
        } // End foreach test_id, print "test_[name]_forward\t";

        if( m_parameters.testViterbi ) {
          for( test_id = 0; test_id <= last_test_id; test_id++ ) {
            if( tests[ test_id ].isRun ) {
              tab_stream << "test_" << tests[ test_id ].name << "_viterbi\t";
              if( m_parameters.coutViterbi && tests[ test_id ].isCout ) {
                cout << "test_" << tests[ test_id ].name << "_viterbi ";
              }
            }
          } // End foreach test_id, print "test_[name]_viterbi\t";

        } // End if testViterbi

        if( m_parameters.testTruepath ) {
          for( test_id = 0; test_id <= last_test_id; test_id++ ) {
            if( tests[ test_id ].isRun ) {
              tab_stream << "test_" << tests[ test_id ].name << "_truepath\t";
              if( m_parameters.coutTruepath && tests[ test_id ].isCout ) {
                cout << "test_" << tests[ test_id ].name << "_truepath ";
              }
            }
          } // End foreach test_id, print "test_[name]_truepath\t";

        } // End if testTruepath

        tab_stream << endl;
        cout << endl;
      } // End if m_parameters.saveResultsToFile

      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun and
       * params.calculateSymmeterizedKullbackLeiblerDistancesTo? are true, then
       * selfEntropyPositions[ test_num ] will hold the self-entropy of the
       * positions of the profile corresponding to test test_num.
       */
      double selfEntropyPositions[ last_test_id + 1 ];

      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun and
       * params.calculateSymmeterizedKullbackLeiblerDistancesTo? are true, then
       * selfEntropyExceptPositions[ test_num ] will hold the self-entropy of
       * the the non-position distributions of the profile corresponding to
       * test test_num.
       */
      double selfEntropyExceptPositions[ last_test_id + 1 ];

      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun and
       * params.calculateSymmeterizedKullbackLeiblerDistancesToTrue are true,
       * then distanceSKLPositions_true[ test_num ] will hold the symmetrized
       * Kullback-Leibler divergence between the positions of the profile
       * corresponding to test test_num and the positions of the true profile.
       */
      double distanceSKLPositions_true[ last_test_id + 1 ];
  
      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun and
       * params.calculateSymmeterizedKullbackLeiblerDistancesToTrue are true,
       * then distanceSKLExceptPositions_true[ test_num ] will hold the
       * symmetrized Kullback-Leibler divergence between the non-position
       * distributions of the profile corresponding to test test_num and those
       * distributions of the true profile.
       */
      double distanceSKLExceptPositions_true[ last_test_id + 1 ];
  
      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun and
       * params.calculateProfileProfileAlignments and
       * params.calculateSymmeterizedKullbackLeiblerDistancesToTrue are true,
       * then distanceSKLPositions_aligned_true[ test_num ] will hold the
       * symmetrized Kullback-Leibler divergence between just the aligned
       * positions of the profile corresponding to test test_num and those
       * positions of the true profile.
       */
      double distanceSKLPositions_aligned_true[ last_test_id + 1 ];
  
      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun and
       * params.calculateProfileProfileAlignments are true, then
       * numProfileProfileAlignedPositions_true[ test_num ] will hold the
       * number of aligned positions in the SKL profile-profile alignment
       * between the profile corresponding to test test_num and the true
       * profile.  See also params.profileProfileIndelOpenCost and
       * params.profileProfileIndelExtensionCost.
       */
      uint32_t numProfileProfileAlignedPositions_true[ last_test_id + 1 ];
  
      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun and
       * params.calculateSymmeterizedKullbackLeiblerDistancesTo? are true, then
       * distanceSKLPositions_starting[ test_num ] will hold the symmetrized
       * Kullback-Leibler divergence between the positions of the profile
       * corresponding to test test_num and the positions of the starting
       * profile.
       */
      double distanceSKLPositions_starting[ last_test_id + 1 ];
  
      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun and
       * params.calculateSymmeterizedKullbackLeiblerDistancesTo? are true, then
       * distanceSKLExceptPositions_starting[ test_num ] will hold the
       * symmetrized Kullback-Leibler divergence between the non-position
       * distributions of the profile corresponding to test test_num and those
       * distributions the starting profile.
       */
      double distanceSKLExceptPositions_starting[ last_test_id + 1 ];
  
      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true, then testIters[ test_num ] will
       * hold the number of iterations of training (or Gibbs) performed (in the
       * case of Gibbs, this includes unused iterations such as burn-in and
       * thinned-away iters too).
       */
      uint32_t testIters[ last_test_id + 1 ];
  
      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true, then testCPUtime[ test_num ] will
       * hold the total CPU time required (in the
       * case of Gibbs, this includes unused iterations such as burn-in and
       * thinned-away iters too).
       */
      double testCPUtime[ last_test_id + 1 ];
  
      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true, then testScore_training_forward[
       * test_num ] will hold the forward score for the corresponding test
       * number for the training set.
       */
      ScoreType testScore_training_forward[ last_test_id + 1 ];
  
      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true and tests[ test_num ].isGibbs is
       * true, then testScore_training_mean_forward[ test_num ] will hold the
       * forward score for mean profile for the corresponding test number for
       * the training set.  Note that testScore_training_forward[ test_num ]
       * will hold whichever is higher, this or
       * testScore_training_mode_forward[ test_num ].
       *
       * Only used if parameters.reportGibbs_mean is true.
       */
      ScoreType testScore_training_mean_forward[ last_test_id + 1 ];

      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true and tests[ test_num ].isGibbs is
       * true, then testScore_training_mode_forward[ test_num ] will hold the
       * forward score for mode profile for the corresponding test number for
       * the training set.  Note that testScore_training_forward[ test_num ]
       * will hold whichever is higher, this or
       * testScore_training_mean_forward[ test_num ].
       *
       * Only used if parameters.reportGibbs_mode is true.
       */
      ScoreType testScore_training_mode_forward[ last_test_id + 1 ];

      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true, then testScore_test_forward[
       * test_num ] will hold the forward score for the corresponding test
       * number for the test set.
       */
      ScoreType testScore_test_forward[ last_test_id + 1 ];
  
      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true and tests[ test_num ].isGibbs is
       * true, then testScore_test_mean_forward[ test_num ] will hold the
       * forward score for mean profile for the corresponding test number for
       * the test set.  Note that testScore_test_forward[ test_num ]
       * will hold whichever is higher, this or
       * testScore_test_mode_forward[ test_num ].
       *
       * Only used if parameters.reportGibbs_mean is true.
       */
      ScoreType testScore_test_mean_forward[ last_test_id + 1 ];

      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true and tests[ test_num ].isGibbs is
       * true, then testScore_test_mode_forward[ test_num ] will hold the
       * forward score for mode profile for the corresponding test number for
       * the test set.  Note that testScore_test_forward[ test_num ]
       * will hold whichever is higher, this or
       * testScore_test_mean_forward[ test_num ].
       *
       * Only used if parameters.reportGibbs_mode is true.
       */
      ScoreType testScore_test_mode_forward[ last_test_id + 1 ];

      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true, then testScore_training_viterbi[
       * test_num ] will hold the viterbi score for the corresponding test
       * number for the training set.
       */
      ScoreType testScore_training_viterbi[ last_test_id + 1 ];
  
      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true, then testScore_test_viterbi[
       * test_num ] will hold the viterbi score for the corresponding test
       * number for the test set.
       */
      ScoreType testScore_test_viterbi[ last_test_id + 1 ];

      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true, then testScore_training_truepath[
       * test_num ] will hold the truepath score for the corresponding test
       * number for the training set.
       */
      ScoreType testScore_training_truepath[ last_test_id + 1 ];
  
      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true, then testScore_test_truepath[
       * test_num ] will hold the truepath score for the corresponding test
       * number for the test set.
       */
      ScoreType testScore_test_truepath[ last_test_id + 1 ];

      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true, then testProfileTree[ test_num ]
       * will hold the profile tree (after training) for the corresponding test
       * number.
       */
      ProfileTreeType testProfileTree[ last_test_id + 1 ];
  
      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true and tests[ test_num ].isGibbs is
       * true, then testProfileTree_mean[ test_num ] will hold the mean profile
       * tree (after sampling) for the corresponding test number.
       *
       * Only used if parameters.reportGibbs_mean is true.
       */
      ProfileTreeType testProfileTree_mean[ last_test_id + 1 ];
  
      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true and tests[ test_num ].isGibbs is
       * true, then testProfileTree_mode[ test_num ] will hold the mode profile
       * tree (after sampling) for the corresponding test number.
       *
       * Only used if parameters.reportGibbs_mode is true.
       */
      ProfileTreeType testProfileTree_mode[ last_test_id + 1 ];
  
      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true, then testProfileSequences[
       * test_num ] will hold the vector of vectors of SequenceTypes
       * corresponding to the segmentation of the sequences into the profiles
       * trained in test test_num.
       */
      //vector<vector<SequenceType> > testProfileSequences[ last_test_id + 1 ];
      //for( test_id = 0; test_id <= last_test_id; test_id++ ) {
      //  testProfileSequences[ test_id ] =
      //    vector<vector<SequenceType> >( m_parameters.numProfiles );
      //}

      DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType> dp;
      MultinomialDistribution<ResidueType,ProbabilityType> residue_dist;
      Fasta<SequenceResidueType>  pattern_sequences( m_parameters.numProfiles );
      vector<Fasta<SequenceResidueType> > training_fastas( m_parameters.numProfiles );
      Fasta<SequenceResidueType>  training_fasta; // All together now
      vector<Fasta<SequenceResidueType> > testing_fastas( m_parameters.numProfiles );
      Fasta<SequenceResidueType>  testing_fasta; // All together now
                    
      RootTypeMultipleAlignment training_root_ma;
      InternalNodeTypeMultipleAlignment training_internal_node_ma;
      RootTypeMultipleAlignment testing_root_ma;
      InternalNodeTypeMultipleAlignment testing_internal_node_ma;

      map<Test *, ProfileTreeType *> testProfileTreePtrMap;
      for( test_id = 0; test_id <= last_test_id; test_id++ ) {
        testProfileTreePtrMap.insert( pair<Test *, ProfileTreeType *>( &tests[ test_id ], &testProfileTree[ test_id ] ) );
      } // End foreach test_id, set up test ptr -> root ptr map

      // TODO: Add a parameter (to the trainer) for using matrices vs. anchors.
      RowVector training_forward_rows_1;
      RowVector training_forward_rows_2;
      RowVector testing_forward_rows_1;
      RowVector testing_forward_rows_2;
      // OLD:
      //SequentialAccessContainer training_forward_matrices;
      //uint32_t longest_training_sequence_length, last_longest_training_sequence_length;
      //SequentialAccessContainer testing_forward_matrices;
      //uint32_t longest_testing_sequence_length, last_longest_testing_sequence_length;

      clock_t start_time, end_time; // For calculating the testCPUtime.
              
      for( uint32_t conservation_rate_i = 0;
           conservation_rate_i < ( m_parameters.conservationRates.size() >0 ? m_parameters.conservationRates.size() : 10U );
           conservation_rate_i++
      ) {
        double conservation_rate =
          ( m_parameters.conservationRates.size()>0 ?
            m_parameters.conservationRates[ conservation_rate_i ] :
            ( .1 * ( conservation_rate_i + 1 ) ) );
        for( uint32_t profile_length_i = 0;
             profile_length_i < ( m_parameters.profileLengths.size()>0 ? m_parameters.profileLengths.size() : 10U );
             profile_length_i++
        ) {
          uint32_t profile_length =
            ( m_parameters.profileLengths.size()>0 ?
              m_parameters.profileLengths[ profile_length_i ] :
              ( 10 * ( profile_length_i + 1 ) ) );

          // Do additional setup of the transition priors for those transitions
          // observed many times (since the number of transitions depends on
          // the profile length).
          if( m_parameters.usePriors || m_parameters.startWithGlobalsDrawnFromPrior ) {
            globalPrior[ 0 ][ Transition::fromMatch ][ TransitionFromMatch::toMatch ] = ( profile_length * m_parameters.priorStrength_internal_transitions * m_parameters.priorMtoM );
            globalPrior[ 0 ][ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] = ( profile_length * m_parameters.priorStrength_internal_transitions * m_parameters.priorMtoD );
#ifdef USE_DEL_IN_DEL_OUT
            // TODO: Create training_parameters_template.priorMtoW
            globalPrior[ 0 ][ Transition::fromMatch ][ TransitionFromMatch::toDeletionOut ] = ( profile_length * m_parameters.priorStrength_internal_transitions * .5 );
#endif // USE_DEL_IN_DEL_OUT
            globalPrior[ 0 ][ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] = ( profile_length * m_parameters.priorStrength_internal_transitions * m_parameters.priorMtoI );
            globalPrior[ 0 ][ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ] = ( profile_length * m_parameters.priorStrength_internal_transitions * m_parameters.priorItoM );
            globalPrior[ 0 ][ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] = ( profile_length * m_parameters.priorStrength_internal_transitions * m_parameters.priorItoI );
            globalPrior[ 0 ][ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ] = ( profile_length * m_parameters.priorStrength_internal_transitions * m_parameters.priorDtoM );
            globalPrior[ 0 ][ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ] = ( profile_length * m_parameters.priorStrength_internal_transitions * m_parameters.priorDtoD );

//            globalPrior[ 0 ][ Transition::fromMatch ][ TransitionFromMatch::toMatch ] = 1.0f + ( m_parameters.priorStrength * m_parameters.priorMtoM );
//            globalPrior[ 0 ][ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] = 1.0f + ( m_parameters.priorStrength * m_parameters.priorMtoD );
//            globalPrior[ 0 ][ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] = 1.0f + ( m_parameters.priorStrength * m_parameters.priorMtoI );
//            globalPrior[ 0 ][ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ] = 1.0f + ( m_parameters.priorStrength * m_parameters.priorItoM );
//            globalPrior[ 0 ][ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] = 1.0f + ( m_parameters.priorStrength * m_parameters.priorItoI );
//            globalPrior[ 0 ][ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ] = 1.0f + ( m_parameters.priorStrength * m_parameters.priorDtoM );
//            globalPrior[ 0 ][ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ] = 1.0f + ( m_parameters.priorStrength * m_parameters.priorDtoD );
          } // End if usePriors
          if( m_parameters.usePriors ) {
            m_parameters.matchEmissionPrior = &matchEmissionPrior;
            m_parameters.globalPrior = &globalPrior;
          }
          
          vector<uint32_t> profile_profile_alignment;
          double profile_profile_alignment_cost;
          for( uint32_t num_training_sequences_per_profile_i = 0;
               num_training_sequences_per_profile_i < ( m_parameters.numTrainingSequencesPerProfiles.size()>0 ? m_parameters.numTrainingSequencesPerProfiles.size() : 2U );
               num_training_sequences_per_profile_i++
          ) {
            uint32_t num_training_sequences_per_profile =
              ( m_parameters.numTrainingSequencesPerProfiles.size()>0 ?
                m_parameters.numTrainingSequencesPerProfiles[ num_training_sequences_per_profile_i ] :
                ( uint32_t )std::pow( 10.0, ( int )( num_training_sequences_per_profile_i + 1 ) ) );
            for( uint32_t expected_deletions_count_i = 0;
                 expected_deletions_count_i < ( m_parameters.expectedDeletionsCounts.size() > 0 ? m_parameters.expectedDeletionsCounts.size() : 1U );
                 expected_deletions_count_i++
            ) {
              double expected_deletions_count =
                ( m_parameters.expectedDeletionsCounts.size() > 0 ?
                  m_parameters.expectedDeletionsCounts[ expected_deletions_count_i ] :
                  1.0 );
              for( uint32_t expected_insertions_count_i = 0;
                   ( m_parameters.useDeletionsForInsertionsParameters ? ( expected_insertions_count_i == 0 ) : ( expected_insertions_count_i < ( m_parameters.expectedInsertionsCounts.size() > 0 ? m_parameters.expectedInsertionsCounts.size() : 1U ) ) );
                   expected_insertions_count_i++
              ) {
                double expected_insertions_count =
                  ( m_parameters.useDeletionsForInsertionsParameters ?
                    expected_deletions_count :
                    ( ( m_parameters.expectedInsertionsCounts.size() > 0 ?
                        m_parameters.expectedInsertionsCounts[ expected_insertions_count_i ] :
                        1.0 ) ) );
                for( uint32_t expected_deletion_length_as_profile_length_fraction_i = 0;
                     expected_deletion_length_as_profile_length_fraction_i < ( m_parameters.expectedDeletionLengthAsProfileLengthFractions.size() > 0 ? m_parameters.expectedDeletionLengthAsProfileLengthFractions.size() : 1U );
                     expected_deletion_length_as_profile_length_fraction_i++
                ) {
                  double expected_deletion_length_as_profile_length_fraction =
                    ( m_parameters.expectedDeletionLengthAsProfileLengthFractions.size() > 0 ?
                      m_parameters.expectedDeletionLengthAsProfileLengthFractions[ expected_deletion_length_as_profile_length_fraction_i ] :
                      1.0 );
                  for( uint32_t expected_insertion_length_as_profile_length_fraction_i = 0;
                       ( m_parameters.useDeletionsForInsertionsParameters ? ( expected_insertion_length_as_profile_length_fraction_i == 0 ) : ( expected_insertion_length_as_profile_length_fraction_i < (m_parameters.expectedInsertionLengthAsProfileLengthFractions.size() > 0 ? m_parameters.expectedInsertionLengthAsProfileLengthFractions.size() : 1U ) ) );
                       expected_insertion_length_as_profile_length_fraction_i++
                  ) {
                    double expected_insertion_length_as_profile_length_fraction =
                      ( m_parameters.useDeletionsForInsertionsParameters ?
                        expected_deletion_length_as_profile_length_fraction :
                        ( ( m_parameters.expectedInsertionLengthAsProfileLengthFractions.size() > 0 ?
                            m_parameters.expectedInsertionLengthAsProfileLengthFractions[ expected_insertion_length_as_profile_length_fraction_i ] :
                            1.0 ) ) );

                    // OLD:
                    //longest_training_sequence_length = 0;
                    //longest_testing_sequence_length = 0;
                    for( uint32_t true_profile_i = 0; true_profile_i < m_parameters.numTrueProfiles; true_profile_i++ ) {
                      // Generate a pattern sequence from the uniform
                      // distribution, of length profile_length.
                      pattern_sequences.reinitialize(
                        m_parameters.numProfiles
                      );
                      uint32_t which_profile;
                    
                      // TODO: Take a unique prefix?
                      pattern_sequences.m_descriptions[ 0 ] = "Root";
                      for( which_profile = 1; which_profile < m_parameters.numProfiles; which_profile++ ) {
                        pattern_sequences.m_descriptions[ which_profile ] =
                          ( "Child " + boost::lexical_cast<std::string>( which_profile ) );
                      } // End foreach which_profile
              
                        // Okay so if we want to have a shared position rate of x, and we are
                        // allowing y = 1/ResidueType::elementCount of the randomly-drawn positions to
                        // accidentally share the parent, then we need to choose not to draw randomly
                        // at a rate of z, where
                        // (1-z)y + z = x
                        // y - zy + z = x
                        // y + (1-y)z = x
                        // (1-y)z = x - y
                        // z = ( (x-y)/(1-y) )
                        double even_base_prob = ( 1.0 / seqan::ValueSize<ResidueType>::VALUE ); // y = .25
                        double shared_pos_rate_trick = // z
                          ( m_parameters.sharedPositionRate - even_base_prob ) / ( 1.0 - even_base_prob );
                      
                        // TODO: REMOVE
                        //cout << "shared position rate is " << m_parameters.sharedPositionRate << endl;
                        //cout << "shared pos rate trick is " << shared_pos_rate_trick << endl;
                      
                        residue_dist.even();
                        for( uint32_t pos_i = 0; pos_i < profile_length; pos_i++ ) {
                          seqan::appendValue(
                            pattern_sequences[ 0 ],
                            residue_dist.draw( m_random )
                          );
                          for( which_profile = 1; which_profile < m_parameters.numProfiles; which_profile++ ) {
                            if( ( m_parameters.sharedPositionRate == 1.0 ) ||
                                ( m_random.nextUniform() <= shared_pos_rate_trick ) ) {
                              // Shared, so the pattern of the child is the same as that of the root.
                              // TODO: This is assuming that the topology is flat: all children are
                              // children of the root.  If that should change, then this needs to
                              // change, too.
                              seqan::appendValue(
                                pattern_sequences[ which_profile ],
                                pattern_sequences[ 0 ][ pos_i ]
                              );
                            } else {
                              // Note that even here the position could be shared (if we happen to
                              // draw the same base as is in the root pos): see above where we
                              // calculate shared_pos_rate_trick.
                              seqan::appendValue(
                                pattern_sequences[ which_profile ],
                                residue_dist.draw( m_random )
                              );
                            } // End if this position is shared .. else ..
                          } // End foreach child
                        } // End foreach pos_i

                    
                      if( be_verbose ) {
                        cout << "Pattern sequences:" << endl;
                        cout << pattern_sequences << endl;
                      } // End if be_verbose
                    
                      if( m_parameters.saveResultsToFile && m_parameters.savePatternSequences ) {
                        fs::path pattern_sequences_filename =
                          ( m_parameters.resultsFilePrefix +
                            run_unique_id +
                            "." + lexical_cast<string>( conservation_rate * 100 ) +
                            "." + lexical_cast<string>( profile_length ) +
                            "." + lexical_cast<string>( num_training_sequences_per_profile ) +
                            "." + lexical_cast<string>( m_parameters.numTestingSequencesPerProfile ) +
                                "." + lexical_cast<string>( expected_deletions_count ) +
                                "." + lexical_cast<string>( expected_insertions_count ) +
                                "." + lexical_cast<string>( expected_deletion_length_as_profile_length_fraction ) +
                                "." + lexical_cast<string>( expected_insertion_length_as_profile_length_fraction ) +
                            "." + lexical_cast<string>( true_profile_i ) +
                            m_parameters.patternSequencesFileSuffix );
                        std::ofstream pattern_sequences_stream( ( dirname / pattern_sequences_filename ).string().c_str() );
                        assert( pattern_sequences_stream.good() );
                        pattern_sequences_stream << pattern_sequences;
                        pattern_sequences_stream.close();
                      } // End if saveResultsToFile
                    
                      // Make a root profile patterned thereon
                      testProfileTree[ TEST_ID_true ].reinitialize( profile_length );

                      // First add the children.  We want to create the whole topology before
                      // taking any reference to the children because modifying the topology
                      // invalidates references to the children (see ProfileTree.hpp for more on
                      // this).
                      for( which_profile = 1; which_profile < m_parameters.numProfiles; which_profile++ ) {
                        // For now make them all children of the root
                        // TODO: Random (or prescribed) alternate topologies
                        testProfileTree[ TEST_ID_true ].addChildToRoot(); // vertex is which_profile.
                      } // End foreach which_profile...
                    
                      RootType & true_root =
                        ( *testProfileTree[ TEST_ID_true ].getProfileTreeRoot() );
                      // TODO: REMOVE?
                      true_root.ensurePositionsKnowTheirRoot();

                      // "true" global transition params.
        
                      // First calculate the appropriate indel values.
                      ProbabilityType deletion_open =
                        INIT_PROBABILITY(ProbabilityType)( expected_deletions_count / profile_length );
                      ProbabilityType insertion_open =
                        ( m_parameters.useDeletionsForInsertionsParameters ?
                          deletion_open :
                          INIT_PROBABILITY(ProbabilityType)( expected_insertions_count / profile_length ) );
        
                      // [ the EV of a geometric is 1/p, where p is prob of stopping, so if q is the prob of continuing, we want ( 1 - q ) = 1/EV. ]
                      ProbabilityType deletion_extension =
                        INIT_PROBABILITY(ProbabilityType)( 1.0 - min( ( 1.0 / ( expected_deletion_length_as_profile_length_fraction * profile_length ) ), ( 1.0 / m_parameters.minExpectedDeletionLength ) ) );
                      ProbabilityType insertion_extension =
                        ( m_parameters.useDeletionsForInsertionsParameters ? deletion_extension : INIT_PROBABILITY(ProbabilityType)( 1.0 - min( ( 1.0 / ( expected_insertion_length_as_profile_length_fraction * profile_length ) ), ( 1.0 / m_parameters.minExpectedInsertionLength ) ) ) );

#ifdef USE_DEL_IN_DEL_OUT
                      // TODO: Make a parameter for this..
                      ProbabilityType deletion_in_open = .5;
                      // TODO: Make a parameter for this..
                      ProbabilityType deletion_in_extension = .5;
                      // TODO: Make a parameter for this..
                      ProbabilityType deletion_out_open = .5;
                      // TODO: Make a parameter for this..
                      ProbabilityType deletion_out_extension = .5;
#endif // USE_DEL_IN_DEL_OUT
                    
                      // Now set up the profile(s)
                      true_root.even();
#ifndef DISALLOW_FLANKING_TRANSITIONS
                      true_root[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] =
                        INIT_PROBABILITY(ProbabilityType)( m_parameters.preAlignInsertion );
                      true_root[ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ] =
                        INIT_PROBABILITY(ProbabilityType)( 1 ) -
                        true_root[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ];
#endif // !DISALLOW_FLANKING_TRANSITIONS
                      true_root[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ] =
                        deletion_open;
#ifdef USE_DEL_IN_DEL_OUT
                      true_root[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletionIn ] =
                        deletion_in_open;
                      true_root[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toMatch ] =
                        1.0 -
                        (
                          true_root[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletion ] +
                          true_root[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletionIn ]
                        );
                      true_root[ galosh::Transition::fromDeletionIn ][ galosh::TransitionFromDeletionIn::toDeletionIn ] =
                        deletion_in_extension;
                      true_root[ galosh::Transition::fromDeletionIn ][ galosh::TransitionFromDeletionIn::toMatch ] =
                        1.0 -
                        true_root[ galosh::Transition::fromDeletionIn ][ galosh::TransitionFromDeletionIn::toDeletionIn ];
#else
                      true_root[ Transition::fromBegin ][ TransitionFromBegin::toMatch ] =
                        INIT_PROBABILITY(ProbabilityType)( 1 ) -
                        true_root[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ];
#endif // USE_DEL_IN_DEL_OUT .. else ..                        
                      true_root[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] =
                        insertion_open;
                      true_root[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] =
                        deletion_open;
#ifdef USE_DEL_IN_DEL_OUT
                      true_root[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletionOut ] =
                        deletion_out_open;
                      true_root[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toMatch ] =
                        1.0 -
                        (
                          true_root[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toInsertion ] +
                          true_root[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletion ] +
                          true_root[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletionOut ]
                        );
                      true_root[ galosh::Transition::fromDeletionOut ][ galosh::TransitionFromDeletionOut::toDeletionOut ] =
                        deletion_out_extension;
                      true_root[ galosh::Transition::fromDeletionOut ][ galosh::TransitionFromDeletionOut::toEnd ] =
                        1.0 -
                        true_root[ galosh::Transition::fromDeletionOut ][ galosh::TransitionFromDeletionOut::toDeletionOut ];
#else
                      true_root[ Transition::fromMatch ][ TransitionFromMatch::toMatch ] =
                        INIT_PROBABILITY(ProbabilityType)( 1.0 ) -
                        (
                         true_root[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] +
                         true_root[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ]
                        );
#endif // USE_DEL_IN_DEL_OUT .. else ..
                      true_root[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] =
                        insertion_extension;
                      true_root[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ] =
                        INIT_PROBABILITY(ProbabilityType)( 1.0 ) -
                        true_root[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ];
                      true_root[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ] =
                        deletion_extension;
                      true_root[ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ] =
                        INIT_PROBABILITY(ProbabilityType)( 1.0 ) -
                        true_root[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ];
                    
#ifdef USE_END_DISTRIBUTION
                      // For now we don't use the End distribution ..
                      //true_root[ Transition::fromEnd ][ TransitionFromEnd::toPostAlign ] = INIT_PROBABILITY(ProbabilityType)( 1 );
                      //true_root[ Transition::fromEnd ][ TransitionFromEnd::toLoop ] = INIT_PROBABILITY(ProbabilityType)( 0 );
#endif // USE_END_DISTRIBUTION
#ifndef DISALLOW_FLANKING_TRANSITIONS
                      true_root[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] =
                        INIT_PROBABILITY(ProbabilityType)( m_parameters.postAlignInsertion );
                      true_root[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] =
                        INIT_PROBABILITY(ProbabilityType)( 1.0 ) -
                        true_root[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ];
#endif // DISALLOW_FLANKING_TRANSITIONS
              
                        // This is a trick to get the values set correctly:
                        // make an even profile, then set one of them higher than you need it, but
                        // normalize to get the right thing.
                        // NOTE: true_root[ j ][ pattern_sequences[ 0 ][ i ] is the same for all i,j right now.
                  
                        ProbabilityType pattern_trick_value =
                          ( ( conservation_rate == 1.0 ) ? INIT_PROBABILITY(ProbabilityType)( 1.0 ) :
                            ( ( INIT_PROBABILITY(ProbabilityType)( 1.0 ) - true_root[ 0 ][ Emission::Match ][ pattern_sequences[ 0 ][ 0 ] ] ) *
                              ( INIT_PROBABILITY(ProbabilityType)( conservation_rate / ( 1.0 - conservation_rate ) ) ) ) );
                      
                        // r is current remaining value (1 - P(base)), p is target value.
                        //x/(r + x) = p
                        //p( r + x) = x
                        // rp + px = x
                        // x - px = rp
                        // x( 1 - p ) = rp
                        // x = r( p / ( 1 - p ) )
                        //cout << "Pattern trick value " << pattern_trick_value << endl; 
                      
                        for( uint32_t pos_i = 0; pos_i < profile_length; pos_i++ ) {
                          if( conservation_rate == 1.0 ) {
                            true_root[ pos_i ][ Emission::Match ].zero();
                          }
                          true_root[ pos_i ][ Emission::Match ][ pattern_sequences[ 0 ][ pos_i ] ] = pattern_trick_value;
                          if( conservation_rate != 1.0 ) {
                            true_root[ pos_i ][ Emission::Match ].normalize( 0 );
                          }
                      
                          // Now modify the children
                          for( which_profile = 1; which_profile < m_parameters.numProfiles; which_profile++ ) {
                            if( pattern_sequences[ which_profile ][ pos_i ] ==
                                pattern_sequences[ 0 ][ pos_i ] ) {
                              // If the child shares the parent pattern at this position, just "useParentPosition".
                              // TODO: REMOVE
                              //cout << "[child " << which_profile << "]: Using parent position at pos " << pos_i << endl;
                              // TODO: Implement sharing positions
                              //testProfileTree[ TEST_ID_true ].getChild( 0, which_profile ).useParentPosition( pos_i );
                            } else {
                              // Otherwise, set the child pos up separately
                              if( conservation_rate == 1.0 ) {
                                testProfileTree[ TEST_ID_true ].getChild( 0, which_profile )[ pos_i ][ Emission::Match ].zero();
                              }
                              testProfileTree[ TEST_ID_true ].getChild( 0, which_profile )[ pos_i ][ Emission::Match ][ pattern_sequences[ which_profile ][ pos_i ] ] = pattern_trick_value;
                              if( conservation_rate != 1.0 ) {
                                testProfileTree[ TEST_ID_true ].getChild( 0, which_profile )[ pos_i ][ Emission::Match ].normalize( 0 );
                              }
                            } // End if this child shares the root's pattern at this position .. else ..
                          } // End foreach which_profile...
                        } // End foreach position, set it up according to the pattern and conservation_rate.

              
                      if( be_verbose ) {
                        // TODO: REMOVE?
                        cout << "> " << pattern_sequences.m_descriptions[ 0 ] << endl;
                        cout << true_root << endl;
                        for( which_profile = 1; which_profile < m_parameters.numProfiles; which_profile++ ) {
                          cout << "> " << pattern_sequences.m_descriptions[ which_profile ] << endl;
                          cout << testProfileTree[ TEST_ID_true ].getChild( 0, which_profile ) << endl;
                        } // End foreach which_profile ..
                      } // End if be_verbose
                    
                      if( m_parameters.saveResultsToFile && m_parameters.saveTrueProfileTrees ) {
                        fs::path true_profile_tree_filename =
                          ( m_parameters.resultsFilePrefix +
                            run_unique_id +
                            "." + lexical_cast<string>( conservation_rate * 100 ) +
                            "." + lexical_cast<string>( profile_length ) +
                            "." + lexical_cast<string>( num_training_sequences_per_profile ) +
                            "." + lexical_cast<string>( m_parameters.numTestingSequencesPerProfile ) +
                                "." + lexical_cast<string>( expected_deletions_count ) +
                                "." + lexical_cast<string>( expected_insertions_count ) +
                                "." + lexical_cast<string>( expected_deletion_length_as_profile_length_fraction ) +
                                "." + lexical_cast<string>( expected_insertion_length_as_profile_length_fraction ) +
                            "." + lexical_cast<string>( true_profile_i ) +
                            m_parameters.trueProfileTreeFileSuffix );
                        writeXML( testProfileTree[ TEST_ID_true ], ( dirname / true_profile_tree_filename ).string().c_str() );
                      } // End if saveResultsToFile
              
              
                      // Now draw the training sequences for the root profile.

                      // One fasta per profile
                      dp.drawSequences(
                        m_parameters,
                        true_root,
                        num_training_sequences_per_profile,
                        "Root randomly generated training sequence ",
                        m_random,
                        training_fastas[ 0 ],
                        training_root_ma
                      );
                      if( be_verbose ) {
                        cout << "Root alignment paths of training sequences are:" << endl;
                        //cout << training_root_ma << endl;
                        training_root_ma.toPairwiseStream( cout );
                      } // End if be_verbose
                      if( m_parameters.saveResultsToFile && m_parameters.saveTrueTrainingAlignments ) {
                        fs::path training_alignments_filename =
                          ( m_parameters.resultsFilePrefix +
                            run_unique_id +
                            "." + lexical_cast<string>( conservation_rate * 100 ) +
                            "." + lexical_cast<string>( profile_length ) +
                            "." + lexical_cast<string>( num_training_sequences_per_profile ) +
                            "." + lexical_cast<string>( m_parameters.numTestingSequencesPerProfile ) +
                                "." + lexical_cast<string>( expected_deletions_count ) +
                                "." + lexical_cast<string>( expected_insertions_count ) +
                                "." + lexical_cast<string>( expected_deletion_length_as_profile_length_fraction ) +
                                "." + lexical_cast<string>( expected_insertion_length_as_profile_length_fraction ) +
                            "." + lexical_cast<string>( true_profile_i ) +
                            ".true" +
                            ".root" +
                            m_parameters.trainingTrueAlignmentsFileSuffix );
                        std::ofstream training_alignments_stream( ( dirname / training_alignments_filename ).string().c_str() );
                        assert( training_alignments_stream.good() );
                        training_root_ma.toPairwiseStream( training_alignments_stream );
                        training_alignments_stream.close();
                      } // End if saveResultsToFile
                      training_fasta = training_fastas[ 0 ]; // All together now
                    
                      // Now draw the test sequences for the root profile.
                      // One fasta per profile
                      dp.drawSequences(
                        m_parameters,
                        true_root,
                        m_parameters.numTestingSequencesPerProfile,
                        "Root randomly generated test sequence ",
                        m_random,
                        testing_fastas[ 0 ],
                        testing_root_ma
                      );
                      if( be_verbose ) {
                        cout << "Root alignment paths of test sequences are:" << endl;
                        //cout << testing_root_ma << endl;
                        testing_root_ma.toPairwiseStream( cout );
                      } // End if be_verbose
                      if( m_parameters.saveResultsToFile && m_parameters.saveTrueTestingAlignments ) {
                        fs::path test_alignments_filename =
                          ( m_parameters.resultsFilePrefix +
                            run_unique_id +
                            "." + lexical_cast<string>( conservation_rate * 100 ) +
                            "." + lexical_cast<string>( profile_length ) +
                            "." + lexical_cast<string>( num_training_sequences_per_profile ) +
                            "." + lexical_cast<string>( m_parameters.numTestingSequencesPerProfile ) +
                                "." + lexical_cast<string>( expected_deletions_count ) +
                                "." + lexical_cast<string>( expected_insertions_count ) +
                                "." + lexical_cast<string>( expected_deletion_length_as_profile_length_fraction ) +
                                "." + lexical_cast<string>( expected_insertion_length_as_profile_length_fraction ) +
                            "." + lexical_cast<string>( true_profile_i ) +
                            ".true" +
                            ".root" +
                            m_parameters.trueTestingAlignmentsFileSuffix );
                        std::ofstream test_alignments_stream( ( dirname / test_alignments_filename ).string().c_str() );
                        assert( test_alignments_stream.good() );
                        testing_root_ma.toPairwiseStream( test_alignments_stream );
                        test_alignments_stream.close();
                      } // End if saveResultsToFile
                      testing_fasta = testing_fastas[ 0 ]; // All together now
                   
                      // Now draw the training and test sequences for the internal nodes.
                      for( which_profile = 1; which_profile < m_parameters.numProfiles; which_profile++ ) {
                        // TODO: REMOVE?
                        testProfileTree[ TEST_ID_true ].getProfileTreeInternalNode( which_profile ).ensurePositionsKnowTheirRoot();

                        dp.drawSequences(
                          m_parameters,
                          testProfileTree[ TEST_ID_true ].getProfileTreeInternalNode( which_profile ),
                          num_training_sequences_per_profile,
                          ( "Child " + boost::lexical_cast<std::string>( which_profile ) + " randomly generated " ),
                          m_random,
                          training_fastas[ which_profile ],
                          training_internal_node_ma
                        );
                        training_fasta += training_fastas[ which_profile ];
              
                        if( be_verbose ) {
                          cout << "Child " << which_profile << " alignment paths are:" << endl;
                          //cout << training_internal_node_ma << endl;
                          training_internal_node_ma.toPairwiseStream( cout );
                        } // End if be_verbose
                        if( m_parameters.saveResultsToFile && m_parameters.saveTrueTrainingAlignments ) {
                          fs::path training_alignments_filename =
                            ( m_parameters.resultsFilePrefix +
                              run_unique_id +
                              "." + lexical_cast<string>( conservation_rate * 100 ) +
                              "." + lexical_cast<string>( profile_length ) +
                              "." + lexical_cast<string>( num_training_sequences_per_profile ) +
                              "." + lexical_cast<string>( m_parameters.numTestingSequencesPerProfile ) +
                                "." + lexical_cast<string>( expected_deletions_count ) +
                                "." + lexical_cast<string>( expected_insertions_count ) +
                                "." + lexical_cast<string>( expected_deletion_length_as_profile_length_fraction ) +
                                "." + lexical_cast<string>( expected_insertion_length_as_profile_length_fraction ) +
                              "." + lexical_cast<string>( true_profile_i ) +
                              ".true" +
                              ".node" + lexical_cast<string>( which_profile ) +
                              m_parameters.trainingTrueAlignmentsFileSuffix );
                          std::ofstream training_alignments_stream( ( dirname / training_alignments_filename ).string().c_str() );
                          assert( training_alignments_stream.good() );
                          training_internal_node_ma.toPairwiseStream( training_alignments_stream );
                          training_alignments_stream.close();
                        } // End if saveResultsToFile
                    
                        dp.drawSequences(
                          m_parameters,
                          testProfileTree[ TEST_ID_true ].getProfileTreeInternalNode( which_profile ),
                          m_parameters.numTestingSequencesPerProfile,
                          ( "Child " + boost::lexical_cast<std::string>( which_profile ) + " randomly generated " ),
                          m_random,
                          testing_fastas[ which_profile ],
                          testing_internal_node_ma
                        );
                        testing_fasta += testing_fastas[ which_profile ];
                        if( be_verbose ) {
                          cout << "Child " << which_profile << " alignment paths are:" << endl;
                          //cout << testing_internal_node_ma << endl;
                          testing_internal_node_ma.toPairwiseStream( cout );
                        } // End if be_verbose
                        if( m_parameters.saveResultsToFile && m_parameters.saveTrueTestingAlignments ) {
                          fs::path test_alignments_filename =
                            ( m_parameters.resultsFilePrefix +
                              run_unique_id +
                              "." + lexical_cast<string>( conservation_rate * 100 ) +
                              "." + lexical_cast<string>( profile_length ) +
                              "." + lexical_cast<string>( num_training_sequences_per_profile ) +
                              "." + lexical_cast<string>( m_parameters.numTestingSequencesPerProfile ) +
                                "." + lexical_cast<string>( expected_deletions_count ) +
                                "." + lexical_cast<string>( expected_insertions_count ) +
                                "." + lexical_cast<string>( expected_deletion_length_as_profile_length_fraction ) +
                                "." + lexical_cast<string>( expected_insertion_length_as_profile_length_fraction ) +
                              "." + lexical_cast<string>( true_profile_i ) +
                              ".true" +
                              ".node" + lexical_cast<string>( which_profile ) +
                              m_parameters.trueTestingAlignmentsFileSuffix );
                          std::ofstream test_alignments_stream( ( dirname / test_alignments_filename ).string().c_str() );
                          assert( test_alignments_stream.good() );
                          testing_internal_node_ma.toPairwiseStream( test_alignments_stream );
                          test_alignments_stream.close();
                        } // End if saveResultsToFile
                      
                      } // End foreach child profile, draw the training and test sequences.
                    
                      if( be_verbose ) {
                        cout << "Randomly generated training sequences are:" << endl;
                        cout << training_fasta << endl;
                    
                        cout << "Randomly generated test sequences are:" << endl;
                        cout << testing_fasta << endl;
                      } // End if be_verbose
                      if( m_parameters.saveResultsToFile && m_parameters.saveTrainingSequences ) {
                        fs::path training_sequences_filename =
                          ( m_parameters.resultsFilePrefix +
                            run_unique_id +
                            "." + lexical_cast<string>( conservation_rate * 100 ) +
                            "." + lexical_cast<string>( profile_length ) +
                            "." + lexical_cast<string>( num_training_sequences_per_profile ) +
                            "." + lexical_cast<string>( m_parameters.numTestingSequencesPerProfile ) +
                                "." + lexical_cast<string>( expected_deletions_count ) +
                                "." + lexical_cast<string>( expected_insertions_count ) +
                                "." + lexical_cast<string>( expected_deletion_length_as_profile_length_fraction ) +
                                "." + lexical_cast<string>( expected_insertion_length_as_profile_length_fraction ) +
                            "." + lexical_cast<string>( true_profile_i ) +
                            m_parameters.trainingSequencesFileSuffix );
                        std::ofstream training_sequences_stream( ( dirname / training_sequences_filename ).string().c_str() );
                        assert( training_sequences_stream.good() );
                        training_sequences_stream << training_fasta;
                        training_sequences_stream.close();
                      } // End if saveResultsToFile
                      if( m_parameters.saveResultsToFile && m_parameters.saveTestingSequences ) {
                        fs::path test_sequences_filename =
                          ( m_parameters.resultsFilePrefix +
                            run_unique_id +
                            "." + lexical_cast<string>( conservation_rate * 100 ) +
                            "." + lexical_cast<string>( profile_length ) +
                            "." + lexical_cast<string>( num_training_sequences_per_profile ) +
                            "." + lexical_cast<string>( m_parameters.numTestingSequencesPerProfile ) +
                                "." + lexical_cast<string>( expected_deletions_count ) +
                                "." + lexical_cast<string>( expected_insertions_count ) +
                                "." + lexical_cast<string>( expected_deletion_length_as_profile_length_fraction ) +
                                "." + lexical_cast<string>( expected_insertion_length_as_profile_length_fraction ) +
                            "." + lexical_cast<string>( true_profile_i ) +
                            m_parameters.testingSequencesFileSuffix );
                        std::ofstream test_sequences_stream( ( dirname / test_sequences_filename ).string().c_str() );
                        assert( test_sequences_stream.good() );
                        test_sequences_stream << testing_fasta;
                        test_sequences_stream.close();
                      } // End if saveResultsToFile

                      // See if we need to reinitialize the training matrices
                      // OLD. See below.
                      //last_longest_training_sequence_length =
                      //  longest_training_sequence_length;
                      //for( uint32_t seq_i = 0;
                      //     seq_i < training_fasta.numSequences();
                      //     seq_i++ ) {
                      //  if( training_fasta[ seq_i ].length() >
                      //      longest_training_sequence_length ) {
                      //    longest_training_sequence_length =
                      //      training_fasta[ seq_i ].length();
                      //  }
                      //} // End foreach training fasta, see if it is the longest.
                      // TODO: Put back?
                      //if( longest_training_sequence_length >
                      //    last_longest_training_sequence_length ) {
                      //  training_forward_matrices.reinitialize(
                      //    profile_length,
                      //    training_fasta,
                      //    training_fasta.numSequences(),
                      //    longest_training_sequence_length
                      //  );
                      //} // End if longest_training_sequence_length > last_longest_training_sequence_length
                      // NEW:
                      training_forward_rows_1.reinitialize(
                        training_fasta,
                        training_fasta.size()
                      );
                      training_forward_rows_2.reinitialize(
                        training_fasta,
                        training_fasta.size()
                      );

                      // See if we need to reinitialize the testing matrices
                      // OLD.  See below.
                      //last_longest_testing_sequence_length =
                      //  longest_testing_sequence_length;
                      //for( uint32_t seq_i = 0;
                      //     seq_i < testing_fasta.numSequences();
                      //     seq_i++ ) {
                      //  if( testing_fasta[ seq_i ].length() >
                      //      longest_testing_sequence_length ) {
                      //    longest_testing_sequence_length =
                      //      testing_fasta[ seq_i ].length();
                      //  }
                      //} // End foreach testing fasta, see if it is the longest.
                      // TODO: Put back?
                      //if( longest_testing_sequence_length >
                      //    last_longest_testing_sequence_length ) {
                      //  testing_forward_matrices.reinitialize(
                      //    profile_length,
                      //    testing_fasta,
                      //    testing_fasta.numSequences(),
                      //    longest_testing_sequence_length
                      //  );
                      //} // End if longest_testing_sequence_length > last_longest_testing_sequence_length
                      // NEW:
                      testing_forward_rows_1.reinitialize(
                        testing_fasta,
                        testing_fasta.size()
                      );
                      testing_forward_rows_2.reinitialize(
                        testing_fasta,
                        testing_fasta.size()
                      );

                      // TODO: Start with the globals not correctly set.  If we start with even
                      // transition probs, the training fails.
                      //ProfileTreeType training_profile_tree =
                      //  testProfileTree[ TEST_ID_true ];

                        //ProfileTreeType( profile_length );
                    
                      //ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType> tree_trainer =
                      //  ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>(
                      //    training_fasta,
                      //    training_profile_tree
                      //  );
                      //tree_trainer.m_parameters = m_parameters;
                      //
                      //tree_trainer.train();
                    
                      // TODO: REMOVE.  TESTING
                      //cout << tree_trainer.m_parameters;
                    
                      // TODO: REMOVE
                      //if( true ) {
                      //  return 0;
                      //}
              
                      if( be_verbose ) {    
                        cout << "The true profile tree is:" << endl;
                        cout << testProfileTree[ TEST_ID_true ] << endl;
                      } // End if be_verbose
              
              
                      // Get "true profile" forward score for the training sequences.
                      // TODO: Put back -- or make new tree scorer..
                        //dp.forward_score( 
                        //  m_parameters,
                        //  testProfileTree[ TEST_ID_true ],
                        //  training_fastas,
                        //  training_forward_matrices
                        //);
                        dp.forward_score(
                          m_parameters,
                          false, // do not use viterbi
                          *testProfileTree[ TEST_ID_true ].getProfileTreeRoot(),
                          training_fasta,
                          training_fasta.size(),
                          training_forward_rows_1,
                          training_forward_rows_2,
                          testScore_training_forward[ TEST_ID_true ]
                        );
                        if( be_verbose ) {
                          cout << "The \"true profile\" forward score for the training sequences is " << testScore_training_forward[ TEST_ID_true ] << endl;
                        } // End if be_verbose
                        
                      for( uint32_t starting_profile_i = 0; starting_profile_i < m_parameters.numStartingProfiles; starting_profile_i++ ) {

                        // TODO: REMOVE
                        //cout << "STARTING PROFILE index is " << starting_profile_i << endl;
              
                        // The starting profile.  To get the length from
                        // the true profile, reinitialize.
                        testProfileTree[ TEST_ID_starting ].reinitialize( testProfileTree[ TEST_ID_true ] );
              
                        // Start over
                        //testProfileTree[ TEST_ID_starting ].evenPositions();

                        if( m_parameters.trainProfilePositions == false ) {
                          testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->copyPositions( *( testProfileTree[ TEST_ID_true ].getProfileTreeRoot() ) );
                        } else {
                          if(
                             ( !m_parameters.alsoStartWithEvenPositions || ( starting_profile_i > 0 ) ) &&
                            m_parameters.startWithPositionsDrawnFromPrior &&
                             ( !m_parameters.startWithUniformPositions || ( starting_profile_i < ( ( m_parameters.numStartingProfiles / 2 ) + ( m_parameters.alsoStartWithEvenPositions ? 1 : 0 ) ) ) ) // If startWithUniformPositions is true ALSO, then only start the first half from the prior, and the last half with uniform (excluding the 0th, if alsoStartWithEvenPositions).
                          ) {
                            // Draw position params from the prior
                            testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->dirichletMixturePositions(
                              matchEmissionPrior,
                              m_random
                            );
                          } else if(
                             ( !m_parameters.alsoStartWithEvenPositions || ( starting_profile_i > 0 ) ) &&
                             m_parameters.startWithUniformPositions
                          ) {
                            // Set the root of the starting profile tree to
                            // uniformly-distributed values.
                            testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->uniformPositions( m_random );
                          } else {
                            // even
                            testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->evenPositions();
                          } // End if draw from priors .. else if startWithUniformPositions .. else ..
                        } // End if !trainProfilePositions .. else ..

                        // Globals either start as even (default) or are random
                        // (if startWithUniformGlobals is true), in which case
                        // they are either drawn from the prior 
                        // or are set to be uniform up to a maximum --
                        // unless trainProfileGlobals is false, in which case
                        // they start as the true globals (and stay there).
                        if( m_parameters.trainProfileGlobals == false ) {
                          testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->copyExceptPositions( *( testProfileTree[ TEST_ID_true ].getProfileTreeRoot() ) );
                        } else {
                          RootType & starting_root =
                            ( *testProfileTree[ TEST_ID_starting ].getProfileTreeRoot() );
                            
                          // Note that we don't need to set the values of the
                          // internal nodes of the starting profile tree, since
                          // they will be set to their respective parent's
                          // values before they are trained.
                          if(
                            m_parameters.startWithGlobalsDrawnFromPrior
                          ) {
                            // Draw globals from the prior
                            starting_root.dirichletMixtureExceptPositions(
                              globalPrior,
                              m_random
                            );
                          } else if( m_parameters.startWithUniformGlobals ) {
                            // Set the globals to be uniform within the range
                            // 0 to min( params.startWithPriorGlobals_minX, (
                            // startWithUniformGlobals_scalar * the truth) ).
                            // ALAS i have to do this manually!
                            ProbabilityType unif_max;
                            
                            unif_max =
                              min( ( ProbabilityType )m_parameters.startWithUniformGlobals_maxNtoN, ( ( ProbabilityType )( m_parameters.startWithUniformGlobals_scalar ) * true_root[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] ) );
                            starting_root[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] =
                              ( ( ProbabilityType )( m_random.nextUniform() ) * unif_max );
                            starting_root[ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ] =
                              ( ProbabilityType )( 1.0 ) -
                              starting_root[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ];
                            
                            unif_max =
                              min( ( ProbabilityType )m_parameters.startWithUniformGlobals_maxBtoD, ( ( ProbabilityType )( m_parameters.startWithUniformGlobals_scalar ) * true_root[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ] ) );
                            starting_root[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ] =
                              ( ( ProbabilityType )( m_random.nextUniform() ) * unif_max );
                            starting_root[ Transition::fromBegin ][ TransitionFromBegin::toMatch ] =
                              ( ProbabilityType )( 1.0 ) -
                              starting_root[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ];
                            
                            // There is a potential problem: if we're not
                            // careful we get such large M->I and M->D probs
                            // that there's nothing left for M->M.
                            assert( m_parameters.startWithUniformGlobals_maxMtoI + m_parameters.startWithUniformGlobals_maxMtoD < 1.0 );
#ifdef USE_DEL_IN_DEL_OUT
                            // TODO: Create m_parameters.startWithUniformGlobals_maxMtoW
                            assert( m_parameters.startWithUniformGlobals_maxMtoI + m_parameters.startWithUniformGlobals_maxMtoD + .5 < 1.0 );
#endif // USE_DEL_IN_DEL_OUT
                            unif_max =
                              min( ( ProbabilityType )m_parameters.startWithUniformGlobals_maxMtoI, ( ( ProbabilityType )( m_parameters.startWithUniformGlobals_scalar ) * true_root[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] ) );
                            starting_root[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] =
                              ( ( ProbabilityType )( m_random.nextUniform() ) * unif_max );
                            
                            unif_max =
                              min( ( ProbabilityType )m_parameters.startWithUniformGlobals_maxMtoD, ( ( ProbabilityType )( m_parameters.startWithUniformGlobals_scalar ) * true_root[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] ) );
                            starting_root[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] =
                              ( ( ProbabilityType )( m_random.nextUniform() ) * unif_max );
#ifdef USE_DEL_IN_DEL_OUT
                            // TODO: ERE I AM!!!
#else
                            starting_root[ Transition::fromMatch ][ TransitionFromMatch::toMatch ] =
                              ( ProbabilityType )( 1.0 ) -
                              (
                               starting_root[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] +
                               starting_root[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ]
                              );
#endif // USE_DEL_IN_DEL_OUT .. else ..                            
                            unif_max =
                              min( ( ProbabilityType )m_parameters.startWithUniformGlobals_maxItoI, ( ( ProbabilityType )( m_parameters.startWithUniformGlobals_scalar ) * true_root[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] ) );
                            starting_root[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] =
                              ( ( ProbabilityType )( m_random.nextUniform() ) * unif_max );
                            starting_root[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ] =
                              ( ProbabilityType )( 1.0 ) -
                              starting_root[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ];
                            
                            unif_max =
                              min( ( ProbabilityType )m_parameters.startWithUniformGlobals_maxDtoD, ( ( ProbabilityType )( m_parameters.startWithUniformGlobals_scalar ) * true_root[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ] ) );
                            starting_root[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ] =
                              ( ( ProbabilityType )( m_random.nextUniform() ) * unif_max );
                            starting_root[ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ] =
                              ( ProbabilityType )( 1.0 ) -
                              starting_root[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ];
                            
                            // We don't allow looping, so we fix the end->loop
                            // prob at 0
                            // For now we don't use the End distribution ..
                            //starting_root[ Transition::fromEnd ][ TransitionFromEnd::toPostAlign ] = 1.0;
                            //starting_root[ Transition::fromEnd ][ TransitionFromEnd::toLoop ] = 0.0;
                            
                            unif_max =
                              min( ( ProbabilityType )m_parameters.startWithUniformGlobals_maxCtoC, ( ( ProbabilityType )( m_parameters.startWithUniformGlobals_scalar ) * true_root[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] ) );
                            starting_root[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] =
                              ( ( ProbabilityType )( m_random.nextUniform() ) * unif_max );
                            
                            starting_root[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] =
                              ( ProbabilityType )( 1.0 ) -
                              starting_root[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ];
                          } else {
                            // even
                            testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->evenExceptPositions();
                          } // End if startWithGlobalsDrawnFromPrior .. else if startWithUniformGlobals .. else ..

                          // TODO: Incorporate into parameters.

                          // When doing lengthadjust, it is best to start with
                          // the indel opens equally probable.
                          if( m_parameters.testLengthadjust ) {
                            ProbabilityType average_indel_open = 1.0;
                            average_indel_open -= starting_root[ Transition::fromMatch ][ TransitionFromMatch::toMatch ];
                            average_indel_open /= 2;
                            starting_root[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] =
                              average_indel_open;
                            starting_root[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] =
                              average_indel_open;
                          } // End if startWithEqualIndelOpens
                        } // End if trainProfileGlobals
              
                        if( be_verbose ) {    
                          cout << "The starting profile tree is:" << endl;
                          cout << testProfileTree[ TEST_ID_starting ] << endl;
                        } // End if be_verbose
              
                        // Get "starting profile" forward score for the training sequences.
                        // TODO: Put back..
                        //  dp.forward_score( 
                        //    m_parameters,
                        //    testProfileTree[ TEST_ID_starting ],
                        //    training_fastas,
                        //    training_forward_matrices
                        //   );
                        dp.forward_score(
                          m_parameters,
                          false, // do not use viterbi
                          *testProfileTree[ TEST_ID_starting ].getProfileTreeRoot(),
                          training_fasta,
                          training_fasta.size(),
                          training_forward_rows_1,
                          training_forward_rows_2,
                          testScore_training_forward[ TEST_ID_starting ]
                        );
                        if( be_verbose ) {
                          cout << "The \"starting profile\" forward score for the training sequences is " << testScore_training_forward[ TEST_ID_starting ] << endl;
                        } // End if be_verbose
              
                        // TODO: REMOVE
                        //cout << "=======YO=======" << endl;
                        //cout << "last_test_id is " << last_test_id << endl;

                        if( m_parameters.saveResultsToFile && m_parameters.saveStartingProfiles ) {
                          fs::path starting_root_filename =
                            ( m_parameters.resultsFilePrefix +
                              run_unique_id +
                              "." + lexical_cast<string>( conservation_rate * 100 ) +
                              "." + lexical_cast<string>( profile_length ) +
                              "." + lexical_cast<string>( num_training_sequences_per_profile ) +
                              "." + lexical_cast<string>( m_parameters.numTestingSequencesPerProfile ) +
                                "." + lexical_cast<string>( expected_deletions_count ) +
                                "." + lexical_cast<string>( expected_insertions_count ) +
                                "." + lexical_cast<string>( expected_deletion_length_as_profile_length_fraction ) +
                                "." + lexical_cast<string>( expected_insertion_length_as_profile_length_fraction ) +
                              "." + lexical_cast<string>( true_profile_i ) +
                              "." + lexical_cast<string>( starting_profile_i ) +
                              m_parameters.startingProfileTreeFileSuffix );
                          writeXML( testProfileTree[ TEST_ID_starting ], ( dirname / starting_root_filename ).string().c_str() );
                        } // End if saveResultsToFile

                        // TODO: REMOVE
                        //cout << "=======YAY=======" << endl;
                        //cout << "last_test_id is " << last_test_id << endl;

                        /// TODO: Move these initializations out of the loop:

                        // Note that we must set the m_profile value before
                        // training...
                        ProfileGibbs<RootType, ScoreType, MatrixValueType, SequenceResidueType> sampler(
                          &m_random,
                          testProfileTree[ TEST_ID_starting ].getProfileTreeRoot(), // not used except to choose the right template type
                          training_fasta
                        );

                        // Note that we must set the m_profileTree value before
                        // training...
                        ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType> trainer(
                          training_fasta,
                          NULL
                        );

                        for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                          // SymmeterizedKullbackLeibler (SKL) distance/divergence stuff
                          if(
                            (
                              ( test_id == TEST_ID_starting ) ||
                              ( test_id == TEST_ID_true )
                            ) &&
                            (
                              m_parameters.calculateSymmeterizedKullbackLeiblerDistancesToTrue ||
                              m_parameters.calculateSymmeterizedKullbackLeiblerDistancesToStarting
                            )
                          ) {
                            // TODO: Support the whole tree.
                            if(
                              ( test_id == TEST_ID_true ) &&
                              ( m_parameters.numProfiles > 1 )
                            ) {
                              cout << "WARNING: Calculating SKL distances using only the root profile, not the whole tree.  TODO: implement it for the whole tree." << endl;
                            }
                            selfEntropyPositions[ test_id ] =
                              testProfileTree[ test_id ].getProfileTreeRoot()->crossEntropyPositions( *testProfileTree[ test_id ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 );
                            if( be_verbose ) {
                              cout << "The " << tests[ test_id ].name << " positions (self-)entropy is " << selfEntropyPositions[ test_id ] << endl;
                            }
                            selfEntropyExceptPositions[ test_id ] =
                              testProfileTree[ test_id ].getProfileTreeRoot()->crossEntropyExceptPositions( *testProfileTree[ test_id ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 );
                            if( be_verbose ) {
                              cout << "The " << tests[ test_id ].name << " non-position distributions (self-)entropy is " << selfEntropyExceptPositions[ test_id ] << endl;
                            }
                          } // End calculating self entropy (starting or true profiles only)
                          if(
                            ( test_id == TEST_ID_starting ) &&
                            (
                              m_parameters.calculateSymmeterizedKullbackLeiblerDistancesToTrue ||
                              m_parameters.calculateSymmeterizedKullbackLeiblerDistancesToStarting
                            ) 
                          ) {
                            //distanceSKLPositions_starting[ test_id ] = 0.0;
                            //distanceSKLExceptPositions_starting[ test_id ] = 0.0;
                            // *now* calculate the distance to true
                            // (assumption here is that TEST_ID_starting >
                            // TEST_ID_true)
                            if(
                              testProfileTree[ TEST_ID_true ].getProfileTreeRoot()->length() !=
                              testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->length()
                            ) {
                              distanceSKLPositions_starting[ TEST_ID_true ] =
                                numeric_limits<double>::max();
                            } else {
                              distanceSKLPositions_starting[ TEST_ID_true ] =
                                ( 
                                  ( testProfileTree[ TEST_ID_true ].getProfileTreeRoot()->crossEntropyPositions( *testProfileTree[ TEST_ID_starting ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) - selfEntropyPositions[ TEST_ID_true ] ) +
                                  ( testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->crossEntropyPositions( *testProfileTree[ TEST_ID_true ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) - selfEntropyPositions[ TEST_ID_starting ] )
                                );
                            } // End if lengths are different .. else ..
                            distanceSKLPositions_true[ TEST_ID_starting ] =
                              distanceSKLPositions_starting[ TEST_ID_true ];
                            if( be_verbose ) {
                              if(
                                testProfileTree[ TEST_ID_true ].getProfileTreeRoot()->length() ==
                                testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->length()
                              ) {
                                cout << "The true profile (positions) cross-entropy with the starting profile is " << testProfileTree[ TEST_ID_true ].getProfileTreeRoot()->crossEntropyPositions( *testProfileTree[ TEST_ID_starting ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) << endl;
                                cout << "The starting profile (positions) cross-entropy with the true profile is " << testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->crossEntropyPositions( *testProfileTree[ TEST_ID_true ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) << endl;
                              } // End if lengths are the same
                              cout << "The starting profile (positions) symmeterized Kullback-Leibler divergence with the true profile is " << distanceSKLPositions_starting[ TEST_ID_true ]  << endl;
                            } // End if be_verbose
                            if(
                              testProfileTree[ TEST_ID_true ].getProfileTreeRoot()->length() !=
                              testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->length()
                            ) {
                            distanceSKLExceptPositions_starting[ TEST_ID_true ] =
                              numeric_limits<double>::max();
                            } else {
                              distanceSKLExceptPositions_starting[ TEST_ID_true ] =
                                (
                                  ( testProfileTree[ TEST_ID_true ].getProfileTreeRoot()->crossEntropyExceptPositions( *testProfileTree[ TEST_ID_starting ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) - selfEntropyExceptPositions[ TEST_ID_true ] ) +
                                  ( testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->crossEntropyExceptPositions( *testProfileTree[ TEST_ID_true ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) - selfEntropyExceptPositions[ TEST_ID_starting ] )
                                );
                            } // End if lengths are different .. else ..
                            distanceSKLExceptPositions_true[ TEST_ID_starting ] =
                              distanceSKLExceptPositions_starting[ TEST_ID_true ];
                            if( be_verbose ) {
                              if(
                                testProfileTree[ TEST_ID_true ].getProfileTreeRoot()->length() ==
                                testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->length()
                              ) {
                                cout << "The true profile (non-positions) cross-entropy with the starting profile is " << testProfileTree[ TEST_ID_true ].getProfileTreeRoot()->crossEntropyExceptPositions( *testProfileTree[ TEST_ID_starting ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) << endl;
                                cout << "The starting profile (non-positions) cross-entropy with the true profile is " << testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->crossEntropyExceptPositions( *testProfileTree[ TEST_ID_true ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) << endl;
                              } // End if lengths are the same
                              cout << "The starting profile (non-positions) symmeterized Kullback-Leibler divergence with the true profile is " << distanceSKLExceptPositions_starting[ TEST_ID_true ]  << endl;
                            } // End if be_verbose
                          } // End if TEST_ID_starting && ( calculateSymmeterizedKullbackLeiblerDistancesToStarting || calculateSymmeterizedKullbackLeiblerDistancesToTrue )

                          if( ( test_id == TEST_ID_true ) ||
                              ( test_id == TEST_ID_starting ) ) {
                            // testProfileSequences[ test_id ] = training_fasta
                            //for( uint32_t node_i = 0; node_i < training_fastas.size(); node_i++ ) {
                            //  testProfileSequences[ test_id ][ node_i ] =
                            //    training_fastas[ node_i ];
                            //}
                            continue; // Don't train true and starting profiles.
                          }
              
                          if( !tests[ test_id ].isRun ) {
                            continue;
                          }
              
                          // TODO: REMOVE
                          //cout << "Running test id " << test_id << endl;
              
                          if( tests[ test_id ].startingPositionsTest == NULL ) {
                            // Use the "starting" test
                            testProfileTree[ test_id ] = testProfileTree[ TEST_ID_starting ];
                          } else {
                            if( tests[ test_id ].startingPositionsTest == &tests[ test_id ] ) {
                              // Use the "true" test
                              testProfileTree[ test_id ] = testProfileTree[ TEST_ID_true ];
                            } else {
                              testProfileTree[ test_id ] =
                                *( testProfileTreePtrMap[ tests[ test_id ].startingPositionsTest ] );
                            }
                          } // End if tests[ test_id ].startingPositionsTest == NULL .. else ..
              
                          if( tests[ test_id ].startingGlobalsTest == NULL ) {
                            // Use the "starting" test
                            if( tests[ test_id ].startingPositionsTest != NULL ) {
                              testProfileTree[ test_id ].getProfileTreeRoot()->copyExceptPositions( ( *testProfileTree[ TEST_ID_starting ].getProfileTreeRoot() ) );
                            }
                          } else {
                            if( tests[ test_id ].startingGlobalsTest == &tests[ test_id ] ) {
                              // Use the "true" test
                              if( tests[ test_id ].startingPositionsTest != &tests[ test_id ] ) {
                                testProfileTree[ test_id ].getProfileTreeRoot()->copyExceptPositions( ( *testProfileTree[ TEST_ID_true ].getProfileTreeRoot() ) );
                              }
                            } else {
                              if( tests[ test_id ].startingGlobalsTest !=
                                  tests[ test_id ].startingPositionsTest ) {
                                testProfileTree[ test_id ].getProfileTreeRoot()->copyExceptPositions(
                                  *( ( testProfileTreePtrMap[ tests[ test_id ].startingGlobalsTest ] )->getProfileTreeRoot() )
                                );
                              }
                            }
                          } // End if tests[ test_id ].startingGlobalsTest == NULL .. else ..
                          if( tests[ test_id ].isGibbs ) {

                            if( !m_parameters.sampleProfilePositions ) {
                              // Make sure that the profile starts with the
                              // right position values, since we're expecting
                              // them to be the true values.
                              testProfileTree[ test_id ].getProfileTreeRoot()->copyPositions(
                                *( testProfileTree[ TEST_ID_true ].getProfileTreeRoot() )
                              );
                            }
                            
                            sampler.m_profile =
                              testProfileTree[ test_id ].getProfileTreeRoot();
                            sampler.m_parameters = m_parameters;
                            
                            tests[ test_id ].parametersModifier.applyModifications( sampler.m_parameters );
                            
                            // TODO: REMOVE
                            //cout << "sampler verbosity level is " << sampler.m_parameters.verbosity << endl;
                            //if( tests[ test_id ].parametersModifier.isModified_verbosity ) {
                            //  cout << "it should be modified to " << tests[ test_id ].parametersModifier.parameters.verbosity << endl;
                            //} else {
                            //  cout << "it should be unmodified (" << m_parameters.verbosity << ")" << endl;
                            //}
                            
                            if( be_verbose ) {
                              cout << "Now (before sampling the " << tests[ test_id ].name << " profile), the profile tree is:" << endl;
                              cout << testProfileTree[ test_id ] << endl;
                            } // End if be_verbose

                            testScore_training_forward[ test_id ] =
                              sampler.sample();
                            testIters[ test_id ] = sampler.m_totalIterations;

                            // The sampler puts the best profile it found in
                            // m_samplingProfile.  Note that this might not be the mode or the overall mean (it could be the mean of one of the chains).
                            testProfileTree[ test_id ].getProfileTreeRoot()->copyFrom( sampler.m_samplingProfile );

                            // Also separately store the overall mean and mode
                            if( m_parameters.reportGibbsMean ) {
                              testScore_training_mean_forward[ test_id ] =
                                sampler.m_averageProfileScore;
                              testProfileTree_mean[ test_id ].getProfileTreeRoot()->copyFrom( sampler.m_averageProfile );
                            }
                            if( m_parameters.reportGibbsMode ) {
                              testScore_training_mode_forward[ test_id ] =
                                sampler.m_bestProfileScore;
                              testProfileTree_mode[ test_id ].getProfileTreeRoot()->copyFrom( sampler.m_bestProfile );
                            }
                            if( be_verbose ) {
                              cout << "Now (after training/sampling the " << tests[ test_id ].name << " profile), the score is " << testScore_training_forward[ test_id ] << ", and the profile tree is:" << endl;
                              cout << testProfileTree[ test_id ] << endl;
                            } // End if be_verbose

                          } else { // if isGibbs .. else ..
                            
                            trainer.m_profileTree = &testProfileTree[ test_id ];
                            trainer.m_parameters = m_parameters;
                            
                            tests[ test_id ].parametersModifier.applyModifications( trainer.m_parameters );
                            
                            // TODO: REMOVE
                            //cout << "trainer verbosity level is " << trainer.m_parameters.verbosity << endl;
                            //if( tests[ test_id ].parametersModifier.isModified_verbosity ) {
                            //  cout << "it should be modified to " << tests[ test_id ].parametersModifier.parameters.verbosity << endl;
                            //} else {
                            //  cout << "it should be unmodified (" << m_parameters.verbosity << ")" << endl;
                            //}
                            
                            if( be_verbose && 1 ) {
                              cout << "Now (before training the " << tests[ test_id ].name << " profile), the profile tree is:" << endl;
                              cout << testProfileTree[ test_id ] << endl;
                            } // End if be_verbose
                            if( be_verbose && false ) {
                              cout << "Now (before training the " << tests[ test_id ].name << " profile), the profile trainer parameters are:" << endl;
                              cout << trainer.m_parameters << endl;
                            } // End if be_verbose
                            start_time = clock();
                            testScore_training_forward[ test_id ] = trainer.train();
                            end_time = clock();
                            testCPUtime[ test_id ] = ( static_cast<double>( end_time - start_time ) ) / CLOCKS_PER_SEC;
                            testIters[ test_id ] = trainer.m_totalIterations;

                            if( be_verbose ) {
                              cout << "Trained profile length is " << testProfileTree[ test_id ].getProfileTreeRoot()->length() << endl;
                            }

                            // If lengths have changed, we may need to resize
                            // the forward matrices.  They have to be greater
                            // than or equal to the right size...
                            // TODO: Put back..
                            //if( training_forward_matrices.size() < ( testProfileTree[ test_id ].getProfileTreeRoot()->length() + 1 ) ) {
                            //  training_forward_matrices.reinitialize(
                            //    testProfileTree[ test_id ].getProfileTreeRoot()->length(),
                            //    training_fasta,
                            //    training_fasta.numSequences(),
                            //    longest_training_sequence_length
                            //  );
                            //} // End if we need to reinitialize the training forward matrices (due to a change in profile length)
                            // TODO: Put back..
                            //if( testing_forward_matrices.size() < ( testProfileTree[ test_id ].getProfileTreeRoot()->length() + 1 ) ) {
                            //  testing_forward_matrices.reinitialize(
                            //    testProfileTree[ test_id ].getProfileTreeRoot()->length(),
                            //    testing_fasta,
                            //    testing_fasta.numSequences(),
                            //    longest_testing_sequence_length
                            //  );
                            //} // End if we need to reinitialize the testing forward matrices (due to a change in profile length)

                            if( be_verbose ) {
                              cout << "Now (after training the " << tests[ test_id ].name << " profile), the score is " << testScore_training_forward[ test_id ] << ", and the profile tree is:" << endl;
                              cout << testProfileTree[ test_id ] << endl;
                            } // End if be_verbose
                          } // End if isGibbs .. else ..

                          if( m_parameters.saveResultsToFile && m_parameters.saveTestProfiles ) {
                            fs::path test_root_filename =
                              ( m_parameters.resultsFilePrefix +
                                run_unique_id +
                                "." + lexical_cast<string>( conservation_rate * 100 ) +
                                "." + lexical_cast<string>( profile_length ) +
                                "." + lexical_cast<string>( num_training_sequences_per_profile ) +
                                "." + lexical_cast<string>( m_parameters.numTestingSequencesPerProfile ) +
                                "." + lexical_cast<string>( expected_deletions_count ) +
                                "." + lexical_cast<string>( expected_insertions_count ) +
                                "." + lexical_cast<string>( expected_deletion_length_as_profile_length_fraction ) +
                                "." + lexical_cast<string>( expected_insertion_length_as_profile_length_fraction ) +
                                "." + lexical_cast<string>( true_profile_i ) +
                                "." + lexical_cast<string>( starting_profile_i ) +
                                "." + tests[ test_id ].name +
                                m_parameters.testProfileTreeFileSuffix );
                            writeXML( testProfileTree[ test_id ], ( dirname / test_root_filename ).string().c_str() );
                          } // End if saveResultsToFile
                          //testProfileSequences[ test_id ] =
                          //  trainer.m_profile_sequences;

                          // SymmeterizedKullbackLeibler (SKL) distance/divergence stuff
                          if(
                            tests[ test_id ].isRun &&
                            (
                              m_parameters.calculateSymmeterizedKullbackLeiblerDistancesToTrue ||
                              m_parameters.calculateSymmeterizedKullbackLeiblerDistancesToStarting
                            )
                          ) {
                            selfEntropyPositions[ test_id ] =
                              testProfileTree[ test_id ].getProfileTreeRoot()->crossEntropyPositions( *testProfileTree[ test_id ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 );
                            if( be_verbose ) {
                              cout << "The " << tests[ test_id ].name << " positions (self-)entropy is " << selfEntropyPositions[ test_id ] << endl;
                            }
                            selfEntropyExceptPositions[ test_id ] =
                              testProfileTree[ test_id ].getProfileTreeRoot()->crossEntropyExceptPositions( *testProfileTree[ test_id ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 );
                            if( be_verbose ) {
                              cout << "The " << tests[ test_id ].name << " non-position distributions (self-)entropy is " << selfEntropyExceptPositions[ test_id ] << endl;
                            }
                          } // End calculating self entropy

                          if( m_parameters.calculateSymmeterizedKullbackLeiblerDistancesToTrue ) {
                            if(
                              testProfileTree[ TEST_ID_true ].getProfileTreeRoot()->length() !=
                              testProfileTree[ test_id ].getProfileTreeRoot()->length()
                            ) {
                            distanceSKLPositions_true[ test_id ] =
                              numeric_limits<double>::max();
                            } else {
                              distanceSKLPositions_true[ test_id ] =
                                (
                                  ( testProfileTree[ test_id ].getProfileTreeRoot()->crossEntropyPositions( *testProfileTree[ TEST_ID_true ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) - selfEntropyPositions[ test_id ] ) +
                                  ( testProfileTree[ TEST_ID_true ].getProfileTreeRoot()->crossEntropyPositions( *testProfileTree[ test_id ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) - selfEntropyPositions[ TEST_ID_true ] )
                                );
                            } // End if lengths are different .. else ..
                            if( be_verbose ) {
                              if(
                                testProfileTree[ TEST_ID_true ].getProfileTreeRoot()->length() ==
                                testProfileTree[ test_id ].getProfileTreeRoot()->length()
                              ) {
                                cout << "The " << tests[ test_id ].name << " (positions) cross-entropy with the true profile is " << testProfileTree[ test_id ].getProfileTreeRoot()->crossEntropyPositions( *testProfileTree[ TEST_ID_true ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) << endl;
                                cout << "The true profile (positions) cross-entropy with the " << tests[ test_id ].name << " profile is " << testProfileTree[ TEST_ID_true ].getProfileTreeRoot()->crossEntropyPositions( *testProfileTree[ test_id ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) << endl;
                              } // End if lengths are the same
                              cout << "The " << tests[ test_id ].name << " profile (positions) symmeterized Kullback-Leibler divergence with the true profile is " << distanceSKLPositions_true[ test_id ]  << endl;
                            } // End if be_verbose
                            if(
                              testProfileTree[ TEST_ID_true ].getProfileTreeRoot()->length() !=
                              testProfileTree[ test_id ].getProfileTreeRoot()->length()
                            ) {
                              distanceSKLExceptPositions_true[ test_id ] =
                                numeric_limits<double>::max();
                            } else {
                              distanceSKLExceptPositions_true[ test_id ] =
                                (
                                  ( testProfileTree[ test_id ].getProfileTreeRoot()->crossEntropyExceptPositions( *testProfileTree[ TEST_ID_true ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) - selfEntropyExceptPositions[ test_id ] ) +
                                  ( testProfileTree[ TEST_ID_true ].getProfileTreeRoot()->crossEntropyExceptPositions( *testProfileTree[ test_id ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) - selfEntropyExceptPositions[ TEST_ID_true ] )
                                );
                            } // End if lengths are different .. else ..
                            if( be_verbose ) {
                              if(
                                testProfileTree[ TEST_ID_true ].getProfileTreeRoot()->length() ==
                                testProfileTree[ test_id ].getProfileTreeRoot()->length()
                              ) {
                                cout << "The " << tests[ test_id ].name << " (non-positions) cross-entropy with the true profile is " << testProfileTree[ test_id ].getProfileTreeRoot()->crossEntropyExceptPositions( *testProfileTree[ TEST_ID_true ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) << endl;
                                cout << "The true profile (non-positions) cross-entropy with the " << tests[ test_id ].name << " profile is " << testProfileTree[ TEST_ID_true ].getProfileTreeRoot()->crossEntropyExceptPositions( *testProfileTree[ test_id ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) << endl;
                              } // End if lengths are the same
                              cout << "The " << tests[ test_id ].name << " profile (non-positions) symmeterized Kullback-Leibler divergence with the true profile is " << distanceSKLExceptPositions_true[ test_id ]  << endl;
                            }
                          } // End if calculateSymmeterizedKullbackLeiblerDistancesToTrue
                          if( m_parameters.calculateSymmeterizedKullbackLeiblerDistancesToStarting ) {
                            if(
                              testProfileTree[ test_id ].getProfileTreeRoot()->length() !=
                              testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->length()
                            ) {
                              distanceSKLPositions_starting[ test_id ] =
                                numeric_limits<double>::max();
                            } else {
                              distanceSKLPositions_starting[ test_id ] =
                                (
                                  ( testProfileTree[ test_id ].getProfileTreeRoot()->crossEntropyPositions( *testProfileTree[ TEST_ID_starting ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) - selfEntropyPositions[ test_id ] ) +
                                  ( testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->crossEntropyPositions( *testProfileTree[ test_id ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) - selfEntropyPositions[ TEST_ID_starting ] )
                                );
                            } // End if lengths are different .. else ..
                            if( be_verbose ) {
                              if(
                                testProfileTree[ test_id ].getProfileTreeRoot()->length() ==
                                testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->length()
                              ) {
                                cout << "The " << tests[ test_id ].name << " (positions) cross-entropy with the starting profile is " << testProfileTree[ test_id ].getProfileTreeRoot()->crossEntropyPositions( *testProfileTree[ TEST_ID_starting ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) << endl;
                                cout << "The starting profile (positions) cross-entropy with the " << tests[ test_id ].name << " profile is " << testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->crossEntropyPositions( *testProfileTree[ test_id ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) << endl;
                              } // End if lengths are the same
                              cout << "The " << tests[ test_id ].name << " profile (positions) symmeterized Kullback-Leibler divergence with the starting profile is " << distanceSKLPositions_starting[ test_id ]  << endl;
                            } // End if be_verbose
                            if(
                              testProfileTree[ test_id ].getProfileTreeRoot()->length() !=
                              testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->length()
                            ) {
                              distanceSKLExceptPositions_starting[ test_id ] =
                                numeric_limits<double>::max();
                            } else {
                              distanceSKLExceptPositions_starting[ test_id ] =
                                (
                                  ( testProfileTree[ test_id ].getProfileTreeRoot()->crossEntropyExceptPositions( *testProfileTree[ TEST_ID_starting ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) - selfEntropyExceptPositions[ test_id ] ) +
                                  ( testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->crossEntropyExceptPositions( *testProfileTree[ test_id ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) - selfEntropyExceptPositions[ TEST_ID_starting ] )
                                );
                            } // End if lengths are different .. else .. 
                            if( be_verbose ) {
                              if(
                                testProfileTree[ test_id ].getProfileTreeRoot()->length() ==
                                testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->length()
                              ) {
                                cout << "The " << tests[ test_id ].name << " (non-positions) cross-entropy with the starting profile is " << testProfileTree[ test_id ].getProfileTreeRoot()->crossEntropyExceptPositions( *testProfileTree[ TEST_ID_starting ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) << endl;
                                cout << "The starting profile (non-positions) cross-entropy with the " << tests[ test_id ].name << " profile is " << testProfileTree[ TEST_ID_starting ].getProfileTreeRoot()->crossEntropyExceptPositions( *testProfileTree[ test_id ].getProfileTreeRoot(), ( ProfileTreeRoot<ResidueType, double> const * )0 ) << endl;
                              } // End if lengths are the same
                              cout << "The " << tests[ test_id ].name << " profile (non-positions) symmeterized Kullback-Leibler divergence with the starting profile is " << distanceSKLExceptPositions_starting[ test_id ]  << endl;
                            } // End if be_verbose
                          } // End if calculateSymmeterizedKullbackLeiblerDistancesToStarting

                          if( m_parameters.calculateProfileProfileAlignments ) {
                            profile_profile_alignment_cost =
                              dp.profileProfile_align_SKL(
                                m_parameters,
                                *testProfileTree[ TEST_ID_true ].getProfileTreeRoot(),
                                *testProfileTree[ test_id ].getProfileTreeRoot(),
                                m_parameters.profileProfileIndelOpenCost,
                                m_parameters.profileProfileIndelExtensionCost,
                                profile_profile_alignment//,
                                // TODO: REMOVE?  Yes, it makes everything align as indel-indel 'match'es.  Bogus.
                                //( 2.0 * m_parameters.profileProfileIndelOpenCost )
                              );
                            uint32_t num_opens = 0;
                            uint32_t num_extensions = 0;
                            numProfileProfileAlignedPositions_true[ test_id ] = 0;
                            if( profile_profile_alignment[ 0 ] > 0 ) {
                              num_opens = 1;
                              num_extensions = profile_profile_alignment[ 0 ] - 1;
                            }
                            for( uint32_t i = 1; i < profile_profile_alignment.size(); i++ ) {
                              if( profile_profile_alignment[ i ] == 0 ) {
                                if( profile_profile_alignment[ i - 1 ] == 0 ) {
                                  num_extensions += 1;
                                } else {
                                  num_opens += 1;
                                }
                              } else {
                                numProfileProfileAlignedPositions_true[ test_id ] += 1;
                                if( profile_profile_alignment[ i ] > 1 ) {
                                  num_opens += 1;
                                  num_extensions += ( profile_profile_alignment[ i ] - 2 );
                                }
                              }
                            }  // End foreach alignment position
                            distanceSKLPositions_aligned_true[ test_id ] =
                              ( profile_profile_alignment_cost - ( num_opens * m_parameters.profileProfileIndelOpenCost ) - ( num_extensions * m_parameters.profileProfileIndelExtensionCost ) );
                            if( be_verbose ) {
                              cout << "After aligning the two profiles, got alignment cost: " << profile_profile_alignment_cost << endl;
                              cout << "The alignment is ( " << profile_profile_alignment[ 0 ];
                              for( uint32_t i = 1; i < profile_profile_alignment.size(); i++ ) {
                                cout << ", " << profile_profile_alignment[ i ];
                              }  // End foreach alignment position
                              cout << " )" << endl;
                              cout << "There were " << numProfileProfileAlignedPositions_true[ test_id ] << " aligned positions." << endl;
                              cout << "The symmeterized KL divergence of the aligned positions is " << distanceSKLPositions_aligned_true[ test_id ] << endl;
                              cout << "The average symmeterized KL divergence over the aligned positions is " << ( distanceSKLPositions_aligned_true[ test_id ] / numProfileProfileAlignedPositions_true[ test_id ] ) << endl;
                            } // End if be_verbose
                          } // End if calculateProfileProfileAlignments
                        } // End foreach test_id, train the root and calculate the forward
                          // score and SKL distances and alignments
              
                        // TODO: REMOVE
                        //cout << "=======YUP=======" << endl;

                        // Training sequences score using viterbi
                        if( m_parameters.testViterbi ) {
                          for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                            if( !tests[ test_id ].isRun ) {
                              continue;
                            }
                          
                            // TODO: REMOVE
                            //cout << "Viterbi-ing test id " << test_id << endl;
                            //cout << "sequences are { " << endl;
                            //for( uint32_t node_i = 0; node_i < testProfileTree[ test_id ].nodeCount(); node_i++ ) {
                            //  cout << "\t" << node_i << ":" << endl;
                            //  for( uint32_t seq_i = 0; seq_i < testProfileSequences[ test_id ][ node_i ].size(); seq_i++ ) {
                            //    cout << testProfileSequences[ test_id ][ node_i ][ seq_i ] << endl;
                            //  }
                            //}
                            //cout << " }" << endl;

                            // TODO: Put back
                            //  dp.forward_score_viterbi(
                            //    m_parameters,
                            //    testProfileTree[ test_id ],
                            //    testProfileSequences[ test_id ],
                            //    training_forward_matrices
                            //  );
                            dp.forward_score(
                              m_parameters,
                              true, // do use viterbi
                              *testProfileTree[ test_id ].getProfileTreeRoot(),
                              training_fasta,
                              training_fasta.size(),
                              training_forward_rows_1,
                              training_forward_rows_2,
                              testScore_training_viterbi[ test_id ]
                            );
                            if( be_verbose ) {    
                              cout << "The " << tests[ test_id ].name << " total score for all training sequences, using viterbi, is: " << testScore_training_viterbi[ test_id ] << endl;
                              //cout << "The corresponding viterbi forward matrices for the training sequences are: " << endl;
                              //cout << training_forward_matrices << endl;
                              //ProfileTreeTypeMultipleAlignment ma =
                              //  dp.forward_viterbiAlign(
                              //    m_parameters,
                              //    testProfileTree[ test_id ],
                              //    training_fasta,
                              //    trainer.m_sequence_count,
                              //    training_forward_matrices
                              //  );
                              //cout << "Conditional profile Multiple Alignment for training sequences is:" << endl;
                              ////cout << training_conditional_ma << endl;
                              //training_conditional_ma.toPairwiseStream( cout );
                            } // End if be_verbose
                          } // End foreach test_id, calculate the viterbi score for the
                          // training sequences.
                        } // End if testViterbi
                        
                        // Training sequences score using truepath
                        if( m_parameters.testTruepath ) {
                          for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                            if( !tests[ test_id ].isRun ) {
                              continue;
                            }
                          
                            // TODO: REMOVE
                            //cout << "Truepath-ing test id " << test_id << endl;
                            //cout << "sequences are { " << endl;
                            //for( uint32_t node_i = 0; node_i < testProfileTree[ test_id ].nodeCount(); node_i++ ) {
                            //  cout << "\t" << node_i << ":" << endl;
                            //  for( uint32_t seq_i = 0; seq_i < testProfileSequences[ test_id ][ node_i ].size(); seq_i++ ) {
                            //    cout << testProfileSequences[ test_id ][ node_i ][ seq_i ] << endl;
                            //  }
                            //}
                            //cout << " }" << endl;

                            if(
                              testProfileTree[ test_id ].getProfileTreeRoot()->length() !=
                              testProfileTree[ TEST_ID_true ].getProfileTreeRoot()->length()
                            ) {
                              testScore_training_truepath[ test_id ] = 0;
                            } else {
                              testScore_training_truepath[ test_id ] =
                                training_root_ma.calculateScore(
                                  *( testProfileTree[ test_id ].getProfileTreeRoot() )
                                );
                            } // End if the lengths are different .. else ..
                            if( be_verbose ) {    
                              cout << "The " << tests[ test_id ].name << " total score for all training sequences, using truepath, is: " << testScore_training_truepath[ test_id ] << endl;
                            } // End if be_verbose
                          } // End foreach test_id, calculate the truepath score for the
                          // training sequences.
                        } // End if testTruepath

                        for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                          if( !tests[ test_id ].isRun ) {
                            continue;
                          }
              
                          // TODO: REMOVE
                          //cout << "Calculating testing-set forward score for test id " << test_id << endl;
                          //cout << "sequences are { " << endl;
                          //for( uint32_t node_i = 0; node_i < testing_fastas.size(); node_i++ ) {
                          //  cout << "\t" << node_i << ":" << endl;
                          //  for( uint32_t seq_i = 0; seq_i < testProfileSequences[ test_id ][ node_i ].size(); seq_i++ ) {
                          //    cout << testing_fastas[ node_i ][ seq_i ] << endl;
                          //  }
                          //}
                          //cout << " }" << endl;
                          //cout << "profile tree is" << endl;
                          //cout << testProfileTree[ test_id ] << endl;
                          //cout << "profile length is " << testProfileTree[ test_id ].getProfileTreeRoot()->length() << endl;
                          //cout << "(test) forward matrices size is " << testing_forward_matrices.size() << endl;
                          // TODO: Put back
                          //  dp.forward_score(
                          //    m_parameters,
                          //    testProfileTree[ test_id ],
                          //    testing_fastas, // TODO: retrain just seq ids.
                          //    testing_forward_matrices
                          //  );
                          dp.forward_score(
                            m_parameters,
                            false, // do not use viterbi
                            *testProfileTree[ test_id ].getProfileTreeRoot(),
                            testing_fasta,
                            testing_fasta.size(),
                            testing_forward_rows_1,
                            testing_forward_rows_2,
                            testScore_test_forward[ test_id ]
                          );
                          if( be_verbose ) {    
                            cout << "The " << tests[ test_id ].name << " total score for all test sequences, using forward, is: " << testScore_test_forward[ test_id ] << endl;
                            //cout << "The corresponding forward matrices for the test sequences are: " << endl;
                            //cout << testing_forward_matrices << endl;
                          } // End if be_verbose
                          if( tests[ test_id ].isGibbs ) {
                            // First mean, then mode
                            if( m_parameters.reportGibbsMean ) {
                              // TODO: Put back
                              //testScore_test_mean_forward[ test_id ] =
                              //  dp.forward_score(
                              //    m_parameters,
                              //    testProfileTree_mean[ test_id ],
                              //    testing_fastas,
                              //    testing_forward_matrices
                              //  );
                              dp.forward_score(
                                m_parameters,
                                false, // do not use viterbi
                                *testProfileTree_mean[ test_id ].getProfileTreeRoot(),
                                testing_fasta,
                                testing_fasta.size(),
                                testing_forward_rows_1,
                                testing_forward_rows_2,
                                testScore_test_mean_forward[ test_id ]
                              );
                              if( be_verbose ) {    
                                cout << "The " << tests[ test_id ].name << " mean profile total score for all test sequences, using forward, is: " << testScore_test_mean_forward[ test_id ] << endl;
                                //cout << "The corresponding forward matrices for the test sequences are: " << endl;
                                //cout << testing_forward_matrices << endl;
                              } // End if be_verbose
                            } // End if reportGibbsMean

                            // mode
                            if( m_parameters.reportGibbsMode ) {
                              // TODO: Put back
                              //testScore_test_mode_forward[ test_id ] =
                              //  dp.forward_score(
                              //    m_parameters,
                              //    testProfileTree_mode[ test_id ],
                              //    testing_fastas,
                              //    testing_forward_matrices
                              //  );
                              dp.forward_score(
                                m_parameters,
                                false, // do not use viterbi
                                *testProfileTree_mode[ test_id ].getProfileTreeRoot(),
                                testing_fasta,
                                testing_fasta.size(),
                                testing_forward_rows_1,
                                testing_forward_rows_2,
                                testScore_test_mode_forward[ test_id ]
                              );
                              if( be_verbose ) {    
                                cout << "The " << tests[ test_id ].name << " mode profile total score for all test sequences, using forward, is: " << testScore_test_mode_forward[ test_id ] << endl;
                                //cout << "The corresponding forward matrices for the test sequences are: " << endl;
                                //cout << testing_forward_matrices << endl;
                              } // End if be_verbose
                            } // End if reportGibbsMode
                          } // End if isGibbs
                        } // End foreach test_id, calculate the forward score for the test
                          // sequences.
                        
                        // Test sequences score using viterbi
                        if( m_parameters.testViterbi ) {
                          for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                            if( !tests[ test_id ].isRun ) {
                              continue;
                            }
                          
                            // TODO: REMOVE
                            //cout << "Calculating testing-set viterbi score for test id " << test_id << endl;
                            // TODO: Put back
                            //testScore_test_viterbi[ test_id ] =
                            //  dp.forward_score_viterbi(
                            //    m_parameters,
                            //    testProfileTree[ test_id ],
                            //    testing_fastas, // TODO: retrain just seq ids.
                            //    testing_forward_matrices
                            //  );
                            dp.forward_score(
                              m_parameters,
                              true, // do use viterbi
                              *testProfileTree[ test_id ].getProfileTreeRoot(),
                              testing_fasta,
                              testing_fasta.size(),
                              testing_forward_rows_1,
                              testing_forward_rows_2,
                              testScore_test_viterbi[ test_id ]
                            );

                            if( be_verbose ) {    
                              cout << "The " << tests[ test_id ].name << " total score for all test sequences, using viterbi, is: " << testScore_test_viterbi[ test_id ] << endl;
                              //cout << "The corresponding viterbi forward matrices for the test sequences are: " << endl;
                              //cout << testing_forward_matrices << endl;
                              //ProfileTreeTypeMultipleAlignment ma =
                              //  dp.forward_viterbiAlign(
                              //    m_parameters,
                              //    testProfileTree[ test_id ],
                              //    testing_fasta,
                              //    trainer.m_sequence_count,
                              //    testing_forward_matrices
                              //  );
                              //cout << "Conditional profile Multiple Alignment for test sequences is:" << endl;
                              ////cout << test_conditional_ma << endl;
                              //test_conditional_ma.toPairwiseStream( cout );
                            } // End if be_verbose
                          } // End foreach test_id, calculate the viterbi score for the
                          // test sequences.
                        } // End if testTruepath

                        // Test sequences score using truepath
                        if( m_parameters.testTruepath ) {
                          for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                            if( !tests[ test_id ].isRun ) {
                              continue;
                            }
                          
                            // TODO: REMOVE
                            //cout << "Calculating testing-set truepath score for test id " << test_id << endl;
                          
                            if(
                              testProfileTree[ test_id ].getProfileTreeRoot()->length() !=
                              testProfileTree[ TEST_ID_true ].getProfileTreeRoot()->length()
                            ) {
                              testScore_test_truepath[ test_id ] = 0;
                            } else {
                              testScore_test_truepath[ test_id ] =
                                testing_root_ma.calculateScore(
                                  *( testProfileTree[ test_id ].getProfileTreeRoot() )
                                );
                            } // End if lengths are different .. else ..
                            if( be_verbose ) {    
                              cout << "The " << tests[ test_id ].name << " total score for all test sequences, using truepath, is: " << testScore_test_truepath[ test_id ] << endl;
                            } // End if be_verbose
                          } // End foreach test_id, calculate the truepath score for the
                          // test sequences.
                        } // End if testTruepath
                        
                        // Summary:
                        if( be_verbose ) {
                          cout << endl << "Summary for starting profile #" << starting_profile_i << endl;
                          if( m_parameters.testViterbi ) {
                            for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                              if( tests[ test_id ].isRun ) {
                                cout << "The " << tests[ test_id ].name << " total score for all training sequences, using viterbi, is: " << testScore_training_viterbi[ test_id ] << endl;
                              }
                            } // End foreach test_id, cout the training viterbi score.
                            cout << endl;
                            for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                              if( tests[ test_id ].isRun ) {
                                cout << "The " << tests[ test_id ].name << " total score for all test sequences, using viterbi, is: " << testScore_test_viterbi[ test_id ] << endl;
                              }
                            } // End foreach test_id, cout the test viterbi score.
                            cout << endl;
                          } // End if testViterbi

                          if( m_parameters.testTruepath ) {
                            for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                              if( tests[ test_id ].isRun ) {
                                cout << "The " << tests[ test_id ].name << " total score for all training sequences, using truepath, is: " << testScore_training_truepath[ test_id ] << endl;
                              }
                            } // End foreach test_id, cout the training truepath score.
                            cout << endl;
                            for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                              if( tests[ test_id ].isRun ) {
                                cout << "The " << tests[ test_id ].name << " total score for all test sequences, using truepath, is: " << testScore_test_truepath[ test_id ] << endl;
                              }
                            } // End foreach test_id, cout the test truepath score.
                            cout << endl;
                          } // End if testTruepath

                          for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                            if( tests[ test_id ].isRun ) {
                              cout << "The " << tests[ test_id ].name << " total score for all training sequences, using forward, is: " << testScore_training_forward[ test_id ] << endl;
                            }
                          } // End foreach test_id, cout the training forward score.
                          cout << endl;
                          for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                            if( tests[ test_id ].isRun ) {
                              cout << "The " << tests[ test_id ].name << " total score for all test sequences, using forward, is: " << testScore_test_forward[ test_id ] << endl;
                            }
                          } // End foreach test_id, cout the test forward score.
                          for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                            if(
                              ( test_id == TEST_ID_true ) ||
                              ( test_id == TEST_ID_starting )
                            ) {
                              continue;
                            }
                            if( tests[ test_id ].isRun ) {
                              cout << "The " << tests[ test_id ].name << " took " << testIters[ test_id ] << " iterations (" << testCPUtime[ test_id ] << " CPU secs)." << endl;
                            }
                          } // End foreach test_id, cout the test iters.
                          cout << endl;
                        } // End if be_verbose

                        if( m_parameters.saveResultsToFile ) {
                          // print the tab results for this starting profile
                          tab_stream << conservation_rate << "\t";
                          tab_stream << profile_length << "\t";
                          tab_stream << num_training_sequences_per_profile << "\t";
                          tab_stream << expected_deletions_count << "\t";
                          tab_stream << expected_insertions_count << "\t";
                          tab_stream << expected_deletion_length_as_profile_length_fraction << "\t";
                          tab_stream << expected_insertion_length_as_profile_length_fraction << "\t";
                          if( m_parameters.numTrueProfiles > 1 ) {
                            tab_stream << true_profile_i << "\t";
                          }
                          cout << conservation_rate << " ";
                          cout << profile_length << " ";
                          cout << num_training_sequences_per_profile << " ";
                          cout << expected_deletions_count << " ";
                          cout << expected_insertions_count << " ";
                          cout << expected_deletion_length_as_profile_length_fraction << " ";
                          cout << expected_insertion_length_as_profile_length_fraction << " ";
                          if( m_parameters.numTrueProfiles > 1 ) {
                            cout << true_profile_i << " ";
                          }

                          for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                            if( tests[ test_id ].isRun ) {
                              if( m_parameters.calculateSymmeterizedKullbackLeiblerDistancesToTrue ) {
                                tab_stream << distanceSKLPositions_true[ test_id ] << "\t";
                                if( ( test_id != TEST_ID_true ) && m_parameters.coutDistances && tests[ test_id ].isCout ) {
                                  cout << distanceSKLPositions_true[ test_id ] << " ";
                                }
                                tab_stream << distanceSKLExceptPositions_true[ test_id ] << "\t";
                                if( ( test_id != TEST_ID_true ) && m_parameters.coutDistances && tests[ test_id ].isCout ) {
                                  cout << distanceSKLExceptPositions_true[ test_id ] << " ";
                                }
                              } // End if calculateSymmeterizedKullbackLeiblerDistancesToTrue
                              if( m_parameters.calculateSymmeterizedKullbackLeiblerDistancesToStarting ) {
                                tab_stream << distanceSKLPositions_starting[ test_id ] << "\t";
                                if( ( !m_parameters.calculateSymmeterizedKullbackLeiblerDistancesToTrue || ( test_id != TEST_ID_true ) ) && ( test_id != TEST_ID_starting ) && m_parameters.coutDistances && tests[ test_id ].isCout ) {
                                  cout << distanceSKLPositions_starting[ test_id ] << " ";
                                }
                                tab_stream << distanceSKLExceptPositions_starting[ test_id ] << "\t";
                                if( ( !m_parameters.calculateSymmeterizedKullbackLeiblerDistancesToTrue || ( test_id != TEST_ID_true ) ) && ( test_id != TEST_ID_starting ) && m_parameters.coutDistances && tests[ test_id ].isCout ) {
                                  cout << distanceSKLExceptPositions_starting[ test_id ] << " ";
                                }
                              } // End if calculateSymmeterizedKullbackLeiblerDistancesToStarting
                            } // End if isRun
                          } // End foreach test_id, print distances
                  
                          for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                            if(
                              ( test_id == TEST_ID_true ) ||
                              ( test_id == TEST_ID_starting )
                            ) {
                              continue;
                            }
                            if( tests[ test_id ].isRun ) {
                              if( m_parameters.calculateProfileProfileAlignments ) {
                                tab_stream << distanceSKLPositions_aligned_true[ test_id ] << "\t";
                                if( tests[ test_id ].isCout ) {
                                  cout << distanceSKLPositions_aligned_true[ test_id ] << " ";
                                }
                                tab_stream << numProfileProfileAlignedPositions_true[ test_id ] << "\t";
                                if( tests[ test_id ].isCout ) {
                                  cout << numProfileProfileAlignedPositions_true[ test_id ] << " ";
                                }
                              } // End if calculateProfileProfileAlignments
                              if( true ) {
                                tab_stream << testProfileTree[ test_id ].getProfileTreeRoot()->length() << "\t";
                                if( tests[ test_id ].isCout ) {
                                  cout << testProfileTree[ test_id ].getProfileTreeRoot()->length() << " ";
                                }
                              } // End if true
                            } // End if isRun
                          } // End foreach test_id, print profile-profile alignment stats and profile lengths

                          for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                            if( tests[ test_id ].isRun ) {
                              if( !( ( test_id == TEST_ID_true ) || ( test_id == TEST_ID_starting ) ) ) {
                                tab_stream << testIters[ test_id ] << "\t";
                                if( tests[ test_id ].isCout ) {
                                  cout << testIters[ test_id ];
                                  cout << " ";
                                }
                                tab_stream << testCPUtime[ test_id ] << "\t";
                                if( tests[ test_id ].isCout ) {
                                  cout << testCPUtime[ test_id ];
                                  cout << " ";
                                }
                              } // End if !( TEST_ID_true || TEST_ID_starting )
                              if( convert_tab_output_to_log_double ) {
                                tab_stream << toLogDouble( testScore_training_forward[ test_id ] ) << "\t";
                              } else {
                                tab_stream << testScore_training_forward[ test_id ] << "\t";
                              }
                              if( tests[ test_id ].isCout ) {
                                cout << tests[ test_id ].coutLeftBrace;
                                cout << testScore_training_forward[ test_id ];
                                cout << tests[ test_id ].coutRightBrace;
                                cout << " ";
                              }
                              if( tests[ test_id ].isGibbs ) {
                                // Mean, then mode
                                if( m_parameters.reportGibbsMean ) {
                                  if( convert_tab_output_to_log_double ) {
                                    tab_stream << toLogDouble( testScore_training_mean_forward[ test_id ] ) << "\t";
                                  } else {
                                    tab_stream << testScore_training_mean_forward[ test_id ] << "\t";
                                  }
                                  if( tests[ test_id ].isCout ) {
                                    cout << tests[ test_id ].coutLeftBrace;
                                    cout << testScore_training_mean_forward[ test_id ];
                                    cout << tests[ test_id ].coutRightBrace;
                                    cout << " ";
                                  }
                                } // End if reportGibbsMean
                                // mode
                                if( m_parameters.reportGibbsMode ) {
                                  if( convert_tab_output_to_log_double ) {
                                    tab_stream << toLogDouble( testScore_training_mode_forward[ test_id ] ) << "\t";
                                  } else {
                                    tab_stream << testScore_training_mode_forward[ test_id ] << "\t";
                                  }
                                  if( tests[ test_id ].isCout ) {
                                    cout << tests[ test_id ].coutLeftBrace;
                                    cout << testScore_training_mode_forward[ test_id ];
                                    cout << tests[ test_id ].coutRightBrace;
                                    cout << " ";
                                  }
                                } // End if reportGibbsMode
                              } // End if isGibbs
                            } // End if isRun
                          } // End foreach test_id, print testScore_training_forward[ test_id ]
                  
                          if( m_parameters.testViterbi ) {
                            for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                              if( tests[ test_id ].isRun ) {
                                if( convert_tab_output_to_log_double ) {
                                  tab_stream << toLogDouble( testScore_training_viterbi[ test_id ] ) << "\t";
                                } else {
                                  tab_stream << testScore_training_viterbi[ test_id ] << "\t";
                                }
                                if( m_parameters.coutViterbi && tests[ test_id ].isCout ) {
                                  //if( test_id == TEST_ID_conditional ) {
                                  //  cout << "<";
                                  //} else if( test_id == TEST_ID_unconditional ) {
                                  //  cout << "(";
                                  //}
                                  cout << tests[ test_id ].coutLeftBrace;
                                  cout << testScore_training_viterbi[ test_id ];
                                  cout << tests[ test_id ].coutRightBrace;
                                  //if( test_id == TEST_ID_conditional ) {
                                  //  cout << ">";
                                  //} else if( test_id == TEST_ID_unconditional ) {
                                  //  cout << ")";
                                  //}
                                  cout << " ";
                                }
                              }
                            } // End foreach test_id, print testScore_training_viterbi[ test_id ]
                          } // End if testViterbi

                          if( m_parameters.testTruepath ) {
                            for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                              if( tests[ test_id ].isRun ) {
                                if( convert_tab_output_to_log_double ) {
                                  tab_stream << toLogDouble( testScore_training_truepath[ test_id ] ) << "\t";
                                } else {
                                  tab_stream << testScore_training_truepath[ test_id ] << "\t";
                                }
                                if( m_parameters.coutTruepath && tests[ test_id ].isCout ) {
                                  //if( test_id == TEST_ID_conditional ) {
                                  //  cout << "<";
                                  //} else if( test_id == TEST_ID_unconditional ) {
                                  //  cout << "(";
                                  //}
                                  cout << tests[ test_id ].coutLeftBrace;
                                  cout << testScore_training_truepath[ test_id ];
                                  cout << tests[ test_id ].coutRightBrace;
                                  //if( test_id == TEST_ID_conditional ) {
                                  //  cout << ">";
                                  //} else if( test_id == TEST_ID_unconditional ) {
                                  //  cout << ")";
                                  //}
                                  cout << " ";
                                }
                              }
                            } // End foreach test_id, print testScore_training_truepath[ test_id ]
                          } // End if testTruepath

                          // After this mark are the test set scores (as opposed to training set)
                          cout << "<==--#\t" << endl << "#--==>\t";
                          for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                            if( tests[ test_id ].isRun ) {
                              if( convert_tab_output_to_log_double ) {
                                tab_stream << toLogDouble( testScore_test_forward[ test_id ] ) << "\t";
                              } else {
                                tab_stream << testScore_test_forward[ test_id ] << "\t";
                              }
                              if( tests[ test_id ].isCout ) {
                                cout << tests[ test_id ].coutLeftBrace;
                                cout << tests[ test_id ].coutLeftBrace;
                                cout << testScore_test_forward[ test_id ];
                                cout << tests[ test_id ].coutRightBrace;
                                cout << tests[ test_id ].coutRightBrace;
                                cout << " ";
                              }
                              if( tests[ test_id ].isGibbs ) {
                                // Mean, then mode
                                if( m_parameters.reportGibbsMean ) {
                                  if( convert_tab_output_to_log_double ) {
                                    tab_stream << toLogDouble( testScore_test_mean_forward[ test_id ] ) << "\t";
                                  } else {
                                    tab_stream << testScore_test_mean_forward[ test_id ] << "\t";
                                  }
                                  if( tests[ test_id ].isCout ) {
                                    cout << tests[ test_id ].coutLeftBrace;
                                    cout << tests[ test_id ].coutLeftBrace;
                                    cout << testScore_test_mean_forward[ test_id ];
                                    cout << tests[ test_id ].coutRightBrace;
                                    cout << tests[ test_id ].coutRightBrace;
                                    cout << " ";
                                  }
                                } // End if reportGibbsMean
                                // mode
                                if( m_parameters.reportGibbsMode ) {
                                  if( convert_tab_output_to_log_double ) {
                                    tab_stream << toLogDouble( testScore_test_mode_forward[ test_id ] ) << "\t";
                                  } else {
                                    tab_stream << testScore_test_mode_forward[ test_id ] << "\t";
                                  }
                                  if( tests[ test_id ].isCout ) {
                                    cout << tests[ test_id ].coutLeftBrace;
                                    cout << tests[ test_id ].coutLeftBrace;
                                    cout << testScore_test_mode_forward[ test_id ];
                                    cout << tests[ test_id ].coutRightBrace;
                                    cout << tests[ test_id ].coutRightBrace;
                                    cout << " ";
                                  }
                                } // End if reportGibbsMode
                              } // End if isGibbs
                            } // End if isRun
                          } // End foreach test_id, print testScore_test_forward[ test_id ]
              
                          if( m_parameters.testViterbi ) {
                            for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                              if( tests[ test_id ].isRun ) {
                                if( convert_tab_output_to_log_double ) {
                                  tab_stream << toLogDouble( testScore_test_viterbi[ test_id ] ) << "\t";
                                } else {
                                  tab_stream << testScore_test_viterbi[ test_id ] << "\t";
                                }
                                if( m_parameters.coutViterbi && tests[ test_id ].isCout ) {
                                  //if( test_id == TEST_ID_conditional ) {
                                  //  cout << "<<";
                                  //} else if( test_id == TEST_ID_unconditional ) {
                                  //  cout << "((";
                                  //}
                                  cout << tests[ test_id ].coutLeftBrace;
                                  cout << tests[ test_id ].coutLeftBrace;
                                  cout << testScore_test_viterbi[ test_id ];
                                  cout << tests[ test_id ].coutRightBrace;
                                  cout << tests[ test_id ].coutRightBrace;
                                  //if( test_id == TEST_ID_conditional ) {
                                  //  cout << ">>";
                                  //} else if( test_id == TEST_ID_unconditional ) {
                                  //  cout << "))";
                                  //}
                                  cout << " ";
                                }
                              }
                            } // End foreach test_id, print testScore_test_viterbi[ test_id ]
                          } // End if testViterbi

                          if( m_parameters.testTruepath ) {
                            for( test_id = 0; test_id <= last_test_id; test_id++ ) {
                              if( tests[ test_id ].isRun ) {
                                if( convert_tab_output_to_log_double ) {
                                  tab_stream << toLogDouble( testScore_test_truepath[ test_id ] ) << "\t";
                                } else {
                                  tab_stream << testScore_test_truepath[ test_id ] << "\t";
                                }
                                if( m_parameters.coutTruepath && tests[ test_id ].isCout ) {
                                  //if( test_id == TEST_ID_conditional ) {
                                  //  cout << "<<";
                                  //} else if( test_id == TEST_ID_unconditional ) {
                                  //  cout << "((";
                                  //}
                                  cout << tests[ test_id ].coutLeftBrace;
                                  cout << tests[ test_id ].coutLeftBrace;
                                  cout << testScore_test_truepath[ test_id ];
                                  cout << tests[ test_id ].coutRightBrace;
                                  cout << tests[ test_id ].coutRightBrace;
                                  //if( test_id == TEST_ID_conditional ) {
                                  //  cout << ">>";
                                  //} else if( test_id == TEST_ID_unconditional ) {
                                  //  cout << "))";
                                  //}
                                  cout << " ";
                                }
                              }
                            } // End foreach test_id, print testScore_test_truepath[ test_id ]
                          } // End if testTruepath
                
                          tab_stream << endl;
                          cout << endl;
                        } // End if m_parameters.saveResultsToFile
                    
                      } // End foreach starting_profile_i
                    
                    } // End foreach true_profile_i
                  } // End foreach expected_insertion_length_as_profile_length_fraction_i
                } // End foreach expected_deletion_length_as_profile_length_fraction_i
              } // End foreach expected_insertions_count_i
            } // End foreach expected_deletions_count_i
          } // End foreach num_training_sequences_per_profile_i
        } // End foreach profile_length_i
      } // End foreach conservation_rate_i

      if( m_parameters.saveResultsToFile ) {
        // End and close the file.
        //tab_stream << "</profuse_test_results>" << endl;
        tab_stream.close();
      } // End if m_parameters.saveResultsToFile

    } // start()
} // End namespace galosh

#endif // __GALOSH_PROFUSETEST_HPP__
