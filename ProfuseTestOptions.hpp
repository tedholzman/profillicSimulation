/**
 * \file ProfuseTestOptions.hpp
 * \author Ted Holzman (tholzman@scharp.org)
 * \par Library:
 * Galosh ProfuseTest
 * \brief Options for ProfuseTest
 * \copyright &copy; 2008, 2011, 2012 by Paul T. Edlefsen, Fred Hutchinson Cancer
 *    Research Center.
 *  All rights reserved.
 *****************************************************************************/



///Useful for modifying GALOSH_DEF_OPT behavior.  See vector-valued options below
#ifndef TMP_EXTRA_STUFF
#define TMP_EXTRA_STUFF BOOST_PP_EMPTY()
#endif

/**
 * \note Note the lack of "protective" \#define symbols in this file.  We may
 * want to include it several times, redefining GALOSH_DEF_OPT as needed.
 */

/**
 * The configuration file contains the non-default, non-command-line options
 *
 */
GALOSH_DEF_OPT(configFile,string,"ProfuseTest.cfg","File path for the configuration file");

/** Save file version
 * 1 was the beginning
 * 2 is after I modified the conditional_then_unconditional_root stuff to use the globals from the conditional, but the position-specific values from the starting profile.
 * 3 is after I added unconditional_with_fixed_starting_globals
 * 4 is after I added unconditional_with_fixed_starting_globals_then_with_fixed_positions
 * 5 is after I added startWithUniformGlobals
 * 6 is after I added startWithUniformGlobals_scalar
 * 7 is after I added cout[Test] params
 * 8 is after I added startWithUniformPositions, startWithPositionsDrawnFromPrior, etc.
 * 9 is after I added convert_tab_output_to_log_double
 *10 is after I added CPU time
 *11 TAH 7/12 moved options to BOOST program options
 */
GALOSH_DEF_OPT(saveFileVersion,uint32_t,11U,"Version of the save file");

GALOSH_DEF_OPT(saveResultsToFile,bool,true,"Should we save the results to a file?");
GALOSH_DEF_OPT(saveResultsParentDirectory,string,".","Parent directory name");



///File prefixes and suffixes, and output options
GALOSH_DEF_OPT(resultsFilePrefix,string,"ProfuseTest.","Prefix for main results");
GALOSH_DEF_OPT(tabFileSuffix,string,".tab","Suffix for tab-delimited output");
GALOSH_DEF_OPT(parametersFileSuffix,string,".Parameters.xml","Suffix for parameters file");
GALOSH_DEF_OPT(saveTrueProfileTrees,bool,true,"Shall we save 'true' Profile Trees?");
GALOSH_DEF_OPT(trueProfileTreeFileSuffix,string,".true.ProfileTree.xml","Suffix for profile tree file");
GALOSH_DEF_OPT(saveStartingProfiles,bool,true,"Shall we save starting profiles?");
GALOSH_DEF_OPT(startingProfileTreeFileSuffix,string,".starting.ProfileTree.xml","Starting profile file suffix.");
GALOSH_DEF_OPT(saveTestProfiles,bool,true,"Shall we save test profiles?");
GALOSH_DEF_OPT(testProfileTreeFileSuffix,string,".ProfileTree.xml","Suffix of profile tree file.");
GALOSH_DEF_OPT(savePatternSequences,bool,true,"Shall we save pattern files?");
GALOSH_DEF_OPT(patternSequencesFileSuffix,string,".pattern_sequences.fasta","Pattern fasta file suffix");
GALOSH_DEF_OPT(saveTests,bool,true,"Shall we save test information?");
GALOSH_DEF_OPT(testsFileSuffix,string,".tests","Tests file suffix.");
GALOSH_DEF_OPT(saveTrainingSequences,bool,true,"Shall we save fasta training sequences?");
GALOSH_DEF_OPT(trainingSequencesFileSuffix,string,".training_sequences.fasta","Suffix for training sequence fasta files");
GALOSH_DEF_OPT(saveTestingSequences,bool,true,"Shall we save testing sequences?");
GALOSH_DEF_OPT(testingSequencesFileSuffix,string,".testing_sequences.fasta","Test sequence fasta file suffix");
GALOSH_DEF_OPT(saveTrueTrainingAlignments,bool,true,"Shall we save true training alignments?");
GALOSH_DEF_OPT(trainingTrueAlignmentsFileSuffix,string,".training_alignments.fasta","True training alignments fasta file suffix");
GALOSH_DEF_OPT(saveTrueTestingAlignments,bool,true,"Shall we save true test alignments?");
GALOSH_DEF_OPT(trueTestingAlignmentsFileSuffix,string,".testing_alignments.fasta","True test alignment fasta file suffix");

/**
 * Number of times we start from the very beginning (with a new pattern sequence).
 */
GALOSH_DEF_OPT(numProfiles,uint32_t,1,"Number of root profiles to run");
/**
 * For each "root" profile, do the whole simulation process a number of different times, starting over with
 * a new pattern sequence.
 */
GALOSH_DEF_OPT(numTrueProfiles,uint32_t,4,"We do the whole thing a number of different times, starting over with a new pattern sequence.");
/**
  * For each true root profile, we run the trainers from a number of
  * different starting profiles.  This is that number.
  */
GALOSH_DEF_OPT(numStartingProfiles,uint32_t,4,"Number of starting profiles per root profile.");
GALOSH_DEF_OPT(numTestingSequencesPerProfile,uint32_t,20,"Number of test sequences per profile");

/**
 * When numProfiles > 1, what fraction of the positions of each child
 * should be shared with its parent?
 */
GALOSH_DEF_OPT(sharedPositionRate,double,1.0,"Fraction of positions of each child to be shared with the parent");

/**
 * Lock the indel parameters of the true profile to be the same for
 * insertions as for deletions?  This makes the expectedInsertionsCounts,
 * expectedInsertionLengthAsProfileLengthFractions, and
 * minExpectedInsertionLength unused, since the corresponding deletion
 * values will be used instead.  It also reduces the number of tests by
 * reducing the number of possible combinations (since deletions and
 * insertions will go in lock step).
 */
GALOSH_DEF_OPT(useDeletionsForInsertionsParameters,bool,true,"Should insertion prob be forced equal to deletion prob?");

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
GALOSH_DEF_OPT(minExpectedDeletionLength,double,1.25,"The deletionExtension value of the true profile will be the minimum of  ( 1.0 / ( expectedDeletionLengthAsProfileLengthFraction * profileLength ) ) and ( 1.0 / minExpectedDeletionLength )");

 /**
  * If useDeletionsForInsertionsParameters is false, the
  * insertionExtension value of the true profile will be the minimum of (
  * 1.0 / ( expectedInsertionLengthAsProfileLengthFraction * profileLength
  * ) ) and ( 1.0 / minExpectedInsertionLength ).
  *
  * @see useDeletionsForInsertionsParameters
  */
GALOSH_DEF_OPT(minExpectedInsertionLength,double,1.25,"If useDeletionsForInsertionsParameters is false, the insertionExtension value of the true profile will be the minimum of (1.0 / ( expectedInsertionLengthAsProfileLengthFraction * profileLength) ) and ( 1.0 / minExpectedInsertionLength ).");

 /**
  * The preAlignInsertion value of the true profile.
  */
GALOSH_DEF_OPT(preAlignInsertion,double,.01,"The preAlignInsertion value of the true profile");

 /**
  * The postAlignInsertion value of the true profile.
  */
GALOSH_DEF_OPT(postAlignInsertion,double,.01,"The postAlignInsertion value of the true profile");

 /**
  * The effective number of sequences "observed" a priori.  Note that we
  * use a different prior strength for main-model transitions: see
  * priorStrength_internal_transitions.
  */
GALOSH_DEF_OPT(priorStrength,float,1.0f,"The effective number of sequences \"observed\" a priori");

 /**
  * The effective number of sequences "observed" a priori, for main-model
  * transitions.
  */
GALOSH_DEF_OPT(priorStrength_internal_transitions,float,10.0f,"The effective number of sequences \"observed\" a priori, for main-model transitions");

 /**
  * The prior contribution (per "a priori sequence": see priorStrength) of
  * M->M transitions.  This will be multiplied by the profile length and
  * by the priorStrength when setting up the global prior.
  */
GALOSH_DEF_OPT(priorMtoM,float,.95f,"The prior contribution (per \"a priori sequence\": see priorStrength) of M->M transitions");

 /**
  * The prior contribution (per "a priori sequence": see priorStrength) of
  * M->I transitions.  This will be multiplied by the profile length and
  * by the priorStrength when setting up the global prior.
  */
GALOSH_DEF_OPT(priorMtoI,float,.025f,"The prior contribution (per \"a priori sequence\": see priorStrength) of M->I transitions");

 /**
  * The prior contribution (per "a priori sequence": see priorStrength) of
  * M->D transitions.  This will be multiplied by the profile length and
  * by the priorStrength when setting up the global prior.
  */
GALOSH_DEF_OPT(priorMtoD,float,.025f,"The prior contribution (per \"a priori sequence\": see priorStrength) of M->D transitions");

 /**
  * The prior contribution (per "a priori sequence": see priorStrength) of
  * I->M transitions.  This will be multiplied by the profile length and
  * by the priorStrength when setting up the global prior.
  */
GALOSH_DEF_OPT(priorItoM,float,.05f,"The prior contribution (per \"a priori sequence\": see priorStrength) of I->M transitions");

 /**
  * The prior contribution (per "a priori sequence": see priorStrength) of
  * I->I transitions.  This will be multiplied by the profile length and
  * by the priorStrength when setting up the global prior.
  */
GALOSH_DEF_OPT(priorItoI,float,.95f,"The prior contribution (per \"a priori sequence\": see priorStrength) of I->I transitions");

 /**
  * The prior contribution (per "a priori sequence": see priorStrength) of
  * D->M transitions.  This will be multiplied by the profile length and
  * by the priorStrength when setting up the global prior.
  */
GALOSH_DEF_OPT(priorDtoM,float,.95f," The prior contribution (per \"a priori sequence\": see priorStrength) of D->M transitions");

 /**
  * The prior contribution (per "a priori sequence": see priorStrength) of
  * D->D transitions.  This will be multiplied by the profile length and
  * by the priorStrength when setting up the global prior.
  */
GALOSH_DEF_OPT(priorDtoD,float,.05f,"The prior contribution (per \"a priori sequence\": see priorStrength) of D->D transitions");

 /**
  * Additionally report the overall mean of all chains found while
  * performing Gibbs sampling?  The best profile is always reported, which
  * may be the overall mean, the mode, or the mean of one of the chains.
  */
GALOSH_DEF_OPT(reportGibbsMean,bool,false,"Additionally report the overall mean of all chains found while performing Gibbs sampling?");

 /**
  * Additionally report the mode found while
  * performing Gibbs sampling?  The best profile is always reported, which
  * may be the overall mean, the mode, or the mean of one of the chains.
  *
  * Note that it takes some extra time to store the mode (this turns on
  * saveGibbsMode in the ProfileGibbs class).
  */
GALOSH_DEF_OPT(reportGibbsMode,bool,false,"Additionally report the mode found while performing Gibbs sampling?");

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
GALOSH_DEF_OPT(startWithUniformGlobals,bool,false,"If startWithGlobalsDrawnFromPrior is not true, and if startWithUniformGlobals is true, then we set the global values of the startingProfile to random values between 0 and min(startWithUniformGlobals_scalar times the true values,startWithUniformGlobals_maxXtoY).  If it is false, we start with the known, true globals.");

 /**
  * @see startWithUniformGlobals
  */
GALOSH_DEF_OPT(startWithUniformGlobals_scalar,double,2.0,"See \"startWithUniformGlobals\"");

 /**
  * @see startWithUniformGlobals
  */
GALOSH_DEF_OPT(startWithUniformGlobals_maxNtoN,double,.2,"See \"startWithUniformGlobals\"");

 /**
  * @see startWithUniformGlobals
  */
GALOSH_DEF_OPT(startWithUniformGlobals_maxBtoD,double,.2,"See \"startWithUniformGlobals\"");

 /**
  * @see startWithUniformGlobals
  */
GALOSH_DEF_OPT(startWithUniformGlobals_maxMtoI,double,.2,"See \"startWithUniformGlobals\"");

 /**
  * @see startWithUniformGlobals
  */
GALOSH_DEF_OPT(startWithUniformGlobals_maxMtoD,double,.2,"See \"startWithUniformGlobals\"");

 /**
  * @see startWithUniformGlobals
  */
GALOSH_DEF_OPT(startWithUniformGlobals_maxItoI,double,.5,"See \"startWithUniformGlobals\"");

 /**
  * @see startWithUniformGlobals
  */
GALOSH_DEF_OPT(startWithUniformGlobals_maxDtoD,double,.5,"See \"startWithUniformGlobals\"");

 /**
  * @see startWithUniformGlobals
  */
GALOSH_DEF_OPT(startWithUniformGlobals_maxCtoC,double,.2,"See \"startWithUniformGlobals\"");

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
GALOSH_DEF_OPT(startWithUniformPositions,bool,false,"If startWithUniformPositions is true, then we set the position-specific values of the startingProfile to random values between 0 and 1.  If it is false, we start with the known, true parameter values.");

 /**
  * If startWithGlobalsDrawnFromPrior is true, the
  * global values of the starting profile will be drawn from the prior.
  *
  * @see startWithUniformGlobals
  */
GALOSH_DEF_OPT(startWithGlobalsDrawnFromPrior,bool,false,"If startWithGlobalsDrawnFromPrior is true, the global values of the starting profile will be drawn from the prior");

 /**
  * If startWithPositionsDrawnFromPrior is true, the
  * position-specific values of the starting profile will be drawn from
  * the prior... but see the notes in startWithUniformPositions.
  *
  * @see startWithUniformPositions
  * @see alsoStartWithEvenPositions
  */
GALOSH_DEF_OPT(startWithPositionsDrawnFromPrior,bool,false,"If startWithPositionsDrawnFromPrior is true, the position-specific values of the starting profile will be drawn from the prior");

 /**
  * Calculate the viterbi scores after training each profile?  Note that
  * this is in addition to the forward scores, which we always calculate.
  */
GALOSH_DEF_OPT(testViterbi,bool,true,"Calculate the viterbi scores after training each profile?");

 /**
  * Write the viterbi scores to STDOUT?
  */
GALOSH_DEF_OPT(coutViterbi,bool,false,"Write the viterbi scores to STDOUT?");

 /**
  * Calculate the truepath scores after training each profile?  Note that
  * this is in addition to the forward scores, which we always calculate.
  */
GALOSH_DEF_OPT(testTruepath,bool,true,"Calculate the truepath scores after training each profile?");

 /**
  * Write the truepath scores to STDOUT?
  */
GALOSH_DEF_OPT(coutTruepath,bool,false,"Write the truepath scores to STDOUT?");

 /**
  * Calculate the SKL distance between the training profiles and the
  * true profile?
  */
GALOSH_DEF_OPT(calculateSymmeterizedKullbackLeiblerDistancesToTrue,bool,true,"Calculate the SKL distance between the training profiles and the true profile?");

 /**
  * Calculate the SKL distance between the training profiles and the
  * starting profile?
  */
GALOSH_DEF_OPT(calculateSymmeterizedKullbackLeiblerDistancesToStarting,bool,false,"Calculate the SKL distance between the training profiles and the starting profile?");

 /**
  * Calculate the SKL distance between the training profiles and the
  * conditional profile?
  */
GALOSH_DEF_OPT(coutDistances,bool,true,"Calculate the SKL distance between the training profiles and the conditional profile?");

 /**
  * Calculate SKL profile-profile alignments between each trained profile
  * and the true profile?
  */
GALOSH_DEF_OPT(calculateProfileProfileAlignments,bool,true,"Calculate SKL profile-profile alignments between each trained profile and the true profile?");

 /**
  * The cost of a gap open when performing SKL profile-profile alignements.
  */
GALOSH_DEF_OPT(profileProfileIndelOpenCost,double,.25,"The cost of a gap open when performing SKL profile-profile alignements");

 /**
  * The cost of a gap extension when performing SKL profile-profile alignments.
  */
GALOSH_DEF_OPT(profileProfileIndelExtensionCost,double,.25,"The cost of a gap extension when performing SKL profile-profile alignment");

 /**
  * Calculate forward (etc) scores for the true profile, too?
  */
GALOSH_DEF_OPT(testTrueProfile,bool,true,"Calculate forward (etc) scores for the true profile");

 /**
  * Write forward (etc) scores for the true profile to STDOUT
  * during training?
  */
GALOSH_DEF_OPT(coutTrueProfile,bool,true,"Write forward (etc) scores for the true profile to STDOUT during training?");

 /**
  * Calculate forward (etc) scores for the pre-training profile, too?
  */
GALOSH_DEF_OPT(testStartingProfile,bool,true,"Calculate forward (etc) scores for the pre-training profile");

 /**
  * Write forward (etc) scores for the pre-training profile to STDOUT
  * during training?
  */
GALOSH_DEF_OPT(coutStartingProfile,bool,true,"Write forward (etc) scores for the pre-training profile to STDOUT during training?");

 /**
  * Also train the unconditional profile, and calculate forward (etc)
  * scores using it?
  */
GALOSH_DEF_OPT(testUnconditionalProfile,bool,true,"Also train the unconditional profile, and calculate forward (etc) scores using it?");

 /**
  * Write the forward (etc) scores for the unconditional profile to STDOUT
  * during training?
  */
GALOSH_DEF_OPT(coutUnconditionalProfile,bool,true,"Write the forward (etc) scores for the unconditional profile to STDOUT during training?");

 /**
  * Also train the unconditional (with fixed starting globals) profile,
  * and calculate forward (etc) scores using it?
  */
GALOSH_DEF_OPT(testUnconditionalWithFixedStartingGlobalsProfile,bool,true,"Also train the unconditional (with fixed starting globals) profile and calculate forward (etc) scores using it?");

 /**
  * Write the forward (etc) scores for the unconditional (with fixed
  * starting globals) profile to STDOUT during training?
  */
GALOSH_DEF_OPT(coutUnconditionalWithFixedStartingGlobalsProfile,bool,false,"Write the forward (etc) scores for the unconditional (with starting globals) profile to STDOUT during training?");

 /**
  * Also train the unconditional (with fixed true globals) profile,
  * and calculate forward (etc) scores using it?
  */
GALOSH_DEF_OPT(testUnconditionalWithFixedTrueGlobalsProfile,bool,true,"Train the unconditional (with fixed true globals) profile, and calculate forward (etc) scores using it?");

 /**
  * Write the forward (etc) scores for the unconditional (with fixed
  * true globals) profile to STDOUT during training?
  */
GALOSH_DEF_OPT(coutUnconditionalWithFixedTrueGlobalsProfile,bool,false,"Write the forward (etc) scores for the unconditional (with fixed true globals) profile to STDOUT during training?");

 /**
  * Also train the "conditional, then unconditional" profile, and
  * calculate forward (etc) scores using it?
  */
GALOSH_DEF_OPT(testConditionalThenUnconditionalProfile,bool,true,"train the \"conditional, then unconditional\" profile, and calculate forward (etc) scores using it?");

 /**
  * Write the forward (etc) scores for the "conditional, then
  * unconditional" profile to STDOUT during training?
  */
GALOSH_DEF_OPT(coutConditionalThenUnconditionalProfile,bool,false,"Write the forward (etc) scores for the \"conditional, then unconditional\" profile to STDOUT during training?");

 /**
  * Also train the "unconditional, then conditional" profile, and
  * calculate forward (etc) scores using it?
  */
GALOSH_DEF_OPT(testUnconditionalThenConditionalProfile,bool,true,"Train the \"unconditional, then conditional\" profile, and calculate forward (etc) scores using it?");

 /**
  * Write the forward (etc) scores for the "unconditional, then
  * conditional" profile to STDOUT during training?
  */
GALOSH_DEF_OPT(coutUnconditionalThenConditionalProfile,bool,false,"Write the forward (etc) scores for the \"unconditional, then conditional\" profile to STDOUT during training?");

 /**
  * Also train the "unconditional (first with fixed starting globals, then
  * with fixed positions)" profile, and calculate forward (etc) scores
  * using it?
  */
GALOSH_DEF_OPT(testUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile,bool,true,"Train the \"unconditional (first with fixed starting globals, then with fixed positions)\" profile, and calculate forward (etc) scores using it?");

 /**
  * Write the forward (etc) scores for the unconditional (first with fixed
  * starting globals, then with fixed positions) profile to STDOUT during
  * training?
  */
GALOSH_DEF_OPT(coutUnconditionalWithFixedStartingGlobalsThenWithFixedPositionsProfile,bool,false,"Write the forward (etc) scores for the unconditional (first with fixed starting globals, then with fixed positions) profile to STDOUT during training?");

 /**
  * Also train the "unconditional (first with fixed true globals, then
  * with fixed positions)" profile, and calculate forward (etc) scores
  * using it?
  */
GALOSH_DEF_OPT(testUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile,bool,true,"Train the \"unconditional (first with fixed true globals, then with fixed positions)\" profile, and calculate forward (etc) scores using it?");

 /**
  * Write the forward (etc) scores for the unconditional (first with fixed
  * true globals, then with fixed positions) profile to STDOUT during
  * training?
  */
GALOSH_DEF_OPT(coutUnconditionalWithFixedTrueGlobalsThenWithFixedPositionsProfile,bool,false,"Write the forward (etc) scores for the unconditional (first with fixed true globals, then with fixed positions) profile to STDOUT during training?");

 /**
  * Also use the best profile found after Gibbs sampling using Conditional
  * ("Per Position") Gibbs and calculate forward (etc) scores using it?
  */
GALOSH_DEF_OPT(testConditionalGibbsProfile,bool,true,"Use the best profile found after Gibbs sampling using Conditional (\"Per Position\") Gibbs and calculate forward (etc) scores using it?");

 /**
  * Write the forward (etc) scores for the best profile found (for
  * Conditional Gibbs) to STDOUT during training?
  */
GALOSH_DEF_OPT(coutConditionalGibbsProfile,bool,true,"Write the forward (etc) scores for the best profile found (for Conditional Gibbs) to STDOUT during training?");

 /**
  * Also use the best profile found after Gibbs sampling using
  * Unconditional ("Simple") Gibbs and calculate forward (etc) scores
  * using it?
  */
GALOSH_DEF_OPT(testUnconditionalGibbsProfile,bool,true,"Use the best profile found after Gibbs sampling using Unconditional (\"Simple\") Gibbs and calculate forward (etc) scores using it?");

 /**
  * Write the forward (etc) scores for the best profile found profile (for
  * Unconditional Gibbs) to STDOUT during training?
  */
GALOSH_DEF_OPT(coutUnconditionalGibbsProfile,bool,true,"Write the forward (etc) scores for the best profile found profile (for Unconditional Gibbs) to STDOUT during training?");

 /**
  * Also use lengthadjust to train the profiles, and calculate forward (etc) scores using it?
  */
GALOSH_DEF_OPT(testLengthadjust,bool,true,"Use lengthadjust to train the profiles, and calculate forward (etc) scores using it?");

 /**
  * Also use Baldi to train the profiles, and calculate forward (etc) scores using it?
  */
GALOSH_DEF_OPT(testBaldi,bool,false,"Use Baldi to train the profiles, and calculate forward (etc) scores using it?");

 /**
  * Also use Baldi / Siegel to train the profiles, and calculate forward (etc) scores using it?
  */
GALOSH_DEF_OPT(testBaldiSiegel,bool,true,"Use Baldi / Siegel to train the profiles, and calculate forward (etc) scores using it?");

 /**
  * If startWithUniformPositions or startWithPositionsDrawnFromPrior are
  * true, then under normal circumstances numStartingProfiles profiles
  * would be randomly generated; if alsoStartWithEvenPositions is *also*
  * true, then the first starting profile (index 0) will have even
  * positions rather than randomly-generated position emission values.
  */
GALOSH_DEF_OPT(alsoStartWithEvenPositions,bool,false,"If startWithUniformPositions or startWithPositionsDrawnFromPrior are true, then under normal circumstances numStartingProfiles profiles would be randomly generated; if alsoStartWithEvenPositions is *also* true, then the first starting profile (index 0) will have even positions rather than randomly-generated position emission values.");
// TODO: REMOVE
GALOSH_DEF_OPT(useRabinerScaling,bool,false,"Deprecated?");

/// This is how we're handling vectors.  It is a work-around because vectors are handled specially
/// by boost::program_options.  It allows the command line to look something like
///
/// --profileLengths 10 20 30
///
/// The TMP_EXTRA_STUFF must be set to include (at least) the ->multitoken() thing.
/// It should also be unset at the bottom of the vector initializations.

#undef TMP_EXTRA_STUFF
#define TMP_EXTRA_STUFF ->multitoken()

GALOSH_DEF_OPT(profileLengths,myVector<int>,myVector<int>(1,100) BOOST_PP_COMMA() string("100"),"Lengths of the profiles/fasta seqs");

/**
  * Use this number of sequences when training.
  *
  * UPDATE: This is now a pointer to a vector.  Tests will be run foreach
  * num_training_sequences_per_profile_i in
  * numTrainingSequencesPerProfiles.  If it is NULL, { 10, 100 } will be
  * used (this is the default).
  */

GALOSH_DEF_OPT(numTrainingSequencesPerProfiles,myVector<uint32_t>,myVector<uint32_t>(1,100) BOOST_PP_COMMA() string("100"),"Number of training sequences for each profile");

#ifdef PROFUSETEST_DEFAULT_TMP_ARRAY_TO_VECTOR
myVector<double> tmp_default_conservation_rates; for(double x=0.1; x<=1.0; x+=0.1){tmp_default_conservation_rates.push_back(x);}
#endif
/**
 * When making the true root profile from the pattern sequence, use this
 * probability for the pattern sequence base at each position, and divide
 * the remaining probability evenly among the remaining bases.
 *
 * UPDATE: This is now a pointer to a vector of rates.  Tests will be run
 * foreach conservation_rate in conservationRates.  If it is NULL,
 * { .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0 }
 * will be used (this is the default).
 *
 * \note Note: It is hard to initialize a vector with different values.  This
 * might become simpler with C++11 or with boost::array stuff.  Currently we're
 * doing it as above, and initializing with .begin() and .end() vector iterators
 *
 */
GALOSH_DEF_OPT(conservationRates,myVector<double>,myVector<double>(tmp_default_conservation_rates) BOOST_PP_COMMA() string("0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0"),"Iterate through each of these conservation rates");

/// Make the expected number of deletions be .5 or 1.0 per sequence.
/// Note that this will apply to the insertions, too, unless
/// m_parameters.useDeletionsForInsertionsParameters is set to false.
GALOSH_DEF_OPT(expectedDeletionsCounts,myVector<double>,myVector<double>(1,1.0) BOOST_PP_COMMA() string("0.5"),"Iterate through this set of expected deletions.  Note that deletions-counts=insertion counts unless useDeletionsForInsertionsParameters is false");

/**
  * If useDeletionsForInsertionsParameters is false, the insertionOpen
  * value of the true profile will be set to ( expectedInsertionsCount /
  * profileLength ).
  *
  * UPDATE: This is now a pointer to a vector.  Tests will be run foreach
  * expected_insertions_count in expectedInsertionCounts.  If it is unspecified,
  * { 1.0 } will be used.
  * @see useDeletionsForInsertionsParameters
  */
GALOSH_DEF_OPT(expectedInsertionsCounts,myVector<double>,myVector<double>(1,0.5) BOOST_PP_COMMA() string("0.5"),"Iterate through this series of expected insertions.  useDeletionsForInsertionsParameters must be false or this parameter is ignored.");

/**
   * The deletionExtension value of the true profile will be the minimum of
   * ( 1.0 / ( expectedDeletionLengthAsProfileLengthFraction *
   * profileLength ) ) and ( 1.0 / minExpectedDeletionLength ).  If
   * useDeletionsForInsertionsParameters is true, the insertionExtension
   * value of the true profile will also be set to be the minimum of ( 1.0
   * / ( expectedDeletionLengthAsProfileLengthFraction * profileLength ) )
   * and ( 1.0 / minExpectedDeletionLength ).
   *
   * UPDATE: This is now a pointer to a vector.  Tests will be run foreach
   * expected_deletion_length_as_profile_length_fraction in
   * expectedDeletionLengthAsProfileLengthFraction.  If it is NULL, { 0.1 }
   * will be used (this is the default).
   *
   * @see useDeletionsForInsertionsParameters
   */

// Make the expected length of each deletion be ( profile_length / 20 ) or (
// profile_length / 10 )...
// Note that this will apply to the insertions, too, unless
// m_parameters.useDeletionsForInsertionsParameters is set to false.
GALOSH_DEF_OPT(expectedDeletionLengthAsProfileLengthFractions,myVector<double>,myVector<double>(1,0.0125) BOOST_PP_COMMA() string("0.0125"),"Expected lengths of deletions. Iterate through all lengths in this list.");

// Make the expected length of each insertion be ( profile_length / 20 ) or (
// profile_length / 10 )...
// Note that this is not used unless
// m_parameters.useDeletionsForInsertionsParameters is set to false.
// ..(or 1.25, whichever is larger).
/**
* If useDeletionsForInsertionsParameters is false, the
* insertionExtension value of the true profile will be the minimum of (
* 1.0 / ( expectedInsertionLengthAsProfileLengthFraction * profileLength
* ) ) and ( 1.0 / minExpectedInsertionLength ).
*
* UPDATE: This is now a pointer to a vector.  Tests will be run foreach
* expected_insertion_length_as_profile_length_fraction in
* expectedInsertionLengthAsProfileLengthFraction.  If it is NULL, { 0.1 }
* will be used (this is the default).
*
* @see useDeletionsForInsertionsParameters
*/
GALOSH_DEF_OPT(expectedInsertionLengthAsProfileLengthFractions,myVector<double>,myVector<double>(1,0.1) BOOST_PP_COMMA() string("0.0125"),"Expected lengths of deletions. Iterate through all lengths in this list.");


/** do this after the vector definition section */
#undef TMP_EXTRA_STUFF
#define TMP_EXTRA_STUFF BOOST_PP_EMPTY()
