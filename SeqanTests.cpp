#include "Ambiguous.hpp"
#include "Algebra.hpp"
#include "Sequence.hpp"
#include "MultinomialDistribution.hpp"
#include "ProfileHMM.hpp"
#include "Profile.hpp"
#include "Fasta.hpp"
#include "Random.hpp"
#include "DynamicProgramming.hpp"
#include "ProfileTrainer.hpp"
#include "ProfileTreeTrainer.hpp"
//#include "ProfileGibbs.hpp" // TODO: Add tests for gibbs...
#include "ProfuseTest.hpp"
#include "AminoAcid20.hpp"

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/find_motif.h>

#ifdef __HAVE_MUSCLE
#include "muscle/distfunc.h"
#include "muscle/clustsetdf.h"
#include "muscle/clust.h"
#include "muscle/tree.h"
#include "muscle/textfile.h"
#endif // __HAVE_MUSCLE

#include <boost/lexical_cast.hpp>

#ifdef __HAVE_MUSCLE
int g_argc;
char **g_argv;
#endif // __HAVE_MUSCLE

using namespace seqan;

// Uncomment this to switch to testing the simple basic stuff (only)
#define TEST_BASICS

#ifdef TEST_BASICS
template <typename T>
static bool isTrue ( T const & tag ) { return false; }
static bool isTrue ( seqan::True const & tag ) { return true; }

void
test_basics ()
{
  // All basic (non-dp) tests can be turned off by setting test_non_dp to false.
  const bool test_non_dp = true;
  const bool test_alphabets = true;
  const bool test_ambiguous = true; // Note this can be used in conjunction with test_multinomials and test_trainer too
  const bool test_sequences = true;
  const bool test_fasta = true;
  const bool test_multinomials = true; // see also test_ambiguous
  const bool test_profile_hmm_states = true;
  const bool test_profiles = true;

  // All dp tests can be turned off using test_dynamic_programming = false.
  const bool test_dynamic_programming = true;
  const bool test_random_sequences = false;
  const bool test_dp_basics = false;
  const bool test_alignment_profiles = true;//false;
  const bool test_priors = false;
#ifdef ALLOW_BOLTZMANN_GIBBS
  const bool test_boltzmann_gibbs = true;
#endif // ALLOW_BOLTZMANN_GIBBS
  const bool test_cross_entropy = false;
  const bool test_profile_profile_alignment = false;
  const bool test_path_cooccurrence = false;

  const bool test_trainer = false;//true; // see also test_ambiguous

///The typical alphabet is convertible to $char$.
///Note that a conversion of a $char$ into another alphabet and back can change the value of the $char$.
  if( test_non_dp ) {
    if( test_alphabets ) {
      Dna a = 'a';
      std::cout << a; //output: 'A'
  
      Dna5 b = 'f'; //'f' is unknown character
      std::cout << b; //output: 'N'
  
      ///Many SeqAn alphabet classes can be converted into each other.
      b = a;
      std::cout << b; //output: 'A'
  
      Iupac c = b;
      std::cout << c; //output: 'A'
  
      std::cout << std::endl;
    } // End if test_alphabets

    if( test_ambiguous ) {
      std::cout << "Dna is " << ( isTrue( galosh::IsAmbiguous<seqan::Dna, seqan::Dna5>::Type() ) ? "" : "not " ) << "ambiguous over Dna5." << std::endl;
      std::cout << "Dna5 is " << ( isTrue( galosh::IsAmbiguous<seqan::Dna5, seqan::Dna>::Type() ) ? "" : "not " ) << "ambiguous over Dna." << std::endl;

      Dna dna_residue;
      Dna5 a = 'a';
      std::cout << "Dna5's element '" << a << "' is ambiguous over " << galosh::ambiguousCount( a, Dna() ) << " Dna elements." << std::endl;
      galosh::ambiguousAssign( dna_residue, a, 0 );
      std::cout << "A possible value of A: " << dna_residue << endl;
      Dna5 n = 'n';
      std::cout << "Dna5's element '" << n << "' is ambiguous over " << galosh::ambiguousCount( n, Dna() ) << " Dna elements." << std::endl;
      size_t num_elements = galosh::ambiguousCount( n, Dna() );
      for( int i = 0; i < num_elements; i++ ) {
        galosh::ambiguousAssign( dna_residue, n, i );
        std::cout << "A possible value of N: " << dna_residue << endl;
      }
    } // End if test_ambiguous
  
    if( test_sequences ) {
      seqan::Peptide prot = "anypeptide";
      std::cout << length(prot) << std::endl;  //output: 10
  
      std::cout << prot[ 9 ] << std::endl;
  
      prot += "anewend";
      std::cout << prot << std::endl;          //ouput: "ANYPEPTIDEANEWEND"
      std::cout << "\tlength( prot ) is " << seqan::length( prot ) << std::endl;
  
      /////
      galosh::Sequence<Dna> dna_seq = "acgt";
      std::cout << dna_seq << std::endl;
      std::cout << "\tlength( dna_seq ) is " << seqan::length( dna_seq ) << std::endl;
    } // End if test_sequences
  
    if( test_fasta ) {
      ///////
      //// Fasta with no gaps
      galosh::Fasta<Dna> fasta;
      fasta.fromFile( "fasta/dna.10.10.fasta" );
      std::cout << fasta << std::endl;
  
      ///////
      //// Fasta: with gaps.  Use "char" since basic seqan alphabets don't support gap chars.
      galosh::Fasta<char> aligned_fasta;
      aligned_fasta.fromFile( "fasta/21U.fa.muscle.first10" );
      std::cout << aligned_fasta << std::endl;

      // TODO:  I'm thinking of just scrapping Fasta altogether, using
      // StringSet instead, with _loadSequences from
      // seqan/graph_utils/utility_parsing.h, as in seqan_tcoffee.
      // Eventually.
    } // End if test_fasta
  
    if( test_multinomials ) {
      Dna c_2( 1 );
      std::cout << c_2 << std::endl; //output: 'C'
  
      galosh::MultinomialDistribution<Dna, realspace> dna_dist;
      dna_dist[ c_2 ] = .4;
      std::cout << dna_dist << std::endl;
  
      galosh::MultinomialDistribution<galosh::StateLabel, float> state_dist;
      std::cout << state_dist << std::endl;

      if( test_ambiguous ) {
        Dna5 a = 'a';
        Dna5 n = 'n';
        std::cout << "The dna_dist, accessed using Dna5 'a', returns " << dna_dist.ambiguousSum( a ) << endl;
        //std::cout << "The dna_dist, accessed using Dna5 'n', returns " << dna_dist.ambiguousSum( n ) << endl;
        std::cout << "The dna_dist, accessed using Dna5 'n', returns " << dna_dist[ n ] << endl;

        // Uncomment this, and it should fail to compile, since AminoAcid is
        // not ambiguous over Dna.
        //AminoAcid aminoacid_n = 'n';
        //std::cout << "The dna_dist, accessed using AminoAcid 'n', returns " << dna_dist[ aminoacid_n ] << endl;
        
        dna_dist.ambiguousIncrement( n, 1.0 );
        std::cout << "After performing an ambiguousIncrement( n, 1.0 ), the dna_dist is " << dna_dist << endl;
        galosh::MultinomialDistribution<Dna, realspace>::AmbiguousValue<Dna5> ambiguous_value = dna_dist[ n ];
        ambiguous_value += 1.0;
        std::cout << "After performing an ( ambiguous_value = dna_dist[ n ] ) += 1.0, the dna_dist is " << dna_dist << ", and the ambiguous_value is " << ambiguous_value << endl;
        ambiguous_value = 1.0;
        std::cout << "After performing an ambiguous_value = 1.0, the dna_dist is " << dna_dist << ", and the ambiguous_value is " << ambiguous_value << endl;
        dna_dist[ n ] += 1.0;
        std::cout << "After performing a dna_dist[ n ] += 1.0, the dna_dist is " << dna_dist << endl;
      } // End if test_ambiguous (and test_multinomials)
  
      // TODO: Use seqan::FrequencyDistribution?
      //FrequencyDistribution<Dna, float> dna_dist_seqan;
      //dna_dist_seqan[ a ] = .5;
      //std::cout << c_2 << endl;
      //std::cout << ordValue( c_2 ) << endl;
      //dna_dist_seqan[ c_2 ] = .5;
      //dna_dist_seqan[ Dna( 'g' ) ] = .5;
      //dna_dist_seqan[ Dna( 'T' ) ] = .5;
      //normalize( dna_dist_seqan );
      //std::cout << dna_dist_seqan << std::endl;
      //// It looks like I could replace my MultinomialDistribution, or make it
      //// derive from FrequencyDistribution.
      //// Note also their "profile" is a StringSet of FrequencyDistributions.
    } // End if test_multinomials
  
    if( test_profile_hmm_states ) {
      ///////
      // ProfileHMM states
      // Start state
      std::cout << "StateLabelId<StartStateLabel> is " << galosh::StateLabelId<galosh::StartStateLabel>::VALUE << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::StartStateLabel>::VALUE ) << " is " << ( isTrue( IsSimple<galosh::StartStateLabel>::Type() ) ? "" : "not " ) << "simple." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::StartStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsEmitting<galosh::StartStateLabel>::Type() ) ? "" : "not " ) << "emitting." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::StartStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsAssociatedWithPosition<galosh::StartStateLabel>::Type() ) ? "" : "not " ) << "associated with an ancestral sequence position." << endl;
      galosh::MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::StartStateLabel, galosh::Plan7>::Type, float> Start_state_dist;
      std::cout << Start_state_dist << std::endl;
  
      // PreAlign state
      std::cout << "StateLabelId<PreAlignStateLabel> is " << galosh::StateLabelId<galosh::PreAlignStateLabel>::VALUE << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::PreAlignStateLabel>::VALUE ) << " is " << ( isTrue( IsSimple<galosh::PreAlignStateLabel>::Type() ) ? "" : "not " ) << "simple." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::PreAlignStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsEmitting<galosh::PreAlignStateLabel>::Type() ) ? "" : "not " ) << "emitting." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::PreAlignStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsAssociatedWithPosition<galosh::PreAlignStateLabel>::Type() ) ? "" : "not " ) << "associated with an ancestral sequence position." << endl;
      galosh::MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::PreAlignStateLabel, galosh::Plan7>::Type, float> PreAlign_state_dist;
      std::cout << PreAlign_state_dist << std::endl;
  
      // Begin state
      std::cout << "StateLabelId<BeginStateLabel> is " << galosh::StateLabelId<galosh::BeginStateLabel>::VALUE << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::BeginStateLabel>::VALUE ) << " is " << ( isTrue( IsSimple<galosh::BeginStateLabel>::Type() ) ? "" : "not " ) << "simple." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::BeginStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsEmitting<galosh::BeginStateLabel>::Type() ) ? "" : "not " ) << "emitting." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::BeginStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsAssociatedWithPosition<galosh::BeginStateLabel>::Type() ) ? "" : "not " ) << "associated with an ancestral sequence position." << endl;
      galosh::MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::BeginStateLabel, galosh::Plan7>::Type, float> Begin_state_dist;
      std::cout << Begin_state_dist << std::endl;
  
      // Match state
      std::cout << "StateLabelId<MatchStateLabel> is " << galosh::StateLabelId<galosh::MatchStateLabel>::VALUE << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::MatchStateLabel>::VALUE ) << " is " << ( isTrue( IsSimple<galosh::MatchStateLabel>::Type() ) ? "" : "not " ) << "simple." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::MatchStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsEmitting<galosh::MatchStateLabel>::Type() ) ? "" : "not " ) << "emitting." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::MatchStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsAssociatedWithPosition<galosh::MatchStateLabel>::Type() ) ? "" : "not " ) << "associated with an ancestral sequence position." << endl;
      galosh::MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, float> Match_state_dist;
      std::cout << Match_state_dist << std::endl;
  
      // Insertion state (both Plan 7 and Plan 9)
      std::cout << "StateLabelId<InsertionStateLabel> is " << galosh::StateLabelId<galosh::InsertionStateLabel>::VALUE << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::InsertionStateLabel>::VALUE ) << " is " << ( isTrue( IsSimple<galosh::InsertionStateLabel>::Type() ) ? "" : "not " ) << "simple." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::InsertionStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsEmitting<galosh::InsertionStateLabel>::Type() ) ? "" : "not " ) << "emitting." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::InsertionStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsAssociatedWithPosition<galosh::InsertionStateLabel>::Type() ) ? "" : "not " ) << "associated with an ancestral sequence position." << endl;
      galosh::MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::InsertionStateLabel, galosh::Plan9>::Type, float> Plan9_Insertion_state_dist;
      std::cout << Plan9_Insertion_state_dist << std::endl;
      galosh::MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::InsertionStateLabel, galosh::Plan7>::Type, float> Plan7_Insertion_state_dist;
      std::cout << Plan7_Insertion_state_dist << std::endl;
  
      // Deletion state (both Plan 7 and Plan 9)
      std::cout << "StateLabelId<DeletionStateLabel> is " << galosh::StateLabelId<galosh::DeletionStateLabel>::VALUE << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::DeletionStateLabel>::VALUE ) << " is " << ( isTrue( IsSimple<galosh::DeletionStateLabel>::Type() ) ? "" : "not " ) << "simple." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::DeletionStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsEmitting<galosh::DeletionStateLabel>::Type() ) ? "" : "not " ) << "emitting." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::DeletionStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsAssociatedWithPosition<galosh::DeletionStateLabel>::Type() ) ? "" : "not " ) << "associated with an ancestral sequence position." << endl;
      galosh::MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::DeletionStateLabel, galosh::Plan9>::Type, float> Plan9_Deletion_state_dist;
      std::cout << Plan9_Deletion_state_dist << std::endl;
      galosh::MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::DeletionStateLabel, galosh::Plan7>::Type, float> Plan7_Deletion_state_dist;
      std::cout << Plan7_Deletion_state_dist << std::endl;
  
      // End state
      std::cout << "StateLabelId<EndStateLabel> is " << galosh::StateLabelId<galosh::EndStateLabel>::VALUE << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::EndStateLabel>::VALUE ) << " is " << ( isTrue( IsSimple<galosh::EndStateLabel>::Type() ) ? "" : "not " ) << "simple." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::EndStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsEmitting<galosh::EndStateLabel>::Type() ) ? "" : "not " ) << "emitting." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::EndStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsAssociatedWithPosition<galosh::EndStateLabel>::Type() ) ? "" : "not " ) << "associated with an ancestral sequence position." << endl;
      galosh::MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::EndStateLabel, galosh::Plan7>::Type, float> End_state_dist;
      std::cout << End_state_dist << std::endl;
  
      // Loop state
      std::cout << "StateLabelId<LoopStateLabel> is " << galosh::StateLabelId<galosh::LoopStateLabel>::VALUE << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::LoopStateLabel>::VALUE ) << " is " << ( isTrue( IsSimple<galosh::LoopStateLabel>::Type() ) ? "" : "not " ) << "simple." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::LoopStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsEmitting<galosh::LoopStateLabel>::Type() ) ? "" : "not " ) << "emitting." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::LoopStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsAssociatedWithPosition<galosh::LoopStateLabel>::Type() ) ? "" : "not " ) << "associated with an ancestral sequence position." << endl;
      galosh::MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::LoopStateLabel, galosh::Plan7>::Type, float> Loop_state_dist;
      std::cout << Loop_state_dist << std::endl;
  
      // PostAlign state
      std::cout << "StateLabelId<PostAlignStateLabel> is " << galosh::StateLabelId<galosh::PostAlignStateLabel>::VALUE << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::PostAlignStateLabel>::VALUE ) << " is " << ( isTrue( IsSimple<galosh::PostAlignStateLabel>::Type() ) ? "" : "not " ) << "simple." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::PostAlignStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsEmitting<galosh::PostAlignStateLabel>::Type() ) ? "" : "not " ) << "emitting." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::PostAlignStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsAssociatedWithPosition<galosh::PostAlignStateLabel>::Type() ) ? "" : "not " ) << "associated with an ancestral sequence position." << endl;
      galosh::MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::PostAlignStateLabel, galosh::Plan7>::Type, float> PostAlign_state_dist;
      std::cout << PostAlign_state_dist << std::endl;
  
      // Terminal state
      std::cout << "StateLabelId<TerminalStateLabel> is " << galosh::StateLabelId<galosh::TerminalStateLabel>::VALUE << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::TerminalStateLabel>::VALUE ) << " is " << ( isTrue( IsSimple<galosh::TerminalStateLabel>::Type() ) ? "" : "not " ) << "simple." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::TerminalStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsEmitting<galosh::TerminalStateLabel>::Type() ) ? "" : "not " ) << "emitting." << endl;
      std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::TerminalStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsAssociatedWithPosition<galosh::TerminalStateLabel>::Type() ) ? "" : "not " ) << "associated with an ancestral sequence position." << endl;
  
      // End ProfileHMM states
      ///////
    } // End if test_profile_hmm_states
  } // End if test_non_dp

  if(
    ( test_non_dp && test_profiles ) ||
    test_dynamic_programming ||
    test_trainer
  ) {
    //typedef float ProbabilityType;
    //typedef float ScoreType;
    //typedef float MatrixValueType;
    typedef realspace ProbabilityType;
    typedef realspace ScoreType;
    typedef realspace MatrixValueType;
    typedef ProfileTreeRoot<Dna, ProbabilityType> ProfileType;

    Dna residue;
    galosh::Sequence<Dna5> dna_seq_a = "a";
    galosh::Sequence<Dna5> dna_seq_b = "ag";
    galosh::Sequence<Dna5> dna_seq_c = "aga";
    vector<galosh::Sequence<Dna5> > dna_seqs( 3 );
    dna_seqs[ 0 ] = dna_seq_a;
    dna_seqs[ 1 ] = dna_seq_b;
    dna_seqs[ 2 ] = dna_seq_c;
    int num_sequences_to_use = dna_seqs.size();
  
#ifdef USE_DEL_IN_DEL_OUT
    // We need at least 3 positions to be able to test extensions of del-ins and del-outs
    galosh::ProfileTreeRoot<Dna, ProbabilityType> dna_profile( 3 );
#else
    galosh::ProfileTreeRoot<Dna, ProbabilityType> dna_profile( 2 );
#endif // USE_DEL_IN_DEL_OUT .. else ..
    galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer forward_matrices(
      dna_profile,
      dna_seqs,
      num_sequences_to_use
    );
    galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer::iterator forward_matrices_iterator =
      forward_matrices.begin();
    galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer backward_matrices(
      dna_profile,
      dna_seqs,
      num_sequences_to_use
    );
    galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer::reverse_iterator backward_matrices_iterator =
      backward_matrices.rbegin();
    galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::PositionSpecificSequenceScoreCoefficientsVector coefficients_vector( 3 );
    galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::PositionEntente position_entente;
    galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::PositionEntente position_entente_unscaled;
    galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::PositionEntente position_entente_backup;
 
    ///////
    /// Profile
    if( test_non_dp && test_profiles ) {
      std::cout << "The dna profile is:" << std::endl;
      std::cout << dna_profile;
      std::cout << endl;

      std::cout << "A default Iupac profile is:" << std::endl;
      galosh::ProfileTreeRoot<Iupac, ProbabilityType> iupac_profile;
      std::cout << iupac_profile;
      std::cout << endl;

      std::cout << "Reading profile from file 'seqantest.DATA.profile'" << std::endl;
      galosh::ProfileTreeRoot<Dna, ProbabilityType> dna_profile_from_file;
      dna_profile_from_file.fromFile( "seqantest.DATA.profile" );
      std::cout << "\tgot:" << std::endl;
      std::cout << dna_profile_from_file;
      std::cout << endl;
    }

#ifndef DISALLOW_FLANKING_TRANSITIONS
    dna_profile[ galosh::Transition::fromPreAlign ][ galosh::TransitionFromPreAlign::toPreAlign ] =
      .1;
    dna_profile[ galosh::Transition::fromPreAlign ][ galosh::TransitionFromPreAlign::toBegin ] =
      1.0 -
      dna_profile[ galosh::Transition::fromPreAlign ][ galosh::TransitionFromPreAlign::toPreAlign ];
#endif // !DISALLOW_FLANKING_TRANSITIONS

    dna_profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletion ] =
      .01;
#ifdef USE_DEL_IN_DEL_OUT
    dna_profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletionIn ] =
      .5;
    dna_profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toMatch ] =
      1.0 -
      (
        dna_profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletion ] +
        dna_profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletionIn ]
      );
#ifndef USE_SWENTRY_SWEXIT
    dna_profile[ galosh::Transition::fromDeletionIn ][ galosh::TransitionFromDeletionIn::toDeletionIn ] =
      .5;
    dna_profile[ galosh::Transition::fromDeletionIn ][ galosh::TransitionFromDeletionIn::toMatch ] =
      1.0 -
      dna_profile[ galosh::Transition::fromDeletionIn ][ galosh::TransitionFromDeletionIn::toDeletionIn ];
#endif // !USE_SWENTRY_SWEXIT
#else
    dna_profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toMatch ] =
      1.0 -
      dna_profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletion ];
#endif // USE_DEL_IN_DEL_OUT .. else ..
    
    dna_profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toInsertion ] =
      .02;
    dna_profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletion ] =
      .01;
#ifdef USE_DEL_IN_DEL_OUT
    dna_profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletionOut ] =
      .5;
    dna_profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toMatch ] =
      1.0 -
      (
        dna_profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toInsertion ] +
        dna_profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletion ] +
        dna_profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletionOut ]
      );
#ifndef USE_SWENTRY_SWEXIT
    dna_profile[ galosh::Transition::fromDeletionOut ][ galosh::TransitionFromDeletionOut::toDeletionOut ] =
      .5;
    dna_profile[ galosh::Transition::fromDeletionOut ][ galosh::TransitionFromDeletionOut::toEnd ] =
      1.0 -
      dna_profile[ galosh::Transition::fromDeletionOut ][ galosh::TransitionFromDeletionOut::toDeletionOut ];
#endif // !USE_SWENTRY_SWEXIT
#else
    dna_profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toMatch ] =
      1.0 -
      (
        dna_profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toInsertion ] +
        dna_profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletion ]
      );
#endif // USE_DEL_IN_DEL_OUT .. else ..
    dna_profile[ 0 ][ galosh::Transition::fromInsertion ][ galosh::TransitionFromInsertion::toInsertion ] =
      .1;
    dna_profile[ 0 ][ galosh::Transition::fromInsertion ][ galosh::TransitionFromInsertion::toMatch ] =
      1.0 -
      dna_profile[ 0 ][ galosh::Transition::fromInsertion ][ galosh::TransitionFromInsertion::toInsertion ];
    dna_profile[ 0 ][ galosh::Transition::fromDeletion ][ galosh::TransitionFromDeletion::toDeletion ] =
      .2;
    dna_profile[ 0 ][ galosh::Transition::fromDeletion ][ galosh::TransitionFromDeletion::toMatch ] =
      1.0 -
      dna_profile[ 0 ][ galosh::Transition::fromDeletion ][ galosh::TransitionFromDeletion::toDeletion ];
#ifndef DISALLOW_FLANKING_TRANSITIONS
    dna_profile[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toPostAlign ] =
      .1;
    dna_profile[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toTerminal ] =
      1.0 -
      dna_profile[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toPostAlign ];
#endif // !DISALLOW_FLANKING_TRANSITIONS
    std::cout << "Profile is:" << std::endl;
    dna_profile.writeExceptPositions( std::cout );

    dna_profile[ 0 ][ galosh::Emission::Match ][ ( residue = 'a' ) ] =
      .7;
    dna_profile[ 0 ][ galosh::Emission::Match ][ ( residue = 'c' ) ] =
      .1;
    dna_profile[ 0 ][ galosh::Emission::Match ][ ( residue = 'g' ) ] =
      .1;
    dna_profile[ 0 ][ galosh::Emission::Match ][ ( residue = 't' ) ] =
      .1;
    dna_profile[ 1 ][ galosh::Emission::Match ][ ( residue = 'a' ) ] =
      .4;
    dna_profile[ 1 ][ galosh::Emission::Match ][ ( residue = 'c' ) ] =
      .1;
    dna_profile[ 1 ][ galosh::Emission::Match ][ ( residue = 'g' ) ] =
      .3;
    dna_profile[ 1 ][ galosh::Emission::Match ][ ( residue = 't' ) ] =
      .2;
#ifdef USE_DEL_IN_DEL_OUT
    dna_profile[ 2 ][ galosh::Emission::Match ][ ( residue = 'a' ) ] =
      .1;
    dna_profile[ 2 ][ galosh::Emission::Match ][ ( residue = 'c' ) ] =
      .1;
    dna_profile[ 2 ][ galosh::Emission::Match ][ ( residue = 'g' ) ] =
      .3;
    dna_profile[ 2 ][ galosh::Emission::Match ][ ( residue = 't' ) ] =
      .5;
#endif // USE_DEL_IN_DEL_OUT

    dna_profile.writePositions( std::cout );
    std::cout << endl;

    if( test_dynamic_programming || test_trainer ) {
      ///////
      //// DynamicProgramming
      galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType> dp;
          
      galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Parameters dp_parameters;
      Random random;
      ScoreType score;
      if( test_dynamic_programming ) {
        if( test_random_sequences ) {
          random.setSeed( static_cast<uint32_t>( std::time( NULL ) ) );
          cout << "After reseeding with seed " << random.getSeed() << ":" << endl;
          // Now try drawing from the dna_profile before training.
          galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::MultipleAlignment<ProfileType, Dna5> ma_before;

          // TODO: REMOVE .. just use the regular dna_profile..
          galosh::ProfileTreeRoot<Dna, ProbabilityType> another_dna_profile( 10 );
          another_dna_profile.copyExceptPositions( dna_profile );
          cout << "Drawing from profile:\n" << another_dna_profile << endl;

          galosh::Fasta<Dna5> random_fasta_before;
          uint32_t zero_length_sequences_drawn =
            dp.drawSequences(
              dp_parameters,
              //dna_profile,
              another_dna_profile,
              5,
              "Randomly generated ",
              random,
              random_fasta_before,
              ma_before
            );
          cout << "Randomly generated " << ( 5 + zero_length_sequences_drawn ) << ", of which " << zero_length_sequences_drawn << " had zero length." << endl;
          //cout << "Randomly generated sequences are:" << endl;
          //cout << random_fasta_before << endl;
          cout << "Alignment paths are:" << endl;
          cout << ma_before << endl;
          ma_before.toPairwiseStream( cout );
          
        } // End if test_random_sequences
    
        if( test_dp_basics || test_alignment_profiles || test_priors ) {
          if( test_dp_basics ) {
            dp_parameters.rabinerScaling_useMaximumValue = false;
            dp_parameters.matrixRowScaleFactor = 1;
            // Test using rabiner scaling at first, or not: (note it keeps getting changed -- scroll down)
            dp_parameters.useRabinerScaling = false;//true;

            ScoreType rabiner_cumulative_inverse_scalar = 1;
                
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::Row forward_row_a_1( 2 );
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::Row forward_row_a_2( 2 );
            // Keep the inverse scalars for future reference
#ifdef USE_DEL_IN_DEL_OUT
            vector<MatrixValueType> rabiner_inverse_scalars_a( 4 );
#else
            vector<MatrixValueType> rabiner_inverse_scalars_a( 3 );
#endif // USE_DEL_IN_DEL_OUT .. else ..
            dp.forward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_a,
              0,
              forward_row_a_2,
              forward_row_a_1
            );
                
            if( dp_parameters.useRabinerScaling ) {
              // Save the inverse scalar for use with the backward stuff..
              rabiner_inverse_scalars_a[ 0 ] = forward_row_a_1.m_rabinerInverseScalar;
              // Start over with the cumulative inverse scalar
              rabiner_cumulative_inverse_scalar = forward_row_a_1.m_rabinerInverseScalar;
              // With rabiner scaling:
              cout << "Rabiner Scaled:" << endl;
              cout << forward_row_a_1 << endl;
              //[ 0.97561, 0, 0 ], [ 0.0243902, 0, 0 ]
              // USE_DEL_IN_DEL_OUT:
              // 
              cout << "rabiner sum is " << forward_row_a_1.m_rabinerInverseScalar << endl;
              //1.025 = exp( 0.0246926 )
              cout << "rabiner cumulative inverse scalar is " << rabiner_cumulative_inverse_scalar << endl;
              //1.025 = exp( 0.0246926 )
              // Without rabiner scaling:
              forward_row_a_2 = forward_row_a_1;
              forward_row_a_2 *= rabiner_cumulative_inverse_scalar;
              cout << "Unscaled:" << endl;
              cout << forward_row_a_2 << endl;
              //[ 1, 0, 0 ], [ 0.025, 0, 0 ]
            } else {
              cout << forward_row_a_1 << endl;
              // Without rabiner scaling:
              //[ 1, 0, 0 ], [ 0.025, 0, 0 ]
              // USE_DEL_IN_DEL_OUT:
              //[ 1, 0, 0, 0, 0 ], [ 0.025, 0, 0, 0, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 1, 0, 0 ], [ 0, 0, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1, deletionIn: 0, deletionOut: 0
            }
            
            dp.forward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_a,
              1,
              forward_row_a_1,
              forward_row_a_2
            );
            if( dp_parameters.useRabinerScaling ) {
              // Save the inverse scalar for use with the backward stuff..
              rabiner_inverse_scalars_a[ 1 ] = forward_row_a_2.m_rabinerInverseScalar;
              rabiner_cumulative_inverse_scalar *= forward_row_a_2.m_rabinerInverseScalar;
              // With rabiner scaling:
              cout << "Rabiner Scaled:" << endl;
              cout << forward_row_a_2 << endl;
              //[ 0, 0, 0.0142197 ], [ 0.985425, 0, 0.000355492 ]
              cout << "rabiner sum is " << forward_row_a_2.m_rabinerInverseScalar << endl;
              //0.617488 = exp(-0.482096)
              cout << "rabiner cumulative inverse scalar is " << rabiner_cumulative_inverse_scalar << endl;
              //0.632925 = exp(-0.457403)
              // Without rabiner scaling:
              forward_row_a_1 = forward_row_a_2;
              forward_row_a_1 *= rabiner_cumulative_inverse_scalar;
              cout << "Unscaled:" << endl;
              cout << forward_row_a_1 << endl;
              // Without rabiner scaling:
              //[ 0, 0, 0.009 ], [ 0.6237, 0, 0.000225 ]
            } else {
              cout << forward_row_a_2 << endl;
              // Without rabiner scaling:
              //[ 0, 0, 0.009 ], [ 0.6237, 0, 0.000225 ]
              // USE_DEL_IN_DEL_OUT:
              //[ 0, 0, 0.009, 0.45, 0 ], [ 0.3087, 0, 0.000225, 0.01125, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 0, 0, 0.01 ], [ 0.343, 0, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1, deletionIn: 0.5, deletionOut: 0
            }
            
            dp.forward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_a,
              2,
              forward_row_a_2,
              forward_row_a_1
            );
            if( dp_parameters.useRabinerScaling ) {
              // Save the inverse scalar for use with the backward stuff..
              rabiner_inverse_scalars_a[ 2 ] = forward_row_a_1.m_rabinerInverseScalar;
              rabiner_cumulative_inverse_scalar *= forward_row_a_1.m_rabinerInverseScalar;
              // With rabiner scaling:
              cout << "Rabiner Scaled:" << endl;
              cout << forward_row_a_1 << endl;
              //[ 0.217226, 0, 0 ], [ 1.11111, 0, 0 ]
              cout << "rabiner sum is " << forward_row_a_1.m_rabinerInverseScalar << endl;
              //0.0130921 = exp(-4.33575)
              cout << "rabiner cumulative inverse scalar is " << rabiner_cumulative_inverse_scalar << endl;
              //0.008286314
              // Without rabiner scaling:
              forward_row_a_2 = forward_row_a_1;
              forward_row_a_2 *= rabiner_cumulative_inverse_scalar;
              cout << "Unscaled:" << endl;
              cout << forward_row_a_2 << endl;
              // Without rabiner scaling:
              //[ 0.0018, 0, 0 ], [ 0.009207, 0, 0 ]
            } else {
              cout << forward_row_a_1 << endl;
              // Without rabiner scaling:
              //[ 0.0018, 0, 0 ], [ 0.009207, 0, 0 ]
              // USE_DEL_IN_DEL_OUT:
              //[ 0, 0, 0.0018, 0.225, 0 ], [ 0.09288, 0, 0.00544725, 0.005625, 0.0385875 ], inverseScalar: 1, cumulativeInverseScalar: 1
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 0, 0, 0.002 ], [ 0.1032, 0, 0.0060025 ], inverseScalar: 1, cumulativeInverseScalar: 1, deletionIn: 0.25, deletionOut: 0.042875
            }
            
#ifdef USE_DEL_IN_DEL_OUT
            dp.forward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_a,
              3,
              forward_row_a_1,
              forward_row_a_2
            );
            if( dp_parameters.useRabinerScaling ) {
              // Save the inverse scalar for use with the backward stuff..
              rabiner_inverse_scalars_a[ 3 ] = forward_row_a_2.m_rabinerInverseScalar;
              rabiner_cumulative_inverse_scalar *= forward_row_a_2.m_rabinerInverseScalar;
              // With rabiner scaling:
              cout << "Rabiner Scaled:" << endl;
              cout << forward_row_a_2 << endl;

              cout << "rabiner sum is " << forward_row_a_2.m_rabinerInverseScalar << endl;

              cout << "rabiner cumulative inverse scalar is " << rabiner_cumulative_inverse_scalar << endl;

              // Without rabiner scaling:
              forward_row_a_1 = forward_row_a_2;
              forward_row_a_1 *= rabiner_cumulative_inverse_scalar;
              cout << "Unscaled:" << endl;
              cout << forward_row_a_1 << endl;
            } else {
              cout << forward_row_a_2 << endl;
              // Without rabiner scaling:
              // USE_DEL_IN_DEL_OUT:
              //[ 0.00036, 0, 0, 0, 0 ], [ 0.01388565, 0, 0, 0, 0.0618075 ], inverseScalar: 1, cumulativeInverseScalar: 1
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:0.013176
              //[ 0.0004, 0, 0 ], [ 0.0154085, 0, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1, deletionIn: 0, deletionOut: 0.068675
            }
#endif // USE_DEL_IN_DEL_OUT

//[ 1, 0, 0 ], [ 0, 0, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1, deletionIn: 0, deletionOut: 0
//[ 0, 0, 0.01 ], [ 0.343, 0, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1, deletionIn: 0.5, deletionOut: 0
//[ 0, 0, 0.002 ], [ 0.1032, 0, 0.0060025 ], inverseScalar: 1, cumulativeInverseScalar: 1, deletionIn: 0.25, deletionOut: 0.042875
//[ 0.0004, 0, 0 ], [ 0.0154085, 0, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1, deletionIn: 0, deletionOut: 0.068675

// .0060025 * .2 = .0012005
// .1032 * (.01*(1-.25)/(1-.5)) = .1032 * .015 = .001548
//.0154085 - ( .0012005 + .001548 ) = .0154085 - .0027485 = .01266
// .01266 / .1 = .1266
// .002 * .8 = .0016
// .1266 - .0016 = .125
// .125 / .25 = .5 (closing the del-in).
// col 1 Z->M->C Match contribution: 0.0125
// col 1 C->C Match contribution: 0
// col 1 M->M->C Match contribution: 0
// using prev_row_match_distribution == (M=0.705,I=0.03,D=0.015,W=0.25)
// col 1 M->D->C Match contribution: 0.001548
// col 1 I->M->C Match contribution: 0
// col 1 D->D->C Match contribution: 0.0012005
// col 1 D->M->C Match contribution: 0.00016

            // TODO: REMOVE
            //exit( 0 );

            bool forward_swap =
              dp.forward_score( 
                dp_parameters,
                dna_profile,
                dna_seq_a,
                forward_row_a_1,
                forward_row_a_2,
                score
              );
            cout << "Score for sequence a: " << score << endl;
            //0.0082863
            // USE_DEL_IN_DEL_OUT:
            //0.06812383
            // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
            //0.0840835
            score =
              dp.forward_score( 
                dp_parameters,
                dna_profile,
                dna_seq_a
              );
            cout << "Again, the score for sequence a: " << score << endl;
            
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::Row forward_row_b_1( 3 );
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::Row forward_row_b_2( 3 );
            
            // For these we assume we're *not* using Rabiner Scaling.
            dp_parameters.useRabinerScaling = false;
            
            dp.forward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_b,
              0,
              forward_row_b_2,
              forward_row_b_1
            );
            cout << forward_row_b_1 << endl;
            //[ 1, 0, 0 ], [ 0.025, 0, 0 ], [ 0.000625, 0, 0 ]
            // USE_DEL_IN_DEL_OUT:
            //[ 1, 0, 0, 0, 0 ], [ 0.025, 0, 0, 0, 0 ], [ 0.000625, 0, 0, 0, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1
            // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
            //[ 1, 0, 0 ], [ 0, 0, 0 ], [ 0, 0, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1, deletionIn: 0, deletionOut: 0
            dp.forward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_b,
              1,
              forward_row_b_1,
              forward_row_b_2
            );
            cout << forward_row_b_2 << endl;
            //[ 0, 0, 0.009 ], [ 0.6237, 0, 0.000225 ], [ 0.0022275, 0.0031185, 5.625e-06 ]
            // USE_DEL_IN_DEL_OUT:
            //[ 0, 0, 0.009, 0.45, 0 ], [ 0.3087, 0, 0.000225, 0.01125, 0 ], [ 0.0011025, 0.002701125, 5.625e-06, 0.00028125, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1
            // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
            //[ 0, 0, 0.01 ], [ 0.343, 0, 0 ], [ 0, 0.00300125, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1, deletionIn: 0.5, deletionOut: 0
            dp.forward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_b,
              2,
              forward_row_b_2,
              forward_row_b_1
            );
            cout << forward_row_b_1 << endl;
            //[ 0.0018, 0, 0 ], [ 0.009207, 0, 0 ], [ 0.181804, 0, 0 ]
            // USE_DEL_IN_DEL_OUT:
            //[ 0, 0, 0.0018, 0.225, 0 ], [ 0.09288, 0, 0.00544725, 0.005625, 0.0385875 ], [ 0.07791322, 0.0006966, 2.041875e-05, 0.000140625, 0.0001378125 ], inverseScalar: 1, cumulativeInverseScalar: 1
            // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
            //[ 0, 0, 0.002 ], [ 0.1032, 0, 0.0060025 ], [ 0.08463525, 0.000774, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1, deletionIn: 0.25, deletionOut: 0
#ifdef USE_DEL_IN_DEL_OUT
            dp.forward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_b,
              3,
              forward_row_b_1,
              forward_row_b_2
            );
            cout << forward_row_b_2 << endl;
            // USE_DEL_IN_DEL_OUT:
            //[ 0.00036, 0, 0, 0, 0 ], [ 0.01388565, 0, 0, 0, 0.0618075 ], [ 0.02486032, 0, 0, 0, 0.01961612 ], inverseScalar: 1, cumulativeInverseScalar: 1
            // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
            //[ 0.0004, 0, 0 ], [ 0.0154085, 0, 0 ], [ 0.02453693, 0, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1, deletionIn: 0.125, deletionOut: 0.02115881
#endif // USE_DEL_IN_DEL_OUT .. else ..

            forward_swap =
              dp.forward_score( 
                dp_parameters,
                dna_profile,
                dna_seq_b,
                forward_row_b_1,
                forward_row_b_2,
                score
              );
            cout << "Score for sequence b: " << score << endl;
            //0.163624
            // USE_DEL_IN_DEL_OUT:
            //0.0400288
            // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
            //0.04569574
            score =
              dp.forward_score( 
                dp_parameters,
                dna_profile,
                dna_seq_b
              );
            cout << "Again, the score for sequence b: " << score << endl;
            
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::Row forward_row_c_1( 4 );
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::Row forward_row_c_2( 4 );
            
            // Test using rabiner scaling again:
            dp_parameters.useRabinerScaling = false;//true;
            // Keep the inverse scalars for future reference
#ifdef USE_DEL_IN_DEL_OUT
            vector<MatrixValueType> rabiner_inverse_scalars_c( 4 );
#else
            vector<MatrixValueType> rabiner_inverse_scalars_c( 3 );
#endif // USE_DEL_IN_DEL_OUT .. else ..
            
            dp.forward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_c,
              0,
              forward_row_c_2,
              forward_row_c_1
            );
            if( dp_parameters.useRabinerScaling ) {
              // Save the inverse scalar for use with the backward stuff..
              rabiner_inverse_scalars_c[ 0 ] = forward_row_c_1.m_rabinerInverseScalar;
            
              // Start over with the cumulative inverse scalar
              rabiner_cumulative_inverse_scalar = forward_row_c_1.m_rabinerInverseScalar;
              // With rabiner scaling:
              cout << "Rabiner Scaled:" << endl;
              cout << forward_row_c_1 << endl;
              //[ 0.975, 0, 0 ], [ 0.024375, 0, 0 ], [ 0.000609375, 0, 0 ], [ 1.52344e-05, 0, 0 ]

              // USE_DEL_IN_DEL_OUT:
              //[ 0.493671, 0, 0, 0, 0 ], [ 0.01234178, 0, 0, 0, 0 ], [ 0.0003085444, 0, 0, 0, 0 ], [ 7.713609e-06, 0, 0, 0, 0 ], inverseScalar: 2.025641, cumulativeInverseScalar: 1
              cout << "rabiner sum is " << forward_row_c_1.m_rabinerInverseScalar << endl;
              //1.02564
              cout << "rabiner cumulative inverse scalar is " << rabiner_cumulative_inverse_scalar << endl;
              //1.02564
              // Without rabiner scaling:
              forward_row_c_2 = forward_row_c_1;
              forward_row_c_2 *= rabiner_cumulative_inverse_scalar;
              cout << "Unscaled:" << endl;
              cout << forward_row_c_2 << endl;
              //[ 1, 0, 0 ], [ 0.025, 0, 0 ], [ 0.000625, 0, 0 ], [ 1.5625e-05, 0, 0 ]
              // USE_DEL_IN_DEL_OUT:
              //[ 1, 0, 0, 0, 0 ], [ 0.025, 0, 0, 0, 0 ], [ 0.000625, 0, 0, 0, 0 ], [ 1.5625e-05, 0, 0, 0, 0 ], inverseScalar: 2.025641, cumulativeInverseScalar: 1
            } else {
              cout << forward_row_c_1 << endl;
              // Without rabiner scaling:
              //[ 1, 0, 0 ], [ 0.025, 0, 0 ], [ 0.000625, 0, 0 ], [ 1.5625e-05, 0, 0 ]
            }
            
            dp.forward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_c,
              1,
              forward_row_c_1,
              forward_row_c_2
            );
            if( dp_parameters.useRabinerScaling ) {
              // Save the inverse scalar for use with the backward stuff..
              rabiner_inverse_scalars_c[ 1 ] = forward_row_c_2.m_rabinerInverseScalar;
            
              rabiner_cumulative_inverse_scalar *= forward_row_c_2.m_rabinerInverseScalar;
              // With rabiner scaling:
              cout << "Rabiner Scaled:" << endl;
              cout << forward_row_c_2 << endl;
              //[ 0, 0, 0.0140899 ], [ 0.97643, 0, 0.000352247 ], [ 0.00348725, 0.00488215, 8.80618e-06 ], [ 0.000610269, 0.00013949, 2.20155e-07 ]
              // USE_DEL_IN_DEL_OUT:

              cout << "rabiner sum is " << forward_row_c_2.m_rabinerInverseScalar << endl;
              //0.622787
              cout << "rabiner cumulative inverse scalar is " << rabiner_cumulative_inverse_scalar << endl;
              //0.638756
              // Without rabiner scaling:
              forward_row_c_1 = forward_row_c_2;
              forward_row_c_1 *= rabiner_cumulative_inverse_scalar;
              cout << "Unscaled:" << endl;
              cout << forward_row_c_1 << endl;
              // Without rabiner scaling:
              //[ 0, 0, 0.009 ], [ 0.6237, 0, 0.000225 ], [ 0.0022275, 0.0031185, 5.625e-06 ], [ 0.000389813, 8.91e-05, 1.40625e-07 ]
              // USE_DEL_IN_DEL_OUT:
              //[ 0, 0, 0.009, 0.45, 0 ], [ 0.3087, 0, 0.000225, 0.01125, 0 ], [ 0.0011025, 0.002701125, 5.625e-06, 0.00028125, 0 ], [ 0.0001929375, 7.7175e-05, 1.40625e-07, 7.03125e-06, 0 ], inverseScalar: 1.386812, cumulativeInverseScalar: 1
            } else {
              cout << forward_row_c_2 << endl;
              // Without rabiner scaling:
              //[ 0, 0, 0.009 ], [ 0.6237, 0, 0.000225 ], [ 0.0022275, 0.0031185, 5.625e-06 ], [ 0.000389813, 8.91e-05, 1.40625e-07 ]
            }
            
            dp.forward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_c,
              2,
              forward_row_c_2,
              forward_row_c_1
            );
            if( dp_parameters.useRabinerScaling ) {
              // Save the inverse scalar for use with the backward stuff..
              rabiner_inverse_scalars_c[ 2 ] = forward_row_c_1.m_rabinerInverseScalar;
            
              // Rabiner_sum is 1.0 since this is the last row, so this is redundant:
              rabiner_cumulative_inverse_scalar *= forward_row_c_1.m_rabinerInverseScalar;
              // With rabiner scaling:
              cout << "Rabiner Scaled:" << endl;
              cout << forward_row_c_1 << endl;
              //[ 0.305915, 0, 0 ], [ 1.56476, 0, 0 ], [ 30.8981, 0, 0 ], [ 1.11111, 0, 0 ]
              // USE_DEL_IN_DEL_OUT:

              cout << "rabiner sum is " << forward_row_c_1.m_rabinerInverseScalar << endl;
              //1
              cout << "rabiner cumulative inverse scalar is " << rabiner_cumulative_inverse_scalar << endl;
              //0.638756
              // Without rabiner scaling:
              forward_row_c_2 = forward_row_c_1;
              forward_row_c_2 *= rabiner_cumulative_inverse_scalar;
              cout << "Unscaled:" << endl;
              cout << forward_row_c_2 << endl;
              // Without rabiner scaling:
              //[ 0.0018, 0, 0 ], [ 0.009207, 0, 0 ], [ 0.181804, 0, 0 ], [ 0.00653776, 0, 0 ]
              // USE_DEL_IN_DEL_OUT:
              //[ 0, 0, 0.0018, 0.225, 0 ], [ 0.09288, 0, 0.00544725, 0.005625, 0.0385875 ], [ 0.07791322, 0.0006966, 2.041875e-05, 0.000140625, 0.0001378125 ], [ 0.001393178, 0.0006017642, 3.404531e-06, 3.515625e-06, 2.411719e-05 ], inverseScalar: 1.160287, cumulativeInverseScalar: 1
            } else {
              cout << forward_row_c_1 << endl;
              // Without rabiner scaling:
              //[ 0.0018, 0, 0 ], [ 0.009207, 0, 0 ], [ 0.181804, 0, 0 ], [ 0.00653776, 0, 0 ]
            }
            
#ifdef USE_DEL_IN_DEL_OUT
            dp.forward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_c,
              3,
              forward_row_c_1,
              forward_row_c_2
            );
            if( dp_parameters.useRabinerScaling ) {
              // Save the inverse scalar for use with the backward stuff..
              rabiner_inverse_scalars_c[ 3 ] = forward_row_c_2.m_rabinerInverseScalar;
            
              rabiner_cumulative_inverse_scalar *= forward_row_c_2.m_rabinerInverseScalar;
              // With rabiner scaling:
              cout << "Rabiner Scaled:" << endl;
              cout << forward_row_c_2 << endl;
              // USE_DEL_IN_DEL_OUT:
              // [ 0.05657595, 0, 0, 0, 0 ], [ 2.182205, 0, 0, 0, 9.713383 ], [ 3.906934, 0, 0, 0, 3.082779 ], [ 1.052585, 0, 0, 0, 0.05852649 ], inverseScalar: 0.001952204, cumulativeInverseScalar: 1
              cout << "rabiner sum is " << forward_row_c_2.m_rabinerInverseScalar << endl;

              cout << "rabiner cumulative inverse scalar is " << rabiner_cumulative_inverse_scalar << endl;

              // Without rabiner scaling:
              forward_row_c_1 = forward_row_c_2;
              forward_row_c_1 *= rabiner_cumulative_inverse_scalar;
              cout << "Unscaled:" << endl;
              cout << forward_row_c_1 << endl;
              // USE_DEL_IN_DEL_OUT:
              //[ 0.00036, 0, 0, 0, 0 ], [ 0.01388565, 0, 0, 0, 0.0618075 ], [ 0.02486032, 0, 0, 0, 0.01961612 ], [ 0.006697731, 0, 0, 0, 0.0003724116 ], inverseScalar: 0.001952204, cumulativeInverseScalar: 1
            } else {
              cout << forward_row_c_2 << endl;
              // Without rabiner scaling:

            }
#endif // USE_DEL_IN_DEL_OUT

            forward_swap =
              dp.forward_score(
                dp_parameters,
                dna_profile,
                dna_seq_c,
                forward_row_c_1,
                forward_row_c_2,
                score
              );
            cout << "Score for sequence c: " << score << endl;
            //0.00588399
            // USE_DEL_IN_DEL_OUT:
            //0.006363128
            // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
            //0.006322764
            score =
              dp.forward_score( 
                dp_parameters,
                dna_profile,
                dna_seq_c
              );
            cout << "Again, the score for sequence c: " << score << endl;
            
            // (1-.01)*.7 * ( 1.0 - .1 ) :0.6237000
            // .9*.01*.025:0.0002250
            // .009*.2:0.0018000
            // .0018*.1*.25:0.0000450
            // .6237*.01:0.0062370
            // .000225*.2:0.0000450
            // .009*.8*.4:0.0028800
            // .0018*.1*.25+.6237*.01+.000225*.2+.009*.8*.4:0.0092070
            // 
            // .9*.99*.1*.025:0.0022275
            // .6237*.02*.25:0.0031185
            // .000625*.9*.01:0.0000056
            // 
            // .01*.0022275:0.0000223
            // .000625*.9*.01*.2:0.0000011
            // .1*.25*.009207:0.0002302
            // .6237*.97*.3:0.1814967
            // .000225*.8*.3:0.0000540
            // .01*.0022275+.000625*.9*.01*.2+.1*.25*.009207+.6237*.97*.3+.000225*.8*.3:0.1818043
            
            // Now calculate them all together
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::RowVector forward_rows_1( dna_seqs, 3 );
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::RowVector forward_rows_2( dna_seqs, 3 );
            forward_rows_1[ 0 ] = forward_row_a_1;
            forward_rows_2[ 0 ] = forward_row_a_2;
            forward_rows_1[ 1 ] = forward_row_b_1;
            forward_rows_2[ 1 ] = forward_row_b_2;
            forward_rows_1[ 2 ] = forward_row_c_1;
            forward_rows_2[ 2 ] = forward_row_c_2;
            
            // Expecting:
            // 0.005883984*0.163624*0.0082863 = 7.97772646e-06
            bool forward_swap_all =
              dp.forward_score(
                dp_parameters,
                false, // Don't use viterbi
                dna_profile,
                dna_seqs,
                num_sequences_to_use,
                forward_rows_1,
                forward_rows_2,
                score
              );
            cout << "Total (product) score for all sequences: " << score << endl;
            // Total (product) score for all sequences: 7.97772e-06
            // USE_DEL_IN_DEL_OUT:
            // Total (product) score for all sequences: 1.735171e-05
            // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
            // Total (product) score for all sequences: 2.429369e-05
            score =
              dp.forward_score(
                dp_parameters,
                dna_profile,
                dna_seqs,
                num_sequences_to_use
              );
            cout << "Again, the total score for all sequences: " << score << endl;
            //Again, the total score for all sequences: 7.97772e-06
            // USE_DEL_IN_DEL_OUT:
            //Again, the total score for all sequences: 1.735171e-05


            // Now test another forward_score method, and the forward_reverseCalculateRow(..) method.
            //dp_parameters.useRabinerScaling = false; // to make evaluation easier.
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer anchor_columns(
              dna_profile,
              dna_seqs,
              num_sequences_to_use
            );
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer anchor_rows(
              dna_profile,
              dna_seqs,
              num_sequences_to_use
            );
            // TODO: MAGIC # (store_every_Nth_column)
            uint32_t store_every_nth_column = 5;
            if( store_every_nth_column > 0 ) {
              anchor_columns.reinitialize(
                dna_profile.length(),
                dna_seqs,
                num_sequences_to_use,
                store_every_nth_column,
                ( ( store_every_nth_column == 0 ) ? 0 : numeric_limits<uint32_t>::max() )
              );
            }
            // TODO: MAGIC # (store_every_Nth_row)
            // TODO: Make this a parameter
            anchor_rows.reinitialize(
              dna_profile.length(),
              dna_seqs,
              num_sequences_to_use,
              numeric_limits<uint32_t>::max(),
              3//1
            );
            vector<ScoreType> sequence_scores( num_sequences_to_use );
            ScoreType largest_sequence_score;
            if( dp_parameters.useRabinerScaling == false ) {
              for( uint32_t seq_i = 0; seq_i < num_sequences_to_use; seq_i++ ) {
                forward_rows_1[ seq_i ].m_rabinerInverseScalar = 1;
                forward_rows_1[ seq_i ].m_rabinerCumulativeInverseScalar = 1;
                forward_rows_2[ seq_i ].m_rabinerInverseScalar = 1;
                forward_rows_2[ seq_i ].m_rabinerCumulativeInverseScalar = 1;
              }
            }
            forward_swap_all =
              dp.forward_score(
                dp_parameters,
                false, // don't use viterbi
                dna_profile,
                dna_seqs,
                num_sequences_to_use,
                &anchor_columns,
                &anchor_rows,
                forward_rows_1,
                forward_rows_2,
                &sequence_scores,
                &largest_sequence_score,
                score
              );
            cout << "Finally, the total score for all sequences: " << score << endl;
            // USE_DEL_IN_DEL_OUT:
            // Finally, the total score for all sequences: 1.735171e-05
            cout << "The sequence scores are:\n";
            for( uint32_t seq_i = 0; seq_i < num_sequences_to_use; seq_i++ ) {
              cout << "sequence " << seq_i << ": " << sequence_scores[ seq_i ] << endl;
            }
            cout << "The largest sequence score is: " << largest_sequence_score << endl;

            // Now use those anchor rows & cols to do reverseCalculateRow.
            // Unfortunately we can only do this if the profile has length > 2.
            // When USE_DEL_IN_DEL_OUT, our profile has length 3..
#ifdef USE_DEL_IN_DEL_OUT
            // Note that we haven't implemented forward_reverseCalculateRow for
            // the case in which USE_DEL_IN_DEL_OUT but
            // !DISALLOW_FLANKING_TRANSITIONS
#ifdef DISALLOW_FLANKING_TRANSITIONS
            // First do it without any anchor cols.
            cout << "reverse-calculating forward row 1 for seq c with no anchor cols." << endl;
            cout << "\tthe subsequent row is:" << endl <<
              ( forward_swap_all ? forward_rows_2[ 2 ] : forward_rows_1[ 2 ] ) << endl;
            if( anchor_columns.m_storeEveryNthColumn > 0 ) {
              forward_matrices_iterator = anchor_columns.begin();
              forward_matrices_iterator++; // for row 1

              // Restore just the final column
              ( forward_swap_all ? forward_rows_1[ 2 ] : forward_rows_2[ 2 ] )[ dna_seqs[ 2 ].length() ] =
                ( *forward_matrices_iterator )[ 2 ][ ( uint32_t )( dna_seqs[ 2 ].length() / anchor_columns.m_storeEveryNthColumn ) + 1 ];
              // And the delIn and delOut values..
              ( forward_swap_all ? forward_rows_1[ 2 ] : forward_rows_2[ 2 ] ).m_deletionIn =
                ( *forward_matrices_iterator )[ 2 ].m_deletionIn;
              ( forward_swap_all ? forward_rows_1[ 2 ] : forward_rows_2[ 2 ] ).m_deletionOut =
                ( *forward_matrices_iterator )[ 2 ].m_deletionOut;

              if( dp_parameters.useRabinerScaling ) {
                // Copy the rabiner scalars.
                ( forward_swap_all ? forward_rows_1[ 2 ] : forward_rows_2[ 2 ] ).m_rabinerInverseScalar =
                  ( *forward_matrices_iterator )[ 2 ].m_rabinerInverseScalar;
                ( forward_swap_all ? forward_rows_1[ 2 ] : forward_rows_2[ 2 ] ).m_rabinerCumulativeInverseScalar =
                  ( *forward_matrices_iterator )[ 2 ].m_rabinerCumulativeInverseScalar;
              } // End if useRabinerScaling
            } // End if m_anchor_columns.m_storeEveryNthColumn > 0
            dp.forward_reverseCalculateRow(
              dp_parameters,
              dna_profile,
              dna_seqs[ 2 ],
              1,
              0,
              ( forward_swap_all ? forward_rows_2[ 2 ] : forward_rows_1[ 2 ] ),
              ( forward_swap_all ? forward_rows_1[ 2 ] : forward_rows_2[ 2 ] )
            );
            forward_row_c_2 =
              ( forward_swap_all ? forward_rows_1[ 2 ] : forward_rows_2[ 2 ] );
            forward_row_c_2.m_rabinerInverseScalar =
              ( forward_swap_all ? forward_rows_1[ 2 ] : forward_rows_2[ 2 ] ).m_rabinerInverseScalar;
            forward_row_c_2.m_rabinerCumulativeInverseScalar =
              ( forward_swap_all ? forward_rows_1[ 2 ] : forward_rows_2[ 2 ] ).m_rabinerCumulativeInverseScalar;
            if( dp_parameters.useRabinerScaling ) {
              // With rabiner scaling:
              cout << "Rabiner Scaled:" << endl;
              cout << forward_row_c_2 << endl;
              //[ 0, 0, 0.0140899 ], [ 0.97643, 0, 0.000352247 ], [ 0.00348725, 0.00488215, 8.80618e-06 ], [ 0.000610269, 0.00013949, 2.20155e-07 ]
              // USE_DEL_IN_DEL_OUT:

              cout << "rabiner sum is " << forward_row_c_2.m_rabinerInverseScalar << endl;
              //0.622787
              cout << "rabiner cumulative inverse scalar is " << rabiner_cumulative_inverse_scalar << endl;
              //0.638756
              // Without rabiner scaling:
              forward_row_c_1 = forward_row_c_2;
              forward_row_c_1 *= forward_row_c_2.m_rabinerCumulativeInverseScalar;
              cout << "Unscaled:" << endl;
              cout << forward_row_c_1 << endl;
              // Without rabiner scaling:
              //[ 0, 0, 0.009 ], [ 0.6237, 0, 0.000225 ], [ 0.0022275, 0.0031185, 5.625e-06 ], [ 0.000389813, 8.91e-05, 1.40625e-07 ]
              // USE_DEL_IN_DEL_OUT:
              //[ 0, 0, 0.009, 0.45, 0 ], [ 0.3087, 0, 0.000225, 0.01125, 0 ], [ 0.0011025, 0.002701125, 5.625e-06, 0.00028125, 0 ], [ 0.0001929375, 7.7175e-05, 1.40625e-07, 7.03125e-06, 0 ], inverseScalar: 1.386812, cumulativeInverseScalar: 1
            } else {
              cout << forward_row_c_2 << endl;
              // Without rabiner scaling:
              //[ 0, 0, 0.009 ], [ 0.6237, 0, 0.000225 ], [ 0.0022275, 0.0031185, 5.625e-06 ], [ 0.000389813, 8.91e-05, 1.40625e-07 ]
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 0, 0, 0.01 ], [ 0.343, 0, 0 ], [ 0, 0.00300125, 0 ], [ 0, 7.503125e-05, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1, deletionIn: 0.5, deletionOut: 0
            }

            cout << "reverse-calculating forward row 1 for seq c with " << anchor_columns.m_storeEveryNthColumn << " anchor cols." << endl;
            cout << "\tthe subsequent row is:" << endl <<
              ( forward_swap_all ? forward_rows_2[ 2 ] : forward_rows_1[ 2 ] ) << endl;
            if( anchor_columns.m_storeEveryNthColumn > 0 ) {
              forward_matrices_iterator = anchor_columns.begin();
              forward_matrices_iterator++; // for row 1
              ( *forward_matrices_iterator )[ 2 ].restoreColumns(
                anchor_columns.m_storeEveryNthColumn,
                ( forward_swap_all ? forward_rows_1[ 2 ] : forward_rows_2[ 2 ] )
              );
                          
              if( dp_parameters.useRabinerScaling ) {
                // Copy the rabiner scalars.
                ( forward_swap_all ? forward_rows_1[ 2 ] : forward_rows_2[ 2 ] ).m_rabinerInverseScalar =
                  ( *forward_matrices_iterator )[ 2 ].m_rabinerInverseScalar;
                ( forward_swap_all ? forward_rows_1[ 2 ] : forward_rows_2[ 2 ] ).m_rabinerCumulativeInverseScalar =
                  ( *forward_matrices_iterator )[ 2 ].m_rabinerCumulativeInverseScalar;
              } // End if useRabinerScaling
            } // End if m_anchor_columns.m_storeEveryNthColumn > 0
            dp.forward_reverseCalculateRow(
              dp_parameters,
              dna_profile,
              dna_seqs[ 2 ],
              1,
              anchor_columns.m_storeEveryNthColumn,
              ( forward_swap_all ? forward_rows_2[ 2 ] : forward_rows_1[ 2 ] ),
              ( forward_swap_all ? forward_rows_1[ 2 ] : forward_rows_2[ 2 ] )
            );
            forward_row_c_2 =
              ( forward_swap_all ? forward_rows_1[ 2 ] : forward_rows_2[ 2 ] );
            forward_row_c_2.m_rabinerInverseScalar =
              ( forward_swap_all ? forward_rows_1[ 2 ] : forward_rows_2[ 2 ] ).m_rabinerInverseScalar;
            forward_row_c_2.m_rabinerCumulativeInverseScalar =
              ( forward_swap_all ? forward_rows_1[ 2 ] : forward_rows_2[ 2 ] ).m_rabinerCumulativeInverseScalar;
            if( dp_parameters.useRabinerScaling ) {
              // With rabiner scaling:
              cout << "Rabiner Scaled:" << endl;
              cout << forward_row_c_2 << endl;
              //[ 0, 0, 0.0140899 ], [ 0.97643, 0, 0.000352247 ], [ 0.00348725, 0.00488215, 8.80618e-06 ], [ 0.000610269, 0.00013949, 2.20155e-07 ]
              // USE_DEL_IN_DEL_OUT:

              cout << "rabiner sum is " << forward_row_c_2.m_rabinerInverseScalar << endl;
              //0.622787
              cout << "rabiner cumulative inverse scalar is " << rabiner_cumulative_inverse_scalar << endl;
              //0.638756
              // Without rabiner scaling:
              forward_row_c_1 = forward_row_c_2;
              forward_row_c_1 *= forward_row_c_2.m_rabinerCumulativeInverseScalar;
              cout << "Unscaled:" << endl;
              cout << forward_row_c_1 << endl;
              // Without rabiner scaling:
              //[ 0, 0, 0.009 ], [ 0.6237, 0, 0.000225 ], [ 0.0022275, 0.0031185, 5.625e-06 ], [ 0.000389813, 8.91e-05, 1.40625e-07 ]
              // USE_DEL_IN_DEL_OUT:
              //[ 0, 0, 0.009, 0.45, 0 ], [ 0.3087, 0, 0.000225, 0.01125, 0 ], [ 0.0011025, 0.002701125, 5.625e-06, 0.00028125, 0 ], [ 0.0001929375, 7.7175e-05, 1.40625e-07, 7.03125e-06, 0 ], inverseScalar: 1.386812, cumulativeInverseScalar: 1
            } else {
              cout << forward_row_c_2 << endl;
              // Without rabiner scaling:
              //[ 0, 0, 0.009 ], [ 0.6237, 0, 0.000225 ], [ 0.0022275, 0.0031185, 5.625e-06 ], [ 0.000389813, 8.91e-05, 1.40625e-07 ]
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 0, 0, 0.01 ], [ 0.343, 0, 0 ], [ 0, 0.00300125, 0 ], [ 0, 7.503125e-05, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1, deletionIn: 0.5, deletionOut: 0
            }
#endif // DISALLOW_FLANKING_TRANSITIONS
#endif // USE_DEL_IN_DEL_OUT

            //// backwards ////

            // Do this part without scaling..
            dp_parameters.useRabinerScaling = false;

            // If we are using Rabiner scaling, we'll need to keep track of the
            // cumulative backward scale factor.
            ScoreType rabiner_backward_cumulative_inverse_scalar = 1.0;

            // a
            rabiner_backward_cumulative_inverse_scalar = 1.0;
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::Row backward_row_a_1( 2 );
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::Row backward_row_a_2( 2 );

#ifdef USE_DEL_IN_DEL_OUT
            if( dp_parameters.useRabinerScaling ) {
              backward_row_a_2.m_rabinerInverseScalar =
                rabiner_inverse_scalars_a[ 3 ];
            } else {
              assert( backward_row_a_2.m_rabinerInverseScalar == 1 );
            }
            assert( backward_row_a_2.m_rabinerInverseScalar != 0 );
            dp.backward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_a,
              3,
              backward_row_a_1,
              backward_row_a_2
            );
            if( dp_parameters.useRabinerScaling ) {
              rabiner_backward_cumulative_inverse_scalar *=
                backward_row_a_2.m_rabinerInverseScalar;
              // With rabiner scaling:
              cout << "Rabiner Scaled:" << endl;
              cout << backward_row_a_2 << endl;
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 0, 0, 0 ], [ 1, 0, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1, deletionIn: 0, deletionOut: .5
              cout << "rabiner sum is " << backward_row_a_2.m_rabinerInverseScalar << endl;
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              // 1
              cout << "rabiner backward cumulative inverse scalar is " << rabiner_backward_cumulative_inverse_scalar << endl;
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              // 1
              // Without rabiner scaling:
              backward_row_a_1 = backward_row_a_2;
              backward_row_a_1 *= rabiner_backward_cumulative_inverse_scalar;
              cout << "Unscaled:" << endl;
              cout << backward_row_a_1 << endl;
              // Without rabiner scaling:
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 0, 0, 0 ], [ 1, 0, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1, deletionIn: 0, deletionOut: .5
            } else {
              cout << backward_row_a_2 << endl;
              // Without rabiner scaling:

            }
#endif // USE_DEL_IN_DEL_OUT

            if( dp_parameters.useRabinerScaling ) {
              backward_row_a_1.m_rabinerInverseScalar =
                rabiner_inverse_scalars_a[ 2 ];
            } else {
              assert( backward_row_a_1.m_rabinerInverseScalar == 1 );
            }
            assert( backward_row_a_1.m_rabinerInverseScalar != 0 );
            dp.backward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_a,
              2,
              backward_row_a_2,
              backward_row_a_1
            );
            if( dp_parameters.useRabinerScaling ) {
              rabiner_backward_cumulative_inverse_scalar *=
                backward_row_a_1.m_rabinerInverseScalar;
              // With rabiner scaling:
              cout << "Rabiner Scaled:" << endl;
              cout << backward_row_a_1 << endl;
              //[ 1, 0, 0 ], [ 40, 0, 0 ], inverseScalar: 0.0225, cumulativeInverseScalar: 1
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 0.09992913, 0.1275691, 0.1133948 ], [ 0.3756201, 0, 0.2834869 ], inverseScalar: 0.7055, cumulativeInverseScalar: 1, deletionIn: 0.07087173, deletionOut: 0.7087172

              // Note that this depends on the flag in backward_calculateRow(..) called calculate_rabiner_inverse_scalar_here:
              cout << "rabiner sum is " << backward_row_a_1.m_rabinerInverseScalar << endl;
              //0.0225
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //0.7055 (with calculate_rabiner_inverse_scalar_here=true)
              cout << "rabiner backward cumulative inverse scalar is " << rabiner_backward_cumulative_inverse_scalar << endl;
              //0.0225
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //0.7055
              // Without rabiner scaling:
              backward_row_a_2 = backward_row_a_1;
              backward_row_a_2 *= rabiner_backward_cumulative_inverse_scalar;
              cout << "Unscaled:" << endl;
              cout << backward_row_a_2 << endl;
              //[ 0.0225, 0, 0 ], [ 0.9, 0, 0 ], inverseScalar: 0.0225, cumulativeInverseScalar: 1
              // Without rabiner scaling:
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 0.0705, 0, 0.08 ], [ 0.265, 0, 0.2 ], inverseScalar: 0.7055, cumulativeInverseScalar: 1, deletionIn: 0.05, deletionOut: 0.25
            } else {
              cout << backward_row_a_1 << endl;
              // Without rabiner scaling:

            }
            
            if( dp_parameters.useRabinerScaling ) {
              backward_row_a_2.m_rabinerInverseScalar =
                rabiner_inverse_scalars_a[ 1 ];
            } else {
              assert( backward_row_a_2.m_rabinerInverseScalar == 1 );
            }
            dp.backward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_a,
              1,
              backward_row_a_1,
              backward_row_a_2
            );
            if( dp_parameters.useRabinerScaling ) {
              rabiner_backward_cumulative_inverse_scalar *=
                backward_row_a_2.m_rabinerInverseScalar;
              // With rabiner scaling:
              cout << "Rabiner Scaled:" << endl;
              cout << backward_row_a_2 << endl;
              //[ 0.3025521, 0.2805377, 0.2532632 ], [ 0.007792714, 0, 0.1558543 ], inverseScalar: 51.33, cumulativeInverseScalar: 1
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 0.2266847, 0.244124, 0.2579423 ], [ 0.1688908, 0, 0.1023581 ], inverseScalar: 0.5539121, cumulativeInverseScalar: 1, deletionIn: 0.1995982, deletionOut: 0.639738
              cout << "rabiner sum is " << backward_row_a_2.m_rabinerInverseScalar << endl;
              //51.33
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //0.5539121
              cout << "rabiner backward cumulative inverse scalar is " << rabiner_backward_cumulative_inverse_scalar << endl;
              //1.154925
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //0.390785
              // Without rabiner scaling:
              backward_row_a_1 = backward_row_a_2;
              backward_row_a_1 *= rabiner_backward_cumulative_inverse_scalar;
              cout << "Unscaled:" << endl;
              cout << backward_row_a_1 << endl;
              // Without rabiner scaling:
              //[ 0.349425, 0.324, 0.2925 ], [ 0.009, 0, 0.18 ], inverseScalar: 51.33, cumulativeInverseScalar: 1
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 0.088585, 0.0954, 0.1008 ], [ 0.066, 0, 0.04 ], inverseScalar: 0.5539121, cumulativeInverseScalar: 1, deletionIn: 0.078, deletionOut: 0.25
            } else {
              cout << backward_row_a_2 << endl;
              // Without rabiner scaling:

            }
            
            if( dp_parameters.useRabinerScaling ) {
              backward_row_a_1.m_rabinerInverseScalar =
                rabiner_inverse_scalars_a[ 0 ];
            }
            dp.backward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_a,
              0,
              backward_row_a_2,
              backward_row_a_1
            );
            if( dp_parameters.useRabinerScaling ) {
              rabiner_backward_cumulative_inverse_scalar *=
                backward_row_a_1.m_rabinerInverseScalar;
              // With rabiner scaling:
              cout << "Rabiner Scaled:" << endl;
              cout << backward_row_a_1 << endl;
              //[ 0.8364677, 0, 0 ], [ 0.1635323, 0, 0 ], inverseScalar: 0.00857744, cumulativeInverseScalar: 1
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 0.9936554, 0, 0 ], [ 0.006344574, 0, 0 ], inverseScalar: 0.1613317, cumulativeInverseScalar: 1, deletionIn: 0, deletionOut: 0
              cout << "rabiner sum is " << backward_row_a_1.m_rabinerInverseScalar << endl;
              //0.00857744
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //0.1613317
              cout << "rabiner backward cumulative inverse scalar is " << rabiner_backward_cumulative_inverse_scalar << endl;
              //0.0099063
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //0.063046
              // Without rabiner scaling:
              backward_row_a_2 = backward_row_a_1;
              backward_row_a_2 *= rabiner_backward_cumulative_inverse_scalar;
              cout << "Unscaled:" << endl;
              cout << backward_row_a_2 << endl;
              // Without rabiner scaling:
              //[ 0.0082863, 0, 0 ], [ 0.00162, 0, 0 ], inverseScalar: 0.00857744, cumulativeInverseScalar: 1
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 0.062646, 0, 0 ], [ 0.0004, 0, 0 ], inverseScalar: 0.1613317, cumulativeInverseScalar: 1, deletionIn: 0, deletionOut: 0
            } else {
              cout << backward_row_a_1 << endl;
              // Without rabiner scaling:

            }
            
            bool backward_swap =
              dp.backward_score( 
                dp_parameters,
                dna_profile,
                dna_seq_a,
                ( dp_parameters.useRabinerScaling ? &rabiner_inverse_scalars_c : NULL ),
                backward_row_a_1,
                backward_row_a_2,
                score
              );
            cout << "Backwards, score for sequence a: " << score << endl;
            //0.0082863

            score =
              dp.backward_score( 
                dp_parameters,
                dna_profile,
                dna_seq_a,
                ( dp_parameters.useRabinerScaling ? &rabiner_inverse_scalars_a : NULL )
              );
            cout << "Again backwards, the score for sequence a: " << score << endl;

            // For this, do use rabiner scaling
            dp_parameters.useRabinerScaling = false;//true;
            // c
            rabiner_backward_cumulative_inverse_scalar = 1.0;
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::Row backward_row_c_1( 4 );
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::Row backward_row_c_2( 4 );

#ifdef USE_DEL_IN_DEL_OUT
            if( dp_parameters.useRabinerScaling ) {
              backward_row_c_2.m_rabinerInverseScalar =
                rabiner_inverse_scalars_c[ 3 ];
            } else {
              assert( backward_row_c_2.m_rabinerInverseScalar == 1 );
            }
            assert( backward_row_c_2.m_rabinerInverseScalar != 0 );
            dp.backward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_c,
              3,
              backward_row_c_1,
              backward_row_c_2
            );
            if( dp_parameters.useRabinerScaling ) {
              rabiner_backward_cumulative_inverse_scalar *=
                backward_row_c_2.m_rabinerInverseScalar;
              // With rabiner scaling:
              cout << "Rabiner Scaled:" << endl;
              cout << backward_row_c_2 << endl;
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 0, 0, 0 ], [ 0, 0, 0 ], [ 0, 0, 0 ], [ 1, 0, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1, deletionIn: 0, deletionOut: 1              
              cout << "rabiner sum is " << backward_row_c_2.m_rabinerInverseScalar << endl;
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //1
              cout << "rabiner backward cumulative inverse scalar is " << rabiner_backward_cumulative_inverse_scalar << endl;
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //1
              // Without rabiner scaling:
              backward_row_c_1 = backward_row_c_2;
              backward_row_c_1 *= rabiner_backward_cumulative_inverse_scalar;
              cout << "Unscaled:" << endl;
              cout << backward_row_c_1 << endl;
              // Without rabiner scaling:
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 0, 0, 0 ], [ 0, 0, 0 ], [ 0, 0, 0 ], [ 1, 0, 0 ], inverseScalar: 1, cumulativeInverseScalar: 1, deletionIn: 0, deletionOut: 1
            } else {
              cout << backward_row_c_2 << endl;
              // Without rabiner scaling:

            }
#endif // USE_DEL_IN_DEL_OUT

            if( dp_parameters.useRabinerScaling ) {
              backward_row_c_1.m_rabinerInverseScalar =
                rabiner_inverse_scalars_c[ 2 ];
            } else {
              assert( backward_row_c_1.m_rabinerInverseScalar == 1 );
            }
            assert( backward_row_c_1.m_rabinerInverseScalar != 0 );
            dp.backward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_c,
              2,
              backward_row_c_2,
              backward_row_c_1
            );
            if( dp_parameters.useRabinerScaling ) {
              rabiner_backward_cumulative_inverse_scalar *=
                backward_row_c_1.m_rabinerInverseScalar;
              // With rabiner scaling:
              cout << "Rabiner Scaled:" << endl;
              cout << backward_row_c_1 << endl;
              //[ 0.0015266, 0, 0 ], [ 0.0610641, 0, 0 ], [ 2.44256, 0, 0 ], [ 97.7025, 0, 0 ]
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 2.381799e-05, 7.939329e-05, 0 ], [ 0.0009527195, 0.003175732, 0 ], [ 0.09950626, 0.1270293, 0.1129149 ], [ 0.3740306, 0, 0.2822873 ], inverseScalar: 0.7084981, cumulativeInverseScalar: 1, deletionIn: 0, deletionOut: 0.7057182
              // Note that this depends on the flag in backward_calculateRow(..) called calculate_rabiner_inverse_scalar_here:
              cout << "rabiner sum is " << backward_row_c_1.m_rabinerInverseScalar << endl;
              // 0.002208375 (with calculate_rabiner_inverse_scalar_here=false)
              // 0.02307656  (with calculate_rabiner_inverse_scalar_here=true)
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              // 0.7084981 (with calculate_rabiner_inverse_scalar_here=true)
              cout << "rabiner backward cumulative inverse scalar is " << rabiner_backward_cumulative_inverse_scalar << endl;
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              // 0.7084981
              // Without rabiner scaling:
              backward_row_c_2 = backward_row_c_1;
              backward_row_c_2 *= rabiner_backward_cumulative_inverse_scalar;
              cout << "Unscaled:" << endl;
              cout << backward_row_c_2 << endl;
              // Without rabiner scaling:
              //[ 1.40625e-05, 0, 0 ], [ 0.0005625, 0, 0 ], [ 0.0225, 0, 0 ], [ 0.9, 0, 0 ]
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 1.6875e-05, 5.625e-05, 0 ], [ 0.000675, 0.00225, 0 ], [ 0.0705, 0.09, 0.08 ], [ 0.265, 0, 0.2 ], inverseScalar: 0.7084981, cumulativeInverseScalar: 1, deletionIn: 0, deletionOut: 0.5
            } else {
              cout << backward_row_c_1 << endl;
              // Without rabiner scaling:
              //[ 1.40625e-05, 0, 0 ], [ 0.0005625, 0, 0 ], [ 0.0225, 0, 0 ], [ 0.9, 0, 0 ]
            }
            
            backward_row_c_2.m_rabinerInverseScalar =
              rabiner_inverse_scalars_c[ 1 ];
            dp.backward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_c,
              1,
              backward_row_c_1,
              backward_row_c_2
            );
            if( dp_parameters.useRabinerScaling ) {
              rabiner_backward_cumulative_inverse_scalar *=
                backward_row_c_2.m_rabinerInverseScalar;
              // With rabiner scaling:
              cout << "Rabiner Scaled:" << endl;
              cout << backward_row_c_2 << endl;
              //[ 0.050422, 0.0970692, 0.0318661 ], [ 1.42466, 2.47085, 0.960887 ], [ 60.9084, 56.4766, 50.9858 ], [ 1.56879, 0, 31.3759 ]
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 0.0009125153, 0.001734782, 0.0004813267 ], [ 0.04062448, 0.04773157, 0.03770393 ], [ 0.1973997, 0.212586, 0.2246191 ], [ 0.147072, 0, 0.08913458 ], inverseScalar: 0.6333957, cumulativeInverseScalar: 1, deletionIn: 0.0003008292, deletionOut: 0.5570911
              cout << "rabiner sum is " << backward_row_c_2.m_rabinerInverseScalar << endl;
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //0.6333957
              cout << "rabiner backward cumulative inverse scalar is " << rabiner_backward_cumulative_inverse_scalar << endl;
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //0.4487596
              // Without rabiner scaling:
              backward_row_c_1 = backward_row_c_2;
              backward_row_c_1 *= rabiner_backward_cumulative_inverse_scalar;
              cout << "Unscaled:" << endl;
              cout << backward_row_c_1 << endl;
              // Without rabiner scaling:
              //[ 0.000289266, 0.000556875, 0.000182813 ], [ 0.00817313, 0.014175, 0.0055125 ], [ 0.349425, 0.324, 0.2925 ], [ 0.009, 0, 0.18 ]
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 0.0004095, 0.0007785, 0.000216 ], [ 0.01823062, 0.02142, 0.01692 ], [ 0.088585, 0.0954, 0.1008 ], [ 0.066, 0, 0.04 ], inverseScalar: 0.6333957, cumulativeInverseScalar: 1, deletionIn: 0.000135, deletionOut: 0.25
            } else {
              cout << backward_row_c_2 << endl;
              // Without rabiner scaling:
              //[ 0.000289266, 0.000556875, 0.000182813 ], [ 0.00817313, 0.014175, 0.0055125 ], [ 0.349425, 0.324, 0.2925 ], [ 0.009, 0, 0.18 ]
            }
            
            backward_row_c_1.m_rabinerInverseScalar =
              rabiner_inverse_scalars_c[ 0 ];
            
            dp.backward_calculateRow(
              dp_parameters,
              dna_profile,
              dna_seq_c,
              0,
              backward_row_c_2,
              backward_row_c_1
            );
            if( dp_parameters.useRabinerScaling ) {
              rabiner_backward_cumulative_inverse_scalar *=
                backward_row_c_1.m_rabinerInverseScalar;
              // With rabiner scaling:
              cout << "Rabiner Scaled:" << endl;
              cout << backward_row_c_1 << endl;
              //[ 1, 0, 0 ], [ 5.33491, 0, 0 ], [ 1.40828, 0, 0 ], [ 0.275324, 0, 0 ]
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 0.181279, 0, 0 ], [ 0.1293017, 0, 0 ], [ 0.677951, 0, 0 ], [ 0.01146834, 0, 0 ], inverseScalar: 0.0777223, cumulativeInverseScalar: 1, deletionIn: 0, deletionOut: 0
              cout << "rabiner sum is " << backward_row_c_1.m_rabinerInverseScalar << endl;
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //0.0777223
              cout << "rabiner backward cumulative inverse scalar is " << rabiner_backward_cumulative_inverse_scalar << endl;
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //0.03487863
              // Without rabiner scaling:
              backward_row_c_2 = backward_row_c_1;
              backward_row_c_2 *= rabiner_backward_cumulative_inverse_scalar;
              cout << "Unscaled:" << endl;
              cout << backward_row_c_2 << endl;
              // Without rabiner scaling:
              //[ 0.00588399, 0, 0 ], [ 0.0313905, 0, 0 ], [ 0.0082863, 0, 0 ], [ 0.00162, 0, 0 ]
              // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
              //[ 0.006322764, 0, 0 ], [ 0.004509865, 0, 0 ], [ 0.023646, 0, 0 ], [ 0.0004, 0, 0 ], inverseScalar: 0.0777223, cumulativeInverseScalar: 1, deletionIn: 0, deletionOut: 0
            } else {
              cout << backward_row_c_1 << endl;
              // Without rabiner scaling:
              //[ 0.00588399, 0, 0 ], [ 0.0313905, 0, 0 ], [ 0.0082863, 0, 0 ], [ 0.00162, 0, 0 ]
            }
            
            backward_swap =
              dp.backward_score( 
                dp_parameters,
                dna_profile,
                dna_seq_c,
                ( dp_parameters.useRabinerScaling ? &rabiner_inverse_scalars_c : NULL ),
                backward_row_c_1,
                backward_row_c_2,
                score
              );
            cout << "Backwards, score for sequence c: " << score << endl;
            //0.00588399
            // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
            //0.006322764
            score =
              dp.backward_score( 
                dp_parameters,
                dna_profile,
                dna_seq_c,
                ( dp_parameters.useRabinerScaling ? &rabiner_inverse_scalars_c : NULL )
              );
            cout << "Again backwards, the score for sequence c: " << score << endl;
            //0.00588399
            
            //.025 * .9 *.025 * .025:0.000014063
            //.0225*.2:0.004500000
            //.9*.4*.8:0.288000000
            //.0225*.2+.9*.4*.8:0.292500000
            //.9*.4*.9:0.324000000
            //.97*.4*.9:0.349200000
            //.0225*.01:0.000225000
            //.0225*.8*.3:0.005400000
            //.0005625*.2:0.000112500
            //.0225*.9*.3:0.006075000
            //.324*.1*.25:0.008100000
            //.0225*.97*.3:0.006547500
            //.324*.02*.25:0.001620000
            //.0005625*.01:0.000005625
            //.0225*.97*.3+.324*.02*.25+.0005625*.01:0.008173125
            //.025 * .9 *.025 * .025*.2:0.000002813
            //.0005625*.8*.4:0.000180000
            //.0005625*.9*.4:0.000202500
            //.014175*.1*.25:0.000354375
            //.0005625*.97*.4:0.000218250
            //.014175*.02*.25:0.000070875
            //.025 * .9 *.025 * .025*.01:0.000000141
            //.0005625*.97*.4+.014175*.02*.25+.025 * .9 *.025 * .025*.01:0.000289266
            //.18*.9*.01:0.001620000
            //.00162*.1*.25:0.000040500
            //.009*.9*.99*.7:0.005613300
            //.2925*.9*.01:0.002632500
            //.00162*.1*.25+.009*.9*.99*.7+.2925*.9*.01:0.008286300
            //
            //0.00653776*.9:0.005883984
            
          } // End if test_dp_basics
        
          score =
            dp.forward_score(
              dp_parameters,
              dna_profile,
              dna_seqs,
              num_sequences_to_use,
              forward_matrices
            );
          if( test_dp_basics ) {
            cout << "And yet again (forward), the total score for all sequences: " << score << endl;
          } // End if test_dp_basics
          if( dp_parameters.useRabinerScaling ) {
            backward_matrices.copyRabinerInverseScalarsFrom( forward_matrices );
          }
          score =
            dp.backward_score(
              dp_parameters,
              dna_profile,
              dna_seqs,
              num_sequences_to_use,
              backward_matrices
            );
          assert( score != 0 );
          assert( !isinf( score ) );
          if( test_dp_basics ) {
            cout << "Yes, again (backward), the total score for all sequences: " << score << endl;
            // USE_DEL_IN_DEL_OUT:
            // 2.429369e-05
            // USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS:
            // 1.735171e-05

            // These may have been previously calculated without using rabiner
            // scaling, so we need to recalculate them with the scaling...
            dp_parameters.useRabinerScaling = false;//true;

            score =
              dp.forward_score( 
                dp_parameters,
                dna_profile,
                dna_seqs,
                num_sequences_to_use,
                forward_matrices
              );
            cout << "And yet again (forward, with scaling), the total score for all sequences: " << score << endl;
            backward_matrices.copyRabinerInverseScalarsFrom( forward_matrices );
            score =
              dp.backward_score( 
                dp_parameters,
                dna_profile,
                dna_seqs,
                num_sequences_to_use,
                backward_matrices
              );
            cout << "Yes, again (backward, with scaling), the total score for all sequences: " << score << endl;

            // Test the position-specific coefficients.
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::PositionSpecificSequenceScoreCoefficients coefficients;
            cout << "Coefficients (before) are " << coefficients << endl;
            //[ (A=0.25,C=0.25,G=0.25,T=0.25), constant:-2.63542e-228 ]
      
            //dp_parameters.debug = DEBUG_All;
            // We use this when useRabinerScaling is true:
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::PositionSpecificSequenceScoreCoefficients unscaled_coefficients;
      
#ifdef USE_DEL_IN_DEL_OUT
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
            // TODO: REMOVE
            cout << "calling calculatePositionSpecificSequenceScoreCoefficients(..) with\nprev_forward_row = " << ( *forward_matrices_iterator )[ 2 ] << "\nbackward_row = " << ( *backward_matrices_iterator )[ 2 ] << endl;
            dp.calculatePositionSpecificSequenceScoreCoefficients(
              dp_parameters,
              dna_profile,
              dna_seq_c,
              3,
              //forward_matrices[ 2 ][ 2 ],
              ( *forward_matrices_iterator )[ 2 ],
              //backward_matrices[ 3 ][ 2 ],
              ( *backward_matrices_iterator )[ 2 ],
              coefficients
            );
            if( dp_parameters.useRabinerScaling ) {
              cout << "Scaled coefficients (after) are " << coefficients << endl;

              cout << "Unscaled coefficients (after) are " << coefficients.createUnscaledCopy() << endl;
              // Without Rabiner Scaling:

            } else {
              cout << "Coefficients (after) are " << coefficients << endl;
              // Without Rabiner Scaling:

            }
      
            score =
              coefficients.calculateScore( dna_profile[ 2 ] );
            cout << "And the score for sequence c, using the coefficients for row 3, is c: " << score << endl;
            // 0.000015296 is the missing!
            // TODO: REMOVE
            //exit( 0 );
#endif // USE_DEL_IN_DEL_OUT

#ifdef USE_DEL_IN_DEL_OUT
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
#else
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
#endif // USE_DEL_IN_DEL_OUT .. else ..
            // TODO: REMOVE
            cout << "calling calculatePositionSpecificSequenceScoreCoefficients(..) with\nprev_forward_row = " << ( *forward_matrices_iterator )[ 2 ] << "\nbackward_row = " << ( *backward_matrices_iterator )[ 2 ] << endl;
            dp.calculatePositionSpecificSequenceScoreCoefficients(
              dp_parameters,
              dna_profile,
              dna_seq_c,
              2,
              //forward_matrices[ 1 ][ 2 ],
              ( *forward_matrices_iterator )[ 2 ],
              //backward_matrices[ 2 ][ 2 ],
              ( *backward_matrices_iterator )[ 2 ],
              coefficients
            );
            if( dp_parameters.useRabinerScaling ) {
              cout << "Scaled coefficients (after) are " << coefficients << endl;
              //[ (A=0.761166,C=0,G=2.31413,T=0), constant:0.00129488, inverseScalar:0.00588399 ]
              cout << "Unscaled coefficients (after) are " << coefficients.createUnscaledCopy() << endl;
              // Without Rabiner Scaling:
              //[ (A=0.00447869,C=0,G=0.0136163,T=0), constant:7.61906e-06 ]
            } else {
              cout << "Coefficients (after) are " << coefficients << endl;
              // Without Rabiner Scaling:
              //[ (A=0.00447869,C=0,G=0.0136163,T=0), constant:7.61906e-06 ]
            }
      
            score =
              coefficients.calculateScore( dna_profile[ 1 ] );
            cout << "And the score for sequence c, using the coefficients for row 2, is c: " << score << endl;

            // TODO: REMOVE
            //exit( 0 );
        
            cout << "Coefficients (before) are " << coefficients << endl;
        
#ifdef USE_DEL_IN_DEL_OUT
            forward_matrices_iterator = forward_matrices.begin();
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
            backward_matrices_iterator++;
#else
            forward_matrices_iterator = forward_matrices.begin();
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT .. else ..
            //dp_parameters.debug = DEBUG_All;
            dp.calculatePositionSpecificSequenceScoreCoefficients(
              dp_parameters,
              dna_profile,
              dna_seq_c,
              1,
              //forward_matrices[ 0 ][ 2 ],
              ( *forward_matrices_iterator )[ 2 ],
              //backward_matrices[ 1 ][ 2 ],
              ( *backward_matrices_iterator )[ 2 ],
              coefficients
            );
            if( dp_parameters.useRabinerScaling ) {
              cout << "Scaled coefficients (after) are " << coefficients << endl;
              //[ (A=1.23849,C=0,G=1.32282,T=0), constant:0.000774347, inverseScalar:0.00588399 ]
              cout << "Unscaled coefficients (after) are " << coefficients.createUnscaledCopy() << endl;
              // Without Rabiner Scaling:
              //[ (A=0.00728727,C=0,G=0.00778344,T=0), constant:4.55625e-06 ]
            } else {
              cout << "Coefficients (after) are " << coefficients << endl;
              // Without Rabiner Scaling:
              //[ (A=0.00728727,C=0,G=0.00778344,T=0), constant:4.55625e-06 ]
            }
      
            score =
              coefficients.calculateScore( dna_profile[ 0 ] );
            cout << "And the score for sequence c, using the coefficients for row 1, is c: " << score << endl;
        
            // TODO: REMOVE
            //exit( 0 );
        
            // Now for seq a, jic.
#ifdef USE_DEL_IN_DEL_OUT
            cout << "Coefficients (before) are " << coefficients << endl;
            //[ (A=0,C=0,G=0,T=0), constant:0.00588399 ]

            //dp_parameters.debug = DEBUG_All;
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();

            dp.calculatePositionSpecificSequenceScoreCoefficients(
              dp_parameters,
              dna_profile,
              dna_seq_a,
              3,
              //forward_matrices[ 2 ][ 0 ],
              ( *forward_matrices_iterator )[ 0 ],
              //backward_matrices[ 3 ][ 0 ],
              ( *backward_matrices_iterator )[ 0 ],
              coefficients
            );
            if( dp_parameters.useRabinerScaling ) {
              cout << "Scaled coefficients (after) are " << coefficients << endl;

              cout << "Unscaled coefficients (after) are " << coefficients.createUnscaledCopy() << endl;
              // Without Rabiner Scaling:

            } else {
              cout << "Coefficients (after) are " << coefficients << endl;
              // Without Rabiner Scaling:

            }
      
            score =
              coefficients.calculateScore( dna_profile[ 2 ] );
            cout << "And the score for sequence a, using the coefficients for row 3, is a: " << score << endl;

            // TODO: REMOVE
            //exit( 0 );
#endif // USE_DEL_IN_DEL_OUT

            cout << "Coefficients (before) are " << coefficients << endl;
            //[ (A=0,C=0,G=0,T=0), constant:0.00588399 ]
      
            //dp_parameters.debug = DEBUG_All;
#ifdef USE_DEL_IN_DEL_OUT
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
#else
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
#endif // USE_DEL_IN_DEL_OUT .. else ..
            dp.calculatePositionSpecificSequenceScoreCoefficients(
              dp_parameters,
              dna_profile,
              dna_seq_a,
              2,
              //forward_matrices[ 1 ][ 0 ],
              ( *forward_matrices_iterator )[ 0 ],
              //backward_matrices[ 2 ][ 0 ],
              ( *backward_matrices_iterator )[ 0 ],
              coefficients
            );
            if( dp_parameters.useRabinerScaling ) {
              cout << "Scaled coefficients (after) are " << coefficients << endl;
              //[ (A=0.782014,C=0,G=0,T=0), constant:0.687195, inverseScalar:0.0082863 ]
              cout << "Unscaled coefficients (after) are " << coefficients.createUnscaledCopy() << endl;
              // Without Rabiner Scaling:
              //[ (A=0.00648,C=0,G=0,T=0), constant:0.0056943 ]
            } else {
              cout << "Coefficients (after) are " << coefficients << endl;
              // Without Rabiner Scaling:
              //[ (A=0.00648,C=0,G=0,T=0), constant:0.0056943 ]
            }
      
            score =
              coefficients.calculateScore( dna_profile[ 1 ] );
            cout << "And the score for sequence a, using the coefficients for row 2, is a: " << score << endl;
        
            // TODO: REMOVE
            //exit( 0 );
        
            cout << "Coefficients (before) are " << coefficients << endl;

#ifdef USE_DEL_IN_DEL_OUT
            forward_matrices_iterator = forward_matrices.begin();
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
            backward_matrices_iterator++;
#else
            forward_matrices_iterator = forward_matrices.begin();
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT .. else ..
            //dp_parameters.debug = DEBUG_All;
            dp.calculatePositionSpecificSequenceScoreCoefficients(
              dp_parameters,
              dna_profile,
              dna_seq_a,
              1,
              //forward_matrices[ 0 ][ 0 ],
              ( *forward_matrices_iterator )[ 0 ],
              //backward_matrices[ 1 ][ 0 ],
              ( *backward_matrices_iterator )[ 0 ],
              coefficients
            );
            if( dp_parameters.useRabinerScaling ) {
              cout << "Scaled coefficients (after) are " << coefficients << endl;
              //[ (A=0.967742,C=0,G=0,T=0), constant:0.322581, inverseScalar:0.0082863 ]
              cout << "Unscaled coefficients (after) are " << coefficients.createUnscaledCopy() << endl;
              // Without Rabiner Scaling:
              //[ (A=0.008019,C=0,G=0,T=0), constant:0.002673 ]
            } else {
              cout << "Coefficients (after) are " << coefficients << endl;
              // Without Rabiner Scaling:
              //[ (A=0.008019,C=0,G=0,T=0), constant:0.002673 ]
            }
      
            score =
              coefficients.calculateScore( dna_profile[ 0 ] );
            cout << "And the score for sequence a, using the coefficients for row 1, is a: " << score << endl;

            // TODO: REMOVE
            //exit( 0 );
        
            // Okay now test getting the Baum-Welch update using the coefficients.
            // Do it without rabiner scaling for now.
            dp_parameters.useRabinerScaling = false; // TODO: REMOVE
            // Have to recalculate these matrices since useRabinerScaling is now false:
            if( dp_parameters.useRabinerScaling == false ) {
              forward_matrices.reinitialize(
                dna_profile.length(),
                dna_seqs,
                num_sequences_to_use
              );
              assert( forward_matrices.size() == ( dna_profile.length() + 1  ) );
              score =
                dp.forward_score( 
                  dp_parameters,
                  dna_profile,
                  dna_seqs,
                  num_sequences_to_use,
                  forward_matrices
                );
              cout << "And yet again (forward), the total score for all sequences: " << score << endl;
              backward_matrices.reinitialize(
                dna_profile,
                dna_seqs,
                num_sequences_to_use
              );
              backward_matrices.copyRabinerInverseScalarsFrom( forward_matrices );
              score =
                dp.backward_score( 
                  dp_parameters,
                  dna_profile,
                  dna_seqs,
                  num_sequences_to_use,
                  backward_matrices
                );
              cout << "Yes, again (backward), the total score for all sequences: " << score << endl;
            } // End if useRabinerScaling == false
          } // End if test_dp_basics
          if( test_dp_basics || test_priors ) {
            position_entente.zero();

#ifdef USE_DEL_IN_DEL_OUT
            forward_matrices_iterator = forward_matrices.begin();
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
            backward_matrices_iterator++;
#else
            forward_matrices_iterator = forward_matrices.begin();
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT .. else ..
            dp.calculatePositionSpecificSequenceScoreCoefficients(
              dp_parameters,
              dna_profile,
              dna_seq_a,
              1,
              //forward_matrices[ 0 ][ 0 ],
              ( *forward_matrices_iterator )[ 0 ],
              //backward_matrices[ 1 ][ 0 ],
              ( *backward_matrices_iterator )[ 0 ],
              coefficients_vector[ 0 ]
            );
            cout << "Scaled coefficients for sequence a, row 1, are " << coefficients_vector[ 0 ] << endl;
      
            cout << "Unscaled coefficients for sequence a, row 1, are " << coefficients_vector[ 0 ].createUnscaledCopy() << endl;
            //[ (A=0.008019,C=0,G=0,T=0), constant:0.002673, inverseScalar:1 ]
            score =
              dp.updatePositionEntenteForSequence(
                dp_parameters,
                dna_profile[ 0 ],
                coefficients_vector[ 0 ],
                NULL,
                position_entente
              );
            cout << "The score, calculated via updatePositionEntenteForSequence at row 1, for sequence a: " << score << endl;
            cout << "The position entente, after seq a, is " << position_entente << endl;
            position_entente_unscaled = position_entente;
            position_entente_unscaled.unscale();
            cout << "Unscaled, the position entente, after seq a, is " << position_entente_unscaled << endl;
            //[ M_dna:(A=0.677419,C=0,G=0,T=0) ]
      
#ifdef USE_DEL_IN_DEL_OUT
            forward_matrices_iterator = forward_matrices.begin();
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
            backward_matrices_iterator++;
#else
            forward_matrices_iterator = forward_matrices.begin();
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT .. else ..
            dp.calculatePositionSpecificSequenceScoreCoefficients(
              dp_parameters,
              dna_profile,
              dna_seq_b,
              1,
              //forward_matrices[ 0 ][ 1 ],
              ( *forward_matrices_iterator )[ 1 ],
              //backward_matrices[ 1 ][ 1 ],
              ( *backward_matrices_iterator )[ 1 ],
              coefficients_vector[ 1 ]
            );
            cout << "Scaled coefficients for sequence b, row 1, are " << coefficients_vector[ 1 ] << endl;
      
            cout << "Unscaled coefficients for sequence b, row 1, are " << coefficients_vector[ 1 ].createUnscaledCopy() << endl;
            // [ (A=0.233553,C=0,G=0.000200475,T=0), constant:0.000116438, inverseScalar:1 ]
            // [ (A=-1.45434,C=-inf,G=-8.51482,T=-inf), constant:-9.05816, inverseScalar:0 ]
            score =
              dp.updatePositionEntenteForSequence(
                dp_parameters,
                dna_profile[ 0 ],
                coefficients_vector[ 1 ],
                NULL,
                position_entente
              );
            cout << "The score, calculated via updatePositionEntenteForSequence at row 1, for sequence b: " << score << endl;
            cout << "The position entente, after seq b, is " << position_entente << endl;
            position_entente_unscaled = position_entente;
            position_entente_unscaled.unscale();
            cout << "Unscaled, the position entente, after seq b, is " << position_entente_unscaled << endl;
            //[ M_dna:(A=1.67659,C=0,G=0.000122522,T=0) ]
      
#ifdef USE_DEL_IN_DEL_OUT
            forward_matrices_iterator = forward_matrices.begin();
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
            backward_matrices_iterator++;
#else
            forward_matrices_iterator = forward_matrices.begin();
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT .. else ..
            dp.calculatePositionSpecificSequenceScoreCoefficients(
              dp_parameters,
              dna_profile,
              dna_seq_c,
              1,
              //forward_matrices[ 0 ][ 2 ],
              ( *forward_matrices_iterator )[ 2 ],
              //backward_matrices[ 1 ][ 2 ],
              ( *backward_matrices_iterator )[ 2 ],
              coefficients_vector[ 2 ]
            );
            cout << "Scaled coefficients for sequence c, row 1, are " << coefficients_vector[ 2 ] << endl;
      
            cout << "Unscaled coefficients for sequence c, row 1, are " << coefficients_vector[ 2 ].createUnscaledCopy() << endl;
            //[ (A=0.00728727,C=0,G=0.00778344,T=0), constant:4.55625e-06, inverseScalar:1 ]
            score =
              dp.updatePositionEntenteForSequence(
                dp_parameters,
                dna_profile[ 0 ],
                coefficients_vector[ 2 ],
                NULL,
                position_entente
              );
            cout << "The score, calculated via updatePositionEntenteForSequence at row 1, for sequence c: " << score << endl;
            cout << "The position entente, after seq c, is " << position_entente << endl;
            position_entente_unscaled = position_entente;
            position_entente_unscaled.unscale();
            cout << "Unscaled, the position entente, after seq c, is " << position_entente_unscaled << endl;
            //[ M_dna:(A=2.54353,C=0,G=0.132404,T=0) ]
      
            position_entente_backup = position_entente;
            position_entente.normalize( 0 );
            cout << "The position entente, normalized (min=0), is " << position_entente << endl;
            //[ M_dna:(A=0.95052,C=0,G=0.0494797,T=0) ]
            position_entente.normalize( 1e-5 );
            cout << "The position entente, normalized (min=1e-5), is " << position_entente << endl;
            //[ M_dna:(A=0.950492,C=1e-05,G=0.0494877,T=1e-05) ]
            position_entente = position_entente_backup;

            // TODO: ERE I AM.  TODO: Add globalentente tests.. (For now: see below, in test_alignment_profiles)
          } // End if test_dp_basics || test_priors
  
          if( test_alignment_profiles ) {
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer::iterator prev_forward_matrices_iterator;
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer::reverse_iterator next_backward_matrices_iterator;
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile ap_c( dna_profile.length() + 1 );
            // Pos 0
#ifdef DISALLOW_FLANKING_TRANSITIONS
            ap_c[ 0 ][
              galosh::Transition::fromPreAlign
            ].zero();
            ap_c[ 0 ][
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
            backward_matrices_iterator++;
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            next_backward_matrices_iterator = backward_matrices_iterator;
            next_backward_matrices_iterator--;
            score =
              dp.calculateAlignmentProfilePosition(
                dp_parameters,
                dna_profile,
                dna_seq_c,
                0,
                //forward_matrices[ 0 ][ 2 ], // ignored
                ( *forward_matrices_iterator )[ 2 ], // ignored
                //forward_matrices[ 0 ][ 2 ],
                ( *forward_matrices_iterator )[ 2 ],
                //backward_matrices[ 0 ][ 2 ],
                ( *backward_matrices_iterator )[ 2 ],
                //backward_matrices[ 1 ][ 2 ],
                ( *next_backward_matrices_iterator )[ 2 ],
                NULL,
                ap_c[ 0 ],
                NULL
              );
            cout << "Score for sequence c, calculated using calculateAlignmentProfilePosition(..) for position 0, is c: " << score << endl;
            cout << "Alignment Profile Position 0, for seq c, is " << ap_c[ 0 ] << endl;
            //Alignment Profile Position 0, for seq c, is [ M:(A=0,C=0,G=0,T=0), M->(M=0.1322818,I=0.03472614,D=0.0002796255), I->(M=0,I=0), D->(M=0,D=0), I:(A=0.03475004,C=0,G=0.0371161,T=0), N->(N=0.03714,B=1), B->(M=0.9992257,D=0.0007743474), C->(C=0,T=0), scalar:0.0491615 ]
            if( dp_parameters.useRabinerScaling == false ) {
              ap_c[ 0 ].unscale();
              cout << "Unscaled, Alignment Profile Position 0, for seq c, is " << ap_c[ 0 ] << endl;
              //Unscaled, Alignment Profile Position 0, for seq c, is [ M:(A=0,C=0,G=0,T=0), M->(M=2.690759,I=0.7063685,D=0.005687895), I->(M=0,I=0), D->(M=0,D=0), I:(A=0.7068546,C=0,G=0.754983,T=0), N->(N=0.7554691,B=20.34112), B->(M=20.32537,D=0.01575109), C->(C=0,T=0) ]
            }

            // TODO: REMOVE
            //exit( 0 );

            // Pos 1
#ifdef DISALLOW_FLANKING_TRANSITIONS
            ap_c[ 1 ][
              galosh::Transition::fromPreAlign
            ].zero();
            ap_c[ 1 ][
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            prev_forward_matrices_iterator = forward_matrices_iterator;
            prev_forward_matrices_iterator--;
            next_backward_matrices_iterator = backward_matrices_iterator;
            next_backward_matrices_iterator--;
            score =
              dp.calculateAlignmentProfilePosition(
                dp_parameters,
                dna_profile,
                dna_seq_c,
                1,
                //forward_matrices[ 0 ][ 2 ],
                ( *prev_forward_matrices_iterator )[ 2 ],
                //forward_matrices[ 1 ][ 2 ],
                ( *forward_matrices_iterator )[ 2 ],
                //backward_matrices[ 1 ][ 2 ],
                ( *backward_matrices_iterator )[ 2 ],
                //backward_matrices[ 2 ][ 2 ],
                ( *next_backward_matrices_iterator )[ 2 ],
                NULL,
                ap_c[ 1 ],
                &( coefficients_vector[ 2 ] )
              );
            cout << "Score for sequence c, calculated using calculateAlignmentProfilePosition(..) for position 1, is c: " << score << endl;
            cout << "Scaled coefficients for sequence c, row 1, are " << coefficients_vector[ 2 ] << endl;
            cout << "Unscaled coefficients for sequence c, row 1, are " << coefficients_vector[ 2 ].createUnscaledCopy() << endl;
            cout << "Alignment Profile Position 1, for seq c, is " << ap_c[ 1 ] << endl;
            // Alignment Profile Position 1, for seq c, is [ M_DNA:(A=0.866944,C=0,G=0.132282,T=0), M->(M->M=0.826229,M->I=0.171719,M->D=0.00127767), I->(I->M=0.171719,I->I=0), D->(D->M=0.00075714,D->D=1.72077e-05), I_DNA:(A=0,C=0,G=0.171719,T=0), N->(N->N=0,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=0,E->J=0), C->(C->C=0,C->T=0), C_DNA:(A=0,C=0,G=0,T=0) ]
            if( dp_parameters.useRabinerScaling == false ) {
              ap_c[ 1 ].unscale();
              cout << "Unscaled, Alignment Profile Position 1, for seq c, is " << ap_c[ 1 ] << endl;
            }

            // TODO: REMOVE
            //exit( 0 );

            // Pos 2
#ifdef DISALLOW_FLANKING_TRANSITIONS
            ap_c[ 2 ][
              galosh::Transition::fromPreAlign
            ].zero();
            ap_c[ 2 ][
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            prev_forward_matrices_iterator = forward_matrices_iterator;
            prev_forward_matrices_iterator--;
            next_backward_matrices_iterator = backward_matrices_iterator;
#ifdef USE_DEL_IN_DEL_OUT
            next_backward_matrices_iterator--;
#endif // USE_DEL_IN_DEL_OUT
            score =
              dp.calculateAlignmentProfilePosition(
                dp_parameters,
                dna_profile,
                dna_seq_c,
                2,
                //forward_matrices[ 1 ][ 2 ],
                ( *prev_forward_matrices_iterator )[ 2 ],
                //forward_matrices[ 2 ][ 2 ],
                ( *forward_matrices_iterator )[ 2 ],
                //backward_matrices[ 2 ][ 2 ],
                ( *backward_matrices_iterator )[ 2 ],
                //backward_matrices[ 3 ][ 2 ], // ignored unless USE_DEL_IN_DEL_OUT
                ( *next_backward_matrices_iterator )[ 2 ], // ignored unless USE_DEL_IN_DEL_OUT
                NULL,
                ap_c[ 2 ],
                &( coefficients_vector[ 2 ] )
              );
            cout << "Score for sequence c, calculated using calculateAlignmentProfilePosition(..) for position 2, is c: " << score << endl;
            cout << "Scaled coefficients for sequence c, row 2, are " << coefficients_vector[ 2 ] << endl;
            cout << "Unscaled coefficients for sequence c, row 2, are " << coefficients_vector[ 2 ].createUnscaledCopy() << endl;
            cout << "Alignment Profile Position 2, for seq c, is " << ap_c[ 2 ] << endl;
            // Alignment Profile Position 2, for seq c, is [ M_DNA:(A=0.304467,C=0,G=0.694239,T=0), M->(M->M=0,M->I=0,M->D=0), I->(I->M=0,I->I=0), D->(D->M=0,D->D=0), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=1,E->J=0), C->(C->C=0.696093,C->T=1), C_DNA:(A=0.695213,C=0,G=0.000880175,T=0) ]
            if( dp_parameters.useRabinerScaling == false ) {
              ap_c[ 2 ].unscale();
              cout << "Unscaled, Alignment Profile Position 2, for seq c, is " << ap_c[ 2 ] << endl;
            }

#ifdef USE_DEL_IN_DEL_OUT
            // Pos 3
#ifdef DISALLOW_FLANKING_TRANSITIONS
            ap_c[ 3 ][
              galosh::Transition::fromPreAlign
            ].zero();
            ap_c[ 3 ][
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
            prev_forward_matrices_iterator = forward_matrices_iterator;
            prev_forward_matrices_iterator--;
            score =
              dp.calculateAlignmentProfilePosition(
                dp_parameters,
                dna_profile,
                dna_seq_c,
                3,
                //forward_matrices[ 2 ][ 2 ],
                ( *prev_forward_matrices_iterator )[ 2 ],
                //forward_matrices[ 3 ][ 2 ],
                ( *forward_matrices_iterator )[ 2 ],
                //backward_matrices[ 3 ][ 2 ],
                ( *backward_matrices_iterator )[ 2 ],
                //backward_matrices[ 3 ][ 2 ], // ignored
                ( *backward_matrices_iterator )[ 2 ], // ignored
                NULL,
                ap_c[ 3 ],
                &( coefficients_vector[ 2 ] )
              );
            cout << "Score for sequence c, calculated using calculateAlignmentProfilePosition(..) for position 3, is c: " << score << endl;
            cout << "Scaled coefficients for sequence c, row 3, are " << coefficients_vector[ 2 ] << endl;
            cout << "Unscaled coefficients for sequence c, row 3, are " << coefficients_vector[ 2 ].createUnscaledCopy() << endl;
            cout << "Alignment Profile Position 3, for seq c, is " << ap_c[ 3 ] << endl;

            if( dp_parameters.useRabinerScaling == false ) {
              ap_c[ 3 ].unscale();
              cout << "Unscaled, Alignment Profile Position 3, for seq c, is " << ap_c[ 3 ] << endl;
            }
#endif // USE_DEL_IN_DEL_OUT
      
            cout << "The alignment profile for seq c prints itself out like so: " << endl;
            cout << ap_c << endl;
            // [ M_DNA:(A=0,C=0,G=0,T=0), M->(M->M=0,M->I=0,M->D=0), I->(I->M=0,I->I=0), D->(D->M=0,D->D=0), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0.134257,N->B=1), B->(B->M=0.999226,B->D=0.000774347), N_DNA:(A=0.133377,C=0,G=0.000880175,T=0), E->(E->C=0,E->J=0), C->(C->C=0,C->T=0), C_DNA:(A=0,C=0,G=0,T=0) ]
            // [ M_DNA:(A=0.866944,C=0,G=0.132282,T=0), M->(M->M=0.826229,M->I=0.171719,M->D=0.00127767), I->(I->M=0.171719,I->I=0), D->(D->M=0.00075714,D->D=1.72077e-05), I_DNA:(A=0,C=0,G=0.171719,T=0), N->(N->N=0,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=0,E->J=0), C->(C->C=0,C->T=0), C_DNA:(A=0,C=0,G=0,T=0) ]
            // [ M_DNA:(A=0.304467,C=0,G=0.694239,T=0), M->(M->M=0,M->I=0,M->D=0), I->(I->M=0,I->I=0), D->(D->M=0,D->D=0), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=1,E->J=0), C->(C->C=0.696093,C->T=1), C_DNA:(A=0.695213,C=0,G=0.000880175,T=0) ]

            // TODO: REMOVE
            //exit( 0 );

            // For comparison, calculate the position ententes.
            forward_matrices_iterator = forward_matrices.begin();
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            dp.calculatePositionSpecificSequenceScoreCoefficients(
              dp_parameters,
              dna_profile,
              dna_seq_c,
              1,
              //forward_matrices[ 0 ][ 2 ],
              ( *forward_matrices_iterator )[ 2 ],
              //backward_matrices[ 1 ][ 2 ],
              ( *backward_matrices_iterator )[ 2 ],
              coefficients_vector[ 2 ]
            );
            position_entente.zero();
            dp.updatePositionEntenteForSequence(
              dp_parameters,
              dna_profile[ 0 ],
              coefficients_vector[ 2 ],
              NULL,
              position_entente
            );
            position_entente.unscale();
            cout << "The position entente at Pos 1, for seq c, is " << position_entente << endl;
            // The position entente at Pos 1, for seq c, is [ M_DNA:(A=0.866944,C=0,G=0.132282,T=0) ]

            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            dp.calculatePositionSpecificSequenceScoreCoefficients(
              dp_parameters,
              dna_profile,
              dna_seq_c,
              2,
              //forward_matrices[ 1 ][ 2 ],
              ( *forward_matrices_iterator )[ 2 ],
              //backward_matrices[ 2 ][ 2 ],
              ( *backward_matrices_iterator )[ 2 ],
              coefficients_vector[ 2 ]
            );
            position_entente.zero();
            dp.updatePositionEntenteForSequence(
              dp_parameters,
              dna_profile[ 1 ],
              coefficients_vector[ 2 ],
              NULL,
              position_entente
            );
            position_entente.unscale();
            cout << "The position entente at Pos 2, for seq c, is " << position_entente << endl;
            // The position entente at Pos 2, for seq c, is [ M_DNA:(A=0.304467,C=0,G=0.694239,T=0) ]

            // TODO: REMOVE
            //exit( 0 );

#ifdef USE_DEL_IN_DEL_OUT
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
            dp.calculatePositionSpecificSequenceScoreCoefficients(
              dp_parameters,
              dna_profile,
              dna_seq_c,
              3,
              //forward_matrices[ 2 ][ 2 ],
              ( *forward_matrices_iterator )[ 2 ],
              //backward_matrices[ 3 ][ 2 ],
              ( *backward_matrices_iterator )[ 2 ],
              coefficients_vector[ 2 ]
            );
            position_entente.zero();
            dp.updatePositionEntenteForSequence(
              dp_parameters,
              dna_profile[ 2 ],
              coefficients_vector[ 2 ],
              NULL,
              position_entente
            );
            position_entente.unscale();
            cout << "The position entente at Pos 3, for seq c, is " << position_entente << endl;

            // TODO: REMOVE
            //exit( 0 );
#endif // USE_DEL_IN_DEL_OUT

            // For comparison, calculate the global ententes.
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::GlobalEntente global_entente;
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::GlobalEntente global_entente_backup;
            global_entente.zero();
#ifdef DISALLOW_FLANKING_TRANSITIONS
            global_entente[
              galosh::Transition::fromPreAlign
            ].zero();
            global_entente[
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
            backward_matrices_iterator++;
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            score =
              dp.updateGlobalEntenteForSequence(
                dp_parameters,
                dna_profile,
                dna_seq_c,
                0,
                //forward_matrices[ 0 ][ 2 ], // ignored
                ( *forward_matrices_iterator )[ 2 ], // ignored
                //forward_matrices[ 0 ][ 2 ],
                ( *forward_matrices_iterator )[ 2 ],
                //backward_matrices[ 0 ][ 2 ],
                ( *backward_matrices_iterator )[ 2 ],
                NULL,
                global_entente
              );
            cout << "Score for sequence c, calculated using updateGlobalEntenteForSequence(..) for position 0, is c: " << score << endl;
            global_entente.unscale();
            cout << "The Global Entente at Pos 0, for seq c, is " << global_entente << endl;
            // The Global Entente at Pos 0, for seq c, is [ M->(M->M=0,M->I=0,M->D=0), I->(I->M=0,I->I=0), D->(D->M=0,D->D=0), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0.134257,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0.133377,C=0,G=0.000880175,T=0), E->(E->C=0,E->J=0), C->(C->C=0,C->T=0), C_DNA:(A=0,C=0,G=0,T=0) ]

            // TODO: REMOVE
            //exit( 0 );

            global_entente.zero();
#ifdef DISALLOW_FLANKING_TRANSITIONS
            global_entente[
              galosh::Transition::fromPreAlign
            ].zero();
            global_entente[
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            prev_forward_matrices_iterator = forward_matrices_iterator;
            prev_forward_matrices_iterator--;
            score =
              dp.updateGlobalEntenteForSequence(
                dp_parameters,
                dna_profile,
                dna_seq_c,
                1,
                //forward_matrices[ 0 ][ 2 ],
                ( *prev_forward_matrices_iterator )[ 2 ],
                //forward_matrices[ 1 ][ 2 ],
                ( *forward_matrices_iterator )[ 2 ],
                //backward_matrices[ 1 ][ 2 ],
                ( *backward_matrices_iterator )[ 2 ],
                NULL,
                global_entente
              );
            cout << "Score for sequence c, calculated using updateGlobalEntenteForSequence(..) for position 1, is c: " << score << endl;
            global_entente.unscale();
            cout << "The Global Entente at Pos 1, for seq c, is " << global_entente << endl;
            // The Global Entente at Pos 1, for seq c, is [ M->(M->M=0,M->I=0.171719,M->D=0), I->(I->M=0,I->I=0), D->(D->M=0,D->D=0), I_DNA:(A=0,C=0,G=0.171719,T=0), N->(N->N=0,N->B=1), B->(B->M=0.999226,B->D=0.000774347), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=0,E->J=0), C->(C->C=0,C->T=0), C_DNA:(A=0,C=0,G=0,T=0) ]

            global_entente.zero();
#ifdef DISALLOW_FLANKING_TRANSITIONS
            global_entente[
              galosh::Transition::fromPreAlign
            ].zero();
            global_entente[
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            prev_forward_matrices_iterator = forward_matrices_iterator;
            prev_forward_matrices_iterator--;
            score =
              dp.updateGlobalEntenteForSequence(
                dp_parameters,
                dna_profile,
                dna_seq_c,
                2,
                //forward_matrices[ 1 ][ 2 ],
                ( *prev_forward_matrices_iterator )[ 2 ],
                //forward_matrices[ 2 ][ 2 ],
                ( *forward_matrices_iterator )[ 2 ],
                //backward_matrices[ 2 ][ 2 ],
                ( *backward_matrices_iterator )[ 2 ],
                NULL,
                global_entente
              );
            cout << "Score for sequence c, calculated using updateGlobalEntenteForSequence(..) for position 2, is c: " << score << endl;
            global_entente.unscale();
            cout << "The Global Entente at Pos 2, for seq c, is " << global_entente << endl;
            // The Global Entente at Pos 2, for seq c, is [ M->(M->M=0.826229,M->I=0,M->D=0.00127767), I->(I->M=0.171719,I->I=0), D->(D->M=0.00075714,D->D=1.72077e-05), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=1,E->J=0), C->(C->C=0.696093,C->T=1), C_DNA:(A=0.695213,C=0,G=0.000880175,T=0) ]

#ifdef USE_DEL_IN_DEL_OUT
            global_entente.zero();
#ifdef DISALLOW_FLANKING_TRANSITIONS
            global_entente[
              galosh::Transition::fromPreAlign
            ].zero();
            global_entente[
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
            prev_forward_matrices_iterator = forward_matrices_iterator;
            prev_forward_matrices_iterator--;
            score =
              dp.updateGlobalEntenteForSequence(
                dp_parameters,
                dna_profile,
                dna_seq_c,
                3,
                //forward_matrices[ 2 ][ 2 ],
                ( *prev_forward_matrices_iterator )[ 2 ],
                //forward_matrices[ 3 ][ 2 ],
                ( *forward_matrices_iterator )[ 2 ],
                //backward_matrices[ 3 ][ 2 ],
                ( *backward_matrices_iterator )[ 2 ],
                NULL,
                global_entente
              );
            cout << "Score for sequence c, calculated using updateGlobalEntenteForSequence(..) for position 3, is c: " << score << endl;
            global_entente.unscale();
            cout << "The Global Entente at Pos 3, for seq c, is " << global_entente << endl;

            // TODO: REMOVE
            //exit( 0 );
      
#endif // USE_DEL_IN_DEL_OUT
            cout << endl;

            // Ok now do seq b
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile ap_b( ( dna_profile.length() + 1 ) );
            // Pos 0
#ifdef DISALLOW_FLANKING_TRANSITIONS
            ap_b[ 0 ][
              galosh::Transition::fromPreAlign
            ].zero();
            ap_b[ 0 ][
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
            backward_matrices_iterator++;
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            next_backward_matrices_iterator = backward_matrices_iterator;
            next_backward_matrices_iterator--;
            score =
              dp.calculateAlignmentProfilePosition(
                dp_parameters,
                dna_profile,
                dna_seq_b,
                0,
                //forward_matrices[ 0 ][ 1 ], // ignored
                ( *forward_matrices_iterator )[ 1 ], // ignored
                //forward_matrices[ 0 ][ 1 ],
                ( *forward_matrices_iterator )[ 1 ],
                //backward_matrices[ 0 ][ 1 ],
                ( *backward_matrices_iterator )[ 1 ],
                //backward_matrices[ 1 ][ 1 ],
                ( *next_backward_matrices_iterator )[ 1 ],
                NULL,
                ap_b[ 0 ],
                NULL
              );
            cout << "Score for sequence b, calculated using calculateAlignmentProfilePosition(..) for position 0, is b: " << score << endl;
            //cout << "Alignment Profile Position 0, for seq b, is " << ap_b[ 0 ] << endl;
            if( dp_parameters.useRabinerScaling == false ) {
              ap_b[ 0 ].unscale();
            //  cout << "Unscaled, Alignment Profile Position 0, for seq b, is " << ap_b[ 0 ] << endl;
            }

            // Pos 1
#ifdef DISALLOW_FLANKING_TRANSITIONS
            ap_b[ 1 ][
              galosh::Transition::fromPreAlign
            ].zero();
            ap_b[ 1 ][
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            prev_forward_matrices_iterator = forward_matrices_iterator;
            prev_forward_matrices_iterator--;
            next_backward_matrices_iterator = backward_matrices_iterator;
            next_backward_matrices_iterator--;
            score =
              dp.calculateAlignmentProfilePosition(
                dp_parameters,
                dna_profile,
                dna_seq_b,
                1,
                //forward_matrices[ 0 ][ 1 ],
                ( *prev_forward_matrices_iterator )[ 1 ],
                //forward_matrices[ 1 ][ 1 ],
                ( *forward_matrices_iterator )[ 1 ],
                //backward_matrices[ 1 ][ 1 ],
                ( *backward_matrices_iterator )[ 1 ],
                //backward_matrices[ 2 ][ 1 ],
                ( *next_backward_matrices_iterator )[ 1 ],
                NULL,
                ap_b[ 1 ],
                NULL
              );
            cout << "Score for sequence b, calculated using calculateAlignmentProfilePosition(..) for position 1, is b: " << score << endl;
            //cout << "Alignment Profile Position 1, for seq b, is " << ap_b[ 1 ] << endl;
            if( dp_parameters.useRabinerScaling == false ) {
              ap_b[ 1 ].unscale();
            //  cout << "Unscaled, Alignment Profile Position 1, for seq b, is " << ap_b[ 1 ] << endl;
            }

            // Pos 2
#ifdef DISALLOW_FLANKING_TRANSITIONS
            ap_b[ 2 ][
              galosh::Transition::fromPreAlign
            ].zero();
            ap_b[ 2 ][
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            prev_forward_matrices_iterator = forward_matrices_iterator;
            prev_forward_matrices_iterator--;
            next_backward_matrices_iterator = backward_matrices_iterator;
#ifdef USE_DEL_IN_DEL_OUT
            next_backward_matrices_iterator--;
#endif // USE_DEL_IN_DEL_OUT
            score =
              dp.calculateAlignmentProfilePosition(
                dp_parameters,
                dna_profile,
                dna_seq_b,
                2,
                //forward_matrices[ 1 ][ 1 ],
                ( *prev_forward_matrices_iterator )[ 1 ],
                //forward_matrices[ 2 ][ 1 ],
                ( *forward_matrices_iterator )[ 1 ],
                //backward_matrices[ 2 ][ 1 ],
                ( *backward_matrices_iterator )[ 1 ],
                //backward_matrices[ 2 ][ 1 ], // ignored
                ( *next_backward_matrices_iterator )[ 1 ], // ignored
                NULL,
                ap_b[ 2 ],
                NULL
              );
            cout << "Score for sequence b, calculated using calculateAlignmentProfilePosition(..) for position 2, is b: " << score << endl;
            //cout << "Alignment Profile Position 2, for seq b, is " << ap_b[ 2 ] << endl;
            if( dp_parameters.useRabinerScaling == false ) {
              ap_b[ 2 ].unscale();
            //  cout << "Unscaled, Alignment Profile Position 2, for seq b, is " << ap_b[ 2 ] << endl;
            }

#ifdef USE_DEL_IN_DEL_OUT
            // Pos 3
#ifdef DISALLOW_FLANKING_TRANSITIONS
            ap_b[ 3 ][
              galosh::Transition::fromPreAlign
            ].zero();
            ap_b[ 3 ][
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
            prev_forward_matrices_iterator = forward_matrices_iterator;
            prev_forward_matrices_iterator--;
            score =
              dp.calculateAlignmentProfilePosition(
                dp_parameters,
                dna_profile,
                dna_seq_b,
                3,
                //forward_matrices[ 2 ][ 1 ],
                ( *prev_forward_matrices_iterator )[ 1 ],
                //forward_matrices[ 3 ][ 1 ],
                ( *forward_matrices_iterator )[ 1 ],
                //backward_matrices[ 3 ][ 1 ],
                ( *backward_matrices_iterator )[ 1 ],
                //backward_matrices[ 3 ][ 1 ], // ignored
                ( *backward_matrices_iterator )[ 1 ], // ignored
                NULL,
                ap_b[ 3 ],
                NULL
              );
            cout << "Score for sequence b, calculated using calculateAlignmentProfilePosition(..) for position 3, is b: " << score << endl;
            //cout << "Alignment Profile Position 3, for seq b, is " << ap_b[ 3 ] << endl;
            if( dp_parameters.useRabinerScaling == false ) {
              ap_b[ 3 ].unscale();
            //  cout << "Unscaled, Alignment Profile Position 3, for seq b, is " << ap_b[ 3 ] << endl;
            }
#endif // USE_DEL_IN_DEL_OUT

            cout << "The alignment profile for seq b prints itself out like so: " << endl;
            cout << ap_b << endl;
            // [ M_DNA:(A=0,C=0,G=0,T=0), M->(M->M=0,M->I=0,M->D=0), I->(I->M=0,I->I=0), D->(D->M=0,D->D=0), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0.000438109,N->B=1), B->(B->M=0.999288,B->D=0.000711617), N_DNA:(A=0.000431921,C=0,G=6.18797e-06,T=0), E->(E->C=0,E->J=0), C->(C->C=0,C->T=0), C_DNA:(A=0,C=0,G=0,T=0) ]
            // [ M_DNA:(A=0.999166,C=0,G=0.000122522,T=0), M->(M->M=0.998308,M->I=0,M->D=0.000980175), I->(I->M=0,I->I=0), D->(D->M=0.000693053,D->D=1.85639e-05), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=0,E->J=0), C->(C->C=0,C->T=0), C_DNA:(A=0,C=0,G=0,T=0) ]
            // [ M_DNA:(A=0.00039603,C=0,G=0.998605,T=0), M->(M->M=0,M->I=0,M->D=0), I->(I->M=0,I->I=0), D->(D->M=0,D->D=0), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=1,E->J=0), C->(C->C=0.00127225,C->T=1), C_DNA:(A=6.18797e-06,C=0,G=0.00126606,T=0) ]

            // TODO: REMOVE
            //exit( 0 );
      
            // For comparison, calculate the position ententes.
            forward_matrices_iterator = forward_matrices.begin();
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            dp.calculatePositionSpecificSequenceScoreCoefficients(
              dp_parameters,
              dna_profile,
              dna_seq_b,
              1,
              //forward_matrices[ 0 ][ 1 ],
              ( *forward_matrices_iterator )[ 1 ],
              //backward_matrices[ 1 ][ 1 ],
              ( *backward_matrices_iterator )[ 1 ],
              coefficients_vector[ 1 ]
            );
            position_entente.zero();
            dp.updatePositionEntenteForSequence(
              dp_parameters,
              dna_profile[ 0 ],
              coefficients_vector[ 1 ],
              NULL,
              position_entente
            );
            position_entente.unscale();
            cout << "The position entente at Pos 1, for seq b, is " << position_entente << endl;
            // The position entente at Pos 1, for seq b, is [ M_DNA:(A=0.999166,C=0,G=0.000122522,T=0) ]
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            dp.calculatePositionSpecificSequenceScoreCoefficients(
              dp_parameters,
              dna_profile,
              dna_seq_b,
              2,
              //forward_matrices[ 1 ][ 1 ],
              ( *forward_matrices_iterator )[ 1 ],
              //backward_matrices[ 2 ][ 1 ],
              ( *backward_matrices_iterator )[ 1 ],
              coefficients_vector[ 1 ]
            );
            position_entente.zero();
            dp.updatePositionEntenteForSequence(
              dp_parameters,
              dna_profile[ 1 ],
              coefficients_vector[ 1 ],
              NULL,
              position_entente
            );
            position_entente.unscale();
            cout << "The position entente at Pos 2, for seq b, is " << position_entente << endl;
            // The position entente at Pos 2, for seq b, is [ M_DNA:(A=0.00039603,C=0,G=0.998605,T=0) ]
#ifdef USE_DEL_IN_DEL_OUT
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
            dp.calculatePositionSpecificSequenceScoreCoefficients(
              dp_parameters,
              dna_profile,
              dna_seq_b,
              3,
              //forward_matrices[ 2 ][ 1 ],
              ( *forward_matrices_iterator )[ 1 ],
              //backward_matrices[ 3 ][ 1 ],
              ( *backward_matrices_iterator )[ 1 ],
              coefficients_vector[ 1 ]
            );
            position_entente.zero();
            dp.updatePositionEntenteForSequence(
              dp_parameters,
              dna_profile[ 2 ],
              coefficients_vector[ 1 ],
              NULL,
              position_entente
            );
            position_entente.unscale();
            cout << "The position entente at Pos 3, for seq b, is " << position_entente << endl;

            // TODO: REMOVE
            //exit( 0 );
#endif // USE_DEL_IN_DEL_OUT
      
            // For comparison, calculate the global ententes.
            global_entente.zero();
#ifdef DISALLOW_FLANKING_TRANSITIONS
            global_entente[
              galosh::Transition::fromPreAlign
            ].zero();
            global_entente[
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
            backward_matrices_iterator++;
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            score =
              dp.updateGlobalEntenteForSequence(
                dp_parameters,
                dna_profile,
                dna_seq_b,
                0,
                //forward_matrices[ 0 ][ 1 ], // ignored
                ( *forward_matrices_iterator )[ 1 ], // ignored
                //forward_matrices[ 0 ][ 1 ],
                ( *forward_matrices_iterator )[ 1 ],
                //backward_matrices[ 0 ][ 1 ],
                ( *backward_matrices_iterator )[ 1 ],
                NULL,
                global_entente
              );
            cout << "Score for sequence b, calculated using updateGlobalEntenteForSequence(..) for position 0, is b: " << score << endl;
            global_entente.unscale();
            cout << "The Global Entente at Pos 0, for seq b, is " << global_entente << endl;
            // The Global Entente at Pos 0, for seq b, is [ M->(M->M=0,M->I=0,M->D=0), I->(I->M=0,I->I=0), D->(D->M=0,D->D=0), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0.000438109,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0.000431921,C=0,G=6.18797e-06,T=0), E->(E->C=0,E->J=0), C->(C->C=0,C->T=0), C_DNA:(A=0,C=0,G=0,T=0) ]
            global_entente.zero();
#ifdef DISALLOW_FLANKING_TRANSITIONS
            global_entente[
              galosh::Transition::fromPreAlign
            ].zero();
            global_entente[
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            prev_forward_matrices_iterator = forward_matrices_iterator;
            prev_forward_matrices_iterator--;
            score =
              dp.updateGlobalEntenteForSequence(
                dp_parameters,
                dna_profile,
                dna_seq_b,
                1,
                //forward_matrices[ 0 ][ 1 ],
                ( *prev_forward_matrices_iterator )[ 1 ],
                //forward_matrices[ 1 ][ 1 ],
                ( *forward_matrices_iterator )[ 1 ],
                //backward_matrices[ 1 ][ 1 ],
                ( *backward_matrices_iterator )[ 1 ],
                NULL,
                global_entente
              );
            cout << "Score for sequence b, calculated using updateGlobalEntenteForSequence(..) for position 1, is b: " << score << endl;
            global_entente.unscale();
            cout << "The Global Entente at Pos 1, for seq b, is " << global_entente << endl;
            // The Global Entente at Pos 1, for seq b, is [ M->(M->M=0,M->I=0,M->D=0), I->(I->M=0,I->I=0), D->(D->M=0,D->D=0), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0,N->B=1), B->(B->M=0.999288,B->D=0.000711617), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=0,E->J=0), C->(C->C=0,C->T=0), C_DNA:(A=0,C=0,G=0,T=0) ]
            global_entente.zero();
#ifdef DISALLOW_FLANKING_TRANSITIONS
            global_entente[
              galosh::Transition::fromPreAlign
            ].zero();
            global_entente[
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            prev_forward_matrices_iterator = forward_matrices_iterator;
            prev_forward_matrices_iterator--;
            score =
              dp.updateGlobalEntenteForSequence(
                dp_parameters,
                dna_profile,
                dna_seq_b,
                2,
                //forward_matrices[ 1 ][ 1 ],
                ( *prev_forward_matrices_iterator )[ 1 ],
                //forward_matrices[ 2 ][ 1 ],
                ( *forward_matrices_iterator )[ 1 ],
                //backward_matrices[ 2 ][ 1 ],
                ( *backward_matrices_iterator )[ 1 ],
                NULL,
                global_entente
              );
            cout << "Score for sequence b, calculated using updateGlobalEntenteForSequence(..) for position 2, is b: " << score << endl;
            global_entente.unscale();
            cout << "The Global Entente at Pos 2, for seq b, is " << global_entente << endl;
            // The Global Entente at Pos 2, for seq b, is [ M->(M->M=0.998308,M->I=0,M->D=0.000980175), I->(I->M=0,I->I=0), D->(D->M=0.000693053,D->D=1.85639e-05), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=1,E->J=0), C->(C->C=0.00127225,C->T=1), C_DNA:(A=6.18797e-06,C=0,G=0.00126606,T=0) ]

#ifdef USE_DEL_IN_DEL_OUT
            global_entente.zero();
#ifdef DISALLOW_FLANKING_TRANSITIONS
            global_entente[
              galosh::Transition::fromPreAlign
            ].zero();
            global_entente[
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
            prev_forward_matrices_iterator = forward_matrices_iterator;
            prev_forward_matrices_iterator--;
            score =
              dp.updateGlobalEntenteForSequence(
                dp_parameters,
                dna_profile,
                dna_seq_b,
                3,
                //forward_matrices[ 2 ][ 1 ],
                ( *prev_forward_matrices_iterator )[ 1 ],
                //forward_matrices[ 3 ][ 1 ],
                ( *forward_matrices_iterator )[ 1 ],
                //backward_matrices[ 3 ][ 1 ],
                ( *backward_matrices_iterator )[ 1 ],
                NULL,
                global_entente
              );
            cout << "Score for sequence b, calculated using updateGlobalEntenteForSequence(..) for position 3, is b: " << score << endl;
            global_entente.unscale();
            cout << "The Global Entente at Pos 3, for seq b, is " << global_entente << endl;
#endif // USE_DEL_IN_DEL_OUT
            cout << endl;

            // TODO: REMOVE
            //exit( 0 );

            // Ok now do seq a
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile ap_a( ( dna_profile.length() + 1 ) );
            // Pos 0
#ifdef DISALLOW_FLANKING_TRANSITIONS
            ap_a[ 0 ][
              galosh::Transition::fromPreAlign
            ].zero();
            ap_a[ 0 ][
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
            backward_matrices_iterator++;
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            next_backward_matrices_iterator = backward_matrices_iterator;
            next_backward_matrices_iterator--;
            score =
              dp.calculateAlignmentProfilePosition(
                dp_parameters,
                dna_profile,
                dna_seq_a,
                0,
                //forward_matrices[ 0 ][ 0 ], // ignored
                ( *forward_matrices_iterator )[ 0 ], // ignored
                //forward_matrices[ 0 ][ 0 ],
                ( *forward_matrices_iterator )[ 0 ],
                //backward_matrices[ 0 ][ 0 ],
                ( *backward_matrices_iterator )[ 0 ],
                //backward_matrices[ 1 ][ 0 ],
                ( *next_backward_matrices_iterator )[ 0 ],
                NULL,
                ap_a[ 0 ],
                NULL
              );
            cout << "Score for sequence a, calculated using calculateAlignmentProfilePosition(..) for position 0, is a: " << score << endl;
            //cout << "Alignment Profile Position 0, for seq a, is " << ap_a[ 0 ] << endl;
            if( dp_parameters.useRabinerScaling == false ) {
              ap_a[ 0 ].unscale();
            //  cout << "Unscaled, Alignment Profile Position 0, for seq a, is " << ap_a[ 0 ] << endl;
            }

            // Pos 1
#ifdef DISALLOW_FLANKING_TRANSITIONS
            ap_a[ 1 ][
              galosh::Transition::fromPreAlign
            ].zero();
            ap_a[ 1 ][
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            prev_forward_matrices_iterator = forward_matrices_iterator;
            prev_forward_matrices_iterator--;
            next_backward_matrices_iterator = backward_matrices_iterator;
            next_backward_matrices_iterator--;
            score =
              dp.calculateAlignmentProfilePosition(
                dp_parameters,
                dna_profile,
                dna_seq_a,
                1,
                //forward_matrices[ 0 ][ 0 ],
                ( *prev_forward_matrices_iterator )[ 0 ],
                //forward_matrices[ 1 ][ 0 ],
                ( *forward_matrices_iterator )[ 0 ],
                //backward_matrices[ 1 ][ 0 ],
                ( *backward_matrices_iterator )[ 0 ],
                //backward_matrices[ 2 ][ 0 ],
                ( *next_backward_matrices_iterator )[ 0 ],
                NULL,
                ap_a[ 1 ],
                NULL
              );
            cout << "Score for sequence a, calculated using calculateAlignmentProfilePosition(..) for position 1, is a: " << score << endl;
            //cout << "Alignment Profile Position 1, for seq a, is " << ap_a[ 1 ] << endl;
            if( dp_parameters.useRabinerScaling == false ) {
              ap_a[ 1 ].unscale();
            //  cout << "Unscaled, Alignment Profile Position 1, for seq a, is " << ap_a[ 1 ] << endl;
            }

            // Pos 2
#ifdef DISALLOW_FLANKING_TRANSITIONS
            ap_a[ 2 ][
              galosh::Transition::fromPreAlign
            ].zero();
            ap_a[ 2 ][
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            prev_forward_matrices_iterator = forward_matrices_iterator;
            prev_forward_matrices_iterator--;
            next_backward_matrices_iterator = backward_matrices_iterator;
#ifdef USE_DEL_IN_DEL_OUT
            next_backward_matrices_iterator--;
#endif // USE_DEL_IN_DEL_OUT
            score =
              dp.calculateAlignmentProfilePosition(
                dp_parameters,
                dna_profile,
                dna_seq_a,
                2,
                //forward_matrices[ 1 ][ 0 ],
                ( *prev_forward_matrices_iterator )[ 0 ],
                //forward_matrices[ 2 ][ 0 ],
                ( *forward_matrices_iterator )[ 0 ],
                //backward_matrices[ 2 ][ 0 ],
                ( *backward_matrices_iterator )[ 0 ],
                //backward_matrices[ 2 ][ 0 ], // ignored
                ( *next_backward_matrices_iterator )[ 0 ], // ignored
                NULL,
                ap_a[ 2 ],
                NULL
              );
            cout << "Score for sequence a, calculated using calculateAlignmentProfilePosition(..) for position 2, is a: " << score << endl;
            //cout << "Alignment Profile Position 2, for seq a, is " << ap_a[ 2 ] << endl;
            if( dp_parameters.useRabinerScaling == false ) {
              ap_a[ 2 ].unscale();
            //  cout << "Unscaled, Alignment Profile Position 2, for seq a, is " << ap_a[ 2 ] << endl;
            }

#ifdef USE_DEL_IN_DEL_OUT
            // Pos 3
#ifdef DISALLOW_FLANKING_TRANSITIONS
            ap_a[ 3 ][
              galosh::Transition::fromPreAlign
            ].zero();
            ap_a[ 3 ][
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
            prev_forward_matrices_iterator = forward_matrices_iterator;
            prev_forward_matrices_iterator--;
            score =
              dp.calculateAlignmentProfilePosition(
                dp_parameters,
                dna_profile,
                dna_seq_a,
                3,
                //forward_matrices[ 2 ][ 0 ],
                ( *prev_forward_matrices_iterator )[ 0 ],
                //forward_matrices[ 3 ][ 0 ],
                ( *forward_matrices_iterator )[ 0 ],
                //backward_matrices[ 3 ][ 0 ],
                ( *backward_matrices_iterator )[ 0 ],
                //backward_matrices[ 3 ][ 0 ], // ignored
                ( *backward_matrices_iterator )[ 0 ], // ignored
                NULL,
                ap_a[ 3 ],
                NULL
              );
            cout << "Score for sequence a, calculated using calculateAlignmentProfilePosition(..) for position 3, is a: " << score << endl;
            //cout << "Alignment Profile Position 3, for seq a, is " << ap_a[ 3 ] << endl;
            if( dp_parameters.useRabinerScaling == false ) {
              ap_a[ 3 ].unscale();
            //  cout << "Unscaled, Alignment Profile Position 3, for seq a, is " << ap_a[ 3 ] << endl;
            }
#endif // USE_DEL_IN_DEL_OUT

            cout << "The alignment profile for seq b prints itself out like so: " << endl;
            cout << ap_a << endl;
            // [ M_DNA:(A=0,C=0,G=0,T=0), M->(M->M=0,M->I=0,M->D=0), I->(I->M=0,I->I=0), D->(D->M=0,D->D=0), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0.00488759,N->B=1), B->(B->M=0.677419,B->D=0.322581), N_DNA:(A=0.00488759,C=0,G=0,T=0), E->(E->C=0,E->J=0), C->(C->C=0,C->T=0), C_DNA:(A=0,C=0,G=0,T=0) ]
            // [ M_DNA:(A=0.677419,C=0,G=0,T=0), M->(M->M=0,M->I=0,M->D=0.677419), I->(I->M=0,I->I=0), D->(D->M=0.312805,D->D=0.00977517), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=0,E->J=0), C->(C->C=0,C->T=0), C_DNA:(A=0,C=0,G=0,T=0) ]
            // [ M_DNA:(A=0.312805,C=0,G=0,T=0), M->(M->M=0,M->I=0,M->D=0), I->(I->M=0,I->I=0), D->(D->M=0,D->D=0), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=1,E->J=0), C->(C->C=0.00488759,C->T=1), C_DNA:(A=0.00488759,C=0,G=0,T=0) ]

            // TODO: REMOVE
            //exit( 0 );
      
            // For comparison, calculate the position ententes.
            forward_matrices_iterator = forward_matrices.begin();
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            dp.calculatePositionSpecificSequenceScoreCoefficients(
              dp_parameters,
              dna_profile,
              dna_seq_a,
              1,
              //forward_matrices[ 0 ][ 0 ],
              ( *forward_matrices_iterator )[ 0 ],
              //backward_matrices[ 1 ][ 0 ],
              ( *backward_matrices_iterator )[ 0 ],
              coefficients_vector[ 0 ]
            );
            position_entente.zero();
            dp.updatePositionEntenteForSequence(
              dp_parameters,
              dna_profile[ 0 ],
              coefficients_vector[ 0 ],
              NULL,
              position_entente
            );
            position_entente.unscale();
            cout << "The position entente at Pos 1, for seq a, is " << position_entente << endl;
            // The position entente at Pos 1, for seq a, is [ M_DNA:(A=0.999166,C=0,G=0.000122522,T=0) ]
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            dp.calculatePositionSpecificSequenceScoreCoefficients(
              dp_parameters,
              dna_profile,
              dna_seq_a,
              2,
              //forward_matrices[ 1 ][ 0 ],
              ( *forward_matrices_iterator )[ 0 ],
              //backward_matrices[ 2 ][ 0 ],
              ( *backward_matrices_iterator )[ 0 ],
              coefficients_vector[ 0 ]
            );
            position_entente.zero();
            dp.updatePositionEntenteForSequence(
              dp_parameters,
              dna_profile[ 1 ],
              coefficients_vector[ 0 ],
              NULL,
              position_entente
            );
            position_entente.unscale();
            cout << "The position entente at Pos 2, for seq a, is " << position_entente << endl;
            // The position entente at Pos 2, for seq a, is [ M_DNA:(A=0.00039603,C=0,G=0.998605,T=0) ]
#ifdef USE_DEL_IN_DEL_OUT
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
            dp.calculatePositionSpecificSequenceScoreCoefficients(
              dp_parameters,
              dna_profile,
              dna_seq_a,
              3,
              //forward_matrices[ 2 ][ 0 ],
              ( *forward_matrices_iterator )[ 0 ],
              //backward_matrices[ 3 ][ 0 ],
              ( *backward_matrices_iterator )[ 0 ],
              coefficients_vector[ 0 ]
            );
            position_entente.zero();
            dp.updatePositionEntenteForSequence(
              dp_parameters,
              dna_profile[ 2 ],
              coefficients_vector[ 0 ],
              NULL,
              position_entente
            );
            position_entente.unscale();
            cout << "The position entente at Pos 3, for seq a, is " << position_entente << endl;

            // TODO: REMOVE
            //exit( 0 );
#endif // USE_DEL_IN_DEL_OUT
      
            // For comparison, calculate the global ententes.
            global_entente.zero();
#ifdef DISALLOW_FLANKING_TRANSITIONS
            global_entente[
              galosh::Transition::fromPreAlign
            ].zero();
            global_entente[
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
            backward_matrices_iterator++;
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            score =
              dp.updateGlobalEntenteForSequence(
                dp_parameters,
                dna_profile,
                dna_seq_a,
                0,
                //forward_matrices[ 0 ][ 0 ], // ignored
                ( *forward_matrices_iterator )[ 0 ], // ignored
                //forward_matrices[ 0 ][ 0 ],
                ( *forward_matrices_iterator )[ 0 ],
                //backward_matrices[ 0 ][ 0 ],
                ( *backward_matrices_iterator )[ 0 ],
                NULL,
                global_entente
              );
            cout << "Score for sequence a, calculated using updateGlobalEntenteForSequence(..) for position 0, is a: " << score << endl;
            global_entente.unscale();
            cout << "The Global Entente at Pos 0, for seq a, is " << global_entente << endl;
            // The Global Entente at Pos 0, for seq a, is [ M->(M->M=0,M->I=0,M->D=0), I->(I->M=0,I->I=0), D->(D->M=0,D->D=0), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0.000438109,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0.000431921,C=0,G=6.18797e-06,T=0), E->(E->C=0,E->J=0), C->(C->C=0,C->T=0), C_DNA:(A=0,C=0,G=0,T=0) ]
            global_entente.zero();
#ifdef DISALLOW_FLANKING_TRANSITIONS
            global_entente[
              galosh::Transition::fromPreAlign
            ].zero();
            global_entente[
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
            backward_matrices_iterator++;
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            prev_forward_matrices_iterator = forward_matrices_iterator;
            prev_forward_matrices_iterator--;
            score =
              dp.updateGlobalEntenteForSequence(
                dp_parameters,
                dna_profile,
                dna_seq_a,
                1,
                //forward_matrices[ 0 ][ 0 ],
                ( *prev_forward_matrices_iterator )[ 0 ],
                //forward_matrices[ 1 ][ 0 ],
                ( *forward_matrices_iterator )[ 0 ],
                //backward_matrices[ 1 ][ 0 ],
                ( *backward_matrices_iterator )[ 0 ],
                NULL,
                global_entente
              );
            cout << "Score for sequence a, calculated using updateGlobalEntenteForSequence(..) for position 1, is a: " << score << endl;
            global_entente.unscale();
            cout << "The Global Entente at Pos 1, for seq a, is " << global_entente << endl;
            // The Global Entente at Pos 1, for seq a, is [ M->(M->M=0,M->I=0,M->D=0), I->(I->M=0,I->I=0), D->(D->M=0,D->D=0), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0,N->B=1), B->(B->M=0.999288,B->D=0.000711617), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=0,E->J=0), C->(C->C=0,C->T=0), C_DNA:(A=0,C=0,G=0,T=0) ]
            global_entente.zero();
#ifdef DISALLOW_FLANKING_TRANSITIONS
            global_entente[
              galosh::Transition::fromPreAlign
            ].zero();
            global_entente[
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
#ifdef USE_DEL_IN_DEL_OUT
            backward_matrices_iterator++;
#endif // USE_DEL_IN_DEL_OUT
            prev_forward_matrices_iterator = forward_matrices_iterator;
            prev_forward_matrices_iterator--;
            score =
              dp.updateGlobalEntenteForSequence(
                dp_parameters,
                dna_profile,
                dna_seq_a,
                2,
                //forward_matrices[ 1 ][ 0 ],
                ( *prev_forward_matrices_iterator )[ 0 ],
                //forward_matrices[ 2 ][ 0 ],
                ( *forward_matrices_iterator )[ 0 ],
                //backward_matrices[ 2 ][ 0 ],
                ( *backward_matrices_iterator )[ 0 ],
                NULL,
                global_entente
              );
            cout << "Score for sequence a, calculated using updateGlobalEntenteForSequence(..) for position 2, is a: " << score << endl;
            cout << "BEFORE UNSCALING, Global Entente at Pos 2, for seq a, is " << global_entente << endl;
            global_entente.unscale();
            cout << "The Global Entente at Pos 2, for seq a, is " << global_entente << endl;
            // The Global Entente at Pos 2, for seq a, is [ M->(M->M=0.998308,M->I=0,M->D=0.000980175), I->(I->M=0,I->I=0), D->(D->M=0.000693053,D->D=1.85639e-05), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=1,E->J=0), C->(C->C=0.00127225,C->T=1), C_DNA:(A=6.18797e-06,C=0,G=0.00126606,T=0) ]

#ifdef USE_DEL_IN_DEL_OUT
            global_entente.zero();
#ifdef DISALLOW_FLANKING_TRANSITIONS
            global_entente[
              galosh::Transition::fromPreAlign
            ].zero();
            global_entente[
              galosh::Transition::fromPostAlign
            ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            forward_matrices_iterator = forward_matrices.begin();
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            forward_matrices_iterator++;
            backward_matrices_iterator = backward_matrices.rbegin();
            prev_forward_matrices_iterator = forward_matrices_iterator;
            prev_forward_matrices_iterator--;
            score =
              dp.updateGlobalEntenteForSequence(
                dp_parameters,
                dna_profile,
                dna_seq_a,
                3,
                //forward_matrices[ 2 ][ 0 ],
                ( *prev_forward_matrices_iterator )[ 0 ],
                //forward_matrices[ 3 ][ 0 ],
                ( *forward_matrices_iterator )[ 0 ],
                //backward_matrices[ 3 ][ 0 ],
                ( *backward_matrices_iterator )[ 0 ],
                NULL,
                global_entente
              );
            cout << "Score for sequence a, calculated using updateGlobalEntenteForSequence(..) for position 3, is a: " << score << endl;
            global_entente.unscale();
            cout << "The Global Entente at Pos 3, for seq a, is " << global_entente << endl;
#endif // USE_DEL_IN_DEL_OUT
            cout << endl;
      
            // TODO: REMOVE
            exit( 0 );

            double euclidean_dist_ab = ap_b.euclideanDistance( ap_a );
            double euclidean_dist_ac = ap_c.euclideanDistance( ap_a );
            double euclidean_dist_bc = ap_c.euclideanDistance( ap_b );
            cout << "Euclidean distances:" << endl;
            cout << "-\t" << euclidean_dist_ab << "\t" << euclidean_dist_ac << endl;
            cout << euclidean_dist_ab << "\t-\t" << euclidean_dist_bc << endl;
            cout << euclidean_dist_ac << "\t" << euclidean_dist_bc << "\t-" << endl;
            // Euclidean distances:
            // -	1.54792	1.57958
            // 1.54792	-	1.15749
            // 1.57958	1.15749	-
      
            galosh::MultinomialDistribution<Dna, double>::DistanceMatrix distance_matrix;
            //distance_matrix = 1.0;
            //cout << "Uniform distance matrix: " << endl << distance_matrix << endl;
            distance_matrix = 0.0;
            distance_matrix[ ( residue = 'a' ) ][ ( residue = 'a' ) ] = 1.0;
            cout << "Test distance matrix: " << endl << distance_matrix << endl;
            galosh::MultinomialDistribution<Dna, ProbabilityType> expected_distances;
            ap_c[ 1 ][ galosh::Emission::Match ].calculateExpectedDistances(
              distance_matrix,
              expected_distances
            );
            cout << "Expected distances (taken over the Alignment Profile Pos 1 for seq c) are" << endl;
            cout << expected_distances << endl;
            // (A=0.866944,C=0,G=0,T=0)
            ProbabilityType dist_ac_1 =
              ap_a[ 1 ][ galosh::Emission::Match ].calculateExpectedValue( expected_distances );
            cout << "dist_ac_1 = " << dist_ac_1 << endl;
            // dist_ca_1 = 0.587285
            ap_a[ 1 ][ galosh::Emission::Match ].calculateExpectedDistances(
              distance_matrix,
              expected_distances
            );
            cout << "Expected distances (taken over the Alignment Profile Pos 1 for seq a) are" << endl;
            cout << expected_distances << endl;
            // (A=0.677419,C=0,G=0,T=0)
            ProbabilityType dist_ca_1 =
              ap_c[ 1 ][ galosh::Emission::Match ].calculateExpectedValue( expected_distances );
            cout << "dist_ca_1 = " << dist_ca_1 << endl;
            // dist_ca_1 = 0.587285
      
            // Copied from muscle's nucmx.cpp and modified
            galosh::MultinomialDistribution<Dna, ScoreType>::DistanceMatrix
              blastz_score_matrix;
            blastz_score_matrix.isSymmetric = true;
            blastz_score_matrix.isProximityMatrix = true;
            blastz_score_matrix[ ( residue = 'a' ) ][ ( residue = 'a' ) ] = exp(91);
            blastz_score_matrix[ ( residue = 'a' ) ][ ( residue = 'c' ) ] = exp(-114);
            blastz_score_matrix[ ( residue = 'a' ) ][ ( residue = 'g' ) ] = exp(-31);
            blastz_score_matrix[ ( residue = 'a' ) ][ ( residue = 't' ) ] = exp(-123);
            blastz_score_matrix[ ( residue = 'c' ) ][ ( residue = 'a' ) ] = exp(-114);
            blastz_score_matrix[ ( residue = 'c' ) ][ ( residue = 'c' ) ] = exp(100);
            blastz_score_matrix[ ( residue = 'c' ) ][ ( residue = 'g' ) ] = exp(-125);
            blastz_score_matrix[ ( residue = 'c' ) ][ ( residue = 't' ) ] = exp(-31);
            blastz_score_matrix[ ( residue = 'g' ) ][ ( residue = 'a' ) ] = exp(-31);
            blastz_score_matrix[ ( residue = 'g' ) ][ ( residue = 'c' ) ] = exp(-125);
            blastz_score_matrix[ ( residue = 'g' ) ][ ( residue = 'g' ) ] = exp(100);
            blastz_score_matrix[ ( residue = 'g' ) ][ ( residue = 't' ) ] = exp(-114);
            blastz_score_matrix[ ( residue = 't' ) ][ ( residue = 'a' ) ] = exp(-123);
            blastz_score_matrix[ ( residue = 't' ) ][ ( residue = 'c' ) ] = exp(-31);
            blastz_score_matrix[ ( residue = 't' ) ][ ( residue = 'g' ) ] = exp(-114);
            blastz_score_matrix[ ( residue = 't' ) ][ ( residue = 't' ) ] = exp(91);
            cout << "BLASTZ DNA Score Matrix: " << endl << blastz_score_matrix << endl;
            // BLASTZ DNA Score Matrix: 
            // (A=(A=91,C=-114,G=-31,T=-123),C=(A=-114,C=100,G=-125,T=-31),G=(A=-31,C=-125,G=100,T=-114),T=(A=-123,C=-31,G=-114,T=91))
      
            galosh::MultinomialDistribution<Dna, ScoreType> expected_blastz_scores_c;
            ap_c[ 1 ][ galosh::Emission::Match ].calculateExpectedDistances(
              blastz_score_matrix,
              expected_blastz_scores_c
            );
            cout << "Expected proximities (taken over the Alignment Profile Pos 1 for seq c) are" << endl;
            cout << expected_blastz_scores_c << endl;
            // (A=90.8572,C=-114.143,G=97.9772,T=-116.022)
            ScoreType blastz_score_ac_1 =
              ap_a[ 1 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_scores_c );
            cout << "blastz_score_ac_1 = " << blastz_score_ac_1 << endl;
            // blastz_score_ac_1 = 90.4678
            ScoreType blastz_score_bc_1 =
              ap_b[ 1 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_scores_c );
            cout << "blastz_score_bc_1 = " << blastz_score_bc_1 << endl;
            // blastz_score_bc_1 = 90.9975
            galosh::MultinomialDistribution<Dna, ScoreType> expected_blastz_scores_b;
            ap_b[ 1 ][ galosh::Emission::Match ].calculateExpectedDistances(
              blastz_score_matrix,
              expected_blastz_scores_b
            );
            cout << "Expected proximities (taken over the Alignment Profile Pos 1 for seq b) are" << endl;
            cout << expected_blastz_scores_b << endl;
            // (A=90.9992,C=-114.001,G=90.9928,T=-122.311)
            ScoreType blastz_score_ab_1 =
              ap_a[ 1 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_scores_b );
            cout << "blastz_score_ab_1 = " << blastz_score_ab_1 << endl;
            // blastz_score_ab_1 = 90.6097
      
            galosh::MultinomialDistribution<Dna, ScoreType> expected_blastz_scores_a;
            ap_a[ 1 ][ galosh::Emission::Match ].calculateExpectedDistances(
              blastz_score_matrix,
              expected_blastz_scores_a
            );
            cout << "Expected proximities (taken over the Alignment Profile Pos 1 for seq a) are" << endl;
            cout << expected_blastz_scores_a << endl;
            // (A=90.6105,C=-114.389,G=-31.3895,T=-123.389)
            ScoreType blastz_score_ca_1 =
              ap_c[ 1 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_scores_a );
            cout << "blastz_score_ca_1 = " << blastz_score_ca_1 << endl;
            // blastz_score_ca_1 = 90.4678
            ScoreType blastz_score_ba_1 =
              ap_b[ 1 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_scores_a );
            cout << "blastz_score_ba_1 = " << blastz_score_ba_1 << endl;
            // blastz_score_ba_1 = 90.6097
      
            cout << "BLASTZ Score proximities (Pos 1):" << endl;
            cout << "-\t" << blastz_score_ab_1 << "\t" << blastz_score_ac_1 << endl;
            cout << blastz_score_ab_1 << "\t-\t" << blastz_score_bc_1 << endl;
            cout << blastz_score_ac_1 << "\t" << blastz_score_bc_1 << "\t-" << endl;
            // BLASTZ Score proximities (Pos 1):
            // -	90.6097	90.4678
            // 90.6097	-	90.9975
            // 90.4678	90.9975	-
      
            // Try again, using the blastz scores as-is (no exponentiating).
            galosh::MultinomialDistribution<Dna, double>::DistanceMatrix
              blastz_double_matrix;
            blastz_double_matrix.isSymmetric = true;
            blastz_double_matrix.isProximityMatrix = true;
            blastz_double_matrix[ ( residue = 'a' ) ][ ( residue = 'a' ) ] = 91;
            blastz_double_matrix[ ( residue = 'a' ) ][ ( residue = 'c' ) ] = -114;
            blastz_double_matrix[ ( residue = 'a' ) ][ ( residue = 'g' ) ] = -31;
            blastz_double_matrix[ ( residue = 'a' ) ][ ( residue = 't' ) ] = -123;
            blastz_double_matrix[ ( residue = 'c' ) ][ ( residue = 'a' ) ] = -114;
            blastz_double_matrix[ ( residue = 'c' ) ][ ( residue = 'c' ) ] = 100;
            blastz_double_matrix[ ( residue = 'c' ) ][ ( residue = 'g' ) ] = -125;
            blastz_double_matrix[ ( residue = 'c' ) ][ ( residue = 't' ) ] = -31;
            blastz_double_matrix[ ( residue = 'g' ) ][ ( residue = 'a' ) ] = -31;
            blastz_double_matrix[ ( residue = 'g' ) ][ ( residue = 'c' ) ] = -125;
            blastz_double_matrix[ ( residue = 'g' ) ][ ( residue = 'g' ) ] = 100;
            blastz_double_matrix[ ( residue = 'g' ) ][ ( residue = 't' ) ] = -114;
            blastz_double_matrix[ ( residue = 't' ) ][ ( residue = 'a' ) ] = -123;
            blastz_double_matrix[ ( residue = 't' ) ][ ( residue = 'c' ) ] = -31;
            blastz_double_matrix[ ( residue = 't' ) ][ ( residue = 'g' ) ] = -114;
            blastz_double_matrix[ ( residue = 't' ) ][ ( residue = 't' ) ] = 91;
            cout << "BLASTZ DNA double Matrix: " << endl << blastz_double_matrix << endl;
            // (A=(A=91,C=-114,G=-31,T=-123),C=(A=-114,C=100,G=-125,T=-31),G=(A=-31,C=-125,G=100,T=-114),T=(A=-123,C=-31,G=-114,T=91))
      
            galosh::MultinomialDistribution<Dna, double> expected_blastz_doubles_c;
            ap_c[ 0 ][ galosh::Emission::Match ].calculateExpectedDistances(
              blastz_double_matrix,
              expected_blastz_doubles_c
            );
            cout << "Expected proximities (taken over the Alignment Profile Pos 0 for seq c) are" << endl;
            cout << expected_blastz_doubles_c << endl;
            // (A=0,C=0,G=0,T=0)
            double blastz_double_ac_0 =
              ap_a[ 0 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_doubles_c );
            cout << "blastz_double_ac_0 = " << blastz_double_ac_0 << endl;
            // blastz_double_ac_0 = 0
            double blastz_double_bc_0 =
              ap_b[ 0 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_doubles_c );
            cout << "blastz_double_bc_0 = " << blastz_double_bc_0 << endl;
            // blastz_double_bc_0 = 1.89205e-301
            galosh::MultinomialDistribution<Dna, double> expected_blastz_doubles_b;
            ap_b[ 0 ][ galosh::Emission::Match ].calculateExpectedDistances(
              blastz_double_matrix,
              expected_blastz_doubles_b
            );
            cout << "Expected proximities (taken over the Alignment Profile Pos 0 for seq b) are" << endl;
            cout << expected_blastz_doubles_b << endl;
            // (A=0,C=0,G=0,T=0)
            double blastz_double_ab_0 =
              ap_a[ 0 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_doubles_b );
            cout << "blastz_double_ab_0 = " << blastz_double_ab_0 << endl;
            // blastz_double_ab_0 = 1.89934e-301
      
            galosh::MultinomialDistribution<Dna, double> expected_blastz_doubles_a;
            ap_a[ 0 ][ galosh::Emission::Match ].calculateExpectedDistances(
              blastz_double_matrix,
              expected_blastz_doubles_a
            );
            cout << "Expected proximities (taken over the Alignment Profile Pos 0 for seq a) are" << endl;
            cout << expected_blastz_doubles_a << endl;
            // (A=0,C=0,G=0,T=0)
            double blastz_double_ca_0 =
              ap_c[ 0 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_doubles_a );
            cout << "blastz_double_ca_0 = " << blastz_double_ca_0 << endl;
            // blastz_double_ca_0 = 1.88476e-301
            double blastz_double_ba_0 =
              ap_b[ 0 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_doubles_a );
            cout << "blastz_double_ba_0 = " << blastz_double_ba_0 << endl;
            // blastz_double_ba_0 = 1.89205e-301
            cout << "BLASTZ double proximities (Pos 0):" << endl;
            cout << "-\t" << blastz_double_ab_0 << "\t" << blastz_double_ac_0 << endl;
            cout << blastz_double_ab_0 << "\t-\t" << blastz_double_bc_0 << endl;
            cout << blastz_double_ac_0 << "\t" << blastz_double_bc_0 << "\t-" << endl;
            // BLASTZ double proximities (Pos 0):
            // -	1.89934e-301	0
            // 1.89934e-301	-	1.89205e-301
            // 0	1.89205e-301	-
      
            ap_c[ 1 ][ galosh::Emission::Match ].calculateExpectedDistances(
              blastz_double_matrix,
              expected_blastz_doubles_c
            );
            cout << "Expected proximities (taken over the Alignment Profile Pos 1 for seq c) are" << endl;
            cout << expected_blastz_doubles_c << endl;
            // (A=74.7912,C=-115.367,G=-13.6471,T=-121.714)
            double blastz_double_ac_1 =
              ap_a[ 1 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_doubles_c );
            cout << "blastz_double_ac_1 = " << blastz_double_ac_1 << endl;
            // blastz_double_ac_1 = 50.665
            double blastz_double_bc_1 =
              ap_b[ 1 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_doubles_c );
            cout << "blastz_double_bc_1 = " << blastz_double_bc_1 << endl;
            // blastz_double_bc_1 = 74.7271
            ap_b[ 1 ][ galosh::Emission::Match ].calculateExpectedDistances(
              blastz_double_matrix,
              expected_blastz_doubles_b
            );
            cout << "Expected proximities (taken over the Alignment Profile Pos 1 for seq b) are" << endl;
            cout << expected_blastz_doubles_b << endl;
            // (A=90.9203,C=-113.92,G=-30.9619,T=-122.911)
            double blastz_double_ab_1 =
              ap_a[ 1 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_doubles_b );
            cout << "blastz_double_ab_1 = " << blastz_double_ab_1 << endl;
            // blastz_double_ab_1 = 61.5912
      
            ap_a[ 1 ][ galosh::Emission::Match ].calculateExpectedDistances(
              blastz_double_matrix,
              expected_blastz_doubles_a
            );
            cout << "Expected proximities (taken over the Alignment Profile Pos 1 for seq a) are" << endl;
            cout << expected_blastz_doubles_a << endl;
            // (A=61.6452,C=-77.2258,G=-21,T=-83.3226)
            double blastz_double_ca_1 =
              ap_c[ 1 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_doubles_a );
            cout << "blastz_double_ca_1 = " << blastz_double_ca_1 << endl;
            // blastz_double_ca_1 = 50.665
            double blastz_double_ba_1 =
              ap_b[ 1 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_doubles_a );
            cout << "blastz_double_ba_1 = " << blastz_double_ba_1 << endl;
            // blastz_double_ba_1 = 61.5912
      
            cout << "BLASTZ double proximities (Pos 1):" << endl;
            cout << "-\t" << blastz_double_ab_1 << "\t" << blastz_double_ac_1 << endl;
            cout << blastz_double_ab_1 << "\t-\t" << blastz_double_bc_1 << endl;
            cout << blastz_double_ac_1 << "\t" << blastz_double_bc_1 << "\t-" << endl;
            // BLASTZ double proximities (Pos 1):
            // -	61.5912	50.665
            // 61.5912	-	74.7271
            // 50.665	74.7271	-
      
            ap_c[ 2 ][ galosh::Emission::Match ].calculateExpectedDistances(
              blastz_double_matrix,
              expected_blastz_doubles_c
            );
            cout << "Expected proximities (taken over the Alignment Profile Pos 2 for seq c) are" << endl;
            cout << expected_blastz_doubles_c << endl;
            // (A=6.18506,C=-121.489,G=59.9854,T=-116.593)
            double blastz_double_ac_2 =
              ap_a[ 2 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_doubles_c );
            cout << "blastz_double_ac_2 = " << blastz_double_ac_2 << endl;
            // blastz_double_ac_2 = 1.93472
            double blastz_double_bc_2 =
              ap_b[ 2 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_doubles_c );
            cout << "blastz_double_bc_2 = " << blastz_double_bc_2 << endl;
            // blastz_double_bc_2 = 59.9042
            ap_b[ 2 ][ galosh::Emission::Match ].calculateExpectedDistances(
              blastz_double_matrix,
              expected_blastz_doubles_b
            );
            cout << "Expected proximities (taken over the Alignment Profile Pos 2 for seq b) are" << endl;
            cout << expected_blastz_doubles_b << endl;
            // (A=-30.9207,C=-124.871,G=99.8482,T=-113.89)
            double blastz_double_ab_2 =
              ap_a[ 2 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_doubles_b );
            cout << "blastz_double_ab_2 = " << blastz_double_ab_2 << endl;
            // blastz_double_ab_2 = -9.67217
      
            ap_a[ 2 ][ galosh::Emission::Match ].calculateExpectedDistances(
              blastz_double_matrix,
              expected_blastz_doubles_a
            );
            cout << "Expected proximities (taken over the Alignment Profile Pos 2 for seq a) are" << endl;
            cout << expected_blastz_doubles_a << endl;
            // (A=28.4653,C=-35.6598,G=-9.69697,T=-38.4751)
            double blastz_double_ca_2 =
              ap_c[ 2 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_doubles_a );
            cout << "blastz_double_ca_2 = " << blastz_double_ca_2 << endl;
            // blastz_double_ca_2 = 1.93472
            double blastz_double_ba_2 =
              ap_b[ 2 ][ galosh::Emission::Match ].calculateExpectedValue( expected_blastz_doubles_a );
            cout << "blastz_double_ba_2 = " << blastz_double_ba_2 << endl;
            // blastz_double_ba_2 = -9.67217
      
            cout << "BLASTZ double proximities (Pos 2):" << endl;
            cout << "-\t" << blastz_double_ab_2 << "\t" << blastz_double_ac_2 << endl;
            cout << blastz_double_ab_2 << "\t-\t" << blastz_double_bc_2 << endl;
            cout << blastz_double_ac_2 << "\t" << blastz_double_bc_2 << "\t-" << endl;
            // BLASTZ double proximities (Pos 2):
            // -	-9.67217	1.93472
            // -9.67217	-	59.9042
            // 1.93472	59.9042	-
      
            // Summary across the positions
            double blastz_double_ab =
              ( blastz_double_ab_0 + blastz_double_ab_1 + blastz_double_ab_2 );
            double blastz_double_ac =
              ( blastz_double_ac_0 + blastz_double_ac_1 + blastz_double_ac_2 );
            double blastz_double_bc =
              ( blastz_double_bc_0 + blastz_double_bc_1 + blastz_double_bc_2 );
            cout << "BLASTZ double proximities:" << endl;
            cout << "-\t" << blastz_double_ab << "\t" << blastz_double_ac << endl;
            cout << blastz_double_ab << "\t-\t" << blastz_double_bc << endl;
            cout << blastz_double_ac << "\t" << blastz_double_bc << "\t-" << endl;
            // BLASTZ double proximities:
            // -	51.919	52.5997
            // 51.919	-	134.631
            // 52.5997	134.631	-
      
            // Probability-of-co-occurence metric
      
            // Pos 0 has 0 prob of co-occurrence, since there's no Matches here.
            //MatrixValueType cooc_ab_0 =
            //  ap_a[ 0 ][ galosh::Emission::Match ].calculateExpectedValue(
            //    ap_b[ 0 ][ galosh::Emission::Match ]
            //  );
            ////cout << "Probability of co-occurrence for Pos 0: " << endl;
            //cout << "cooc_ab_0 is " << cooc_ab_0 << endl;
            //// Should be the same as:
            ////MatrixValueType cooc_ba_0 =
            ////  ap_b[ 0 ][ galosh::Emission::Match ].calculateExpectedValue(
            ////    ap_a[ 0 ][ galosh::Emission::Match ]
            ////  );
            ////cout << "cooc_ba_0 is " << cooc_ba_0 << endl;
            //MatrixValueType cooc_ac_0 =
            //  ap_a[ 0 ][ galosh::Emission::Match ].calculateExpectedValue(
            //    ap_c[ 0 ][ galosh::Emission::Match ]
            //  );
            //cout << "cooc_ac_0 is " << cooc_ac_0 << endl;
            //MatrixValueType cooc_bc_0 =
            //  ap_b[ 0 ][ galosh::Emission::Match ].calculateExpectedValue(
            //    ap_c[ 0 ][ galosh::Emission::Match ]
            //  );
            //cout << "cooc_bc_0 is " << cooc_bc_0 << endl;
            //
            //cout << "Probability of co-occurrence (Pos 0):" << endl;
            //cout << "-\t" << cooc_ab_0 << "\t" << cooc_ac_0 << endl;
            //cout << cooc_ab_0 << "\t-\t" << cooc_bc_0 << endl;
            //cout << cooc_ac_0 << "\t" << cooc_bc_0 << "\t-" << endl;
      
            MatrixValueType cooc_ab_1 =
              ap_a[ 1 ][ galosh::Emission::Match ].calculateExpectedValue(
                ap_b[ 1 ][ galosh::Emission::Match ]
              );
            //cout << "Probability of co-occurrence for Pos 1: " << endl;
            cout << "cooc_ab_1 is " << cooc_ab_1 << endl;
            // Should be the same as:
            //MatrixValueType cooc_ba_1 =
            //  ap_b[ 1 ][ galosh::Emission::Match ].calculateExpectedValue(
            //    ap_a[ 1 ][ galosh::Emission::Match ]
            //  );
            //cout << "cooc_ba_1 is " << cooc_ba_1 << endl;
            MatrixValueType cooc_ac_1 =
              ap_a[ 1 ][ galosh::Emission::Match ].calculateExpectedValue(
                ap_c[ 1 ][ galosh::Emission::Match ]
              );
            cout << "cooc_ac_1 is " << cooc_ac_1 << endl;
            MatrixValueType cooc_bc_1 =
              ap_b[ 1 ][ galosh::Emission::Match ].calculateExpectedValue(
                ap_c[ 1 ][ galosh::Emission::Match ]
              );
            cout << "cooc_bc_1 is " << cooc_bc_1 << endl;
      
            cout << "Probability of co-occurrence (Pos 1):" << endl;
            cout << "-\t" << cooc_ab_1 << "\t" << cooc_ac_1 << endl;
            cout << cooc_ab_1 << "\t-\t" << cooc_bc_1 << endl;
            cout << cooc_ac_1 << "\t" << cooc_bc_1 << "\t-" << endl;
            // Probability of co-occurrence (Pos 1):
            // -	0.676854	0.587285
            // 0.676854	-	0.866237
            // 0.587285	0.866237	-
      
            MatrixValueType cooc_ab_2 =
              ap_a[ 2 ][ galosh::Emission::Match ].calculateExpectedValue(
                ap_b[ 2 ][ galosh::Emission::Match ]
              );
            //cout << "Probability of co-occurrence for Pos 2: " << endl;
            cout << "cooc_ab_2 is " << cooc_ab_2 << endl;
            // Should be the same as:
            //MatrixValueType cooc_ba_2 =
            //  ap_b[ 2 ][ galosh::Emission::Match ].calculateExpectedValue(
            //    ap_a[ 2 ][ galosh::Emission::Match ]
            //  );
            //cout << "cooc_ba_2 is " << cooc_ba_2 << endl;
            MatrixValueType cooc_ac_2 =
              ap_a[ 2 ][ galosh::Emission::Match ].calculateExpectedValue(
                ap_c[ 2 ][ galosh::Emission::Match ]
              );
            cout << "cooc_ac_2 is " << cooc_ac_2 << endl;
            MatrixValueType cooc_bc_2 =
              ap_b[ 2 ][ galosh::Emission::Match ].calculateExpectedValue(
                ap_c[ 2 ][ galosh::Emission::Match ]
              );
            cout << "cooc_bc_2 is " << cooc_bc_2 << endl;
      
            cout << "Probability of co-occurrence (Pos 2):" << endl;
            cout << "-\t" << cooc_ab_2 << "\t" << cooc_ac_2 << endl;
            cout << cooc_ab_2 << "\t-\t" << cooc_bc_2 << endl;
            cout << cooc_ac_2 << "\t" << cooc_bc_2 << "\t-" << endl;
            // Probability of co-occurrence (Pos 2):
            // -	0.00012388	0.0952388
            // 0.00012388	-	0.693391
            // 0.0952388	0.693391	-
      
            // Summary across the positions (note we exclude pos 0, since no Matches)
            MatrixValueType cooc_ab =
              ( cooc_ab_1 * cooc_ab_2 );
            MatrixValueType cooc_ac =
              ( cooc_ac_1 * cooc_ac_2 );
            MatrixValueType cooc_bc =
              ( cooc_bc_1 * cooc_bc_2 );
            cout << "Probability of co-occurrence:" << endl;
            cout << "-\t" << cooc_ab << "\t" << cooc_ac << endl;
            cout << cooc_ab << "\t-\t" << cooc_bc << endl;
            cout << cooc_ac << "\t" << cooc_bc << "\t-" << endl;
            // Probability of co-occurrence:
            // -	8.3849e-05	0.0559323
            // 8.3849e-05	-	0.600641
            // 0.0559323	0.600641	-
      
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::DistanceMatrix<MatrixValueType> cooc_matrix(
              3,
              true,
              true
            );
            cooc_matrix( 0, 1 ) = cooc_ab;
            cooc_matrix( 0, 2 ) = cooc_ac;
            cooc_matrix( 1, 2 ) = cooc_bc;
            cout << "Probability of co-occurrence as distance matrix:" << endl;
            cout << cooc_matrix << endl;
            // [ 0, 8.3849e-05, 0.0559323 ]
            // [ 8.3849e-05, 0, 0.600641 ]
            // [ 0.0559323, 0.600641, 0 ]
      
            // Ok now try the newfangled way
            vector<galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile> aps( 3 );
            aps[ 0 ].reinitialize( ( dna_profile.length() + 1 ) );
            aps[ 1 ].reinitialize( ( dna_profile.length() + 1 ) );
            aps[ 2 ].reinitialize( ( dna_profile.length() + 1 ) );
            dp.calculateAlignmentProfiles(
              dp_parameters,
              dna_profile,
              dna_seqs,
              num_sequences_to_use,
              forward_matrices,
              aps
            );
            cout << "After running calculateAlignmentProfiles(..), got these aps:" << endl;
            cout << "aps[ 0 ]:" << endl;
            cout << aps[ 0 ] << endl;
            cout << "aps[ 1 ]:" << endl;
            cout << aps[ 1 ] << endl;
            cout << "aps[ 2 ]:" << endl;
            cout << aps[ 2 ] << endl;
            // aps[ 0 ]:
            // [ M_DNA:(A=0,C=0,G=0,T=0), M->(M->M=0,M->I=0,M->D=0), I->(I->M=0,I->I=0), D->(D->M=0,D->D=0), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0.00488759,N->B=1), B->(B->M=0.677419,B->D=0.322581), N_DNA:(A=0.00488759,C=0,G=0,T=0), E->(E->C=0,E->J=0), C->(C->C=0,C->T=0), C_DNA:(A=0,C=0,G=0,T=0) ]
            // [ M_DNA:(A=0.677419,C=0,G=0,T=0), M->(M->M=0,M->I=0,M->D=0.677419), I->(I->M=0,I->I=0), D->(D->M=0.312805,D->D=0.00977517), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=0,E->J=0), C->(C->C=0,C->T=0), C_DNA:(A=0,C=0,G=0,T=0) ]
            // [ M_DNA:(A=0.312805,C=0,G=0,T=0), M->(M->M=0,M->I=0,M->D=0), I->(I->M=0,I->I=0), D->(D->M=0,D->D=0), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=1,E->J=0), C->(C->C=0.00488759,C->T=1), C_DNA:(A=0.00488759,C=0,G=0,T=0) ]
            // 
            // aps[ 1 ]:
            // [ M_DNA:(A=0,C=0,G=0,T=0), M->(M->M=0,M->I=0,M->D=0), I->(I->M=0,I->I=0), D->(D->M=0,D->D=0), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0.000438109,N->B=1), B->(B->M=0.999288,B->D=0.000711617), N_DNA:(A=0.000431921,C=0,G=6.18797e-06,T=0), E->(E->C=0,E->J=0), C->(C->C=0,C->T=0), C_DNA:(A=0,C=0,G=0,T=0) ]
            // [ M_DNA:(A=0.999166,C=0,G=0.000122522,T=0), M->(M->M=0.998308,M->I=0,M->D=0.000980175), I->(I->M=0,I->I=0), D->(D->M=0.000693053,D->D=1.85639e-05), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=0,E->J=0), C->(C->C=0,C->T=0), C_DNA:(A=0,C=0,G=0,T=0) ]
            // [ M_DNA:(A=0.00039603,C=0,G=0.998605,T=0), M->(M->M=0,M->I=0,M->D=0), I->(I->M=0,I->I=0), D->(D->M=0,D->D=0), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=1,E->J=0), C->(C->C=0.00127225,C->T=1), C_DNA:(A=6.18797e-06,C=0,G=0.00126606,T=0) ]
            // 
            // aps[ 2 ]:
            // [ M_DNA:(A=0,C=0,G=0,T=0), M->(M->M=0,M->I=0,M->D=0), I->(I->M=0,I->I=0), D->(D->M=0,D->D=0), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0.134257,N->B=1), B->(B->M=0.999226,B->D=0.000774347), N_DNA:(A=0.133377,C=0,G=0.000880175,T=0), E->(E->C=0,E->J=0), C->(C->C=0,C->T=0), C_DNA:(A=0,C=0,G=0,T=0) ]
            // [ M_DNA:(A=0.866944,C=0,G=0.132282,T=0), M->(M->M=0.826229,M->I=0.171719,M->D=0.00127767), I->(I->M=0.171719,I->I=0), D->(D->M=0.00075714,D->D=1.72077e-05), I_DNA:(A=0,C=0,G=0.171719,T=0), N->(N->N=0,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=0,E->J=0), C->(C->C=0,C->T=0), C_DNA:(A=0,C=0,G=0,T=0) ]
            // [ M_DNA:(A=0.304467,C=0,G=0.694239,T=0), M->(M->M=0,M->I=0,M->D=0), I->(I->M=0,I->I=0), D->(D->M=0,D->D=0), I_DNA:(A=0,C=0,G=0,T=0), N->(N->N=0,N->B=0), B->(B->M=0,B->D=0), N_DNA:(A=0,C=0,G=0,T=0), E->(E->C=1,E->J=0), C->(C->C=0.696093,C->T=1), C_DNA:(A=0.695213,C=0,G=0.000880175,T=0) ]
      
            // MatchEmissionCoocucurrenceProbability
            MatrixValueType cooc_01 =
              aps[ 0 ].calculateMatchEmissionCooccurrenceProbability( aps[ 1 ] );
            cout << "cooc_01 is " << cooc_01 << endl;
            cout << "should be the same as cooc_ab == " << cooc_ab << endl;
            // cooc_01 is 8.3849e-05
            // should be the same as cooc_ab == 8.3849e-05
      
            cooc_matrix = ( MatrixValueType )0;
            cooc_matrix.setToMatchEmissionCooccurrenceProbabilities( aps );
            cout << "After setting using setToMatchEmissionCooccurrenceProbabilities( aps ), cooc_matrix is " << endl;
            cout << cooc_matrix << endl;
            // [ 0.0449018, 8.3849e-05, 0.0559323 ]
            // [ 8.3849e-05, 0.99555, 0.600641 ]
            // [ 0.0559323, 0.600641, 0.441971 ]
      
#ifdef __HAVE_MUSCLE
            cout << "Creating muscle dist func." << endl;
            DistFunc muscle_dist_func;
            muscle_dist_func.SetCount( 3 );
            muscle_dist_func.SetName( 0, "a" );
            muscle_dist_func.SetName( 1, "b" );
            muscle_dist_func.SetName( 2, "c" );
            muscle_dist_func.SetId( 0, 0 );
            muscle_dist_func.SetId( 1, 1 );
            muscle_dist_func.SetId( 2, 2 );
            muscle_dist_func.SetDist( 0, 1, ( float )( 1.0 / galosh::toDouble( cooc_ab ) ) );
            muscle_dist_func.SetDist( 0, 2, ( float )( 1.0 / galosh::toDouble( cooc_ac ) ) );
            muscle_dist_func.SetDist( 1, 2, ( float )( 1.0 / galosh::toDouble( cooc_bc ) ) );
            cout << "done." << endl;
      
            cout << "Creating muscle_clust." << endl;
            ClustSetDF muscle_clust_set( muscle_dist_func );
            Clust muscle_clust;
            muscle_clust.Create( muscle_clust_set, CLUSTER_UPGMB );
            cout << "done." << endl;
      
            cout << "Creating muscle_tree." << endl;
            Tree muscle_tree;
            muscle_tree.FromClust( muscle_clust );
            cout << "done." << endl;
      
            cout << "Writing tree to file test_muscle_tree.dnd" << endl;
            TextFile f( "test_muscle_tree.dnd", true );
            muscle_tree.ToFile( f );
            cout << "done." << endl;
#endif // __HAVE_MUSCLE
      
            // Try cross-entropy.
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::DistanceMatrix<double> entropy_matrix;
            entropy_matrix.setToCrossEntropies( aps );
            cout << "After setting using setToCrossEntropies( aps ), entropy_matrix is " << endl;
            cout << entropy_matrix << endl;
            // [ 0.0449018, 8.3849e-05, 0.0559323 ]
            // [ 8.3849e-05, 0.99555, 0.600641 ]
            // [ 0.0559323, 0.600641, 0.441971 ]
      
#ifdef __HAVE_MUSCLE
            cout << "Creating muscle dist func." << endl;
            //DistFunc muscle_dist_func;
            //muscle_dist_func.SetCount( 3 );
            //muscle_dist_func.SetName( 0, "a" );
            //muscle_dist_func.SetName( 1, "b" );
            //muscle_dist_func.SetName( 2, "c" );
            //muscle_dist_func.SetId( 0, 0 );
            //muscle_dist_func.SetId( 1, 1 );
            //muscle_dist_func.SetId( 2, 2 );
            for( uint32_t from_i; from_i < 3; from_i++ ) {
              for( uint32_t to_i; to_i < 3; to_i++ ) {
                muscle_dist_func.SetDist( from_i, to_i, ( float )entropy_matrix( from_i, to_i ) );
              }
            }
            cout << "done." << endl;
      
            cout << "Creating muscle_clust_entropy." << endl;
            ClustSetDF muscle_clust_set_entropy( muscle_dist_func );
            Clust muscle_clust_entropy;
            muscle_clust_entropy.Create( muscle_clust_set_entropy, CLUSTER_UPGMB );
            cout << "done." << endl;
      
            cout << "Creating muscle_tree_entropy." << endl;
            Tree muscle_tree_entropy;
            muscle_tree_entropy.FromClust( muscle_clust_entropy );
            cout << "done." << endl;
      
            cout << "Writing tree to file test_muscle_tree_entropy.dnd" << endl;
            TextFile f_entropy( "test_muscle_tree_entropy.dnd", true );
            muscle_tree_entropy.ToFile( f_entropy );
            cout << "done." << endl;
#endif // __HAVE_MUSCLE
      
            // Try KL Distance
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::DistanceMatrix<double> kl_matrix;
            kl_matrix.setToKullbackLeiblerDivergences( aps );
            cout << "After setting using setToKullbackLeiblerDivergences( aps ), kl_matrix is " << endl;
            cout << kl_matrix << endl;
            // [ 0.0449018, 8.3849e-05, 0.0559323 ]
            // [ 8.3849e-05, 0.99555, 0.600641 ]
            // [ 0.0559323, 0.600641, 0.441971 ]
      
#ifdef __HAVE_MUSCLE
            cout << "Creating muscle dist func." << endl;
            //DistFunc muscle_dist_func;
            //muscle_dist_func.SetCount( 3 );
            //muscle_dist_func.SetName( 0, "a" );
            //muscle_dist_func.SetName( 1, "b" );
            //muscle_dist_func.SetName( 2, "c" );
            //muscle_dist_func.SetId( 0, 0 );
            //muscle_dist_func.SetId( 1, 1 );
            //muscle_dist_func.SetId( 2, 2 );
            for( uint32_t from_i; from_i < 3; from_i++ ) {
              for( uint32_t to_i; to_i < 3; to_i++ ) {
                muscle_dist_func.SetDist( from_i, to_i, ( float )kl_matrix( from_i, to_i ) );
              }
            }
            cout << "done." << endl;
      
            cout << "Creating muscle_clust_kl." << endl;
            ClustSetDF muscle_clust_set_kl( muscle_dist_func );
            Clust muscle_clust_kl;
            muscle_clust_kl.Create( muscle_clust_set_kl, CLUSTER_UPGMB );
            cout << "done." << endl;
      
            cout << "Creating muscle_tree_kl." << endl;
            Tree muscle_tree_kl;
            muscle_tree_kl.FromClust( muscle_clust_kl );
            cout << "done." << endl;
      
            cout << "Writing tree to file test_muscle_tree_kl.dnd" << endl;
            TextFile f_kl( "test_muscle_tree_kl.dnd", true );
            muscle_tree_kl.ToFile( f_kl );
            cout << "done." << endl;
#endif // __HAVE_MUSCLE
      
            // mark
            // Try SKL.
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::DistanceMatrix<double> skl_matrix;
            skl_matrix.setToSymmeterizedKullbackLeiblerDivergences( aps );
            cout << "After setting using setToSymmeterizedKullbackLeiblerDivergences( aps ), skl_matrix is " << endl;
            cout << skl_matrix << endl;
            // [ 0.0449018, 8.3849e-05, 0.0559323 ]
            // [ 8.3849e-05, 0.99555, 0.600641 ]
            // [ 0.0559323, 0.600641, 0.441971 ]
      
#ifdef __HAVE_MUSCLE
            cout << "Creating muscle dist func." << endl;
            //DistFunc muscle_dist_func;
            //muscle_dist_func.SetCount( 3 );
            //muscle_dist_func.SetName( 0, "a" );
            //muscle_dist_func.SetName( 1, "b" );
            //muscle_dist_func.SetName( 2, "c" );
            //muscle_dist_func.SetId( 0, 0 );
            //muscle_dist_func.SetId( 1, 1 );
            //muscle_dist_func.SetId( 2, 2 );
            for( uint32_t from_i; from_i < 3; from_i++ ) {
              for( uint32_t to_i; to_i < 3; to_i++ ) {
                muscle_dist_func.SetDist( from_i, to_i, ( float )skl_matrix( from_i, to_i ) );
              }
            }
            cout << "done." << endl;
      
            cout << "Creating muscle_clust_skl." << endl;
            ClustSetDF muscle_clust_set_skl( muscle_dist_func );
            Clust muscle_clust_skl;
            muscle_clust_skl.Create( muscle_clust_set_skl, CLUSTER_UPGMB );
            cout << "done." << endl;
      
            cout << "Creating muscle_tree_skl." << endl;
            Tree muscle_tree_skl;
            muscle_tree_skl.FromClust( muscle_clust_skl );
            cout << "done." << endl;
      
            cout << "Writing tree to file test_muscle_tree_skl.dnd" << endl;
            TextFile f_skl( "test_muscle_tree_skl.dnd", true );
            muscle_tree_skl.ToFile( f_skl );
            cout << "done." << endl;
#endif // __HAVE_MUSCLE
            // endmark
      
            // TODO: REMOVE
            exit( 0 );
          } // End if test_alignment_profiles
      
          // Test using priors
          if( test_priors ) {
            galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::DirichletMixtureMatchEmissionPrior<float> emission_prior( 1 );
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 'a' ) ] = 1.0;
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 'c' ) ] = 1.0;
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 'g' ) ] = 1.0;
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 't' ) ] = 1.0;
  
            position_entente = position_entente_backup;
            cout << "The position entente, before incorporating any priors, is " << position_entente << endl;
            // The position entente, before incorporating any priors, is [ M_DNA:(A=1,C=0,G=0.0520553,T=0), scalar:-0.933553 ]
            emission_prior.incorporatePrior( position_entente );
            cout << "The position entente, after incorporating a simple Laplace prior, is " << position_entente << endl;
            // The position entente, after incorporating a simple Laplace prior, is [ M_DNA:(A=1.39315,C=0.393155,G=0.44521,T=0.393155), scalar:-0.933553 ]
            position_entente_unscaled = position_entente;
            position_entente_unscaled.unscale();
            cout << "Unscaled, the position entente is " << position_entente_unscaled << endl;
            // Unscaled, the position entente is [ M_DNA:(A=3.54353,C=1,G=1.1324,T=1) ]
            position_entente.normalize( 0 );
            cout << "The position entente, normalized (min=0), is " << position_entente << endl;
            //The position entente, normalized (min=0), is [ M_dna:(A=0.530792,C=0.149792,G=0.169625,T=0.149792) ]
            // Revert
            position_entente = position_entente_backup;
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 'a' ) ] = 10.0;
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 'c' ) ] = 10.0;
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 'g' ) ] = 10.0;
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 't' ) ] = 10.0;
            emission_prior.incorporatePrior( position_entente );
            cout << "The position entente, after incorporating a stronger prior (10.0 on each residue), is " << position_entente << endl;
            // The position entente, after incorporating a stronger prior (10.0 on each residue), is [ M_DNA:(A=4.93155,C=3.93155,G=3.9836,T=3.93155), scalar:-0.933553 ]
            position_entente_unscaled = position_entente;
            position_entente_unscaled.unscale();
            cout << "Unscaled, the position entente is " << position_entente_unscaled << endl;
            // Unscaled, the position entente is [ M_DNA:(A=12.5435,C=10,G=10.1324,T=10) ]
            position_entente.normalize( 0 );
            cout << "The position entente, normalized (min=0), is " << position_entente << endl;
            //The position entente, normalized (min=0), is [ M_dna:(A=0.293925,C=0.234324,G=0.237427,T=0.234324) ]
            // Revert
            position_entente = position_entente_backup;
            emission_prior.reinitialize( 2 );
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 'a' ) ] = 1.0;
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 'c' ) ] = 1.0;
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 'g' ) ] = 1.0;
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 't' ) ] = 1.0;
            emission_prior[ 1 ][ galosh::Emission::Match ][ ( residue = 'a' ) ] = 1.0;
            emission_prior[ 1 ][ galosh::Emission::Match ][ ( residue = 'c' ) ] = 1.0;
            emission_prior[ 1 ][ galosh::Emission::Match ][ ( residue = 'g' ) ] = 1.0;
            emission_prior[ 1 ][ galosh::Emission::Match ][ ( residue = 't' ) ] = 1.0;
            emission_prior.incorporatePrior( position_entente );
            cout << "The position entente, after incorporating an even mixture of the Laplace prior with itself, is " << position_entente << endl;
            // The position entente, after incorporating an even mixture of the Laplace prior with itself, is [ M_DNA:(A=0.208683,C=0.0588913,G=0.0666888,T=0.0588913), scalar:-0.933553 ]
            position_entente_unscaled = position_entente;
            position_entente_unscaled.unscale();
            cout << "Unscaled, the position entente is " << position_entente_unscaled << endl;
            // Unscaled, the position entente is [ M_DNA:(A=0.530792,C=0.149792,G=0.169625,T=0.149792) ]
            position_entente.normalize( 0 );
            cout << "The position entente, normalized (min=0), is " << position_entente << endl;
            // The position entente, normalized (min=0), is [ M_DNA:(A=0.530792,C=0.149792,G=0.169625,T=0.149792) ]
            // Revert
            position_entente = position_entente_backup;
          
            emission_prior[ 1 ][ galosh::Emission::Match ][ ( residue = 'a' ) ] = 10.0;
            emission_prior[ 1 ][ galosh::Emission::Match ][ ( residue = 'c' ) ] = 10.0;
            emission_prior[ 1 ][ galosh::Emission::Match ][ ( residue = 'g' ) ] = 10.0;
            emission_prior[ 1 ][ galosh::Emission::Match ][ ( residue = 't' ) ] = 10.0;
            emission_prior.incorporatePrior( position_entente );
            cout << "The position entente, after incorporating an even mixture of the first 2 priors, is " << position_entente << endl;
            //The position entente, after incorporating an even mixture of the first 2 priors, is [ M_dna:(A=0.175345,C=0.070789,G=0.0762317,T=0.070789), scalar:-0.933553 ]
            position_entente_unscaled = position_entente;
            position_entente_unscaled.unscale();
            cout << "Unscaled, the position entente is " << position_entente_unscaled << endl;
            //Unscaled, the position entente is [ M_dna:(A=0.445995,C=0.180054,G=0.193897,T=0.180054) ]
            position_entente.normalize( 0 );
            cout << "The position entente, normalized (min=0), is " << position_entente << endl;
            //The position entente, normalized (min=0), is [ M_dna:(A=0.445995,C=0.180054,G=0.193897,T=0.180054) ]
            // Revert
            position_entente = position_entente_backup;
          
            emission_prior.m_mixingProbs[ 0 ] = INIT_PROBABILITY(ProbabilityType)( .9 );
            emission_prior.m_mixingProbs[ 1 ] = INIT_PROBABILITY(ProbabilityType)( .1 );
            emission_prior.incorporatePrior( position_entente );
            cout << "The position entente, after incorporating an UNeven mixture of those 2 priors (90% on the Laplace prior), is " << position_entente << endl;
            //The position entente, after incorporating an UNeven mixture of those 2 priors (90% on the Laplace prior), is [ M_dna:(A=0.20325,C=0.0608303,G=0.068244,T=0.0608303), scalar:-0.933553 ]
            position_entente_unscaled = position_entente;
            position_entente_unscaled.unscale();
            cout << "Unscaled, the position entente is " << position_entente_unscaled << endl;
            //Unscaled, the position entente is [ M_dna:(A=0.516972,C=0.154724,G=0.173581,T=0.154724) ]
            position_entente.normalize( 0 );
            cout << "The position entente, normalized (min=0), is " << position_entente << endl;
            //The position entente, normalized (min=0), is [ M_dna:(A=0.516972,C=0.154724,G=0.173581,T=0.154724) ]
            // Revert
            position_entente = position_entente_backup;
            // emission_prior[ 0 ] is A/G rich, [ 1 ] is C/T rich.
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 'a' ) ] = 10.0;
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 'c' ) ] = 1.0;
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 'g' ) ] = 10.0;
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 't' ) ] = 1.0;
            emission_prior[ 1 ][ galosh::Emission::Match ][ ( residue = 'a' ) ] = 1.0;
            emission_prior[ 1 ][ galosh::Emission::Match ][ ( residue = 'c' ) ] = 10.0;
            emission_prior[ 1 ][ galosh::Emission::Match ][ ( residue = 'g' ) ] = 1.0;
            emission_prior[ 1 ][ galosh::Emission::Match ][ ( residue = 't' ) ] = 10.0;
            emission_prior.m_mixingProbs[ 0 ] = INIT_PROBABILITY(ProbabilityType)( .5 );
            emission_prior.m_mixingProbs[ 1 ] = INIT_PROBABILITY(ProbabilityType)( .5 );
            emission_prior.incorporatePrior( position_entente );
            cout << "The position entente, after incorporating an even mixture of an A/G rich prior with a C/T rich prior, is " << position_entente << endl;
            //The position entente, after incorporating an even mixture of an A/G rich prior with a C/T rich prior, is [ M_dna:(A=0.19976,C=0.0160247,G=0.161345,T=0.0160247), scalar:-0.933553 ]
            position_entente_unscaled = position_entente;
            position_entente_unscaled.unscale();
            cout << "Unscaled, the position entente is " << position_entente_unscaled << endl;
            //Unscaled, the position entente is [ M_dna:(A=0.508096,C=0.0407594,G=0.410385,T=0.0407594) ]
            position_entente.normalize( 0 );
            cout << "The position entente, normalized (min=0), is " << position_entente << endl;
            //The position entente, normalized (min=0), is [ M_dna:(A=0.508096,C=0.0407594,G=0.410385,T=0.0407594) ]
            // Revert
            position_entente = position_entente_backup;
            emission_prior.m_mixingProbs[ 0 ] = INIT_PROBABILITY(ProbabilityType)( .1 );
            emission_prior.m_mixingProbs[ 1 ] = INIT_PROBABILITY(ProbabilityType)( .9 );
            emission_prior.incorporatePrior( position_entente );
            cout << "The position entente, after incorporating an UNeven mixture of those 2 priors (90% on the C/T rich prior), is " << position_entente << endl;
            //The position entente, after incorporating an UNeven mixture of those 2 priors (90% on the C/T rich prior), is [ M_dna:(A=0.192762,C=0.023023,G=0.154346,T=0.023023), scalar:-0.933553 ]
            position_entente_unscaled = position_entente;
            position_entente_unscaled.unscale();
            cout << "Unscaled, the position entente is " << position_entente_unscaled << endl;
            //Unscaled, the position entente is [ M_dna:(A=0.490296,C=0.0585596,G=0.392585,T=0.0585596) ]
            position_entente.normalize( 0 );
            cout << "The position entente, normalized (min=0), is " << position_entente << endl;
            //The position entente, normalized (min=0), is [ M_dna:(A=0.490296,C=0.0585596,G=0.392585,T=0.0585596) ]
            // Revert
            position_entente = position_entente_backup;
            // emission_prior[ 0 ] is A/T rich, [ 1 ] is C/G rich.
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 'a' ) ] = 10.0;
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 'c' ) ] = 1.0;
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 'g' ) ] = 1.0;
            emission_prior[ 0 ][ galosh::Emission::Match ][ ( residue = 't' ) ] = 10.0;
            emission_prior[ 1 ][ galosh::Emission::Match ][ ( residue = 'a' ) ] = 1.0;
            emission_prior[ 1 ][ galosh::Emission::Match ][ ( residue = 'c' ) ] = 10.0;
            emission_prior[ 1 ][ galosh::Emission::Match ][ ( residue = 'g' ) ] = 10.0;
            emission_prior[ 1 ][ galosh::Emission::Match ][ ( residue = 't' ) ] = 1.0;
            emission_prior.m_mixingProbs[ 0 ] = INIT_PROBABILITY(ProbabilityType)( .5 );
            emission_prior.m_mixingProbs[ 1 ] = INIT_PROBABILITY(ProbabilityType)( .5 );
            emission_prior.incorporatePrior( position_entente );
            cout << "The position entente, after incorporating an even mixture of an A/T rich prior with a C/G rich prior, is " << position_entente << endl;
            //The position entente, after incorporating an even mixture of an A/T rich prior with a C/G rich prior, is [ M_dna:(A=0.198163,C=0.0176218,G=0.0197314,T=0.157638), scalar:-0.933553 ]
            position_entente_unscaled = position_entente;
            position_entente_unscaled.unscale();
            cout << "Unscaled, the position entente is " << position_entente_unscaled << endl;
            //Unscaled, the position entente is [ M_dna:(A=0.504034,C=0.0448216,G=0.0501873,T=0.400957) ]
            position_entente.normalize( 0 );
            cout << "The position entente, normalized (min=0), is " << position_entente << endl;
            //The position entente, normalized (min=0), is [ M_dna:(A=0.504034,C=0.0448216,G=0.0501873,T=0.400957) ]
            // Revert
            position_entente = position_entente_backup;
            emission_prior.m_mixingProbs[ 0 ] = INIT_PROBABILITY(ProbabilityType)( .001 );
            emission_prior.m_mixingProbs[ 1 ] = INIT_PROBABILITY(ProbabilityType)( .999 );
            emission_prior.incorporatePrior( position_entente );
            cout << "The position entente, after incorporating an UNeven mixture of those 2 priors (99.9% on the C/G rich prior), is " << position_entente << endl;
            //The position entente, after incorporating an UNeven mixture of those 2 priors (99.9% on the C/G rich prior), is [ M_dna:(A=0.0675671,C=0.148218,G=0.150328,T=0.0270417), scalar:-0.933553 ]
            position_entente_unscaled = position_entente;
            position_entente_unscaled.unscale();
            cout << "Unscaled, the position entente is " << position_entente_unscaled << endl;
            //Unscaled, the position entente is [ M_dna:(A=0.171859,C=0.376997,G=0.382363,T=0.0687815) ]
            position_entente.normalize( 0 );
            cout << "The position entente, normalized (min=0), is " << position_entente << endl;
            //The position entente, normalized (min=0), is [ M_dna:(A=0.171859,C=0.376997,G=0.382363,T=0.0687815) ]
          
            // ERE I AM.  Testing transition priors
            //      galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::DirichletMixtureTransitionPrior<float> transition_prior( 1 );
            //      //transition_prior[ 0 ][ Transition::DNA ][ ( residue = 'a' ) ] = 1.0;
            //      //transition_prior[ 0 ][ Transition::DNA ][ ( residue = 'c' ) ] = 1.0;
            //      //transition_prior[ 0 ][ Transition::DNA ][ ( residue = 'g' ) ] = 1.0;
            //      //transition_prior[ 0 ][ Transition::DNA ][ ( residue = 't' ) ] = 1.0;
            //      cout << "The global entente, before incorporating any priors, is " << global_entente << endl;
            //      transition_prior.incorporatePrior( global_entente );
            //      cout << "The global entente, after incorporating a simple Laplace prior, is " << global_entente << endl;
            //      //The global entente, after incorporating a simple Laplace prior, is [ M_dna:(A=0.208683,C=0.0588913,G=0.0666888,T=0.0588913), scalar:-0.933553 ]
            //      global_entente_unscaled = global_entente;
            //      global_entente_unscaled.unscale();
            //      cout << "Unscaled, the global entente is " << global_entente_unscaled << endl;
            //      //Unscaled, the global entente is [ M_dna:(A=0.530792,C=0.149792,G=0.169625,T=0.149792) ]
            //      global_entente.normalize( 0 );
            //      cout << "The global entente, normalized (min=0), is " << global_entente << endl;
            //      //The global entente, normalized (min=0), is [ M_dna:(A=0.530792,C=0.149792,G=0.169625,T=0.149792) ]
            //      // Revert
            //      global_entente = global_entente_backup;
          
          } // End if test_priors
        } // End if test_dp_basics || test_alignment_profiles || test_priors
  
#ifdef ALLOW_BOLTZMANN_GIBBS
        // Test Boltzmann-Gibbs conversion
        if( test_boltzmann_gibbs ) {
          //DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::PositionBoltzmannGibbs pos_boltzmann_gibbs( dna_profile[ 0 ] );
          galosh::PositionBoltzmannGibbs<Dna, ScoreType> pos_boltzmann_gibbs( dna_profile[ 0 ] );
          galosh::PositionSpecificParameters<Dna, ProbabilityType> pos_unboltzmann_gibbs;
          pos_boltzmann_gibbs.toProfilePosition( pos_unboltzmann_gibbs );
          cout << "The first profile position, " << dna_profile[ 0 ] << ", represented in Boltzmann-Gibbs form with default temp (1.0), is " << pos_boltzmann_gibbs << ".  Converted back to regular representation, it is " << pos_unboltzmann_gibbs << "." << endl;
          // The first profile position, [ M_DNA:(A=0.7,C=0.1,G=0.1,T=0.1) ], represented in Boltzmann-Gibbs form with default temp (1.0), is [ M_DNA:(A=-0.356675,C=-2.30259,G=-2.30259,T=-2.30259) ].  Converted back to regular representation, it is [ M_DNA:(A=0.7,C=0.1,G=0.1,T=0.1) ].
          //DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::PositionBoltzmannGibbs pos_boltzmann_gibbs_temp100( dna_profile[ 0 ], 100.0 );
          galosh::PositionBoltzmannGibbs<Dna, ScoreType> pos_boltzmann_gibbs_temp100( dna_profile[ 0 ], 100.0 );
          pos_boltzmann_gibbs_temp100.toProfilePosition( pos_unboltzmann_gibbs );
          cout << "The first profile position, " << dna_profile[ 0 ] << ", represented in Boltzmann-Gibbs form with temperature 100.0, is " << pos_boltzmann_gibbs_temp100 << ".  Converted back to regular representation, it is " << pos_unboltzmann_gibbs << "." << endl;
          // The first profile position, [ M_DNA:(A=0.7,C=0.1,G=0.1,T=0.1) ], represented in Boltzmann-Gibbs form with temperature 100.0, is [ M_DNA:(A=-0.00356675,C=-0.0230259,G=-0.0230259,T=-0.0230259) ].  Converted back to regular representation, it is [ M_DNA:(A=0.7,C=0.1,G=0.1,T=0.1) ].
        } // End if test_boltzmann_gibbs
#endif // ALLOW_BOLTZMANN_GIBBS
      
        if( test_cross_entropy ) {
          // First, on one position only
          double self_entropy_0 =
            dna_profile[ 0 ].crossEntropy( dna_profile[ 0 ] );
          cout << "The first profile position's (self-)entropy is " << self_entropy_0 << endl;
      
          double self_entropy_1 =
            dna_profile[ 1 ].crossEntropy( dna_profile[ 1 ] );
          cout << "The second profile position's (self-)entropy is " << self_entropy_1 << endl;
      
          double cross_entropy_01 =
            dna_profile[ 0 ].crossEntropy( dna_profile[ 1 ] );
          cout << "The first profile position's cross-entropy with the second position is " << cross_entropy_01 << endl;
          cout << "The first profile position's KL divergence with the second position is " << ( cross_entropy_01 - self_entropy_0 ) << endl;
      
          double cross_entropy_10 =
            dna_profile[ 1 ].crossEntropy( dna_profile[ 0 ] );
          cout << "The second profile position's cross-entropy with the first position is " << cross_entropy_10 << endl;
          cout << "The second profile position's KL divergence with the first position is " << ( cross_entropy_10 - self_entropy_1 ) << endl;
      
          cout << "The symmeterized KL divergence between the two positions is " << ( ( cross_entropy_01 - self_entropy_0 ) + ( cross_entropy_10 - self_entropy_1 ) ) << endl;
      
          galosh::ProfilePosition<Dna, ProbabilityType> evenpos( &dna_profile );
          evenpos.even();
          double self_entropy_evenpos =
            evenpos.crossEntropy( evenpos );
          cout << "The even position's (self-)entropy is " << self_entropy_evenpos << endl;
      
          double cross_entropy_0evenpos =
            dna_profile[ 0 ].crossEntropy( evenpos );
          cout << "The first profile position's cross-entropy with the even position is " << cross_entropy_0evenpos << endl;
          cout << "The first profile position's KL divergence with the even position is " << ( cross_entropy_0evenpos - self_entropy_0 ) << endl;
      
          double cross_entropy_evenpos0 =
            evenpos.crossEntropy( dna_profile[ 0 ] );
          cout << "The even position's cross-entropy with the first position is " << cross_entropy_evenpos0 << endl;
          cout << "The even position's KL divergence with the first position is " << ( cross_entropy_evenpos0 - self_entropy_evenpos ) << endl;
      
          cout << "The symmeterized KL divergence between the two positions is " << ( ( cross_entropy_0evenpos - self_entropy_0 ) + ( cross_entropy_evenpos0 - self_entropy_evenpos ) ) << endl;
      
          galosh::ProfilePosition<Dna, ProbabilityType> A3pos( &dna_profile );
          A3pos[ galosh::Emission::Match ][ ( residue = 'a' ) ] = .3;
          A3pos[ galosh::Emission::Match ][ ( residue = 'c' ) ] = ((1.0-.3)/3);
          A3pos[ galosh::Emission::Match ][ ( residue = 'g' ) ] = ((1.0-.3)/3);
          A3pos[ galosh::Emission::Match ][ ( residue = 't' ) ] = ((1.0-.3)/3);
          A3pos.normalize( 0.0 );
          double self_entropy_A3pos =
            A3pos.crossEntropy( A3pos );
          cout << "The A3 position's (self-)entropy is " << self_entropy_A3pos << endl;
      
          double cross_entropy_A3posevenpos =
            A3pos.crossEntropy( evenpos );
          cout << "The A3 position's cross-entropy with the even position is " << cross_entropy_A3posevenpos << endl;
          cout << "The A3 position's KL divergence with the even position is " << ( cross_entropy_A3posevenpos - self_entropy_A3pos ) << endl;
      
          double cross_entropy_evenposA3pos =
            evenpos.crossEntropy( A3pos );
          cout << "The even position's cross-entropy with the A3 position is " << cross_entropy_evenposA3pos << endl;
          cout << "The even position's KL divergence with the A3 position is " << ( cross_entropy_evenposA3pos - self_entropy_evenpos ) << endl;
      
          cout << "The symmeterized KL divergence between the two positions is " << ( ( cross_entropy_A3posevenpos - self_entropy_A3pos ) + ( cross_entropy_evenposA3pos - self_entropy_evenpos ) ) << endl;
      
          galosh::ProfilePosition<Dna, ProbabilityType> T3pos( &dna_profile );
          T3pos[ galosh::Emission::Match ][ ( residue = 't' ) ] = .3;
          T3pos[ galosh::Emission::Match ][ ( residue = 'a' ) ] = ((1.0-.3)/3);
          T3pos[ galosh::Emission::Match ][ ( residue = 'c' ) ] = ((1.0-.3)/3);
          T3pos[ galosh::Emission::Match ][ ( residue = 'g' ) ] = ((1.0-.3)/3);
          T3pos.normalize( 0.0 );
          double self_entropy_T3pos =
            T3pos.crossEntropy( T3pos );
          cout << "The T3 position's (self-)entropy is " << self_entropy_T3pos << endl;
      
          double cross_entropy_A3posT3pos =
            A3pos.crossEntropy( T3pos );
          cout << "The A3 position's cross-entropy with the T3 position is " << cross_entropy_A3posT3pos << endl;
          cout << "The A3 position's KL divergence with the T3 position is " << ( cross_entropy_A3posT3pos - self_entropy_A3pos ) << endl;
      
          double cross_entropy_T3posA3pos =
            T3pos.crossEntropy( A3pos );
          cout << "The T3 position's cross-entropy with the A3 position is " << cross_entropy_T3posA3pos << endl;
          cout << "The T3 position's KL divergence with the A3 position is " << ( cross_entropy_T3posA3pos - self_entropy_T3pos ) << endl;
      
          cout << "The symmeterized KL divergence between the two positions is " << ( ( cross_entropy_A3posT3pos - self_entropy_A3pos ) + ( cross_entropy_T3posA3pos - self_entropy_T3pos ) ) << endl;
      
          galosh::ProfilePosition<Dna, ProbabilityType> T4pos( &dna_profile );
          T4pos[ galosh::Emission::Match ][ ( residue = 't' ) ] = .4;
          T4pos[ galosh::Emission::Match ][ ( residue = 'a' ) ] = ((1.0-.4)/3);
          T4pos[ galosh::Emission::Match ][ ( residue = 'c' ) ] = ((1.0-.4)/3);
          T4pos[ galosh::Emission::Match ][ ( residue = 'g' ) ] = ((1.0-.4)/3);
          T4pos.normalize( 0.0 );
          double self_entropy_T4pos =
            T4pos.crossEntropy( T4pos );
          cout << "The T4 position's (self-)entropy is " << self_entropy_T4pos << endl;
      
          double cross_entropy_A3posT4pos =
            A3pos.crossEntropy( T4pos );
          cout << "The A3 position's cross-entropy with the T4 position is " << cross_entropy_A3posT4pos << endl;
          cout << "The A3 position's KL divergence with the T4 position is " << ( cross_entropy_A3posT4pos - self_entropy_A3pos ) << endl;
      
          double cross_entropy_T4posA3pos =
            T4pos.crossEntropy( A3pos );
          cout << "The T4 position's cross-entropy with the A3 position is " << cross_entropy_T4posA3pos << endl;
          cout << "The T4 position's KL divergence with the A3 position is " << ( cross_entropy_T4posA3pos - self_entropy_T4pos ) << endl;
      
          cout << "The symmeterized KL divergence between the two positions is " << ( ( cross_entropy_A3posT4pos - self_entropy_A3pos ) + ( cross_entropy_T4posA3pos - self_entropy_T4pos ) ) << endl;
      
          galosh::ProfilePosition<Dna, ProbabilityType> T5pos( &dna_profile );
          T5pos[ galosh::Emission::Match ][ ( residue = 't' ) ] = .5;
          T5pos[ galosh::Emission::Match ][ ( residue = 'a' ) ] = ((1.0-.5)/3);
          T5pos[ galosh::Emission::Match ][ ( residue = 'c' ) ] = ((1.0-.5)/3);
          T5pos[ galosh::Emission::Match ][ ( residue = 'g' ) ] = ((1.0-.5)/3);
          T5pos.normalize( 0.0 );
          double self_entropy_T5pos =
            T5pos.crossEntropy( T5pos );
          cout << "The T5 position's (self-)entropy is " << self_entropy_T5pos << endl;
      
          double cross_entropy_A3posT5pos =
            A3pos.crossEntropy( T5pos );
          cout << "The A3 position's cross-entropy with the T5 position is " << cross_entropy_A3posT5pos << endl;
          cout << "The A3 position's KL divergence with the T5 position is " << ( cross_entropy_A3posT5pos - self_entropy_A3pos ) << endl;
      
          double cross_entropy_T5posA3pos =
            T5pos.crossEntropy( A3pos );
          cout << "The T5 position's cross-entropy with the A3 position is " << cross_entropy_T5posA3pos << endl;
          cout << "The T5 position's KL divergence with the A3 position is " << ( cross_entropy_T5posA3pos - self_entropy_T5pos ) << endl;
      
          cout << "The symmeterized KL divergence between the two positions is " << ( ( cross_entropy_A3posT5pos - self_entropy_A3pos ) + ( cross_entropy_T5posA3pos - self_entropy_T4pos ) ) << endl;
      
          galosh::ProfilePosition<Dna, ProbabilityType> T6pos( &dna_profile );
          T6pos[ galosh::Emission::Match ][ ( residue = 't' ) ] = .6;
          T6pos[ galosh::Emission::Match ][ ( residue = 'a' ) ] = ((1.0-.6)/3);
          T6pos[ galosh::Emission::Match ][ ( residue = 'c' ) ] = ((1.0-.6)/3);
          T6pos[ galosh::Emission::Match ][ ( residue = 'g' ) ] = ((1.0-.6)/3);
          T6pos.normalize( 0.0 );
          double self_entropy_T6pos =
            T6pos.crossEntropy( T6pos );
          cout << "The T6 position's (self-)entropy is " << self_entropy_T6pos << endl;
      
          double cross_entropy_A3posT6pos =
            A3pos.crossEntropy( T6pos );
          cout << "The A3 position's cross-entropy with the T6 position is " << cross_entropy_A3posT6pos << endl;
          cout << "The A3 position's KL divergence with the T6 position is " << ( cross_entropy_A3posT6pos - self_entropy_A3pos ) << endl;
      
          double cross_entropy_T6posA3pos =
            T6pos.crossEntropy( A3pos );
          cout << "The T6 position's cross-entropy with the A3 position is " << cross_entropy_T6posA3pos << endl;
          cout << "The T6 position's KL divergence with the A3 position is " << ( cross_entropy_T6posA3pos - self_entropy_T6pos ) << endl;
      
          cout << "The symmeterized KL divergence between the two positions is " << ( ( cross_entropy_A3posT6pos - self_entropy_A3pos ) + ( cross_entropy_T6posA3pos - self_entropy_T6pos ) ) << endl;
      
          galosh::ProfilePosition<Dna, ProbabilityType> A4pos( &dna_profile );
          A4pos[ galosh::Emission::Match ][ ( residue = 'a' ) ] = .4;
          A4pos[ galosh::Emission::Match ][ ( residue = 'c' ) ] = ((1.0-.4)/3);
          A4pos[ galosh::Emission::Match ][ ( residue = 'g' ) ] = ((1.0-.4)/3);
          A4pos[ galosh::Emission::Match ][ ( residue = 't' ) ] = ((1.0-.4)/3);
          A4pos.normalize( 0.0 );
          double self_entropy_A4pos =
            A4pos.crossEntropy( A4pos );
          cout << "The A4 position's (self-)entropy is " << self_entropy_A4pos << endl;
      
          double cross_entropy_A3posA4pos =
            A3pos.crossEntropy( A4pos );
          cout << "The A3 position's cross-entropy with the A4 position is " << cross_entropy_A3posA4pos << endl;
          cout << "The A3 position's KL divergence with the A4 position is " << ( cross_entropy_A3posA4pos - self_entropy_A3pos ) << endl;
      
          double cross_entropy_A4posA3pos =
            A4pos.crossEntropy( A3pos );
          cout << "The A4 position's cross-entropy with the A3 position is " << cross_entropy_A4posA3pos << endl;
          cout << "The A4 position's KL divergence with the A3 position is " << ( cross_entropy_A4posA3pos - self_entropy_A4pos ) << endl;
      
          cout << "The symmeterized KL divergence between the two positions is " << ( ( cross_entropy_A3posA4pos - self_entropy_A3pos ) + ( cross_entropy_A4posA3pos - self_entropy_A4pos ) ) << endl;
      
          galosh::ProfilePosition<Dna, ProbabilityType> A5pos( &dna_profile );
          A5pos[ galosh::Emission::Match ][ ( residue = 'a' ) ] = .5;
          A5pos[ galosh::Emission::Match ][ ( residue = 'c' ) ] = ((1.0-.5)/3);
          A5pos[ galosh::Emission::Match ][ ( residue = 'g' ) ] = ((1.0-.5)/3);
          A5pos[ galosh::Emission::Match ][ ( residue = 't' ) ] = ((1.0-.5)/3);
          A5pos.normalize( 0.0 );
          double self_entropy_A5pos =
            A5pos.crossEntropy( A5pos );
          cout << "The A5 position's (self-)entropy is " << self_entropy_A5pos << endl;
      
          double cross_entropy_A3posA5pos =
            A3pos.crossEntropy( A5pos );
          cout << "The A3 position's cross-entropy with the A5 position is " << cross_entropy_A3posA5pos << endl;
          cout << "The A3 position's KL divergence with the A5 position is " << ( cross_entropy_A3posA5pos - self_entropy_A3pos ) << endl;
      
          double cross_entropy_A5posA3pos =
            A5pos.crossEntropy( A3pos );
          cout << "The A5 position's cross-entropy with the A3 position is " << cross_entropy_A5posA3pos << endl;
          cout << "The A5 position's KL divergence with the A3 position is " << ( cross_entropy_A5posA3pos - self_entropy_A5pos ) << endl;
      
          cout << "The symmeterized KL divergence between the two positions is " << ( ( cross_entropy_A3posA5pos - self_entropy_A3pos ) + ( cross_entropy_A5posA3pos - self_entropy_A5pos ) ) << endl;
      
          galosh::ProfilePosition<Dna, ProbabilityType> A6pos( &dna_profile );
          A6pos[ galosh::Emission::Match ][ ( residue = 'a' ) ] = .6;
          A6pos[ galosh::Emission::Match ][ ( residue = 'c' ) ] = ((1.0-.6)/3);
          A6pos[ galosh::Emission::Match ][ ( residue = 'g' ) ] = ((1.0-.6)/3);
          A6pos[ galosh::Emission::Match ][ ( residue = 't' ) ] = ((1.0-.6)/3);
          A6pos.normalize( 0.0 );
          double self_entropy_A6pos =
            A6pos.crossEntropy( A6pos );
          cout << "The A6 position's (self-)entropy is " << self_entropy_A6pos << endl;
      
          double cross_entropy_A3posA6pos =
            A3pos.crossEntropy( A6pos );
          cout << "The A3 position's cross-entropy with the A6 position is " << cross_entropy_A3posA6pos << endl;
          cout << "The A3 position's KL divergence with the A6 position is " << ( cross_entropy_A3posA6pos - self_entropy_A3pos ) << endl;
      
          double cross_entropy_A6posA3pos =
            A6pos.crossEntropy( A3pos );
          cout << "The A6 position's cross-entropy with the A3 position is " << cross_entropy_A6posA3pos << endl;
          cout << "The A6 position's KL divergence with the A3 position is " << ( cross_entropy_A6posA3pos - self_entropy_A6pos ) << endl;
      
          cout << "The symmeterized KL divergence between the two positions is " << ( ( cross_entropy_A3posA6pos - self_entropy_A3pos ) + ( cross_entropy_A6posA3pos - self_entropy_A6pos ) ) << endl;
      
          double cross_entropy_A6posT3pos =
            A6pos.crossEntropy( T3pos );
          cout << "The A6 position's cross-entropy with the T3 position is " << cross_entropy_A6posT3pos << endl;
          cout << "The A6 position's KL divergence with the T3 position is " << ( cross_entropy_A6posT3pos - self_entropy_A6pos ) << endl;
      
          double cross_entropy_T3posA6pos =
            T3pos.crossEntropy( A6pos );
          cout << "The T3 position's cross-entropy with the A6 position is " << cross_entropy_T3posA6pos << endl;
          cout << "The T3 position's KL divergence with the A6 position is " << ( cross_entropy_T3posA6pos - self_entropy_T3pos ) << endl;
      
          cout << "The symmeterized KL divergence between the two positions is " << ( ( cross_entropy_A6posT3pos - self_entropy_A6pos ) + ( cross_entropy_T3posA6pos - self_entropy_T3pos ) ) << endl;
      
          double cross_entropy_A6posT4pos =
            A6pos.crossEntropy( T4pos );
          cout << "The A6 position's cross-entropy with the T4 position is " << cross_entropy_A6posT4pos << endl;
          cout << "The A6 position's KL divergence with the T4 position is " << ( cross_entropy_A6posT4pos - self_entropy_A6pos ) << endl;
      
          double cross_entropy_T4posA6pos =
            T4pos.crossEntropy( A6pos );
          cout << "The T4 position's cross-entropy with the A6 position is " << cross_entropy_T4posA6pos << endl;
          cout << "The T4 position's KL divergence with the A6 position is " << ( cross_entropy_T4posA6pos - self_entropy_T4pos ) << endl;
      
          cout << "The symmeterized KL divergence between the two positions is " << ( ( cross_entropy_A6posT4pos - self_entropy_A6pos ) + ( cross_entropy_T4posA6pos - self_entropy_T4pos ) ) << endl;
      
          double cross_entropy_A6posT5pos =
            A6pos.crossEntropy( T5pos );
          cout << "The A6 position's cross-entropy with the T5 position is " << cross_entropy_A6posT5pos << endl;
          cout << "The A6 position's KL divergence with the T5 position is " << ( cross_entropy_A6posT5pos - self_entropy_A6pos ) << endl;
      
          double cross_entropy_T5posA6pos =
            T5pos.crossEntropy( A6pos );
          cout << "The T5 position's cross-entropy with the A6 position is " << cross_entropy_T5posA6pos << endl;
          cout << "The T5 position's KL divergence with the A6 position is " << ( cross_entropy_T5posA6pos - self_entropy_T5pos ) << endl;
      
          cout << "The symmeterized KL divergence between the two positions is " << ( ( cross_entropy_A6posT5pos - self_entropy_A6pos ) + ( cross_entropy_T5posA6pos - self_entropy_T4pos ) ) << endl;
      
          double cross_entropy_A6posT6pos =
            A6pos.crossEntropy( T6pos );
          cout << "The A6 position's cross-entropy with the T6 position is " << cross_entropy_A6posT6pos << endl;
          cout << "The A6 position's KL divergence with the T6 position is " << ( cross_entropy_A6posT6pos - self_entropy_A6pos ) << endl;
      
          double cross_entropy_T6posA6pos =
            T6pos.crossEntropy( A6pos );
          cout << "The T6 position's cross-entropy with the A6 position is " << cross_entropy_T6posA6pos << endl;
          cout << "The T6 position's KL divergence with the A6 position is " << ( cross_entropy_T6posA6pos - self_entropy_T6pos ) << endl;
      
          cout << "The symmeterized KL divergence between the two positions is " << ( ( cross_entropy_A6posT6pos - self_entropy_A6pos ) + ( cross_entropy_T6posA6pos - self_entropy_T6pos ) ) << endl;
      
      
          // Now try weighting it.  With this kind of weights, the "entropy" and
          // "KL divergence" don't have their usual properties.  To calculate
          // weighted entropy or KL divergence, weights within any Multinomial
          // distribution should be the same.
          galosh::MatchEmissionParameters<Dna, double> mep_weights;
          mep_weights[ galosh::Emission::Match ][ ( residue = 'a' ) ] = 1.0;
          mep_weights[ galosh::Emission::Match ][ ( residue = 'c' ) ] = 0.0;
          mep_weights[ galosh::Emission::Match ][ ( residue = 'g' ) ] = 0.0;
          mep_weights[ galosh::Emission::Match ][ ( residue = 't' ) ] = 0.0;
          double weighted_self_entropy_0 =
            dna_profile[ 0 ].crossEntropy( dna_profile[ 0 ], &mep_weights );
          cout << "The first profile position's (self-)entropy *of the \"A\" probability only* is " << weighted_self_entropy_0 << endl;
      
          double weighted_self_entropy_1 =
            dna_profile[ 1 ].crossEntropy( dna_profile[ 1 ], &mep_weights );
          cout << "The second profile position's (self-)entropy *of the \"A\" probability only* is " << weighted_self_entropy_1 << endl;
      
          double weighted_cross_entropy_01 =
            dna_profile[ 0 ].crossEntropy( dna_profile[ 1 ], &mep_weights );
          cout << "The first profile position's cross-entropy with the second position *of the \"A\" probability only* is " << weighted_cross_entropy_01 << endl;
          //cout << "The first profile position's KL divergence with the second position *of the \"A\" probability only* is " << ( weighted_cross_entropy_01 - weighted_self_entropy_0 ) << endl;
      
          double weighted_cross_entropy_10 =
            dna_profile[ 1 ].crossEntropy( dna_profile[ 0 ], &mep_weights );
          cout << "The second profile position's cross-entropy with the first position *of the \"A\" probability only* is " << weighted_cross_entropy_10 << endl;
          //cout << "The second profile position's KL divergence with the first position *of the \"A\" probability only* is " << ( weighted_cross_entropy_10 - weighted_self_entropy_1 ) << endl;
      
          //cout << "The symmeterized KL divergence between the two positions *of the \"A\" probability only* is " << ( ( weighted_cross_entropy_01 - weighted_self_entropy_0 ) + ( weighted_cross_entropy_10 - weighted_self_entropy_1 ) ) << endl;
      
          // Now on a whole profile
          double self_entropy_p =
            dna_profile.crossEntropy( dna_profile );
          cout << "The profile's (self-)entropy is " << self_entropy_p << endl;
      
          ProfileType even_profile = dna_profile;
          even_profile.even();
          double self_entropy_e =
            even_profile.crossEntropy( even_profile );
          cout << "The even profile's (self-)entropy is " << self_entropy_e << endl;
      
          double cross_entropy_pe =
            dna_profile.crossEntropy( even_profile );
          cout << "The profile's cross-entropy with the even profile is " << cross_entropy_pe << endl;
          cout << "The profile's KL divergence with the even profile is " << ( cross_entropy_pe - self_entropy_p ) << endl;
      
          double cross_entropy_ep =
            even_profile.crossEntropy( dna_profile );
          cout << "The even profile's cross-entropy with the profile is " << cross_entropy_ep << endl;
          cout << "The even profile's KL divergence with the profile is " << ( cross_entropy_ep - self_entropy_e ) << endl;
      
          cout << "The symmeterized KL divergence between the two profiles is " << ( ( cross_entropy_pe - self_entropy_p ) + ( cross_entropy_ep - self_entropy_e ) ) << endl;
      
          // Now try weighting it.
          galosh::ProfileTreeRoot<Dna, double> prof_weights( dna_profile.length() );
          prof_weights.zero();
          prof_weights[ 0 ][ galosh::Emission::Match ][ ( residue = 'a' ) ] = 1.0;
          prof_weights[ 0 ][ galosh::Emission::Match ][ ( residue = 'c' ) ] = 1.0;
          prof_weights[ 0 ][ galosh::Emission::Match ][ ( residue = 'g' ) ] = 1.0;
          prof_weights[ 0 ][ galosh::Emission::Match ][ ( residue = 't' ) ] = 1.0;
          double weighted_self_entropy_p =
            dna_profile.crossEntropy( dna_profile, &prof_weights );
          cout << "The profile's (self-)entropy *at the first position only* is " << weighted_self_entropy_p << endl;
      
          double weighted_self_entropy_e =
            even_profile.crossEntropy( even_profile, &prof_weights );
          cout << "The even profile's (self-)entropy *at the first position only* is " << weighted_self_entropy_e << endl;
      
          double weighted_cross_entropy_pe =
            dna_profile.crossEntropy( even_profile, &prof_weights );
          cout << "The profile's cross-entropy with the even profile *at the first position only* is " << weighted_cross_entropy_pe << endl;
          cout << "The profile's KL divergence with the even profile *at the first position only* is " << ( weighted_cross_entropy_pe - weighted_self_entropy_p ) << endl;
      
          double weighted_cross_entropy_ep =
            even_profile.crossEntropy( dna_profile, &prof_weights );
          cout << "The even profile's cross-entropy with the profile *at the first position only* is " << weighted_cross_entropy_ep << endl;
          cout << "The even profile's KL divergence with the profile *at the first position only* is " << ( weighted_cross_entropy_ep - weighted_self_entropy_e ) << endl;
      
          cout << "The symmeterized KL divergence between the two profiles *at the first position only* is " << ( ( weighted_cross_entropy_pe - weighted_self_entropy_p ) + ( weighted_cross_entropy_ep - weighted_self_entropy_e ) ) << endl;
      
          // TODO: REMOVE!
          exit( 0 );
        } // End if test_cross_entropy
      
        if( test_profile_profile_alignment ) {
          random.setSeed( static_cast<uint32_t>( std::time( NULL ) ) );
          cout << "After reseeding with seed " << random.getSeed() << ":" << endl;
      
          ProfileType profile1( 10 );
          profile1.uniform( random );
          ProfileType profile2( 10 );
          profile2.uniform( random );
          for( uint32_t i = 0; i < 5; i++ ) {
            profile2[ i ] = profile1[ i ];
          }
          for( uint32_t i = 8; i < 10; i++ ) {
            profile2[ i ] = profile1[ i - 2 ];
          }
      
          vector<uint32_t> alignment;
          double cost;
      
          double indel_open_cost;
          double indel_extension_cost;
      
          for( indel_open_cost = 0; indel_open_cost <= 0.5; indel_open_cost += .0625 ) {
            indel_extension_cost = indel_open_cost;
            cout << "indel_open_cost = " << indel_open_cost << endl;
            cout << "indel_extension_cost = " << indel_extension_cost << endl;
            cost =
              dp.profileProfile_align_SKL(
                dp_parameters,
                profile1,
                profile2,
                indel_open_cost,
                indel_extension_cost,
                alignment
              );
            cout << "After aligning the two profiles, got alignment cost: " << cost << endl;
            cout << "The alignment is ( " << alignment[ 0 ];
            uint32_t num_opens = 0;
            uint32_t num_extensions = 0;
            uint32_t num_aligned_positions = 0;
            if( alignment[ 0 ] > 0 ) {
              num_opens = 1;
              num_extensions = alignment[ 0 ] - 1;
            }
            for( uint32_t i = 1; i < alignment.size(); i++ ) {
              cout << ", " << alignment[ i ];
              if( alignment[ i ] == 0 ) {
                if( alignment[ i - 1 ] == 0 ) {
                  num_extensions += 1;
                } else {
                  num_opens += 1;
                }
              } else {
                num_aligned_positions += 1;
                if( alignment[ i ] > 1 ) {
                  num_opens += 1;
                  num_extensions += ( alignment[ i ] - 2 );
                }
              }
            }  // End foreach alignment position
            cout << " )" << endl;
            cout << "There were " << num_opens << " gap opens." << endl;
            cout << "There were " << num_extensions << " gap extensions." << endl;
            cout << "There were " << num_aligned_positions << " aligned positions." << endl;
            cout << "The symmeterized KL divergence of the aligned positions is " << 
              ( cost - ( num_opens * indel_open_cost ) - ( num_extensions * indel_extension_cost ) ) << endl;
            cout << "The average symmeterized KL divergence over the aligned positions is " << 
              ( ( cost - ( num_opens * indel_open_cost ) - ( num_extensions * indel_extension_cost ) ) / num_aligned_positions ) << endl;
          } // End foreach indel_open_cost
      
          // TODO: REMOVE!
          exit( 0 );
        } // End if test_profile_profile_alignment
      
        if( test_path_cooccurrence ) {
          ScoreType cooc_prob =
            dp.calculatePathCooccurrenceProbability(
              dp_parameters,
              dna_profile,
              dna_profile
            );
          cout << "The profile's path cooccurrence probability (with itself) is " << cooc_prob << endl;
      
          // Make a profile with exactly one path.
          ProfileType another_profile( 2 );
#ifndef DISALLOW_FLANKING_TRANSITIONS
          another_profile[ galosh::Transition::fromPreAlign ][ galosh::TransitionFromPreAlign::toPreAlign ] = 0;
          another_profile[ galosh::Transition::fromPreAlign ][ galosh::TransitionFromPreAlign::toBegin ] = 1;
#endif // !DISALLOW_FLANKING_TRANSITIONS
          another_profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletion ] = 0;
          another_profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toMatch ] = 1;
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
          another_profile[ galosh::Emission::PreAlignInsertion ].even();
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
          
          another_profile[ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toInsertion ] = 0;
          another_profile[ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletion ] = 0;
          another_profile[ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toMatch ] = 1;
      
          another_profile[ galosh::Transition::fromInsertion ][ galosh::TransitionFromInsertion::toInsertion ] = 0;
          another_profile[ galosh::Transition::fromInsertion ][ galosh::TransitionFromInsertion::toMatch ] = 1;
          another_profile[ galosh::Transition::fromDeletion ][ galosh::TransitionFromDeletion::toDeletion ] = 0;
          another_profile[ galosh::Transition::fromDeletion ][ galosh::TransitionFromDeletion::toMatch ] = 1;
          
          another_profile[ galosh::Emission::Insertion ].even();
      
          another_profile[ 0 ][ galosh::Emission::Match ][ ( residue = 'a' ) ] = 1;
          another_profile[ 0 ][ galosh::Emission::Match ][ ( residue = 'c' ) ] = 0;
          another_profile[ 0 ][ galosh::Emission::Match ][ ( residue = 'g' ) ] = 0;
          another_profile[ 0 ][ galosh::Emission::Match ][ ( residue = 't' ) ] = 0;
      
          another_profile[ 1 ][ galosh::Emission::Match ][ ( residue = 'a' ) ] = 0;
          another_profile[ 1 ][ galosh::Emission::Match ][ ( residue = 'c' ) ] = 1;
          another_profile[ 1 ][ galosh::Emission::Match ][ ( residue = 'g' ) ] = 0;
          another_profile[ 1 ][ galosh::Emission::Match ][ ( residue = 't' ) ] = 0;
      
#ifdef USE_END_DISTRIBUTION
          // For now we we never allow loops, so the End distribution is constant.
          another_profile[ galosh::Transition::fromEnd ][ galosh::TransitionFromEnd::toPostAlign ] = 1;
          another_profile[ galosh::Transition::fromEnd ][ galosh::TransitionFromEnd::toLoop ] = 0;
#endif // USE_END_DISTRIBUTION
#ifndef DISALLOW_FLANKING_TRANSITIONS
          another_profile[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toPostAlign ] = 0;
          another_profile[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toTerminal ] = 1;
#endif // !DISALLOW_FLANKING_TRANSITIONS
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
          another_profile[ galosh::Emission::PostAlignInsertion ].even();
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
      
          cooc_prob =
            dp.calculatePathCooccurrenceProbability(
              dp_parameters,
              another_profile,
              another_profile
            );
          cout << "The one-path profile's path cooccurrence probability (with itself) is " << cooc_prob << endl;
          cooc_prob = 1;
          cout << "This should be " << cooc_prob << endl;
      
          cooc_prob =
            dp.calculatePathCooccurrenceProbability(
              dp_parameters,
              dna_profile,
              another_profile
            );
          cout << "The one-path profile's path cooccurrence probability with the other profile is " << cooc_prob << endl;
      
          another_profile[ 1 ].copyFrom( dna_profile[ 1 ] );
          cooc_prob =
            dp.calculatePathCooccurrenceProbability(
              dp_parameters,
              dna_profile,
              another_profile
            );
          cout << "After setting their second position's match emission probabilties to be equal, the one-path profile's path cooccurrence probability with the other profile is " << cooc_prob << endl;
      
          another_profile[ 0 ].copyFrom( dna_profile[ 0 ] );
          cooc_prob =
            dp.calculatePathCooccurrenceProbability(
              dp_parameters,
              dna_profile,
              another_profile
            );
          cout << "After setting both of their positions' match emission probabilties to be equal, the one-path profile's path cooccurrence probability with the other profile is " << cooc_prob << endl;
      
          // TODO: REMOVE!
          exit( 0 );
        } // End if test_path_cooccurrence
      } // End if test_dynamic_programming

      if( test_trainer ) {

        if( test_ambiguous ) {
          dna_seqs[ 0 ] = "a";
          dna_seqs[ 1 ] = "an";
          dna_seqs[ 2 ] = "ng";
        } // End if test_ambiguous && test_trainer

        galosh::ProfileTrainer<ProfileType, ScoreType, MatrixValueType, Dna5> trainer(
          &dna_profile,
          dna_seqs,
          num_sequences_to_use
        );
        random.setSeed( 1 );
        //random.setSeed( std::time( NULL ) );

        dna_profile.uniformPositions( random );

        trainer.m_parameters.usePriors = false;
        trainer.m_parameters.alwaysAccept = false;

        //trainer.m_parameters.useUnconditionalBaumWelch = true;
        //training_parameters_template.unconditionalIsolatesGlobals = false;
        //trainer.m_parameters.useRabinerScaling = false;
        //trainer.m_parameters.trainProfileGlobals = false;
        //trainer.m_parameters.trainProfilePositions = false;
        //trainer.m_parameters.profileValueMinimum = 0;
      
        // For testing Baldi-style gradient ascent
#ifdef ALLOW_BOLTZMANN_GIBBS
        // NOTE: This has NO EFFECT unless  ALLOW_BOLTZMANN_GIBBS is #defined and true 
        trainer.m_parameters.baldiLearningRate = 1;
        trainer.m_parameters.baldiTemperature = 1;
        trainer.m_parameters.baldiHybrid = false;
        trainer.m_parameters.siegelMaxFindingThePeakAttempts_positions = 1000;
        trainer.m_parameters.siegelEpsilonScaleFactor = 1.5;
        trainer.m_parameters.minBaumWelchInverseScalar = 0;
        trainer.m_parameters.maxBaumWelchInverseScalar = 0;
        //trainer.m_parameters.maxPositionCycles = 10;
#endif //ALLOW_BOLTZMANN_GIBBS

        // TODO: REMOVE
        //trainer.m_parameters.debug = galosh::DEBUG_All;
        //trainer.m_parameters.verbosity = galosh::VERBOSITY_All;
        trainer.m_parameters.verbosity = galosh::VERBOSITY_Low;
      
        cout << "The profile (before) is:" << endl;
        cout << dna_profile << endl;
        score = trainer.train();
        cout << "Now (after training), the score is " << score << ", and the profile is:" << endl;
        cout << *trainer.m_profile << endl;
        // Without rabiner scaling (and no priors), I get:
        // The profile (before) is:
        // [ M->(M=0.97,I=0.02,D=0.01), I->(M=0.9,I=0.1), D->(M=0.8,D=0.2), I:(A=0.25,C=0.25,G=0.25,T=0.25), N->(N=0.1,B=0.9), B->(M=0.99,D=0.01), C->(C=0.1,T=0.9) ]
        // [ M:(A=0.1359667,C=0.3251242,G=0.2348561,T=0.3040529) ]
        // [ M:(A=8.000848e-05,C=0.08962184,G=0.2114788,T=0.6988194) ]
        // 
        // The total forward score before training is 1.483109e-08
        // Training..{Iteration#1[Positions:1#1+.+.score:0.0001223387;distance:0.0959217][Positions:1#2+.+.score:0.0001226138;distance:6.38681e-08][Globals:1#1..+.score:0.01538392;distance:0.136444][Globals:1#2..+.score:0.01557982;distance:0.000375724][Globals:1#3..+.score:0.015604;distance:1.07417e-07]distance:0.187433}{Iteration#2[Positions:2#1+.+.score:0.01561316;distance:1.07302e-08][Globals:2#1..+.score:0.01561923;distance:4.6678e-09]distance:1.36475e-08}done.
        // Done training.  Final score is 0.01561923.  Starting score was 1.483109e-08.
        // Now (after training), the score is 0.01561923, and the profile is:
        // [ M->(M=0.6666868,I=1e-05,D=0.3333032), I->(M=0.99999,I=1e-05), D->(M=1e-05,D=0.99999), I:(A=0.9999538,C=1e-05,G=2.61701e-05,T=1e-05), N->(N=1e-05,B=0.99999), B->(M=0.9999467,D=5.331392e-05), C->(C=0.250033,T=0.749967) ]
        // [ M:(A=0.99997,C=1e-05,G=1.000003e-05,T=1e-05) ]
        // [ M:(A=1e-05,C=1e-05,G=0.99997,T=1e-05) ]

        // For -DDISALLOW_FLANKING_TRANSITIONS
        // Without rabiner scaling (and no priors), I get:
        // The profile (before) is:
        // [ M->(M=0.97,I=0.02,D=0.01), I->(M=0.9,I=0.1), D->(M=0.8,D=0.2), I:(A=0.25,C=0.25,G=0.25,T=0.25), N->(N=0,B=1), B->(M=0.99,D=0.01), C->(C=0,T=1) ]
        // [ M:(A=0.1359667,C=0.3251242,G=0.2348561,T=0.3040529) ]
        // [ M:(A=8.000848e-05,C=0.08962184,G=0.2114788,T=0.6988194) ]
        // 
        // The total forward score before training is 1.802175e-12
        // Training..{Iteration#1[Positions:1#1+.+.score:1.486565e-05;distance:0.0891845][Positions:1#2+.+.score:1.5159e-05;distance:0.000563073][Positions:1#3+.-.score:1.516069e-05;distance:3.14723e-06][Globals:1#1..+.score:0.007852757;distance:0.15347][Globals:1#2..+.score:0.00827672;distance:0.000213343][Globals:1#3..+.score:0.008582542;distance:0.000107642][Globals:1#4..+.score:0.008784385;distance:4.60037e-05]distance:0.180136}{Iteration#2[Positions:2#1+.-.score:0.00890306;distance:0.000455815][Positions:2#2+.-.score:0.0089031;distance:1.54114e-07][Globals:2#1..+.score:0.009067931;distance:2.27764e-05][Globals:2#2..+.score:0.009157112;distance:6.61909e-06]distance:0.000503561}{Iteration#3[Positions:3#1+.-.score:0.009159103;distance:7.68023e-06][Globals:3#1..+.score:0.009207456;distance:1.88029e-06]distance:8.75468e-06}done.
        // Done training.  Final score is 0.009207456.  Starting score was 1.802175e-12.
        // Now (after training), the score is 0.009207456, and the profile is:
        // [ M->(M=0.3345434,I=0.3345434,D=0.3309132), I->(M=0.99999,I=1e-05), D->(M=0.99999,D=1e-05), I:(A=1e-05,C=1e-05,G=0.99997,T=1e-05), N->(N=0,B=1), B->(M=0.996383,D=0.003617042), C->(C=0,T=1) ]
        // [ M:(A=0.99997,C=1e-05,G=1e-05,T=1e-05) ]
        // [ M:(A=0.5027267,C=1e-05,G=0.4972533,T=1e-05) ]
        
        // For -DUSE_FLANKING_EMISSION_DISTRIBUTIONS:
        // The profile (before) is:
        // [ M->(M=0.97,I=0.02,D=0.01), I->(M=0.9,I=0.1), D->(M=0.8,D=0.2), I:(A=0.25,C=0.25,G=0.25,T=0.25), N->(N=0.1,B=0.9), N:(A=0.25,C=0.25,G=0.25,T=0.25), B->(M=0.99,D=0.01), C->(C=0.1,T=0.9), C:(A=0.25,C=0.25,G=0.25,T=0.25) ]
        // [ M:(A=0.1359667,C=0.3251242,G=0.2348561,T=0.3040529) ]
        // [ M:(A=8.000848e-05,C=0.08962184,G=0.2114788,T=0.6988194) ]
        //
        // Without rabiner scaling (and no priors), I get:
        // The total forward score before training is 1.483109e-08
        // Training..{Iteration#1[Positions:1#1+.+.score:0.0001223387;distance:0.0697613][Positions:1#2+.+.score:0.0001226138;distance:4.64496e-08][Globals:1#1..+.score:0.01538752;distance:0.165987][Globals:1#2..+.score:0.01558073;distance:0.000336033][Globals:1#3..+.score:0.01560423;distance:2.99469e-05][Globals:1#4..+.score:0.01561012;distance:4.5483e-09]distance:0.198815}{Iteration#2[Positions:2#1+.+.score:0.01561928;distance:7.80326e-09][Globals:2#1..+.score:0.01562076;distance:1.77951e-10]distance:7.93268e-09}done.
        // Done training.  Final score is 0.01562076.  Starting score was 1.483109e-08.
        // Now (after training), the score is 0.01562076, and the profile is:
        // [ M->(M=0.6666638,I=1e-05,D=0.3333262), I->(M=0.99999,I=1e-05), D->(M=1e-05,D=0.99999), I:(A=1e-05,C=1e-05,G=0.99997,T=1e-05), N->(N=1e-05,B=0.99999), N:(A=0.9999683,C=1e-05,G=1.165196e-05,T=1e-05), B->(M=0.9999866,D=1.333183e-05), C->(C=0.2500085,T=0.7499915), C:(A=0.9999647,C=1e-05,G=1.533362e-05,T=1e-05) ]
        // [ M:(A=0.99997,C=1e-05,G=1.000003e-05,T=1e-05) ]
        // [ M:(A=1.000909e-05,C=1e-05,G=0.99997,T=1e-05) ]

        // ==ambiguous== For -DUSE_FLANKING_EMISSION_DISTRIBUTIONS:
        // ==ambiguous== Without rabiner scaling (and no priors), when test_ambiguous is true, I get:
        // ==ambiguous== The total forward score before training is 2.056723e-05
        // ==ambiguous== Training..{Iteration#1[Positions:1#1+.+.score:0.002554536;distance:0.0376974][Positions:1#2-.-.score:0.002554536;distance:0][Globals:1#1..+.score:0.06126719;distance:0.0688172][Globals:1#2..+.score:0.06769019;distance:0.0253723][Globals:1#3..+.score:0.06866747;distance:0.00184733][Globals:1#4..+.score:0.06881911;distance:0.000241211]distance:0.0860183}{Iteration#2[Positions:2#1-.-.score:0.06881911;distance:0][Globals:2#1..+.score:0.0688425;distance:2.41402e-05][Globals:2#2..+.score:0.06884427;distance:2.23149e-06]distance:2.98419e-05}{Iteration#3[Positions:3#1-.-.score:0.06884427;distance:0][Globals:3#1..+.score:0.06884427;distance:2.02621e-07]distance:1.47361e-07}done.
        // ==ambiguous== Done training.  Final score is 0.06884427.  Starting score was 2.056723e-05.
        // ==ambiguous== Now (after training), the score is 0.06884427, and the profile is:
        // ==ambiguous== [ M->(M=0.6666523,I=1e-05,D=0.3333377), I->(M=0.5,I=0.5), D->(M=0.99999,D=1e-05), I:(A=0.25,C=0.25,G=0.25,T=0.25), N->(N=1e-05,B=0.99999), N:(A=0.953267,C=0.01557767,G=0.01557767,T=0.01557767), B->(M=0.99999,D=1e-05), C->(C=1e-05,T=0.99999), C:(A=0.125174,C=0.125174,G=0.6244779,T=0.125174) ]
        // ==ambiguous== [ M:(A=0.8625903,C=0.04577977,G=0.04585019,T=0.04577977) ]
        // ==ambiguous== [ M:(A=0.1252795,C=0.1250575,G=0.6246055,T=0.1250575) ]


        // TODO: Make this depend on a parameter in ProfileTrainer.hpp governing whether or not the forward matrices will be initialized.
        if( trainer.m_forward_matrices.size() != ( trainer.m_profile->length() + 1 ) ) {
          trainer.m_forward_matrices.reinitialize(
            trainer.m_profile->length(),
            trainer.m_sequences,
            trainer.m_sequence_count
          );
        } // End if the forward matrices are not initialized.
      
        score =
          dp.forward_score_viterbi(
            trainer.m_parameters,
            *trainer.m_profile,
            trainer.m_sequences,
            trainer.m_sequence_count,
            trainer.m_forward_matrices
          );
        cout << "The total score for all sequences, using viterbi, is: " << score << endl;
        //  cout << "The viterbi forward matrices are: " << endl;
        //  cout << forward_matrices << endl;
        DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::MultipleAlignment<ProfileType, Dna5> ma(
          trainer.m_profile,
          &( trainer.m_sequences ),
          trainer.m_sequence_count
        );
        dp.forward_viterbiAlign(
          trainer.m_parameters,
          trainer.m_forward_matrices,
          ma
        );
        // TODO: Put back.
        cout << "Viterbi score, calculated using ma.calculateScore(), is " << ma.calculateScore() << endl;
        cout << "Multiple Alignment is:" << endl;
        cout << ma << endl;
        ma.toPairwiseStream( cout );
        // Without rabiner scaling (and with no priors), I get:
        // The total score for all sequences, using viterbi, is: 0.0156185
        // Viterbi score, calculated using ma.calculateScore(), is 0.0156185
        // Multiple Alignment is:
        // Sequence 0:
        // ( 0, 1, 0, 0 )
        // Sequence 1:
        // ( 0, 1, 1, 0 )
        // Sequence 2:
        // ( 0, 1, 1, 1 )
        // 
        // > Sequence 0
        // 
        // Model    1 AG 2
        //            |v
        // Sequence 1 A- 1
        // 
        // 
        // 
        // > Sequence 1
        // 
        // Model    1 AG 2
        //            ||
        // Sequence 1 AG 2
        // 
        // 
        // 
        // > Sequence 2
        // 
        // Model    1 AG- 2
        //            ||^
        // Sequence 1 AGA 3

        // For -DDISALLOW_FLANKING_TRANSITIONS
        // Without rabiner scaling (and no priors), I get:
        // The total score for all sequences, using viterbi, is: 0.009156954
        // Viterbi score, calculated using ma.calculateScore(), is 9.157138e-08
        // Multiple Alignment is:
        // Sequence 0:
        // ( 0, 1, 0, 0 )
        // Sequence 1:
        // ( 0, 1, 1, 0 )
        // Sequence 2:
        // ( 0, 2, 1, 0 )
        // 
        // > Sequence 0
        // 
        // Model    1 AA 2
        //            |v
        // Sequence 1 A- 1
        // 
        // 
        // 
        // > Sequence 1
        // 
        // Model    1 AA 2
        //            | 
        // Sequence 1 AG 2
        // 
        // 
        // 
        // > Sequence 2
        // 
        // Model    1 A-A 2
        //            |^|
        // Sequence 1 AGA 3

        // For -DUSE_FLANKING_EMISSION_DISTRIBUTIONS:
        // Without rabiner scaling (and with no priors), I get:
        // The total score for all sequences, using viterbi, is: 0.01562054
        // Viterbi score, calculated using ma.calculateScore(), is 0.01562054
        // Multiple Alignment is:
        // Sequence 0:
        // ( 0, 1, 0, 0 )
        // Sequence 1:
        // ( 0, 1, 1, 0 )
        // Sequence 2:
        // ( 0, 1, 1, 1 )
        // 
        // > Sequence 0
        // 
        // Model    1 AG 2
        //            |v
        // Sequence 1 A- 1
        // 
        // 
        // 
        // > Sequence 1
        // 
        // Model    1 AG 2
        //            ||
        // Sequence 1 AG 2
        // 
        // 
        // 
        // > Sequence 2
        // 
        // Model    1 AG- 2
        //            ||^
        // Sequence 1 AGA 3

        // ==ambiguous== For -DUSE_FLANKING_EMISSION_DISTRIBUTIONS:
        // ==ambiguous== Without rabiner scaling (and no priors), when test_ambiguous is true, I get:
        // ==ambiguous== The total score for all sequences, using viterbi, is: 0.06884288
        // ==ambiguous== Viterbi score, calculated using ma.calculateScore(), is 0.06884288
        // ==ambiguous== Multiple Alignment is:
        // ==ambiguous== Sequence 0:
        // ==ambiguous== ( 0, 1, 0, 0 )
        // ==ambiguous== Sequence 1:
        // ==ambiguous== ( 0, 1, 1, 0 )
        // ==ambiguous== Sequence 2:
        // ==ambiguous== ( 0, 1, 1, 0 )
        // ==ambiguous== 
        // ==ambiguous== > Sequence 0
        // ==ambiguous== 
        // ==ambiguous== Model    1 AG 2
        // ==ambiguous==            |v
        // ==ambiguous== Sequence 1 A- 1
        // ==ambiguous== 
        // ==ambiguous== 
        // ==ambiguous== 
        // ==ambiguous== > Sequence 1
        // ==ambiguous== 
        // ==ambiguous== Model    1 AG 2
        // ==ambiguous==            | 
        // ==ambiguous== Sequence 1 AN 2
        // ==ambiguous== 
        // ==ambiguous== 
        // ==ambiguous== 
        // ==ambiguous== > Sequence 2
        // ==ambiguous== 
        // ==ambiguous== Model    1 AG 2
        // ==ambiguous==            ||
        // ==ambiguous== Sequence 1 NG 2

      } // End if test_trainer
    } // End if test_dynamic_programming || test_trainer
  } // End if test_profiles || test_dynamic_programming || test_trainer

  ///////
  // Note: is seeds something implementing BLAST?  Could Robert use it instead of wublast?
  return;
} // test_basics ()
#endif // TEST_BASICS

namespace galosh {
template <typename ProbabilityType,
          typename ScoreType,
          typename MatrixValueType>
class Tests {
public:
  template <typename ResidueType, typename SequenceResidueType>
  int
  test ( int const & argc, char const **argv, ResidueType const & residue_type_tag, SequenceResidueType const & sequence_residue_type_tag )
  {
    typedef ProfileTreeRoot<ResidueType, ProbabilityType> ProfileType;
    typedef ProfileTreeRoot<ResidueType, ProbabilityType> InternalNodeType;
  
    const uint32_t train_iterations = 1; //5;
    bool use_lengthadjust = false;//true;
    double lengthadjust_threshold = .5;//.025;
    double lengthadjust_threshold_increment = .0005;

    // Calculate viterbi score after training (root profile only)?
    const bool calculate_viterbi_score = true;
    // Calculate multiple alignment after training (root profile only)?
    const bool calculate_multiple_alignment = false;//true;

    // The "TREE TRAINER" is "PROLIFIC" iff MUSCLE is compiled in (see Makefile).  Otherwise, it's just a regular ProfileTrainer.
    // Use this to switch between the tree trainer and the regular trainer: comment this out to use the nontree trainer.
    #define USE_TREE_TRAINER
    // Comment this out to use the nontree trainer after the tree trainer (only if USE_TREE_TRAINER is not commented out)
    //#define USE_NONTREE_TRAINER_TOO

#ifdef USE_TREE_TRAINER
# define PROFILE_TRAINER tree_trainer
# define PROFILE_TRAINER_IS_TREE_TRAINER true
# define ROOT_PROFILE_PTR ( tree_trainer.m_profileTree->getProfileTreeRoot() )
#else
# define PROFILE_TRAINER nontree_trainer
# define PROFILE_TRAINER_IS_TREE_TRAINER false
# define ROOT_PROFILE_PTR ( nontree_trainer.m_profile )
# undef USE_NONTREE_TRAINER_TOO
#endif // USE_TREE_TRAINER

    Random random = Random();
    DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType> dp;
  
    ScoreType score;
  
    unsigned long seed = std::time( NULL );
    // TODO: REMOVE!  TESTING.  This used to cause us a problem.
    //seed = 1220097217L;
    // seed = 1227122479L;
    //seed = 1227908196L;
    seed = 1303806565L;
    cout << "NOTE: RANDOM SEED IS " << seed << endl;
    random.setSeed( seed );
  
    Fasta<SequenceResidueType> fasta;
     
    int initial_profile_length = 0; // 0 means use the median seq length, 1 means use max
    bool have_true_profile = false;//true;
    ProfileType profile_true;
    //profile_true.fromFile( "test_results/v8_seed1303599137_type/ProfuseTest.v8_seed1303599137_type.90.100.100.100.1.1.0.10000000000000001.0.10000000000000001.1.true.ProfileTree.prof" );
    //profile_true.fromFile( "test_results/v8_seed1303717895_type/ProfuseTest.v8_seed1303717895_type.50.100.100.100.1.1.0.10000000000000001.0.10000000000000001.0.true.ProfileTree.prof" );
    //profile_true.fromFile( "test_results/v8_seed1303842002_type/ProfuseTest.v8_seed1303842002_type.50.100.100.100.1.1.0.10000000000000001.0.10000000000000001.0.true.ProfileTree.prof" );
    //profile_true.fromFile( "test_results/v8_seed1303849236_type/ProfuseTest.v8_seed1303849236_type.30.100.100.100.1.1.0.10000000000000001.0.10000000000000001.0.true.ProfileTree.prof" );
    bool have_starting_profile = false;
    ProfileType profile_starting;
    // TODO: REMOVE
//#ifdef DISALLOW_FLANKING_TRANSITIONS
//    have_starting_profile = true;
//#ifdef USE_DEL_IN_DEL_OUT
//#ifdef USE_SWENTRY_SWEXIT
//    profile_starting.fromFile( "deleteme.swee.prof" );
//#else
//    profile_starting.fromFile( "deleteme.prof" );
//    //profile_starting.fromFile( "deleteme.prof.shortened.3" );
//#endif // USE_SWENTRY_SWEXIT .. else ..
//#else // USE_DEL_IN_DEL_OUT .. else ..
//    profile_starting.fromFile( "deleteme.noflank.prof" );
//#endif // USE_DEL_IN_DEL_OUT
//#endif // DISALLOW_FLANKING_TRANSITIONS
    have_starting_profile = false;//true;
    //profile_starting.fromFile( "test_results_hiatus_DNA/v9_seed1304287230_type/ProfuseTest.v9_seed1304287230_type.50.100.100.100.1.1.0.10000000000000001.0.10000000000000001.1.1.starting.ProfileTree.prof" );
    //profile_starting.fromFile( "test_results/v8_seed1303599137_type/ProfuseTest.v8_seed1303599137_type.90.100.100.100.1.1.0.10000000000000001.0.10000000000000001.1.2.starting.ProfileTree.prof" );
    //profile_starting.fromFile( "test_results/v8_seed1303717895_type/ProfuseTest.v8_seed1303717895_type.50.100.100.100.1.1.0.10000000000000001.0.10000000000000001.0.2.starting.ProfileTree.prof" );
    //profile_starting.fromFile( "test_results/v8_seed1303842002_type/ProfuseTest.v8_seed1303842002_type.50.100.100.100.1.1.0.10000000000000001.0.10000000000000001.0.1.starting.ProfileTree.prof" );
    //profile_starting.fromFile( "test_results/v8_seed1303849236_type/ProfuseTest.v8_seed1303849236_type.30.100.100.100.1.1.0.10000000000000001.0.10000000000000001.0.0.starting.ProfileTree.prof" );
    bool have_trained_profile = false;
    bool start_with_trained_profile = false;
    ProfileType profile_trained;
    // TODO: REMOVE
    //profile_trained.fromFile( "/Users/pedlefse/src/prefab4/in/1a02N_1a3qA.profuse.CQA_LA-max.train.prof" );

    // TODO: REMOVE.  ERE I AM.  TESTING.
    ProfileType other_starting_profile;
    //other_starting_profile.fromFile( "test_results/v8_seed1303599137_type/ProfuseTest.v8_seed1303599137_type.90.100.100.100.1.1.0.10000000000000001.0.10000000000000001.1.0.starting.ProfileTree.prof" );
    //other_starting_profile.fromFile( "test_results/v8_seed1303717895_type/ProfuseTest.v8_seed1303717895_type.50.100.100.100.1.1.0.10000000000000001.0.10000000000000001.0.3.starting.ProfileTree.prof" );
    // This is the test: 
    //profile_starting.copyExceptPositions( other_starting_profile );
    // MORE TESTING
     
//                                ProbabilityType average_indel_open = 1.0;
//                                average_indel_open -= profile_starting[ Transition::fromMatch ][ TransitionFromMatch::toMatch ];
//                                average_indel_open /= 2;
//                                profile_starting[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] =
//                                  average_indel_open;
//                                profile_starting[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] =
//                                  average_indel_open;
//                                // TODO: REMOVE?!
//                                profile_starting[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] = average_indel_open;
//                                profile_starting[ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ] = ( 1.0 - average_indel_open );
//                                profile_starting[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] = average_indel_open;
//                                profile_starting[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] = ( 1.0 - average_indel_open );

    if( argc > 1 ) {
      // At least one argument: fasta filename
      fasta.fromFile( argv[ 1 ] );
      if( argc > 2 ) {
        // At least 2 arguments: 2nd is initial_profile_length
        try {
          initial_profile_length = boost::lexical_cast<int>( argv[ 2 ] );
        } catch( boost::bad_lexical_cast & ) {
          std::cerr << "Unable to interpret the argument '" << argv[ 2 ] << "' as an integer for use as the " << ( use_lengthadjust ? "initial " : "" ) << "profile length." << std::endl;
          exit( 1 );
        } // End try .. catch block for lexical_cast
        if( argc > 3 ) {
          // At least 3 arguments: 3rd is lengthadjust threshold (or 0 for no lenthadjust)
          try {
            lengthadjust_threshold = boost::lexical_cast<double>( argv[ 3 ] );
          } catch( boost::bad_lexical_cast & ) {
            std::cerr << "Unable to interpret the argument '" << argv[ 3 ] << "' as a real value for use as the lengthadjust threshold (or as 0 to indicate no lengthadjust)." << std::endl;
            exit( 1 );
          } // End try .. catch block for lexical_cast
          if( lengthadjust_threshold < 0 ) {
            std::cerr << "The given lengthadjust threshold value, " << lengthadjust_threshold << ", is negative.  You must supply a value between 0 and 1 (or 0 to indicate no lengthadjust)." << std::endl;
            exit( 1 );
          }
          if( lengthadjust_threshold > 1 ) {
            std::cerr << "The given lengthadjust threshold value, " << lengthadjust_threshold << ", is greater than 1.  You must supply a value between 0 and 1 (or 0 to indicate no lengthadjust)." << std::endl;
            exit( 1 );
          }
          if( lengthadjust_threshold == 0 ) {
            std::cout << "Lengthadjust disabled." << endl;
            use_lengthadjust = false;
          } else {
            use_lengthadjust = true;
          }
          if( use_lengthadjust && ( argc > 4 ) ) {
            // At least 4 arguments: 4th is lengthadjust threshold increment
            try {
              lengthadjust_threshold_increment = boost::lexical_cast<double>( argv[ 4 ] );
            } catch( boost::bad_lexical_cast & ) {
              std::cerr << "Unable to interpret the argument '" << argv[ 4 ] << "' as a real value for use as the lengthadjust threshold increment." << std::endl;
              exit( 1 );
            } // End try .. catch block for lexical_cast
            if( lengthadjust_threshold_increment <= 0 ) {
              std::cerr << "The given lengthadjust threshold increment value, " << lengthadjust_threshold << ", is negative.  You must supply a value between 0 and 1." << std::endl;
              exit( 1 );
            }
            if( lengthadjust_threshold_increment > 1 ) {
              std::cerr << "The given lengthadjust threshold increment value, " << lengthadjust_threshold << ", is greater than 1.  You must supply a value between 0 and 1." << std::endl;
              exit( 1 );
            }
            if( lengthadjust_threshold_increment == 0 ) {
              std::cerr << "WARNING: Lengthadjust threshold incrementing disabled." << endl;
            }
          } // At least 4 arguments: 4th is lengthadjust threshold increment
        } // At least 3 arguments: 3rd is lengthadjust threshold (or 0 for no lenthadjust)
      } // At least 2 arguments: 2nd is initial_profile_length

      // At least one argument: fasta filename
    } else { // No arguments, use defaults
      int default_initial_profile_length;
      getDefaultFasta(
        fasta,
        default_initial_profile_length
      );
      if( initial_profile_length == 0 ) {
        initial_profile_length = default_initial_profile_length;
      }
    }
    if( fasta.size() == 0 ) {
      cout << "Fasta file contains no sequences!" << endl;
      return 1; // TODO: RETURN ERROR CODES ?
    }

    uint32_t num_sequences_to_use = fasta.size();
    // TODO: REMOVE
    //num_sequences_to_use = 50; // Keep it simple
  
    // Set up the profile with some reasonable initial global values
    ProfileType profile;

    if( have_trained_profile && start_with_trained_profile && ( profile_trained.length() > 0 ) ) {
      profile = profile_trained;
      initial_profile_length = profile.length();
    } else if( have_starting_profile && ( profile_starting.length() > 0  ) ) {
      profile = profile_starting;
      initial_profile_length = profile.length();
#ifndef DISALLOW_FLANKING_TRANSITIONS
      // TODO: REMOVE
      //profile[ Transition::fromPreAlign ][ Transition::FromPreAlign::toPreAlign ] = .2;
      //profile[ Transition::fromPreAlign ][ Transition::FromPreAlign::toBegin ] = .8;
      //profile[ Transition::fromPostAlign ][ Transition::FromPostAlign::toPostAlign ] = .2;
      //profile[ Transition::fromPostAlign ][ Transition::FromPostAlign::toTerminal ] = .8;
#endif // !DISALLOW_FLANKING_TRANSITIONS
  
    } else { // if ( have_trained_profile && start_with_trained_profile ) .. else if( have_starting_profile ) .. else ..
    
      if( initial_profile_length == 1 ) {
        // Guess the length.  Use the maximum length of the sequences.
        initial_profile_length = 0;
        uint32_t seq_length;
        for( uint32_t seq_i = 0; seq_i < num_sequences_to_use; seq_i++ ) {
          seq_length = fasta[ seq_i ].length();
          if( seq_length > initial_profile_length ) {
            initial_profile_length = seq_length;
          }
        } // End foreach seq, find the one with the max length
      } else if( initial_profile_length == 0 ) {
        // Guess the length.  Use the median length of the sequences.
        vector<uint32_t> lengths( num_sequences_to_use );
        for( uint32_t seq_i = 0; seq_i < num_sequences_to_use; seq_i++ ) {
          lengths[ seq_i ] = fasta[ seq_i ].length();
        }
        // Sort them.
        std::sort( lengths.begin(), lengths.end() );
        if( ( num_sequences_to_use % 2 ) == 1 ) {
          // Odd number, use the middle value.
          initial_profile_length = lengths[ ( num_sequences_to_use - 1 ) / 2 ];
        } else {
          // Even number, use the average of the two middle values.
          initial_profile_length = lengths[ ( num_sequences_to_use / 2 ) - 1 ];
          initial_profile_length += lengths[ ( num_sequences_to_use / 2 ) ];
          // Round it, rather than truncate it:
          initial_profile_length =
            uint32_t( ( ( float )initial_profile_length / 2.0f ) + .5f );
        }
      } // End if initial_profile_length == 0
      if( initial_profile_length == 1 ) { // We can't handle length 1 profiles.
        initial_profile_length = 2;
      }
      assert( initial_profile_length > 1 );

      profile.reinitialize( ( uint32_t )initial_profile_length );
      // Start from a uniform profile -- except for the globals..
      //profile.uniformPositions( random );
      // NEW WAY: start nearer to even...
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DirichletMixtureMatchEmissionPrior<float> a_matchEmissionPrior;
      a_matchEmissionPrior.reinitializeToEven( 100 );
      cout << "a_matchEmissionPrior: " << a_matchEmissionPrior << endl;
      // Draw position params from the prior
      profile.dirichletMixturePositions(
        a_matchEmissionPrior,
        random
      );
 
      profile.normalizePositions( 1E-5 ); // Ensure values aren't too tiny.

#ifndef DISALLOW_FLANKING_TRANSITIONS
      profile[ galosh::Transition::fromPreAlign ][ galosh::TransitionFromPreAlign::toPreAlign ] =
        .01;
      profile[ galosh::Transition::fromPreAlign ][ galosh::TransitionFromPreAlign::toBegin ] =
        1.0 -
        profile[ galosh::Transition::fromPreAlign ][ galosh::TransitionFromPreAlign::toPreAlign ];
#endif // !DISALLOW_FLANKING_TRANSITIONS
      profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletion ] =
        .01;
#ifdef USE_DEL_IN_DEL_OUT
      profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletionIn ] =
        .5;
      profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toMatch ] =
        1.0 -
        (
          profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletion ] +
          profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletionIn ]
        );
#ifndef USE_SWENTRY_SWEXIT
      profile[ galosh::Transition::fromDeletionIn ][ galosh::TransitionFromDeletionIn::toDeletionIn ] =
        .99;//.5;
      profile[ galosh::Transition::fromDeletionIn ][ galosh::TransitionFromDeletionIn::toMatch ] =
        1.0 -
        profile[ galosh::Transition::fromDeletionIn ][ galosh::TransitionFromDeletionIn::toDeletionIn ];
#endif // !USE_SWENTRY_SWEXIT
#else
      profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toMatch ] =
        1.0 -
        profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletion ];
#endif // USE_DEL_IN_DEL_OUT .. else ..
      
      profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toInsertion ] =
        .01;
      profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletion ] =
        .01;
  #ifdef USE_DEL_IN_DEL_OUT
      profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletionOut ] =
        .5;
      profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toMatch ] =
        1.0 -
        (
          profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toInsertion ] +
          profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletion ] +
          profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletionOut ]
        );
#ifndef USE_SWENTRY_SWEXIT
      profile[ galosh::Transition::fromDeletionOut ][ galosh::TransitionFromDeletionOut::toDeletionOut ] =
        .99;
      profile[ galosh::Transition::fromDeletionOut ][ galosh::TransitionFromDeletionOut::toEnd ] =
        1.0 -
        profile[ galosh::Transition::fromDeletionOut ][ galosh::TransitionFromDeletionOut::toDeletionOut ];
#endif // !USE_SWENTRY_SWEXIT
#else
      profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toMatch ] =
        1.0 -
        (
          profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toInsertion ] +
          profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletion ]
        );
#endif // USE_DEL_IN_DEL_OUT .. else ..
      profile[ 0 ][ galosh::Transition::fromInsertion ][ galosh::TransitionFromInsertion::toInsertion ] =
        .5;
      profile[ 0 ][ galosh::Transition::fromInsertion ][ galosh::TransitionFromInsertion::toMatch ] =
        1.0 -
        profile[ 0 ][ galosh::Transition::fromInsertion ][ galosh::TransitionFromInsertion::toInsertion ];
      profile[ 0 ][ galosh::Transition::fromDeletion ][ galosh::TransitionFromDeletion::toDeletion ] =
        .5;
      profile[ 0 ][ galosh::Transition::fromDeletion ][ galosh::TransitionFromDeletion::toMatch ] =
        1.0 -
        profile[ 0 ][ galosh::Transition::fromDeletion ][ galosh::TransitionFromDeletion::toDeletion ];
#ifndef DISALLOW_FLANKING_TRANSITIONS
      profile[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toPostAlign ] =
        .01;
      profile[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toTerminal ] =
        1.0 -
        profile[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toPostAlign ];
#endif // !DISALLOW_FLANKING_TRANSITIONS
    } // End if have_trained_profile or have_starting_profile .. else ..

    if( use_lengthadjust ) {
      cout << "Initial profile length: " << initial_profile_length << endl;
      cout << "Lengthadjust threshold: " << lengthadjust_threshold << endl;
      cout << "Lengthadjust threshold increment: " << lengthadjust_threshold_increment << endl;
    } else {
      cout << "Profile length: " << initial_profile_length << endl;
    }

    // TODO: REMOVE
    //num_sequences_to_use = 25;

    //cout << "FASTA from file:" << endl;
    //cout << fasta << endl;

#ifdef USE_TREE_TRAINER
    ProfileTree<ResidueType, ProbabilityType, InternalNodeType> profile_tree( initial_profile_length );
    profile_tree.getProfileTreeRoot()->copyFrom( profile );

    ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType> tree_trainer(
      fasta,
      num_sequences_to_use,
      &profile_tree
    );
#endif // USE_TREE_TRAINER
  
    // Used if PROFILE_TRAINER is nontree_trainer (in which case PROFILE_TRAINER_IS_TREE_TRAINER should be false)
    ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType> nontree_trainer( &profile, fasta, num_sequences_to_use );
  
    // To match ProfuseTest.cpp params, we set them as below:
    // We need this to get the right (default) params for the priorMtoM, etc.
    typename ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters training_parameters_template;
    {
      ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType> profuse_test;
      training_parameters_template =
        profuse_test.m_parameters;
    }

    // TODO: REMOVE!
    //training_parameters_template.profileValueMinimum = 0;

        
    // TODO: Try with the lengthadjust: maxBaumWelchInverseScalar > 0
    training_parameters_template.minIterations = 1;
    training_parameters_template.maxIterations = 1000;
    // Having maxPositionCycles_globals > 1 seems ok; takes about the same
    // number of iterations, converges to roughly the same place; takes
    // longer by virtue of having more pos cycles per iteration of course.
    training_parameters_template.maxPositionCycles = 1;
    // Having maxPositionCycles_globals > 1 seems to make convergence way
    // slower when lengthadjust is on.  Length keeps adjusting..
    training_parameters_template.maxPositionCycles_globals = 1;
    training_parameters_template.minBaumWelchInverseScalar = 0; // Straight-up bw.
    training_parameters_template.maxBaumWelchInverseScalar = 0; // Straight-up bw.
    training_parameters_template.minBaumWelchInverseScalar_globals = 0; // Straight-up bw.
    training_parameters_template.maxBaumWelchInverseScalar_globals = 0; // Straight-up bw.
    training_parameters_template.scorePercentChangeMinimum_position_cycle = 1;
    training_parameters_template.scorePercentChangeMinimum_iteration = .01;
    training_parameters_template.usePriors = false;//true;
    // TODO: Debug alwaysAccept..?  It doesn't *always*.. ok?
    training_parameters_template.alwaysAccept = false;//true;
  
    training_parameters_template.proposeProfileLengthChanges = use_lengthadjust;
    training_parameters_template.useAlignmentProfiles = true;
    training_parameters_template.numIterationsBetweenLengthChanges = 0;
    training_parameters_template.proposeDeletingThreshold =
      lengthadjust_threshold; //.01;//.025;//.1;
    training_parameters_template.proposeDeletingThreshold_increment =
      lengthadjust_threshold_increment; //.0005;//.00005;//.0005;//5E-5;//.0003125;//.00625;//.025;
    training_parameters_template.proposeInsertingThreshold =
      training_parameters_template.proposeDeletingThreshold;// / 4;//seqan::ValueSize<ResidueType>::VALUE; // TODO: Figure this out...
    training_parameters_template.proposeInsertingPreAlignThreshold = //.35; //.5;
      training_parameters_template.proposeInsertingThreshold;
    training_parameters_template.proposeInsertingPostAlignThreshold = //.35;//.5;
      training_parameters_template.proposeInsertingThreshold;
    training_parameters_template.proposeInsertingThreshold_increment =
      training_parameters_template.proposeDeletingThreshold_increment;// TODO: ERE I AM! PUT BACK? / seqan::ValueSize<ResidueType>::VALUE;
  
    if( have_trained_profile && start_with_trained_profile ) {
      // When we start with the trained profile, we need to get past the
      // length modification wait time (which is just one iteration).
      training_parameters_template.minIterations =
        max( training_parameters_template.minIterations, ( uint32_t )2 );
    }
  
    // Use rabiner scaling? (default true)
    // NOTE: You must change the MatrixValueType to logspace or bfloat iff this is false!\
    // TODO: I don't think I've tested this in a long long time.  Probably safest to disable it for now.  At least make the default false.
    training_parameters_template.useRabinerScaling = false;

    // Train globals first?
    //training_parameters_template.trainGlobalsFirst = true; // note: breaks UBW.

    // Use Ubw?
    training_parameters_template.useUnconditionalBaumWelch = true;
    training_parameters_template.unconditionalIsolatesGlobals = false;
  
    // Use Baldi?
    // For testing Baldi-style gradient ascent
#ifdef ALLOW_BOLTZMANN_GIBBS
    // NOTE about priors:  since globals are presently not updated using Baldi, you can still usePriors and they will affect the globals *but not the positions*.
    training_parameters_template.baldiLearningRate = 1; // 0 means noBaldi!
    training_parameters_template.baldiTemperature = 1;
    training_parameters_template.baldiHybrid = false;
    training_parameters_template.siegelMaxFindingThePeakAttempts_positions = 1000; // 0 means Baldi not Baldi / Siegel !!!
    training_parameters_template.siegelEpsilonScaleFactor = 1.5;
    training_parameters_template.siegelMaxRefiningThePeakSteps_positions = 1000;
    training_parameters_template.siegelRefiningThePeakStepsConvergenceThreshold = 1E-5;
    training_parameters_template.minBaumWelchInverseScalar = 0;
    training_parameters_template.maxBaumWelchInverseScalar = 0;
    //training_parameters_template.maxPositionCycles = 10;
#endif //ALLOW_BOLTZMANN_GIBBS
  
    // This is copied from ProfuseTest.hpp, to test with the same priors used there.
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DirichletMixtureMatchEmissionPrior<float> matchEmissionPrior;
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DirichletMixtureGlobalPrior<float> globalPrior;
    if( training_parameters_template.usePriors ) {
      double priorStrength = 1;
      double priorStrength_internal_transitions = 1;
      double priorStrength_flanking_self_transitions = 1;
      double priorStrength_flanking_other_transitions = 1;
      // Here's how it is set up in ProfuseTest.hpp:
      matchEmissionPrior.reinitializeToEven( priorStrength );
      globalPrior.reinitializeToEven( priorStrength );
      training_parameters_template.matchEmissionPrior = &matchEmissionPrior;
      training_parameters_template.globalPrior = &globalPrior;
      // Put a strong prior on a small C->C and N->N
#ifndef DISALLOW_FLANKING_TRANSITIONS
      globalPrior[ 0 ][ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] = ( priorStrength_flanking_self_transitions * .01 );
      globalPrior[ 0 ][ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ] = ( priorStrength_flanking_other_transitions * .99 );
      globalPrior[ 0 ][ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] = ( priorStrength_flanking_self_transitions * .01 );
      globalPrior[ 0 ][ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] = ( priorStrength_flanking_other_transitions * .99 );
#endif // !DISALLOW_FLANKING_TRANSITIONS

      // Do additional setup of the transition priors for those transtions
      // observed many times (since the number of transitions depends on
      // the profile length).
      globalPrior[ 0 ][ Transition::fromMatch ][ TransitionFromMatch::toMatch ] = ( initial_profile_length * priorStrength_internal_transitions * training_parameters_template.priorMtoM );
      globalPrior[ 0 ][ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] = ( initial_profile_length * priorStrength_internal_transitions * training_parameters_template.priorMtoD );
#ifdef USE_DEL_IN_DEL_OUT
      // TODO: Create training_parameters_template.priorMtoW
      globalPrior[ 0 ][ Transition::fromMatch ][ TransitionFromMatch::toDeletionOut ] = ( initial_profile_length * priorStrength_internal_transitions * .5 );
#endif // USE_DEL_IN_DEL_OUT
      globalPrior[ 0 ][ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] = ( initial_profile_length * priorStrength_internal_transitions * training_parameters_template.priorMtoI );
#ifdef USE_DEL_IN_DEL_OUT
#ifndef USE_SWENTRY_SWEXIT
      // TODO: Create training_parameters_template.priorWtoW
      globalPrior[ 0 ][ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toDeletionOut ] = ( priorStrength_flanking_other_transitions * .99 );
      // TODO: Create training_parameters_template.priorWtoE
      globalPrior[ 0 ][ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toEnd ] = ( priorStrength_flanking_other_transitions * .01 );
#endif // !USE_SWENTRY_SWEXIT
      globalPrior[ 0 ][ Transition::fromBegin ][ TransitionFromBegin::toMatch ] = ( priorStrength_flanking_other_transitions * .49 );
      globalPrior[ 0 ][ Transition::fromBegin ][ TransitionFromBegin::toDeletion ] = ( priorStrength_flanking_other_transitions * .01 );
      // TODO: Create training_parameters_template.priorBtoZ
      globalPrior[ 0 ][ Transition::fromBegin ][ TransitionFromBegin::toDeletionIn ] = ( priorStrength_flanking_other_transitions * .50 );
      // TODO: Create training_parameters_template.priorZtoZ
#ifndef USE_SWENTRY_SWEXIT
      globalPrior[ 0 ][ Transition::fromDeletionIn ][ TransitionFromDeletionIn::toDeletionIn ] = ( priorStrength_flanking_other_transitions * .99 );
      // TODO: Create training_parameters_template.priorZtoM
      globalPrior[ 0 ][ Transition::fromDeletionIn ][ TransitionFromDeletionIn::toMatch ] = ( priorStrength_flanking_other_transitions * .01 );
#endif // !USE_SWENTRY_SWEXIT
#endif // USE_DEL_IN_DEL_OUT
      globalPrior[ 0 ][ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ] = ( initial_profile_length * priorStrength_internal_transitions * training_parameters_template.priorItoM );
      globalPrior[ 0 ][ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] = ( initial_profile_length * priorStrength_internal_transitions * training_parameters_template.priorItoI );
      globalPrior[ 0 ][ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ] = ( initial_profile_length * priorStrength_internal_transitions * training_parameters_template.priorDtoM );
      globalPrior[ 0 ][ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ] = ( initial_profile_length * priorStrength_internal_transitions * training_parameters_template.priorDtoD );
    } // End if usePriors
  
    //training_parameters_template.verbosity = VERBOSITY_All;
    training_parameters_template.verbosity = VERBOSITY_Low;
    //training_parameters_template.verbosity = VERBOSITY_None;
    //training_parameters_template.debug = DEBUG_All;
      
    //training_parameters_template.minBaumWelchInverseScalar = (1.0/4.0);
      
    // TODO: REMOVE.  Testing matrixRowScaleFactor
    //training_parameters_template.matrixRowScaleFactor =
    //  pow( numeric_limits<double>::max(), .25 ) - 1.0;
  
    PROFILE_TRAINER.m_parameters = training_parameters_template;
    // TODO: REMOVE
    //cout << "just before training, trainer parameters are " << PROFILE_TRAINER.m_parameters << endl;
  
    if( true ) {
      cout << "The profile (before) is:" << endl;
      cout << *ROOT_PROFILE_PTR << endl;

    } // End if print the before-profile
    // Do the training
#ifdef USE_TREE_TRAINER
    cout << "Training using the tree trainer." << endl;
    score = tree_trainer.train();
    uint32_t total_iterations =
      tree_trainer.m_totalIterations;

# ifdef USE_NONTREE_TRAINER_TOO
    profile.copyFrom( *profile_tree.getProfileTreeRoot() );
    // The profile is just a bunch of counts right now..
    tree_trainer.normalizeProfileOfCounts( profile );
            
    nontree_trainer.reinitialize(
      &profile
    );
    // Train globals first after first training..
    // TODO: REMOVE
    training_parameters_template.trainGlobalsFirst = true;
    // Train for at least 2 iterations, to make sure we get a chance to
    // do the lengthadjust.
    training_parameters_template.minIterations = 2;

    // start with higher propose*Threshold_* after first training..
    // TODO: REMOVE
    {
      training_parameters_template.proposeDeletingThreshold =
        .5; // TODO: MAGIC #
      training_parameters_template.proposeInsertingThreshold =
        training_parameters_template.proposeDeletingThreshold;
      training_parameters_template.proposeInsertingPreAlignThreshold =
        training_parameters_template.proposeInsertingThreshold;
      training_parameters_template.proposeInsertingPostAlignThreshold =
        training_parameters_template.proposeInsertingThreshold;
      training_parameters_template.proposeDeletingThreshold_increment =
        .005; // TODO: MAGIC #
      training_parameters_template.proposeInsertingThreshold_increment =
        training_parameters_template.proposeDeletingThreshold_increment;
    } // End TODO: REMOVE

    nontree_trainer.m_parameters = training_parameters_template;
    if( false ) {
      cout << "The profile (before retraining using non-tree trainer) is:" << endl;
      cout << profile << endl;
    } // End if print the before-profile

    // Note that this trains profile, not
    // *profile_tree.getProfileTreeRoot().
    score = nontree_trainer.train();
    total_iterations += nontree_trainer.m_iteration;
    cout << "Now (after training once), the score is " << score << ", and the profile length is " << nontree_trainer.m_profile->length() << "." << endl;
# else
    cout << "Now (after training once), the score is " << score << ", and the profile length is " << tree_trainer.m_profileTree->getProfileTreeRoot()->length() << "." << endl;
# endif // USE_NONTREE_TRAINER_TOO .. else ..

    for( uint32_t train_i = 2; train_i <= train_iterations; train_i++ ) {
      // Do it again
# ifndef USE_NONTREE_TRAINER_TOO
      profile.copyFrom( *profile_tree.getProfileTreeRoot() );
      // The profile is just a bunch of counts right now..
      tree_trainer.normalizeProfileOfCounts( profile );
# endif

      profile_tree.reinitialize( initial_profile_length );
      profile_tree.getProfileTreeRoot()->copyFrom( profile );
            
      tree_trainer.reinitialize(
        &profile_tree
      );
      // Train globals first after first training..
      // TODO: REMOVE
      training_parameters_template.trainGlobalsFirst = true;
      // Train for at least 2 iterations, to make sure we get a chance to
      // do the lengthadjust.
      training_parameters_template.minIterations = 2;

      // start with higher propose*Threshold_* after first training..
      // TODO: REMOVE
      {
        training_parameters_template.proposeDeletingThreshold =
          .5; // TODO: MAGIC #
        training_parameters_template.proposeInsertingThreshold =
          training_parameters_template.proposeDeletingThreshold;
        training_parameters_template.proposeInsertingPreAlignThreshold =
          training_parameters_template.proposeInsertingThreshold;
        training_parameters_template.proposeInsertingPostAlignThreshold =
          training_parameters_template.proposeInsertingThreshold;
        training_parameters_template.proposeDeletingThreshold_increment =
          .005; // TODO: MAGIC #
        training_parameters_template.proposeInsertingThreshold_increment =
          training_parameters_template.proposeDeletingThreshold_increment;
      } // End TODO: REMOVE
 
     tree_trainer.m_parameters = training_parameters_template;
      if( true ) {
        cout << "The profile (after training " << train_i << " times) is:" << endl;
        cout << *( tree_trainer.m_profileTree->getProfileTreeRoot() ) << endl;
      } // End if print the before-profile
      score = tree_trainer.train();
      total_iterations += tree_trainer.m_totalIterations;

# ifdef USE_NONTREE_TRAINER_TOO
      profile.copyFrom( *profile_tree.getProfileTreeRoot() );
      // The profile is just a bunch of counts right now..
      tree_trainer.normalizeProfileOfCounts( profile );
            
      nontree_trainer.reinitialize(
        &profile
      );
      nontree_trainer.m_parameters = training_parameters_template;
      if( false ) {
        cout << "The profile (before retraining using non-tree trainer) is:" << endl;
        cout << profile << endl;
      } // End if print the before-profile

      // Note that this trains profile, not
      // *profile_tree.getProfileTreeRoot().
      score = nontree_trainer.train();
      total_iterations += nontree_trainer.m_iteration;
# else
      cout << "Now (after training " << train_i << " times), the score is " << score << ", and the profile length is " << nontree_trainer.m_profile->length() << "." << endl;
# endif // USE_NONTREE_TRAINER_TOO
      cout << "Now (after training " << train_i << " times), the score is " << score << ", and the profile length is " << tree_trainer.m_profileTree->getProfileTreeRoot()->length() << "." << endl;
    } // End for each train_i <= train_iterations

    // Retrain again with lower length adjust thresholds (?)
# ifdef USE_NONTREE_TRAINER_TOO
    // TODO: REMOVE.  TESTING.
    if( false ) {
      // train one last time, with more agressive length adjusting.
      // Already copied and trained.
      //profile.copyFrom( *profile_tree.getProfileTreeRoot() );
      // The profile is just a bunch of counts right now..
      //tree_trainer.normalizeProfileOfCounts( profile );
            
      nontree_trainer.reinitialize(
        &profile
      );
      // Train globals first after first training..
      // TODO: REMOVE
      training_parameters_template.trainGlobalsFirst = true;
      // Train for at least 2 iterations, to make sure we get a chance to
      // do the lengthadjust.
      training_parameters_template.minIterations = 2;

      // redo it with lower propose*Threshold ...
      // TODO: REMOVE
      {
        training_parameters_template.proposeDeletingThreshold =
          .25; // TODO: MAGIC #
        training_parameters_template.proposeInsertingThreshold =
          training_parameters_template.proposeDeletingThreshold;
        training_parameters_template.proposeInsertingPreAlignThreshold =
          training_parameters_template.proposeInsertingThreshold;
        training_parameters_template.proposeInsertingPostAlignThreshold =
          training_parameters_template.proposeInsertingThreshold;
        training_parameters_template.proposeDeletingThreshold_increment =
          .005; // TODO: MAGIC #
        training_parameters_template.proposeInsertingThreshold_increment =
          training_parameters_template.proposeDeletingThreshold_increment;
      } // End TODO: REMOVE
      nontree_trainer.m_parameters = training_parameters_template;
      if( false ) {
        cout << "The profile (before retraining using non-tree trainer) is:" << endl;
        cout << profile << endl;
      } // End if print the before-profile

      // Note that this trains profile, not
      // *profile_tree.getProfileTreeRoot().
      score = nontree_trainer.train();
      total_iterations += nontree_trainer.m_iteration;
    }
# endif // USE_NONTREE_TRAINER_TOO

    // The profile is just a bunch of counts right now..
    tree_trainer.normalizeProfileOfCounts(
      *profile_tree.getProfileTreeRoot()
    );

    tree_trainer.m_totalIterations = total_iterations;
#else
    cout << "Training using the non-tree trainer." << endl;
    score = nontree_trainer.train();
#endif // USE_TREE_TRAINER
    cout << "Now (after training), the score is " << score << ", and the profile (length " << ROOT_PROFILE_PTR->length() << ") is:" << endl;
#ifdef USE_TREE_TRAINER
    // TODO: PUT BACK.
    if( tree_trainer.m_profileTree->nodeCount() < 5 ) {
      cout << *tree_trainer.m_profileTree << endl;
    } else {
      cout << *( tree_trainer.m_profileTree->getProfileTreeRoot() ) << endl;
    }
    //cout << endl;
    cout << "Its length is " << tree_trainer.m_profileTree->getProfileTreeRoot()->length() << endl;
    // TODO: REMOVE
    cout << "The tree is:" << endl;
    tree_trainer.m_profileTree->showGraph( cout );
# ifdef USE_NONTREE_TRAINER_TOO
    cout << "We also used the non-tree trainer afterwards, and got this profile (length " << profile.length() << "):" << endl;
    cout << profile << endl;
# endif // USE_NONTREE_TRAINER_TOO
    cout << "It took " << tree_trainer.m_totalIterations << " iterations." << endl;
#else
    cout << *nontree_trainer.m_profile << endl;
    cout << "Its length is " << nontree_trainer.m_profile->length() << endl;
    cout << "It took " << nontree_trainer.m_iteration << " iterations." << endl;
#endif // USE_TREE_TRAINER
    cout << "That score, again, was " << score << endl;

#ifdef USE_TREE_TRAINER  
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template TreeMultipleAlignment<ProfileTree<ResidueType, ProbabilityType, InternalNodeType>, SequenceResidueType > tree_multiple_alignment;
    if(
      ( tree_trainer.m_profileTree->nodeCount() > 1 )
    ) {
      tree_multiple_alignment.reinitialize(
        tree_trainer.m_profileTree,
        &fasta,
        num_sequences_to_use,
        &tree_trainer.m_profile_sequenceIndices,
        &tree_trainer.m_profileProfileAlignments
      );
  
      dp.alignLeaves(
        PROFILE_TRAINER.m_parameters,
        tree_multiple_alignment
      );
  
      cout << "The multiple alignment prints itself out like so:" << endl;
      //cout << tree_multiple_alignment << endl;
      //tree_multiple_alignment.toPairwiseStream( cout );
      //tree_multiple_alignment.toPileupStream( cout, &fasta.m_descriptions );
      tree_multiple_alignment.toAlignedFastaStream( cout, &fasta.m_descriptions );
    } // End if nodeCount() > 1
#endif // USE_TREE_TRAINER
  
    if( have_true_profile ) {
      // Do profile-profile alignment with true profile
      double indel_open_cost = .25;
      double indel_extension_cost = .25;
      vector<uint32_t> profile_profile_alignment;
      double profile_profile_alignment_cost =
        dp.profileProfile_align_SKL(
          PROFILE_TRAINER.m_parameters,
          profile_true,
          *ROOT_PROFILE_PTR,
          indel_open_cost,
          indel_extension_cost,
          profile_profile_alignment
        );
      double aligned_positions_cost = 0.0;
      uint32_t num_opens = 0;
      uint32_t num_extensions = 0;
      uint32_t num_aligned_positions = 0;
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
          num_aligned_positions += 1;
          if( profile_profile_alignment[ i ] > 1 ) {
            num_opens += 1;
            num_extensions += ( profile_profile_alignment[ i ] - 2 );
          }
        }
      }  // End foreach alignment position
      aligned_positions_cost =
        ( profile_profile_alignment_cost - ( num_opens * indel_open_cost ) - ( num_extensions * indel_extension_cost ) );
      cout << "After aligning the two profiles, got alignment cost: " << profile_profile_alignment_cost << endl;
      cout << "The alignment is ( " << profile_profile_alignment[ 0 ];
      for( uint32_t i = 1; i < profile_profile_alignment.size(); i++ ) {
        cout << ", " << profile_profile_alignment[ i ];
      }  // End foreach alignment position
      cout << " )" << endl;
      cout << "There were " << num_aligned_positions << " aligned positions." << endl;
      cout << "The symmeterized KL divergence of the aligned positions is " << aligned_positions_cost << endl;
      cout << "The average symmeterized KL divergence over the aligned positions is " << ( aligned_positions_cost / num_aligned_positions ) << endl;
    } // End if have_true_profile
  
    if( calculate_viterbi_score || calculate_multiple_alignment ) {
      if( !PROFILE_TRAINER_IS_TREE_TRAINER ) {
        // Free memory first
        nontree_trainer.reinitialize();
      } // TODO: else, free memory for tree trainer

      // Note that (for now) this will use the root profile only.
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer viterbi_matrices(
#ifdef USE_NONTREE_TRAINER_TOO
        *nontree_trainer.m_profile,
#else
        *ROOT_PROFILE_PTR,
#endif // USE_NONTREE_TRAINER_TOO
        fasta,
        num_sequences_to_use
      );
      score =
        dp.forward_score_viterbi(
          PROFILE_TRAINER.m_parameters,
#ifdef USE_NONTREE_TRAINER_TOO
          *nontree_trainer.m_profile,
#else
          *ROOT_PROFILE_PTR,
#endif // USE_NONTREE_TRAINER_TOO
          fasta,
          num_sequences_to_use,
          viterbi_matrices // TODO: ARE THESE MATRICES SIZED CORRECTLY?
        );
      cout << "The total score for all sequences, using viterbi, is: " << score << endl;
      // End calculating viterbi score
  
      if( calculate_multiple_alignment ) {
        // Show multiple alignment
        typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template MultipleAlignment<ProfileType, SequenceResidueType> ma(
#ifdef USE_NONTREE_TRAINER_TOO
          ( nontree_trainer.m_profile ),
#else
          ROOT_PROFILE_PTR,
#endif // USE_NONTREE_TRAINER_TOO
          &fasta,
          num_sequences_to_use
        );
        dp.forward_viterbiAlign(
          PROFILE_TRAINER.m_parameters,
          viterbi_matrices,
          ma
        );
        cout << "Viterbi score, calculated using ma.calculateScore(), is " << ma.calculateScore() << endl;
        cout << "Multiple Alignment is:" << endl;
        //cout << ma << endl;
        //ma.toPairwiseStream( cout, &fasta.m_descriptions );
#ifdef USE_TREE_TRAINER            
        if( // if we've previously called the alignLeaves(..) method..
          ( tree_trainer.m_profileTree->nodeCount() > 1 ) &&
          (
            tree_multiple_alignment.m_sequenceOrder.size() ==
            num_sequences_to_use
          )
        ) {
          for( uint32_t seq_i = 0; seq_i < num_sequences_to_use; seq_i++ ) {
            ma.m_sequenceOrder[ seq_i ] =
              tree_multiple_alignment.m_sequenceOrder[ seq_i ];
            cout << "Seq " << ma.m_sequenceOrder[ seq_i ] << ": ";
            cout << fasta.m_descriptions[ ma.m_sequenceOrder[ seq_i ] ] << endl;
          }
        } // End if we have a sequence order, use it.
#endif // USE_TREE_TRAINER
        ma.toPileupStream( cout, &fasta.m_descriptions );
      } // End if calculate_multiple_alignment
    } // End if calculate_viterbi_score || calculate_multiple_alignment
    
    return 0;
  } // test( int const &, char **, ResidueType const & )

  void
  getDefaultFasta ( Fasta<Dna> & fasta, int & initial_profile_length )
  {
    //initial_profile_length = 100;
    //fasta.fromFile( "fasta/dna.100.100.fasta" );
      
    initial_profile_length = 10;
    fasta.fromFile( "fasta/dna.10.100.fasta" );
    //fasta.fromFile( "fasta/dna.10.10.fasta" );
  
    //initial_profile_length = 10;
    //fasta.fromFile( "fasta/dna.biclustering.test.10.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test2.10.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test3.10.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test4.10.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test5.10.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test6.11.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test7.11.10.fasta" );
          
    //initial_profile_length = 10;
    //fasta.fromFile( "fasta/dna.biclustering.test.15.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test2.15.10.fasta" );
    // This one shares the A block and the G block with a C at pos 2:
    //fasta.fromFile( "fasta/dna.biclustering.test3.15.10.fasta" );
    // This one shares the A block and the T block with a C at pos 2:
    //fasta.fromFile( "fasta/dna.biclustering.test4.15.10.fasta" );
      
    //initial_profile_length = 21;
    //fasta.fromFile( "fasta/21U.fa" );
    //fasta.fromFile( "fasta/21U.last500.fasta" );
    //fasta.fromFile( "fasta/21U.last50.fasta" );
    //fasta.fromFile( "fasta/21U.last10.fasta" );
  
    // Testing handling of "N" ambiguity char
    //initial_profile_length = 3;
    //fasta.fromFile( "fasta/one_dna_seq_with_N.fasta" );
  
    //initial_profile_length = 0; // no idea
    //fasta.fromFile( "/Users/Paul/src/data/RepeatMasker/Tigger2/hg18_chr22.fa.noNs" );
    //fasta.fromFile( "fasta/one_dna_seq_without_N.fasta" );

    //initial_profile_length = 0; // no idea
    //fasta.fromFile( "Rfam/RF00376.seed.fasta.noNs" );

    return;
  } // getDefaultFasta ( Fasta<Dna> &, int & )

  void
  getDefaultFasta ( Fasta<Dna5> & fasta, int & initial_profile_length )
  {
    //initial_profile_length = 100;
    //fasta.fromFile( "fasta/dna.100.100.fasta" );
      
    //initial_profile_length = 10;
    //fasta.fromFile( "fasta/dna.10.100.fasta" );
    //fasta.fromFile( "fasta/dna.10.10.fasta" );
  
    //initial_profile_length = 10;
    //fasta.fromFile( "fasta/dna.biclustering.test.10.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test2.10.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test3.10.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test4.10.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test5.10.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test6.11.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test7.11.10.fasta" );
          
    //initial_profile_length = 10;
    //fasta.fromFile( "fasta/dna.biclustering.test.15.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test2.15.10.fasta" );
    // This one shares the A block and the G block with a C at pos 2:
    //fasta.fromFile( "fasta/dna.biclustering.test3.15.10.fasta" );
    // This one shares the A block and the T block with a C at pos 2:
    //fasta.fromFile( "fasta/dna.biclustering.test4.15.10.fasta" );
      
    //initial_profile_length = 21;
    //fasta.fromFile( "fasta/21U.fa" );
    //fasta.fromFile( "fasta/21U.last500.fasta" );
    //fasta.fromFile( "fasta/21U.last50.fasta" );
    //fasta.fromFile( "fasta/21U.last10.fasta" );
  
    // Testing handling of "N" ambiguity char
    //initial_profile_length = 3;
    //fasta.fromFile( "fasta/one_dna_seq_with_N.fasta" );
  
    //initial_profile_length = 0; // no idea
    //fasta.fromFile( "/Users/Paul/src/data/RepeatMasker/Tigger2/hg18_chr22.fa" );

    initial_profile_length = 0; // no idea
    fasta.fromFile( "Rfam/RF00376.seed.fasta" );
    //fasta.fromFile( "Rfam/RF00376.seed.fasta.one_seq_with_N" );
    //fasta.fromFile( "Rfam/RF00376.seed.fasta.two_seqs_one_has_N" );
    //fasta.fromFile( "Rfam/RF00376.seed.fasta.again" );

    return;
  } // getDefaultFasta ( Fasta<Dna5> &, int & )

  void
  getDefaultFasta ( Fasta<Iupac> & fasta, int & initial_profile_length )
  {
    // #define USE_DEL_IN_DEL_OUT and DISALLOW_FLANKING_TRANSITIONS first for this one:
    //initial_profile_length = 262;
    //fasta.fromFile( "data/RepeatMasker/MIR/hg18_chr22.allMIRs.fa.justMIR.minLength150" );
    //fasta.fromFile( "data/RepeatMasker/MIR/mm9/mm9_chr16-19.allMIRs.lt25pctDiv.fa.justMIR" );
    //fasta.fromFile( "deleteme.fasta" );

    initial_profile_length = 100;
    //fasta.fromFile( "test_results_hiatus_DNA/v9_seed1304287230_type/ProfuseTest.v9_seed1304287230_type.50.100.100.100.1.1.0.10000000000000001.0.10000000000000001.1.training_sequences.fasta" );
    //fasta.fromFile( "test_results/v8_seed1303599137_type/ProfuseTest.v8_seed1303599137_type.90.100.100.100.1.1.0.10000000000000001.0.10000000000000001.1.training_sequences.fasta" );
    //fasta.fromFile( "test_results/v8_seed1303717895_type/ProfuseTest.v8_seed1303717895_type.50.100.100.100.1.1.0.10000000000000001.0.10000000000000001.0.training_sequences.fasta" );
    //fasta.fromFile( "test_results/v8_seed1303842002_type/ProfuseTest.v8_seed1303842002_type.50.100.100.100.1.1.0.10000000000000001.0.10000000000000001.0.training_sequences.fasta" );
    //fasta.fromFile( "test_results/v8_seed1303849236_type/ProfuseTest.v8_seed1303849236_type.30.100.100.100.1.1.0.10000000000000001.0.10000000000000001.0.training_sequences.fasta" );

    //initial_profile_length = 100;
    fasta.fromFile( "dna.100.100.fasta" ); //tah 
    //fasta.fromFile( "fasta/dna.100.100.USE_DEL_IN_DEL_OUT.fasta" );
      
    //initial_profile_length = 10;
    //fasta.fromFile( "fasta/dna.10.100.fasta" );
    //fasta.fromFile( "fasta/dna.10.10.fasta" );
  
    //initial_profile_length = 10;
    //fasta.fromFile( "fasta/dna.biclustering.test.10.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test2.10.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test3.10.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test4.10.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test5.10.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test6.11.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test7.11.10.fasta" );
          
    //initial_profile_length = 10;
    //fasta.fromFile( "fasta/dna.biclustering.test.15.10.fasta" );
    //fasta.fromFile( "fasta/dna.biclustering.test2.15.10.fasta" );
    // This one shares the A block and the G block with a C at pos 2:
    //fasta.fromFile( "fasta/dna.biclustering.test3.15.10.fasta" );
    // This one shares the A block and the T block with a C at pos 2:
    //fasta.fromFile( "fasta/dna.biclustering.test4.15.10.fasta" );
      
    //initial_profile_length = 21;
    //fasta.fromFile( "fasta/21U.fa" );
    //fasta.fromFile( "fasta/21U.last500.fasta" );
    //fasta.fromFile( "fasta/21U.last50.fasta" );
    //fasta.fromFile( "fasta/21U.last10.fasta" );
  
    // Testing handling of "N" ambiguity char
    //initial_profile_length = 3;
    //fasta.fromFile( "fasta/one_dna_seq_with_N.fasta" );
  
    //initial_profile_length = 0; // no idea
    //fasta.fromFile( "/Users/Paul/src/data/RepeatMasker/Tigger2a/hg18_chr22.fa" );

    //initial_profile_length = 0; // no idea
    //initial_profile_length = 84; // this is what it settles on, with lengthadjust on
    //fasta.fromFile( "Rfam/RF00376.seed.fasta" );
    //fasta.fromFile( "Rfam/RF00376.seed.fasta.one_seq_with_N" );
    //fasta.fromFile( "Rfam/RF00376.seed.fasta.two_seqs_one_has_N" );
    //fasta.fromFile( "Rfam/RF00376.seed.fasta.again" );

    //initial_profile_length = 84;
    //fasta.fromFile( "Rfam/RF00376.full.fasta" );
    //fasta.fromFile( "Rfam/RF00376.full.fasta.one_seq_with_W" );

    return;
  } // getDefaultFasta ( Fasta<Iupac> &, int & )

  void
  getDefaultFasta ( Fasta<AminoAcid> & fasta, int & initial_profile_length )
  {
    initial_profile_length = 1; // use max length
    fasta.fromFile( "/Users/pedlefse/src/prefab4/in/1a02N_1a3qA" );

    //initial_profile_length = 100;
    //fasta.fromFile( "test_results/v9_seed1303861634_type/ProfuseTest.v9_seed1303861634_type.90.100.100.100.1.1.0.10000000000000001.0.10000000000000001.0.training_sequences.fasta" );

    //fasta.fromFile( "fasta/hbr.globins.fasta" );
    //initial_profile_length = 141;
  
    //fasta.fromFile( "fasta/sct.globins.fasta" );
    //initial_profile_length = 150;
  
    //fasta.fromFile( "fasta/sct3sdha.globins.fasta" );
    //initial_profile_length = 150;
  
    //fasta.fromFile( "fasta/hl0.globins.fasta" );
    //initial_profile_length = 158;
  
    //fasta.fromFile( "fasta/spg.globins.fasta" );
    //initial_profile_length = 144;
  
    //fasta.fromFile( "fasta/hbrspgi3da.globins.fasta" );
    //initial_profile_length = 146;
  
    //fasta.fromFile( "fasta/hbrspgi3dasct3sdha.globins.fasta" );
    //initial_profile_length = 150;
  
    //fasta.fromFile( "fasta/aminoacid.PREFAB.1b6e.fasta" );
    //initial_profile_length = 119;
  
    //fasta.fromFile( "fasta/aminoacid.biclustering.PREFAB.1b6e_1hyrB.fasta" );     //initial_profile_length = 119;
  
    //fasta.fromFile( "fasta/hemoglobins.fasta" );
    //initial_profile_length = 150;

    return;
  } // getDefaultFasta ( Fasta<AminoAcid> &, int & )
}; // End class Tests

} // End namespace galosh

int
main ( int const argc, char const ** argv )
{
#ifdef TEST_BASICS
  test_basics();
  return 0;
#endif

#ifdef __PROFUSE_USE_AMINOS
  typedef seqan::AminoAcid20 ResidueType;
  typedef seqan::AminoAcid SequenceResidueType;
#else // __PROFUSE_USE_AMINOS .. else
  typedef seqan::Dna ResidueType;
  typedef seqan::Iupac SequenceResidueType;
#endif // __PROFUSE_USE_AMINOS .. else ..

  //typedef bfloat ProbabilityType;
  //typedef logspace ProbabilityType;
  //typedef floatrealspace ProbabilityType;
  typedef doublerealspace ProbabilityType;
  
  typedef bfloat ScoreType; // Preferred
  //typedef logspace ScoreType; // SLOWer than bfloat
  //typedef realspace ScoreType; // Only for very few & small sequences
  //typedef longdoublerealspace ScoreType;
  
  typedef bfloat MatrixValueType;
  //typedef logspace MatrixValueType;
  //typedef doublerealspace MatrixValueType;
  //typedef longdoublerealspace MatrixValueType;
  //typedef floatrealspace MatrixValueType;

  galosh::Tests<ProbabilityType, ScoreType, MatrixValueType> tests;
  return tests.test( argc, argv, ResidueType(), SequenceResidueType() );
} // main (..)
