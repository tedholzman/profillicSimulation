/*---------------------------------------------------------------------------##
##  File:
##      @(#) ProfileTreeTrainer.hpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      Class definition for the ProfileTreeTrainer class, which repeatedly
##      uses ProfileTrainer objects to train an entire profile tree.
##
#******************************************************************************
#*  Copyright (C) 2007    Paul Edlefsen                                       *
#*  All rights reserved.                                                      *
#*****************************************************************************/

#if     _MSC_VER > 1000
#pragma once
#endif

#ifndef __GALOSH_PROFILETREETRAINER_HPP__
#define __GALOSH_PROFILETREETRAINER_HPP__

#include "Parameters.hpp"
using galosh::Parameters;
using galosh::DebugLevel;
using galosh::VerbosityLevel;

#include "Profile.hpp"
using galosh::ProfileTreeRoot;

#include "ProfileTree.hpp"
using galosh::ProfileTree;

#include "Sequence.hpp"
using galosh::Sequence;

#include "ProfileTrainer.hpp"
using galosh::ProfileTrainer;

#include <string>
using std::string;
#include <iostream>
using std::cout;
using std::endl;

#include <math.h> // isnan

#ifdef __HAVE_MUSCLE
#include "muscle/distfunc.h"
#include "muscle/clustsetdf.h"
#include "muscle/clust.h"
#include "muscle/tree.h"
#include "muscle/textfile.h"
#endif // __HAVE_MUSCLE

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

namespace galosh {

  // TODO: Move
  // Copied from seqan's basic_definition.h
  template <typename T>
  struct Tag
  {
  };

  //// Metrics for constructing distance/proximity matrices.
  // COOCExpectedCountDistance indicates the Coocurrence Expected Count proximity metric.
  struct COOCExpectedCountMetric_;
  typedef Tag<COOCExpectedCountMetric_> const COOCExpectedCountMetric;
  // COOCProbabilityMetric indicates the Coocurrence Probability proximity metric.
  struct COOCProbabilityMetric_;
  typedef Tag<COOCProbabilityMetric_> const COOCProbabilityMetric;
  // CrossEntropy indicates the CrossEntropy distance metric.
  struct CrossEntropyMetric_;
  typedef Tag<CrossEntropyMetric_> const CrossEntropyMetric;
  // SKL indicates the Symmeterized Kullback-Leibler Divergence distance metric.
  struct SKLMetric_;
  typedef Tag<SKLMetric_> const SKLMetric;
      
template <class ResidueType,
          class ProbabilityType,
          class ScoreType,
          class MatrixValueType,
          class SequenceResidueType>
  class ProfileTreeTrainer {
  public:
  
    typedef ProfileTreeRoot<ResidueType, ProbabilityType> InternalNodeType;
    typedef Sequence<SequenceResidueType> SequenceType;

    class Parameters :
    public ProfileTrainer<ProfileTreeRoot<ResidueType, ProbabilityType>, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::Parameters
    {
      // Boost serialization
    private:
      typedef typename ProfileTrainer<ProfileTreeRoot<ResidueType, ProbabilityType>,ScoreType,MatrixValueType,SequenceResidueType,InternalNodeType>::Parameters profile_trainer_parameters_t; 
      friend class boost::serialization::access;
      template<class Archive>
      void serialize ( Archive & ar, const unsigned int /* file_version */ )
      {
        // save/load base class information
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( profile_trainer_parameters_t );

        ar & BOOST_SERIALIZATION_NVP( shareProfilePositions );
        ar & BOOST_SERIALIZATION_NVP( shareProfilePositions_percentChangeTo0Threshold );
        ar & BOOST_SERIALIZATION_NVP( childSequenceMixtureThreshold );
      } // serialize( Archive &, const unsigned int )

    public:
  
      /// PARAMETERS
      /**
       * After training a subfamily profile with all positions differing from
       * its parent profile, should we try to retrain after identifying the
       * most closely similar profile positions and forcing the child profile
       * to use the parent profile at those positions?
       *
       * @see shareProfilePositions_percentChangeTo0Threshold
       */
      bool shareProfilePositions;
  #define DEFAULT_shareProfilePositions false

      /**
       * After training a subfamily profile with all positions differing from
       * its parent profile, we try to retrain after identifying the most
       * closely similar profile positions and forcing the child profile to use
       * the parent profile at those positions.  We do this until the score
       * after removing a position drops below some threshold.
       *
       * This is the threshold for the percent by which the score has changed
       * from the original (all-positions-different) profile.
       *
       * @see shareProfilePositions
       */
      double shareProfilePositions_percentChangeTo0Threshold;
  #define DEFAULT_shareProfilePositions_percentChangeTo0Threshold -50

      /**
       * We have to make a decision about whether a sequence belongs to the
       * child profile or the parent profile before recursively breaking the
       * profile into more subfamilies.  If the sequence subfamily mixture
       * parameter exceeds this threshold, we say that the sequence belongs to
       * the child profile.
       */
      double childSequenceMixtureThreshold;
  #define DEFAULT_childSequenceMixtureThreshold .5

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

      virtual void
      copyFrom ( const Parameters & copy_from );
    
      virtual void
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

    }; // End inner class Parameters

    template <class ParametersType>
    class ParametersModifierTemplate :
      public ProfileTrainer<ProfileTreeRoot<ResidueType, ProbabilityType>,ScoreType,MatrixValueType,SequenceResidueType,InternalNodeType>::template ParametersModifierTemplate<ParametersType>
    {
      typedef typename ProfileTrainer<ProfileTreeRoot<ResidueType, ProbabilityType>,ScoreType,MatrixValueType,SequenceResidueType,InternalNodeType>::template ParametersModifierTemplate<ParametersType> base_parameters_modifier_t; 

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
        ar & BOOST_SERIALIZATION_NVP( isModified_shareProfilePositions );
        ar & BOOST_SERIALIZATION_NVP( isModified_shareProfilePositions_percentChangeTo0Threshold );
        ar & BOOST_SERIALIZATION_NVP( isModified_childSequenceMixtureThreshold );
      } // serialize( Archive &, const unsigned int )

    public:
  
      /// isModified flags for Parameters
      /**
       * After training a subfamily profile with all positions differing from
       * its parent profile, should we try to retrain after identifying the
       * most closely similar profile positions and forcing the child profile
       * to use the parent profile at those positions?
       *
       * @see shareProfilePositions_percentChangeTo0Threshold
       */
      bool isModified_shareProfilePositions;

      /**
       * After training a subfamily profile with all positions differing from
       * its parent profile, we try to retrain after identifying the most
       * closely similar profile positions and forcing the child profile to use
       * the parent profile at those positions.  We do this until the score
       * after removing a position drops below some threshold.
       *
       * This is the threshold for the percent by which the score has changed
       * from the original (all-positions-different) profile.
       *
       * @see shareProfilePositions
       */
      bool isModified_shareProfilePositions_percentChangeTo0Threshold;

      /**
       * We have to make a decision about whether a sequence belongs to the
       * child profile or the parent profile before recursively breaking the
       * profile into more subfamilies.  If the sequence subfamily mixture
       * parameter exceeds this threshold, we say that the sequence belongs to
       * the child profile.
       */
      bool isModified_childSequenceMixtureThreshold;

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
        writeParametersModifier( os );

        return os;
      } // friend operator<< ( basic_ostream &, ParametersModifierTemplate const& )

      template<class CharT, class Traits>
      void
      writeParametersModifier (
        std::basic_ostream<CharT,Traits>& os
      );

      template<class AnyParameters>
      void
      applyModifications ( AnyParameters & target_parameters );

    }; // End inner class ParametersModifierTemplate

    typedef ParametersModifierTemplate<typename ProfileTreeTrainer::Parameters> ParametersModifier;

    /**
     * The parameters currently being used.
     */
    Parameters m_parameters;

  ProfileTree<ResidueType, ProbabilityType, InternalNodeType> * m_profileTree;

    vector<SequenceType> const & m_sequences;

    /**
     * The number of sequences to use in training (the index one greater than
     * that of the last one to be used; must be <= m_sequences.size().
     */
    uint32_t m_sequence_count;

    // TODO: Make this a type, rather than a class.
    DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType> m_dynamic_programming;

    /**
     * After each round of subfamily training, we break the sequences into
     * groups, assigning them to either the parent or the child profile, by
     * comparing the sequenceSubfamilyIdentifier parameters to
     * parameters.childSequenceMixtureThreshold.  This is a vector of vectors,
     * one for each profile (including the root).  Each contained vector has a
     * unique set of sequence indices, and every sequence index is in exactly
     * one of the contained vectors.
     */
    vector<vector<uint32_t > > m_profile_sequenceIndices;

    /**
     * Per-node scores.
     */
    vector<ScoreType> m_profile_scores;

    /**
     * Per-node alignments among children (assuming nodes have 2 or no
     * children).  Will be zero-length if the node has no children.
     */
    vector<vector<uint32_t> > m_profileProfileAlignments;

    /**
     * The total number of trainer iterations elapsed in all ProfileTrainers
     * used by this ProfileTreeTrainer's train() method.
     */
    uint32_t m_totalIterations;
    
    /**
     * This is the starting score.  It is not changed once it is set.
     */
    // Not used
    //ScoreType m_scoreBeforeTraining;
    
    /** 
     * The score after the performing the current or most recent training
     * step.
     */
    // Not used
    //ScoreType m_endingScore;

    /**
     * Construct a profile tree trainer with the given tree and sequences.
     */  
    ProfileTreeTrainer (
      vector<SequenceType> const & sequences,
      ProfileTree<ResidueType, ProbabilityType, InternalNodeType> * profile_tree
    ) :
      m_profileTree( profile_tree ),
      m_sequences( sequences ),
      m_sequence_count( sequences.size() ),
      m_profile_sequenceIndices(), // see restart()
      m_profileProfileAlignments(), // see restart()
      m_totalIterations( 0 )
    {
      if( m_parameters.debug >= DEBUG_All ) {
        cout << "[debug] ProfileTreeTrainer::<init>( vector<SequenceType> const &, ProfileTree<ResidueType, ProbabilityType, InternalNodeType> & )" << endl;
      } // End if DEBUG_All
      // Do nothing else
    } // <init>( vector<SequenceType> const &, ProfileTree<ResidueType, ProbabilityType, InternalNodeType> & )

    /**
     * Construct a profile trainer with the given tree and sequences, and the
     * number of sequences to use (use the first num_sequences_to_use sequences
     * only).
     */  
    ProfileTreeTrainer (
      vector<SequenceType> const & sequences,
      uint32_t const & num_sequences_to_use, // use only the first X sequences...
      ProfileTree<ResidueType, ProbabilityType, InternalNodeType> * profile_tree
    ) :
      m_profileTree( profile_tree ),
      m_sequences( sequences ),
      m_sequence_count( min( ( size_t )num_sequences_to_use, sequences.size() ) ),
      m_profile_sequenceIndices(), // see restart()
      m_profileProfileAlignments(), // see restart()
      m_totalIterations( 0 )
    {
      if( m_parameters.debug >= DEBUG_All ) {
        cout << "[debug] ProfileTreeTrainer::<init>( vector<SequenceType> const &, uint32_t const &, ProfileTree<ResidueType, ProbabilityType, InternalNodeType> & )" << endl;
      } // End if DEBUG_All
      // Do nothing else
    } // <init>( vector<SequenceType> const &, uint32_t const &, ProfileTree<ResidueType, ProbabilityType, InternalNodeType> & )

    /**
     * Reinitialize this profile tree trainer with the given tree.
     */  
    void
    reinitialize (
      ProfileTree<ResidueType, ProbabilityType, InternalNodeType> * profile_tree
    )
    {
      m_profileTree = profile_tree;
    } // reinitialize( ProfileTree<ResidueType, ProbabilityType, InternalNodeType> & )

    /**
     * Utility method to convert the sort of profile that we get at each node
     * after running train_usingAlignmentProfiles(..) to a standard profile.
     * After running train_usingAlignmentProfiles(..) the profiles have counts
     * at each position's Match emission distribution parameter, rather than
     * probabilities.  We could simply normalize, but we want positions with
     * fewer counts to have more entropy.  The optional uint32_t argument
     * (maximum_count) is the number of sequences used to train the profile
     * (which is the maximum number of counts possible at any one position) by
     * default [Note that if you call it with a 0 maximum_count argument, the
     * default will be used instead].  The optional double argument
     * (minimum_value) is m_parameters.profileValueMinimum by default [Note
     * that if you call it with a negative minimum_value argument, the default
     * will be used instead].
     */
    void
    normalizeProfileOfCounts (
      InternalNodeType & profile,
      uint32_t maximum_count = 0,
      double minimum_value = -1
    ) const
    {
      if( maximum_count == 0 ) {
        maximum_count = m_sequences.size();
      }
      if( minimum_value < 0 ) {
        minimum_value = m_parameters.profileValueMinimum;
      }
      // The profile is just a bunch of counts right now..
      // We want positions with fewer counts to have increased entropy.
      // We make all positions have the same counts by adding the missing
      // counts to all positions equally. (TODO: Use priors?)
      float missing_count;
      for( uint32_t pos_i = 0; pos_i < profile.length(); pos_i++ ) {
        missing_count =
          ( maximum_count - ( float )toDouble( profile[ pos_i ][ Emission::Match ].total() ) );
        profile[ pos_i ][ Emission::Match ] +=
          (
            missing_count /
            ( float )( seqan::ValueSize<ResidueType>::VALUE )
          );
        // TODO: REMOVE
        if( abs( toDouble( profile[ pos_i ][ Emission::Match ].total() ) - maximum_count ) >= 1E-5 ) {
          cout << "UH-OH: The total Match emission count at position " << pos_i << " should be " << maximum_count << ", but it is " << toDouble( profile[ pos_i ][ Emission::Match ].total() ) << endl;
          cout << "\t missing_count is " << missing_count << endl;
          assert( abs( toDouble( profile[ pos_i ][ Emission::Match ].total() ) - maximum_count ) < 1E-5 );
        }
      } // End foreach pos, make the counts all add to maximum_count
      profile.normalize( minimum_value );

      return;
    } // normalizeProfileOfCounts( InternalNodeType & profile [, uint32_t [, double ]] ) const

    /**
     * Set this trainer to its initial values.
     */
    void
    restart ()
    {
      if( m_parameters.debug >= DEBUG_All ) {
        cout << "[debug] ProfileTreeTrainer::restart()" << endl;
      } // End if DEBUG_All

      // TODO: REMOVE
      //cout << "SequenceTypes are: " << endl;
      //for( uint32_t seq_i = 0; seq_i < m_sequence_count; seq_i++ ) {
      //  cout << m_sequences[ seq_i ] << endl;
      //}

      // Initialize the vector of vectors of subfamily id probs
      // There is one vector for each node except for the root.
      uint32_t tree_size = m_profileTree->nodeCount();
      // Initialize the vector of vectors of sequence indices
      // There is one vector for each profile (including the root)
      if( m_profile_sequenceIndices.size() !=
          tree_size ) {
        m_profile_sequenceIndices.resize( tree_size );
      }
      // Set up the sequenceIndices for the first node to include all of the sequences (for now).
      m_profile_sequenceIndices[ 0 ].resize( m_sequence_count );
      for( uint32_t seq_i = 0; seq_i < m_sequence_count; seq_i++ ) {
        m_profile_sequenceIndices[ 0 ][ seq_i ] = seq_i;
      }
      // Each other vector is of length 0, for now.
      for( uint32_t node_i = 1; node_i < tree_size; node_i++ ) {
        if( m_profile_sequenceIndices[ node_i ].size() != 0 ) {
          m_profile_sequenceIndices[ node_i ].resize( 0 );
        }
      }

      m_totalIterations = 0;
    } // restart()

    ScoreType
    train ()
    {
      // TODO: Make this a parameter!
      return
        train_usingAlignmentProfiles();
    } // train()

    /**
     * Divide and conquer using AlignmentProfiles, ala Muscle, sort-of.
     */
    ScoreType
    train_usingAlignmentProfiles ()
    {

      // First train everything together.
      if( m_parameters.debug >= DEBUG_All ) {
        cout << "[debug] ProfileTreeTrainer::train_usingAlignmentProfiles()" << endl;
      } // End if DEBUG_All

      restart();

      // TODO: REMOVE?
      static const bool show_profiles = false;
      static const bool be_extra_verbose = false;

      //uint32_t seq_i;
    
      ProfileTrainer<ProfileTreeRoot<ResidueType, ProbabilityType>, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType> trainer(
          m_profileTree->getProfileTreeRoot(),
          m_sequences,
          m_sequence_count
        );
      // Use our own special parameters
      trainer.m_parameters = m_parameters;

      if( show_profiles ) {
        cout << "The root profile (before) is:" << endl;
        cout << ( *m_profileTree->getProfileTreeRoot() ) << endl;
      } // End if show_profiles

      //trainer.m_parameters.verbosity = VERBOSITY_All;
      //trainer.m_parameters.debug = DEBUG_Medium;

      // Do it.
      ScoreType onefamily_score;

      // TODO: Make a parameter?
      bool have_trained_profile = false;
      if( have_trained_profile ) {
        if( be_extra_verbose ) {
        cout << "Calculating forward score.." << endl;
        }
        // Then don't retrain it.
        onefamily_score = m_dynamic_programming.forward_score(
          trainer.m_parameters,
          ( *m_profileTree->getProfileTreeRoot() ),
          m_sequences,
          m_sequence_count
        );
        m_totalIterations = 0;
      } else {
        onefamily_score = trainer.train();
        m_totalIterations = trainer.m_iteration;
      }

      // TODO: REMOVE.  For profusetest.  TODO: something better.
      //return onefamily_score;

      // TODO: REMOVE?
      if( m_profileTree->nodeCount() > 1 ) {
        cout << "Now (after training everything as one family), the score is " << onefamily_score << "." << endl;
      }

      if( show_profiles ) {
        cout << "The root profile (after training everything as one family) is:" << endl;
        cout << *trainer.m_profile << endl;
      }

      // TODO: Get it working without relying on Muscle..
#ifndef __HAVE_MUSCLE
      return onefamily_score;
#else // __HAVE_MUSCLE

      m_profile_sequenceIndices.resize( 1 );
      m_profile_sequenceIndices[ 0 ].resize( m_sequence_count );
      for( uint32_t seq_i = 0; seq_i < m_sequence_count; seq_i++ ) {
        m_profile_sequenceIndices[ 0 ][ seq_i ] = seq_i;
      }
      m_profile_scores.resize( 1 );
      m_profile_scores[ 0 ] = onefamily_score;

      // TODO: REMOVE
      ScoreType score_before = onefamily_score;
      ScoreType score_after;
      ScoreType overall_score = onefamily_score;
      if( be_extra_verbose ) {
        cout << "(Onefamily) Overall score is now " << overall_score << endl;
      }
      uint32_t node_to_split = 0;
      // If node_to_split is not the root, which child of its parent is it? (min 1)
      uint32_t node_to_split_which_child;
      uint32_t node_to_split_parent;

      // TODO: Make these parameters

      // ERE I AM.  Trying to figure out what the right value is here.  Note
      // that for DNA, the maximum possible skl distance between two emission
      // distributions is 1416.79. -- if there's prob 1 on different residues.
      const double indel_open_cost = 1;//1400.0;//20.0;//10.0;//.25;
      const double indel_extension_cost = 1;//1400.0;//20.0;//10.0;//.25;

      // If true, do a full (re)training at every split.
      static const bool retrain_children = false;

      //double aligned_positions_cost
      //uint32_t num_opens;
      //uint32_t num_extensions;
      uint32_t num_aligned_positions;
      vector<uint32_t> profile_profile_alignment;

      uint32_t total_positions_in_parent;
      uint32_t parent_position;
      double child_1_weight;
      uint32_t child_2_pos;

      InternalNodeType normalized_child_1;
      InternalNodeType normalized_child_2;

      Tree muscle_tree;
      unsigned muscle_node_id = 0;
      vector<SequenceType> child_sequences;
      uint32_t right_vertex, left_vertex;
      // Do a depth-first traversal
      do {
        score_before = score_after = m_profile_scores[ node_to_split ];

        // Don't try to split a node if there's only one sequence in it.
        if( m_profile_sequenceIndices[ node_to_split ].size() == 1 ) {
          if( be_extra_verbose ) {
            cout << "Can't split leaf node " << node_to_split << endl;
          }
          // Rejected.  Don't decend further here.
          while( ( node_to_split != 0 ) && ( node_to_split_which_child == m_profileTree->childCount( node_to_split_parent ) ) ) {
            if( be_extra_verbose ) {
              cout << "This is the last child of its parent.  Ascending tree." << endl;
            }
            if( be_extra_verbose ) {
              cout << "node_to_split is currently " << node_to_split << endl;
              cout << "node_to_split_parent is currently " << node_to_split_parent << endl;
              cout << "node_to_split_which_child is currently " << node_to_split_which_child << endl;
            } // End if be_extra_verbose

            const uint32_t child_count =
              m_profileTree->childCount( node_to_split_parent );
            // TODO: For now we are assuming that child_count is 2!
            assert( child_count == 2 );

            // ERE I AM.  Can we assume that there are exactly two children?  It seems that we should do the alignment allowing D->I transitions, but make a distinction between true matches and gap-matches: only true matches are to become positions in the new parent profile.  The problem is that right now we lose that information (we "allow" those transitions via a hack, and our return vector communicates only advances, which just isn't rich enough a representation).  I could modify it to optionally also return gap-gap match indicators.  Instead I could just make it include even the insertions in the new parent profile.  The easiest right now is probably to just treat gap-gap matches like matches but when a gap is just a gap, don't include it in the parent profile.

            InternalNodeType & parent =
              m_profileTree->getProfileTreeInternalNode( node_to_split_parent );
            InternalNodeType & child_1 =
              m_profileTree->getChild( node_to_split_parent, 1 );
            InternalNodeType & child_2 =
              m_profileTree->getChild( node_to_split_parent, 2 );

            // Since we are doing SKL-style alignment, we need to normalize the
            // children before aligning them.  But we don't actually want to
            // normalize them just yet, since we blend the children into the
            // parent unnormalized.
            normalized_child_1.copyFrom( child_1 );
            // TODO: REMOVE? Testing.
            normalizeProfileOfCounts(
              normalized_child_1,
              0, // Use the total number of sequences, effectively making leaves more entropic than internal nodes.
              //m_profile_sequenceIndices[ m_profileTree->getChildVertex( node_to_split_parent, 1 ) ].size(), // The maximum number of counts is the number of sequences assigned to child_1
              0 // normalize with 0 minimum value
            );
            // TODO: PUT BACK?
            //normalized_child_1.normalize( 0 );
            normalized_child_2.copyFrom( child_2 );
            // TODO: REMOVE? Testing.
            normalizeProfileOfCounts(
              normalized_child_2,
              0, // Use the total number of sequences, effectively making leaves more entropic than internal nodes.
              //m_profile_sequenceIndices[ m_profileTree->getChildVertex( node_to_split_parent, 2 ) ].size(), // The maximum number of counts is the number of sequences assigned to child_2
              0 // normalize with 0 minimum value
            );
            // TODO: PUT BACK?
            //normalized_child_2.normalize( 0 );

            // TODO: REMOVE
            //cout << "About to call profileProfile_align_SKL( [vertex " << m_profileTree->getChildVertex( node_to_split_parent, 1 ) << ", with " << m_profile_sequenceIndices[ m_profileTree->getChildVertex( node_to_split_parent, 1 ) ].size() << " sequences], [vertex " << m_profileTree->getChildVertex( node_to_split_parent, 2 ) << ", with " << m_profile_sequenceIndices[ m_profileTree->getChildVertex( node_to_split_parent, 2 ) ].size() << " sequences] )" << endl;
            //cout << "\tvertex " << m_profileTree->getChildVertex( node_to_split_parent, 1 ) << " is:" << endl << child_1 << "\t\tnormalized:" << endl << normalized_child_1;
            //cout << "\tvertex " << m_profileTree->getChildVertex( node_to_split_parent, 2 ) << " is:" << endl << child_2 << "\t\tnormalized:" << endl << normalized_child_2 << endl;

            // Get the profile-profile alignment.
            m_dynamic_programming.profileProfile_align_SKL(
              m_parameters,
              normalized_child_1,
              normalized_child_2,
              indel_open_cost,
              indel_extension_cost,
              profile_profile_alignment,
              ( indel_open_cost * 2 ) // do allow "gap matches"
            );
            // TODO: REMOVE
            //exit( 0 );
            if( be_extra_verbose ) {
              cout << "The alignment between the children of node " << node_to_split_parent << " is ( " << profile_profile_alignment[ 0 ];
              for( uint32_t i = 1; i < profile_profile_alignment.size(); i++ ) {
                cout << ", " << profile_profile_alignment[ i ];
              }  // End foreach alignment position
              cout << " )" << endl;
            } // End if be_extra_verbose

            // Save it.
            if( m_profileProfileAlignments.size() < ( node_to_split_parent + 1 ) ) {
              m_profileProfileAlignments.resize( node_to_split_parent + 1 );
            }
            m_profileProfileAlignments[ node_to_split_parent ] =
              profile_profile_alignment;

            // Okay now we need to get the number of matches in the alignment.
            num_aligned_positions = 0;
            total_positions_in_parent = profile_profile_alignment[ 0 ];
            for( uint32_t i = 1; i < profile_profile_alignment.size(); i++ ) {
              if( profile_profile_alignment[ i ] == 0 ) {
                total_positions_in_parent += 1;
              } else {
                num_aligned_positions += 1;
                total_positions_in_parent += profile_profile_alignment[ i ];
              }
            }  // End foreach alignment position
            if( be_extra_verbose ) {
              cout << "There are " << num_aligned_positions << " positions that align between the children of the profile at node " << node_to_split_parent << "." << endl;
              cout << "\tThere are " << total_positions_in_parent << " total positions in the parent, including parts deleted in one child." << endl;
            }
            // Set the length of the parent..
            parent.reinitialize( total_positions_in_parent );

            // Get the indel values from the root, I guess.
            if( node_to_split_parent == 0 ) {
              // If this *is* the root, restore globals from a child.
              parent.copyExceptPositions( child_1 );
            } else {
              parent.copyExceptPositions( *( m_profileTree->getProfileTreeRoot() ) );
            }

            // Make each match position a weighted average of the positions of
            // the two children.

            // The weights are relative to the number of sequences assigned to
            // each branch.
            assert( m_profile_sequenceIndices[ child_2.getProfileTreeVertex() ].size() > 0 );
            //child_1_weight =
            //  (
            //    ( double )m_profile_sequenceIndices[ child_1.getProfileTreeVertex() ].size() /
            //    ( double )m_profile_sequenceIndices[ child_2.getProfileTreeVertex() ].size()
            //  );
            //if( be_extra_verbose ) {
            //  cout << "The weight on the first child's values is " << child_1_weight << " (the weight on the second child's values is 1.0)." << endl;
            //}
            child_2_pos = profile_profile_alignment[ 0 ];
            for( parent_position = 0;
                 parent_position < child_2_pos;
                 parent_position++ ) {
              // These are from child_2 only
              parent[ parent_position ].copyFrom(
                child_2[ parent_position ]
              );
            }
            for( uint32_t child_1_pos_plus1 = 1;
                 child_1_pos_plus1 < profile_profile_alignment.size();
                 child_1_pos_plus1++ ) {
              if( profile_profile_alignment[ child_1_pos_plus1 ] == 0 ) {
                // From child_1 only.
                parent[ parent_position ].copyFrom(
                  child_1[ child_1_pos_plus1 - 1 ]
                );
                parent_position += 1;
              } else {
                parent[ parent_position ].copyFrom(
                  child_1[ child_1_pos_plus1 - 1 ]
                );
                //parent[ parent_position ] *= child_1_weight;
                parent[ parent_position ] +=
                  child_2[ child_2_pos ];
                //parent[ parent_position ].normalize(
                //  0 //m_parameters.profileValueMinimum
                //);

                // Now advance these..
                child_2_pos += 1;
                parent_position += 1;

                // The insertions are unique to child 2
                for( uint32_t tmp_i = 2;
                     tmp_i <= profile_profile_alignment[ child_1_pos_plus1 ];
                     tmp_i++
                ) {
                  // These are from child_2 only
                  parent[ parent_position ].copyFrom(
                    child_2[ child_2_pos ]
                  );
                  parent_position += 1;
                  child_2_pos += 1;
                }

              } // End if there are no child_2 advances .. else there are ..
            }  // End foreach alignment position, set the parent profile
               // values.

            assert( parent_position == total_positions_in_parent );

            if( be_extra_verbose ) {
              if( show_profiles ) {
                cout << "The new profile at node " << node_to_split_parent << " is:" << endl;
                cout << parent << endl;
              }
            } // End if be_extra_verbose

            // Before moving up, normalize the children nodes (now that we
            // won't need them unnormalized anymore).
            child_1.normalize( m_parameters.profileValueMinimum );
            child_2.normalize( m_parameters.profileValueMinimum );

            // We're at the last child of its parent.  Move up.
            node_to_split = node_to_split_parent;
            if( node_to_split != 0 ) {
              node_to_split_parent = m_profileTree->getParentVertex( node_to_split );
              node_to_split_which_child =
                m_profileTree->getChildIndexInParent( node_to_split_parent, node_to_split );
            }
            if( !retrain_children ) {
              muscle_node_id = muscle_tree.GetParent( muscle_node_id );
              if( be_extra_verbose ) {
                cout << "muscle_node_id is now " << muscle_node_id << endl;
              } // End if be_extra_verbose
            } // End if !retrain_children
          } // While we're at the last child of our parent, move up.
          if( node_to_split == 0 ) {
            // Make sure that the m_profileProfileAlignments vector is long
            // enough (leaves at the end won't have been added to it).
            if(
              m_profileProfileAlignments.size() !=
              m_profile_sequenceIndices.size()
            ) {
              m_profileProfileAlignments.resize(
                m_profile_sequenceIndices.size()
              );
            }

            // Can't move up.  Done.
            return overall_score;
          }
          if( be_extra_verbose ) {
            cout << "Moving to next sibling." << endl;
          }
          // Move along to the next child.
          node_to_split_which_child += 1;
          node_to_split =
            m_profileTree->getChildVertex( node_to_split_parent, node_to_split_which_child );
          if( be_extra_verbose ) {
            cout << "node_to_split is now " << node_to_split << endl;
          } // End if be_extra_verbose
          if( !retrain_children ) {
            assert( node_to_split_which_child == 2 );
            // left child
            muscle_node_id =
              muscle_tree.GetLeft( muscle_tree.GetParent( muscle_node_id ) );
            if( be_extra_verbose ) {
              cout << "muscle_node_id is now " << muscle_node_id << endl;
            } // End if be_extra_verbose
          } // End if !retrain_children
        } else {
          if( be_extra_verbose ) {
            cout << "Splitting node " << node_to_split << "." << endl;
          }
          if( retrain_children || ( node_to_split == 0 ) ) {
            // We only need to calculate the tree once if we're not retraining
            // the children.
            if( be_extra_verbose ) {
              cout << "Creating muscle tree." << endl;
            }
            createMuscleTree( node_to_split, muscle_tree );
            muscle_node_id =
              muscle_tree.GetRootNodeIndex();
            if( be_extra_verbose ) {
              cout << "Its root node index is " << muscle_node_id << "." << endl;
            }
          }
          if( be_extra_verbose ) {
            cout << "Calling splitNode( " << node_to_split << ", muscle_tree, " << muscle_node_id << ", " << retrain_children << " )" << endl;
          }
          splitNode(
            node_to_split,
            muscle_tree,
            muscle_node_id,
            retrain_children
          );

          InternalNodeType & right_profile =
            m_profileTree->getChild( node_to_split, 1 );
          right_vertex =
            right_profile.getProfileTreeVertex();
          // Retrain only for leaves unless retrain_children is true.
          if(
            retrain_children ||
            ( m_profile_sequenceIndices[ right_vertex ].size() == 1 )
          ) {
            getSequencesFromIndices(
              m_profile_sequenceIndices[ right_vertex ],
              child_sequences
            );
            m_profile_scores[ right_vertex ] =
              createProfileFromSequences(
                child_sequences,
                right_profile
              );
            if( be_extra_verbose ) {
              cout << "The right alone score is " << m_profile_scores[ right_vertex ] << "." << endl;
            } // End if be_extra_verbose
            if( show_profiles ) {
              cout << "The right profile (after training) is:" << endl;
              cout << right_profile << endl;
            }
          } else {
            m_profile_scores[ right_vertex ] = 
              m_profile_scores[ node_to_split ];
          } // End if retrain_children or right node is a leaf, retrain right
            // node .. else ..

          InternalNodeType & left_profile =
            m_profileTree->getChild( node_to_split, 2 );
          left_vertex =
            left_profile.getProfileTreeVertex();
          // Retrain only for leaves unless retrain_children is true.
          if(
            retrain_children ||
            ( m_profile_sequenceIndices[ left_vertex ].size() == 1 )
          ) {
            getSequencesFromIndices(
              m_profile_sequenceIndices[ left_vertex ],
              child_sequences
            );
            m_profile_scores[ left_vertex ] =
              createProfileFromSequences(
                child_sequences,
                left_profile
              );

            if( be_extra_verbose ) {
              cout << "The left alone score is " << m_profile_scores[ left_vertex ] << "." << endl;
            } // End if be_extra_verbose

            if( show_profiles ) {
              cout << "The left profile (after training) is:" << endl;
              cout << left_profile << endl;
            }
          } else {
            // ? In this case the scores are meaningless, I guess.
            m_profile_scores[ left_vertex ] = 
              m_profile_scores[ node_to_split ];
          } // End if retrain_children or left node is a leaf, retrain left
            // node .. else ..

          score_after =
            (
              m_profile_scores[ right_vertex ] *
              m_profile_scores[ left_vertex ]
            );

          if( be_extra_verbose ) {
            cout << "The new score is " << score_after << " (the old score was " << score_before << ")." << endl;
          }
          overall_score /= score_before;
          overall_score *= score_after;
          if( be_extra_verbose ) {
            cout << "Overall score is now " << overall_score << endl;
          }

          // Walk down.
          // Try splitting the first newly-created child.
          node_to_split_which_child = 1;
          node_to_split_parent = node_to_split;
          node_to_split =
            m_profileTree->getChildVertex( node_to_split, node_to_split_which_child );
          assert( node_to_split != 0 ); // getChildVertex will return 0 if it gets confused.
          if( be_extra_verbose ) {
            cout << "node_to_split is now " << node_to_split << endl;
          } // End if be_extra_verbose
          if( !retrain_children ) {
            // the first child is the right child (we're moving down)
            muscle_node_id = muscle_tree.GetRight( muscle_node_id );
            if( be_extra_verbose ) {
              cout << "muscle_node_id is now " << muscle_node_id << endl;
            } // End if be_extra_verbose
          } // End if !retrain_children
        } // If the score didn't improve .. else ..
      } while( true ); // We break when done.
      assert( false /* we should never get to this point. */ );
#endif // __HAVE_MUSCLE

    } // train_usingAlignmentProfiles()

    /**
     * Fill the given vector (reference) to contain the sequences with the
     * given indices.
     */
    void
    getSequencesFromIndices (
      vector<uint32_t> const & indices,
      vector<SequenceType> & sequences
    ) const {
      uint32_t num_seqs = indices.size();
      sequences.resize( num_seqs );
      for( uint32_t seq_i = 0; seq_i < num_seqs; seq_i++ ) {
        sequences[ seq_i ] =
          m_sequences[ indices[ seq_i ] ];
      }
      return;
    } // getSequencesFromIndices( vector<uint32_t> const &, vector<SequenceType> & ) const

    /**
     * Calculate and return the (unscaled but unnormalized) alignment profiles
     * between the sequences under node node_id and the profile at that node.
     * Note that the Match probabilities do not sum to 1.  They sum to the
     * probability of a Match at the corresponding position.  This is great for
     * passing to the calculateMatchEmissionCooccurrenceExpectedCount() and
     * calculateMatchEmissionCooccurrenceProbability() methods, but the
     * crossEntropy and KL-distance metrics require all distributions to be
     * proper distributions (that is, they must sum to 1.0).
     */
    // TODO: ERE I AM.  I have just discovered that all the rigamarole about calculating the deletion fraction etc is wasted effort: the sum of the (unscaled) Match emissions tells the whole story, even in the last row.  The insertions may still need some fancy calculation in ProfileTrainer, but we can do away with any fanciness for the deletions -- ie we don't need to wait until the next-to-last row to do last-row deletions.
    // TODO: I'm presently going to check that in ProfileTrainer we don't assume that the AlignmentProfiles are weighted by the sequence score, because they are not! -- Actually they are.  They are the conditional expected counts, so if all is right then the sum of the match emissions is the conditional probability of a match, given that the sequence is created (that is, divided by the sequence score).  Anyway, it's all good!
    void
    calculateAlignmentProfiles (
      uint32_t const & node_id,
      vector<typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile> & aps
    ) const {
      static const bool be_extra_verbose = true;

      uint32_t const node_sequence_count =
        m_profile_sequenceIndices[ node_id ].size();

      uint32_t node_length;
      if( node_id == 0 ) {
        node_length = m_profileTree->getProfileTreeRoot()->length();
      } else {
        node_length = m_profileTree->getProfileTreeInternalNode( node_id ).length();
      }

      aps.resize( node_sequence_count );
      for( uint32_t seq_i = 0; seq_i < node_sequence_count; seq_i++ ) {
        aps[ seq_i ].reinitialize(
          (
            node_length + 1
          )
        );
      }
      // TODO: This should depend on the parameter we'll add to ProfileTrainer to determine whether to use forward matrices.
      vector<SequenceType> sequences;
      getSequencesFromIndices(
        m_profile_sequenceIndices[ node_id ],
        sequences
      );
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer forward_matrices(
        node_length,
        sequences,
        node_sequence_count
      );
      // TODO: REMOVE?  We can pass NULLs instead...
      vector<ScoreType> sequence_scores( node_sequence_count );
      ScoreType score =
        m_dynamic_programming.forward_score( 
          m_parameters,
          (
            ( node_id == 0 ) ?
            *( m_profileTree->getProfileTreeRoot() ) :
            m_profileTree->getProfileTreeInternalNode( node_id )
          ),
          sequences,
          node_sequence_count,
          forward_matrices,
          &sequence_scores,
          NULL
        );
      if( be_extra_verbose ) {
        cout << "Total sequence score is " << score << endl;
        //cout << "Sequence scores are " << endl;
        //for( uint32_t seq_i = 0; seq_i < node_sequence_count; seq_i++ ) {
        //  cout << sequence_scores[ seq_i ] << endl; 
        //}
        //cout << endl;
      } // End if be_extra_verbose

      // Allocate some temporary backward matrices...
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::RowVector backward_rows_1(
        sequences,
        node_sequence_count
      );
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::RowVector backward_rows_2(
        sequences,
        node_sequence_count
      );
      if( be_extra_verbose ) {
        cout << "Calculating alignment profiles." << endl;
      } // End if be_extra_verbose
      m_dynamic_programming.calculateAlignmentProfiles(
        m_parameters,
        (
          ( node_id == 0 ) ?
          *( m_profileTree->getProfileTreeRoot() ) :
          m_profileTree->getProfileTreeInternalNode( node_id )
        ),
        sequences,
        node_sequence_count,
        NULL,         //&sequence_scores,
        forward_matrices,         //trainer.m_forward_matrices,
        backward_rows_1,
        backward_rows_2,
        aps
      );

      if( be_extra_verbose ) {
        cout << "Unscaling the alignment profiles." << endl;
      } // End if be_extra_verbose
      // Now unscale them.
      for( uint32_t seq_i = 0; seq_i < node_sequence_count; seq_i++ ) {
        aps[ seq_i ].unscale();
      }
      if( be_extra_verbose ) {
        cout << "Returning from calculateAlignmentProfiles(..)." << endl;
      } // End if be_extra_verbose

      return;
    } // calculateAlignmentProfiles( uint32_t const &, vector<AlignmentProfile> & ) const

    // COOCExpectedCountMetric
    void
    calculateDistanceMatrix (
      uint32_t const & node_id,
      vector<typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile> & aps,
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DistanceMatrix<MatrixValueType> & cooc_matrix,
      Tag<COOCExpectedCountMetric_> const tag
    ) const {
      static const bool be_extra_verbose = true;

      // Cooccurrence expected count
      // The second argument sez count cooccurrence of deletions, too.
      cooc_matrix.setToMatchEmissionCooccurrenceExpectedCounts( aps, true );
      if( be_extra_verbose ) {
        cout << "After setting using setToMatchEmissionCooccurrenceExpectedCounts( aps ), cooc_matrix is " << endl;
        cout << cooc_matrix << endl;
      }

      assert( cooc_matrix.isProximityMatrix() );
      return;
    } // calculateDistanceMatrix( uint32_t const &, vector<AlignmentProfile> &, DistanceMatrix<MatrixValueType> &, Tag<COOCExpectedCountMetric_> const ) const

    // COOCProbabilityMetric
    void
    calculateDistanceMatrix (
      uint32_t const & node_id,
      vector<typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile> & aps,
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DistanceMatrix<MatrixValueType> & cooc_matrix,
      Tag<COOCProbabilityMetric_> const tag
    ) const {
      static const bool be_extra_verbose = true;

      uint32_t const node_sequence_count =
        m_profile_sequenceIndices[ node_id ].size();

      // Cooccurrence probability
      // The second argument sez include cooccurrence of deletions, too.
      cooc_matrix.setToMatchEmissionCooccurrenceProbabilities( aps, true );
      cout << "After setting using setToMatchEmissionCooccurrenceProbabilities( aps ), cooc_matrix is " << endl;
      cout << cooc_matrix << endl;

      assert( cooc_matrix.isProximityMatrix() );
      return;
    } // calculateDistanceMatrix( uint32_t const &, vector<AlignmentProfile> &, DistanceMatrix<MatrixValueType> &, Tag<COOCProbabilityMetric_> const ) const

    // CrossEntropyMetric
    void
    calculateDistanceMatrix (
      uint32_t const & node_id,
      vector<typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile> & aps,
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DistanceMatrix<MatrixValueType> & entropy_matrix,
      Tag<CrossEntropyMetric_> const tag
    ) const {
      static const bool be_extra_verbose = true;

      // Note that the cross-entropy metric requires that the alignment profile
      // be normalized first.  Thus we can't downweight by the probability of
      // deletion.  Instead we can use entropy weights.
      // TODO: Decide: Use weights or no?
      static const bool use_entropy_weights = true;

      uint32_t node_length;
      if( node_id == 0 ) {
        node_length = m_profileTree->getProfileTreeRoot()->length();
      } else {
        node_length = m_profileTree->getProfileTreeInternalNode( node_id ).length();
      }

      if( use_entropy_weights ) {
        // We'll count only the MatchEmissionParameters, weighed by the
        // corresponding deletion fraction.
        typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile entropy_weights(
          (
            node_length + 1
          )
        );

        // TODO: Right now it's not actually weighing by the
        // deletion fraction, which it really ought to.  For this we must be
        // allowed to pass different weights for each alignment profile.

        entropy_weights.zero();
        // Note that we leave the last and first positions with a zero weight,
        // since there's no matches there.
        for( uint32_t pos_i = 1; pos_i < entropy_weights.length(); pos_i++ ) {
          entropy_weights[ pos_i ][ Emission::Match ] = 1.0;
        } // End foreach pos_i, set the weights to 1.0.

        // We *must* normalize first, or else the SKL distances have no
        // meaning.
        for( uint32_t seq_i = 0; seq_i < aps.size(); seq_i++ ) {
          aps[ seq_i ].normalize( 0 );
        }

        entropy_matrix.setToCrossEntropies( aps, &entropy_weights );
      } else { // if use_entropy_weights .. else ..
        // We *must* normalize first, or else the SKL distances have no
        // meaning.
        for( uint32_t seq_i = 0; seq_i < aps.size(); seq_i++ ) {
          aps[ seq_i ].normalize( 0 );
        }

        entropy_matrix.setToCrossEntropies( aps );
      } // End if use_entropy_weights

      if( be_extra_verbose ) {
        cout << "After setting using setToCrossEntropies( aps ), entropy_matrix is " << endl;
        cout << entropy_matrix << endl;
      } // be_extra_verbose

      assert( !entropy_matrix.isProximityMatrix() );
      return;
    } // calculateDistanceMatrix( uint32_t const &, vector<AlignmentProfile> &, DistanceMatrix<MatrixValueType> &, Tag<CrossEntropyMetric_> const ) const

    // SKLMetric
    void
    calculateDistanceMatrix (
      uint32_t const & node_id,
      vector<typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile> & aps,
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DistanceMatrix<MatrixValueType> & skl_matrix,
      Tag<SKLMetric_> const tag
    ) const {
      static const bool be_extra_verbose = true;

      // Note that this metric requires that the alignment profile
      // be normalized first.  Thus we can't downweight by the probability of
      // deletion.  Instead we can use entropy weights.
      // TODO: Decide: Use weights or no?
      static const bool use_entropy_weights = true;

      uint32_t node_length;
      if( node_id == 0 ) {
        node_length = m_profileTree->getProfileTreeRoot()->length();
      } else {
        node_length = m_profileTree->getProfileTreeInternalNode( node_id ).length();
      }

      if( use_entropy_weights ) {
        // We'll count only the MatchEmissionParameters, weighed by the
        // corresponding deletion fraction.
        typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile entropy_weights(
          (
            node_length + 1
          )
        );

        // TODO: Right now it's not actually weighing by the
        // deletion fraction, which it really ought to.  For this we must be
        // allowed to pass different weights for each alignment profile.

        entropy_weights.zero();
        // Note that we leave the last and first positions with a zero weight,
        // since there's no matches there.
        for( uint32_t pos_i = 1; pos_i < entropy_weights.length(); pos_i++ ) {
          entropy_weights[ pos_i ][ Emission::Match ] = 1.0;
        } // End foreach pos_i, set the weights to 1.0.

        // We *must* normalize first, or else the SKL distances have no
        // meaning.
        for( uint32_t seq_i = 0; seq_i < aps.size(); seq_i++ ) {
          aps[ seq_i ].normalize( 0 );
        }

        skl_matrix.setToSymmeterizedKullbackLeiblerDivergences( aps, &entropy_weights );
      } else { // if use_entropy_weights .. else ..
        // We *must* normalize first, or else the SKL distances have no
        // meaning.
        for( uint32_t seq_i = 0; seq_i < aps.size(); seq_i++ ) {
          aps[ seq_i ].normalize( 0 );
        }

        skl_matrix.setToSymmeterizedKullbackLeiblerDivergences( aps );
      } // End if use_entropy_weights

      if( be_extra_verbose ) {
        cout << "After setting using setToSymmeterizedKullbackLeiblerDivergences( aps ), skl_matrix is " << endl;
        cout << skl_matrix << endl;
      } // End if be_extra_verbose

      assert( !skl_matrix.isProximityMatrix() );
      return;
    } // calculateDistanceMatrix( uint32_t const &, vector<AlignmentProfile> &, DistanceMatrix<MatrixValueType> &, Tag<SKLMetric_> const ) const

#ifdef __HAVE_MUSCLE
    void
    createMuscleTree (
      uint32_t const & node_id,
      Tree & muscle_tree
    )
    {
      // TODO: REMOVE?
      static const bool show_profiles = false;
      static const bool be_extra_verbose = true; // false;
      static const bool write_tree_to_file = true; // false;

      if( node_id != 0 ) {
        // TODO: Why is this necessary?!
        // I believe one problem is that we store the nodes in a vector in
        // ProfileTree, so when we grow it, all the data is copied -- then I
        // dunno why the default copy can't get the m_root pointers correct (it
        // seems to be doing a memcpy).
        // TODO: FIX!
        m_profileTree->getProfileTreeInternalNode( node_id ).ensurePositionsKnowTheirRoot();
      }
      uint32_t const node_sequence_count =
        m_profile_sequenceIndices[ node_id ].size();

      uint32_t node_length;
      if( node_id == 0 ) {
        node_length = m_profileTree->getProfileTreeRoot()->length();
      } else {
        node_length = m_profileTree->getProfileTreeInternalNode( node_id ).length();
      }

      if( be_extra_verbose ) {
        cout << "This node's length is: " << node_length << endl;
        cout << "This node's sequences are: " << endl;
        for( uint32_t seq_i = 0; seq_i < node_sequence_count; seq_i++ ) {
          cout << "> " << ( static_cast<Fasta<SequenceResidueType> const *>( &m_sequences ) )->m_descriptions[ m_profile_sequenceIndices[ node_id ][ seq_i ] ] << " (index " << m_profile_sequenceIndices[ node_id ][ seq_i ] << ")" << endl;
        }
      } // End if be_extra_verbose

      if( be_extra_verbose && show_profiles ) {
        cout << "The node is: " << endl;
        if( node_id == 0 ) {
          cout << *m_profileTree->getProfileTreeRoot();
        } else {
          cout << m_profileTree->getProfileTreeInternalNode( node_id );
        }
      } // End if be_extra_verbose

      // Calculate a distance matrix, based on the alignment profiles.
      vector<typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile> aps;
      calculateAlignmentProfiles( node_id, aps );

      // Calculate a distance matrix based on the alignment profiles.  Note
      // that the alignment profiles may be weighted and/or normalized first
      // (altering them), depending on the method.  Also the matrix might be a
      // proximity matrix (test with dist_or_prox_matrix.isProximityMatrix()).
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DistanceMatrix<MatrixValueType> dist_or_prox_matrix;
      // TODO: This class should accept the metric as an argument
      calculateDistanceMatrix( node_id, aps, dist_or_prox_matrix, Tag<COOCExpectedCountMetric_>() );
      //calculateDistanceMatrix( node_id, aps, dist_or_prox_matrix, Tag<COOCProbabilityMetric_>() );
      //calculateDistanceMatrix( node_id, aps, dist_or_prox_matrix, Tag<CrossEntropyMetric_>() );
      //calculateDistanceMatrix( node_id, aps, dist_or_prox_matrix, Tag<SKLMetric_>() );

      //#ifdef __HAVE_MUSCLE
      DistFunc muscle_dist_func;
      muscle_dist_func.SetCount( node_sequence_count );
      for( uint32_t seq_i = 0; seq_i < node_sequence_count; seq_i++ ) {
        muscle_dist_func.SetName( seq_i, ( static_cast<Fasta<SequenceResidueType> const *>( &m_sequences ) )->m_descriptions[ m_profile_sequenceIndices[ node_id ][ seq_i ] ].c_str() );
        muscle_dist_func.SetId( seq_i, seq_i );
      } // End foreach seq_i
        
      float prox;
      float dist;

      if( be_extra_verbose ) {
        cout << "Creating muscle dist func." << endl;
      }
      if( dist_or_prox_matrix.isProximityMatrix() ) {
        MatrixValueType max_prox = dist_or_prox_matrix.maximumValue();
        for( uint32_t seq_i = 0; seq_i < node_sequence_count; seq_i++ ) {
          for( uint32_t other_seq_i = 0; other_seq_i < seq_i; other_seq_i++ ) {
            // TODO: Dehackify!!  Trying to do 0 - x, but negative numbers aren't ok.  Should do maxValue - x.  aps[ 0 ].length() is the maximum theoretical value that the COOCExpectedCount metric can take (if each position had expected count = 1).  The problem is that floating point numbers don't represent the tiny differences among these values very well.
            // TODO: ERE I AM.  The cooc matrix itself contains 0s.  Something is lost along the way.. ?!
            muscle_dist_func.SetDist( seq_i, other_seq_i, toDouble( max_prox - dist_or_prox_matrix( seq_i, other_seq_i ) ) );
            muscle_dist_func.SetDist( other_seq_i, seq_i, toDouble( max_prox - dist_or_prox_matrix( other_seq_i, seq_i ) ) );
            //if( dist_or_prox_matrix( seq_i, other_seq_i ) == 0.0 ) {
            //  muscle_dist_func.SetDist( seq_i, other_seq_i, pow( numeric_limits<float>::max(), .25F ) );
            //} else {
            //  muscle_dist_func.SetDist( seq_i, other_seq_i, ( float )toDouble( 1.0 / dist_or_prox_matrix( seq_i, other_seq_i ) ) );
            //}
            //if( prox == dist_or_prox_matrix( other_seq_i, seq_i ) ) {
            //  muscle_dist_func.SetDist( other_seq_i, seq_i, pow( numeric_limits<float>::max(), .5F ) );
            //} else {
            //  muscle_dist_func.SetDist( other_seq_i, seq_i, ( float )toDouble( 1.0 / dist_or_prox_matrix( other_seq_i, seq_i ) ) );
            //}
          } // End foreach other_seq_i
        } // End foreach seq_i
      } else { 
        // It's a distance matrix
        for( uint32_t seq_i = 0; seq_i < node_sequence_count; seq_i++ ) {
          for( uint32_t other_seq_i = 0; other_seq_i < seq_i; other_seq_i++ ) {
            dist = ( float )toDouble( dist_or_prox_matrix( seq_i, other_seq_i ) );
            muscle_dist_func.SetDist( seq_i, other_seq_i, dist );
            dist = ( float )toDouble( dist_or_prox_matrix( other_seq_i, seq_i ) );
            muscle_dist_func.SetDist( other_seq_i, seq_i, dist );
          } // End foreach other_seq_i
        } // End foreach seq_i
      } // End if dist_or_prox_matrix.isProximityMatrix() .. else ..
      if( be_extra_verbose ) {
        cout << "done." << endl;
      } // End if be_extra_verbose

      if( be_extra_verbose ) {
        cout << "Creating muscle_clust_set." << endl;
      } // End if be_extra_verbose
      ClustSetDF muscle_clust_set( muscle_dist_func );
      // TODO: REMOVE
      assert( muscle_clust_set.GetLeafCount() == node_sequence_count );
      if( be_extra_verbose ) {
        cout << "done." << endl;
      } // End if be_extra_verbose
      if( be_extra_verbose ) {
        cout << "Creating muscle_clust." << endl;
      } // End if be_extra_verbose
      Clust muscle_clust;
      //muscle_clust.Create( muscle_clust_set, CLUSTER_UPGMB );
      //muscle_clust.Create( muscle_clust_set, CLUSTER_UPGMA );
      //muscle_clust.Create( muscle_clust_set, CLUSTER_UPGMAMin );
      //muscle_clust.Create( muscle_clust_set, CLUSTER_UPGMAMax );
      muscle_clust.Create( muscle_clust_set, CLUSTER_NeighborJoining );
      if( be_extra_verbose ) {
        cout << "done." << endl;
      } // End if be_extra_verbose

      if( be_extra_verbose ) {
        cout << "Creating muscle_tree." << endl;
      } // End if be_extra_verbose
      muscle_tree.FromClust( muscle_clust );
      if( be_extra_verbose ) {
        cout << "done." << endl;
      } // End if be_extra_verbose

      // Fix the root first
      if( be_extra_verbose ) {
        //cout << "Fixing root using ROOT_MinAvgLeafDist" << endl; // TODO: See below -- debug (segfaults)
        cout << "Fixing root using ROOT_MidLongestSpan" << endl;
        if( true ) {
          cout << "Before fixing it, root is " << muscle_tree.GetRootNodeIndex() << endl;
          cout << "There are " << muscle_tree.GetNodeCount() << " nodes in the tree." << endl;
        }
      } // End if be_extra_verbose
      //FixRoot( muscle_tree, ROOT_MinAvgLeafDist ); // TODO: Causes an error in phy3.cpp: "RootByMinAvgLeafDist, internal error 6"
      FixRoot( muscle_tree, ROOT_MidLongestSpan );

      const unsigned uRootNodeIndex = muscle_tree.GetRootNodeIndex();
      if( be_extra_verbose ) {
        cout << "done." << endl;

        // TODO: REMOVE
        if( true ) {
          cout << "uRootNodeIndex is " << uRootNodeIndex << endl;
          cout << "m_profile_sequenceIndices[ " << node_id << " ] are ( ";
          for( uint32_t i = 0; i < m_profile_sequenceIndices[ node_id ].size(); i++ ) {
            if( i > 0 ) {
              cout << ", ";
            }
            cout << m_profile_sequenceIndices[ node_id ][ i ];
          }
          cout << " )" << endl;
        } // End if true

      } // End if be_extra_verbose

      if( write_tree_to_file ) {
        if( be_extra_verbose ) {
          cout << "Writing tree to file test_muscle_tree.dnd" << endl;
        } // End if be_extra_verbose
        TextFile f( "test_muscle_tree.dnd", true );
        muscle_tree.ToFile( f );
        if( be_extra_verbose ) {
          cout << "done." << endl;
        } // End if be_extra_verbose
      } // End if write_tree_to_file

      return;
    } // createMuscleTree( uint32_t const &, Tree & )
#endif // __HAVE_MUSCLE

    /**
     * Given a (reference to a) vector of sequences, train the given profile
     * (reference), and return the score.  If the vector has length 1,
     * delegates to createProfileFromSequence(..).  Otherwise the profile will
     * first be initialized as a copy of its parent in the tree, and then it
     * will be trained with for at least 2 iterations.  Note that this will
     * increment m_totalIterations with the number of training iterations.
     */
    ScoreType
    createProfileFromSequences (
      vector<SequenceType> const & sequences,
      InternalNodeType & profile
    ) {
      static const bool be_extra_verbose = true;

      assert( sequences.size() != 0 );
      if( sequences.size() == 1 ) {
        return createProfileFromSequence( sequences[ 0 ], profile );
      }

      // copyFrom() clobbers profileTreeVertex, so we save & restore it.
      uint32_t profile_vertex = profile.getProfileTreeVertex();
      // TODO: REMOVE
      //cout << "[createProfileFromSequences] Saved profile tree vertex: " << profile_vertex << endl;
      profile.copyFrom( m_profileTree->getParent( profile_vertex ) );
      ProfileTrainer<InternalNodeType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType> trainer(
        &profile,
        sequences,
        sequences.size()
      );
      // Use our own special parameters
      trainer.m_parameters = m_parameters;
  
      if( be_extra_verbose ) {
        trainer.m_parameters.verbosity = VERBOSITY_Low;
      } else {
        trainer.m_parameters.verbosity = VERBOSITY_None;
      }
  
      // Let's also NOT train the profile globals.
      // TODO: ?
      //trainer.m_parameters.trainProfileGlobals = false;
  
      // Make sure it gets to adjust the length.
      trainer.m_parameters.minIterations = 2;
  
      // Do it.
      ScoreType score = trainer.train();
      m_totalIterations += trainer.m_iteration;

      profile.setProfileTreeVertex( profile_vertex );

      return score;
    } // createProfileFromSequences( vector<SequenceType> const &, InternalNodeType & )

    /**
     * Given a (reference to a) sequence, reinitialize the given profile
     * (reference) to represent that sequence exactly, and return 1.0.
     */
    ScoreType
    createProfileFromSequence (
      SequenceType const & sequence,
      InternalNodeType & profile
    ) const {
      // fromSequence(..) clobbers the vertex, so we save & restore it.
      uint32_t profile_vertex = profile.getProfileTreeVertex();
      // TODO: REMOVE
      //cout << "[createProfileFromSequence] Saved profile tree vertex: " << profile_vertex << endl;
      profile.fromSequence( sequence );
      profile.setProfileTreeVertex( profile_vertex );
      return 1.0;
    } // createProfileFromSequence( SequenceType const &, InternalNodeType & ) const

#ifdef __HAVE_MUSCLE
    void
    splitNode (
      uint32_t const & node_id,
      Tree const & muscle_tree,
      unsigned const & muscle_node_id,
      bool const & retrain_children
    )
    {
      // TODO: REMOVE?
      static const bool be_extra_verbose = true;

      uint32_t const node_sequence_count =
        m_profile_sequenceIndices[ node_id ].size();

      if( node_sequence_count == 1 ) {
        if( be_extra_verbose ) {
          cout << "Not splitting node " << node_id << ", since there is only one sequence assigned to it." << endl;
        }
        return;
      } // End if node_sequence_count == 1

      uint32_t node_length;
      if( node_id == 0 ) {
        node_length = m_profileTree->getProfileTreeRoot()->length();
      } else {
        node_length = m_profileTree->getProfileTreeInternalNode( node_id ).length();
      }

      if( be_extra_verbose ) {
        cout << "Splitting based on the muscle tree node " << muscle_node_id << ", which corresponds to our tree node " << node_id << "." << endl;
        cout << "\tThere are " << node_sequence_count << " sequences assigned to this node." << endl;
      } // End if be_extra_verbose

      unsigned *left_leaves = new unsigned[node_sequence_count];
      unsigned *right_leaves = new unsigned[node_sequence_count];
      unsigned left_leaves_count;
      unsigned right_leaves_count;

      unsigned muscle_left_id =
        muscle_tree.GetLeft( muscle_node_id );
      assert( NULL_NEIGHBOR != muscle_left_id );
      unsigned muscle_right_id =
        muscle_tree.GetRight( muscle_node_id );
      assert( NULL_NEIGHBOR != muscle_right_id );
    
      GetLeavesExcluding(muscle_tree, muscle_left_id, muscle_right_id, left_leaves, &left_leaves_count);
      GetLeavesExcluding(muscle_tree, muscle_right_id, muscle_left_id,
      right_leaves, &right_leaves_count);
      if( be_extra_verbose ) {
        cout << endl;
        cout << "Split at muscle node " << muscle_node_id << ':' << endl;
        cout << "left (muscle node " << muscle_left_id << ") = {" << endl;
        for (unsigned n = 0; n < left_leaves_count; ++n) {
          if( retrain_children ) {
            cout << '\t' << m_profile_sequenceIndices[ node_id ][ left_leaves[ n ] ] << '\t' << muscle_tree.GetName( left_leaves[ n ] ) << endl;
          } else {
            // indices are relative to the root node
            cout << '\t' << m_profile_sequenceIndices[ 0 ][ left_leaves[ n ] ] << '\t' << muscle_tree.GetName( left_leaves[ n ] ) << endl;
          }
        }
        cout << " }" << endl;
        cout << "right (muscle node " << muscle_right_id << ") = {" << endl;
        for (unsigned n = 0; n < right_leaves_count; ++n) {
          if( retrain_children ) {
            cout << '\t' << m_profile_sequenceIndices[ node_id ][ right_leaves[ n ] ] << '\t' << muscle_tree.GetName( right_leaves[ n ] ) << endl;
          } else {
            // indices are relative to the root node
            cout << '\t' << m_profile_sequenceIndices[ 0 ][ right_leaves[ n ] ] << '\t' << muscle_tree.GetName( right_leaves[ n ] ) << endl;
          }
        }
        cout << " }" << endl;
        cout << endl;
      } // End if be_extra_verbose

      // Try splitting the profile again...
      // For scoring, we'll first need to divide the sequences into right and
      // left groups, build profiles separately (starting from where we'd left
      // off), so that we can recalculate the score at the end as the product..
      if( be_extra_verbose ) {
        cout << "Assigning sequences to left and right profiles:" << endl;
      } // End if be_extra_verbose
      vector<uint32_t> right_sequence_indices;
      vector<uint32_t> left_sequence_indices;

      // TODO: Inefficient.  Use a list?
      left_sequence_indices.resize( left_leaves_count );
      for (unsigned n = 0; n < left_leaves_count; ++n) {
        if( retrain_children ) {
          left_sequence_indices[ n ] =
            m_profile_sequenceIndices[ node_id ][ left_leaves[ n ] ];
        } else {
          // indices are relative to the root node
          left_sequence_indices[ n ] =
            m_profile_sequenceIndices[ 0 ][ left_leaves[ n ] ];
        }
      }
      // TODO: Inefficient.  Use a list?
      right_sequence_indices.resize( right_leaves_count );
      for (unsigned n = 0; n < right_leaves_count; ++n) {
        if( retrain_children ) {
          right_sequence_indices[ n ] =
            m_profile_sequenceIndices[ node_id ][ right_leaves[ n ] ];
        } else {
          // indices are relative to the root node
          right_sequence_indices[ n ] =
            m_profile_sequenceIndices[ 0 ][ right_leaves[ n ] ];
        }
      }

      if( be_extra_verbose ) {
        cout << "There are " << right_sequence_indices.size() << " sequences assigned to the right profile." << endl;
        cout << "There are " << left_sequence_indices.size() << " sequences assigned to the left profile." << endl;
      } // End if be_extra_verbose

      // Initialize the vector of vectors of sequence indices
      // Set up the profile tree with what we learned.
      m_profileTree->addChild( node_id ); // right 
      uint32_t right_child_vertex = 
        m_profileTree->getChildVertex( node_id, 1 );

      m_profileTree->addChild( node_id ); // left
      uint32_t left_child_vertex = 
        m_profileTree->getChildVertex( node_id, 2 );

      // TODO: REMOVE
      assert( m_profileTree->getChild( node_id, 1 ).getProfileTreeVertex() == right_child_vertex );
      assert( m_profileTree->getChild( node_id, 2 ).getProfileTreeVertex() == left_child_vertex );

      if( be_extra_verbose ) {
        cout << "Right child vertex: " << right_child_vertex << endl;
        cout << "Left child vertex: " << left_child_vertex << endl;
      }

      // TODO: This is expensive, potentially.  Use a list?
      m_profile_sequenceIndices.resize( m_profile_sequenceIndices.size() + 2 );
      assert( right_child_vertex == ( m_profile_sequenceIndices.size() - 2 ) );
      assert( left_child_vertex == ( m_profile_sequenceIndices.size() - 1 ) );
      m_profile_sequenceIndices[ right_child_vertex ] = right_sequence_indices;
      m_profile_sequenceIndices[ left_child_vertex ] = left_sequence_indices;

      // TODO: This is expensive, potentially.  Use a list?
      m_profile_scores.resize( m_profile_scores.size() + 2 );
      assert( right_child_vertex == ( m_profile_scores.size() - 2 ) );
      assert( left_child_vertex == ( m_profile_scores.size() - 1 ) );
      // Note that if retrain_children is true, this will be changed later.
      m_profile_scores[ left_child_vertex ] = 
        m_profile_scores[ right_child_vertex ] = m_profile_scores[ node_id ];

      // TODO: This is expensive.  Use a stringbuffer?
      string right_name;
      if( m_profile_sequenceIndices[  right_child_vertex ].size() > 1 ) {
        right_name = "[";
      }
      right_name +=
        boost::lexical_cast<std::string>(
          m_profile_sequenceIndices[ right_child_vertex ][ 0 ]
        );
      for( uint32_t seq_i = 1; seq_i < m_profile_sequenceIndices[  right_child_vertex ].size(); seq_i++ ) {
        right_name += ",";
        right_name +=
          boost::lexical_cast<std::string>(
            m_profile_sequenceIndices[ right_child_vertex ][ seq_i ]
          );
      }
      if( m_profile_sequenceIndices[  right_child_vertex ].size() > 1 ) {
        right_name += "]";
      }
      m_profileTree->setName(
        ( right_child_vertex ),
        right_name
      );
      string left_name;
      if( m_profile_sequenceIndices[  left_child_vertex ].size() > 1 ) {
        left_name = "[";
      }
      left_name +=
        boost::lexical_cast<std::string>(
          m_profile_sequenceIndices[ left_child_vertex ][ 0 ]
        );
      for( uint32_t seq_i = 1; seq_i < m_profile_sequenceIndices[  left_child_vertex ].size(); seq_i++ ) {
        left_name += ", ";
        left_name +=
          boost::lexical_cast<std::string>(
            m_profile_sequenceIndices[ left_child_vertex ][ seq_i ]
          );
      }
      if( m_profile_sequenceIndices[  left_child_vertex ].size() > 1 ) {
        left_name += "]";
      }
      m_profileTree->setName(
        ( left_child_vertex ),
        left_name
      );

      delete[] left_leaves;
      delete[] right_leaves;

      return;
    } // splitNode( uint32_t const &, Tree const &, bool const & )
#endif // __HAVE_MUSCLE
  }; // End class ProfileTreeTrainer

  //======//// potentially non-inline implementations ////========//

  ////// Class galosh::ProfileTreeTrainer::Parameters ////
  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_INIT
  ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
        Parameters ()
      {
        if( DEFAULT_debug >= DEBUG_All ) {
          cout << "[debug] ProfileTreeTrainer::Parameters::<init>()" << endl;
        } // End if DEBUG_All
        resetToDefaults();
      } // <init>()

      // Copy constructor
  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
      template <class AnyParameters>
  GALOSH_INLINE_INIT
  ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      Parameters ( const AnyParameters & copy_from )
      {
        if( copy_from.debug >= DEBUG_All ) {
          cout << "[debug] ProfileTreeTrainer::Parameters::<init>( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
      } // <init>( AnyParameters const & )

      // Copy constructor/operator
  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
      template <class AnyParameters>
  GALOSH_INLINE_TRIVIAL
  typename ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters &
  ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      operator= (
        const AnyParameters & copy_from
      )
      {
        if( copy_from.debug >= DEBUG_All ) {
          cout << "[debug] ProfileTreeTrainer::Parameters::operator=( copy_from )" << endl;
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
  ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      copyFromNonVirtual (
        AnyParameters const & copy_from
      )
      {
        ProfileTrainer<ProfileTreeRoot<ResidueType, ProbabilityType>, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::Parameters::copyFromNonVirtual( copy_from );
        //if( copy_from.debug >= DEBUG_All ) {
        //  cout << "[debug] ProfileTreeTrainer::Parameters::copyFromNonVirtual( copy_from )" << endl;
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
  ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      copyFromNonVirtualDontDelegate (
        AnyParameters const & copy_from
      )
      {
        shareProfilePositions = copy_from.shareProfilePositions;
        shareProfilePositions_percentChangeTo0Threshold = copy_from.shareProfilePositions_percentChangeTo0Threshold;
        childSequenceMixtureThreshold = copy_from.childSequenceMixtureThreshold;
      } // copyFromNonVirtualDontDelegate( AnyParameters const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_TRIVIAL
  void
  ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
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
  ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      resetToDefaults ()
      {
        ProfileTrainer<ProfileTreeRoot<ResidueType, ProbabilityType>, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::Parameters::resetToDefaults();
        // TODO: Why isn't the compiler finding "debug" in galosh::Parameters?
        //if( debug >= DEBUG_All ) {
        //  cout << "[debug] ProfileTreeTrainer::Parameters::resetToDefaults()" << endl;
        //} // End if DEBUG_All
        shareProfilePositions = DEFAULT_shareProfilePositions;
        shareProfilePositions_percentChangeTo0Threshold = DEFAULT_shareProfilePositions_percentChangeTo0Threshold;
        childSequenceMixtureThreshold = DEFAULT_childSequenceMixtureThreshold;
      } // resetToDefaults()

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
      template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
      void
  ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      writeParameters (
        std::basic_ostream<CharT,Traits>& os
      ) const
      {
        //os << static_cast<typename ProfileTrainer<ProfileTreeRoot<ResidueType, ProbabilityType>, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::Parameters>( parameters ) << endl;
        ProfileTrainer<ProfileTreeRoot<ResidueType, ProbabilityType>, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::Parameters::writeParameters( os );
        os << endl;
        
        os << "[ProfileTreeTrainer]" << endl;
        os << "shareProfilePositions = " << shareProfilePositions << endl;
        os << "shareProfilePositions_percentChangeTo0Threshold = " << shareProfilePositions_percentChangeTo0Threshold << endl;
        os << "childSequenceMixtureThreshold = " << childSequenceMixtureThreshold << endl;
      } // writeParameters( basic_ostream & ) const

  ////// Class galosh::ProfileTreeTrainer::ParametersModifierTemplate ////
  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
    template <class ParametersType>
  GALOSH_INLINE_INIT
  ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      ParametersModifierTemplate ()
      {
        if( base_parameters_modifier_t::parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProfileTreeTrainer::ParametersModifierTemplate::<init>()" << endl;
        } // End if DEBUG_All
        isModified_reset();
      } // <init>()

      // Copy constructor
  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
    template <class ParametersType>
      template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_INIT
  ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      ParametersModifierTemplate ( const AnyParametersModifierTemplate & copy_from )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProfileTreeTrainer::ParametersModifierTemplate::<init>( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
      } // <init>( AnyParametersModifierTemplate const & )

      // Copy constructor/operator
  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
    template <class ParametersType>
      template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_TRIVIAL
  typename ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::template ParametersModifierTemplate<ParametersType> &
  ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      operator= (
        const AnyParametersModifierTemplate & copy_from
      )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProfileTreeTrainer::ParametersModifierTemplate::operator=( copy_from )" << endl;
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
  GALOSH_INLINE_TRIVIAL
      void
  ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProfileTreeTrainer::ParametersModifierTemplate::copyFromNonVirtual( copy_from )" << endl;
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
  ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      isModified_copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      )
      {
        base_parameters_modifier_t::isModified_copyFromNonVirtual( copy_from );
        isModified_shareProfilePositions = copy_from.isModified_shareProfilePositions;
        isModified_shareProfilePositions_percentChangeTo0Threshold = copy_from.isModified_shareProfilePositions_percentChangeTo0Threshold;
        isModified_childSequenceMixtureThreshold = copy_from.isModified_childSequenceMixtureThreshold;
      } // isModified_copyFromNonVirtual( AnyParametersModifierTemplate const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
    template <class ParametersType>
  GALOSH_INLINE_TRIVIAL
      void
  ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
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
  ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      isModified_reset ()
      {
        base_parameters_modifier_t::isModified_reset();
        isModified_shareProfilePositions = false;
        isModified_shareProfilePositions_percentChangeTo0Threshold = false;
        isModified_childSequenceMixtureThreshold = false;
      } // isModified_reset()

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
    template <class ParametersType>
      template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
      void
  ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      writeParametersModifier (
        std::basic_ostream<CharT,Traits>& os
      )
      {
        //base_parameters_modifier_t::operator<<( os, parameters_modifier );
        base_parameters_modifier_t::writeParametersModifier( os );
        os << endl;

        os << "[ProfileTreeTrainer]" << endl;
        if( isModified_shareProfilePositions ) {
          os << "shareProfilePositions = " << base_parameters_modifier_t::parameters.shareProfilePositions << endl;
        }
        if( isModified_shareProfilePositions_percentChangeTo0Threshold ) {
          os << "shareProfilePositions_percentChangeTo0Threshold = " << base_parameters_modifier_t::parameters.shareProfilePositions_percentChangeTo0Threshold << endl;
        }
        if( isModified_childSequenceMixtureThreshold ) {
          os << "childSequenceMixtureThreshold = " << base_parameters_modifier_t::parameters.childSequenceMixtureThreshold << endl;
        }
      } // writeParametersModifier( basic_ostream & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
    template <class ParametersType>
      template<class AnyParameters>
  GALOSH_INLINE_PARAMETERSMODIFIER_APPLY_MODIFICATIONS
      void
  ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      applyModifications ( AnyParameters & target_parameters )
      {
        base_parameters_modifier_t::applyModifications( target_parameters );

        if( isModified_shareProfilePositions ) {
          target_parameters.shareProfilePositions =
            base_parameters_modifier_t::parameters.shareProfilePositions;
        }
        if( isModified_shareProfilePositions_percentChangeTo0Threshold ) {
          target_parameters.shareProfilePositions_percentChangeTo0Threshold =
            base_parameters_modifier_t::parameters.shareProfilePositions_percentChangeTo0Threshold;
        }
        if( isModified_childSequenceMixtureThreshold ) {
          target_parameters.childSequenceMixtureThreshold =
            base_parameters_modifier_t::parameters.childSequenceMixtureThreshold;
        }
      } // applyModifications( AnyParameters & )

} // End namespace galosh

#endif // __GALOSH_PROFILETREETRAINER_HPP__
