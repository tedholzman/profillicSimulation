/*---------------------------------------------------------------------------##
##  File:
##      @(#) ProfileGibbs.hpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      Class definition for the Galosh Profile HMM Gibbs Sampler class.
##
#******************************************************************************
#*  Copyright (C) 2008    Paul Edlefsen                                       *
#*  All rights reserved.                                                      *
#*****************************************************************************/

#if     _MSC_VER > 1000
#pragma once
#endif

#ifndef __GALOSH_PROFILEGIBBS_HPP__
#define __GALOSH_PROFILEGIBBS_HPP__

#include "Galosh.hpp"

#include "Parameters.hpp"
using galosh::Parameters;
using galosh::DebugLevel;
using galosh::VerbosityLevel;

#include "Profile.hpp"
using galosh::ProfileTreeRoot;
using galosh::ProfileTreeInternalNode;

#include "ProfileTreeTrainer.hpp"
using galosh::ProfileTreeTrainer;

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

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

#ifndef HMMOC_BFLOAT_ALGEBRA_HPP
#include "Algebra.hpp"           //necessary for doublerealspace defn TAH
#endif

namespace galosh {

  // TODO: Move somewhere else?
  double
  calculateAutocorrelation (
    vector<double> const & data,
    double const mean_of_highest_indices,
    double const mean_of_lowest_indices,
    double const stdev_of_highest_indices,
    double const stdev_of_lowest_indices,
    uint32_t const first_index = 0,
    uint32_t last_index = 0, // 0 means use ( data.size() - 1 )
    uint32_t const lag = 1
  )
  {
    if( data.size() == 0 ) {
      return 0;
    }
    if( last_index == 0 ) {
      last_index = ( data.size() - 1 );
    }
    double covariance = 0.0;
    for( uint32_t high_i = ( first_index + lag ); high_i <= last_index; high_i++ ) {
      covariance +=
        (
          ( data[ high_i ] - mean_of_highest_indices ) *
          ( data[ high_i - lag ] - mean_of_lowest_indices )
        );
    } // End foreach high index
    covariance /= ( ( last_index - ( first_index + lag ) ) - 1 );
    // Now it's the covariance.
  
    // TODO: REMOVE
    //cout << "Autocovariance is " << covariance << endl;
  
    return ( covariance / ( stdev_of_highest_indices * stdev_of_lowest_indices ) );
  } // calculateAutocorrelation(..)

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  class ProfileGibbs {
  public:
    typedef typename profile_traits<ProfileType>::ResidueType ResidueType;
    typedef typename profile_traits<ProfileType>::ProbabilityType ProbabilityType;
    typedef Sequence<SequenceResidueType> SequenceType;

    class Parameters :
      public ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters
    {
      // Boost serialization
    private:
      typedef typename ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters profile_tree_trainer_parameters_t;

      friend class boost::serialization::access;
      template<class Archive>
      void serialize ( Archive & ar, const unsigned int /* file_version */ )
      {
        // save/load base class information
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( profile_tree_trainer_parameters_t );

        ar & BOOST_SERIALIZATION_NVP( sampleProfileGlobals );
        ar & BOOST_SERIALIZATION_NVP( sampleProfilePositions );
        ar & BOOST_SERIALIZATION_NVP( sampleGlobalsFirst );
        ar & BOOST_SERIALIZATION_NVP( numChains );
        ar & BOOST_SERIALIZATION_NVP( minGibbsIterations );
        ar & BOOST_SERIALIZATION_NVP( maxGibbsIterations );
        ar & BOOST_SERIALIZATION_NVP( gibbsIterationsIncrement );
        ar & BOOST_SERIALIZATION_NVP( minRHat );
        ar & BOOST_SERIALIZATION_NVP( burnInFraction );
        ar & BOOST_SERIALIZATION_NVP( useUnconditionalGibbs );
        ar & BOOST_SERIALIZATION_NVP( saveGibbsMode );
        ar & BOOST_SERIALIZATION_NVP( positionShouldBeSampled );
      } // serialize( Archive &, const unsigned int )

    public:
  
      /// PARAMETERS
      /**
       * Should we sample the global parameters at all, or just do positions?
       * If true, sample the globals.  If false, keep them fixed at their
       * starting values.
       */
      // TODO: Implement globals sampling.
      bool sampleProfileGlobals;
  #define DEFAULT_sampleProfileGlobals true
  
      /**
       * Should we sample the position parameters at all?  If true, sample the
       * position-specific parameters.  If false, keep them fixed at their
       * starting values.
       */
      bool sampleProfilePositions;
  #define DEFAULT_sampleProfilePositions true
  
      /**
       * Should we start by sampling the global parameters, then do positions?
       * If true, start with globals.  If false, start with position params,
       * then do globals.
       */
      // TODO: Implement globals sampling.
      bool sampleGlobalsFirst;
  #define DEFAULT_sampleGlobalsFirst false
  
      /**
       * How many chains should we run?  For reliability of the convergence
       * criteria, it is strongly recommended that this be at least 4.
       * @see minRHat
       * Note that this must be at least 1.
       */
      uint32_t numChains;
  #define DEFAULT_numChains 4

      /**
       * How many times at minimum should we iterate ( one iteration includes
       * both positions and globals ), before checking convergence criteria and
       * potentially iteration more?
       * @see gibbsIterationsIncrement
       * @see burnInFraction
       * Note that this must be at least 1.
       */
      uint32_t minGibbsIterations; // at least 1
  #define DEFAULT_minGibbsIterations 30
  
      /**
       * How many times at maximum should we iterate (even if the chains
       * haven't converged)?
       *
       * Note that this must be at least 1.
       */
      uint32_t maxGibbsIterations; // at least 1
  #define DEFAULT_maxGibbsIterations 1000
  
      /**
       * If the chains have not yet converged, how many more iterations should
       * we do before we check again?
       * @see minRHat
       * If this is 0, will use minGibbsIterations instead.
       */
      uint32_t gibbsIterationsIncrement;
  #define DEFAULT_gibbsIterationsIncrement 100

      /**
       * For convergence detection we use the "R Hat" statistic (see the
       * Gelman, et al. book, page 297).
       */
      double minRHat;
      #define DEFAULT_minRHat 1.1

      /**
       * What fraction of the iterations should we discard at the beginning
       * (because they are draws from the not-yet-converged chain)?
       */
      double burnInFraction;
  #define DEFAULT_burnInFraction .2

      /**
       * Use the simpler kind of Gibbs, in which parameter updates are
       * simultaneous?
       * If true, updates all parameters at once.
       * Otherwise, use block-gibbs (position-specific parameters are blocks).
       */
      bool useUnconditionalGibbs;
  #define DEFAULT_useUnconditionalGibbs false

      /**
       * Save the mode?
       * If true, keep track of the highest-scoring profile found while
       * sampling (non-burn-in iterations).  It will be stored in
       * m_bestProfile, and its score in m_bestProfileScore.
       */
      bool saveGibbsMode;
  #define DEFAULT_saveGibbsMode false
    
      /**
       * Minimum profile distribution value.  Applies to everything except the
       * E (End) state.  Given as a double, despite the Profile ProbabilityType.
       */
      double profileValueMinimum;
  #define DEFAULT_profileValueMinimum 1E-5

      /**
       * It is possible to turn off sampling on a position-by-position basis.
       * If this vector is non-null, it must be of the same length as the
       * profile, and the truth value of the i^th position determines whether
       * the i^th profile position should be sampled (true means yes, sample
       * that position; false means keep it fixed at its initial value).
       */
      vector<bool> * positionShouldBeSampled;
  #define DEFAULT_positionShouldBeSampled NULL

      Parameters ();
    
      // Copy constructor
      template <class AnyParameters>
      Parameters ( const AnyParameters & copy_from );
    
      // Copy constructor/operator
      template <class AnyParameters>
      Parameters & operator= (
        AnyParameters const& copy_from
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
      public ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::template ParametersModifierTemplate<ParametersType>
    {
      typedef typename ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::template ParametersModifierTemplate<ParametersType> base_parameters_modifier_t; 

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
        ar & BOOST_SERIALIZATION_NVP( isModified_sampleProfileGlobals );
        ar & BOOST_SERIALIZATION_NVP( isModified_sampleProfilePositions );
        ar & BOOST_SERIALIZATION_NVP( isModified_sampleGlobalsFirst );
        ar & BOOST_SERIALIZATION_NVP( isModified_numChains );
        ar & BOOST_SERIALIZATION_NVP( isModified_minGibbsIterations );
        ar & BOOST_SERIALIZATION_NVP( isModified_maxGibbsIterations );
        ar & BOOST_SERIALIZATION_NVP( isModified_gibbsIterationsIncrement );
        ar & BOOST_SERIALIZATION_NVP( isModified_minRHat );
        ar & BOOST_SERIALIZATION_NVP( isModified_burnInFraction );
        ar & BOOST_SERIALIZATION_NVP( isModified_useUnconditionalGibbs );
        ar & BOOST_SERIALIZATION_NVP( isModified_saveGibbsMode );
        ar & BOOST_SERIALIZATION_NVP( isModified_profileValueMinimum );
        ar & BOOST_SERIALIZATION_NVP( isModified_positionShouldBeSampled );
      } // serialize( Archive &, const unsigned int )

    public:
  
      /// isModified flags for Parameters
      bool isModified_sampleProfileGlobals;
  
      bool isModified_sampleProfilePositions;
  
      bool isModified_sampleGlobalsFirst;
  
      bool isModified_numChains;

      bool isModified_minGibbsIterations;
  
      bool isModified_maxGibbsIterations;
  
      bool isModified_gibbsIterationsIncrement;

      bool isModified_minRHat;

      bool isModified_burnInFraction;

      bool isModified_useUnconditionalGibbs;

      bool isModified_saveGibbsMode;
    
      bool isModified_profileValueMinimum;

      bool isModified_positionShouldBeSampled;

      ParametersModifierTemplate ();
    
      // Copy constructor
      template <class AnyParametersModifierTemplate>
      ParametersModifierTemplate ( const AnyParametersModifierTemplate & copy_from );
    
      // Copy constructor/operator
      template <class AnyParametersModifierTemplate>
      ParametersModifierTemplate & operator= (
        AnyParametersModifierTemplate const& copy_from
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

      void reset ();

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
      );

      template<class AnyParameters>
      void
      applyModifications ( AnyParameters & target_parameters );

      void
      applyModifications ( Parameters & target_parameters );

    }; // End inner class ParametersModifierTemplate

    typedef ParametersModifierTemplate<typename ProfileGibbs::Parameters> ParametersModifier;

    static int const samplingPhaseCount = 2;
    typedef uint8_t SamplingPhase;
    static int const SAMPLING_PHASE_Positions = 0;
    static int const SAMPLING_PHASE_Globals = 1;

    Random * m_random;
    DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType> m_dynamic_programming;
    Parameters m_parameters;
    ProfileType * m_profile;
    vector<SequenceType> const & m_sequences;

    Parameters m_samplingParameters;
    
    /**
     * This is the starting score.  It is not changed once it is set.
     */
    ScoreType m_scoreBeforeSampling;
    
    /** 
     * The score after the performing the current or most recent sampling
     * step.
     */
    ScoreType m_endingScore;
    
    /** 
     * The score at the beginning of the current iteration.
     */
    ScoreType m_startingScore_iteration;

    /** 
     * The current position of the profile (from the end, backwards)
     */
    uint32_t m_row_i;

    /**
     * The index of the current sequence
     */
    uint32_t m_seq_i;

    /**
     * The number of sequences to use in sampling (the index one greater than
     * that of the last one to be used; must be <= m_sequences.size().
     */
    // TODO: What size should this be?
    const uint32_t m_sequence_count;

    /**
     * The current Gibbs iteration.
     */
    uint32_t m_iteration;

    /**
     * The total number of Gibbs iterations, including all chains, burn-in,
     * thinned iters, etc.
     */
    uint32_t m_totalIterations;

    SamplingPhase m_samplingPhase;

    /**
     * The profile being altered by gibbs sampling.  After training, this will
     * be the highest-scoring profile found overall.
     */
    ProfileType m_samplingProfile;

    /**
     * The overall mean profile after calling sample().
     */
    ProfileType m_averageProfile;

    /**
     * The score of the overall mean profile after calling sample().
     */
    ScoreType m_averageProfileScore;

    /**
     * The mode: the best sampled profile (after calling sample()).
     * Used only if parameters.saveGibbsMode is true.
     */
    ProfileType m_bestProfile;

    /**
     * The score at the mode: the score of the best sampled profile (after
     * calling sample()).
     * Used only if parameters.saveGibbsMode is true.
     */
    ScoreType m_bestProfileScore;

    /**
     * The forward rows, indexed first by row, then by sequence.
     */
    vector<typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer> m_forward_matrices;

    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::PositionEntente m_position_entente;

    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::GlobalEntente m_global_entente;

    /**
     * This is the prior we use for the Match-state emission probabilities.
     */
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DirichletMixtureMatchEmissionPrior<float> m_matchEmissionPrior;

    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DirichletMixtureGlobalPrior<float> m_globalPrior;

    /**
     * Construct a profile trainer with the given profile and sequences.
     */  
    ProfileGibbs (
      Random * random,
      ProfileType * profile,
      vector<SequenceType> const & sequences
    );

    /**
     * Construct a profile Gibbs sampler with the given profile and sequences,
     * and the number of sequences to use (use the first num_sequences_to_use
     * sequences only).
     */  
    ProfileGibbs (
      Random * random,
      ProfileType * profile,
      vector<SequenceType> const & sequences,
      uint32_t const & num_sequences_to_use // use only the first X sequences...
    );

    /**
     * Store the given parameters, and starting profile, and set the sampler to
     * its initial values.
     *
     * Note that (for now) the profile must be the same length as the original
     * profile.
     */
    void
    restart (
      const Parameters & parameters,
      ProfileType * profile
    );
  
    ScoreType
    sample ();

  }; // End class ProfileGibbs

  //======//// potentially non-inline implementations ////========//

  ////// Class galosh::ProfileGibbs::Parameters ////
  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_INIT
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
  Parameters ()
  {
    if( DEFAULT_debug >= DEBUG_All ) {
      cout << "[debug] ProfileGibbs::Parameters::<init>()" << endl;
    } // End if DEBUG_All
    this->resetToDefaults();
  } // galosh::ProfileGibbs::Parameters::<init>()

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class AnyParameters>
  GALOSH_INLINE_INIT
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      // Copy constructor
      Parameters ( const AnyParameters & copy_from )
      {
        //if( copy_from.debug >= DEBUG_All ) {
        //  cout << "[debug] ProfileGibbs::Parameters::<init>( copy_from )" << endl;
        //} // End if DEBUG_All
        copyFromNonVirtual( copy_from );
      } // <init>( AnyParameters const & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class AnyParameters>
  GALOSH_INLINE_COPY
  typename ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters &
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      // Copy constructor/operator
      operator= (
        AnyParameters const& copy_from
      )
      {
        if( copy_from.debug >= DEBUG_All ) {
          cout << "[debug] ProfileGibbs::Parameters::operator=( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
        return *this;
      } // operator=( AnyParameters const & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class AnyParameters>
  GALOSH_INLINE_COPY
  void
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      copyFromNonVirtual (
        AnyParameters const & copy_from
      )
      {
        ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::copyFromNonVirtual( copy_from );
        //if( copy_from.debug >= DEBUG_All ) {
        //  cout << "[debug] ProfileGibbs::Parameters::copyFromNonVirtual( copy_from )" << endl;
        //} // End if DEBUG_All
        copyFromNonVirtualDontDelegate( copy_from );
      } // copyFromNonVirtual( AnyParameters const & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class AnyParameters>
  GALOSH_INLINE_COPY
  void
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
  copyFromNonVirtualDontDelegate (
    AnyParameters const & copy_from
  )
  {
    sampleProfileGlobals = copy_from.sampleProfileGlobals;
    sampleProfilePositions = copy_from.sampleProfilePositions;
    sampleGlobalsFirst = copy_from.sampleGlobalsFirst;
    numChains = copy_from.numChains;
    minGibbsIterations = copy_from.minGibbsIterations;
    maxGibbsIterations = copy_from.maxGibbsIterations;
    gibbsIterationsIncrement = copy_from.gibbsIterationsIncrement;
    minRHat = copy_from.minRHat;
    burnInFraction = copy_from.burnInFraction;
    useUnconditionalGibbs = copy_from.useUnconditionalGibbs;
    saveGibbsMode = copy_from.saveGibbsMode;
    profileValueMinimum = copy_from.profileValueMinimum;
    positionShouldBeSampled = copy_from.positionShouldBeSampled;
  } // copyFromNonVirtualDontDelegate( AnyParameters const & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_COPY
  void
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      copyFrom ( const Parameters & copy_from )
      {
        copyFromNonVirtual( copy_from );
      } // copyFrom( Parameters const & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_REINITIALIZE
  void
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      resetToDefaults ()
      {
        ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::resetToDefaults();
        // TODO: Why isn't the compiler finding "debug" in galosh::Parameters?
        //if( debug >= DEBUG_All ) {
        //  cout << "[debug] ProfileGibbs::Parameters::resetToDefaults()" << endl;
        //} // End if DEBUG_All
        sampleProfileGlobals = DEFAULT_sampleProfileGlobals;
        sampleProfilePositions = DEFAULT_sampleProfilePositions;
        sampleGlobalsFirst = DEFAULT_sampleGlobalsFirst;
        numChains = DEFAULT_numChains;
        minGibbsIterations = DEFAULT_minGibbsIterations;
        maxGibbsIterations = DEFAULT_maxGibbsIterations;
        gibbsIterationsIncrement = DEFAULT_gibbsIterationsIncrement;
        minRHat = DEFAULT_minRHat;
        burnInFraction = DEFAULT_burnInFraction;
        useUnconditionalGibbs = DEFAULT_useUnconditionalGibbs;
        saveGibbsMode = DEFAULT_saveGibbsMode;
        profileValueMinimum = DEFAULT_profileValueMinimum;
        positionShouldBeSampled = DEFAULT_positionShouldBeSampled;
      } // resetToDefaults()
    
  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
      template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
      void
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      writeParameters ( 
        std::basic_ostream<CharT,Traits>& os
      ) const
      {
        ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::writeParameters( os );
        os << endl;

        os << "[ProfileGibbs]" << endl;
        os << "sampleProfileGlobals = " << sampleProfileGlobals << endl;
        os << "sampleProfilePositions = " << sampleProfilePositions << endl;
        os << "sampleGlobalsFirst = " << sampleGlobalsFirst << endl;
        os << "numChains = " << numChains << endl;
        os << "minGibbsIterations = " << minGibbsIterations << endl;
        os << "maxGibbsIterations = " << maxGibbsIterations << endl;
        os << "gibbsIterationsIncrement = " << gibbsIterationsIncrement << endl;
        os << "minRHat = " << minRHat << endl;
        os << "burnInFraction = " << burnInFraction << endl;
        os << "useUnconditionalGibbs = " << useUnconditionalGibbs << endl;
        os << "saveGibbsMode = " << saveGibbsMode << endl;
        os << "profileValueMinimum = " << profileValueMinimum << endl;
        os << "positionShouldBeSampled = " << positionShouldBeSampled << endl;

        return os;
      } // writeParameters( basic_ostream & ) const

  ////// Class galosh::ProfileGibbs::ParametersModifierTemplate ////
  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class ParametersType>
  GALOSH_INLINE_INIT
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      ParametersModifierTemplate ()
      {
        if( base_parameters_modifier_t::parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProfileGibbs::ParametersModifierTemplate::<init>()" << endl;
        } // End if DEBUG_All
        isModified_reset();
      } // <init>()

      // Copy constructor
  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_INIT
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      ParametersModifierTemplate ( const AnyParametersModifierTemplate & copy_from )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProfileGibbs::ParametersModifierTemplate::<init>( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
      } // <init>( AnyParametersModifierTemplate const & )
    
  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_COPY
      // Copy constructor/operator
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType> &
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      operator= (
        AnyParametersModifierTemplate const& copy_from
      )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProfileGibbs::ParametersModifierTemplate::operator=( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
        return *this;
      } // operator=( AnyParametersModifierTemplate const & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_COPY
      void
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      isModified_copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      )
      {
        base_parameters_modifier_t::isModified_copyFromNonVirtual( copy_from );
        isModified_sampleProfileGlobals = copy_from.isModified_sampleProfileGlobals;
        isModified_sampleProfilePositions = copy_from.isModified_sampleProfilePositions;
        isModified_sampleGlobalsFirst = copy_from.isModified_sampleGlobalsFirst;
        isModified_numChains = copy_from.isModified_numChains;
        isModified_minGibbsIterations = copy_from.isModified_minGibbsIterations;
        isModified_maxGibbsIterations = copy_from.isModified_maxGibbsIterations;
        isModified_gibbsIterationsIncrement = copy_from.isModified_gibbsIterationsIncrement;
        isModified_minRHat = copy_from.isModified_minRHat;
        isModified_burnInFraction = copy_from.isModified_burnInFraction;
        isModified_useUnconditionalGibbs = copy_from.isModified_useUnconditionalGibbs;
        isModified_saveGibbsMode = copy_from.isModified_saveGibbsMode;
        isModified_profileValueMinimum = copy_from.isModified_profileValueMinimum;
        isModified_positionShouldBeSampled = copy_from.isModified_positionShouldBeSampled;
      } // isModified_copyFromNonVirtual ( AnyParametersModifierTemplate const & )


  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class ParametersType>
  GALOSH_INLINE_REINITIALIZE
  void
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
  reset ()
      {
        isModified_reset();
        base_parameters_modifier_t::parameters.resetToDefaults();
      } // reset()

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class ParametersType>
  GALOSH_INLINE_REINITIALIZE
      void
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      isModified_reset ()
      {
        base_parameters_modifier_t::isModified_reset();

        isModified_sampleProfileGlobals = false;
        isModified_sampleProfilePositions = false;
        isModified_sampleGlobalsFirst = false;
        isModified_numChains = false;
        isModified_minGibbsIterations = false;
        isModified_maxGibbsIterations = false;
        isModified_gibbsIterationsIncrement = false;
        isModified_minRHat = false;
        isModified_burnInFraction = false;
        isModified_useUnconditionalGibbs = false;
        isModified_saveGibbsMode = false;
        isModified_profileValueMinimum = false;
        isModified_positionShouldBeSampled = false;
      } // isModified_reset()

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class ParametersType>
  template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
      void
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      writeParametersModifier (
        std::basic_ostream<CharT,Traits>& os
      )
      {
        DynamicProgramming<ResidueType, ProbabilityType,ScoreType,MatrixValueType>::template ParametersModifierTemplate<ParametersType>::writeParametersModifier( os );
        os << endl;

        os << "[ProfileGibbs]" << endl;
        if( isModified_sampleProfileGlobals ) {
          os << "sampleProfileGlobals = " << base_parameters_modifier_t::parameters.sampleProfileGlobals << endl;
        }
        if( isModified_sampleProfilePositions ) {
          os << "sampleProfilePositions = " << base_parameters_modifier_t::parameters.sampleProfilePositions << endl;
        }
        if( isModified_sampleGlobalsFirst ) {
          os << "sampleGlobalsFirst = " << base_parameters_modifier_t::parameters.sampleGlobalsFirst << endl;
        }
        if( isModified_numChains ) {
          os << "numChains = " << base_parameters_modifier_t::parameters.numChains << endl;
        }
        if( isModified_minGibbsIterations ) {
          os << "minGibbsIterations = " << base_parameters_modifier_t::parameters.minGibbsIterations << endl;
        }
        if( isModified_maxGibbsIterations ) {
          os << "maxGibbsIterations = " << base_parameters_modifier_t::parameters.maxGibbsIterations << endl;
        }
        if( isModified_gibbsIterationsIncrement ) {
          os << "gibbsIterationsIncrement = " << base_parameters_modifier_t::parameters.gibbsIterationsIncrement << endl;
        }
        if( isModified_minRHat ) {
          os << "minRHat = " << base_parameters_modifier_t::parameters.minRHat << endl;
        }
        if( isModified_burnInFraction ) {
          os << "burnInFraction = " << base_parameters_modifier_t::parameters.burnInFraction << endl;
        }
        if( isModified_useUnconditionalGibbs ) {
          os << "useUnconditionalGibbs = " << base_parameters_modifier_t::parameters.useUnconditionalGibbs << endl;
        }
        if( isModified_saveGibbsMode ) {
          os << "saveGibbsMode = " << base_parameters_modifier_t::parameters.saveGibbsMode << endl;
        }
        if( isModified_profileValueMinimum ) {
          os << "profileValueMinimum = " << base_parameters_modifier_t::parameters.profileValueMinimum << endl;
        }
        if( isModified_positionShouldBeSampled ) {
          os << "positionShouldBeSampled = " << base_parameters_modifier_t::parameters.positionShouldBeSampled << endl;
        }
      } // writeParametersModifier ( basic_ostream & ) const

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class ParametersType>
  template <class AnyParametersExceptGibbsParameters>
  GALOSH_INLINE_PARAMETERSMODIFIER_APPLY_MODIFICATIONS
      void
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      applyModifications ( AnyParametersExceptGibbsParameters & target_parameters )
      {
        base_parameters_modifier_t::applyModifications( target_parameters );
      } // applyModifications( AnyParametersExceptGibbsParameters & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class ParametersType>
  GALOSH_INLINE_PARAMETERSMODIFIER_APPLY_MODIFICATIONS
      void
  // Note this is explicitly defined to only override when the type of Parameters being applied to is a ProfileGibbs::Parameters type.
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      applyModifications ( typename ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters & target_parameters )
      {
        base_parameters_modifier_t::applyModifications( target_parameters );

        if( isModified_sampleProfileGlobals ) {
          target_parameters.sampleProfileGlobals =
            base_parameters_modifier_t::parameters.sampleProfileGlobals;
        }
        if( isModified_sampleProfilePositions ) {
          target_parameters.sampleProfilePositions =
            base_parameters_modifier_t::parameters.sampleProfilePositions;
        }
        if( isModified_sampleGlobalsFirst ) {
          target_parameters.sampleGlobalsFirst =
            base_parameters_modifier_t::parameters.sampleGlobalsFirst;
        }
        if( isModified_numChains ) {
          target_parameters.numChains =
            base_parameters_modifier_t::parameters.numChains;
        }
        if( isModified_minGibbsIterations ) {
          target_parameters.minGibbsIterations =
            base_parameters_modifier_t::parameters.minGibbsIterations;
        }
        if( isModified_maxGibbsIterations ) {
          target_parameters.maxGibbsIterations =
            base_parameters_modifier_t::parameters.maxGibbsIterations;
        }
        if( isModified_gibbsIterationsIncrement ) {
          target_parameters.gibbsIterationsIncrement =
            base_parameters_modifier_t::parameters.gibbsIterationsIncrement;
        }
        if( isModified_minRHat ) {
          target_parameters.minRHat =
            base_parameters_modifier_t::parameters.minRHat;
        }
        if( isModified_burnInFraction ) {
          target_parameters.burnInFraction =
            base_parameters_modifier_t::parameters.burnInFraction;
        }
        if( isModified_useUnconditionalGibbs ) {
          target_parameters.useUnconditionalGibbs =
            base_parameters_modifier_t::parameters.useUnconditionalGibbs;
        }
        if( isModified_saveGibbsMode ) {
          target_parameters.saveGibbsMode =
            base_parameters_modifier_t::parameters.saveGibbsMode;
        }
        if( isModified_profileValueMinimum ) {
          target_parameters.profileValueMinimum =
            base_parameters_modifier_t::parameters.profileValueMinimum;
        }
        if( isModified_positionShouldBeSampled ) {
          target_parameters.positionShouldBeSampled =
            base_parameters_modifier_t::parameters.positionShouldBeSampled;
        }

      } // applyModifications( AnyParameters & )

  ////// Class galosh::ProfileGibbs ////
  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_INIT
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::
    /**
     * Construct a profile Gibbs sampler with the given profile and sequences.
     */  
    ProfileGibbs (
      Random * random,
      ProfileType * profile,
      vector<SequenceType> const & sequences
    ) :
      m_random( random ),
      m_profile( profile ),
      m_sequences( sequences ),
      m_sequence_count( sequences.size() ),
      m_forward_matrices(), // see restart()
      m_dynamic_programming(),
      m_global_entente(), // see restart()
      m_position_entente(), // see restart()
      m_matchEmissionPrior(), // see restart()
      m_globalPrior(), // see restart()
      m_averageProfile(),
      m_bestProfile()
    {
      if( m_parameters.debug >= DEBUG_All ) {
        cout << "[debug] ProfileGibbs::<init>( ProfileType, vector<SequenceType> )" << endl;
      } // End if DEBUG_All
      // Do nothing else
    } // <init>( ProfileType *, vector<SequenceType> const & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_INIT
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::
    /**
     * Construct a profile Gibbs sampler with the given profile and sequences,
     * and the number of sequences to use (use the first num_sequences_to_use
     * sequences only).
     */  
    ProfileGibbs (
      Random * random,
      ProfileType * profile,
      vector<SequenceType> const & sequences,
      uint32_t const & num_sequences_to_use // use only the first X sequences...
    ) :
      m_random( random ),
      m_profile( profile ),
      m_sequences( sequences ),
      m_sequence_count( min( ( size_t )( ( num_sequences_to_use > 0 ) ? num_sequences_to_use : sequences.size() ), sequences.size() ) ),
      m_forward_matrices(), // see restart()
      m_dynamic_programming(),
      m_global_entente(), // see restart()
      m_position_entente(), // see restart()
      m_matchEmissionPrior(), // see restart()
      m_globalPrior(), // see restart()
      m_averageProfile(),
      m_bestProfile()
    {
      if( m_parameters.debug >= DEBUG_All ) {
        cout << "[debug] ProfileGibbs::<init>( ProfileType, vector<SequenceType> )" << endl;
      } // End if DEBUG_All
      // Do nothing else
    } // <init>( ProfileType *, vector<SequenceType> const &, uint32_t const & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_REINITIALIZE
  void
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::
    /**
     * Store the given parameters, and starting profile, and set the trainer to
     * its initial values.
     *
     * Note that (for now) the profile must be the same length as the original
     * profile.
     */
    restart (
      const Parameters & parameters,
      ProfileType * profile
    )
    {
      if( parameters.debug >= DEBUG_All ) {
        cout << "[debug] ProfileGibbs::restart( parameters, profile )" << endl;
      } // End if DEBUG_All
      m_samplingParameters = parameters;

      m_profile = profile;
      m_samplingProfile = *m_profile;
      if( m_samplingParameters.saveGibbsMode ) {
        m_bestProfile = *m_profile;
        m_bestProfileScore = 0;
      }

      // Start out sampling position-specific params.  Globals after.  Or the
      // other way 'round, depending on the parameters.
      // TODO: USE THIS
      m_samplingPhase =
        ( m_samplingParameters.sampleGlobalsFirst ?
          SAMPLING_PHASE_Globals :
          SAMPLING_PHASE_Positions );

      m_position_entente.reinitialize();
      m_global_entente.reinitialize();

      if( m_samplingParameters.usePriors ) {
        if( m_samplingParameters.matchEmissionPrior != NULL ) {
          m_matchEmissionPrior =
            *( m_samplingParameters.matchEmissionPrior );
        } else {
          m_matchEmissionPrior.reinitializeToLaplace();
        }

        // TODO: REMOVE
        //cout << "MatchEmission prior is " << m_matchEmissionPrior << endl;

        if( m_samplingParameters.globalPrior != NULL ) {
          m_globalPrior =
            *( m_samplingParameters.globalPrior );
        } else {
          m_globalPrior.reinitializeToLaplace();
        }

        // TODO: REMOVE
        //cout << "Global prior is " << m_globalPrior << endl;

      } else {
        m_matchEmissionPrior.reinitializeToEven( 0.0f );
        m_globalPrior.reinitializeToEven( 0.0f );
      } // End if usePriors .. else ..
    } // restart( const Parameters &, const ProfileType & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_GIBBS
    ScoreType
  ProfileGibbs<ProfileType, ScoreType, MatrixValueType, SequenceResidueType>::
    sample ()
    {
      // TODO: REMOVE
      bool be_extra_verbose = false;
      bool dont_zero_ententes_hack = false;

      if( be_extra_verbose || ( m_parameters.debug >= DEBUG_All ) ) {
        cout << "[debug] ProfileGibbs::sample()" << endl;
      } // End if DEBUG_All
      // TODO: REMOVE
      //cout << "BEFORE RESTARTING, PROFILE IS " << *m_profile << endl;
      restart( m_parameters, m_profile );

      // MARK
      ScoreType score;

      // TODO: Make some of these member vars so we don't keep reallocating the memory.
      // TODO: Don't allocate the whole counts matrix if the type is PerPosition.

      // TODO: Make an Algebra method for uin32_t type?
      //vector<ProfileTreeRoot<ResidueType, uint32_t> > counts;
      vector<ProfileTreeRoot<ResidueType, doublerealspace> > counts;

        // TODO: REMOVE.  Only used for the dont_zero_ententes_hack
        vector<typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::PositionEntente > position_ententes;
        vector<typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::GlobalEntente > global_ententes;

        typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer::reverse_iterator forward_matrices_for_chain_i_iter;

        // For Per Position Gibbs:
        // Every other position, we must swap which of the SamplerState vectors
        // is current.
        vector<typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::SamplerState> sampler_states_1( m_sequence_count );
        vector<typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::SamplerState> sampler_states_2( m_sequence_count );
        bool swap = true;

        vector<ProfileType> starting_profiles;
        vector<ProfileType> profiles;
        vector<ProfileType> average_profiles;
        vector<double> average_log_scores;
        vector<double> log_score_variances;
        vector<uint32_t> log_score_autocorrelation_first_lag_below_thresholds;

        vector<ScoreType> best_scores;
        //vector<ProfileType> best_profiles;
        vector<vector<double> > log_scores;

        double average_log_score;
        double log_score_variance;
        double average_within_chain_log_score_variance;
        double between_chain_log_score_variance;
        double var_hat_plus;
        double V_hat;
        double R_hat;
        double log_score_stdev;

        ScoreType best_score_overall = INIT_PROBABILITY(ScoreType)( 0.0 );
        // Note that, to speed things up, I commented out the storing of the best profile found during sampling, since the average profile is almost always the best.
        ProfileType best_profile_overall = m_samplingProfile;
        double log_score;

        bool done;
        uint32_t chain_i, iter_i, iter_j;
        uint32_t last_row = m_samplingProfile.length();
        uint32_t row_i;

        double autocorrelation;
        uint32_t lag;

        // For Per Position Gibbs:
        uint32_t seq_i;

        uint32_t gibbs_type; // There's two kinds : 0 = Simple, 1 = Per Pos
        vector<string> gibbs_type_label = vector<string>( 2 );
        gibbs_type_label[ 0 ] = "Simple"; // unconditionalGibbs
        gibbs_type_label[ 1 ] = "Per Position";



      uint32_t max_iters;
      uint32_t iterations_increment =
        (
          ( m_samplingParameters.gibbsIterationsIncrement == 0 ) ?
          m_samplingParameters.minGibbsIterations :
          m_samplingParameters.gibbsIterationsIncrement
        );
      uint32_t burn_in_iters;
      uint32_t max_lag;
      double autocorrelation_threshold;

        ///// PARAMETERS ////
        //m_samplingParameters.numChains = 10;
        max_iters = m_samplingParameters.minGibbsIterations;
        burn_in_iters =
          max(
            ( uint32_t )1,
            ( uint32_t )(
              m_samplingParameters.burnInFraction *
              max_iters
            )
          );

        //m_samplingParameters.minRHat = 1.1;
        //m_samplingParameters.gibbsIterationsIncrement = ( max_iters - burn_in_iters );//100;

        // 0 means don't calculate autocorrelations at all.
        max_lag = 0;//100;
        autocorrelation_threshold = .1;

        // TODO: Let the user supply the starting profiles
        // Set up the starting profiles
        starting_profiles.resize( m_samplingParameters.numChains );
        for( chain_i = 0; chain_i < m_samplingParameters.numChains; chain_i++ ) {
          // Argh
          //starting_profiles[ chain_i ] = m_samplingProfile;
          starting_profiles[ chain_i ].copyFrom( m_samplingProfile );

          if( m_samplingParameters.sampleProfilePositions ) {
            if( m_samplingParameters.usePriors ) {
              // Draw the params from the prior.
              starting_profiles[ chain_i ].dirichletMixturePositions(
                m_matchEmissionPrior,
                *m_random
              );
            } else { // if usePriors .. else .. 
              starting_profiles[ chain_i ].uniformPositions( *m_random );
            } // End if usePriors .. else .. 
          } // End if sampleProfilePositions
          if( m_samplingParameters.sampleProfileGlobals ) {
            if( m_samplingParameters.usePriors ) {
              starting_profiles[ chain_i ].dirichletMixtureExceptPositions(
                m_globalPrior,
                *m_random
              );
            } else { // if usePriors .. else ..
              starting_profiles[ chain_i ].uniformExceptPositions( *m_random );
            }  // End if usePriors .. else ..
          } // End if sampleProfileGlobals

          // TODO: REMOVE
         //// open the archive
         // string filename;
         // if( chain_i == 0 ) {
         //   filename = "test_results/v7_seed1205808921_typeD/ProfuseTest.v7_seed1205808921_typeD.40.100.0.100.2.0.starting.ProfileTree.xml";
         // } else if( chain_i == 1 ) {
         //   filename = "test_results/v7_seed1205808921_typeD/ProfuseTest.v7_seed1205808921_typeD.40.100.0.100.2.1.starting.ProfileTree.xml";
         // } else if( chain_i == 2 ) {
         //   filename = "test_results/v7_seed1205808921_typeD/ProfuseTest.v7_seed1205808921_typeD.40.100.0.100.2.2.starting.ProfileTree.xml";
         // } else if( chain_i == 3 ) {
         //   filename = "test_results/v7_seed1205808921_typeD/ProfuseTest.v7_seed1205808921_typeD.40.100.0.100.2.3.starting.ProfileTree.xml";
         // }
         // std::ifstream ifs( filename.c_str() );
         //assert( ifs.good() );
         //boost::archive::xml_iarchive ia( ifs );
         //// restore the profile from the archive
         //ia >> BOOST_SERIALIZATION_NVP( starting_profiles[ chain_i ] );
         //starting_profiles[ chain_i ].copyExceptPositions( m_samplingProfile );

          if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
            cout << "[Chain " << chain_i << "] Starting profile: " << starting_profiles[ chain_i ] << endl;
          }

          // TODO: REMOVE? Train each starting profile separately.
          // m_samplingProfile.copyFrom( starting_profiles[ chain_i ] );
          // m_samplingProfile.normalize( 1E-5 ); // Ensure values aren't too tiny.
          // m_samplingParameters.trainProfileGlobals = false;
          // m_samplingParameters.verbosity = VERBOSITY_Low;
          // //cout << "The profile (before) is:" << endl;
          // //cout << m_samplingProfile << endl;
          // score = trainer2.train();
          // cout << "Now (after training), the score is " << score << ", and the profile is:" << endl;
          // cout << *trainer2.m_profile << endl;

        } // End foreach chain, draw the starting profile.
        
        //for( gibbs_type = 0; gibbs_type < 2; gibbs_type++ ) {
        gibbs_type = ( m_samplingParameters.useUnconditionalGibbs ? 0 : 1 );
        if( !m_samplingParameters.sampleProfilePositions ) {
          // If we're not sampling positions, then the types are the same.
          gibbs_type = 0;
        }
        if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
          if( m_samplingParameters.sampleProfilePositions ) {
            cout << "Using (" << gibbs_type_label[ gibbs_type ] << ") Gibbs:" << endl;
          } else {
            cout << "Using (Globals-only) Gibbs:" << endl;
          }
        }

        // Setup the first time through, or if m_samplingParameters.numChains depends on
        // gibbs_type.
          if( profiles.size() != m_samplingParameters.numChains ) {
            counts.resize( m_samplingParameters.numChains );
            if( dont_zero_ententes_hack ) {
              position_ententes.resize( m_samplingParameters.numChains );
              global_ententes.resize( m_samplingParameters.numChains );
            }
            profiles.resize( m_samplingParameters.numChains );
            average_profiles.resize( m_samplingParameters.numChains );
            //best_profiles.resize( m_samplingParameters.numChains );
            average_log_scores.resize( m_samplingParameters.numChains );
            log_score_variances.resize( m_samplingParameters.numChains );
            log_score_autocorrelation_first_lag_below_thresholds.resize( m_samplingParameters.numChains );
            best_scores.resize( m_samplingParameters.numChains );
            log_scores.resize( m_samplingParameters.numChains );
            m_forward_matrices.resize( m_samplingParameters.numChains );
          } // ( profiles.size() != m_samplingParameters.numChains )

          // Chain-specific setup
          for( chain_i = 0; chain_i < m_samplingParameters.numChains; chain_i++ ) {
            // TODO: REMOVE
            //cout << "About to call counts[ chain_i ].reinitialize(..)" << endl;
            counts[ chain_i ].reinitialize(
              m_samplingProfile.length()
            );
            if( dont_zero_ententes_hack ) {
              position_ententes[ chain_i ].reinitialize();
              global_ententes[ chain_i ].reinitialize();
            } // End if dont_zero_ententes_hack
            // Argh.
            //profiles[ chain_i ] = starting_profiles[ chain_i ];
            profiles[ chain_i ].copyFrom( starting_profiles[ chain_i ] );
            average_profiles[ chain_i ].reinitialize(
              profiles[ chain_i ].length()
            );
            average_profiles[ chain_i ].zero();
            m_averageProfile.reinitialize(
              m_samplingProfile.length()
            );
 
            m_forward_matrices[ chain_i ].reinitialize(
              profiles[ chain_i ].length(),
              m_sequences,
              m_sequence_count
            );
            // And fill the forward matrices..
            score =
              m_dynamic_programming.forward_score( 
                m_samplingParameters,
                profiles[ chain_i ],
                m_sequences,
                m_sequence_count,
                m_forward_matrices[ chain_i ]
              );
            if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
              cout << "\tScore: " << score << endl;
            }
            m_dynamic_programming.drawPaths(
              m_samplingParameters,
              profiles[ chain_i ],
              m_sequences,
              m_sequence_count,
              m_forward_matrices[ chain_i ],
              *m_random,
              &( counts[ chain_i ] ),
              ( typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template MultipleAlignment<ProfileType, SequenceResidueType> * )NULL
            );
          
            best_scores[ chain_i ] = score;
            if(
              m_samplingParameters.saveGibbsMode &&
              ( m_bestProfileScore < score )
            ) {
              m_bestProfileScore = score;
              m_bestProfile.copyFrom( profiles[ chain_i ] );
            }
            average_log_scores[ chain_i ] = 0.0;

            // TODO: REMOVE
            //cout << " ABOUT TO SET UP best_profiles[ " << chain_i << " ]" << endl;
            // Argh.
            //best_profiles[ chain_i ] = profiles[ chain_i ];
            //best_profiles[ chain_i ].copyFrom( profiles[ chain_i ] );
          } // End foreach chain_i, do setup.

          // Iterate until the chains are sufficiently converged (We need to
          // gather a large enough sample to ensure that all R_hat measurements
          // are less than 1.1).
          done = false;
          iter_i = 0;
          m_totalIterations = 0;
          do { // while( !done )
            if( true || be_extra_verbose ) {
              cout << "max_iters is now " << max_iters << endl;
              cout << "burn_in_iters is now " << burn_in_iters << endl;
            }

            // Setup for a new max_iters
            for( chain_i = 0; chain_i < m_samplingParameters.numChains; chain_i++ ) {
              log_scores[ chain_i ].resize( max_iters - burn_in_iters );
            } // End foreach chain_i, do more setup.

            for( ; iter_i < max_iters; iter_i++, m_totalIterations++ ) {
              if( false && be_extra_verbose ) {
                cout << "iteration " << iter_i << ":" << endl;
              }
              for( chain_i = 0; chain_i < m_samplingParameters.numChains; chain_i++ ) {
                if( false && be_extra_verbose ) {
                  cout << "chain " << chain_i << ":" << endl;
                }
                if( gibbs_type == 1 ) {
                  counts[ chain_i ].zero();
                } // End if using Per Position Gibbs
                // TODO: REMOVE.
                //else {
                //  cout << "counts[ " << chain_i << " ] is " << counts[ chain_i ] << endl;
                //}
                // When using Simple Gibbs, start at row_i = profile_length - 1;
                row_i = ( ( gibbs_type == 0 ) ? ( last_row - 1 ) : last_row );
                forward_matrices_for_chain_i_iter = m_forward_matrices[ chain_i ].rbegin();
                do { // Foreach row_i downto (and including) 0
                  if( false && be_extra_verbose ) {
                    cout << "row " << row_i << ":" << endl;
                  }
                  // First calculate the position_entente.
                  if( row_i != last_row  ) {
                    if( gibbs_type == 1 ) {
                      // Every other time, we must swap which of the
                      // SamplerStates is current.
                      swap = !swap;
                    
                      for( seq_i = 0; seq_i < m_sequence_count; seq_i++ ) {
                        m_dynamic_programming.drawPartialPathForPosition(
                          m_samplingParameters,
                          profiles[ chain_i ],
                          m_sequences[ seq_i ],
                          row_i,
                          ( *forward_matrices_for_chain_i_iter )[ seq_i ],
                          ( swap ? sampler_states_2[ seq_i ] : sampler_states_1[ seq_i ] ), // ignored if row_i == last_row
                          ( swap ? sampler_states_1[ seq_i ] : sampler_states_2[ seq_i ] ),
                          *m_random,
                          // TODO: Make an Algebra method for uint32_t type?
                          //( ProfileTreeRoot<ResidueType, uint32_t> * )NULL,
                          ( ProfileTreeRoot<ResidueType, doublerealspace> * )NULL,
                          &( counts[ chain_i ][ row_i ] ),
                          ( vector<uint32_t> * )NULL
                        );
                      } // End foreach seq_i
                    } // End if using Per Position Gibbs .. else ..
                    if( m_samplingParameters.sampleProfilePositions ) {
                      if( dont_zero_ententes_hack ) { //&& ( iter_i >= burn_in_iters ) ) {
                        position_ententes[ chain_i ] += counts[ chain_i ][ row_i ];
                        m_position_entente = position_ententes[ chain_i ];
                      } else {
                        // Argh!
                        //m_position_entente = counts[ chain_i ][ row_i ];
                        m_position_entente.copyFrom( counts[ chain_i ][ row_i ] );
                      }
                      if( m_samplingParameters.usePriors ) {
                        // TODO: REMOVE
                        //cout << "The position entente, before incorporating any priors, is " << m_position_entente << endl;
                        // Do it.
                        m_matchEmissionPrior.incorporatePrior( m_position_entente );
                        // TODO: REMOVE
                        //cout << "The position entente, after incorporating a simple Laplace prior, is " << m_position_entente << endl;
                      } // End if usePriors
                      
                      m_position_entente.normalize( m_samplingParameters.profileValueMinimum );
                      // Try drawing from a Dirichlet with those counts for the positions:
                      profiles[ chain_i ][ row_i ].dirichlet( m_position_entente, *m_random );
                      // TODO: REMOVE
                      //cout << "After drawing from a dirichlet( " << m_position_entente << " ), profile[ " << row_i << " ] is " << profiles[ chain_i ][ row_i ] << endl;
                    } // End if sampleProfilePositions
                  } // End if row_i != last_row

                  if( gibbs_type == 1 ) {
                    // Now draw the next bit of partial path...
                    if( row_i != last_row ) {
                      counts[ chain_i ][ row_i ].zero();
                    }
                    // TODO: REMOVE
                    //else {
                    //  cout << "counts[ " << chain_i << " ][ " << row_i << " ] are " << counts[ chain_i ][ row_i ] << endl;
                    //}
                    for( seq_i = 0; seq_i < m_sequence_count; seq_i++ ) {
                      if( false && be_extra_verbose ) {
                        cout << "Drawing partial path for seq " << seq_i << ".." << endl;
                      }

                      // Update counts for this position and sequence
                      m_dynamic_programming.drawPartialPathForPosition(
                        m_samplingParameters,
                        profiles[ chain_i ],
                        m_sequences[ seq_i ],
                        row_i,
                        ( *forward_matrices_for_chain_i_iter )[ seq_i ],
                        ( swap ? sampler_states_2[ seq_i ] : sampler_states_1[ seq_i ] ), // ignored if row_i == last_row
                        ( swap ? sampler_states_1[ seq_i ] : sampler_states_2[ seq_i ] ),
                        *m_random,
                        ( ( row_i == last_row ) ? NULL : &( counts[ chain_i ] ) ),
                        ( ( row_i == last_row ) ? NULL : &( counts[ chain_i ][ row_i ] ) ),
                        ( vector<uint32_t> * )NULL
                      );
                    } // End foreach seq_i
                  } // End if using Per Position Gibbs, update counts

                  ++forward_matrices_for_chain_i_iter;
                } while( row_i-- > 0 ); // End for each forward row, downto 0 (incl. 0)
            
                // Globals. TODO: use sampleGlobalsFirst
                if( m_samplingParameters.sampleProfileGlobals ) {
                  if( dont_zero_ententes_hack ) { //&& ( iter_i >= burn_in_iters ) ) {
                    global_ententes[ chain_i ] += counts[ chain_i ];
                    m_global_entente = global_ententes[ chain_i ];
                  } else {
                    m_global_entente = counts[ chain_i ];
                  }
                  if( m_samplingParameters.usePriors ) {
                    // TODO: REMOVE
                    //cout << "The global entente, before incorporating any priors, is " << m_global_entente << endl;
                    m_globalPrior.incorporatePrior( m_global_entente );
                    // TODO: REMOVE
                    //cout << "The global entente, after incorporating the prior, is " << m_global_entente << endl;
                  } // End if usePriors
                  m_global_entente.normalize( m_samplingParameters.profileValueMinimum );
                  // Try drawing from a Dirichlet with those counts for Match transitions:
                  profiles[ chain_i ].dirichletExceptPositions( m_global_entente, *m_random );
                  //cout << "After drawing from a dirichlet( " << m_global_entente << " ), profile globals are ";
                  //profile.writeExceptPositions( cout );
                  //cout << endl;
                } // End if m_samplingParameters.sampleProfileGlobals

                score =
                  m_dynamic_programming.forward_score( 
                    m_samplingParameters,
                    profiles[ chain_i ],
                    m_sequences,
                    m_sequence_count,
                    m_forward_matrices[ chain_i ]
                  );
                if( false || be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
                  cout << "[" << gibbs_type_label[ gibbs_type ] << ", Iteration " << iter_i << ", Chain " << chain_i << "] Score: " << score << endl;
                }
                log_score = toLogDouble( score );
            
                // Profile has changed since the last iteration.
                if( iter_i >= burn_in_iters ) {
                  average_profiles[ chain_i ] += profiles[ chain_i ];
                  average_log_scores[ chain_i ] += log_score;
                  log_scores[ chain_i ][ iter_i - burn_in_iters ] = log_score;
                  if( score > best_scores[ chain_i ] ) {
                    best_scores[ chain_i ] = score;
                    // Argh.
                    //best_profiles[ chain_i ] = profiles[ chain_i ];
                    //best_profiles[ chain_i ].copyFrom( profiles[ chain_i ] );
                  }
                  if(
                    m_samplingParameters.saveGibbsMode &&
                    ( m_bestProfileScore < score )
                  ) {
                    m_bestProfileScore = score;
                    m_bestProfile.copyFrom( profiles[ chain_i ] );
                  }
                } // End if we're beyond the burn-in.
            
                if( gibbs_type == 0 ) { // if using Simple Gibbs
                  m_dynamic_programming.drawPaths(
                    m_samplingParameters,
                    profiles[ chain_i ],
                    m_sequences,
                    m_sequence_count,
                    m_forward_matrices[ chain_i ],
                    *m_random,
                    &counts[ chain_i ],
                    ( typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template MultipleAlignment<ProfileType, SequenceResidueType> * )NULL
                  );
                } // End if using Simple Gibbs
              } // End foreach chain_i
            } // End foreach iter_i

            // Report per-chain best profiles
            for( chain_i = 0; chain_i < m_samplingParameters.numChains; chain_i++ ) {
              if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
                //cout << "[" << gibbs_type_label[ gibbs_type ] << ", Chain " << chain_i << "] The best profile is " << best_profiles[ chain_i ];
                //cout << "\tIts score is " << best_scores[ chain_i ] << endl;
                cout << "\tBest score is " << best_scores[ chain_i ] << endl;
              }
              //if( best_score_overall < best_scores[ chain_i ] ) {
              //  best_score_overall = best_scores[ chain_i ];
              //  // Argh.
              //  //best_profile_overall = best_profiles[ chain_i ];
              //  best_profile_overall.copyFrom( best_profiles[ chain_i ] );
              //}
            } // End foreach chain_i, report best profiles

            // Calculate per-chain averages
            m_averageProfile.zero();
            average_log_score = 0.0;
            for( chain_i = 0; chain_i < m_samplingParameters.numChains; chain_i++ ) {
              // TODO: REMOVE!
              average_profiles[ chain_i ] /= ( max_iters - burn_in_iters );
              //average_profiles[ chain_i ].normalize( 0 );
              if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
                cout << "[" << gibbs_type_label[ gibbs_type ] << ", Chain " << chain_i << "] The average profile is " << average_profiles[ chain_i ];
              }
              score =
                m_dynamic_programming.forward_score( 
                  m_samplingParameters,
                  average_profiles[ chain_i ],
                  m_sequences,
                  m_sequence_count,
                  m_forward_matrices[ chain_i ]
                );
              if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
                cout << "\tIts score is " << score << endl;
              }
              if( best_score_overall < score ) {
                best_score_overall = score;
                // Argh.
                //best_profile_overall = average_profiles[ chain_i ];
                best_profile_overall.copyFrom( average_profiles[ chain_i ] );
              }
              m_averageProfile += average_profiles[ chain_i ];
            
              average_log_scores[ chain_i ] /= ( max_iters - burn_in_iters );
              if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
                cout << "[" << gibbs_type_label[ gibbs_type ] << ", Chain " << chain_i << "] The average log score is " << average_log_scores[ chain_i ] << endl;
                cout << endl;
              }
              average_log_score += average_log_scores[ chain_i ];
            } // end foreach chain_i, calculate per-chain averages

            // Calculate the overall average profile and its score
            // TODO: REMOVE!
            m_averageProfile /= m_samplingParameters.numChains;
            //m_averageProfile.normalize( 0 );
            if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
              cout << "[" << gibbs_type_label[ gibbs_type ] << "] The overall average profile is " << m_averageProfile;
            }
            m_averageProfileScore =
              m_dynamic_programming.forward_score( 
                m_samplingParameters,
                m_averageProfile,
                m_sequences,
                m_sequence_count,
                m_forward_matrices[ 0 ]
              );
            if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
              cout << "\tIts score is " << m_averageProfileScore << endl;
            }
            if( best_score_overall < m_averageProfileScore ) {
              best_score_overall = m_averageProfileScore;
              best_profile_overall.copyFrom( m_averageProfile );
            }
            
            average_log_score /= m_samplingParameters.numChains;
            if( m_samplingParameters.numChains > 1 ) {
              if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
                cout << "[" << gibbs_type_label[ gibbs_type ] << "] The overall average log score is " << average_log_score << endl;
                cout << endl;
              }
            }
            
            // Calculate score variances
            if( ( max_iters - burn_in_iters ) > 0 ) {
              if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
                cout << endl;
              }
              log_score_variance = 0.0;
              average_within_chain_log_score_variance = 0.0;
              between_chain_log_score_variance = 0.0;
              for( chain_i = 0; chain_i < m_samplingParameters.numChains; chain_i++ ) {
                log_score_variances[ chain_i ] = 0.0;
                for( iter_j = burn_in_iters; iter_j < max_iters; iter_j++ ) {
                  log_score_variances[ chain_i ] +=
                    std::pow( ( log_scores[ chain_i ][ iter_j - burn_in_iters ] - average_log_scores[ chain_i ] ), 2 );
                  log_score_variance +=
                    std::pow( ( log_scores[ chain_i ][ iter_j - burn_in_iters ] - average_log_score ), 2 );
                }
                // To get an unbiased estimator of the true variance, we need to
                // divide by one fewer iteration:
                if( ( max_iters - burn_in_iters ) > 1 ) {
                  log_score_variances[ chain_i ] /= ( ( max_iters - burn_in_iters ) - 1 );
                  if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
                    cout << "[" << gibbs_type_label[ gibbs_type ] << ", Chain " << chain_i << "] The log-score variance is " << log_score_variances[ chain_i ] << endl;
                    cout << "\tLog-score stdev is " << std::pow( log_score_variances[ chain_i ], .5 ) << endl;
                    cout << endl;
                  }

                  average_within_chain_log_score_variance +=
                    log_score_variances[ chain_i ];

                  log_score_stdev =
                    std::pow( log_score_variances[ chain_i ], .5 );
                  if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
                    cout << endl;
                  }
                } // End if ( max_iters - burn_in_iters > 1 )

                between_chain_log_score_variance +=
                  std::pow( ( average_log_scores[ chain_i ] - average_log_score ), 2 );
            
              } // End foreach chain_i, calculate variances
            
              // To get an unbiased estimator of the true variance, we need to
              // divide by one fewer iteration:
              if( ( max_iters - burn_in_iters ) > 1 ) {
                log_score_variance /=
                  ( ( m_samplingParameters.numChains * ( max_iters - burn_in_iters ) ) - 1 );
                if( m_samplingParameters.numChains > 1 ) {
                  if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
                    cout << "[" << gibbs_type_label[ gibbs_type ] << "] The overall log-score variance is " << log_score_variance << endl;
                    cout << "\tLog-score stdev is " << std::pow( log_score_variance, .5 ) << endl;
                    cout << endl;
                  }
              
                  average_within_chain_log_score_variance /= m_samplingParameters.numChains;
                  if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
                    cout << "[" << gibbs_type_label[ gibbs_type ] << "] The average within-chain log-score variance is " << average_within_chain_log_score_variance << endl;
                    //cout << "\tIts stdev is " << std::pow( average_within_chain_log_score_variance, .5 ) << endl;
                    cout << endl;
                  }

                  between_chain_log_score_variance /=
                    ( m_samplingParameters.numChains - 1 );
                  if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
                    cout << "[" << gibbs_type_label[ gibbs_type ] << "] The between-chain log-score variance is " << between_chain_log_score_variance << endl;
                    //cout << "\tIts stdev is " << std::pow( between_chain_log_score_variance, .5 ) << endl;
                    cout << endl;
                  }

                  // See Gelman, et al. page 296, equation 11.3.  This is an
                  // overestimate of the true marginal posterior variance of the
                  // estimand, as opposed to the
                  // average_within_chain_log_score_variance, which is an
                  // underestimate (both approach the true variance in the limit of
                  // many iterations).
                  var_hat_plus =
                    (
                      ( ( ( max_iters - burn_in_iters - 1.0 ) / ( double )( max_iters - burn_in_iters ) ) * average_within_chain_log_score_variance ) +
                      between_chain_log_score_variance
                    );
                  if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
                    cout << "[" << gibbs_type_label[ gibbs_type ] << "] The R-hat numerator (an overestimate of the true variance), is " << var_hat_plus << endl;
                    //cout << "\tIts stdev is " << std::pow( var_hat_plus, .5 ) << endl;
                    cout << endl;
                  }

                  // V_hat is discussed in "General methods for monitoring
                  // convergence of iterative simulations" by Brooks and Gelman
                  // (1998) (it is not in the Gelman book cited above).  It is
                  // the "pooled posterior variance estimate."
                  V_hat =
                    var_hat_plus + ( between_chain_log_score_variance / m_samplingParameters.numChains );
                  if( true || be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
                    cout << "[" << gibbs_type_label[ gibbs_type ] << "] V-hat, the R-hat numerator (the pooled posterior variance estimate), is " << V_hat << endl;
                    //cout << "\tIts stdev is " << std::pow( V_hat, .5 ) << endl;
                    cout << endl;
                  }

                  // R_hat is the ratio of the overestimated variance and the
                  // underestimated variance.  See Gelman et al. page 297.  Or
                  // see Brooks and Gelman 1998, equation 1.1.
                  R_hat =
                    ( V_hat / average_within_chain_log_score_variance );
                  if( true || be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
                    cout << "[" << gibbs_type_label[ gibbs_type ] << "] R-hat, which is the ratio of the pooled posterior variance estimate to the average within-chain variance, is " << R_hat << endl;
                  }
                  // Same calculation, a different way:
                  //R_hat =
                  //  (
                  //    ( ( ( m_samplingParameters.numChains + 1.0 ) / m_samplingParameters.numChains ) * ( var_hat_plus / average_within_chain_log_score_variance ) ) -
                  //    ( ( max_iters - burn_in_iters - 1.0 ) / ( double )( m_samplingParameters.numChains * ( max_iters - burn_in_iters ) ) )
                  //  );
                  //cout << "[" << gibbs_type_label[ gibbs_type ] << "] R-hat, which is the ratio of the pooled posterior variance estimate to the average within-chain variance, is " << R_hat << endl;
                  //cout << endl;

                  if( R_hat < m_samplingParameters.minRHat ) {
                    done = true;
                  } else {
                    done = false;
                  }
                } // End if m_samplingParameters.numChains > 1
              } else { // if ( max_iters - burn_in_iters ) > 1 .. else ..
                done = false;
              } // End if ( max_iters - burn_in_iters ) > 1 .. else ..
            } else { // if ( max_iters - burn_in_iters ) > 0 .. else ..
              done = false;
            } // End if ( max_iters - burn_in_iters ) > 0 .. else ..

            if( !done && ( max_iters == m_samplingParameters.maxGibbsIterations ) ) {
              // Uh-oh.  We have to stop now.
              // TODO: ?
              //if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
                cout << "Maximum iterations reached.  Stopping even though the chains have not yet converged." << endl;
                //}
              done = true;
            }

            if( !done ) {

              // Unless we store every profile, we have a problem with changing
              // the number of burn in iterations.  For now what we do is call
              // all iterations up to this point 'burn-in' iterations.

              // Start over...
              for( chain_i = 0; chain_i < m_samplingParameters.numChains; chain_i++ ) {
                average_profiles[ chain_i ].zero();
                average_log_scores[ chain_i ] = 0.0;
                if( dont_zero_ententes_hack ) {
                  position_ententes[ chain_i ].zero();
                  global_ententes[ chain_i ].zero();
                }
              } // End foreach chain_i, start over.
              
              // Now increment max_iters.
              max_iters += iterations_increment;
              if( max_iters > m_samplingParameters.maxGibbsIterations ) {
                max_iters = m_samplingParameters.maxGibbsIterations;
              }
              burn_in_iters = iter_i;

            } // If we're not done, undo some stuff...

          } while( !done );

          // Calculate the effective number of independent draws (eq. 11.4 of
          // Gelman et. al).
          double n_eff =
            (
              ( m_samplingParameters.numChains * ( max_iters - burn_in_iters ) * var_hat_plus ) /
              between_chain_log_score_variance
            );
          if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
            cout << endl;
            cout << "[" << gibbs_type_label[ gibbs_type ] << "] Number of draws is " << ( m_samplingParameters.numChains * ( max_iters - burn_in_iters ) ) << endl;
            cout << "[" << gibbs_type_label[ gibbs_type ] << "] The effective number of independent draws is " << n_eff << endl;
          }

          if( max_lag > 0 ) {
            // Calculate autocorrelations
            for( chain_i = 0; chain_i < m_samplingParameters.numChains; chain_i++ ) {
              log_score_autocorrelation_first_lag_below_thresholds[ chain_i ] = 0;
              if( ( max_iters - burn_in_iters ) > 1 ) {
                for( lag = 1; lag <= max_lag; lag++ ) {
                  autocorrelation =
                    calculateAutocorrelation(
                      log_scores[ chain_i ],
                      average_log_scores[ chain_i ],
                      average_log_scores[ chain_i ],
                      log_score_stdev,
                      log_score_stdev,
                      0,
                      ( max_iters - burn_in_iters - 1 ),
                      lag
                    );
                  if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
                    cout << "[" << gibbs_type_label[ gibbs_type ] << ", Chain " << chain_i << "] The log-score lag-" << lag << "-autocorrelation is " << autocorrelation << endl;
                  }
                  if( autocorrelation < autocorrelation_threshold ) {
                    if( log_score_autocorrelation_first_lag_below_thresholds[ chain_i ] == 0 ) {
                      log_score_autocorrelation_first_lag_below_thresholds[ chain_i ] = lag;
                    }
                    // TODO: REMOVE to see all of the autocorrelations up to
                    // max_lag.
                    break;
                  }
                } // End foreach lag
              } // End if( ( max_iters - burn_in_iters ) > 1 )
            } // End foreach chain_i

            if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
              for( chain_i = 0; chain_i < m_samplingParameters.numChains; chain_i++ ) {
                cout << "[" << gibbs_type_label[ gibbs_type ] << ", Chain " << chain_i << "] First lag at which the autocorrelation is < " << autocorrelation_threshold << ": " << log_score_autocorrelation_first_lag_below_thresholds[ chain_i ] << endl;
              } // End foreach chain_i
            }
          } // End if max_lag > 0
          // TODO: Thin?

       //} // End foreach gibbs_type (Simple == 0 and Per Position == 1)

        if(
          m_samplingParameters.saveGibbsMode &&
          ( best_score_overall < m_bestProfileScore )
        ) {
          best_score_overall = m_bestProfileScore;
          best_profile_overall.copyFrom( m_bestProfile );
        }

        if( be_extra_verbose || ( m_samplingParameters.verbosity >= VERBOSITY_Low ) ){
          cout << endl;
          cout << "The best profile overall is " << best_profile_overall;
          cout << "\tIts score is " << best_score_overall << endl;
        }

        // ARGH!  This is not working as it should (it isn't calling the
        // operator= that I defined, so it doesn't set the m_root pointers of
        // the ProfilePositions properly).  I've added fixing this to the TODO
        // Bugs list.
        // m_samplingProfile = best_profile_overall;
        m_samplingProfile.copyFrom( best_profile_overall );

      // ENDMARK
      return best_score_overall;
    } // sample()

} // End namespace galosh

#endif // __GALOSH_PROFILEGIBBS_HPP__
