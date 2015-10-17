/*---------------------------------------------------------------------------##
##  Library:
##      galosh::profillicSimulation
##  File:
##      ProfillicSimulation.hpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      Class definition for the ProfillicSimulation class, which implements the
##      profillic simulations for evaluating the methods.
##
#******************************************************************************
#*
#*    This file is part of profillic, a suite of programs for estimating parameters of
#*    Profile HMMs.  Please see the document CITING, which should have been
#*    included with this file.  You may use at will, subject to the license
#*    (Apache v2.0), but *please cite the relevant papers* in your documentation
#*    and publications associated with uses of this library.  Thank you!
#*
#*    Copyright (C) 2009, 2011 by Paul T. Edlefsen, Fred Hutchinson Cancer
#*    Research Center.
#*
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 *    
 *        http://www.apache.org/licenses/LICENSE-2.0
 *    
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
#*****************************************************************************/


#if     _MSC_VER > 1000
#pragma once
#endif

#ifndef __GALOSH_PROFUSETEST_HPP__
#define __GALOSH_PROFUSETEST_HPP__

#include "Profillic.hpp"

#include "Parameters.hpp"
using galosh::Parameters;
using galosh::DebugLevel;
using galosh::VerbosityLevel;

#include "DynamicProgramming.hpp"
using galosh::DynamicProgramming;

#include "ProlificParameters.hpp"
using galosh::ProlificParameters;

#include "Profile.hpp"
using galosh::ProfileTreeRoot;
using galosh::ProfileTreeInternalNode;

#include "ProfileTree.hpp"
using galosh::ProfileTree;

#include "ProfileTreeTrainer.hpp"
using galosh::ProfileTreeTrainer;

//#include "ProfileGibbs.hpp"
//using galosh::ProfileGibbs;

#include "ProfileTrainer.hpp"
using galosh::ProfileTrainer;

#include "Sequence.hpp"
using galosh::Sequence;

#include "Random.hpp"
using galosh::Random;

#include <string>
using std::string;
#include <iostream>
using std::cout;
using std::endl;
//#include <cmath> // for isnan( double )
//using std::isnan; // doesn't work !  why?
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

#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;
#include "CommandlineParameters.hpp"

#include <boost/lexical_cast.hpp>

/**
 * \class ProfillicSimulation
 * \brief Public / commandline interface to ProfillicSimulation.
 *
 */
namespace galosh {

template <class ResidueType,
          class ProbabilityType,
          class ScoreType,
          class MatrixValueType,
          class SequenceResidueType>
  class ProfillicSimulation {
  public:
    typedef Sequence<SequenceResidueType> SequenceType;

    class Parameters :
      //public ProfileGibbs<ProfileTreeRoot<ResidueType, ProbabilityType>,ScoreType,MatrixValueType,SequenceResidueType>::Parameters
        public ProfileTreeTrainer<ResidueType,ProbabilityType,ScoreType,MatrixValueType,SequenceResidueType>::Parameters
    {
    public:
       po::options_description m_profusetest_options;

    private:
      //typedef typename ProfileGibbs<ProfileTreeRoot<ResidueType, ProbabilityType>,ScoreType,MatrixValueType,SequenceResidueType>::Parameters profile_gibbs_parameters_t;
      typedef typename ProfileTreeTrainer<ResidueType,ProbabilityType,ScoreType,MatrixValueType,SequenceResidueType>::Parameters profile_tree_trainer_parameters_t;
      // Boost serialization
      friend class boost::serialization::access;
      template<class Archive>
      void serialize ( Archive & ar, const unsigned int /* file_version */ )
      {
        // save/load base class information
//        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( profile_gibbs_parameters_t );
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( profile_tree_trainer_parameters_t );

        /**
         * ProfillicSimulation Members to be serialized
         *   TAH 9/13
         **/
        #undef GALOSH_DEF_OPT
        #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) ar & BOOST_SERIALIZATION_NVP( NAME )
        #include "ProfillicSimulationOptions.hpp"  /// serialize ProfillicSimulation parameters

      } // serialize( Archive &, const unsigned int )

    public:
  
      /**
       * Define ProfillicSimulation::Parameters "members".  These are tightly tied to the options.
       *    TAH 9/13
       **/
      #undef GALOSH_DEF_OPT
      #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) TYPE NAME
      #include "ProfillicSimulationOptions.hpp"  /// declare Parameters members specific to ProfillicSimulation

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
       #include "ProfillicSimulationOptions.hpp"
#undef GALOSH_DEF_OPT

    }; // End inner class Parameters

    template <class ParametersType>
    class ParametersModifierTemplate :
    //public ProfileGibbs<ProfileTreeRoot<ResidueType, ProbabilityType>,ScoreType,MatrixValueType,SequenceResidueType>::template ParametersModifierTemplate<ParametersType>
public ProfileTreeTrainer<ResidueType,ProbabilityType,ScoreType,MatrixValueType,SequenceResidueType>::template ParametersModifierTemplate<ParametersType>
    {
      //typedef typename ProfileGibbs<ProfileTreeRoot<ResidueType, ProbabilityType>,ScoreType,MatrixValueType,SequenceResidueType>::template ParametersModifierTemplate<ParametersType> base_parameters_modifier_t; 
      typedef typename ProfileTreeTrainer<ResidueType,ProbabilityType,ScoreType,MatrixValueType,SequenceResidueType>::template ParametersModifierTemplate<ParametersType> base_parameters_modifier_t; 

      // Boost serialization
    private:
      friend class boost::serialization::access;
      template<class Archive>
      void serialize ( Archive & ar, const unsigned int /* file_version */ )
      {
        // save/load base class information.  This will serialize the
        // parameters too.
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( base_parameters_modifier_t );

        // Serialize the ProfillicSimulation::ParameterModifierTemplate specific isModified_<member>s TAH 9/13
        #undef GALOSH_DEF_OPT
        #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) ar & BOOST_SERIALIZATION_NVP(isModified_##NAME)
        #include "ProfillicSimulationOptions.hpp"  //archive ProfillicSimulation::ParameterModifierTemplate isModified_<member>s

      } // serialize( Archive &, const unsigned int )

    public:
  
      /// isModified flags for Parameters
      /**
       * Declare isModified_<member>s for isModified_<member>s of ProfileTrainer::ParameterModifierTemplate
       * TAH 9/13
       **/
      #undef GALOSH_DEF_OPT
      #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) bool isModified_##NAME
      #include "ProfillicSimulationOptions.hpp"  // Declare isModified<member>s for ProfillicSimulation::ParametersModifierTemplate

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

    typedef ParametersModifierTemplate<typename ProfillicSimulation::Parameters> ParametersModifier;

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
      //typename ProfileGibbs<RootType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifier parametersModifier;
      typename ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifier parametersModifier;

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
    ProfillicSimulation ();

    ///TAH 6/12 constructor from commandline options
   ProfillicSimulation ( int argc, char **argv );

    /**
     * Construct a profuse test object, using the provided seed.
     */  
    ProfillicSimulation (
      uint32_t const seed
    );

    void
    start ();

  }; // End class ProfillicSimulation

  //======//// potentially non-inline implementations ////========//


  ////// Class galosh::ProfillicSimulation::Parameters ////
  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_INIT
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      Parameters ()
      {
         #ifdef DEBUG
         cout << "[debug] ProfillicSimulation::Parameters::<init>()" << endl;
         cout << "[debug] using ProfillicSimulationOptions.hpp" << endl;
         #endif
         /**
          *  Describe all options/parameters to the options_description object.  In the main
          *  routine (whatever it may be) these descriptions will be parsed from the commandline
          *  and possibly other sources.  This Parameters initializer calls it's parent to get
          *  the ProlificParameter::Parameters options.
          **/
         #undef GALOSH_DEF_OPT
#define PROFUSETEST_DEFAULT_TMP_ARRAY_TO_VECTOR
         #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP)          \
         this->galosh::Parameters::m_galosh_options_description.add_options()(#NAME,po::value<TYPE>(&NAME)->default_value(DEFAULTVAL) TMP_EXTRA_STUFF,HELP)
         #include "ProfillicSimulationOptions.hpp"  /// define all the commandline options for this module

      } // <init>()

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class AnyParameters>
  GALOSH_INLINE_INIT
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      // Copy constructor
      Parameters ( const AnyParameters & copy_from )
      {
        //if( static_cast<galosh::Parameters>( copy_from ).debug >= DEBUG_All ) {
        //  cout << "[debug] ProfillicSimulation::Parameters::<init>( copy_from )" << endl;
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
  typename ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters &
      // Copy constructor/operator
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      operator= (
        const AnyParameters & copy_from
      )
      {
        if( copy_from.debug >= DEBUG_All ) {
          cout << "[debug] ProfillicSimulation::Parameters::operator=( copy_from )" << endl;
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
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      copyFromNonVirtual (
        AnyParameters const & copy_from
      )
      {
        //ProfileGibbs<ProfileTreeRoot<ResidueType, ProbabilityType>,ScoreType,MatrixValueType,SequenceResidueType>::Parameters::copyFromNonVirtual( copy_from );
        ProfileTreeTrainer<ResidueType,ProbabilityType,ScoreType,MatrixValueType,SequenceResidueType>::Parameters::copyFromNonVirtual( copy_from );
        //if( copy_from.debug >= DEBUG_All ) {
        //  cout << "[debug] ProfillicSimulation::Parameters::copyFromNonVirtual( copy_from )" << endl;
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
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      copyFromNonVirtualDontDelegate (
        AnyParameters const & copy_from
      )
      {
    /// TAH 9/13
    #undef GALOSH_DEF_OPT
    #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) NAME = copy_from. NAME
    #include "ProfillicSimulationOptions.hpp"  /// copy all ProfillicSimulation::Parameters members

      } // copyFromNonVirtualDontDelegate( AnyParameters const & )


  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_TRIVIAL
      void
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
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
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      resetToDefaults ()
      {
        //ProfileGibbs<ProfileTreeRoot<ResidueType, ProbabilityType>,ScoreType,MatrixValueType,SequenceResidueType>::Parameters::resetToDefaults();
        ProfileTreeTrainer<ResidueType,ProbabilityType,ScoreType,MatrixValueType,SequenceResidueType>::Parameters::resetToDefaults();
        #ifdef DEBUG
          cout << "[debug] ProfillicSimulation::Parameters::resetToDefaults()" << endl;
        #endif

        /// TAH 9/13
        #undef GALOSH_DEF_OPT
        #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) NAME = this->Parameters::m_galosh_options_map[#NAME].template as<TYPE>()
        #include "ProfillicSimulationOptions.hpp"  /// reset all Parameters members
                                              /// (ProfillicSimulation and through inheritance tree)

      } // resetToDefaults()

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
  void
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters::
      writeParameters (
        std::basic_ostream<CharT,Traits>& os
      ) const
      {
        //ProfileGibbs<ProfileTreeRoot<ResidueType, ProbabilityType>,ScoreType,MatrixValueType,SequenceResidueType>::Parameters::writeParameters( os );
        ProfileTreeTrainer<ResidueType,ProbabilityType,ScoreType,MatrixValueType,SequenceResidueType>::Parameters::writeParameters( os );
        os << endl;

        //Note: we must comment out [ProfillicSimulation] because it means something special to
        //in configuration files sensu program_options
        os << "#[ProfillicSimulation]" << endl;

        /**
         * write out all ProfillicSimulation specific parameters in the style of a configuration
         * file, so that program_options parsers can read it back in
         *   TAH 9/13
         **/
        #undef GALOSH_DEF_OPT
        #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) os << #NAME << " = " << lexical_cast<string>(NAME) << endl
        #include "ProfillicSimulationOptions.hpp"  /// write all ProfillicSimulation::Parameters members to os

      } // writeParameters ( basic_ostream & )

  ////// Class galosh::ProfillicSimulation::ParametersModifierTemplate ////
  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class ParametersType>
  GALOSH_INLINE_INIT
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      ParametersModifierTemplate ()
      {
        if( base_parameters_modifier_t::parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProfillicSimulation::ParametersModifierTemplate::<init>()" << endl;
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
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      // Copy constructor
      ParametersModifierTemplate ( const AnyParametersModifierTemplate & copy_from )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProfillicSimulation::ParametersModifierTemplate::<init>( copy_from )" << endl;
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
  typename ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::template ParametersModifierTemplate<ParametersType> &
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      // Copy constructor/operator
  operator= (
        const AnyParametersModifierTemplate & copy_from
      )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProfillicSimulation::ParametersModifierTemplate::operator=( copy_from )" << endl;
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
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProfillicSimulation::ParametersModifierTemplate::copyFromNonVirtual( copy_from )" << endl;
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
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      isModified_copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      )
      {
        /// Copy all parent isModified_<member>s
        base_parameters_modifier_t::isModified_copyFromNonVirtual( copy_from );

        /// TAH 9/13 copy all ProfillicSimulation::isModified_<member>s
        #undef GALOSH_DEF_OPT
        #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) isModified_##NAME = copy_from.isModified_##NAME
        #include "ProfillicSimulationOptions.hpp"  // Copy ProfillicSimulation::ParametersModifierTemplate::isModified_<member>s

      } // isModified_copyFromNonVirtual( AnyParametersModifierTemplate const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  template <class ParametersType>
  GALOSH_INLINE_TRIVIAL
  void
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
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
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      isModified_reset ()
      {
        /// Set all parent ParametersModifierTemplate::isModified_<member>s to false
        base_parameters_modifier_t::isModified_reset();
        /// TAH 9/13 set all ProfillicSimulation::ParametersModifierTemplate::isModified_<member>s to false
        #undef GALOSH_DEF_OPT
        #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) isModified_##NAME = false;
        #include "ProfillicSimulationOptions.hpp" // Reset ProfillicSimulation::ParametersModifierTemplate::isModified_<member>s to false

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
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      writeParametersModifier (
        std::basic_ostream<CharT,Traits>& os
      ) const
      {
        base_parameters_modifier_t::writeParametersModifier( os );
        os << endl;

        /// TAH 9/13 must comment out tags in square braces for program_options config file parser
        os << "#[ProfillicSimulation]" << endl;
        /**
         * write out ProfillicSimulation::parameters iff they've been modified
         *   TAH 9/13
         **/
        #undef GALOSH_DEF_OPT
        #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) if( isModified_##NAME ) os << #NAME << " = " << lexical_cast<string>(galosh::ParametersModifierTemplate<ParametersType>::parameters. NAME)
        #include "ProfillicSimulationOptions.hpp" // write out changed ProfillicSimulationParameters::ParametersModifierTemplate parameters

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
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::ParametersModifierTemplate<ParametersType>::
      applyModifications ( AnyParameters & target_parameters )
      {
        /// Set the parameters of another object iff they've been changed in this one
        base_parameters_modifier_t::applyModifications( target_parameters );

        /**
         * Set the parameters of a foreign Parameters object to this Parameter object's values
         * iff they have changed.
         *    TAH 9/13
         **/
        #undef GALOSH_DEF_OPT
        #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) if( isModified_##NAME ) target_parameters. NAME = this->parameters. NAME
        #include "ProfillicSimulationOptions.hpp" // copy changed parameters

      } // applyModifications( Parameters & )

  ////// Class galosh::ProfillicSimulation ////
  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_INIT
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::
    ProfillicSimulation (
    ) :
      m_parameters(),
      m_random( static_cast<uint32_t>( std::time( NULL ) ) )
    {
      if( m_parameters.debug >= DEBUG_All ) {
        cout << "[debug] ProfillicSimulation::<init>()" << endl;
      } // End if DEBUG_All
      // Do nothing else
    } // <init>()

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_INIT
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::
  /**
   * Construct a profuse test object, using the given command-line arguments.
   */  
  ProfillicSimulation ( const int argc, char ** const argv ) :
      m_parameters(),
      m_random( static_cast<uint32_t>( std::time( NULL ) ) )
  {
    try {
      string config_file;
      // Declare a group of options that will be 
      // allowed only on command line
      po::options_description generic( "Generic options" );
      generic.add_options()
        ( "version,v", "print version string" )
        ( "help,h", "produce help message" )
        ( "config,c", po::value<string>( &config_file )->default_value( "ProfillicSimulation.cfg" ),
          "name of a file of a configuration." )
        ;
      
      // Declare a group of options that will be 
      // allowed both on command line and in
      // config file
      //po::options_description config( "Configuration" );
      //config.add_options()
      //  ;
      
      typedef ProfileTreeRoot<ResidueType, ProbabilityType> ProfileType;
      typedef ProfileTreeRoot<ResidueType, ProbabilityType> InternalNodeType;

      po::options_description cmdline_options;
      //cmdline_options.add( generic ).add( m_parameters.m_galosh_options_description ).add( config ).add( lengthadjust_opts );
      cmdline_options.add( generic ).add( m_parameters.m_galosh_options_description );
      
      po::options_description config_file_options;
      //config_file_options.add( m_parameters.m_galosh_options_description ).add( config ).add( lengthadjust_opts );
      config_file_options.add( m_parameters.m_galosh_options_description );
      
      po::options_description visible( "Basic options" );
      //visible.add( generic ).add( config );
      visible.add( generic );
      
      po::positional_options_description p;
      p.add( "config", 1 );
          
      store( po::command_line_parser( argc, argv ).options( cmdline_options ).positional( p ).run(), m_parameters.m_galosh_options_map );
      notify( m_parameters.m_galosh_options_map );
      
      // TODO: REMOVE
      //cout << m_parameters << endl;
      //cout << endl;
      
#define USAGE() " " << argv[ 0 ] << " [options] [<config file>]"
      
          
//      cout << "Usage: " << argv[ 0 ] << " <config file>" << endl;
      
      // Read in the config file.
      if( config_file.length() > 0 ) {
      // TODO: REMOVE
        cout << "Config file is: " << m_parameters.m_galosh_options_map["config"].template as<string>() << endl;
        ifstream ifs( config_file.c_str() );
        if( !ifs ) {
          //if(!m_parameters.m_galosh_options_map["config"].defaulted()) {         //TAH 3/13 don't choke if config file was defaulted and is missing
             cout << "Can't open the config file named \"" << config_file << "\"\n";
             return;
             //} 
        } else {
          store( parse_config_file( ifs, config_file_options ), m_parameters.m_galosh_options_map );
          notify( m_parameters.m_galosh_options_map );
        }
      }
      
      if( m_parameters.m_galosh_options_map.count( "help" ) > 0 ) {
        cout << "Usage: " << USAGE() << endl;
        cout << visible << "\n";
        return;
      }
      
      if( m_parameters.m_galosh_options_map.count( "version" ) ) {
        cout << "ProfillicSimulation, version 1.2\n";
        return;
      }
      
      //if( m_parameters.m_galosh_options_map.count( "debug" ) ) {
      //  cout << "[DEBUGGING]\n";
      //  return;
      //}
      
      // Required options
      if( m_parameters.m_galosh_options_map.count( "config" ) == 0 ) {
        cout << "Usage: " << USAGE() << endl;
        return;
      }
      
      if( m_parameters.m_galosh_options_map.count( "seed" ) ) {
        if( m_parameters.seed != 0 ) {
          m_random.setSeed( m_parameters.seed );
        }
      }
      // TODO: REMOVE
      //cout << "Hi, ok, done constructing" << endl;

    } catch( std::exception& e ) { /// exceptions thrown by boost stuff (etc)
      cerr << "error: " << e.what() << endl;
    } catch( string &err ) {      /// exceptions thrown as strings
      cerr << "error: " << err << endl;
    } catch( ... ) {               /// anything else
      cerr << "Strange unknown exception" << endl;
    }

  } // <init>( const int argc, char ** const argv )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_INIT
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::
    /**
     * Construct a profuse test object, using the provided seed.
     */  
    ProfillicSimulation (
      uint32_t const seed
    ) :
      m_parameters(),
      m_random( seed )
    {
      if( m_parameters.debug >= DEBUG_All ) {
        cout << "[debug] ProfillicSimulation::<init>( uint32_t )" << endl;
      } // End if DEBUG_All
      // Do nothing else
    } // <init>( uint32_t )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType>
  GALOSH_INLINE_TRAIN
  void
  ProfillicSimulation<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::
    start ()
    {
      bool be_verbose = ( m_parameters.verbosity > VERBOSITY_Low );

      // TODO: REMOVE "true"; DEBUGGING.
      if( true || ( m_parameters.verbosity > VERBOSITY_High ) ) {
        // Print out the parameters
        for( const auto& it : ( m_parameters.m_galosh_options_map ) ) {
          std::cout << it.first.c_str() << " = ";
          auto& value = it.second.value();
          if( auto v = boost::any_cast<bool>(&value) ) {
            std::cout << *v;
          } else if( auto v = boost::any_cast<double>(&value) ) {
            std::cout << *v;
          } else if( auto v = boost::any_cast<float>(&value) ) {
            std::cout << *v;
          } else if( auto v = boost::any_cast<unsigned int>(&value) ) {
            std::cout << *v;
          } else if( auto v = boost::any_cast<int>(&value) ) {
            std::cout << *v;
          } else if( auto v = boost::any_cast<std::string>(&value) ) {
            std::cout << *v;
          } else if( auto v = boost::any_cast<myVector<double> >(&value) ) {
            std::cout << *v;
          } else if( auto v = boost::any_cast<myVector<uint32_t> >(&value) ) {
            std::cout << *v;
          } else {
            std::cout << "error";
          }
          std::cout << endl;
        }
        std::cout << endl;
      } // End if we should print out all of the parameters

      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DirichletMixtureMatchEmissionPrior<float> matchEmissionPrior;

      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DirichletMixtureGlobalPrior<float> globalPrior;

      // Set up the priors
      if( m_parameters.usePriors || m_parameters.startWithPositionsDrawnFromPrior ) {
        matchEmissionPrior.reinitializeToEven( m_parameters.priorStrength );
      } // End if m_parameters.usePriors || m_parameters.startWithPositionsDrawnFromPrior
      if( m_parameters.usePriors || m_parameters.startWithGlobalsDrawnFromPrior ) {
        globalPrior.reinitializeToEven( m_parameters.priorStrength );
        // NOTE: We will do additional set-up of the global prior for each
        // change in the profile length, since we need to adjust for profile
        // length.

        // TODO: Make these parameters, DEHACKIFY
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

      // TODO: PUT BACK WHEN ProfileGibbs IS READY.
      //if( m_parameters.testConditionalGibbsProfile ) {
      //  tests[ TEST_ID_gibbs_conditional ].name = "cGibbs";
      //  tests[ TEST_ID_gibbs_conditional ].isRun = true;
      //} else {
      //  tests[ TEST_ID_gibbs_conditional ].isRun = false;
      //}
      //if( m_parameters.coutConditionalGibbsProfile ) {
      //  tests[ TEST_ID_gibbs_conditional ].isCout = true;
      //} else {
      //  tests[ TEST_ID_gibbs_conditional ].isCout = false;
      //}
      //tests[ TEST_ID_gibbs_conditional ].isGibbs = true;
      //tests[ TEST_ID_gibbs_conditional ].coutLeftBrace = "`";
      //tests[ TEST_ID_gibbs_conditional ].coutRightBrace = "'";
      //// Start with the starting profile globals and positions
      //tests[ TEST_ID_gibbs_conditional ].startingGlobalsTest =
      //  NULL;
      //tests[ TEST_ID_gibbs_conditional ].startingPositionsTest =
      //  NULL;
      //// Modifiers for conditional gibbs test
      //// inherit these
      ////tests[ TEST_ID_gibbs_conditional ].parametersModifier.parameters.sampleProfilePositions = true;
      ////tests[ TEST_ID_gibbs_conditional ].parametersModifier.isModified_sampleProfilePositions = true;
      ////if( m_parameters.trainProfileGlobals ) {
      ////  tests[ TEST_ID_gibbs_conditional ].parametersModifier.parameters.sampleProfileGlobals = true;
      ////  tests[ TEST_ID_gibbs_conditional ].parametersModifier.isModified_sampleProfileGlobals = true;
      ////} else {
      ////  tests[ TEST_ID_gibbs_conditional ].parametersModifier.parameters.sampleProfileGlobals = false;
      ////  tests[ TEST_ID_gibbs_conditional ].parametersModifier.isModified_sampleProfileGlobals = false;
      ////}
      //
      //tests[ TEST_ID_gibbs_conditional ].parametersModifier.parameters.useUnconditionalGibbs = false;
      //tests[ TEST_ID_gibbs_conditional ].parametersModifier.isModified_useUnconditionalGibbs = true;
      //if( m_parameters.reportGibbsMode ) {
      //  tests[ TEST_ID_gibbs_conditional ].parametersModifier.parameters.saveGibbsMode = true;
      //  tests[ TEST_ID_gibbs_conditional ].parametersModifier.isModified_saveGibbsMode = true;
      //} // End if reportGibbsMode
      //
      //if( m_parameters.testUnconditionalGibbsProfile ) {
      //  tests[ TEST_ID_gibbs_unconditional ].name = "uGibbs";
      //  tests[ TEST_ID_gibbs_unconditional ].isRun = true;
      //} else {
      //  tests[ TEST_ID_gibbs_unconditional ].isRun = false;
      //}
      //if( m_parameters.coutUnconditionalGibbsProfile ) {
      //  tests[ TEST_ID_gibbs_unconditional ].isCout = true;
      //} else {
      //  tests[ TEST_ID_gibbs_unconditional ].isCout = false;
      //}
      //tests[ TEST_ID_gibbs_unconditional ].isGibbs = true;
      //tests[ TEST_ID_gibbs_unconditional ].coutLeftBrace = "*";
      //tests[ TEST_ID_gibbs_unconditional ].coutRightBrace = "*";
      //// Start with the starting profile globals and positions
      //tests[ TEST_ID_gibbs_unconditional ].startingGlobalsTest =
      //  // TODO: Put back: NULL;
      //  &tests[ TEST_ID_gibbs_unconditional ]; // true
      //tests[ TEST_ID_gibbs_unconditional ].startingPositionsTest = NULL;
      //// Modifiers for unconditional gibbs test
      //// inherit these, don't modify them.
      ////tests[ TEST_ID_gibbs_unconditional ].parametersModifier.parameters.sampleProfilePositions = true;
      ////tests[ TEST_ID_gibbs_unconditional ].parametersModifier.isModified_sampleProfilePositions = true;
      ////if( m_parameters.trainProfileGlobals ) {
      ////  tests[ TEST_ID_gibbs_unconditional ].parametersModifier.parameters.sampleProfileGlobals = true;
      ////  tests[ TEST_ID_gibbs_unconditional ].parametersModifier.isModified_sampleProfileGlobals = true;
      ////} else {
      ////  tests[ TEST_ID_gibbs_unconditional ].parametersModifier.parameters.sampleProfileGlobals = false;
      ////  tests[ TEST_ID_gibbs_unconditional ].parametersModifier.isModified_sampleProfileGlobals = false;
      ////}
      //tests[ TEST_ID_gibbs_unconditional ].parametersModifier.parameters.useUnconditionalGibbs = true;
      //tests[ TEST_ID_gibbs_unconditional ].parametersModifier.isModified_useUnconditionalGibbs = true;
      //if( m_parameters.reportGibbsMode ) {
      //  tests[ TEST_ID_gibbs_unconditional ].parametersModifier.parameters.saveGibbsMode = true;
      //  tests[ TEST_ID_gibbs_unconditional ].parametersModifier.isModified_saveGibbsMode = true;
      //} // End if reportGibbsMode

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


      uint32_t test_id;
    
      // Now go back through the tests and turn on tests[ test_num ].isRun for
      // any test that is required for another test with tests[ .. ].isRun ==
      // true.
      bool cycle_run_tests_again = true;
      while ( cycle_run_tests_again ) {
        cycle_run_tests_again = false;
        for( test_id = 0; test_id <= LAST_TEST_ID; test_id++ ) {
          // TODO: REMOVE
          //cout << "test_id is " << test_id << endl;
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

      if( m_parameters.verbosity >= VERBOSITY_Meta ) {
        cout << "Training with seed " << m_random.getSeed() << ":" << endl;
      } // End if VERBOSITY_Meta

      // Results
      string run_unique_id =
        ( "v" + boost::lexical_cast<string>( m_parameters.saveFileVersion ) + "_seed" + boost::lexical_cast<string>( m_random.getSeed() ) + "_ProfillicSimulation" );
      fs::path dirname =
        ( static_cast<fs::path>( m_parameters.saveResultsParentDirectory ) /
          run_unique_id );
      std::ofstream tab_stream;
      if( m_parameters.saveResultsToFile ) {
        if( !fs::exists( m_parameters.saveResultsParentDirectory ) ) {
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
        // Parameters file
        fs::path parameters_filename =
           ( m_parameters.resultsFilePrefix +
             run_unique_id +
             m_parameters.parametersFileSuffix );
         ofstream parameters_stream( ( dirname / parameters_filename ).string().c_str() );
        assert( parameters_stream.good() );

        // Print out the parameters
        for( const auto& it : ( m_parameters.m_galosh_options_map ) ) {
          parameters_stream << it.first.c_str() << " = ";
          auto& value = it.second.value();
          if( auto v = boost::any_cast<bool>(&value) ) {
            parameters_stream << *v;
          } else if( auto v = boost::any_cast<double>(&value) ) {
            parameters_stream << *v;
          } else if( auto v = boost::any_cast<float>(&value) ) {
            parameters_stream << *v;
          } else if( auto v = boost::any_cast<unsigned int>(&value) ) {
            parameters_stream << *v;
          } else if( auto v = boost::any_cast<int>(&value) ) {
            parameters_stream << *v;
          } else if( auto v = boost::any_cast<std::string>(&value) ) {
            parameters_stream << *v;
          } else if( auto v = boost::any_cast<myVector<double> >(&value) ) {
            parameters_stream << *v;
          } else if( auto v = boost::any_cast<myVector<uint32_t> >(&value) ) {
            parameters_stream << *v;
          } else {
            parameters_stream << "error";
          }
          parameters_stream << endl;
        }
        parameters_stream << endl;
        parameters_stream.close();

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
      vector<ScoreType> testScore_training_forward( last_test_id + 1 );
  
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
      vector<ScoreType> testScore_training_mean_forward( last_test_id + 1 );

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
      vector<ScoreType> testScore_training_mode_forward( last_test_id + 1 );

      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true, then testScore_test_forward[
       * test_num ] will hold the forward score for the corresponding test
       * number for the test set.
       */
      vector<ScoreType> testScore_test_forward( last_test_id + 1 );
  
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
      vector<ScoreType> testScore_test_mean_forward( last_test_id + 1 );

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
      vector<ScoreType> testScore_test_mode_forward( last_test_id + 1 );

      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true, then testScore_training_viterbi[
       * test_num ] will hold the viterbi score for the corresponding test
       * number for the training set.
       */
      vector<ScoreType> testScore_training_viterbi( last_test_id + 1 );
  
      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true, then testScore_test_viterbi[
       * test_num ] will hold the viterbi score for the corresponding test
       * number for the test set.
       */
      vector<ScoreType> testScore_test_viterbi( last_test_id + 1 );

      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true, then testScore_training_truepath[
       * test_num ] will hold the truepath score for the corresponding test
       * number for the training set.
       */
      vector<ScoreType> testScore_training_truepath( last_test_id + 1 );
  
      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true, then testScore_test_truepath[
       * test_num ] will hold the truepath score for the corresponding test
       * number for the test set.
       */
      vector<ScoreType> testScore_test_truepath( last_test_id + 1 );

      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true, then testProfileTree[ test_num ]
       * will hold the profile tree (after training) for the corresponding test
       * number.
       */
      vector<ProfileTreeType> testProfileTree( last_test_id + 1 );
  
      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true and tests[ test_num ].isGibbs is
       * true, then testProfileTree_mean[ test_num ] will hold the mean profile
       * tree (after sampling) for the corresponding test number.
       *
       * Only used if parameters.reportGibbs_mean is true.
       */
      vector<ProfileTreeType> testProfileTree_mean( last_test_id + 1 );
  
      /**
       * For convenience we number the tests from 0 to last_test_id.
       *
       * If tests[ test_num ].isRun is true and tests[ test_num ].isGibbs is
       * true, then testProfileTree_mode[ test_num ] will hold the mode profile
       * tree (after sampling) for the corresponding test number.
       *
       * Only used if parameters.reportGibbs_mode is true.
       */
      vector<ProfileTreeType> testProfileTree_mode( last_test_id + 1 );
  
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
           conservation_rate_i < ( m_parameters.conservationRates.size() > 0 ? m_parameters.conservationRates.size() : 3U );
           conservation_rate_i++
      ) {
        double conservation_rate =
          ( m_parameters.conservationRates.size()>0 ?
            m_parameters.conservationRates[ conservation_rate_i ] :
            ( .25 * ( conservation_rate_i + 1 ) ) );
        for( uint32_t profile_length_i = 0;
             profile_length_i < ( m_parameters.profileLengths.size() > 0 ? m_parameters.profileLengths.size() : 10U );
             profile_length_i++
        ) {
          uint32_t profile_length =
            ( m_parameters.profileLengths.size() > 0 ?
              m_parameters.profileLengths[ profile_length_i ] :
              ( 100 * ( profile_length_i + 1 ) ) );

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
               num_training_sequences_per_profile_i < ( m_parameters.numTrainingSequencesPerProfiles.size() > 0 ? m_parameters.numTrainingSequencesPerProfiles.size() : 1U );
               num_training_sequences_per_profile_i++
          ) {
            uint32_t num_training_sequences_per_profile =
              ( m_parameters.numTrainingSequencesPerProfiles.size() > 0 ?
                m_parameters.numTrainingSequencesPerProfiles[ num_training_sequences_per_profile_i ] :
                ( uint32_t )std::pow( 100.0, ( int )( num_training_sequences_per_profile_i + 1 ) ) );
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
                      0.1 );
                  for( uint32_t expected_insertion_length_as_profile_length_fraction_i = 0;
                       ( m_parameters.useDeletionsForInsertionsParameters ? ( expected_insertion_length_as_profile_length_fraction_i == 0 ) : ( expected_insertion_length_as_profile_length_fraction_i < (m_parameters.expectedInsertionLengthAsProfileLengthFractions.size() > 0 ? m_parameters.expectedInsertionLengthAsProfileLengthFractions.size() : 1U ) ) );
                       expected_insertion_length_as_profile_length_fraction_i++
                  ) {
                    double expected_insertion_length_as_profile_length_fraction =
                      ( m_parameters.useDeletionsForInsertionsParameters ?
                        expected_deletion_length_as_profile_length_fraction :
                        ( ( m_parameters.expectedInsertionLengthAsProfileLengthFractions.size() > 0 ?
                            m_parameters.expectedInsertionLengthAsProfileLengthFractions[ expected_insertion_length_as_profile_length_fraction_i ] :
                            0.1 ) ) );

                    // Pattern sequences
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
                      // TODO: REMOVE?  This may now be unnecessary.
                      true_root.ensurePositionsKnowTheirRoot();

                      // "true" global transition params.
        
                      // First calculate the appropriate indel values.
                      ProbabilityType deletion_open =
                        ( expected_deletions_count / profile_length );
                      ProbabilityType insertion_open =
                        ( m_parameters.useDeletionsForInsertionsParameters ?
                          deletion_open :
                          ( expected_insertions_count / profile_length ) );
        
                      // [ the EV of a geometric is 1/p, where p is prob of stopping, so if q is the prob of continuing, we want ( 1 - q ) = 1/EV. ]
                      ProbabilityType deletion_extension =
                        ( 1.0 - min( ( 1.0 / ( expected_deletion_length_as_profile_length_fraction * profile_length ) ), ( 1.0 / m_parameters.minExpectedDeletionLength ) ) );
                      ProbabilityType insertion_extension =
                        ( m_parameters.useDeletionsForInsertionsParameters ? deletion_extension : ( 1.0 - min( ( 1.0 / ( expected_insertion_length_as_profile_length_fraction * profile_length ) ), ( 1.0 / m_parameters.minExpectedInsertionLength ) ) ) );

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
                        ( m_parameters.preAlignInsertion );
                      true_root[ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ] =
                        ( 1 ) -
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
                        ( 1 ) -
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
                        ( 1.0 ) -
                        (
                         true_root[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] +
                         true_root[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ]
                        );
#endif // USE_DEL_IN_DEL_OUT .. else ..
                      true_root[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] =
                        insertion_extension;
                      true_root[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ] =
                        ( 1.0 ) -
                        true_root[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ];
                      true_root[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ] =
                        deletion_extension;
                      true_root[ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ] =
                        ( 1.0 ) -
                        true_root[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ];
                    
#ifdef USE_END_DISTRIBUTION
                      // For now we don't use the End distribution ..
                      //true_root[ Transition::fromEnd ][ TransitionFromEnd::toPostAlign ] = ( 1 );
                      //true_root[ Transition::fromEnd ][ TransitionFromEnd::toLoop ] = ( 0 );
#endif // USE_END_DISTRIBUTION
#ifndef DISALLOW_FLANKING_TRANSITIONS
                      true_root[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] =
                        ( m_parameters.postAlignInsertion );
                      true_root[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] =
                        ( 1.0 ) -
                        true_root[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ];
#endif // DISALLOW_FLANKING_TRANSITIONS
              
                        // This is a trick to get the values set correctly:
                        // make an even profile, then set one of them higher than you need it, but
                        // normalize to get the right thing.
                        // NOTE: true_root[ j ][ pattern_sequences[ 0 ][ i ] is the same for all i,j right now.
                  
                        ProbabilityType pattern_trick_value =
                          ( ( conservation_rate == 1.0 ) ? ( 1.0 ) :
                            ( ( ( 1.0 ) - true_root[ 0 ][ Emission::Match ][ pattern_sequences[ 0 ][ 0 ] ] ) *
                              ( ( conservation_rate / ( 1.0 - conservation_rate ) ) ) ) );
                      
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

                      training_forward_rows_1.reinitialize(
                        training_fasta,
                        training_fasta.size()
                      );
                      training_forward_rows_2.reinitialize(
                        training_fasta,
                        training_fasta.size()
                      );

                      testing_forward_rows_1.reinitialize(
                        testing_fasta,
                        testing_fasta.size()
                      );
                      testing_forward_rows_2.reinitialize(
                        testing_fasta,
                        testing_fasta.size()
                      );

                      //ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType> tree_trainer =
                      //  ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>(
                      //    training_fasta,
                      //    training_profile_tree
                      //  );
                      //tree_trainer.m_parameters = m_parameters;
                      //
                      //tree_trainer.train();
                    
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
                            // TODO: ERE I AM IN FIXING DEL_IN_DEL_OUT
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
                          // the indel opens (M->I and M->D) equally probable.
                          // Note for consistency I'm doing this ALWAYS.
                          //if( m_parameters.testLengthadjust ) {
                            ProbabilityType average_indel_open = 1.0;
                            average_indel_open -= starting_root[ Transition::fromMatch ][ TransitionFromMatch::toMatch ];
                            average_indel_open /= 2;
                            starting_root[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] =
                              average_indel_open;
                            starting_root[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] =
                              average_indel_open;
                           //} // End if startWithEqualIndelOpens
                        } // End if trainProfileGlobals
              
                        if( be_verbose ) {    
                          cout << "The starting profile tree is:" << endl;
                          cout << testProfileTree[ TEST_ID_starting ] << endl;
                        } // End if be_verbose
              
                        // Get "starting profile" forward score for the training sequences.
                        // TODO: Put back option to use matrices?
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
                        // TODO: PUT BACK WHEN ProfileGibbs IS READY.
                        //ProfileGibbs<RootType, ScoreType, MatrixValueType, SequenceResidueType> sampler(
                        //  &m_random,
                        //  testProfileTree[ TEST_ID_starting ].getProfileTreeRoot(), // not used except to choose the right template type
                        //  training_fasta
                        //);

                        // Note that we must set the m_profileTree value before
                        // training...
                        ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType> trainer(
                          training_fasta,
                          NULL
                        );
                        // Use our own special parameters; see below where it gets modified, too.
                        trainer.m_parameters = m_parameters;

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

                            // TODO: PUT BACK WHEN ProfileGibbs IS READY.
                            //if( !m_parameters.sampleProfilePositions ) {
                            //  // Make sure that the profile starts with the
                            //  // right position values, since we're expecting
                            //  // them to be the true values.
                            //  testProfileTree[ test_id ].getProfileTreeRoot()->copyPositions(
                            //    *( testProfileTree[ TEST_ID_true ].getProfileTreeRoot() )
                            //  );
                            //}
                            //
                            //sampler.m_profile =
                            //  testProfileTree[ test_id ].getProfileTreeRoot();
                            //sampler.m_parameters = m_parameters;
                            //
                            //tests[ test_id ].parametersModifier.applyModifications( sampler.m_parameters );
                            //
                            //// TODO: REMOVE
                            ////cout << "sampler verbosity level is " << sampler.m_parameters.verbosity << endl;
                            ////if( tests[ test_id ].parametersModifier.isModified_verbosity ) {
                            ////  cout << "it should be modified to " << tests[ test_id ].parametersModifier.parameters.verbosity << endl;
                            ////} else {
                            ////  cout << "it should be unmodified (" << m_parameters.verbosity << ")" << endl;
                            ////}
                            //
                            //if( be_verbose ) {
                            //  cout << "Now (before sampling the " << tests[ test_id ].name << " profile), the profile tree is:" << endl;
                            //  cout << testProfileTree[ test_id ] << endl;
                            //} // End if be_verbose
                            //
                            //testScore_training_forward[ test_id ] =
                            //  sampler.sample();
                            //testIters[ test_id ] = sampler.m_totalIterations;
                            //
                            //// The sampler puts the best profile it found in
                            //// m_samplingProfile.  Note that this might not be the mode or the overall mean (it could be the mean of one of the chains).
                            //testProfileTree[ test_id ].getProfileTreeRoot()->copyFrom( sampler.m_samplingProfile );
                            //
                            //// Also separately store the overall mean and mode
                            //if( m_parameters.reportGibbsMean ) {
                            //  testScore_training_mean_forward[ test_id ] =
                            //    sampler.m_averageProfileScore;
                            //  testProfileTree_mean[ test_id ].getProfileTreeRoot()->copyFrom( sampler.m_averageProfile );
                            //}
                            //if( m_parameters.reportGibbsMode ) {
                            //  testScore_training_mode_forward[ test_id ] =
                            //    sampler.m_bestProfileScore;
                            //  testProfileTree_mode[ test_id ].getProfileTreeRoot()->copyFrom( sampler.m_bestProfile );
                            //}
                            //if( be_verbose ) {
                            //  cout << "Now (after training/sampling the " << tests[ test_id ].name << " profile), the score is " << testScore_training_forward[ test_id ] << ", and the profile tree is:" << endl;
                            //  cout << testProfileTree[ test_id ] << endl;
                            //} // End if be_verbose

                          } else { // if isGibbs .. else ..
                            
                            trainer.m_profileTree = &testProfileTree[ test_id ];

                            // Copy the parameters and then set them again.
                            trainer.m_parameters = m_parameters;
                            tests[ test_id ].parametersModifier.applyModifications( trainer.m_parameters );
                            
                            // TODO: REMOVE
                            //cout << "trainer debug level is " << trainer.m_parameters.debug << endl;
                            //if( tests[ test_id ].parametersModifier.isModified_debug ) {
                            //  cout << "it should be modified to " << tests[ test_id ].parametersModifier.parameters.debug << endl;
                            //} else {
                            //  cout << "it should be unmodified (" << m_parameters.debug << ")" << endl;
                            //}
                            //cout << "trainer verbosity level is " << trainer.m_parameters.verbosity << endl;
                            //if( tests[ test_id ].parametersModifier.isModified_verbosity ) {
                            //  cout << "it should be modified to " << tests[ test_id ].parametersModifier.parameters.verbosity << endl;
                            //} else {
                            //  cout << "it should be unmodified (" << m_parameters.verbosity << ")" << endl;
                            //}
                            
                            if( be_verbose && true ) {
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
