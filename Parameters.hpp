/*---------------------------------------------------------------------------##
##  Library:
##      galosh::prolific
##  File:
##      Parameters.hpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      Class definition for the Galosh Parameters class.
##
#******************************************************************************
#*
#*    This file is part of prolific, a library of useful C++ classes for
#*    working with genomic sequence data and Profile HMMs.  Please see the
#*    document CITING, which should have been included with this file.  You may
#*    use at will, subject to the license (Apache v2.0), but *please cite the
#*    relevant papers* in your documentation and publications associated with
#*    uses of this library.  Thank you!
#*
#*    Copyright (C) 2006, 2008, 2011 by Paul T. Edlefsen, Fred Hutchinson Cancer
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

#ifndef __GALOSH_PARAMETERS_HPP__
#define __GALOSH_PARAMETERS_HPP__

#include "Prolific.hpp"

#include <string>
using std::string;
#include <iostream>
using namespace std;
//using std::ostream;
//using std::iosbase;
//using std::cout;
//using std::endl;

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/version.hpp>

namespace galosh {

  enum DebugLevel {
    DEBUG_None = 0,
    DEBUG_Special = 1,
    DEBUG_Medium = 500,
    DEBUG_All = 1000
  };

  enum VerbosityLevel {
    VERBOSITY_None = 0,
    VERBOSITY_Meta = 5, // lower than low: show only info exterior to training
    VERBOSITY_Low = 10,
    VERBOSITY_Medium = 50,
    VERBOSITY_High = 100,
    VERBOSITY_All = 1000
  };
  
  class Parameters {
    // Boost serialization
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize ( Archive & ar, const unsigned int /* file_version */ )
    {
      ar & BOOST_SERIALIZATION_NVP( debug );
      ar & BOOST_SERIALIZATION_NVP( verbosity );
    } // serialize( Archive &, const unsigned int )

  public:
  
    DebugLevel debug;
#define DEFAULT_debug DEBUG_None

    VerbosityLevel verbosity;
#define DEFAULT_verbosity VERBOSITY_None
  
    Parameters ()
    {
      if( DEFAULT_debug ) {
        cout << "[debug] Parameters::<init>()" << endl;
      }
      resetToDefaults();
    } // <init>()
  
    virtual ~Parameters () {};

    // Copy constructor
    Parameters ( const Parameters & copy_from )
    {
      if( copy_from.debug >= DEBUG_All ) {
        cout << "[debug] Parameters::<init>( copy_from )" << endl;
      }
      copyFrom( copy_from );
    } // <init>( Parameters & )
  
    // Copy constructor/operator
    Parameters & operator= (
      const Parameters & copy_from
    )
    {
      if( copy_from.debug >= DEBUG_All ) {
        cout << "[debug] Parameters::operator=( copy_from )" << endl;
      }
      copyFrom( copy_from );
      return *this;
    } // operator=( const Parameters & )

    void
    copyFromNonVirtual ( const Parameters & copy_from )
    {
      if( copy_from.debug >= DEBUG_All ) {
        cout << "[debug] Parameters::copyFromNonVirtual( copy_from )" << endl;
      }
      debug = copy_from.debug;
      verbosity = copy_from.verbosity;
    } // copyFromNonVirtual( Parameters const & )
  
    virtual void
    copyFrom ( const Parameters & copy_from )
    {
      copyFromNonVirtual( copy_from );
    } // copyFrom( Parameters & )
  
    virtual void
    resetToDefaults ()
    {
      debug = DEFAULT_debug;
      if( debug ) {
        cout << "[debug] Parameters::resetToDefaults()" << endl;
      }
      verbosity = DEFAULT_verbosity;
    } // resetToDefaults()

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
    ) const
    {
      os << "[Parameters]" << endl;
      os << "debug = " << debug << endl;
      os << "verbosity = " << verbosity << endl;
    } // writeParameters( basic_ostream & ) const

  }; // End class Parameters

  template <class ParametersType>
  class ParametersModifierTemplate {
    // Boost serialization
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize ( Archive & ar, const unsigned int /* file_version */ )
    {
      // Serialize the isModified_ stuff
      ar & BOOST_SERIALIZATION_NVP( isModified_debug );
      ar & BOOST_SERIALIZATION_NVP( isModified_verbosity );
      // Serialize the parameters
      ar & BOOST_SERIALIZATION_NVP( parameters );
    } // serialize( Archive &, const unsigned int )

  public:

    ParametersType parameters;
    bool isModified_debug;
    bool isModified_verbosity;
  
    ParametersModifierTemplate () :
      parameters()
    {
      if( parameters.debug ) {
        cout << "[debug] ParametersModifierTemplate::<init>()" << endl;
      }
      isModified_reset();
    } // <init>()
  
    // Copy constructor
    template <class AnyParametersModifierTemplate>
    ParametersModifierTemplate ( const AnyParametersModifierTemplate & copy_from ) :
      parameters()
    {
      if( copy_from.parameters.debug >= DEBUG_All ) {
        cout << "[debug] ParametersModifierTemplate::<init>( copy_from )" << endl;
      }
      copyFromNonVirtual( copy_from );
    } // <init>( AnyParametersModifierTemplate const & )
  
    // Copy constructor/operator
    template <class AnyParametersModifierTemplate>
    ParametersModifierTemplate & operator= (
      const AnyParametersModifierTemplate & copy_from
    )
    {
      if( copy_from.parameters.debug >= DEBUG_All ) {
        cout << "[debug] ParametersModifierTemplate::operator=( copy_from )" << endl;
      }
      copyFromNonVirtual( copy_from );
      return *this;
    } // operator=( AnyParametersModifierTemplate const & )
  
    template <class AnyParametersModifierTemplate>
    void
    copyFromNonVirtual ( const AnyParametersModifierTemplate & copy_from )
    {
      if( copy_from.parameters.debug >= DEBUG_All ) {
        cout << "[debug] ParametersModifierTemplate::copyFromNonVirtual( copy_from )" << endl;
      }
      isModified_copyFromNonVirtual( copy_from );

      parameters.copyFromNonVirtual( copy_from.parameters );
    } // copyFromNonVirtual( AnyParametersModifierTemplate const & )

    template <class AnyParametersModifierTemplate>
    void
    isModified_copyFromNonVirtual ( const AnyParametersModifierTemplate & copy_from )
    {
      isModified_debug = copy_from.isModified_debug;
      isModified_verbosity = copy_from.isModified_verbosity;
    } // isModified_copyFromNonVirtual( AnyParametersModifierTemplate const & )

    void reset ()
    {
      isModified_reset();
      parameters.resetToDefaults();
    } // reset()

    void
    isModified_reset ()
    {
      if( parameters.debug ) {
        cout << "[debug] ParametersModifierTemplate::isModified_reset()" << endl;
      }
      isModified_debug = false;
      isModified_verbosity = false;
    } // isModified_reset()

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
    ) const
    {
      os << "[ParametersModifierTemplate]" << endl;
      if( isModified_debug ) {
        os << "debug = " << parameters.debug << endl;
      }
      if( isModified_verbosity ) {
        os << "verbosity = " << parameters.verbosity << endl;
      }
    } // writeParametersModifier ( basic_ostream & ) const

    template<class AnyParameters>
    void
    applyModifications ( AnyParameters & target_parameters ) const
    {
      if( isModified_debug ) {
        target_parameters.debug = parameters.debug;
      }
      if( isModified_verbosity ) {
        target_parameters.verbosity = parameters.verbosity;
      }
    } // applyModifications( AnyParameters & ) const

  }; // End class ParametersModifierTemplate

  typedef ParametersModifierTemplate<galosh::Parameters> ParametersModifier;

} // End namespace galosh

#endif // __GALOSH_PARAMETERS_HPP__
