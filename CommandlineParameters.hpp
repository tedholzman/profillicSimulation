/**
 * @file CommandlineParameters
 * @author  Ted Holzman <tholzman@scharp.org>
 * @version 0.01
 * @date 7/2012
 * @section LICENSE
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * \brief macros for describing and accessing commandline parameters
 * \dependencies Depends on BOOST program_options library
 */

#ifndef COMMANDLINEPARAMETERS_HPP_
#define COMMANDLINEPARAMETERS_HPP_

#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>

#include <boost/preprocessor/facilities/empty.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/map.hpp>

/// myVector is a trivial wrapper for std::vector.  It is convenient to use this instead of
/// standard vectors for command line options
template <typename T>
class myVector : public std::vector<T>
{
   public:
	  typedef std::vector<T> parent_t;
	  friend class boost::serialization::access;
	  myVector () {};
      myVector<T> (int n, const T& value)
      : std::vector<T>(n, value)
      {
      }
   private:
      template<class Archive>
      void serialize (Archive & ar, const unsigned int version)
      {
          ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( parent_t );
      }
};

///This function is an overload of the boost::program_options::validate
///The intention is to allow multiple values on a line in the config file
template <typename T>
void
validate (boost::any& v,
          const std::vector<std::string>& values,
          myVector<T>*,
          int)
{
    using namespace boost::program_options;
    myVector<T> tvalues;
    // Make sure no previous assignment to 'a' was made.
    validators::check_first_occurrence(v);
    for(vector<string>::const_iterator it = values.begin(); it!=values.end(); ++it) {
    	stringstream ss(*it);
    	copy(istream_iterator<T>(ss),istream_iterator<T>(),back_inserter(tvalues));
    }
#ifdef DEBUG
    cerr << "tvalues vector is " << tvalues.size() << " long, and contains:" << endl;
    for(int i = 0; i<tvalues.size(); i++) {cerr << tvalues[i] << endl;}
#endif
    v = tvalues;

} // validate(...myVector...)

/// Here we are trying to overload some routines in boost::serialize so it will
/// serialize program_options data, in particular the variables_map
namespace boost {
   namespace serialization {
      // The variables_map is basically a map of variable_value.  It's convenient
      // to be able to serialize both.
      template<class Archive>
      inline void save_construct_data (
         Archive & ar,
         const boost::program_options::variable_value* t,
         const unsigned int file_version
      ) {
         // save data required to construct instance
         ar << t->defaulted();
         ar << t->value();
      }

      template<class Archive>
      inline void load_construct_data (
         Archive & ar,
         boost::program_options::variable_value* t,
         const unsigned int file_version
      ) {
         // retrieve data from archive required to construct new instance
         bool defaulted;
         boost::any value;
         ar >> defaulted;
         ar >> value;
         // TAH - Whacky syntax (borrowed code) new(*anyclass)anyclass-constructor ==
         // new(sizeof(anyclass),*anyclass) -- argh!
         ::new(t)boost::program_options::variable_value(value, defaulted);
      }

      template<class Archive>
      inline void save (
         Archive & ar,
         const boost::program_options::variable_value &t,
         const unsigned int /* file_version */
      ) {
         ar & t.value();
         ar & t.defaulted();
      }

      template<class Archive>
      inline void load (
         Archive & ar,
         boost::program_options::variable_value &t,
         const unsigned int /* file_version */
      ) {
         boost::any v;
         ar & v;
         bool defaulted;
         ar & defaulted;
         boost::program_options::variable_value vv(v, defaulted);
         memcpy(&t, &vv, sizeof(vv));
      }

      /**
       * serialize for variable_value
       */
      template<class Archive>
      inline void serialize (
         Archive & ar,
         boost::program_options::variable_value &t,
         const unsigned int file_version
      ) {
         // we should have nothing to do, as it's all in the constructor, but
         // boost serialization of pairs uses the default constructor and not
         // the above stuff.
         boost::serialization::split_free(ar, t, file_version);
      }

      /** serialize for variables_map
       * just use the superclass's serialize function
       */
      template<class Archive>
      inline void serialize (
         Archive & ar,
         boost::program_options::variables_map &t,
         const unsigned int file_version)
      {
         serialize(
            ar,
            (std::map<std::string,boost::program_options::variable_value> &)t,
            file_version
         );
      }
   }  // serialization
}  // namespace boost

#endif /* COMMANDLINEPARAMETERS_HPP_ */
