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
#include <boost/preprocessor/facilities/empty.hpp>

namespace po = boost::program_options;

#ifndef DEFAULT_OPTIONS_DESCRIPTION
#define DEFAULT_OPTIONS_DESCRIPTION desc
#endif

#ifndef DEFAULT_VARIABLES_MAP
#define DEFAULT_VARIABLES_MAP vm
#endif

#ifndef TMP_EXTRA_STUFF
#define TMP_EXTRA_STUFF BOOST_PP_EMPTY()
#endif

#define GALOSH_DEF_OPT6(DESC,VM,NAME,TYPE,DEFAULTVAL,HELP) \
		  DESC.add_options()(#NAME,po::value<TYPE>()->default_value(DEFAULTVAL) TMP_EXTRA_STUFF,HELP)
/// Note we have to duplicate code here instead of having OPT call OPT6 because on occasion
/// OPT will be given a macro argument which expands to contain commas.  When this happens,
/// the inner macro call will have the wrong number of arguments.  See @ProfuseTestOptions.hpp
/// for the \a profileLengths entry.
#define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP)          \
		  DEFAULT_OPTIONS_DESCRIPTION.add_options()(#NAME,po::value<TYPE>()->default_value(DEFAULTVAL) TMP_EXTRA_STUFF,HELP)

/// myVector is a trivial wrapper for std::vector.  It is convenient to use this instead of
/// standard vectors for command line options
template <typename T>
class myVector : public std::vector<T>
{
   public:
	  myVector () {};

      myVector<T> (int n, const T& value)
      : std::vector<T>(n, value)
      {
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

    cerr << "tvalues vector is " << tvalues.size() << " long, and contains:" << endl;
    for(int i = 0; i<tvalues.size(); i++) {cerr << tvalues[i] << endl;}
    v = tvalues;
}

#endif /* COMMANDLINEPARAMETERS_HPP_ */
