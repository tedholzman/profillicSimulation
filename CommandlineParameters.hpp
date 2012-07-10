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

namespace po = boost::program_options;

#ifndef DEFAULT_OPTIONS_DESCRIPTION
#define DEFAULT_OPTIONS_DESCRIPTION desc
#endif

#ifndef DEFAULT_VARIABLES_MAP
#define DEFAULT_VARIABLES_MAP vm
#endif

#define GALOSH_DEF_OPT6(DESC,VM,NAME,TYPE,DEFAULTVAL,HELP) \
		  DESC.add_options()(#NAME,po::value<TYPE>()->default_value(DEFAULTVAL),HELP)

#define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP)          \
		GALOSH_DEF_OPT6(DEFAULT_OPTIONS_DESCRIPTION,DEFAULT_VARIABLES_MAP,NAME,TYPE,DEFAULTVAL,HELP)

#endif /* COMMANDLINEPARAMETERS_HPP_ */
