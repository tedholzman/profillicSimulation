/*
 * NewSeqanTests.hpp
 *
 *  Created on: Dec 28, 2011
 *      Author: tedholzman
 */

#ifndef NEWSEQANTESTS_HPP_
#define NEWSEQANTESTS_HPP_

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
///\#include "ProfileGibbs.hpp" // \todo Add tests for gibbs...
#include "ProfuseTest.hpp"
#include "AminoAcid20.hpp"

#include <iostream>
#include <string>
#include </usr/include/stdio.h>
#include </usr/include/stdlib.h>
#include </usr/include/string.h>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/find_motif.h>

#ifdef __HAVE_MUSCLE
#include "muscle/distfunc.h"
#include "muscle/clustsetdf.h"
#include "muscle/clust.h"
#include "muscle/tree.h"
#endif // __HAVE_MUSCLE
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#ifdef __HAVE_MUSCLE
int g_argc;
char **g_argv;
#endif // __HAVE_MUSCLE

/**
 * \def FLOAT_POINT_EQUALITY_TOL_PCT
 * When we test if two floating point numbers are equal, they must
 * really fall within FLOAT_POINT_EQUALITY_TOL_PCT \% of each other
 */
#define FLOAT_POINT_EQUALITY_TOL_PCT 0.00000001

using namespace seqan;

/**
 * \brief This function converts a Sequence, or a String to a std::string
 *
 * This is basically a hack, and shouldn't be necessary.  Probably a better solution
 * would be to invent the correct "assign" routine.  If we keep it, it should
 * probably go into the Galosh.hpp file.  --TAH
 *
 * \param source  a reference to a seqan::String
 * \return a std::string, translated from the compressed bits of the seqan::String in the appropriate alphabet
 *
 */
template<typename TSource>
inline std::string toString(TSource & source) {
	typename seqan::Iterator<TSource, seqan::Standard>::Type it = seqan::begin(
			source, seqan::Standard());
	typename seqan::Iterator<TSource, seqan::Standard>::Type it_end =
			seqan::end(source, seqan::Standard());
	std::string retVal = "";
	for (; it < it_end; ++it) {
		typename GetValue<TSource const>::Type val_ = getValue(it);
		retVal += val_;
	}
	return retVal;
} //end of toString()

/**
 * \fn static bool isTrue(T const & tag)
 * \brief Function for determining seqan True values
 * \param tag
 *    reference to any type
 * \return
 *    boolean false if type refers to anything but a seqan True object
 */
template<typename T>
static bool isTrue(T const & tag) {
	return false;
}//end of isTrue(), all types except seqan::True

/**
 * \fn static bool isTrue(seqan::True const & tag)
 * \brief Function for determining seqan True values
 * \param tag
 *    reference to any type
 * \return
 *    boolean true if type refers to seqan True object
 */
static bool isTrue(seqan::True const & tag) {
	return true;
}
/**
 * \typedef AnyStateLabel
 * \brief a boost variant that allows iterating through HMM states
 * \see state_VALUE_generic,state_SIMPLE_generic
 */

/**
 * \note using boost \a variants to iterate through sets of types, such as
 * HMM states
 */
#include <boost/variant.hpp>
typedef
boost::variant<
	galosh::StartStateLabel,
	galosh::PreAlignStateLabel,
	galosh::BeginStateLabel,
	galosh::MatchStateLabel,
	galosh::InsertionStateLabel,
	galosh::DeletionStateLabel,
	galosh::EndStateLabel,
	galosh::LoopStateLabel,
	galosh::PostAlignStateLabel,
	galosh::TerminalStateLabel
> AnyStateLabel;

/**
 * \class state_VALUE_generic
 * \brief utility for iterating through state types
 *
 * This is a \a visitor subclass of the \a boost::variant system.  We're using it to
 * iterate through states.  When used with \a boost::apply_visitor it is basically a
 * get-type generic accessor for the state label numeric VALUE field.  This treats
 * objects of the different classes subsumed within the AnyClassLabel type '
 * (a boost::variant) as if they were objects of the same class.
 *
 * \todo The \ref AnyStateLabel typedef and these accessor classes should probably
 * be moved to the Galosh.hpp (or someplace similar).  They are probably useful
 * in other programs.
 */
class
state_VALUE_generic: public boost::static_visitor<int>
{
public:
   template <typename T>
   int operator()( T & state ) const
    {
	   return galosh::StateLabelId<T>::VALUE;
    }
};
/**
 * \class state_SIMPLE_generic
 * \brief generic accessor for isSimple "field" of state labels
 */
class
state_SIMPLE_generic: public boost::static_visitor<bool>
{
public:
	template <typename T>
	bool operator()( T & state) const
	{
		return isTrue(typename IsSimple<T>::Type() );
	}
};
/**
 * \class state_EMITTING_generic
 * \brief generic accessor for isSimple "field" of state labels
 */
class
state_EMITTING_generic: public boost::static_visitor<bool>
{
public:
	template <typename T>
	bool operator()( T & state) const
	{
		return isTrue(typename galosh::IsEmitting<T>::Type() );
	}
};
/**
 * \class state_ASSOCIATED_generic
 * \brief generic accessor for IsAssociatedWithPosition "field" of state labels
 */
class
state_ASSOCIATED_generic: public boost::static_visitor<bool>
{
public:
	template <typename T>
	bool operator()( T & state) const
	{
		return isTrue(typename galosh::IsAssociatedWithPosition<T>::Type() );
	}
};
/**
 * \class state_DIST7
 * \brief generic accessor for plan9 multinomial distribution associated with state labels
 *
 * It would be marginally better to have this class return a distribution object.
 * Unfortunately the distribution objects are all different classes, and that means
 * this would have to return a \a boost::variant.  I am not certain that will
 * work.  So for the time being I'm taking the easy way out.
 *
 */
using namespace galosh;
#define UNDEFINED_DIST_STRING ""
class
state_DIST7: public boost::static_visitor<std::string>
{
public:
	inline std::string operator()(TerminalStateLabel & state)  const {return UNDEFINED_DIST_STRING;}
	inline std::string operator()(StartStateLabel & state)     const {return MultinomialDistribution< StateLabelTransitionTargets< StartStateLabel,     Plan7 >::Type, float >().toString();}
	inline std::string operator()(PreAlignStateLabel & state)  const {return MultinomialDistribution<StateLabelTransitionTargets<PreAlignStateLabel,  Plan7>::Type, float>().toString();}
	inline std::string operator()(BeginStateLabel & state)     const {return MultinomialDistribution<StateLabelTransitionTargets<BeginStateLabel,     Plan7>::Type, float>().toString();}
	inline std::string operator()(MatchStateLabel & state)     const {return MultinomialDistribution<StateLabelTransitionTargets<MatchStateLabel,     Plan7>::Type, float>().toString();}
	inline std::string operator()(InsertionStateLabel & state) const {return MultinomialDistribution<StateLabelTransitionTargets<InsertionStateLabel, Plan7>::Type, float>().toString();}
	inline std::string operator()(DeletionStateLabel & state)  const {return MultinomialDistribution<StateLabelTransitionTargets<DeletionStateLabel,  Plan7>::Type, float>().toString();}
	inline std::string operator()(EndStateLabel & state)       const {return MultinomialDistribution<StateLabelTransitionTargets<EndStateLabel,       Plan7>::Type, float>().toString();}
	inline std::string operator()(LoopStateLabel & state)      const {return MultinomialDistribution<StateLabelTransitionTargets<LoopStateLabel,      Plan7>::Type, float>().toString();}
	inline std::string operator()(PostAlignStateLabel & state) const {return MultinomialDistribution<StateLabelTransitionTargets<PostAlignStateLabel, Plan7>::Type, float>().toString();}
}; // end of state_DIST7

class
state_DIST9: public boost::static_visitor<std::string>
{
public:
	inline std::string operator()(TerminalStateLabel & state)  const {return UNDEFINED_DIST_STRING;}
	inline std::string operator()(StartStateLabel & state)     const {return UNDEFINED_DIST_STRING;}
	inline std::string operator()(PreAlignStateLabel & state)  const {return UNDEFINED_DIST_STRING;}
	inline std::string operator()(BeginStateLabel & state)     const {return UNDEFINED_DIST_STRING;}
	inline std::string operator()(MatchStateLabel & state)     const {return UNDEFINED_DIST_STRING;}
	inline std::string operator()(InsertionStateLabel & state) const {return MultinomialDistribution<StateLabelTransitionTargets<InsertionStateLabel, Plan9>::Type, float>().toString();}
	inline std::string operator()(DeletionStateLabel & state)  const {return MultinomialDistribution<StateLabelTransitionTargets<DeletionStateLabel,  Plan9>::Type, float>().toString();}
	inline std::string operator()(EndStateLabel & state)       const {return UNDEFINED_DIST_STRING;}
	inline std::string operator()(LoopStateLabel & state)      const {return UNDEFINED_DIST_STRING;}
	inline std::string operator()(PostAlignStateLabel & state) const {return UNDEFINED_DIST_STRING;}
}; // end of state_DIST9


/**
 * \fn std::string numberToString( T Number )
 * \brief Tiny helper function to turn any number into a string.
 */
template <typename T>
inline
std::string numberToString ( T Number )
{
	stringstream ss;
	ss << Number;
	return ss.str();
}

/**
 * \fn std::string stateInfo(AnyStateLabel & sl)
 * \param sl  State label about which to return summary info
 * \return String with state label information
 * \brief This is a helper function for the test_profile_hmm_states unit test
 *
 */
std::string stateInfo(AnyStateLabel & sl)
{
   std::string retVal = "";
   int slNumeric = apply_visitor(state_VALUE_generic(),sl);
   retVal += "label: ";
   retVal += numberToString(slNumeric);
   char code = galosh::StateLabel(slNumeric);
   retVal += "; code: ";
   retVal += code;
   bool simple = apply_visitor(state_SIMPLE_generic(),sl);
   retVal += "; ";
   if(!simple) retVal += "not ";
   retVal += "simple";
   bool emitting = apply_visitor(state_EMITTING_generic(),sl);
   retVal += "; ";
   if(!emitting) retVal += "not ";
   retVal += "emitting";
   bool associated = apply_visitor(state_ASSOCIATED_generic(),sl);
   retVal += "; ";
   if(!associated) retVal += "not ";
   retVal += "associated";
   std::string dist7 = apply_visitor(state_DIST7(),sl);
   if(dist7.compare("") != 0) retVal += "; Plan7 dist=" + dist7;
   std::string dist9 = apply_visitor(state_DIST9(),sl);
   if(dist9.compare("") != 0) retVal += "; Plan9 dist=" + dist9;
   return retVal;
} // of stateInfo

std::map<std::string,AnyStateLabel> stateLabels;
/**
 *  \var stateLabels
 *  \brief a useful map that contains the full text name of the state label,
 *  as the key, and a variant containing a reference to the appropriate state
 *  label object.
 *  \see state_VALUE_generic, state_SIMPLE_generic
 */
/**
 * \def ADD_MAP(x)
 * \brief Tiny macro to insert values into stateLabels.
 */
#define ADD_MAP(x) stateLabels.insert(pair<std::string,AnyStateLabel>(#x,x()))
/**
 * \fn void initialize_globals()
 * \brief initialize complex data structures that are only used for reference
 *
 * Designed for maps, iterators, program options, etc. -- which would clog up
 * the code in other places.
 */
void initialize_globals() {
   ADD_MAP(StartStateLabel);
   ADD_MAP(PreAlignStateLabel);
   ADD_MAP(BeginStateLabel);
   ADD_MAP(MatchStateLabel);
   ADD_MAP(InsertionStateLabel);
   ADD_MAP(DeletionStateLabel);
   ADD_MAP(EndStateLabel);
   ADD_MAP(LoopStateLabel);
   ADD_MAP(PostAlignStateLabel);
   ADD_MAP(TerminalStateLabel);
} // of initialize_globals

int Max_Unit_Name_Width;

#endif /* NEWSEQANTESTS_HPP_ */
