/**
 * @file NewSeqanTests.cpp
 * @author  Ted Holzman <tholzman@scharp.org>
 *  modified from the SeqanTest.cpp file by Paul Edlefsen <pedlefsen@gmail.com>
 * @version 0.01
 * @date 11/18/2011
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
 * \brief Unit tests for galosh prolific functions
 * @section DESCRIPTION
 * This program is a test harness for many of the features of the Galosh::prolific
 * (Copyright &copy; 2011 by Paul T Edlefsen, Fred Hutchinson Cancer Research Center.)
 * a suite of programs for working with Profile HMMs.  It is a modification of the
 * original SeqanTest.cpp program.  This version uses the Boost test harness to
 * create and report on individual tests.  It also uses the Boost options parser
 * to select individual units and various test parameters.  Finally, we're attempting
 * to keep these comments Doxygen compliant, and generate documentation from them.
 *
 * To the greatest extent possible, this code is following
 * the code in Paul's SeqanTest.cpp.  However, instead of having internal
 * consts which turn test elements on and off, this code allows command
 * line parameters which select tests or groups of tests. <p> Test "suites"
 * may be run as a group.  Test "cases" may also be run individually. For
 * example: <p><pre>
 *   newseqantest --run_test=test_non_dp 
 *      or
 *   newseqantest --run_test=test_alphabets</pre>
 *
 * To maintain modularity, we're doing some atypical (i.e. different
 * from the simplest examples in the BOOST documentation) here.  For one thing,
 * both the \a TEST module and the \a PROGRAM_OPTIONS module are parsing the
 * command line.  This way we can use command line options to modify the behavior
 * of the unit tests.  The \a TEST module doesn't give hooks into its own parser.
 * For the same reason, we've DEFINEd NO_MAIN above:  it allows us to do our own
 * initializations along with the \a TEST initializations.
 *
 */
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
#include <boost/algorithm/string/case_conv.hpp>

#ifdef __HAVE_MUSCLE
int g_argc;
char **g_argv;
#endif // __HAVE_MUSCLE
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
}

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
	inline std::string operator()(StartStateLabel & state)     const {return MultinomialDistribution<StateLabelTransitionTargets<StartStateLabel,     Plan7>::Type, float>().toString();}
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
}

/**
 * \brief The next several lines are required for the test harness
 *
 */
#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE galosh_profuse tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>
using boost::test_tools::output_test_stream;

/**
 * \var g_argc
 * keep a global copy of argc  around for subsequent modules to use
 */
int g_argc;
/**
 * \var g_argv
 * keep a global copy of argv around for subsequent modules to use
 */
char **g_argv = 0;
using namespace galosh;

#include <boost/program_options.hpp>
/**
 * \note boost variants are used to iterate through sets of types, such as
 * HMM states
 */

namespace po = boost::program_options;
po::options_description utest_opts("Additional unit test options");
po::variables_map cmdline_opt_map;
po::basic_parsed_options<char> *parsed;
int main(int argc, char* argv[]) {
	g_argc = argc;
	g_argv = (char **) argv;

	po::parsed_options p = po::command_line_parser(argc, argv).options(
			utest_opts).allow_unregistered().run();
	parsed = &p;
	po::store(p, cmdline_opt_map, true);
	po::notify(cmdline_opt_map);
	// load prototype for user's unit test init function
	extern ::boost::unit_test::test_suite* init_unit_test_suite(int argc,
			char* argv[]);
	boost::unit_test::init_unit_test_func init_func = &init_unit_test_suite;
	return ::boost::unit_test::unit_test_main(init_func, argc, argv);
}
/**
 * \fn std::string getOptVal(const std::string optName)
 * \brief fetches the value of program option optName
 * \param optName
 *    The string value of command line parameter you want to fetch
 * \brief This is a kludge, to be repaired when I fully understand the program_options boost module
 * \todo this is not good code.  It's got to be repaired when there's time to
 * understand the nature of the \a allow_unregistered option of the \a basic_command_line_parser.
 * \note Accesses a global variable, \a parsed.  (Shame!)
 * \return Returns the string value of option optName
 */
std::string getOptVal(const std::string optName) {
	std::string retVal = "";
	vector<po::basic_option<char> >::iterator it;
	for (it = parsed->options.begin(); it < parsed->options.end(); it++) {
		if (it->string_key.compare(optName) == 0)
			return it->value[0];
	}
	return retVal;
}
/**
 * \def FLOAT_POINT_EQUALITY_TOL_PCT
 * When we test if two floating point numbers are equal, they must
 * really fall within FLOAT_POINT_EQUALITY_TOL_PCT \% of each other
 */
#define FLOAT_POINT_EQUALITY_TOL_PCT 0.00000001
/**
 * \brief This is a canonical example of a unit test.
 *
 * This test succeeds if \f$ 1 = 1\f$.
 * It also defines commandline option value.
 */
BOOST_AUTO_TEST_CASE( sanity_check ) {
	std::cout << "Sanity check: should always be TRUE" << std::endl;
	BOOST_CHECK_MESSAGE( 1 == 1, "Test failed: 1 doesn't equal 1");
} //sanity_check

/**
 * \namespace test_non_dp
 * \brief suite marker, really a macro call, marking the beginning of
 * a set of tests.
 *
 * The macro parameter, \a test_non_dp is name of suite, and also defines
 * a commandline option to select the/a test_non_dp test(s)
 */
BOOST_AUTO_TEST_SUITE( test_non_dp )

/**
 * \brief Tests of seqan alphabet conversions
 */
BOOST_AUTO_TEST_CASE( test_alphabets )
{
	std::cout << "test_alphabets: Testing alphabet conversions" << std::endl;
	Dna a = 'a';
	BOOST_CHECK_MESSAGE( /// Test conversion of Dna to char
			(char) a == 'A',
			"Dna to string yielded "<<a<< " instead of A"
	);
	Dna5 b = 'f';//'f' is unknown character
	BOOST_CHECK_MESSAGE(///Test conversion of Dna5 to char (and converting unknowns to N's)
			(char) b == 'N',
			"Dna5 yielded " << b << " instead of N"
	);
	//Many SeqAn alphabet classes can be converted into each other.
	b = a;
	BOOST_CHECK_MESSAGE(///Test conversion of Dna5 to Dna to char
			(char) b == 'A',
			"Conversion from Dna to Dna5 yielded " << b << " instead of A"
	);
	Iupac c = a;
	BOOST_CHECK_MESSAGE(///Test conversion of Iupac to Dna to char
			(char) c == 'A',
			"Conversion from Dna to Iupac yielded " << c << " instead of A"
	);
} // end of test_alphabets
bool ambiguity_tested = false;
/**
 * \brief Tests of ambiguous letters in seqan alphabets
 */
BOOST_AUTO_TEST_CASE( test_ambiguous )
{
	std::cout << "test_ambiguity: Testing inter-alphabet ambiguities" << std::endl;
	ambiguity_tested=true;
	BOOST_CHECK_MESSAGE( //Test if Dna is ambiguous over Dna5
			!isTrue( galosh::IsAmbiguous<seqan::Dna, seqan::Dna5>::Type() ),
			"Dna appears to have ambiguities with respect to Dna5.  It shouldn't."
	);
	BOOST_CHECK_MESSAGE(//Test if Dna5 is ambiguous over Dna
			isTrue( galosh::IsAmbiguous<seqan::Dna5, seqan::Dna>::Type() ),
			"Dna5 doesn't appear to have ambiguities with respect to Dna.  It should"
	);
	Dna dna_residue;
	Dna5 a = 'a';
	BOOST_CHECK_MESSAGE(//Test how many elements in Dna5 a match an element in Dna
			galosh::ambiguousCount( a, Dna() ) == 1,
			"'A' in Dna5 matches " << (galosh::ambiguousCount( a, Dna() ) == 1) << " letters in Dna, should only match 1"
	);
	galosh::ambiguousAssign( dna_residue, a, 0 );
	BOOST_CHECK_MESSAGE(//Find a Dna letter matching a Dna5 letter
			(char) dna_residue == 'A',
			"ambiguousAssign error: Dna5 'A' matches Dna " << dna_residue
	);
	Dna5 n = 'n';
	size_t num_elements = galosh::ambiguousCount( n, Dna() );
	BOOST_CHECK_MESSAGE(//Test how many Dna elements match Dna5's letter n
			num_elements == 4,
			"Dna5's N element matches " << num_elements << " Dna elements, when it should match 4"
	);
	string nucleotides ("ATGC");
	// Test that ambiguousAssign matches all four nucleotides against N
	for( int i = 0; i < num_elements; i++ ) {
		galosh::ambiguousAssign( dna_residue, n, i );
		size_t pos = nucleotides.find((char)dna_residue);
		BOOST_CHECK_MESSAGE(
				pos < nucleotides.length(),
				"ambiguousAssign returned " << (char)dna_residue << " which doesn't appear to be a nucleotide!"
		);
		if(pos<nucleotides.length()) nucleotides.erase(pos,1);
	}
	BOOST_CHECK_MESSAGE(
			nucleotides.length() == 0,
			"ambiguousAssign missed some nucleotides: " + nucleotides
	);
} //end of test_ambiguous

/**
 * \brief Test of basic seqan and galosh sequence manipulation operations
 */
BOOST_AUTO_TEST_CASE( test_sequences )
{
	std::cout << "test_sequences: Testing basic seqan and galosh biosequence operations" << std::endl;
	string anypeptide = "anypeptide";
	seqan::Peptide prot = anypeptide;
	BOOST_CHECK_MESSAGE( ///Test if seqan peptide sequence length is correct
			length(prot) == anypeptide.length(),
			"Length of " << toCString(prot) << " should be " << anypeptide.length() <<
			"but ended up being " << length(prot)
	);
	BOOST_CHECK_MESSAGE(///Test array-like indexing of seqan peptide sequence
			(char)prot[ 9 ] == toupper(anypeptide[9]),
			"The 9th character of " << toCString(prot) << " should be " << anypeptide[9] <<
			"but ended up being " << (char)prot[9]
	);
	prot += "anewend"; anypeptide += "anewend";
	string protStr;
	assign(protStr,prot);
	BOOST_CHECK_MESSAGE(///Test concatenation of string to seqan peptide sequence
			boost::to_upper_copy(anypeptide).compare(protStr) == 0,
			"Test sequence should = " << boost::to_upper_copy(anypeptide) <<
			"but ended up = " << toCString(prot)
	);
	BOOST_CHECK_MESSAGE(///Test length of concatenated seqan peptide sequence
			length(prot) == anypeptide.length(),
			"Test sequence should be " << anypeptide.length() << " characters long, " <<
			"but ended up being " << length(prot) << " characters long"
	);

	/// \todo Figure out what is going on with galosh::Sequence<>.<br>
	///It takes a lot of effort to convert it into something that's
	///comparable to a string (either string or char *).  Possibly
	///insert new overloaded "assign" into galosh::Sequence<>.
	string anyDNA = "acgt";
	galosh::Sequence<Dna> dna_seq = anyDNA;
	/// \internal Trying (basically) to convert a galosh::Sequence into std::string
	/// None of the really obvious ways work (cast, convert, assign).  For some reason that
	/// I don't yet understand, a galosh::Sequence, which is a direct subclass of
	/// of a seqan::String, cannot be cast into the String and cannot be converted
	/// the way a seqan::Peptide can, into a std::string.  Currently I've written
	/// a "toString" routine, which does the job along with the ridiculous dynamic
	/// cast.
	seqan::String<Dna,seqan::Alloc<> > *sp = dynamic_cast<seqan::String<Dna,seqan::Alloc<> > *>(&dna_seq);
	string dna_seq_str = toString(*sp);
	string anyDNA_s = (boost::to_upper_copy(anyDNA));
	BOOST_CHECK_MESSAGE(///Test creation of galosh Sequence from string
			anyDNA_s.compare(dna_seq_str) == 0,
			"Test sequence should = " << boost::to_upper_copy(anyDNA) <<
			" but ended up being = " << *toCString(dna_seq)
	);
	BOOST_CHECK_MESSAGE(///Test length of created galosh::Sequence
			length(dna_seq) == anyDNA.length(),
			"Test sequence should be " << anyDNA.length() << " characters long, " <<
			"but ended up being " << length(dna_seq) << " characters long"
	);
} //end of test_sequences

/**
 *  \def definitions for \a default fasta files.
 */
#define DEFAULT_UNGAPPED_FASTA_FILE "fasta/dna.10.10.fasta"
#define DEFAULT_GAPPED_FASTA_FILE "fasta/21U.fa.muscle.first10"
/**
 * \brief Tests of basic fasta file read and conversions.
 *
 * Also, first (currently mediocre) example of defining and using
 * program options other than those defined in the test interface.  We also use
 * the boost \a output_test_stream to determine if a fasta data structure
 * correctly represents the file it was read from.
 */
BOOST_AUTO_TEST_CASE( test_fasta )
{
	std::string ungapped_fasta_file = DEFAULT_UNGAPPED_FASTA_FILE;
	std::string gapped_fasta_file = DEFAULT_GAPPED_FASTA_FILE;
	/**
	 * \todo Try to use the boost program options API "correctly".
	 */
	std::string tmp = getOptVal("ungapped_fasta");
	if(tmp.compare("") != 0) ungapped_fasta_file = tmp;
	tmp = getOptVal("gapped_fasta");
	if(tmp.compare("") != 0) gapped_fasta_file = tmp;
	utest_opts.add_options()
	("ungapped_fasta",po::value<std::string>(),"Path of ungapped fasta file to use in test-fasta")
	("gapped_fasta",po::value<std::string>(),"Path of gapped fasta file to use in test-fasta");
	po::notify(cmdline_opt_map);

	std::cout << "test_fasta:     Testing ungapped and gapped fasta file reads and conversions" << std::endl;

	//
	//// Fasta with no gaps
	galosh::Fasta<Dna> fasta;
	fasta.fromFile(ungapped_fasta_file);
	{
   	   output_test_stream output( ungapped_fasta_file, true );
	   output << fasta;
	   BOOST_CHECK( output.match_pattern() ); /// Ungapped file correctly scanned?
	}
	///////
	//// Fasta: with gaps.  Use "char" since basic seqan alphabets don't support gap chars.

	galosh::Fasta<char> aligned_fasta;
	aligned_fasta.fromFile(gapped_fasta_file);
	/**
	 * \todo currently using unit test module output_test_stream to check the validity
	 * of the fasta files; not just printing them out.  output_test_stream.match_pattern()
	 * kindof sucks:  the output must match the "pattern file" exactly, character
	 * for character.  In the case of fasta files, it might be better
	 * to use a checksum.  In any case, it would be nice if pattern_match used regular
	 * expressions or at least wildcards.
	 */
	{
	   output_test_stream output( gapped_fasta_file, true );
	   output << aligned_fasta;
	   BOOST_CHECK( output.match_pattern() ); /// Gapped file correctly scanned?
	}

	///
	/// \todo  I'm thinking of just scrapping Fasta altogether, using
	/// StringSet instead, with _loadSequences from
	/// seqan/graph_utils/utility_parsing.h, as in seqan_tcoffee.
	/// Eventually.

} //end of test_fasta
/**
 *
 */
BOOST_AUTO_TEST_CASE( test_multinomials )
{
	/// \todo note: ambiguous test \a within multinomials <br>
	///can that be done within this model?  Can we tell if a test (e.g.
	///ambiguous) is going to be executed? Currently using kludge.

	std::cout << "test_multinomials: Testing multinomial distribution generation and manipulation for DNA" << std::endl;

	Dna c_2( 1 );
	BOOST_CHECK_MESSAGE(/// Test conversion of Dna (created from int) to char
			(char) c_2 == 'C',
			"Dna to string yielded "<< c_2 << " instead of C"
	);

	galosh::MultinomialDistribution<Dna, realspace> dna_dist;
	dna_dist[ c_2 ] = .4;
	/**
	 * \note This value is illegal:  doesn't sum to 1
	 */

	BOOST_CHECK_MESSAGE( /// Test modification of a Multinomial distribution object
			dna_dist.m_elementCount == 4 && dna_dist.m_probs[0] == 0.25 &&
			dna_dist.m_probs[1].prob() == 0.4 &&
			dna_dist.m_probs[2].prob() == 0.25 && dna_dist.m_probs[3].prob() == 0.25,
			"dist should be A=T=G=0.25, C=0.4; but instead it is " << dna_dist
	);
    /**
     * \todo it would be nice if all the constants in this module (e.g. 12,0.25)
     * were \#DEFINEd and some reasonable subset of them modifiable via command
     * line options
     */
	galosh::MultinomialDistribution<galosh::StateLabel, float> state_dist;
	bool elValOK = true;
    int n_elements = (sizeof state_dist.m_probs/sizeof state_dist.m_probs[0]);
    float evenDist = 1.0F/(float)n_elements;
	for(int i = 0; i<n_elements;i++) {
	   elValOK = state_dist.m_probs[i] == evenDist;
	   if(!elValOK) break;
	}
	BOOST_CHECK_MESSAGE( /// Test structure of a StateLabel distribution
	   n_elements == 12 && elValOK,
	   "Bad structure for StateLabel distribution " << state_dist
	);

	if( ambiguity_tested )
	{
		std::cout << "test_multinomials:       Testing multinomial distribution with ambiguities (for DNA)" << std::endl;
		Dna5 a = 'a';
		Dna5 n = 'n';
		BOOST_CHECK_MESSAGE( /// Check initial distributions of single Dna5 ambiguous and unambiguous elements
	       dna_dist.ambiguousSum( a ) == 0.25 && dna_dist.ambiguousSum( n ) == 1.15,
           "Initial 'A' distribution should be 0.25 and it's " << dna_dist.ambiguousSum( a ) <<
           "initial 'N' distribution should be 1.15 and it's " << dna_dist[ n ]
		);
        /// \internal
		/// Uncomment this, and it should fail to compile, since AminoAcid is
		/// not ambiguous over Dna.
		/// AminoAcid aminoacid_n = 'n';
		/// std::cout << "The dna_dist, accessed using AminoAcid 'n', returns " << dna_dist[ aminoacid_n ] << endl;
		output_test_stream output;
		dna_dist.ambiguousIncrement( n, 1.0 );
		output << dna_dist;
		BOOST_CHECK_MESSAGE( /// Test ambiguousIncrement of 1.0 to dna distribution
		    output.is_equal("(A=0.5,C=0.65,G=0.5,T=0.5)"),
		    "Distribution should equal (A=0.5,C=0.65,G=0.5,T=0.5), but instead equals " << dna_dist
		);
		galosh::MultinomialDistribution<Dna, realspace>::AmbiguousValue<Dna5> ambiguous_value = dna_dist[ n ];
		ambiguous_value += 1.0;
		output << dna_dist;
		BOOST_CHECK_MESSAGE( /// Test autoincrement (+=)  of ambiguous DNA value (n), effect on distribution
		   output.is_equal("(A=0.75,C=0.9,G=0.75,T=0.75)"),
		   "Distribution should equal (A=0.75,C=0.9,G=0.75,T=0.75), but instead equals " << dna_dist
		);
        BOOST_CHECK_CLOSE(ambiguous_value.prob(),3.15,FLOAT_POINT_EQUALITY_TOL_PCT);  /// Test probability of autoincremented ambiguous element
		ambiguous_value = 1.0;
        output << dna_dist;
		BOOST_CHECK_MESSAGE( /// Test direct assignment of ambiguous DNA value (n), effect on distribution
		   output.is_equal("(A=0.25,C=0.25,G=0.25,T=0.25)"),
		   "Distribution should equal (A=0.25,C=0.25,G=0.25,T=0.25), but instead equals " << dna_dist
		);
        BOOST_CHECK_CLOSE(ambiguous_value.prob(),1.00,FLOAT_POINT_EQUALITY_TOL_PCT);  /// Test probability of autoincremented ambiguous element
		dna_dist[ n ] += 1.0;
		output << dna_dist;
		BOOST_CHECK_MESSAGE( /// Test direct autoincrement,  distrib[n] += 1
		   output.is_equal("(A=0.5,C=0.5,G=0.5,T=0.5)"),
		   "Distribution should equal (A=0.5,C=0.5,G=0.5,T=0.5), but instead equals " << dna_dist
		);
	} // End if test_ambiguous (and test_multinomials)

	/// \todo Use seqan::FrequencyDistribution?
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

} //end of test_multinomials


/**
 * \brief Test properties of profile HMM states.
 */
BOOST_AUTO_TEST_CASE( test_profile_hmm_states )
{
	using namespace galosh;
	/**
	 *  \var stateLabels
	 *  \brief a useful map that contains the full text name of the state label,
	 *  as the key, and a variant containing a reference to the appropriate state
	 *  label object.
	 *  \see state_VALUE_generic, state_SIMPLE_generic
	 */
	std::map<std::string,AnyStateLabel> stateLabels;
	/**
	 * \def ADD_MAP(x)
	 * \brief Tiny macro to insert values into stateLabels.
     *
	 * \todo Move this, stateLabels, AnyStateLabel typedef, and accessor classes
	 * elsewhere
	 */
#define ADD_MAP(x) stateLabels.insert(pair<std::string,AnyStateLabel>(#x,x()))
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

	/**
	 * \var curCorrect
	 * store correct info about each state.
	 * \todo This concept is a little fragile.  When the algebra or methodology changes, the "correct answers" will
	 * change.  Re-think this later on.
	 */
	std::map<std::string,std::string> curCorrect;
#define ADD_MAP_CURCORRECT(x,y) curCorrect.insert(pair<std::string,std::string>(x,y))
    ADD_MAP_CURCORRECT("StartStateLabel","label: 0; code: S; simple; not emitting; not associated; Plan7 dist=(N=1)");
	ADD_MAP_CURCORRECT("PreAlignStateLabel","label: 1; code: N; simple; emitting; not associated; Plan7 dist=(N=0.5,B=0.5)");
	ADD_MAP_CURCORRECT("BeginStateLabel","label: 2; code: B; simple; not emitting; not associated; Plan7 dist=(M=0.5,D=0.5)");
	ADD_MAP_CURCORRECT("MatchStateLabel","label: 3; code: M; simple; emitting; associated; Plan7 dist=(M=0.333333,I=0.333333,D=0.333333)");
	ADD_MAP_CURCORRECT("InsertionStateLabel","label: 4; code: I; simple; emitting; not associated; Plan7 dist=(M=0.5,I=0.5); Plan9 dist=(M=0.333333,I=0.333333,D=0.333333)");
	ADD_MAP_CURCORRECT("DeletionStateLabel","label: 5; code: D; simple; not emitting; associated; Plan7 dist=(M=0.5,D=0.5); Plan9 dist=(M=0.333333,I=0.333333,D=0.333333)");
	ADD_MAP_CURCORRECT("EndStateLabel","label: 6; code: E; simple; not emitting; not associated; Plan7 dist=(C=0.5,J=0.5)");
	ADD_MAP_CURCORRECT("LoopStateLabel","label: 7; code: J; simple; emitting; not associated; Plan7 dist=(J=0.5,B=0.5)");
	ADD_MAP_CURCORRECT("PostAlignStateLabel","label: 8; code: C; simple; emitting; not associated; Plan7 dist=(C=0.5,T=0.5)");
	ADD_MAP_CURCORRECT("TerminalStateLabel","label: 9; code: T; simple; not emitting; not associated");

    /**
     * Iterate through all states, check each one against a known model.
     */
	std::cout << "test_profile_hmm_states: Testing properties and access methods of HMM states" << std::endl;
    for(std::map<std::string,AnyStateLabel>::iterator asl=stateLabels.begin(); asl!=stateLabels.end(); asl++)
       BOOST_CHECK_MESSAGE(
    	      ((std::string)curCorrect(asl->first)).compare(stateInfo(asl->second)) == 0,
    	      asl->first << "should be " << curCorrect(asl->first) << "and instead it's" <<
    	         stateInfo(asl->second)
    	   );

    	std::cout << asl->first << ": " << stateInfo(asl->second) << std::endl;

  } // End test_profile_hmm_states

BOOST_AUTO_TEST_SUITE_END()//end of test_non_db
