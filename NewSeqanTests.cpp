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
//#include "ProfileGibbs.hpp" // \todo Add tests for gibbs...
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
 *
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
 *
 * \brief Function for determining seqan True values
 * \param tag
 *    reference to any type
 * \return 
 *    boolean true if type refers to seqan True object
 */
static bool isTrue(seqan::True const & tag) {
	return true;
}
#include <boost/variant.hpp>
typedef boost::variant<
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
 * iterate through states.
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

class
state_SIMPLE_generic: public boost::static_visitor<bool>
{
public:
	template <typename T>
	bool operator()( T & state) const
	{
		return isTrue(typename IsSimple<T>::Type());
	}
};

template <typename T>
string numberToString ( T Number )
{
	stringstream ss;
	ss << Number;
	return ss.str();
}

/**
 * \fn std::string stateInfo(AnyStateLabel sl)
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

   //std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::StartStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsEmitting<galosh::StartStateLabel>::Type() ) ? "" : "not " ) << "emitting." << endl;
   //std::cout << galosh::StateLabel( (int)galosh::StateLabelId<galosh::StartStateLabel>::VALUE ) << " is " << ( isTrue( galosh::IsAssociatedWithPosition<galosh::StartStateLabel>::Type() ) ? "" : "not " ) << "associated with an ancestral sequence position." << endl;
   //galosh::MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::StartStateLabel, galosh::Plan7>::Type, float> Start_state_dist;
   //std::cout << Start_state_dist << std::endl;

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
 * \section DESCRIPTION
 * To maintain modularity, we're doing some atypical (i.e. different
 * from the simplest examples in the BOOST documentation) here.  For one thing,
 * both the \a TEST module and the \a PROGRAM_OPTIONS module are parsing the
 * command line.  This way we can use command line options to modify the behavior
 * of the unit tests.  The \a TEST module doesn't give hooks into its own parser.
 * For the same reason, we've DEFINEd NO_MAIN above:  it allows us to do our own
 * initializations along with the \a TEST initializations.
 */

/**
 * \var keep a global copy of argc and argv around for subsequent modules to use
 */
int g_argc;
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
 * \brief this is a kludge, to be repaired when I fully understand the program_options boost module
 * \todo this is not good code.  It's got to be repaired when there's time to
 * understand the nature of the \a allow_unregistered option of the \a basic_command_line_parser
 * \param optName string value of the option name
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
 * \def When we test if two floating point numbers are equal, they must
 * really fall within FLOAT_POINT_EQUALITY_TOL_PCT % of each other
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
 * \fn void BOOST_AUTO_TEST_SUITE( string test_non_dp )
 * \brief suite marker, really a macro call
 * \param test_non_dp name of suite, also defines commandline option
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
     * were #DEFINEd and some reasonable subset of them modifiable via command
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
		BOOST_CHECK_MESSAGE( /// Test ambiguousIncrement of 1.0 to dna dist
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
	std::map<std::string,AnyStateLabel> stateLabels;
#define ADD_MAP(x) stateLabels.insert(pair<std::string,AnyStateLabel>(#x,x()))
    ADD_MAP(StartStateLabel);
	ADD_MAP(PreAlignStateLabel);
	ADD_MAP(BeginStateLabel);
	ADD_MAP(MatchStateLabel);
	ADD_MAP(InsertionStateLabel);
	ADD_MAP(EndStateLabel);
	ADD_MAP(LoopStateLabel);
	ADD_MAP(PostAlignStateLabel);
	ADD_MAP(TerminalStateLabel);

	std::cout << "test_profile_hmm_states: Testing properties and access methods of HMM states" << std::endl;
    for(std::map<std::string,AnyStateLabel>::iterator asl=stateLabels.begin(); asl!=stateLabels.end(); asl++)
	   std::cout << asl->first << ":  " << stateInfo(asl->second) << endl;


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
  } // End test_profile_hmm_states

BOOST_AUTO_TEST_SUITE_END()//end of test_non_db
