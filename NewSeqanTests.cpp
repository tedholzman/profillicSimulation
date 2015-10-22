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
#include "NewSeqanTests.hpp"

/**
 * \brief The next several lines are required for the test harness
 *
 */
#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE galosh_profuse_tests
#ifdef DEBUG
#define BOOST_RT_PARAM_DEBUG
#endif
#include <boost/test/included/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>
using boost::test_tools::output_test_stream;
using boost::debug::const_string;

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
//
//namespace po = boost::program_options;
po::options_description galosh_opt_desc("Galosh-specific unit test options");
po::variables_map galosh_opt_map;
po::parsed_options parsed(new po::options_description());

int main(int argc, char* argv[])
{
   g_argc = argc;
   g_argv = (char **) argv;
   /**
    * make sure booleans are displayed as True and False
    */
   std::cout << std::boolalpha;

   /**
    * All about program option processing
    */
   galosh_opt_desc.add_options()
      ("ghelp","Galosh specific help messages")
   ;
   parsed =
      po::command_line_parser(argc, argv).options(galosh_opt_desc).allow_unregistered().run();
   po::store(parsed, galosh_opt_map);
   po::notify(galosh_opt_map);

   initialize_globals();
   extern ::boost::unit_test::test_suite* init_unit_test_suite(int argc,char* argv[]);
   boost::unit_test::init_unit_test_func init_func = &init_unit_test_suite;

   /**
    * \var Max_Unit_Name_Width
    *
    * This is one of several collections of information about the current set of
    * of tests to be run
    */
   Max_Unit_Name_Width = 0;
   for(
      std::_Rb_tree_iterator<
	     std::pair<const long unsigned int,boost::unit_test::test_unit*>
	     > val = s_frk_impl().m_test_units.begin();
	  val != (s_frk_impl().m_test_units.end());
	  val++
   ) {
#ifdef DEBUG
		std::cerr << val->second->p_name.value << " " <<
	    val->first << " " << val->second->p_type_name<< std::endl;
#endif
		if(strcmp(((const_string)(val->second->p_type_name)).begin(),"case") == 0)
		{
           if(val->second->p_name.value.length() > ::Max_Unit_Name_Width)
        	      Max_Unit_Name_Width = val->second->p_name.value.length();
        }
	}
	::boost::unit_test::unit_test_main(init_func, argc, argv);
	return 1;
}
/**
 * \fn std::string getOptVal(po::parsed_options *p,const std::string optName)
 * \brief fetches the value of program option optName
 * \param optName
 *    The string value of command line parameter you want to fetch
 * \param p
 *    A pointer to the boost program option (po::) class that contains the parsed options
 * \brief This is a kludge, to be repaired when I fully understand the program_options boost module
 * \todo this is not good code.  It's got to be repaired when there's time to
 * understand the nature of the \a allow_unregistered option of the \a basic_command_line_parser.
 *
 * \return Returns the string value of option optName
 */
std::string getOptVal(po::parsed_options *p,const std::string optName) {
	std::string retVal = "";
	for (vector<po::basic_option<char> >::iterator it = (*p).options.begin(); it != (*p).options.end(); it++) {
		if (it->string_key.compare(optName) == 0)
			return it->value[0];
	}
	return retVal;
}
/**
 * \fn std::string getOpt(po::parsed_options *p,const std::string optName)
 * \brief tells whether the option, optName has been specified
 * \param optName
 *    The string value of command line parameter you want to fetch
 * \param p
 *    A pointer to the boost program option (po::) class that contains the parsed options
 *
 * \return Returns True if the option is specified (explicitly), False otherwise
 *
 * \todo See if this code works with implicit, default, and options specified
 * other than on the command line.
 */
bool getOpt(po::parsed_options *p, const std::string optName)
{
	for (vector<po::basic_option<char> >::iterator it = p->options.begin(); it != p->options.end(); it++) {
		if (it->string_key.compare(optName) == 0) return true;
	}
	return false;
}

/**
 * \def TEST_DESC(name,description)
 *
 * Macro for formatting minimal output from each test case.
 */
#define TEST_DESC(name,description) (string(name) + string(":") +\
string(::Max_Unit_Name_Width - (string(name)).length() + 1, ' ') +\
string(description))
/**
 * \brief This is a canonical example of a unit test.
 *
 * This test succeeds if \f$ 1 = 1\f$.
 * It also defines commandline option value.
 */
BOOST_AUTO_TEST_CASE( sanity_check ) {
	std::cout << TEST_DESC("Sanity check","should always be TRUE") << std::endl;
	BOOST_CHECK_MESSAGE( 1 == 1, "Expected TRUE, saw " << (1 == 1));
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
	std::cout << TEST_DESC("test_alphabets","Testing alphabet conversions") << std::endl;
	Dna a = 'a';
	BOOST_CHECK_MESSAGE( /// Test conversion of Dna to char
			(char) a == 'A',
			"Dna to string yielded "<<a<< ", expected A"
	);
	Dna5 b = 'f';//'f' is unknown character
	BOOST_CHECK_MESSAGE( ///Test conversion of Dna5 to char (and converting unknowns to N's)
			(char) b == 'N',
			"Dna5 yielded " << b << ", expected N"
	);
	//Many SeqAn alphabet classes can be converted into each other.
	b = a;
	BOOST_CHECK_MESSAGE( ///Test conversion of Dna5 to Dna to char
			(char) b == 'A',
			"Conversion from Dna to Dna5 yielded " << b << ", expected A"
	);
	Iupac c = a;
	BOOST_CHECK_MESSAGE( ///Test conversion of Iupac to Dna to char
			(char) c == 'A',
			"Conversion from Dna to Iupac yielded " << c << ", expected A"
	);
} // end of test_alphabets
bool ambiguity_tested = false;
/**
 * \brief Tests of ambiguous letters in seqan alphabets
 */
BOOST_AUTO_TEST_CASE( test_ambiguous )
{
	std::cout << TEST_DESC("test_ambiguity","Testing inter-alphabet ambiguities") << std::endl;
	ambiguity_tested=true;
	BOOST_CHECK_MESSAGE( ///Test if Dna is ambiguous over Dna5
			!isTrue( galosh::IsAmbiguous<seqan::Dna, seqan::Dna5>::Type() ),
			"Dna should not have ambiguities with respect to Dna5."
	);
	BOOST_CHECK_MESSAGE( ///Test if Dna5 is ambiguous over Dna
			isTrue( galosh::IsAmbiguous<seqan::Dna5, seqan::Dna>::Type() ),
			"Dna5 should have ambiguities with respect to Dna."
	);
	Dna dna_residue;
	Dna5 a = 'a';
	BOOST_CHECK_MESSAGE( ///Test how many elements in Dna5 a match an element in Dna
			galosh::ambiguousCount( a, Dna() ) == 1,
			"'A' in Dna5 matches " << (galosh::ambiguousCount( a, Dna() ) == 1) << " letters in Dna, should match 1"
	);
	galosh::ambiguousAssign( dna_residue, a, 0 );
	BOOST_CHECK_MESSAGE( ///Find a Dna letter matching a Dna5 letter
			(char) dna_residue == 'A',
			"ambiguousAssign - Dna5 'A' matches Dna " << dna_residue << ", expected 'A'"
	);
	Dna5 n = 'n';
	size_t num_elements = galosh::ambiguousCount( n, Dna() );
	BOOST_CHECK_MESSAGE( /// Test how many Dna elements match Dna5's letter n
			num_elements == 4,
			"Dna5's N element matches " << num_elements << " Dna elements, expected 4"
	);
	string nucleotides ("ATGC");
	/// Test that ambiguousAssign matches all four nucleotides against N
	for( int i = 0; i < num_elements; i++ ) {
		galosh::ambiguousAssign( dna_residue, n, i );
		size_t pos = nucleotides.find((char)dna_residue);
		BOOST_CHECK_MESSAGE(
				pos != string::npos,
				"ambiguousAssign returned " << (char)dna_residue << ", expected a unique nucleotide."
		);
		if(pos<nucleotides.length()) nucleotides.erase(pos,1);
	}
	BOOST_CHECK_MESSAGE(
			nucleotides.length() == 0,
			"ambiguousAssign should match all possible nucleotides; missed these: " + nucleotides
	);
} //end of test_ambiguous

/**
 * \brief Test of basic seqan and galosh sequence manipulation operations
 */
BOOST_AUTO_TEST_CASE( test_sequences )
{
	std::cout << TEST_DESC("test_sequences","Testing basic seqan and galosh biosequence operations") << std::endl;
	string anypeptide = "anypeptide";
	seqan::Peptide prot = anypeptide;
	BOOST_CHECK_MESSAGE( ///Test if seqan peptide sequence length is correct
			seqan::length(prot) == anypeptide.length(),
			"Length of '" << (*toCString(&prot)) << "' should be " << anypeptide.length() <<
			"; ended up being " << seqan::length(prot)
	);
	BOOST_CHECK_MESSAGE( ///Test array-like indexing of seqan peptide sequence
			(char)prot[ 9 ] == toupper(anypeptide[9]),
			"The 9th character of " << (*toCString(&prot)) << " should be " << anypeptide[9] <<
			"; ended up being " << (char)prot[9]
	);
	prot += "anewend"; anypeptide += "anewend";
	string protStr;
	seqan::assign(protStr,prot);
	BOOST_CHECK_MESSAGE( ///Test concatenation of string to seqan peptide sequence
			boost::to_upper_copy(anypeptide).compare(protStr) == 0,
			"Test sequence should = " << boost::to_upper_copy(anypeptide) <<
			"; ended up = " << (*toCString(&prot))
	);
	BOOST_CHECK_MESSAGE( ///Test length of concatenated seqan peptide sequence
			seqan::length(prot) == anypeptide.length(),
			"Test sequence should be " << anypeptide.length() << " characters long, " <<
			"; ended up being " << seqan::length(prot) << " characters long"
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
	seqan::String<Dna,seqan::Alloc<> > *sp = static_cast<seqan::String<Dna,seqan::Alloc<> > *>(&dna_seq);
	string dna_seq_str = toString(*sp);
	string anyDNA_s = (boost::to_upper_copy(anyDNA));
	BOOST_CHECK_MESSAGE( ///Test creation of galosh Sequence from string
			anyDNA_s.compare(dna_seq_str) == 0,
			"Test sequence should = " << boost::to_upper_copy(anyDNA) <<
			"; ended up being = " << dna_seq_str
	);
	BOOST_CHECK_MESSAGE( ///Test length of created galosh::Sequence
			seqan::length(dna_seq) == anyDNA.length(),
			"Test sequence should be " << anyDNA.length() << " characters long, " <<
			"; ended up being " << seqan::length(dna_seq) << " characters long"
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
struct test_fasta_options {
	test_fasta_options() {
		std::string tmp = getOptVal(&(::parsed),"ungapped_fasta");
		ungapped_fasta_file = (tmp.compare("") != 0) ? tmp : DEFAULT_UNGAPPED_FASTA_FILE;
		tmp = getOptVal(&(::parsed),"gapped_fasta");
		gapped_fasta_file   = (tmp.compare("") != 0) ? tmp : DEFAULT_GAPPED_FASTA_FILE;
	}
 	std::string ungapped_fasta_file;
	std::string gapped_fasta_file;
};

BOOST_FIXTURE_TEST_CASE( test_fasta, test_fasta_options )
{
	/**
	 * \todo Try to use the boost program options API "correctly".
	 * \todo currently using unit test module output_test_stream to check the validity
	 * of the fasta files; not just printing them out.  output_test_stream.match_pattern()
	 * kindof sucks:  the output must match the "pattern file" exactly, character
	 * for character.  In the case of fasta files, it might be better
	 * to use a checksum.  In any case, it would be nice if pattern_match used regular
	 * expressions or at least wildcards.
	 */

	std::cout << TEST_DESC("test_fasta","Testing ungapped and gapped fasta file reads and conversions") << std::endl;
	//// Fasta with no gaps
	galosh::Fasta<Dna> fasta;
	fasta.fromFile(ungapped_fasta_file);
   	output_test_stream output(ungapped_fasta_file, true );
	output << fasta;
	BOOST_CHECK( output.match_pattern() ); /// Test that the ungapped fasta file is correctly scanned?

	//// Fasta: with gaps.  Use "char" since basic seqan alphabets don't support gap chars.

	galosh::Fasta<char> aligned_fasta;
	aligned_fasta.fromFile(gapped_fasta_file);
    output_test_stream goutput(gapped_fasta_file, true );
    goutput << aligned_fasta;
	BOOST_CHECK( goutput.match_pattern() ); ///Test that Gapped fasta is file correctly scanned?

	/// \todo  I'm thinking of just scrapping Fasta altogether, using
	/// StringSet instead, with _loadSequences from
	/// seqan/graph_utils/utility_parsing.h, as in seqan_tcoffee.
	/// Eventually.  --Paul

} //end of test_fasta

/**
 * \brief Test basic operations with multinomial distributions of DNA
 *
 * Tests initial distributions and ambiguity codes.  Tests various ways of setting
 * distribution elements.  Also demonstrates, in the code proper, various ways
 * of accessing and displaying the elements.
 */
BOOST_AUTO_TEST_CASE( test_multinomials )
{
	/// \todo note: ambiguous test \a within multinomials <br>
	///can that be done within this model?  Can we tell if a test (e.g.
	///ambiguous) is going to be executed? Currently using kludge.

	std::cout << TEST_DESC("test_multinomials","Testing multinomial distribution generation and manipulation for DNA") << std::endl;

	Dna c_2( 1 );
    BOOST_CHECK_MESSAGE( ///Test conversion of Dna (created from int) to char
			(char) c_2 == 'C',
			"Dna to string yielded '"<< c_2 << "'; expecting 'C'"
	);

	galosh::MultinomialDistribution<Dna, realspace> dna_dist;
	dna_dist[ c_2 ] = .4;
	BOOST_CHECK_MESSAGE( ///Test modification of a Multinomial distribution object
			dna_dist.m_elementCount == 4 && dna_dist.m_probs[0] == 0.25 &&
			dna_dist.m_probs[1].prob() == 0.4 &&
			dna_dist.m_probs[2].prob() == 0.25 && dna_dist.m_probs[3].prob() == 0.25,
			"dist should be A=T=G=0.25, C=0.4; yielded " << dna_dist
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
	BOOST_CHECK_MESSAGE( ///Test structure of a StateLabel distribution
	   n_elements == 12 && elValOK,
	   "Expected " << evenDist << " for each of 12 states, saw " << state_dist
	);

	if( ambiguity_tested )
	{
		std::cout << TEST_DESC("test_multinomials","Testing multinomial distribution with ambiguities (for DNA)") << std::endl;
		Dna5 a = 'a';
		Dna5 n = 'n';

		double ambigN = dna_dist.ambiguousSum( n ).prob();
		BOOST_CHECK_MESSAGE( ///Check initial distributions of single Dna5 ambiguous and unambiguous elements
	       dna_dist.ambiguousSum( a ) == 0.25 && ambigN == 1.15,
           "Initial 'A' distribution should be 0.25; yielded " << dna_dist.ambiguousSum( a ) <<
           "; initial 'N' distribution should be 1.15 and it yielded " << dna_dist[ n ]
		);
        /// \internal
		/// Uncomment this, and it should fail to compile, since AminoAcid is
		/// not ambiguous over Dna.
		/// AminoAcid aminoacid_n = 'n';
		/// std::cout << "The dna_dist, accessed using AminoAcid 'n', returns " << dna_dist[ aminoacid_n ] << endl;
		output_test_stream output;
		dna_dist.ambiguousIncrement( n, 1.0 );
		output << dna_dist;
		BOOST_CHECK_MESSAGE( ///Test ambiguousIncrement of 1.0 to dna distribution
		    output.is_equal("(A=0.5,C=0.65,G=0.5,T=0.5)"),
		    "Distribution should equal (A=0.5,C=0.65,G=0.5,T=0.5); turns out to equal " << dna_dist
		);
		galosh::MultinomialDistribution<Dna, realspace>::AmbiguousValue<Dna5> ambiguous_value = dna_dist[ n ];
		ambiguous_value += 1.0;
		output << dna_dist;
		BOOST_CHECK_MESSAGE( ///Test autoincrement (+=)  of ambiguous DNA value (n), effect on distribution
		   output.is_equal("(A=0.75,C=0.9,G=0.75,T=0.75)"),
		   "Distribution should equal (A=0.75,C=0.9,G=0.75,T=0.75); turns out to equal " << dna_dist
		);
		/**
		 * \todo See if BOOST_CHECK_CLOSE works with non-float arithmetic.  It almost certainly won't.
		 * See how easy it is make it work.
		 */
        BOOST_CHECK_CLOSE(ambiguous_value.prob(),3.15,FLOAT_POINT_EQUALITY_TOL_PCT);  /// Test probability of autoincremented ambiguous element
		ambiguous_value = 1.0;
        output << dna_dist;
		BOOST_CHECK_MESSAGE( /// Test direct assignment of ambiguous DNA value (n), effect on distribution
		   output.is_equal("(A=0.25,C=0.25,G=0.25,T=0.25)"),
		   "Distribution should equal (A=0.25,C=0.25,G=0.25,T=0.25); yielded " << dna_dist
		);
        BOOST_CHECK_CLOSE(ambiguous_value.prob(),1.00,FLOAT_POINT_EQUALITY_TOL_PCT);  /// Test probability of autoincremented ambiguous element
		dna_dist[ n ] += 1.0;
		output << dna_dist;
		BOOST_CHECK_MESSAGE( /// Test direct autoincrement,  distrib[n] += 1
		   output.is_equal("(A=0.5,C=0.5,G=0.5,T=0.5)"),
		   "Distribution should equal (A=0.5,C=0.5,G=0.5,T=0.5); yielded " << dna_dist
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
struct test_profile_hmm_states_init {
	/**
	 * \var curCorrect
	 * store correct info about each state.
	 * \todo This concept is a little fragile.  When the algebra or methodology changes, the "correct answers" will
	 * change.  Re-think this later on.  Probably put the right answers into a file.
	 */
	std::map<std::string,std::string> curCorrect;
	test_profile_hmm_states_init ()
	{
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
	}
	~test_profile_hmm_states_init()
	{
	   curCorrect.clear();
	}
   }; // of fixture struct test_profile_hmm_states

   BOOST_FIXTURE_TEST_CASE(test_profile_hmm_states,test_profile_hmm_states_init)
   {
	/**
     * Iterate through all states, check each one against a known model.
     */
	   std::cout << TEST_DESC("test_profile_hmm_states","Testing properties and access methods of HMM states") << std::endl;
       for(std::map<std::string,AnyStateLabel>::iterator asl=stateLabels.begin(); asl!=stateLabels.end(); asl++) {
    	      BOOST_CHECK_MESSAGE(
    	         curCorrect[asl->first].compare(stateInfo(asl->second)) == 0,
    	         asl->first << "should be " << curCorrect[asl->first] << ";\nturns out to be " <<
    	            stateInfo(asl->second)
   	      );
       } //for
   } // End test_profile_hmm_states

BOOST_AUTO_TEST_SUITE_END()//end of test_non_dp

BOOST_AUTO_TEST_SUITE(dynamic_programming_1)

   struct test_profiles_options {
      /**
       * \todo  Can the following typedefs be parameterized?  Moreover, can they
       * be controlled by command-line options instead of being compiled in?
       */
	  typedef realspace ProbabilityType;
      typedef realspace ScoreType;
      typedef realspace MatrixValueType;
      typedef ProfileTreeRoot<Dna, ProbabilityType> ProfileType;

	  bool use_del_in_del_out;
	  output_test_stream output;

	  //constructor
	  test_profiles_options()
	  {
         use_del_in_del_out = getOpt(&(::parsed),"use_del_in_del_out");
	  }
    }; // of fixture struct test_profiles_options

#define GALOSH_SINGLE_VAR_TEST(variable,teststream,answer,label) teststream << variable;\
BOOST_CHECK_MESSAGE(teststream.is_equal(answer),label answer "; equals " << variable);
    BOOST_FIXTURE_TEST_CASE(test_profiles,test_profiles_options)
    {
       std::cout << TEST_DESC("Test Profiles","Basic operations on profiles") << std::endl;
  	   using namespace galosh;

       // We need at least 3 positions to be able to test extensions of del-ins and del-outs
       galosh::ProfileTreeRoot<Dna, ProbabilityType> dna_profile(use_del_in_del_out ? 3 : 2 );
       ///////
       /// Profile
       GALOSH_SINGLE_VAR_TEST( ///Check default DNA profile
          dna_profile,output,
       	  "[ M->(M=0.3333333,I=0.3333333,D=0.3333333), I->(M=0.5,I=0.5), D->(M=0.5,D=0.5), I:(A=0.25,C=0.25,G=0.25,T=0.25), N->(N=0.5,B=0.5), B->(M=0.5,D=0.5), C->(C=0.5,T=0.5) ]\n[ M:(A=0.25,C=0.25,G=0.25,T=0.25) ]\n[ M:(A=0.25,C=0.25,G=0.25,T=0.25) ]\n",
     	  "Default DNA profile: Should equal:\n"
       );

       galosh::ProfileTreeRoot<Iupac, ProbabilityType> iupac_profile;
       GALOSH_SINGLE_VAR_TEST(///Check default IUPAC profile
          iupac_profile,output,
          "[ M->(M=0.3333333,I=0.3333333,D=0.3333333), I->(M=0.5,I=0.5), D->(M=0.5,D=0.5), I:(U=0.0625,T=0.0625,A=0.0625,W=0.0625,C=0.0625,Y=0.0625,M=0.0625,H=0.0625,G=0.0625,K=0.0625,R=0.0625,D=0.0625,S=0.0625,B=0.0625,V=0.0625,N=0.0625), N->(N=0.5,B=0.5), B->(M=0.5,D=0.5), C->(C=0.5,T=0.5) ]",
          "Default IUPAC DNA profile: Should equal:\n"
       );

       std::cout << "Reading profile from file 'seqantest.DATA.profile'" << std::endl;
       galosh::ProfileTreeRoot<Dna, ProbabilityType> dna_profile_from_file;
       dna_profile_from_file.fromFile( "seqantest.DATA.profile" );
       std::cout << "\tgot:" << std::endl;
       std::cout << dna_profile_from_file;
       std::cout << endl;


       Dna residue;
       galosh::Sequence<Dna5> dna_seq_a = "a";
       galosh::Sequence<Dna5> dna_seq_b = "ag";
       galosh::Sequence<Dna5> dna_seq_c = "aga";
       vector<galosh::Sequence<Dna5> > dna_seqs( 3 );
       dna_seqs[ 0 ] = dna_seq_a;
       dna_seqs[ 1 ] = dna_seq_b;
       dna_seqs[ 2 ] = dna_seq_c;
       int num_sequences_to_use = dna_seqs.size();
       /*
       galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer
          forward_matrices(
             dna_profile,
             dna_seqs,
             num_sequences_to_use
          );
       galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer::iterator forward_matrices_iterator =
          forward_matrices.begin();
       galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer backward_matrices(
          dna_profile,
          dna_seqs,
          num_sequences_to_use
       );
       galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer::reverse_iterator backward_matrices_iterator =
          backward_matrices.rbegin();
       galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::PositionSpecificSequenceScoreCoefficientsVector coefficients_vector( 3 );
       galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::PositionEntente position_entente;
       galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::PositionEntente position_entente_unscaled;
       galosh::DynamicProgramming<Dna, ProbabilityType, ScoreType, MatrixValueType>::PositionEntente position_entente_backup;
       */
   }
       
BOOST_AUTO_TEST_SUITE_END()//end of dynamic_programming_1
