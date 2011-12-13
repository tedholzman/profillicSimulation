#include "Algebra.hpp"
#include "Profile.hpp"
#include "Fasta.hpp"

#include <iostream>

#include <seqan/basic.h>

#ifdef __HAVE_MUSCLE
int g_argc;
char **g_argv;
#endif // __HAVE_MUSCLE

using namespace seqan;

namespace galosh {

/////////////
/**
 * Given a pointer to a HMMer profile (created using hmmer::AllocPlan7Shell()),
 * allocate its body and populate it with data from the given galosh Profile.
 * Later, you should free the hmm with FreePlan7( &hmm ).
 */
template <class ProfileType, class ResidueType>
void
profileToConsensus (
  ProfileType const & profile,
  Sequence<ResidueType> & sequence
)
{
  uint32_t profile_length = profile.length();
  sequence.reinitialize( profile_length );
  for( uint32_t pos_i = 0; pos_i < profile_length; pos_i++ ) {
    sequence[ pos_i ] =
      profile[ pos_i ][ Emission::Match ].maximumValueType();
  }
  // That's it.
  return;
} // profileToConsensus( ProfileType const &, Sequence & )

} // End namespace galosh

int
main ( int const argc, char const ** argv )
{
#ifdef __PROFUSE_USE_AMINOS
  typedef seqan::AminoAcid20 ResidueType;
#else // __PROFUSE_USE_AMINOS .. else
  typedef seqan::Dna ResidueType;
#endif // __PROFUSE_USE_AMINOS .. else ..

  if( argc < 2 ) {
    cout << "Usage: " << argv[ 0 ] << " <input (galosh Profile) filename> [<output (Fasta) filename>]" << endl;
    exit( 1 );
  }
  const bool be_verbose = ( argc >= 3 );
  if( be_verbose ) {
    cout << "Reading profile from file '" << argv[ 1 ] << "'" << endl;
  }
  galosh::ProfileTreeRoot<ResidueType, floatrealspace> profile;
  profile.fromFile( argv[ 1 ] );
  if( be_verbose ) {
    cout << "\tgot:" << std::endl;
    cout << profile;
    cout << endl;
  }

  galosh::Fasta<ResidueType> fasta( 1 );
  fasta.m_descriptions[ 0 ] = "Consensus"; // TODO: MAGIC #
  profileToConsensus( profile, fasta[ 0 ] );

  if( argc >= 3 ) {
    if( be_verbose ) {
      cout << "Writing Fasta to file '" << argv[ 2 ] << "'" << endl;
    }
    std::ofstream fasta_stream( argv[ 2 ] );
    assert( fasta_stream.good() );
    fasta_stream << fasta;
    fasta_stream.close();
    if( be_verbose ) {
      cout << "\tdone." << endl;
    }
  } else {
    if( be_verbose ) {
      cout << "Fasta is:" << endl;
    }
    cout << fasta;
    cout << endl;
  }

  exit( 0 );
} // main (..)
