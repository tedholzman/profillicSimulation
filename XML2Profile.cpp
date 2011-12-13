#include "Algebra.hpp"
#include "ProfileTree.hpp"

#include <iostream>

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include "boost/filesystem.hpp"
namespace fs = boost::filesystem;

#include <seqan/basic.h>

#ifdef __HAVE_MUSCLE
int g_argc;
char **g_argv;
#endif // __HAVE_MUSCLE

using namespace galosh;

namespace galosh {

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

} // End namespace galosh

int
main ( int const argc, char const ** argv )
{
  typedef floatrealspace ProbabilityType;
#ifdef __PROFUSE_USE_AMINOS
  typedef seqan::AminoAcid20 ResidueType;
#else // __PROFUSE_USE_AMINOS .. else
  typedef seqan::Dna ResidueType;
#endif // __PROFUSE_USE_AMINOS .. else ..
  typedef ProfileTreeRoot<ResidueType, ProbabilityType> RootType;
  typedef ProfileTreeRoot<ResidueType, ProbabilityType> InternalNodeType;
  typedef ProfileTree<ResidueType, ProbabilityType, InternalNodeType > ProfileTreeType;

  // For now we assume a Dna distribution.  TODO: Generalize.
  if( argc < 2 ) {
    cout << "Usage: " << argv[ 0 ] << " <input ProfileTree.xml filename> [<output (galosh Profile) filename> ]" << endl;
    exit( 1 );
  }
  const bool be_verbose = true; //( argc >= 2 );
  if( be_verbose ) {
    cout << "Reading ProfileTree.xml from file '" << argv[ 1 ] << "'" << endl;
  }
  ProfileTreeType profile_tree;
  readXML( profile_tree, argv[ 1 ] );
  if( profile_tree.getProfileTreeRoot()->length() == 0 ) {
    cout << "No profiles were found in the ProfileTree.xml file '" << argv[ 1 ] << "'" << endl;
    return 1;
  } else if( profile_tree.nodeCount() > 1 ) {
    if( be_verbose ) {
      cout << "WARNING: Using only the root profile in the given ProfileTree.xml file." << endl;
    }
  }
  //if( be_verbose ) {
  //  cout << "\tgot:" << std::endl;
  //  cout << *( profile_tree.getProfileTreeRoot() );
  //  cout << endl;
  //}

  if( argc >= 3 ) {
    if( be_verbose ) {
      cout << "Writing Profile to file '" << argv[ 2 ] << "'" << endl;
    }
    std::ofstream profile_stream( argv[ 2 ] );
    assert( profile_stream.good() );
    profile_stream << *( profile_tree.getProfileTreeRoot() );
    profile_stream.close();
    if( be_verbose ) {
      cout << "\tdone." << endl;
    }
  } else {
    if( be_verbose ) {
      cout << "Profile is:" << endl;
    }
    cout << *( profile_tree.getProfileTreeRoot() );
    cout << endl;
  }

  exit( 0 );
} // main (..)
