#include "Galosh.hpp"
#include "Parameters.hpp"
#include "HMMoC-BFloat-Algebra/Algebra.hpp" //TAH 11/4/2011

#include <sstream> // For stringstream.
#include <ctime> // for std::time

// TODO: REMOVE ..
#include <iostream>
using namespace galosh;

/////////////
// Required for linking with muscle (even though we don't access it from ProfileTrainer.cpp).
#ifdef __HAVE_MUSCLE
int g_argc;
char **g_argv;
#endif // __HAVE_MUSCLE

int
test_algebras ( /*supportedTypesPtr_t supported_types, */ int argc, char **argv )
{
  clock_t start = clock();

  bfloat bfloat0 = 1.0;
  cout << "[0] bfloat0: " << bfloat0 << endl;
  bfloat0 *= 1E-10;
  cout << "[1] bfloat0 after *= 1E-10: " << bfloat0 << endl;
  cout << "[1] log(bfloat0) is now: " << log( bfloat0 ) << endl;
  for( int i = 2; i <= 30; i++ ) {
    bfloat0 *= bfloat0;
    cout << "[" << i << "] bfloat0 after *= bfloat0: " << bfloat0 << endl;
    cout << "[" << i << "] log(bfloat0) is now: " << log( bfloat0 ) << endl;
    cout << "[" << i << "] 1.0 - bfloat0 is now: " << ( 1.0 - bfloat0 ) << endl;
  }
  
  bfloat bfloat1 = 1.0;
  cout << "[0] bfloat1: " << bfloat1 << endl;
  bfloat1 *= 1E-1;
  cout << "[1] bfloat1 after *= 1E-1: " << bfloat1 << endl;
  cout << "[1] log(bfloat1) is now: " << log( bfloat1 ) << endl;
  for( int i = 2; i <= 10; i++ ) {
    bfloat1 *= 1E-1;
    cout << "[" << i << "] bfloat1 after *= 1E-1: " << bfloat1 << endl;
    cout << "[" << i << "] log(bfloat1) is now: " << log( bfloat1 ) << endl;
    cout << "[" << i << "] 1.0 - bfloat1 is now: " << ( 1.0 - bfloat1 ) << endl;
  }

  cout << "largest float: " << numeric_limits<float>::max() << endl;
  cout << "smallest float: " << numeric_limits<float>::min() << endl;

  bfloat bfloat2 = ( bfloat1 * 1E-5 );
  cout << " bfloat1 is " << bfloat1 << endl;
  cout << " bfloat2 is " << bfloat2 << endl;
  cout << " ( bfloat2 - bfloat1 ) = " << ( bfloat2 - bfloat1 ) << endl;
  cout << " ( bfloat1 - bfloat2 ) = " << ( bfloat1 - bfloat2 ) << endl;
  cout << " bfloat0 is " << bfloat0 << endl;
  cout << " ( bfloat0 - bfloat1 ) = " << ( bfloat0 - bfloat1 ) << endl;
  cout << " ( bfloat1 - bfloat0 ) = " << ( bfloat1 - bfloat0 ) << endl;

  clock_t end = clock();
  double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  cout << "Elapsed CPU time: " << cpu_time_used << endl;
  return 0;
} // test_algebras(..)

int
main ( int argc, char **argv )
{
//  // Support only DNA sequences
//  //PolySequence::SupportedTypes supportedTypes( true );
//
//  // Support only AminoAcid sequences
//  PolySequence::SupportedTypes supportedTypes( false, true );

  return //test_algebras( &supportedTypes, argc, argv );
    test_algebras( argc, argv );

} // main( int, char** )


