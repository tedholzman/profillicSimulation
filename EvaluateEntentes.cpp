#include "Algebra.hpp"
#include "Profile.hpp"
#include "Fasta.hpp"
#include "DynamicProgramming.hpp"
#include "ProfileTree.hpp"
#include "ProfileTreeTrainer.hpp"
#include "Random.hpp"

#include <iostream>

#include <seqan/basic.h>

// TODO: REMOVE? Testing
#include <alglib/ap.h>
#include <alglib/statistics.h>

#ifdef __HAVE_MUSCLE
int g_argc;
char **g_argv;
#endif // __HAVE_MUSCLE

using namespace seqan;

namespace galosh {

/////////////
/**
 * Shuffle the values of the array in-place.  Assumptions: the array is of
 * length at least size.  This uses the given Random object to generate random
 * unsigned long integers.
 */
template <class ArrayType, class ValueType>
void
fisher_yates_shuffle (
  Random & random,
  ArrayType & array,
  uint32_t const & size,
  ValueType const & value_tag
) {
  uint32_t j;
  ValueType tmp;
  for( uint32_t i = ( size - 1 ); --i; ) {
    j = random.nextUniform( i );
    if( i == j ) {
      // Then don't move it.
      continue;
    }
    tmp = array[ i ];
    array[ i ] = array[ j ];
    array[ j ] = tmp;
  } // End foreach index i
} // fisher_yates_shuffle(..)

/**
 * Given an iterator into a sorted array and a value (not necessarily in the
 * array), return the lowest index that the value would have if it were added
 * into the array and the array were resorted.  Put another way, count how many
 * times we need to increment the iterator before the given value is less than
 * the iterator value (or within epsilon of it).
 */
template <class IteratorType, class ValueType>
uint32_t
locate_value_in_sorted_array (
  IteratorType it,
  IteratorType const & it_end,
  ValueType const & value,
  ValueType const & epsilon = 0
) {
  uint32_t location;
  // TODO: REMOVE
  //cout << "EPSILON: " << epsilon << endl;
  for( location = 0; it != it_end; ++it, ++location ) {
    if( ( value <= *it ) || ( epsilon && ( ( ( value - *it ) <= epsilon ) ) ) ) {
      // TODO: REMOVE
      //cout << "value: " << value << ", *it: " << *it << ", ( value - *it ): " << ( *it - value ) << endl;
      // TODO: REMOVE
      //exit( 0 );
      return location;
    }
  } // End foreach location
  return location; // This means that it is larger than every value.
} // locate_value_in_sorted_array(..)

} // End namespace galosh

int
main ( int const argc, char const ** argv )
{
#ifdef __PROFUSE_USE_AMINOS
  typedef seqan::AminoAcid20 ResidueType;
  typedef seqan::AminoAcid SequenceResidueType;
#else // __PROFUSE_USE_AMINOS .. else
  typedef seqan::Dna ResidueType;
  typedef seqan::Iupac SequenceResidueType;
#endif // __PROFUSE_USE_AMINOS .. else ..

  if( argc < 3 ) {
    cout << "Usage: " << argv[ 0 ] << " <input (galosh Profile) filename> <input Fasta filename> [<num seqs in group 1> [<number of draws> [<random seed>]]]" << endl;
    exit( 1 );
  }

  //typedef bfloat ProbabilityType;
  //typedef logspace ProbabilityType;
  //typedef floatrealspace ProbabilityType;
  typedef doublerealspace ProbabilityType;
  
  typedef bfloat ScoreType; // Preferred
  //typedef logspace ScoreType; // SLOWer than bfloat
  //typedef realspace ScoreType; // Only for very few & small sequences

  // MatrixValueType should be bfloat or logspace (at least until we fix Rabiner Scaling)
  typedef bfloat MatrixValueType;
  //typedef logspace MatrixValueType;

  // if be_verbose, should we also print out the input profile?
  const bool be_verbose_show_profile = false;
  // if be_verbose, should we also print out the input seqs?
  const bool be_verbose_show_sequences = false;

  const bool be_verbose = ( argc >= 7 );
  if( be_verbose ) {
    cout << "Reading profile from file '" << argv[ 1 ] << "'" << endl;
  }
  typedef galosh::ProfileTreeRoot<ResidueType, ProbabilityType> ProfileType;
  typedef galosh::ProfileTreeRoot<ResidueType, ProbabilityType> InternalNodeType;

  ProfileType profile;
  profile.fromFile( argv[ 1 ] );
  if( be_verbose ) {
    if( be_verbose_show_profile ) {
      cout << "\tgot:" << std::endl;
      cout << profile;
      cout << endl;
    } else{ 
      cout << "\tdone." << endl;
    }
  }

  galosh::Fasta<SequenceResidueType> fasta;
  if( be_verbose ) {
    cout << "Reading sequences from Fasta file '" << argv[ 2 ] << "'" << endl;
  }
  fasta.fromFile( argv[ 2 ] );
  if( be_verbose ) {
    if( be_verbose_show_sequences ) {
      cout << "\tgot:" << endl;
      cout << fasta;
      cout << endl;
    } else {
      cout << "\tdone." << endl;
    }
  } // End if be_verbose

  uint32_t num_seqs_in_group_1 = 0;
  if( argc >= 4 ) {
    try {
      num_seqs_in_group_1 = boost::lexical_cast<uint32_t>( argv[ 3 ] );
    } catch( boost::bad_lexical_cast & ) {
      std::cerr << "Unable to interpret the argument '" << argv[ 3 ] << "' as an unsigned long value for use as the number of sequences to include in group 1." << std::endl;
      exit( 1 );
    } // End try .. catch block for lexical_cast
  } // End if argc >= 4

  uint32_t num_draws = 0;
  if( argc >= 5 ) {
    try {
      num_draws = boost::lexical_cast<uint32_t>( argv[ 4 ] );
    } catch( boost::bad_lexical_cast & ) {
      std::cerr << "Unable to interpret the argument '" << argv[ 4 ] << "' as an unsigned long value for use as the number of sequences to generate." << std::endl;
      exit( 1 );
    } // End try .. catch block for lexical_cast
  } // End if argc >= 5

  uint32_t random_seed = static_cast<uint32_t>( std::time( NULL ) );
  if( argc >= 6 ) {
    try {
      random_seed = boost::lexical_cast<uint32_t>( argv[ 5 ] );
    } catch( boost::bad_lexical_cast & ) {
      std::cerr << "Unable to interpret the argument '" << argv[ 5 ] << "' as an unsigned long value for use as the random seed." << std::endl;
      exit( 1 );
    } // End try .. catch block for lexical_cast
  } // End if argc >= 6
  galosh::Random random( random_seed );

  galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType> dp;
  galosh::ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType>::Parameters parameters;

  // This is now NECESSARY as there is a bug in the rabiner scaling code.  Actually false is now the default so it's here to remind us to FIX RABINER SCALING!
  parameters.useRabinerScaling = false;

  galosh::ProfileTree<ResidueType, ProbabilityType, InternalNodeType> profile_tree( profile.length() );
  profile_tree.getProfileTreeRoot()->copyFrom( profile );
  galosh::ProfileTreeTrainer<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType> tree_trainer(
    fasta,
    &profile_tree
  );
  tree_trainer.m_parameters = parameters;

  tree_trainer.restart();

  // Calculate the alignment profiles
  vector<galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile> aps;
  tree_trainer.calculateAlignmentProfiles( 0, aps ); // Node 0 (it's a stump, not a tree)

  // TODO: REMOVE
  //cout << "Alignment profiles are: " << endl;
  //for( uint32_t ap_i = 0; ap_i < aps.size(); ap_i++ ) {
  //  cout << "\tAlignment profile " << ap_i << ": " << endl;
  //  cout << aps[ ap_i ] << endl;
  //}

  //exit( 0 );

  // TODO: ERE I AM.  Can we compute what the profile would be for each group separately, and compute the entropy of it, and/or the distances to it?
  static const bool do_explorations = false;
  if( do_explorations &&  num_seqs_in_group_1 && ( num_seqs_in_group_1 != aps.size() ) ) {

    uint32_t const num_seqs_in_group_2 = ( aps.size() - num_seqs_in_group_1 );

    // Calculate the group-specific profiles
    ProfileType group_1_profile( profile.length() ); // Will be recalculated for each draw
    ProfileType group_2_profile( profile.length() ); // Will be recalculated for each draw
    ProfileType group_1_profile_real( profile.length() ); // Calc'd once, for "draw 0" (real order)
    ProfileType group_2_profile_real( profile.length() ); // Calc'd once, for "draw 0" (real order)

    vector<uint32_t> shuffled_sequence_index_map( aps.size() );
    for( uint32_t i = 0; i < aps.size(); i++ ) {
      shuffled_sequence_index_map[ i ] = i;
    }

    // Euclidean distances, per pos, then per draw
    vector<vector<double> > group_1_euclidean_dist_to_profile( profile.length() );
    vector<vector<double> > group_2_euclidean_dist_to_profile( profile.length() );
    vector<vector<double> > group_1_euclidean_dist_to_group_2( profile.length() );
    for( uint32_t pos_i = 0; pos_i < profile.length(); pos_i++ ) {
      group_1_euclidean_dist_to_profile[ pos_i ].resize( num_draws + 1 );
      group_2_euclidean_dist_to_profile[ pos_i ].resize( num_draws + 1 );
      group_1_euclidean_dist_to_group_2[ pos_i ].resize( num_draws + 1 );
    } // End foreach pos_i  
    vector<double> group_1_euclidean_dist_to_profile_overall( num_draws + 1 );
    vector<double> group_2_euclidean_dist_to_profile_overall( num_draws + 1 );
    vector<double> group_1_euclidean_dist_to_group_2_overall( num_draws + 1 );
    // Skl distances, per pos, then per draw
    vector<vector<double> > group_1_skl_dist_to_profile( profile.length() );
    vector<vector<double> > group_2_skl_dist_to_profile( profile.length() );
    vector<vector<double> > group_1_skl_dist_to_group_2( profile.length() );
    for( uint32_t pos_i = 0; pos_i < profile.length(); pos_i++ ) {
      group_1_skl_dist_to_profile[ pos_i ].resize( num_draws + 1 );
      group_2_skl_dist_to_profile[ pos_i ].resize( num_draws + 1 );
      group_1_skl_dist_to_group_2[ pos_i ].resize( num_draws + 1 );
    } // End foreach pos_i  
    vector<double> group_1_skl_dist_to_profile_overall( num_draws + 1 );
    vector<double> group_2_skl_dist_to_profile_overall( num_draws + 1 );
    vector<double> group_1_skl_dist_to_group_2_overall( num_draws + 1 );

    galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::DistanceMatrix<doublerealspace> dist_or_prox_matrix;

    vector<galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile> profiles_as_aps( 3 );
    profiles_as_aps[ 0 ].reinitialize( profile.length() );
    profiles_as_aps[ 1 ].reinitialize( group_1_profile.length() );
    profiles_as_aps[ 2 ].reinitialize( group_2_profile.length() );

    // "draw" 0 is the unshuffled values.
    for( uint32_t draw_i = 0; draw_i < ( num_draws + 1 ); draw_i++ ) {
      if( draw_i > 0 ) {
        galosh::fisher_yates_shuffle( random, shuffled_sequence_index_map, aps.size(), ( uint32_t )0 );
        if( 0 ) {
          cout << "[Draw " << draw_i << "] " << "SHUFFLED INDICES: ( " << shuffled_sequence_index_map[ 0 ];
          for( uint32_t i = 1; i < aps.size(); i++ ) {
            cout << ", " << shuffled_sequence_index_map[ i ];
          }
          cout << " )" << endl;
        }
      } // End if draw_i > 0

      group_1_profile.zero();
      group_2_profile.zero();
      for( uint32_t ap_i = 0; ap_i < num_seqs_in_group_1; ap_i++ ) {
        for( uint32_t ap_pos_i = 1; ap_pos_i < aps[ 0 ].length(); ap_pos_i++ ) {
          group_1_profile[ ap_pos_i - 1 ][ galosh::Emission::Match ] +=
            aps[ shuffled_sequence_index_map[ ap_i ] ][ ap_pos_i ][ galosh::Emission::Match ];
        }
      }
      group_1_profile.normalizePositions( 1E-5 ); //0 );
      for( uint32_t ap_i = num_seqs_in_group_1; ap_i < num_seqs_in_group_2; ap_i++ ) {
        for( uint32_t ap_pos_i = 1; ap_pos_i < aps[ 0 ].length(); ap_pos_i++ ) {
          group_2_profile[ ap_pos_i - 1 ][ galosh::Emission::Match ] +=
            aps[ shuffled_sequence_index_map[ ap_i ] ][ ap_pos_i ][ galosh::Emission::Match ];
        }
      }
      group_2_profile.normalizePositions( 1E-5 ); //0 );

      if( draw_i == 0 ) {
        group_1_profile_real.copyFrom( group_1_profile );
        group_2_profile_real.copyFrom( group_2_profile );
      }
  
      // Euclidean distances

      // Now calculate per-position distances 
      for( uint32_t pos_i = 0; pos_i < profile.length(); pos_i++ ) {
        group_1_euclidean_dist_to_profile[ pos_i ][ draw_i ] =
          profile[ pos_i ][ galosh::Emission::Match ].euclideanDistance( 
            group_1_profile[ pos_i ][ galosh::Emission::Match ]
          );
        group_2_euclidean_dist_to_profile[ pos_i ][ draw_i ] =
          profile[ pos_i ][ galosh::Emission::Match ].euclideanDistance( 
            group_2_profile[ pos_i ][ galosh::Emission::Match ]
          );
        group_1_euclidean_dist_to_group_2[ pos_i ][ draw_i ] =
          group_1_profile[ pos_i ][ galosh::Emission::Match ].euclideanDistance( 
            group_2_profile[ pos_i ][ galosh::Emission::Match ]
          );
        if( 0 ) {
          cout << "[Draw " << draw_i << "] " << "[" << pos_i << "] " << "Euclidean distance b/n profile & group 1: " << group_1_euclidean_dist_to_profile[ pos_i ][ draw_i ] << endl;
          cout << "[Draw " << draw_i << "] " << "[" << pos_i << "] " << "Euclidean distance b/n profile & group 2: " << group_2_euclidean_dist_to_profile[ pos_i ][ draw_i ] << endl;
          cout << "[Draw " << draw_i << "] " << "[" << pos_i << "] " << "Euclidean distance b/n group 1 & group 2: " << group_1_euclidean_dist_to_group_2[ pos_i ][ draw_i ] << endl;
        }
      } // End foreach pos_i
      group_1_euclidean_dist_to_profile_overall[ draw_i ] =
        profile.euclideanDistance( group_1_profile );
      group_2_euclidean_dist_to_profile_overall[ draw_i ] =
        profile.euclideanDistance( group_2_profile );
      group_1_euclidean_dist_to_group_2_overall[ draw_i ] =
        group_1_profile.euclideanDistance( group_2_profile );
      if( 0 ) {
        cout << "[Draw " << draw_i << "] " << "[Overall] Euclidean distance b/n profile & group 1: " << group_1_euclidean_dist_to_profile_overall[ draw_i ] << endl;
        cout << "[Draw " << draw_i << "] " << "[Overall] Euclidean distance b/n profile & group 2: " << group_2_euclidean_dist_to_profile_overall[ draw_i ] << endl;
        cout << "[Draw " << draw_i << "] " << "[Overall] Euclidean distance b/n group 1 & group 2: " << group_1_euclidean_dist_to_group_2_overall[ draw_i ] << endl;
      }
  
      ////// --------
      // TODO: Make these work:
      //profiles_as_aps[ 0 ] = profile;
      //profiles_as_aps[ 1 ] = group_1_profile;
      //profiles_as_aps[ 2 ] = group_2_profile;
      // TODO: Is this necessary?
      //profiles_as_aps[ 0 ].zero();
      //profiles_as_aps[ 1 ].zero();
      //profiles_as_aps[ 2 ].zero();
      for( uint32_t pos_i = 0; pos_i < profile.length(); pos_i++ ) {
        profiles_as_aps[ 0 ][ pos_i ][ galosh::Emission::Match ] = profile[ pos_i ][ galosh::Emission::Match ];
        profiles_as_aps[ 1 ][ pos_i ][ galosh::Emission::Match ] = group_1_profile[ pos_i ][ galosh::Emission::Match ];
        profiles_as_aps[ 2 ][ pos_i ][ galosh::Emission::Match ] = group_2_profile[ pos_i ][ galosh::Emission::Match ];
      }
  
      // Skl distances

      // Now calculate per-position distances 
      for( uint32_t pos_i = 0; pos_i < profile.length(); pos_i++ ) {
        dist_or_prox_matrix.setToSymmeterizedKullbackLeiblerDivergences( profiles_as_aps, pos_i );
  
        group_1_skl_dist_to_profile[ pos_i ][ draw_i ] =
          toDouble( dist_or_prox_matrix( 0, 1 ) );
        group_2_skl_dist_to_profile[ pos_i ][ draw_i ] =
          toDouble( dist_or_prox_matrix( 0, 2 ) );
        group_1_skl_dist_to_group_2[ pos_i ][ draw_i ] =
          toDouble( dist_or_prox_matrix( 1, 2 ) );
        if( 0 ) {
          cout << "[Draw " << draw_i << "] " << "[" << pos_i << "] " << "Skl distance b/n profile & group 1: " << group_1_skl_dist_to_profile[ pos_i ][ draw_i ] << endl;
          cout << "[Draw " << draw_i << "] " << "[" << pos_i << "] " << "Skl distance b/n profile & group 2: " << group_2_skl_dist_to_profile[ pos_i ][ draw_i ] << endl;
          cout << "[Draw " << draw_i << "] " << "[" << pos_i << "] " << "Skl distance b/n group 1 & group 2: " << group_1_skl_dist_to_group_2[ pos_i ][ draw_i ] << endl;
        }
      } // End foreach pos_i
      dist_or_prox_matrix.setToSymmeterizedKullbackLeiblerDivergences( profiles_as_aps );
      group_1_skl_dist_to_profile_overall[ draw_i ] =
        toDouble( dist_or_prox_matrix( 0, 1 ) );
      group_2_skl_dist_to_profile_overall[ draw_i ] =
        toDouble( dist_or_prox_matrix( 0, 2 ) );
      group_1_skl_dist_to_group_2_overall[ draw_i ] =
        toDouble( dist_or_prox_matrix( 1, 2 ) );
      if( 0 ) {
        cout << "[Draw " << draw_i << "] " << "[Overall] Skl distance b/n profile & group 1: " << group_1_skl_dist_to_profile_overall[ draw_i ] << endl;
        cout << "[Draw " << draw_i << "] " << "[Overall] Skl distance b/n profile & group 2: " << group_2_skl_dist_to_profile_overall[ draw_i ] << endl;
        cout << "[Draw " << draw_i << "] " << "[Overall] Skl distance b/n group 1 & group 2: " << group_1_skl_dist_to_group_2_overall[ draw_i ] << endl;
      }
    } // End foreach draw_i

    // Sort all of the draws in each array (except for draw 0, which is the real data)
    for( uint32_t pos_i = 0; pos_i < profile.length(); pos_i++ ) {
      std::sort( ( group_1_euclidean_dist_to_profile[ pos_i ].begin() + 1 ), group_1_euclidean_dist_to_profile[ pos_i ].end() );
      std::sort( ( group_2_euclidean_dist_to_profile[ pos_i ].begin() + 1 ), group_2_euclidean_dist_to_profile[ pos_i ].end() );
      std::sort( ( group_1_euclidean_dist_to_group_2[ pos_i ].begin() + 1 ), group_1_euclidean_dist_to_group_2[ pos_i ].end() );
    } // End foreach pos_i
    std::sort( ( group_1_euclidean_dist_to_profile_overall.begin() + 1 ), group_1_euclidean_dist_to_profile_overall.end() );
    std::sort( ( group_2_euclidean_dist_to_profile_overall.begin() + 1 ), group_2_euclidean_dist_to_profile_overall.end() );
    std::sort( ( group_1_euclidean_dist_to_group_2_overall.begin() + 1 ), group_1_euclidean_dist_to_group_2_overall.end() );

    for( uint32_t pos_i = 0; pos_i < profile.length(); pos_i++ ) {
      std::sort( ( group_1_skl_dist_to_profile[ pos_i ].begin() + 1 ), group_1_skl_dist_to_profile[ pos_i ].end() );
      std::sort( ( group_2_skl_dist_to_profile[ pos_i ].begin() + 1 ), group_2_skl_dist_to_profile[ pos_i ].end() );
      std::sort( ( group_1_skl_dist_to_group_2[ pos_i ].begin() + 1 ), group_1_skl_dist_to_group_2[ pos_i ].end() );
    } // End foreach pos_i
    std::sort( ( group_1_skl_dist_to_profile_overall.begin() + 1 ), group_1_skl_dist_to_profile_overall.end() );
    std::sort( ( group_2_skl_dist_to_profile_overall.begin() + 1 ), group_2_skl_dist_to_profile_overall.end() );
    std::sort( ( group_1_skl_dist_to_group_2_overall.begin() + 1 ), group_1_skl_dist_to_group_2_overall.end() );

    uint32_t location;
    double bothtails, lefttail, righttail;

    // Euclidean distance
    double distance_epsilon = 1E-9;
    for( uint32_t pos_i = 0; pos_i < profile.length(); pos_i++ ) {
      if( profile[ pos_i ][ galosh::Emission::Match ] != group_1_profile_real[ pos_i ][ galosh::Emission::Match ] ) {
        location = galosh::locate_value_in_sorted_array( ( group_1_euclidean_dist_to_profile[ pos_i ].begin() + 1 ), group_1_euclidean_dist_to_profile[ pos_i ].end(), group_1_euclidean_dist_to_profile[ pos_i ][ 0 ], distance_epsilon );
        lefttail = ( ( double )( location + 1 ) / num_draws );
        righttail = ( ( double )( num_draws - location ) / num_draws );
        if( lefttail < .5 ) {
          bothtails = ( lefttail * 2.0 );
        } else if( righttail < .5 ) {
          bothtails = ( righttail * 2.0 );
        } else { // == .5!
          bothtails = 1;
        }
        cout << "[" << pos_i << "] " << "Euclidean distance b/n group 1 & profile: location is " << location << endl;
        cout << "[" << pos_i << "] " << "Euclidean distance b/n group 1 & profile: p-value is " << righttail << endl;
        if( righttail < .05 ) {
          cout << "[" << pos_i << "] " << "==== SIGNIFICANT p->1 " << righttail << " *";
          if( righttail < ( .05 / profile.length() ) ) {
            cout << "!";
          }
          if( righttail < ( .05 / ( ( profile.length() * 6 ) + 6 ) ) ) {
            cout << "+";
          }
          cout << endl;
          // TODO: REMOVE
          //cout << " EPSILON: " << numeric_limits<double>::epsilon() << endl;
          cout << "[" << pos_i << "] " << "profile position is " << profile[ pos_i ][ galosh::Emission::Match ] << endl;
          cout << "[" << pos_i << "] " << "group 1 entente is " << group_1_profile_real[ pos_i ][ galosh::Emission::Match ] << endl;
          cout << "[" << pos_i << "] " << "Euclidean distance is " << group_1_euclidean_dist_to_profile[ pos_i ][ 0 ] << endl;
          cout << "Sorted permutation distances are ( " << group_1_euclidean_dist_to_profile[ pos_i ][ 1 ];
          for( uint32_t draw_i = 2; draw_i < ( num_draws + 1 ); draw_i++ ) {
            cout << ", " << group_1_euclidean_dist_to_profile[ pos_i ][ draw_i ];
          }
          cout << ")" << endl;
          cout << "[" << pos_i << "] " << "profile position is " << profile[ pos_i ][ galosh::Emission::Match ] << endl;
          cout << "[" << pos_i << "] " << "group 1 entente is " << group_1_profile_real[ pos_i ][ galosh::Emission::Match ] << endl;
        }
      } // End if the profile and the group_1_entente are different at pos_i

      if( profile[ pos_i ][ galosh::Emission::Match ] != group_2_profile_real[ pos_i ][ galosh::Emission::Match ] ) {
        location = galosh::locate_value_in_sorted_array( ( group_2_euclidean_dist_to_profile[ pos_i ].begin() + 1 ), group_2_euclidean_dist_to_profile[ pos_i ].end(), group_2_euclidean_dist_to_profile[ pos_i ][ 0 ], distance_epsilon );
        lefttail = ( ( double )( location + 1 ) / num_draws );
        righttail = ( ( double )( num_draws - location ) / num_draws );
        if( lefttail < .5 ) {
          bothtails = ( lefttail * 2.0 );
        } else if( righttail < .5 ) {
          bothtails = ( righttail * 2.0 );
        } else { // == .5!
          bothtails = 1;
        }
        cout << "[" << pos_i << "] " << "Euclidean distance b/n group 2 & profile: location is " << location << endl;
        cout << "[" << pos_i << "] " << "Euclidean distance b/n group 2 & profile: p-value is " << righttail << endl;
        if( righttail < .05 ) {
          cout << "[" << pos_i << "] " << "==== SIGNIFICANT p->2 " << righttail << " *";
          if( righttail < ( .05 / profile.length() ) ) {
            cout << "!";
          }
          if( righttail < ( .05 / ( ( profile.length() * 6 ) + 6 ) ) ) {
            cout << "+";
          }
          cout << endl;
          cout << "[" << pos_i << "] " << "profile position is " << profile[ pos_i ][ galosh::Emission::Match ] << endl;
          cout << "[" << pos_i << "] " << "group 2 entente is " << group_2_profile_real[ pos_i ][ galosh::Emission::Match ] << endl;
          cout << "[" << pos_i << "] " << "Euclidean distance is " << group_2_euclidean_dist_to_profile[ pos_i ][ 0 ] << endl;
          cout << "Sorted permutation distances are ( " << group_2_euclidean_dist_to_profile[ pos_i ][ 1 ];
          for( uint32_t draw_i = 2; draw_i < ( num_draws + 1 ); draw_i++ ) {
            cout << ", " << group_2_euclidean_dist_to_profile[ pos_i ][ draw_i ];
          }
          cout << ")" << endl;
        }
      } // End if the profile and the group_2_entente are different at pos_i

      if( group_1_profile_real[ pos_i ][ galosh::Emission::Match ] != group_2_profile_real[ pos_i ][ galosh::Emission::Match ] ) {
        location = galosh::locate_value_in_sorted_array( ( group_1_euclidean_dist_to_group_2[ pos_i ].begin() + 1 ), group_1_euclidean_dist_to_group_2[ pos_i ].end(), group_1_euclidean_dist_to_group_2[ pos_i ][ 0 ], distance_epsilon );
        lefttail = ( ( double )( location + 1 ) / num_draws );
        righttail = ( ( double )( num_draws - location ) / num_draws );
        if( lefttail < .5 ) {
          bothtails = ( lefttail * 2.0 );
        } else if( righttail < .5 ) {
          bothtails = ( righttail * 2.0 );
        } else { // == .5!
          bothtails = 1;
        }
        cout << "[" << pos_i << "] " << "Euclidean distance b/n group 1 & group 2: location is " << location << endl;
        cout << "[" << pos_i << "] " << "Euclidean distance b/n group 1 & group 2: p-value is " << righttail << endl;
        if( righttail < .05 ) {
          cout << "[" << pos_i << "] " << "==== SIGNIFICANT 1->2 " << righttail << " *";
          if( righttail < ( .05 / profile.length() ) ) {
            cout << "!";
          }
          if( righttail < ( .05 / ( ( profile.length() * 6 ) + 6 ) ) ) {
            cout << "+";
          }
          cout << endl;
          cout << "[" << pos_i << "] " << "group 1 entente is " << group_1_profile_real[ pos_i ][ galosh::Emission::Match ] << endl;
          cout << "[" << pos_i << "] " << "group 2 entente is " << group_2_profile_real[ pos_i ][ galosh::Emission::Match ] << endl;
          cout << "[" << pos_i << "] " << "Euclidean distance is " << group_1_euclidean_dist_to_group_2[ pos_i ][ 0 ] << endl;
          cout << "Sorted permutation distances are ( " << group_1_euclidean_dist_to_group_2[ pos_i ][ 1 ];
          for( uint32_t draw_i = 2; draw_i < ( num_draws + 1 ); draw_i++ ) {
            cout << ", " << group_1_euclidean_dist_to_group_2[ pos_i ][ draw_i ];
          }
          cout << ")" << endl;
        }
      } // End if the group_1_entente and the group_2_entente are different at pos_i

    } // End foreach pos_i

    location = galosh::locate_value_in_sorted_array( ( group_1_euclidean_dist_to_profile_overall.begin() + 1 ), group_1_euclidean_dist_to_profile_overall.end(), group_1_euclidean_dist_to_profile_overall[ 0 ], distance_epsilon );
    lefttail = ( ( double )( location + 1 ) / num_draws );
    righttail = ( ( double )( num_draws - location ) / num_draws );
    if( lefttail < .5 ) {
      bothtails = ( lefttail * 2.0 );
    } else if( righttail < .5 ) {
      bothtails = ( righttail * 2.0 );
    } else { // == .5!
      bothtails = 1;
    }
    cout << "[Overall] Euclidean distance b/n group 1 & profile: location is " << location << endl;
    cout << "[Overall] Euclidean distance b/n group 1 & profile: p-value is " << righttail << endl;
    if( righttail < .05 ) {
      cout << "==== SIGNIFICANT p->1 " << righttail << " *";
      if( righttail < ( .05 / ( ( profile.length() * 6 ) + 6 ) ) ) {
        cout << "+";
      }
      cout << endl;
    }

    location = galosh::locate_value_in_sorted_array( ( group_2_euclidean_dist_to_profile_overall.begin() + 1 ), group_2_euclidean_dist_to_profile_overall.end(), group_2_euclidean_dist_to_profile_overall[ 0 ], distance_epsilon );
    lefttail = ( ( double )( location + 1 ) / num_draws );
    righttail = ( ( double )( num_draws - location ) / num_draws );
    if( lefttail < .5 ) {
      bothtails = ( lefttail * 2.0 );
    } else if( righttail < .5 ) {
      bothtails = ( righttail * 2.0 );
    } else { // == .5!
      bothtails = 1;
    }
    cout << "[Overall] Euclidean distance b/n group 2 & profile: location is " << location << endl;
    cout << "[Overall] Euclidean distance b/n group 2 & profile: p-value is " << righttail << endl;
    if( righttail < .05 ) {
      cout << "==== SIGNIFICANT p->2 " << righttail << " *";
      if( righttail < ( .05 / ( ( profile.length() * 6 ) + 6 ) ) ) {
        cout << "+";
      }
      cout << endl;
    }

    location = galosh::locate_value_in_sorted_array( ( group_1_euclidean_dist_to_group_2_overall.begin() + 1 ), group_1_euclidean_dist_to_group_2_overall.end(), group_1_euclidean_dist_to_group_2_overall[ 0 ], distance_epsilon );
    lefttail = ( ( double )( location + 1 ) / num_draws );
    righttail = ( ( double )( num_draws - location ) / num_draws );
    if( lefttail < .5 ) {
      bothtails = ( lefttail * 2.0 );
    } else if( righttail < .5 ) {
      bothtails = ( righttail * 2.0 );
    } else { // == .5!
      bothtails = 1;
    }
    cout << "[Overall] Euclidean distance b/n group 1 & group 2: location is " << location << endl;
    cout << "[Overall] Euclidean distance b/n group 1 & group 2: p-value is " << righttail << endl;
    if( righttail < .05 ) {
      cout << "==== SIGNIFICANT 1->2 " << righttail << " *";
      if( righttail < ( .05 / ( ( profile.length() * 6 ) + 6 ) ) ) {
        cout << "+";
      }
      cout << endl;
    }

    // Skl distance

    for( uint32_t pos_i = 0; pos_i < profile.length(); pos_i++ ) {
      location = galosh::locate_value_in_sorted_array( ( group_1_skl_dist_to_profile[ pos_i ].begin() + 1 ), group_1_skl_dist_to_profile[ pos_i ].end(), group_1_skl_dist_to_profile[ pos_i ][ 0 ], distance_epsilon );
      lefttail = ( ( double )( location + 1 ) / num_draws );
      righttail = ( ( double )( num_draws - location ) / num_draws );
      if( lefttail < .5 ) {
        bothtails = ( lefttail * 2.0 );
      } else if( righttail < .5 ) {
        bothtails = ( righttail * 2.0 );
      } else { // == .5!
        bothtails = 1;
      }
      cout << "[" << pos_i << "] " << "Skl distance b/n group 1 & profile: location is " << location << endl;
      cout << "[" << pos_i << "] " << "Skl distance b/n group 1 & profile: p-value is " << righttail << endl;
      if( righttail < .05 ) {
        cout << "[" << pos_i << "] " << "==== SIGNIFICANT p->1 " << righttail << " *";
        if( righttail < ( .05 / profile.length() ) ) {
          cout << "!";
        }
        if( righttail < ( .05 / ( ( profile.length() * 6 ) + 6 ) ) ) {
          cout << "+";
        }
        cout << endl;
        cout << "[" << pos_i << "] " << "profile position is " << profile[ pos_i ][ galosh::Emission::Match ] << endl;
        cout << "[" << pos_i << "] " << "group 1 entente is " << group_1_profile_real[ pos_i ][ galosh::Emission::Match ] << endl;
      }

      location = galosh::locate_value_in_sorted_array( ( group_2_skl_dist_to_profile[ pos_i ].begin() + 1 ), group_2_skl_dist_to_profile[ pos_i ].end(), group_2_skl_dist_to_profile[ pos_i ][ 0 ], distance_epsilon );
      lefttail = ( ( double )( location + 1 ) / num_draws );
      righttail = ( ( double )( num_draws - location ) / num_draws );
      if( lefttail < .5 ) {
        bothtails = ( lefttail * 2.0 );
      } else if( righttail < .5 ) {
        bothtails = ( righttail * 2.0 );
      } else { // == .5!
        bothtails = 1;
      }
      cout << "[" << pos_i << "] " << "Skl distance b/n group 2 & profile: location is " << location << endl;
      cout << "[" << pos_i << "] " << "Skl distance b/n group 2 & profile: p-value is " << righttail << endl;
      if( righttail < .05 ) {
        cout << "[" << pos_i << "] " << "==== SIGNIFICANT p->2 " << righttail << " *";
        if( righttail < ( .05 / profile.length() ) ) {
          cout << "!";
        }
        if( righttail < ( .05 / ( ( profile.length() * 6 ) + 6 ) ) ) {
          cout << "+";
        }
        cout << endl;
        cout << "[" << pos_i << "] " << "profile position is " << profile[ pos_i ][ galosh::Emission::Match ] << endl;
        cout << "[" << pos_i << "] " << "group 2 entente is " << group_2_profile_real[ pos_i ][ galosh::Emission::Match ] << endl;
      }
  
      location = galosh::locate_value_in_sorted_array( ( group_1_skl_dist_to_group_2[ pos_i ].begin() + 1 ), group_1_skl_dist_to_group_2[ pos_i ].end(), group_1_skl_dist_to_group_2[ pos_i ][ 0 ], distance_epsilon );
      lefttail = ( ( double )( location + 1 ) / num_draws );
      righttail = ( ( double )( num_draws - location ) / num_draws );
      if( lefttail < .5 ) {
        bothtails = ( lefttail * 2.0 );
      } else if( righttail < .5 ) {
        bothtails = ( righttail * 2.0 );
      } else { // == .5!
        bothtails = 1;
      }
      cout << "[" << pos_i << "] " << "Skl distance b/n group 1 & group 2: location is " << location << endl;
      cout << "[" << pos_i << "] " << "Skl distance b/n group 1 & group 2: p-value is " << righttail << endl;
      if( righttail < .05 ) {
        cout << "[" << pos_i << "] " << "==== SIGNIFICANT 1->2 " << righttail << " *";
        if( righttail < ( .05 / profile.length() ) ) {
          cout << "!";
        }
        if( righttail < ( .05 / ( ( profile.length() * 6 ) + 6 ) ) ) {
          cout << "+";
        }
        cout << endl;
        cout << "[" << pos_i << "] " << "group 1 entente is " << group_1_profile_real[ pos_i ][ galosh::Emission::Match ] << endl;
        cout << "[" << pos_i << "] " << "group 2 entente is " << group_2_profile_real[ pos_i ][ galosh::Emission::Match ] << endl;
      }
    } // End foreach pos_i

    location = galosh::locate_value_in_sorted_array( ( group_1_skl_dist_to_profile_overall.begin() + 1 ), group_1_skl_dist_to_profile_overall.end(), group_1_skl_dist_to_profile_overall[ 0 ], distance_epsilon );
    lefttail = ( ( double )( location + 1 ) / num_draws );
    righttail = ( ( double )( num_draws - location ) / num_draws );
    if( lefttail < .5 ) {
      bothtails = ( lefttail * 2.0 );
    } else if( righttail < .5 ) {
      bothtails = ( righttail * 2.0 );
    } else { // == .5!
      bothtails = 1;
    }
    cout << "[Overall] Skl distance b/n group 1 & profile: location is " << location << endl;
    cout << "[Overall] Skl distance b/n group 1 & profile: p-value is " << righttail << endl;
    if( righttail < .05 ) {
      cout << "==== SIGNIFICANT p->1 " << righttail << " *";
      if( righttail < ( .05 / ( ( profile.length() * 6 ) + 6 ) ) ) {
        cout << "+";
      }
      cout << endl;
    }

    location = galosh::locate_value_in_sorted_array( ( group_2_skl_dist_to_profile_overall.begin() + 1 ), group_2_skl_dist_to_profile_overall.end(), group_2_skl_dist_to_profile_overall[ 0 ], distance_epsilon );
    lefttail = ( ( double )( location + 1 ) / num_draws );
    righttail = ( ( double )( num_draws - location ) / num_draws );
    if( lefttail < .5 ) {
      bothtails = ( lefttail * 2.0 );
    } else if( righttail < .5 ) {
      bothtails = ( righttail * 2.0 );
    } else { // == .5!
      bothtails = 1;
    }
    cout << "[Overall] Skl distance b/n group 2 & profile: location is " << location << endl;
    cout << "[Overall] Skl distance b/n group 2 & profile: p-value is " << righttail << endl;
    if( righttail < .05 ) {
      cout << "==== SIGNIFICANT p->2 " << righttail << " *";
      if( righttail < ( .05 / ( ( profile.length() * 6 ) + 6 ) ) ) {
        cout << "+";
      }
      cout << endl;
    }

    location = galosh::locate_value_in_sorted_array( ( group_1_skl_dist_to_group_2_overall.begin() + 1 ), group_1_skl_dist_to_group_2_overall.end(), group_1_skl_dist_to_group_2_overall[ 0 ], distance_epsilon );
    lefttail = ( ( double )( location + 1 ) / num_draws );
    righttail = ( ( double )( num_draws - location ) / num_draws );
    if( lefttail < .5 ) {
      bothtails = ( lefttail * 2.0 );
    } else if( righttail < .5 ) {
      bothtails = ( righttail * 2.0 );
    } else { // == .5!
      bothtails = 1;
    }
    cout << "[Overall] Skl distance b/n group 1 & group 2: location is " << location << endl;
    cout << "[Overall] Skl distance b/n group 1 & group 2: p-value is " << righttail << endl;
    if( righttail < .05 ) {
      cout << "==== SIGNIFICANT 1->2 " << righttail << " *";
      if( righttail < ( .05 / ( ( profile.length() * 6 ) + 6 ) ) ) {
        cout << "+";
      }
      cout << endl;
    }



    if( 0 ) {
      for( uint32_t draw_i = 0; draw_i < ( num_draws + 1 ); draw_i++ ) {
        cout << group_1_skl_dist_to_group_2_overall[ draw_i ] << endl;
      }
      cout << "Location of our value is " << location << endl;
    }

    //
    //for( uint32_t draw_i = 0; draw_i < ( num_draws + 1 ); draw_i++ ) {
    //  cout << group_1_skl_dist_to_group_2[ 0 ][ draw_i ] << endl;
    //}

    // TODO: ...?
    exit( 0 );
  } // End if do_explorations && num_seqs_in_group_1


  if( true ) {
    // Calculate a distance matrix based on the alignment profiles.  Note
    // that the alignment profiles may be weighted and/or normalized first
    // (altering them), depending on the method.  Also the matrix might be a
    // proximity matrix (test with dist_or_prox_matrix.isProximityMatrix()).
    galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::DistanceMatrix<doublerealspace> dist_or_prox_matrix;
        // Cooccurrence expected count
        // The second argument sez count cooccurrence of deletions, too.
    //dist_or_prox_matrix.setToMatchEmissionCooccurrenceExpectedCounts( aps, true );
    // TODO: This class should accept the metric as an argument
        //tree_trainer.calculateDistanceMatrix( 0, aps, dist_or_prox_matrix, galosh::Tag<galosh::COOCExpectedCountMetric_>() );
    //calculateDistanceMatrix( 0, aps, dist_or_prox_matrix, galosh::Tag<galosh::COOCProbabilityMetric_>() );
    //calculateDistanceMatrix( 0, aps, dist_or_prox_matrix, galosh::Tag<galosh::CrossEntropyMetric_>() );
    //calculateDistanceMatrix( 0, aps, dist_or_prox_matrix, galosh::Tag<galosh::SKLMetric_>() );
  
    // Eucliean distances
    // TODO: PUT BACK
    //dist_or_prox_matrix.setToEuclideanDistances( aps );
    //cout << "Euclidean distances among the APs:" << endl;

    for( uint32_t ap_pos_i = 1; ap_pos_i < aps[ 0 ].length(); ap_pos_i++ ) {
    dist_or_prox_matrix.setToEuclideanDistances( aps, ap_pos_i );
    cout << "Euclidean distances among the APs AT POS " << ap_pos_i << endl;
   
    //cout << "[" << ap_pos_i << "] " << "Euclidean Distance Matrix is " << dist_or_prox_matrix << endl;
    cout << "[" << ap_pos_i << "] " << "Max dist: " << dist_or_prox_matrix.maximumOffDiagonalValue() << endl;
    cout << "[" << ap_pos_i << "] " << "Min dist: " << dist_or_prox_matrix.minimumOffDiagonalValue() << endl;
    doublerealspace avg = toDouble( dist_or_prox_matrix.averageOffDiagonal() );
    cout << "[" << ap_pos_i << "] " << "Avg dist: " << avg << endl;
    doublerealspace var = dist_or_prox_matrix.varianceOffDiagonal( avg );
    cout << "[" << ap_pos_i << "] " << "SD of dist: " << pow( var, .5 ) << endl;
    //cout << "[" << ap_pos_i << "] " << "SE of Avg: " << pow( ( var / ( aps.size() * ( aps.size() - 1 ) ) ), .5 ) << endl;
    cout << "[" << ap_pos_i << "] " << "SE of Avg: " << pow( ( var / aps.size() ), .5 ) << endl;
  
    if( num_seqs_in_group_1 && ( num_seqs_in_group_1 != aps.size() ) ) {
      // Compare within and between two groups defined by the first
      // num_seqs_in_group_1 seqs, and the rest.
  
      uint32_t const num_seqs_in_group_2 = ( aps.size() - num_seqs_in_group_1 );
  
      doublerealspace group_1_avg;
      doublerealspace group_1_var;
      doublerealspace group_1_se;
      doublerealspace group_2_avg;
      doublerealspace group_2_var;
      doublerealspace group_2_se;
  
      // TODO: REMOVE.  TESTING
      alglib::real_1d_array group_1_alglib_array;
      // TODO: REMOVE.  Set length.
      group_1_alglib_array.setlength( static_cast<alglib::ae_int_t>( ( num_seqs_in_group_1 * ( num_seqs_in_group_1 - 1 ) ) / 2 ) );
      // TODO: REMOVE.  TESTING
      alglib::real_1d_array group_2_alglib_array;
      // TODO: REMOVE.  Set length.
      group_2_alglib_array.setlength( static_cast<alglib::ae_int_t>( ( num_seqs_in_group_2 * ( num_seqs_in_group_2 - 1 ) ) / 2 ) );
  
      if( num_seqs_in_group_1 > 1 ) {
        // TODO: REMOVE.  TESTING
        alglib::ae_int_t group_1_alglib_array_i = 0;
        for( uint32_t from = 0; from < num_seqs_in_group_1; from++ ) {
          for( uint32_t to = 0; to < from; to++ ) {
            group_1_alglib_array[ group_1_alglib_array_i ] = toDouble( dist_or_prox_matrix( from, to ) );
            ++group_1_alglib_array_i;
          } // End foreach to..
        } // End foreach from..

        cout << "[" << ap_pos_i << "] " << "Euclidean distances among the APs within the " << num_seqs_in_group_1 << " sequences in group 1:" << endl;
        cout << "[" << ap_pos_i << "] " << "Max dist: " << dist_or_prox_matrix.maximumOffDiagonalValue( 0, num_seqs_in_group_1 ) << endl;
        cout << "[" << ap_pos_i << "] " << "Min dist: " << dist_or_prox_matrix.minimumOffDiagonalValue( 0, num_seqs_in_group_1 ) << endl;
        group_1_avg = dist_or_prox_matrix.averageOffDiagonal( 0, num_seqs_in_group_1 );
        cout << "[" << ap_pos_i << "] " << "Avg dist: " << group_1_avg << endl;
        group_1_var = dist_or_prox_matrix.varianceOffDiagonal( 0, num_seqs_in_group_1, group_1_avg );
        cout << "[" << ap_pos_i << "] " << "SD of dist: " << pow( group_1_var, .5 ) << endl;
        cout << "[" << ap_pos_i << "] " << "SE of Avg: " << pow( ( group_1_var / ( num_seqs_in_group_1 * ( num_seqs_in_group_1 - 1 ) ) ), .5 ) << endl;
      } // End if( num_seqs_in_group_1 > 1 )

      if( num_seqs_in_group_2 > 1 ) {
        // TODO: REMOVE.  TESTING
        alglib::ae_int_t group_2_alglib_array_i = 0;
        for( uint32_t from = num_seqs_in_group_1; from < aps.size(); from++ ) {
          for( uint32_t to = num_seqs_in_group_1; to < from; to++ ) {
            group_2_alglib_array[ group_2_alglib_array_i ] = toDouble( dist_or_prox_matrix( from, to ) );
            ++group_2_alglib_array_i;
          } // End foreach to..
        } // End foreach from..
  
        cout << "[" << ap_pos_i << "] " << "Euclidean distances among the APs within the " << num_seqs_in_group_2 << " sequences in group 2:" << endl;
        cout << "[" << ap_pos_i << "] " << "Max dist: " << dist_or_prox_matrix.maximumOffDiagonalValue( num_seqs_in_group_1, num_seqs_in_group_2 ) << endl;
        cout << "[" << ap_pos_i << "] " << "Min dist: " << dist_or_prox_matrix.minimumOffDiagonalValue( num_seqs_in_group_1, num_seqs_in_group_2 ) << endl;
        group_2_avg = dist_or_prox_matrix.averageOffDiagonal( num_seqs_in_group_1, num_seqs_in_group_2 );
        cout << "[" << ap_pos_i << "] " << "Avg dist: " << group_2_avg << endl;
        group_2_var = dist_or_prox_matrix.varianceOffDiagonal( num_seqs_in_group_1, num_seqs_in_group_2, group_2_avg );
        cout << "[" << ap_pos_i << "] " << "SD of dist: " << pow( group_2_var, .5 ) << endl;
        //cout << "[" << ap_pos_i << "] " << "SE of Avg: " << pow( ( group_2_var / ( num_seqs_in_group_2 * ( num_seqs_in_group_2 - 1 ) ) ), .5 ) << endl;
        cout << "[" << ap_pos_i << "] " << "SE of Avg: " << pow( ( group_2_var / num_seqs_in_group_2 ), .5 ) << endl;
      } // End if( num_seqs_in_group_2 > 1 )

      // Now we need to compute those values *across* groups.
      doublerealspace across_max = dist_or_prox_matrix( 0, num_seqs_in_group_1 );
      doublerealspace across_min = dist_or_prox_matrix( 0, num_seqs_in_group_1 );
      doublerealspace across_avg = 0;
      doublerealspace val;
      uint32_t tmp_count = 0;
      for( uint32_t from = 0; from < num_seqs_in_group_1; from++ ) {
        for( uint32_t to = num_seqs_in_group_1; to < aps.size(); to++ ) {
          val = dist_or_prox_matrix( from, to );
          if( val > across_max ) {
            across_max = val;
          }
          if( val < across_min ) {
            across_min = val;
          }
          across_avg += val;
          ++tmp_count;
        } // End foreach to
      } // End foreach from
      across_avg /= tmp_count;

      // Now for the variance (we needed to first compute the average)
      doublerealspace across_var = 0;
      // TODO: REMOVE.  TESTING
      alglib::real_1d_array across_alglib_array;
      // TODO: REMOVE.  Set length.
      across_alglib_array.setlength( static_cast<alglib::ae_int_t>( tmp_count ) );
      tmp_count = 0;
      for( uint32_t from = 0; from < num_seqs_in_group_1; from++ ) {
        for( uint32_t to = num_seqs_in_group_1; to < aps.size(); to++ ) {
          val = dist_or_prox_matrix( from, to );
          // TODO: REMOVE
          across_alglib_array[ static_cast<alglib::ae_int_t>( tmp_count ) ] = toDouble( val );
          // TODO: REMOVE
          //if( val == 0 ) {
          //  cout << "dist is 0?! from " << from << " to " << to << endl;
          //  cout << "\tap[ " << from << " ]: " << aps[ from ][ ap_pos_i ][ galosh::Emission::Match ] << endl;
          //  cout << "\tap[ " << to << " ]: " << aps[ to ][ ap_pos_i ][ galosh::Emission::Match ] << endl;
          //  cout << "\tSequence " << from << " is:\n\t\t" << fasta[ from ] << endl;
          //  cout << "\tSequence " << to << " is:\n\t\t" << fasta[ to ] << endl;
          //} // End if ( val == 0 )
          val -= across_avg;
          across_var += ( val * val );
          ++tmp_count;
        } // End foreach to
      } // End foreach from
      across_var /= ( tmp_count - 1 ); // -1 since we're estimating the mean
      //across_var /= ( aps.size() - 1 ); // -1 since we're estimating the mean -- here I'm using the total number of seqs as the df for estimating the var, since that makes more sense to me (adding one more seq shouldn't double the df!).  TODO: Determine truly correct df.

      cout << "[" << ap_pos_i << "] " << "Euclidean distances among the APs across / between the two groups:" << endl;
      cout << "[" << ap_pos_i << "] " << "Max dist: " << across_max << endl;
      cout << "[" << ap_pos_i << "] " << "Min dist: " << across_min << endl;
      cout << "[" << ap_pos_i << "] " << "Avg dist: " << across_avg << endl;
      cout << "[" << ap_pos_i << "] " << "SD of dist: " << pow( across_var, .5 ) << endl;
      //cout << "[" << ap_pos_i << "] " << "SE of Avg: " << pow( ( across_var / tmp_count ), .5 ) << endl;
      cout << "[" << ap_pos_i << "] " << "SE of Avg: " << pow( ( across_var / aps.size() ), .5 ) << endl;

      // TODO: ERE I AM.  Calc stats?  t-test?  What to what?
      // TODO: REMOVE.  Trying alglib.
      double bothtails, lefttail, righttail;
      if( num_seqs_in_group_1 > 1 ) {
        //unequalvariancettest( across_alglib_array, across_alglib_array.length(), group_1_alglib_array, group_1_alglib_array.length(), bothtails, lefttail, righttail );
        //cout << "[" << ap_pos_i << "] " << "two-sample unequal variance t-test across vs group 1: p-value is " << bothtails << endl;
        mannwhitneyutest( across_alglib_array, across_alglib_array.length(), group_1_alglib_array, group_1_alglib_array.length(), bothtails, lefttail, righttail );
        cout << "[" << ap_pos_i << "] " << "u-test across vs group 1: p-value is " << bothtails << endl;
        if( bothtails < .05 ) {
          cout << "[" << ap_pos_i << "] " << "==== SIGNIFICANT b->1 " << bothtails << " *";
          if( bothtails < ( .05 / ( aps.size() * 2 ) ) ) {
            cout << "!";
          }
          cout << endl;
        }
      }
      if( num_seqs_in_group_2 > 1 ) {
        //unequalvariancettest( across_alglib_array, across_alglib_array.length(), group_2_alglib_array, group_2_alglib_array.length(), bothtails, lefttail, righttail );
        //cout << "[" << ap_pos_i << "] " << "two-sample unequal variance t-test across vs group 2: p-value is " << bothtails << endl;
        mannwhitneyutest( across_alglib_array, across_alglib_array.length(), group_2_alglib_array, group_2_alglib_array.length(), bothtails, lefttail, righttail );
        cout << "[" << ap_pos_i << "] " << "u-test across vs group 2: p-value is " << bothtails << endl;
        if( bothtails < .05 ) {
          cout << "[" << ap_pos_i << "] " << "==== SIGNIFICANT b->2 " << bothtails << " *";
          if( bothtails < ( .05 / ( aps.size() * 2 ) ) ) {
            cout << "!";
          }
          cout << endl;
        }
      }

    } // End if( num_seqs_in_group_1 )
    } // End foreach ap_pos_i ..

    //exit( 0 );
  }

  // Since before normalization, the numbers will be very very tiny, we pass values temporarily through this one (which stores them in MatrixValueType).
  galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile tmp_ap_mvt( 1 ); // TODO: Just use an AlignmentProfilePosition object?

  // By residue, foreach sequence and position:
  // Use doubles, so we can handle negative values. ..
  vector<galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, doublerealspace>::AlignmentProfile> ap_deviations_from_profile( aps.size() );

  for( uint32_t ap_i = 0; ap_i < aps.size(); ap_i++ ) {
    ap_deviations_from_profile[ ap_i ].reinitialize( profile.length() + 1 );
    // TODO: Can we include P(Del) here too?  It's hard to figure out what the profile mean is.  TODO: Estimate it from the seqs?
    // MatrixValueType pos_del_prob;

    // The APs have one more pos than does the profile (just like the dp matrix does) -- so the first entry doesn't correspond to any match emission distribution.
    //for( uint32_t ap_pos_i = 0; ap_pos_i < aps[ ap_i ].length(); ap_pos_i++ ) {
    for( uint32_t ap_pos_i = 1; ap_pos_i < aps[ ap_i ].length(); ap_pos_i++ ) {
      tmp_ap_mvt[ 0 ] = aps[ ap_i ][ ap_pos_i ];
      //pos_del_prob = ( 1.0 - tmp_ap_mvt[ 0 ][ galosh::Emission::Match ].total() );
      //pos_del_prob -= ???
      //ap_deviations_from_profile[ ap_i ][ ap_pos_i ][ galosh::Emission::Match ].normalize( 0 );
      //ap_deviations_from_profile[ ap_i ][ ap_pos_i ][ galosh::Emission::Match ] -= profile[ ap_pos_i ][ galosh::Emission::Match ];

      // Can we do it with all of the params?
      //tmp_ap_mvt[ 0 ].normalize( 0 );
      tmp_ap_mvt[ 0 ][ galosh::Emission::Match ].normalize( 0 );
      // TODO: Make this work:
      //ap_deviations_from_profile[ ap_i ][ ap_pos_i ] = tmp_ap_mvt[ 0 ];
      ap_deviations_from_profile[ ap_i ][ ap_pos_i ][ galosh::Emission::Match ] = tmp_ap_mvt[ 0 ][ galosh::Emission::Match ];
      // Now do the subtraction in doublerealspace type
      // TODO: Make this work:
      //ap_deviations_from_profile[ ap_i ][ ap_pos_i ] -= profile[ ap_pos_i - 1 ];
      ap_deviations_from_profile[ ap_i ][ ap_pos_i ][ galosh::Emission::Match ] -= profile[ ap_pos_i - 1 ][ galosh::Emission::Match ];

      // TODO: ERE I AM
      //cout << "[" << ap_i << "," << ap_pos_i << "] " << aps[ ap_i ][ ap_pos_i ][ galosh::Emission::Match ] << " => " << ap_deviations_from_profile[ ap_i ][ ap_pos_i ][ galosh::Emission::Match ] << endl;
      //cout << "[" << ap_i << "," << ap_pos_i << "] " << tmp_ap_mvt[ 0 ][ galosh::Emission::Match ] << " => " << ap_deviations_from_profile[ ap_i ][ ap_pos_i ][ galosh::Emission::Match ] << endl;
    } // End foreach position ap_pos_i...
  } // End foreach alignment profile ap_i...

  // Calculate a distance matrix based on the alignment profiles.  Note
  // that the alignment profiles may be weighted and/or normalized first
  // (altering them), depending on the method.  Also the matrix might be a
  // proximity matrix (test with dist_or_prox_matrix.isProximityMatrix()).
  galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, doublerealspace>::DistanceMatrix<doublerealspace> dist_or_prox_matrix;
  // Cooccurrence expected count
  //dist_or_prox_matrix.setToMatchEmissionCooccurrenceExpectedCounts( ap_deviations_from_profile, false );
  //cout << "Match emission Cooc EC Matrix is " << dist_or_prox_matrix << endl;
  // SKL distances
  //dist_or_prox_matrix.setToSymmeterizedKullbackLeiblerDivergences( ap_deviations_from_profile );
  //cout << "Match emission SKL Distance Matrix is " << dist_or_prox_matrix << endl;

  // Eucliean distances
  cout << "Euclidean distances for AP deviations from the profile values: " << endl;
  dist_or_prox_matrix.setToEuclideanDistances( ap_deviations_from_profile );
  //cout << "Euclidean Distance Matrix is " << dist_or_prox_matrix << endl;
  cout << "Max dist: " << dist_or_prox_matrix.maximumOffDiagonalValue() << endl;
  cout << "Min dist: " << dist_or_prox_matrix.minimumOffDiagonalValue() << endl;
  doublerealspace avg = dist_or_prox_matrix.averageOffDiagonal();
  cout << "Avg dist: " << avg << endl;
  doublerealspace var = dist_or_prox_matrix.varianceOffDiagonal( avg );
  cout << "SD of dist: " << pow( var, .5 ) << endl;
  // TODO: ERE I AM.  Shouldn't we just use aps.size(), since that is actually the number of data points we have?  It seems like there can't be more dfs than aps.size() - 1.
  //cout << "SE of Avg: " << pow( ( var / ( aps.size() * ( aps.size() - 1 ) ) ), .5 ) << endl;
  cout << "SE of Avg: " << pow( ( var / aps.size() ), .5 ) << endl;

  if( num_seqs_in_group_1 && ( num_seqs_in_group_1 != aps.size() ) ) {
    // Compare within and between two groups defined by the first
    // num_seqs_in_group_1 seqs, and the rest.

    uint32_t const num_seqs_in_group_2 = ( aps.size() - num_seqs_in_group_1 );

    doublerealspace group_1_avg;
    doublerealspace group_1_var;
    doublerealspace group_1_se;
    doublerealspace group_2_avg;
    doublerealspace group_2_var;
    doublerealspace group_2_se;

    // TODO: REMOVE.  TESTING
    alglib::real_1d_array group_1_alglib_array;
    // TODO: REMOVE.  Set length.
    group_1_alglib_array.setlength( static_cast<alglib::ae_int_t>( ( num_seqs_in_group_1 * ( num_seqs_in_group_1 - 1 ) ) / 2 ) );
    // TODO: REMOVE.  TESTING
    alglib::real_1d_array group_2_alglib_array;
    // TODO: REMOVE.  Set length.
    group_2_alglib_array.setlength( static_cast<alglib::ae_int_t>( ( num_seqs_in_group_2 * ( num_seqs_in_group_2 - 1 ) ) / 2 ) );

    if( num_seqs_in_group_1 > 1 ) {
      // TODO: REMOVE.  TESTING
      alglib::ae_int_t group_1_alglib_array_i = 0;
      for( uint32_t from = 0; from < num_seqs_in_group_1; from++ ) {
        for( uint32_t to = 0; to < from; to++ ) {
          group_1_alglib_array[ group_1_alglib_array_i ] = toDouble( dist_or_prox_matrix( from, to ) );
          ++group_1_alglib_array_i;
        } // End foreach to..
      } // End foreach from..

      cout << "Euclidean distances among the APs within the " << num_seqs_in_group_1 << " sequences in group 1:" << endl;
      cout << "Max dist: " << dist_or_prox_matrix.maximumOffDiagonalValue( 0, num_seqs_in_group_1 ) << endl;
      cout << "Min dist: " << dist_or_prox_matrix.minimumOffDiagonalValue( 0, num_seqs_in_group_1 ) << endl;
      group_1_avg = dist_or_prox_matrix.averageOffDiagonal( 0, num_seqs_in_group_1 );
      cout << "Avg dist: " << group_1_avg << endl;
      group_1_var = dist_or_prox_matrix.varianceOffDiagonal(  0, num_seqs_in_group_1, group_1_avg );
      cout << "SD of dist: " << pow( group_1_var, .5 ) << endl;
      //group_1_se = pow( ( group_1_var / ( num_seqs_in_group_1 * ( num_seqs_in_group_1 - 1 ) ) ), .5 );
      group_1_se = pow( ( group_1_var / num_seqs_in_group_1 ), .5 );
      cout << "SE of Avg: " << group_1_se << endl;
    } // End if( num_seqs_in_group_1 > 1 )
    
    if( num_seqs_in_group_2 > 1 ) {
      // TODO: REMOVE.  TESTING
      alglib::ae_int_t group_2_alglib_array_i = 0;
      for( uint32_t from = num_seqs_in_group_1; from < aps.size(); from++ ) {
        for( uint32_t to = num_seqs_in_group_1; to < from; to++ ) {
          group_2_alglib_array[ group_2_alglib_array_i ] = toDouble( dist_or_prox_matrix( from, to ) );
          ++group_2_alglib_array_i;
        } // End foreach to..
      } // End foreach from..

      cout << "Euclidean distances among the APs within the " << num_seqs_in_group_2 << " sequences in group 2:" << endl;
      cout << "Max dist: " << dist_or_prox_matrix.maximumOffDiagonalValue( num_seqs_in_group_1, num_seqs_in_group_2 ) << endl;
      cout << "Min dist: " << dist_or_prox_matrix.minimumOffDiagonalValue( num_seqs_in_group_1, num_seqs_in_group_2 ) << endl;
      group_2_avg = dist_or_prox_matrix.averageOffDiagonal( num_seqs_in_group_1, num_seqs_in_group_2 );
      cout << "Avg dist: " << group_2_avg << endl;
      group_2_var = dist_or_prox_matrix.varianceOffDiagonal(  num_seqs_in_group_1, num_seqs_in_group_2, group_2_avg );
      cout << "SD of dist: " << pow( group_2_var, .5 ) << endl;
      //group_2_se = pow( ( group_2_var / ( num_seqs_in_group_2 * ( num_seqs_in_group_2 - 1 ) ) ), .5 );
      group_2_se = pow( ( group_2_var / num_seqs_in_group_2 ), .5 );
      cout << "SE of Avg: " << group_2_se << endl;
    } // End if( num_seqs_in_group_2 > 1 )
    
    // Now we need to compute those values *across* groups.
    doublerealspace across_max = dist_or_prox_matrix( 0, num_seqs_in_group_1 );
    doublerealspace across_min = dist_or_prox_matrix( 0, num_seqs_in_group_1 );
    doublerealspace across_avg = 0;
    doublerealspace val;
    uint32_t tmp_count = 0;
    for( uint32_t from = 0; from < num_seqs_in_group_1; from++ ) {
      for( uint32_t to = num_seqs_in_group_1; to < aps.size(); to++ ) {
        val = dist_or_prox_matrix( from, to );
        if( val > across_max ) {
          across_max = val;
        }
        if( val < across_min ) {
          across_min = val;
        }
        across_avg += val;
        ++tmp_count;
      } // End foreach to
    } // End foreach from
    across_avg /= tmp_count;
    
    // Now for the variance (we needed to first compute the average)
    doublerealspace across_var = 0;
    // TODO: REMOVE.  TESTING
    alglib::real_1d_array across_alglib_array;
    // TODO: REMOVE.  Set length.  Is this necessary? ?
    across_alglib_array.setlength( static_cast<alglib::ae_int_t>( tmp_count ) );
    tmp_count = 0;
    for( uint32_t from = 0; from < num_seqs_in_group_1; from++ ) {
      for( uint32_t to = num_seqs_in_group_1; to < aps.size(); to++ ) {
        val = dist_or_prox_matrix( from, to );
        // TODO: REMOVE
        across_alglib_array[ static_cast<alglib::ae_int_t>( tmp_count ) ] = toDouble( val );
        // TODO: REMOVE
        //if( val == 0 ) {
        //  cout << "dist is 0?! from " << from << " to " << to << endl;
        //  for( uint32_t ap_pos_i = 0; ap_pos_i < aps[ 0 ].length(); ap_pos_i++ ) {
        //    cout << "AP POS " << ap_pos_i << ":" << endl;
        //    cout << "\tap[ " << from << " ]: " << aps[ from ][ ap_pos_i ][ galosh::Emission::Match ] << endl;
        //    cout << "\tap[ " << to << " ]: " << aps[ to ][ ap_pos_i ][ galosh::Emission::Match ] << endl;
        //  } // End forach ap_pos_i
        //}
        val -= across_avg;
        across_var += ( val * val );
        ++tmp_count;
      } // End foreach to
    } // End foreach from
    across_var /= ( tmp_count - 1 ); // -1 since we estimated the mean
    
    cout << "Euclidean distances among the APs across / between the two groups:" << endl;
    cout << "Max dist: " << across_max << endl;
    cout << "Min dist: " << across_min << endl;
    cout << "Avg dist: " << across_avg << endl;
    cout << "SD of dist: " << pow( across_var, .5 ) << endl;
    doublerealspace across_se = pow( ( across_var / tmp_count ), .5 );
    cout << "SE of Avg: " << across_se << endl;

    // TODO: ERE I AM.  Calc stats?  t-test?  What to what?
    // TODO: REMOVE.  Trying alglib.
    double bothtails, lefttail, righttail;
    // TODO: REMOVE
    cout << "across distances: ( " << across_alglib_array[ 0 ];
    for( alglib::ae_int_t across_alglib_array_i = 1; across_alglib_array_i < across_alglib_array.length(); across_alglib_array_i++ ) {
      cout << ", " << across_alglib_array[ across_alglib_array_i ];
    }
    cout << " )" << endl;
    if( num_seqs_in_group_1 > 1 ) {
      // TODO: REMOVE
      cout << "group 1 distances: ( " << group_1_alglib_array[ 0 ];
      for( alglib::ae_int_t group_1_alglib_array_i = 1; group_1_alglib_array_i < group_1_alglib_array.length(); group_1_alglib_array_i++ ) {
        cout << ", " << group_1_alglib_array[ group_1_alglib_array_i ];
      }
      cout << " )" << endl;
      //unequalvariancettest( across_alglib_array, across_alglib_array.length(), group_1_alglib_array, group_1_alglib_array.length(), bothtails, lefttail, righttail );
      //cout << "two-sample unequal variance t-test across vs group 1: p-value is " << bothtails << endl;
      mannwhitneyutest( across_alglib_array, across_alglib_array.length(), group_1_alglib_array, group_1_alglib_array.length(), bothtails, lefttail, righttail );
      cout << "u-test across vs group 1: p-value is " << bothtails << endl;
      if( bothtails < .05 ) {
        cout << "==== SIGNIFICANT b->1 " << bothtails << " *";
        if( bothtails < ( .05 / 2 ) ) {
          cout << "!";
        }
        cout << endl;
      }
    }
    if( num_seqs_in_group_2 > 1 ) {
      // TODO: REMOVE
      cout << "group 2 distances: ( " << group_2_alglib_array[ 0 ];
      for( alglib::ae_int_t group_2_alglib_array_i = 1; group_2_alglib_array_i < group_2_alglib_array.length(); group_2_alglib_array_i++ ) {
        cout << ", " << group_2_alglib_array[ group_2_alglib_array_i ];
      }
      cout << " )" << endl;
      //unequalvariancettest( across_alglib_array, across_alglib_array.length(), group_2_alglib_array, group_2_alglib_array.length(), bothtails, lefttail, righttail );
      //cout << "two-sample unequal variance t-test across vs group 2: p-value is " << bothtails << endl;
      mannwhitneyutest( across_alglib_array, across_alglib_array.length(), group_2_alglib_array, group_2_alglib_array.length(), bothtails, lefttail, righttail );
      cout << "u-test across vs group 2: p-value is " << bothtails << endl;
      if( bothtails < .05 ) {
        cout << "==== SIGNIFICANT b->2 " << bothtails << " *";
        if( bothtails < ( .05 / 2 ) ) {
          cout << "!";
        }
        cout << endl;
      }
    }
  } // End if( num_seqs_in_group_1 )
  
  exit( 0 );

// / Setup for calculating the other stats

  // When we want to square the deviations, we need to be able to store negative numbers, so.. 
  galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, doublerealspace>::AlignmentProfile tmp_ap_dbl( 1 ); // TODO: Just use an AlignmentProfilePosition object?

  // Per-residue average across sequences
  // These are doubles -- note that they store totals until the end, when they are divided to become averages...
  galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, doublerealspace>::AlignmentProfile ap_deviations_from_profile_average( profile.length() + 1 );

  // Per-residue var across sequences, by position
  // MatrixValueType again, since the squares might be tiny.  Note again that these are sums until the end..
  galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile ap_squared_deviations_from_profile_average( profile.length() + 1 );

  // Per-residue overall average
  // Doubles again
  galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, doublerealspace>::AlignmentProfile ap_deviations_from_profile_average_over_positions( 1 ); // TODO: Just use an AlignmentProfilePosition object?

  // Per-residue overall var
  // MatrixValueTypes ..
  galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile ap_squared_deviations_from_profile_average_over_positions( 1 ); // TODO: Just use an AlignmentProfilePosition object?

  // Total across residues, foreach position, averaged across sequences:
  vector<double> ap_deviations_from_profile_totals_average( profile.length() + 1 );

  // Total across residues, foreach position, cov among sequences:
  galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::DistanceMatrix<doublerealspace> ap_deviations_from_profile_totals_cov( ( profile.length() + 1 ), true, true ); // true ( symmetric ), true ( proximity matrix )

  // Total across residues, foreach position, correlation among sequences:
  galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::DistanceMatrix<doublerealspace> ap_deviations_from_profile_totals_corr( ( profile.length() + 1 ), true, true ); // true ( symmetric ), true ( proximity matrix )

  // Total across residues, foreach sequence, averaged across positions:
  vector<double> ap_deviations_from_profile_totals_average_by_seq( aps.size() );

  // Total across residues, foreach sequence, cov across positions:
  galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::DistanceMatrix<doublerealspace> ap_deviations_from_profile_totals_cov_by_seq( aps.size(), true, true ); // true ( symmetric ), true ( proximity matrix )

  // Total across residues, foreach sequence, correlation across positions:
  galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::DistanceMatrix<doublerealspace> ap_deviations_from_profile_totals_corr_by_seq( aps.size(), true, true ); // true ( symmetric ), true ( proximity matrix )

  ////////// Calc stats:
  double tmp_total, tmp_total_k;
  for( uint32_t ap_i = 0; ap_i < aps.size(); ap_i++ ) {
    // The APs have one more pos than does the profile (just like the dp matrix does) -- so the first entry doesn't correspond to any match emission distribution.
    //for( uint32_t ap_pos_i = 0; ap_pos_i < aps[ ap_i ].length(); ap_pos_i++ ) {
    for( uint32_t ap_pos_i = 1; ap_pos_i < aps[ ap_i ].length(); ap_pos_i++ ) {
      ap_deviations_from_profile_average[ ap_pos_i ][ galosh::Emission::Match ] += ap_deviations_from_profile[ ap_i ][ ap_pos_i ][ galosh::Emission::Match ];
      ap_deviations_from_profile_average_over_positions[ 0 ][ galosh::Emission::Match ] += ap_deviations_from_profile[ ap_i ][ ap_pos_i ][ galosh::Emission::Match ];

      // The squared deviations
      tmp_ap_dbl[ 0 ][ galosh::Emission::Match ] = ap_deviations_from_profile[ ap_i ][ ap_pos_i ][ galosh::Emission::Match ];
      tmp_ap_dbl[ 0 ][ galosh::Emission::Match ] *= ap_deviations_from_profile[ ap_i ][ ap_pos_i ][ galosh::Emission::Match ];
      ap_squared_deviations_from_profile_average[ ap_pos_i ][ galosh::Emission::Match ] += tmp_ap_dbl[ 0 ][ galosh::Emission::Match ];
      ap_squared_deviations_from_profile_average_over_positions[ 0 ][ galosh::Emission::Match ] += tmp_ap_dbl[ 0 ][ galosh::Emission::Match ];

      tmp_total =
        toDouble( ap_deviations_from_profile[ ap_i ][ ap_pos_i ][ galosh::Emission::Match ].total() );
      ap_deviations_from_profile_totals_average[ ap_pos_i ] +=
        tmp_total;
      ap_deviations_from_profile_totals_average_by_seq[ ap_i ] +=
        tmp_total;

      // Covariances among positions
      for( uint32_t ap_pos_k = 1; ap_pos_k <= ap_pos_i; ap_pos_k++ ) {
        tmp_total_k =
          toDouble( ap_deviations_from_profile[ ap_i ][ ap_pos_k ][ galosh::Emission::Match ].total() );
        ap_deviations_from_profile_totals_cov( ap_pos_i, ap_pos_k ) +=
          ( tmp_total * tmp_total_k );
      } // End foreach position ap_pos_k...
      // Covariances among sequences
      for( uint32_t ap_k = 0; ap_k <= ap_i; ap_k++ ) {
        tmp_total_k =
          toDouble( ap_deviations_from_profile[ ap_k ][ ap_pos_i ][ galosh::Emission::Match ].total() );
        ap_deviations_from_profile_totals_cov_by_seq( ap_i, ap_k  ) +=
          ( tmp_total * tmp_total_k );
      } // End foreach position ap_k...
    } // End foreach position ap_pos_i...
    ap_deviations_from_profile_totals_average_by_seq[ ap_i ] /=
      ( profile.length() );  // would be + 1 but the first pos has no match emissions
    cout << "average total deviation per position for seq " << ap_i << ": " <<
      ap_deviations_from_profile_totals_average_by_seq[ ap_i ] << endl;
  } // End foreach alignment profile ap_i...

  // Finish calculating the covariances (divide by the totals)
  vector<double> sds_by_pos( aps[ 0 ].length() );
  vector<double> sds_by_seq( aps.size() );
  for( uint32_t ap_pos_i = 1; ap_pos_i < aps[ 0 ].length(); ap_pos_i++ ) {
    ap_deviations_from_profile_totals_average[ ap_pos_i ] /=
      ( aps.size() );
    cout << "average total deviation per seq for pos " << ap_pos_i << ": " <<
      ap_deviations_from_profile_totals_average[ ap_pos_i ] << endl;

    for( uint32_t ap_pos_k = 1; ap_pos_k <= ap_pos_i; ap_pos_k++ ) {
      ap_deviations_from_profile_totals_cov( ap_pos_i, ap_pos_k ) /=
        aps.size(); // number of seqs (no correction because we know the means)
    } // End foreach position ap_pos_k...
    sds_by_pos[ ap_pos_i ] = toDouble( pow( ap_deviations_from_profile_totals_cov( ap_pos_i, ap_pos_i ), .5 ) ); // sqrt of var
    cout << "SD for pos " << ap_pos_i << ": " << sds_by_pos[ ap_pos_i ] << endl;
  } // End foreach ap_pos_i
  for( uint32_t ap_i = 1; ap_i < aps.size(); ap_i++ ) {
    for( uint32_t ap_k = 1; ap_k <= ap_i; ap_k++ ) {
      ap_deviations_from_profile_totals_cov_by_seq( ap_i, ap_k ) /=
        profile.length(); // number of positions (no correction because we know the means) -- also note that we don't add 1, since the first ap has no match emissions...
    } // End foreach position ap_k...
    sds_by_seq[ ap_i ] = toDouble( pow( ap_deviations_from_profile_totals_cov_by_seq( ap_i, ap_i ), .5 ) ); // sqrt of var
    cout << "SD for seq " << ap_i << ": " << sds_by_seq[ ap_i ] << endl;
  } // End foreach ap_i

  exit( 0 );

  // Calculate the correlations (divide by the sds)
  double average_corr_among_positions = 0;
  double average_corr_among_sequences = 0;
  uint32_t tmp_count = 0;
  cout << "Correlations among positions:" << endl;
  for( uint32_t ap_pos_i = 1; ap_pos_i < aps[ 0 ].length(); ap_pos_i++ ) {
    cout << "Position " << ap_pos_i << ": ";
    for( uint32_t ap_pos_k = 1; ap_pos_k <= ap_pos_i; ap_pos_k++ ) {
      ap_deviations_from_profile_totals_corr( ap_pos_i, ap_pos_k ) =
        ap_deviations_from_profile_totals_cov( ap_pos_i, ap_pos_k );
      ap_deviations_from_profile_totals_corr( ap_pos_i, ap_pos_k ) /=
        ( sds_by_pos[ ap_pos_i ] * sds_by_pos[ ap_pos_k ] );
      cout << ap_deviations_from_profile_totals_corr( ap_pos_i, ap_pos_k ) << "\t";
      if( ap_pos_k < ap_pos_i ) {
        average_corr_among_positions +=
          toDouble( ap_deviations_from_profile_totals_corr( ap_pos_i, ap_pos_k ) );
        tmp_count += 1;
      }
    } // End foreach position ap_pos_k...
    cout << "\n";
  } // End foreach ap_pos_i
  average_corr_among_positions /= tmp_count;
  cout << "Average correlation among positions: " << average_corr_among_positions << endl;
  // ...now for the by_seq corrs:
  tmp_count = 0;
  cout << "Correlations among sequences:" << endl;
  for( uint32_t ap_i = 1; ap_i < aps.size(); ap_i++ ) {
    cout << "Sequence " << ap_i << ": ";
    for( uint32_t ap_k = 1; ap_k <= ap_i; ap_k++ ) {
      ap_deviations_from_profile_totals_corr_by_seq( ap_i, ap_k ) =
        ap_deviations_from_profile_totals_cov_by_seq( ap_i, ap_k );
      ap_deviations_from_profile_totals_corr_by_seq( ap_i, ap_k ) /=
        ( sds_by_seq[ ap_i ] * sds_by_seq[ ap_k ] );
      cout << ap_deviations_from_profile_totals_corr_by_seq( ap_i, ap_k ) << "\t";
      if( ap_k < ap_i ) {
        average_corr_among_sequences +=
          toDouble( ap_deviations_from_profile_totals_corr_by_seq( ap_i, ap_k ) );
        tmp_count += 1;
      }
    } // End foreach sequence ap_k...
    cout << "\n";
  } // End foreach ap_i
  average_corr_among_sequences /= tmp_count;
  cout << "Average correlation among sequences: " << average_corr_among_sequences << endl;

  exit( 0 );

  for( uint32_t ap_pos_i = 1; ap_pos_i < aps[ 0 ].length(); ap_pos_i++ ) {
    // ap_deviations_from_profile_average[ ap_pos_i ][ galosh::Emission::Match ] /=
    //   aps.size();
    // cout << "Average ap_devs [" << ap_pos_i << "] " << ap_deviations_from_profile_average[ ap_pos_i ][ galosh::Emission::Match ] << endl;

    // We know the true mean deviation (0) so we need't subtract 1 to get the var:
    ap_squared_deviations_from_profile_average[ ap_pos_i ][ galosh::Emission::Match ] /=
      aps.size();
    cout << "Var of ap_devs [" << ap_pos_i << "] " << ap_squared_deviations_from_profile_average[ ap_pos_i ][ galosh::Emission::Match ] << endl;
    // Also print expected var
    tmp_ap_dbl[ 0 ][ galosh::Emission::Match ] = 1.0;
    tmp_ap_dbl[ 0 ][ galosh::Emission::Match ] -= profile[ ap_pos_i - 1 ][ galosh::Emission::Match ];
    tmp_ap_dbl[ 0 ][ galosh::Emission::Match ] *= profile[ ap_pos_i - 1 ][ galosh::Emission::Match ];
    //( profile[ ap_pos_i - 1 ] * ( 1.0 - profile[ ap_pos_i - 1 ] ) );
    cout << "Expected vars at pos [" << ap_pos_i << "] " << tmp_ap_dbl[ 0 ][ galosh::Emission::Match ] << endl;
  } // End foreach position ap_pos_i...
  exit( 0 );

  ap_deviations_from_profile_average_over_positions[ 0 ][ galosh::Emission::Match ] /=
    ( aps.size() * ( profile.length() ) ); // Note that we use profile.length(), not plus one, because the 0th entry has no match emissions anyway.
  cout << "Average ap_devs overall, by residue: " << ap_deviations_from_profile_average_over_positions[ 0 ][ galosh::Emission::Match ] << endl;

  // We know the true mean deviation (0) so we need't subtract 1 to get the var:
  MatrixValueType average_squared_deviation_emissions = ap_squared_deviations_from_profile_average_over_positions[ 0 ][ galosh::Emission::Match ].total();
  average_squared_deviation_emissions /=
    ( seqan::ValueSize<ResidueType>::VALUE * aps.size() * ( profile.length() ) );

  // We know the true mean deviation (0) so we need't subtract 1 to get the var:
  ap_squared_deviations_from_profile_average_over_positions[ 0 ][ galosh::Emission::Match ] /=
    ( aps.size() * ( profile.length() ) ); // Note that we use profile.length(), not plus one, because the 0th entry has no match emissions anyway.
  cout << "Var of ap_devs overall, by residue: " << ap_squared_deviations_from_profile_average_over_positions[ 0 ][ galosh::Emission::Match ] << endl;

  MatrixValueType average_deviation_emissions = ap_deviations_from_profile_average_over_positions[ 0 ][ galosh::Emission::Match ].total();
  average_deviation_emissions /= seqan::ValueSize<ResidueType>::VALUE;

  cout << "Average ap_devs overall: " << average_deviation_emissions << endl;
  cout << "Var of ap_devs overall: " << average_squared_deviation_emissions << endl;

  exit( 0 );

} // main (..)
