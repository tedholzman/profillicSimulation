#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) ProfillicSimulationTabProcessor.pl
##  Author:
##      Paul T Edlefsen   Fred Hutchinson Cancer Research Center
##  Description:
##      Script for post-processing ProfillicSimulation results (.tab) files.
##      Output is averages; probabilities are log10-transformed first.
##
#******************************************************************************

use Getopt::Std;

use strict;
use vars qw( $VERSION $DEBUG $VERBOSE
             $opt_D $opt_V $opt_o $opt_O
             $opt_t $opt_a $opt_d $opt_l $opt_i $opt_f $opt_v $opt_T
             $opt_s $opt_S $opt_z $opt_x $opt_X );
$VERSION = '1.2';

# This means -D, -o, -O, and -V are ok, but nothin' else.
# opt_D means print debugging output.
# opt_o is an optional filename to put the output.  Otherwise, STDOUT.
# opt_O is just like opt_o except it'll overwrite the file if it exists.
# opt_V means be verbose.
# opt_t means keep true profiles separate in the output.
# opt_a means show the profile-profile alignment scores in the output.
# opt_d means show the SKL divergences in the output.
# opt_l means show profile lengths in the output.
# opt_i means show training iterations in the output.
# opt_f means show the forward scores in the output.
# opt_v means show the viterbi scores in the output.
# opt_T means show the testing scores only (not the training scores).
# opt_s means also calculate and show standard deviations (in addition to means)
# opt_S means also calculate and show standard errors (in addition to means)
# opt_z means show differences to the true, rather than the actual, log-prob.
# opt_x is the lowest starting profile index to use (default 0)
# opt_X is the highest starting profile index to use (default, the max)
if( not getopts('DVo:O:tadlifvTsSzx:X:') ) {
  usage();
}

$DEBUG = $opt_D;
$VERBOSE = $opt_V;

my $keep_true_profiles_separated = $opt_t;
my $show_PPAlign_scores = $opt_a;
my $show_SKL_divergences = $opt_d;
my $show_profile_lengths = $opt_l;
my $show_training_iters = $opt_i;
my $show_forward_scores = $opt_f;
my $show_viterbi_scores = $opt_v;
my $show_testing_scores = 1;
my $show_training_scores = !$opt_T;

my $show_standard_deviations = $opt_s;
my $show_standard_errors = $opt_S;

my $diffs_to_true = $opt_z;

my $use_starting_profile_indices_low = $opt_x;
my $use_starting_profile_indices_high = $opt_X;

## By default we show both forward and viterbi scores.
unless( $show_forward_scores || $show_viterbi_scores ) {
  $show_forward_scores = 1;
  $show_viterbi_scores = 1;
}

if( $VERBOSE ) {
  $| = 1; # Autoflush.
}

# ProfillicSimulation output can be put through toLogDouble(), or not.  We auto-detect this.
my $output_is_LogDouble = -1; # -1 means we don't know yet

my $output_file = $opt_o || $opt_O;
if( $output_file ) {

  if( !$opt_O && -e $output_file ) {
    print "The file \"$output_file\" exists.  Use -O to overwrite.\n";
    exit;
  }

  if( $VERBOSE ) { print "Opening file \"$output_file\" for writing.."; }

  open OUTPUT_FH, ">$output_file";

  if( $VERBOSE ) { print ".done.\n"; }

} else {

  open OUTPUT_FH, ">&STDOUT";

}

my $tab_file = shift || usage();

if( $VERBOSE ) { print "Opening file \"$tab_file\".."; }

open( TAB_FH, $tab_file ) or
  die "Unable to open ProfillicSimulation .tab file \"$tab_file\": $!";

if( $VERBOSE ) { print ".done.\n"; }

if( $VERBOSE ) { print "Parsing ProfillicSimulation .tab file.."; }

if( $DEBUG ) { print "[DEBUG]\n" };

sub scientificNotationToLogBase10 {
  my $scientific_notation_val = shift;

  my ( $coefficientBeforeDecimal, $coefficientAfterDecimal, $exponent ) = ( $scientific_notation_val =~ /^(\d)\.(\d+)e(-.+)/ );
  unless( defined $exponent ) {
    # Not in scientific notation, clearly.  Just return it as-is.
    return $scientific_notation_val;
  }
  return(  $exponent + ".$coefficientBeforeDecimal$coefficientAfterDecimal" );
} # End sub scientificNotationToLogBase10( $scientific_notation_val )

sub naturalLogToLogBase10 {
  my $natural_log_val = shift;

  if( $natural_log_val !~ /^-/ ) {
    # Not a natrual log value, clearly.  Just return it.
    return $natural_log_val;
  }
  # multiply by log10( e )
  return ( $natural_log_val * log10( exp( 1 ) ) );
} # End sub naturalLogToLogBase10( $natural_log_val )

my @column_headers;
my @data_column_is_true; # true if column corresponds to a "true" profile score
my @data_column_is_PPAlign; # true if column corresponds to a profile-profile alignment score
my @data_column_is_SKL; # true if column corresponds to an SKL distance measure
my @data_column_is_length; # true if column corresponds to a profile length
my @data_column_is_iters; # true if column corresponds to the number of iterations it took
my @data_column_is_viterbi; # true if column corresponds to a viterbi score
my @data_column_is_forward; # true if column corresponds to a forward score
my @data_column_is_training; # true if column corresponds to a training score
my @data_column_is_testing; # true if column corresponds to a test score
my ( $true_profile_id_key_column, $true_profile_id_line_column );
my @data;
my ( @test_descriptors, $test_descriptor, $last_test_descriptor );
my ( $true_profile_id, $log10_true_value, $log10_value );
my @line_values;
my @filtered_values;
my $only_one_true_profile = 1; # True until we see the column 'true_profile_id'
my @true_column_i_for_test_column_i;
my %true_data_columns;
my ( $true_line_value_i, $true_column_i, $newkey );
my $starting_profile_index = 0;
while( <TAB_FH> ) {

  # Chop and Chomp won't remove ^Ms or leading ws.
  ( $_ ) = ( $_ =~ /^\s*(.+?)\s*$/ );

  if( 0 && $DEBUG ) {
    print "TAB LINE: $_\n";
  }

  ## Shouldn't have blank lines.
  #next unless( $_ );

  unless( scalar( @column_headers ) ) {
    # First line:
    @column_headers = split( "\t", $_ );
    # Find where the 'true_profile_id' column is.  As a workaround to
    # a bug in some early .tab files (in which the column headers
    # before 'true_profile_id' were omitted), we compute the number of
    # data columns by treating the first column after the
    # true_profile_id_key_column as the first data column, then in the
    # data rows we get the data by subtracting that number of columns
    # from the end...
    # 
    # Note that when there is only one true profile, the column is
    # omitted from the output.
    $true_profile_id_key_column = 0;
    for( my $column_i; $column_i <= $#column_headers; $column_i++ ) {
      if( $column_headers[ $column_i ] eq 'expected_insertion_length_as_profile_length_fraction' ) {
        ## This saves us if there is only one true profile (and thus no colunmn called 'true_profile_id').
        $true_profile_id_key_column = $column_i;
        next;
      }
      if( $column_headers[ $column_i ] eq 'true_profile_id' ) {
        $true_profile_id_key_column = $column_i;
        $only_one_true_profile = 0; # There must be more than one.
        next;
      }
      if( $column_headers[ $column_i ] =~ m/_true_/ ) {
        # The @data_column_is_true array is indexed from the first
        # column after the $true_profile_id_key_column.
        $data_column_is_true[ $column_i - ( $true_profile_id_key_column + 1 ) ] = 1;
        ## Set up the map between test columns and their corresponding true columns.
        $newkey = $column_headers[ $column_i ];
        $newkey =~ s/_true_/_/;
##        print( "SETTING UP MAPPING FOR column $column_i ($column_headers[ $column_i ]): key is $newkey\n" );
        $true_data_columns{ $newkey } = $column_i;
      }

      ## Note that we assume that the 'true_profile_id' column comes
      ## before any forward or viterbi score column.
      if( $column_headers[ $column_i ] =~ m/PPAlign/ ) { # Must be processed before /SKL/
        # The @data_column_is_PPAlign array is indexed from the first
        # column after the $true_profile_id_key_column.
        $data_column_is_PPAlign[ $column_i - ( $true_profile_id_key_column + 1 ) ] = 1;
      } elsif( $column_headers[ $column_i ] =~ m/SKL/ ) {
        # The @data_column_is_SKL array is indexed from the first
        # column after the $true_profile_id_key_column.
        $data_column_is_SKL[ $column_i - ( $true_profile_id_key_column + 1 ) ] = 1;
      } elsif( $column_headers[ $column_i ] =~ m/length/ ) {
        # The @data_column_is_length array is indexed from the first
        # column after the $true_profile_id_key_column.
        $data_column_is_length[ $column_i - ( $true_profile_id_key_column + 1 ) ] = 1;
      } elsif( $column_headers[ $column_i ] =~ m/iters/ ) {
        # The @data_column_is_iters array is indexed from the first
        # column after the $true_profile_id_key_column.
        $data_column_is_iters[ $column_i - ( $true_profile_id_key_column + 1 ) ] = 1;
      } elsif( $column_headers[ $column_i ] =~ m/viterbi/ ) {
        # The @data_column_is_viterbi array is indexed from the first
        # column after the $true_profile_id_key_column.
        $data_column_is_viterbi[ $column_i - ( $true_profile_id_key_column + 1 ) ] = 1;
      } elsif( $column_headers[ $column_i ] =~ m/forward/ ) {
        # The @data_column_is_forward array is indexed from the first
        # column after the $true_profile_id_key_column.
        $data_column_is_forward[ $column_i - ( $true_profile_id_key_column + 1 ) ] = 1;
      }
      if( $column_headers[ $column_i ] =~ m/training/ ) {
        # The @data_column_is_training array is indexed from the first
        # column after the $true_profile_id_key_column.
        $data_column_is_training[ $column_i - ( $true_profile_id_key_column + 1 ) ] = 1;
      } elsif( $column_headers[ $column_i ] =~ m/test_/ ) {
        # The @data_column_is_testing array is indexed from the first
        # column after the $true_profile_id_key_column.
        $data_column_is_testing[ $column_i - ( $true_profile_id_key_column + 1 ) ] = 1;
      }
      if( !$data_column_is_true[ $column_i ] ) {
        ## Set up the map between other columns and their corresponding true columns.
        $newkey = $column_headers[ $column_i ];
        $newkey =~ s/_(?:un)?(?:conditional|starting)[^_]*_/_/;
##TODO: REMOVE
##        print( "KEY IS $newkey .. TRUE COL FOR $column_headers[ $column_i ] IS $true_data_columns{ $newkey }\n" );
        $true_column_i_for_test_column_i[ $column_i ] = $true_data_columns{ $newkey };
      }
    } # End foreach header column_i...

    next;
  }

  @line_values = split( "\t", $_ );

  unless( defined $true_profile_id_line_column ) {
    ## See note above, where we set the $true_profile_id_key_column.
    if( 0 && $DEBUG ) {
      print( "There are ", scalar( @line_values ), " line values: ", join( ",", @line_values ) );
      print( "There are ", scalar( @column_headers ), " column headers: ", join( ",", @column_headers ) );
    }
    $true_profile_id_line_column =
      ( $#line_values - ( $#column_headers - $true_profile_id_key_column ) );
    if( 0 && $DEBUG ) {
      print( "So \$true_profile_id_line_column is $true_profile_id_line_column\n" );
    }
  }

  if( $only_one_true_profile ) {
    $true_profile_id = 0;
    if( $DEBUG ) {
      print "TRUE PROFILE ID: 0\n(there's only one true profile per settings combination in the tab file)\n";
    }
  } else {
    $true_profile_id = $line_values[ $true_profile_id_line_column ];
    if( 0 && $DEBUG ) {
      print "TRUE PROFILE ID: $true_profile_id\n(parsed from element $true_profile_id_line_column of: (", join( ",", @line_values ), ")\n";
    }
  }
  # We store separately the test descriptors corresponding to each set
  # of test parameters, which is *all* of the text in the columns up
  # to (but not including) the $true_profile_id_line_column, with the
  # tabs put back.
  $test_descriptor =
    join( "\t", @line_values[ 0 .. ( $true_profile_id_line_column - 1 ) ] );
  if( 0 && $DEBUG ) {
    print "TEST DESCRIPTOR: $test_descriptor\n";
  }
  if( $test_descriptor ne $last_test_descriptor ) {
    push( @test_descriptors, $test_descriptor );
  }

  # @data is an array in which each entry corresponds to a set of test
  # parameters.  It contains a reference to an array in which each
  # entry corresponds to a "true" profile, and each of these arrays
  # contains a reference to another array in which each entry
  # corresponds to a different starting profile.  *those* arrays
  # contain the test results (again as an array reference).
  if( $test_descriptor ne $last_test_descriptor ) {
    push( @data, [] );
    $starting_profile_index = 0;
  } else {
    $starting_profile_index += 1;
  }
  $last_test_descriptor = $test_descriptor;

  if( defined( $use_starting_profile_indices_low ) && ( $starting_profile_index < $use_starting_profile_indices_low ) ) {
    ## TODO: REMOVE
#print( "skipping $starting_profile_index < $use_starting_profile_indices_low\n" );
    next;
  }
  if( defined( $use_starting_profile_indices_high ) && ( $starting_profile_index > $use_starting_profile_indices_high ) ) {
    ## TODO: REMOVE
#print( "skipping $starting_profile_index > $use_starting_profile_indices_high\n" );
    next;
  }
## TODO: REMOVE
#  print( "starting_profile_index: $starting_profile_index\n" );
  if( !$diffs_to_true && ( $show_PPAlign_scores && $show_SKL_divergences && $show_profile_lengths && $show_training_iters && $show_viterbi_scores && $show_forward_scores && $show_training_scores && $show_testing_scores ) ) {
    ## Not filtered..
    @filtered_values =
      @line_values[ ( $true_profile_id_line_column + 1 ) .. $#line_values ];
    if( $output_is_LogDouble < 0 ) {
      if( $filtered_values[ 0 ] =~ /^-/ ) {
        $output_is_LogDouble = 1; # log values start with a negative sign.
      } else {
        $output_is_LogDouble = 0;
      }
    } # End if we haven't yet determined whether the output is in log format or in raw (exponential) format aka scientific notation.
    # Convert to log base 10.
    if( !$output_is_LogDouble ) {
      # Then the number is in scientific notation.
      @filtered_values = map { scientificNotationToLogBase10( $_ ) } @filtered_values;
    } else {
      @filtered_values = map { naturalLogToLogBase10( $_ ) } @filtered_values;
    }
  } else {

    @filtered_values = ();
    for( my $line_value_i = ( $true_profile_id_line_column + 1 );
         $line_value_i <= $#line_values;
         $line_value_i++
    ) {
      if(
         ( $show_PPAlign_scores && $data_column_is_PPAlign[ $line_value_i - ( $true_profile_id_line_column + 1 ) ] ) ||
         ( $show_SKL_divergences && $data_column_is_SKL[ $line_value_i - ( $true_profile_id_line_column + 1 ) ] ) ||
         ( $show_profile_lengths && $data_column_is_length[ $line_value_i - ( $true_profile_id_line_column + 1 ) ] ) ||
         ( $show_training_iters && $data_column_is_iters[ $line_value_i - ( $true_profile_id_line_column + 1 ) ] )
      ) {
        if( $output_is_LogDouble < 0 ) {
          if( $line_values[ $line_value_i ] =~ /^-/ ) {
            $output_is_LogDouble = 1; # log values start with a negative sign.
          } else {
            $output_is_LogDouble = 0;
          }
        } # End if we haven't yet determined whether the output is in log format or in raw (exponential) format aka scientific notation.
        # Convert to log base 10.
        if( !$output_is_LogDouble ) {
          # Then the number is in scientific notation.
          $log10_value = scientificNotationToLogBase10( $line_values[ $line_value_i ] );
        } else {
          $log10_value = naturalLogToLogBase10( $line_values[ $line_value_i ] );
        }
        push( @filtered_values, $log10_value );
      } elsif(
          (
           ( $show_viterbi_scores && $data_column_is_viterbi[ $line_value_i - ( $true_profile_id_line_column + 1 ) ] ) ||
            ( $show_forward_scores && $data_column_is_forward[ $line_value_i - ( $true_profile_id_line_column + 1 ) ] )
          ) &&
          ( ( $show_training_scores && $data_column_is_training[ $line_value_i - ( $true_profile_id_line_column + 1 ) ] ) ||
            ( $show_testing_scores && $data_column_is_testing[ $line_value_i - ( $true_profile_id_line_column + 1 ) ] )
          )
               ) {
        if( $output_is_LogDouble < 0 ) {
          if( $line_values[ $line_value_i ] =~ /^-/ ) {
            $output_is_LogDouble = 1; # log values start with a negative sign.
          } else {
            $output_is_LogDouble = 0;
          }
        } # End if we haven't yet determined whether the output is in log format or in raw (exponential) format aka scientific notation.
        # Convert to log base 10.
        if( !$output_is_LogDouble ) {
          # Then the number is in scientific notation.
          $log10_value = scientificNotationToLogBase10( $line_values[ $line_value_i ] );
        } else {
          $log10_value = naturalLogToLogBase10( $line_values[ $line_value_i ] );
        }
        if( $diffs_to_true ) {
          if( $data_column_is_true[ $line_value_i - ( $true_profile_id_line_column + 1 ) ] ) {
            $log10_true_value = $log10_value;
            ## TODO: REMOVE
            #print( "true: $log10_true_value\n" );
          } else {
            # The relevant true column is this colname with "test" replaced by "true".
            $true_column_i =
              $true_column_i_for_test_column_i[ $line_value_i - ( $true_profile_id_line_column + 1 ) ];
            $true_line_value_i = $true_column_i + ( $true_profile_id_line_column + 1 );
            if( 0 && $DEBUG ) {
              print( "for $column_headers[ $line_value_i - ( $true_profile_id_line_column + 1 ) ] we use true column $column_headers[ $true_line_value_i - ( $true_profile_id_line_column + 1 ) ]\n" );
            }
            # Convert to log base 10.
            if( !$output_is_LogDouble ) {
              # Then the number is in scientific notation.
              $log10_true_value = scientificNotationToLogBase10( $line_values[ $true_line_value_i ] );
            } else {
              $log10_true_value = naturalLogToLogBase10( $line_values[ $true_line_value_i ] );
            }
            ## TODO: REMOVE
            #print( "from: $log10_value\t" );
            $log10_value -= $log10_true_value;
            ## TODO: REMOVE
            #print( "to: $log10_value\n" );
          }
        }
        push( @filtered_values, $log10_value );
      }
    } # End foreach $line_value_i
  } # End unless( we should filter the data ) .. else ..

  if( 0 && $DEBUG ) {
    print( "filtered, this is: ", join( ",", @filtered_values ), "\n" );
  }
  push( @{ $data[ $#data ]->[ $true_profile_id ] }, [ @filtered_values ] );

} # End while( <TAB_FH> )

if( $VERBOSE ) { print ".done.\n"; }

if( $VERBOSE ) { print "Closing file \"$tab_file\".."; }
close TAB_FH;
if( $VERBOSE ) { print ".done.\n"; }

### TODO: REMOVE
#use Data::Dumper;
#print Dumper( \@data );

# First calculate the filtered column headers.
my @filtered_column_headers;

  # put the @column_headers (up to and including the true_profile_id
  # column header) in the right place, accounting for the different
  # lengths of the header line and the data lines.
  for( my $blank_column_i = 0;
       $blank_column_i < ( $true_profile_id_line_column - $true_profile_id_key_column );
       $blank_column_i++
   ) {
    push( @filtered_column_headers, "" );
  }
  if( $keep_true_profiles_separated ) {
    push(
      @filtered_column_headers,
      @column_headers[ 0 .. $true_profile_id_key_column ]
    );
  } else {
    if( $true_profile_id_key_column > 0 ) {
      push(
        @filtered_column_headers,
        @column_headers[ 0 .. ( $true_profile_id_key_column - 1 ) ]
      );
    }
  }

### TODO:REMOVE
#print "\$true_profile_id_line_column is $true_profile_id_line_column\n";
#print "\$true_profile_id_key_column is $true_profile_id_key_column\n";
#print "\@column_headers is ( ", join( ", ", @column_headers ), " )\n";
#print "\@column_headers[ 0 .. ( $true_profile_id_key_column - 1 ) ] is ( ", join( ", ", @column_headers[ 0 .. ( $true_profile_id_key_column - 1 ) ] ), " )\n";
#print "\@filtered_column_headers is now ( ", join( ", ", @filtered_column_headers ), " )\n";
for( my $column_i = ( $true_profile_id_key_column + 1 );
     $column_i <= $#column_headers;
     $column_i++
) {
  if(
     ( $show_PPAlign_scores && $data_column_is_PPAlign[ $column_i - ( $true_profile_id_line_column + 1 ) ] ) ||
     ( $show_SKL_divergences && $data_column_is_SKL[ $column_i - ( $true_profile_id_line_column + 1 ) ] ) ||
     ( $show_profile_lengths && $data_column_is_length[ $column_i - ( $true_profile_id_line_column + 1 ) ] ) ||
     ( $show_training_iters && $data_column_is_iters[ $column_i - ( $true_profile_id_line_column + 1 ) ] ) ||
     (
      (
       ( $show_viterbi_scores && $data_column_is_viterbi[ $column_i - ( $true_profile_id_line_column + 1 ) ] ) ||
        ( $show_forward_scores && $data_column_is_forward[ $column_i - ( $true_profile_id_line_column + 1 ) ] )
      ) &&
      ( ( $show_training_scores && $data_column_is_training[ $column_i - ( $true_profile_id_line_column + 1 ) ] ) ||
        ( $show_testing_scores && $data_column_is_testing[ $column_i - ( $true_profile_id_line_column + 1 ) ] )
      )
     )
  ) {
    push( @filtered_column_headers, $column_headers[ $column_i ] );
    if( $show_standard_deviations ) {
      push( @filtered_column_headers, 'sd.' . $column_headers[ $column_i ] );
    }
    if( $show_standard_errors ) {
      push( @filtered_column_headers, 'se.' . $column_headers[ $column_i ] );
    }
  }
} # End foreach $column_i..

if( $VERBOSE ) { print "Calculating averages.."; }

my @test_descriptor_averages;
my @true_profile_averages;
for( my $test_descriptor_i = 0; $test_descriptor_i <= $#test_descriptors; $test_descriptor_i++ ) {
  $test_descriptor_averages[ $test_descriptor_i ] = [];
  $true_profile_averages[ $test_descriptor_i ] = [];
  for( my $true_profile_i = 0; $true_profile_i <= $#{ $data[ $test_descriptor_i ] }; $true_profile_i++ ) {
    $true_profile_averages[ $test_descriptor_i ][ $true_profile_i ] = [];
    for( my $starting_profile_i = 0; $starting_profile_i <= $#{ $data[ $test_descriptor_i ]->[ $true_profile_i ] }; $starting_profile_i++ ) {
      for( my $test_result_i = 0; $test_result_i <= $#{ $data[ $test_descriptor_i ]->[ $true_profile_i ]->[ $starting_profile_i ] }; $test_result_i++ ) {
        $true_profile_averages[ $test_descriptor_i ][ $true_profile_i ][ $test_result_i ] +=
          $data[ $test_descriptor_i ]->[ $true_profile_i ]->[ $starting_profile_i ]->[ $test_result_i ];
  
        if( $starting_profile_i == $#{ $data[ $test_descriptor_i ]->[ $true_profile_i ] } ) {
          # If it is the last of the values to be averaged, divide by the total
          $true_profile_averages[ $test_descriptor_i ][ $true_profile_i ][ $test_result_i ] /=
            scalar( @{ $data[ $test_descriptor_i ]->[ $true_profile_i ] } );
        }
      } # End foreach $test_result_i
    } # End foreach $starting_profile_i
    for( my $test_result_i = 0; $test_result_i < scalar( @{ $true_profile_averages[ $test_descriptor_i ][ $true_profile_i ] } ); $test_result_i++ ) {
      $test_descriptor_averages[ $test_descriptor_i ]->[ $test_result_i ] +=
        $true_profile_averages[ $test_descriptor_i ][ $true_profile_i ][ $test_result_i ];

      if( $true_profile_i == $#{ $data[ $test_descriptor_i ] } ) {
        # If it is the last of the values to be averaged, divide by the total.  Note it's an average of averages.
        $test_descriptor_averages[ $test_descriptor_i ]->[ $test_result_i ] /=
          scalar( @{ $data[ $test_descriptor_i ] } );
      }
    } # End foreach $test_result_i
  } # End foreach true_profile_i
} # End foreach test descriptor

if( $VERBOSE ) { print ".done.\n"; }

my @test_descriptor_stdevs; # At first, sum of squares...
my @true_profile_stdevs; # At first, sum of squares...
my @test_descriptor_stderrs;
my @true_profile_stderrs;
if( $show_standard_deviations || $show_standard_errors ) {

  if( $VERBOSE ) { print "Calculating standard deviations.."; }

  my $tmp;
  for( my $test_descriptor_i = 0; $test_descriptor_i <= $#test_descriptors; $test_descriptor_i++ ) {
    $test_descriptor_stdevs[ $test_descriptor_i ] = [];
    $true_profile_stdevs[ $test_descriptor_i ] = [];
    $test_descriptor_stderrs[ $test_descriptor_i ] = [];
    $true_profile_stderrs[ $test_descriptor_i ] = [];
    for( my $true_profile_i = 0; $true_profile_i <= $#{ $data[ $test_descriptor_i ] }; $true_profile_i++ ) {
      $true_profile_stdevs[ $test_descriptor_i ][ $true_profile_i ] = [];
      $true_profile_stderrs[ $test_descriptor_i ][ $true_profile_i ] = [];
      for( my $starting_profile_i = 0; $starting_profile_i <= $#{ $data[ $test_descriptor_i ]->[ $true_profile_i ] }; $starting_profile_i++ ) {
        for( my $test_result_i = 0; $test_result_i <= $#{ $data[ $test_descriptor_i ]->[ $true_profile_i ]->[ $starting_profile_i ] }; $test_result_i++ ) {
          $true_profile_stdevs[ $test_descriptor_i ][ $true_profile_i ][ $test_result_i ] +=
            ( ( $data[ $test_descriptor_i ]->[ $true_profile_i ]->[ $starting_profile_i ]->[ $test_result_i ] - $true_profile_averages[ $test_descriptor_i ][ $true_profile_i ][ $test_result_i ] ) ** 2 );
          $test_descriptor_stdevs[ $test_descriptor_i ]->[ $test_result_i ] +=
            ( ( $data[ $test_descriptor_i ]->[ $true_profile_i ]->[ $starting_profile_i ]->[ $test_result_i ] - $test_descriptor_averages[ $test_descriptor_i ]->[ $test_result_i ] ) ** 2 );
    
          if( $starting_profile_i == $#{ $data[ $test_descriptor_i ]->[ $true_profile_i ] } ) {
            if( $starting_profile_i == 0 ) {
              # Oh, it's the only starting profile for this true profile... so there's no stdev among the starting profiles..
              $true_profile_stdevs[ $test_descriptor_i ][ $true_profile_i ][ $test_result_i ] = undef;
              $true_profile_stderrs[ $test_descriptor_i ][ $true_profile_i ][ $test_result_i ] = undef;
            } else {
              # If it is the last of the values, divide by the total - 1, then sqrt.
              $true_profile_stdevs[ $test_descriptor_i ][ $true_profile_i ][ $test_result_i ] /=
                ( scalar( @{ $data[ $test_descriptor_i ]->[ $true_profile_i ] } ) - 1 );
              $true_profile_stdevs[ $test_descriptor_i ][ $true_profile_i ][ $test_result_i ] =
                $true_profile_stdevs[ $test_descriptor_i ][ $true_profile_i ][ $test_result_i ] ** .5;
              $true_profile_stderrs[ $test_descriptor_i ][ $true_profile_i ][ $test_result_i ] =
                $true_profile_stdevs[ $test_descriptor_i ][ $true_profile_i ][ $test_result_i ] / ( scalar( @{ $data[ $test_descriptor_i ]->[ $true_profile_i ] } ) ** .5 );
            }
            if( ( $true_profile_i == $#{ $data[ $test_descriptor_i ] } ) &&
                ( ( scalar( @{ $data[ $test_descriptor_i ]->[ $true_profile_i ] } ) *
                    scalar( @{ $data[ $test_descriptor_i ] } ) ) > 1 ) ) {
              $test_descriptor_stdevs[ $test_descriptor_i ]->[ $test_result_i ] /=
                (
                 ## Note that we're assuming that there's the same
                 ## number of starting profiles for each true profile...
                 ( scalar( @{ $data[ $test_descriptor_i ]->[ $true_profile_i ] } ) *
                   scalar( @{ $data[ $test_descriptor_i ] } ) )
                 - 1
                );
              $test_descriptor_stdevs[ $test_descriptor_i ]->[ $test_result_i ] =
                $test_descriptor_stdevs[ $test_descriptor_i ]->[ $test_result_i ] ** .5;
              $test_descriptor_stderrs[ $test_descriptor_i ]->[ $test_result_i ] =
                $test_descriptor_stdevs[ $test_descriptor_i ]->[ $test_result_i ] /
                (
                 ## Note that we're assuming that there's the same
                 ## number of starting profiles for each true profile...
                 ( scalar( @{ $data[ $test_descriptor_i ]->[ $true_profile_i ] } ) *
                   scalar( @{ $data[ $test_descriptor_i ] } ) )
                 ** .5
                );

            } # End if this is the last true profile too
          } # End if this is the last starting profile
        } # End foreach $test_result_i
      } # End foreach $starting_profile_i
    } # End foreach true_profile_i
  } # End foreach test descriptor
  
  if( $VERBOSE ) { print ".done.\n"; }

} # End if $show_standard_deviations || $show_standard_errors

if( $VERBOSE ) { print "Writing processed file.."; }
# Header
print OUTPUT_FH join( "\t", @filtered_column_headers ), "\n";

for( my $test_descriptor_i = 0; $test_descriptor_i <= $#test_descriptors; $test_descriptor_i++ ) {
  unless( $keep_true_profiles_separated ) {
    print OUTPUT_FH $test_descriptors[ $test_descriptor_i ];
  }
  for( my $true_profile_i = 0; $true_profile_i <= $#{ $data[ $test_descriptor_i ] }; $true_profile_i++ ) {
    if( $keep_true_profiles_separated ) {
      print OUTPUT_FH $test_descriptors[ $test_descriptor_i ], "\t", $true_profile_i;
    }
    for( my $starting_profile_i = 0; $starting_profile_i <= $#{ $data[ $test_descriptor_i ]->[ $true_profile_i ] }; $starting_profile_i++ ) {
      for( my $test_result_i = 0; $test_result_i <= $#{ $data[ $test_descriptor_i ]->[ $true_profile_i ]->[ $starting_profile_i ] }; $test_result_i++ ) {
        if( $starting_profile_i == $#{ $data[ $test_descriptor_i ]->[ $true_profile_i ] } ) {
          # .. and print.
          if( $keep_true_profiles_separated ) {
            print OUTPUT_FH "\t", $true_profile_averages[ $test_descriptor_i ][ $true_profile_i ][ $test_result_i ];
            if( $show_standard_deviations ) {
              print OUTPUT_FH "\t", $true_profile_stdevs[ $test_descriptor_i ][ $true_profile_i ][ $test_result_i ];
            }
            if( $show_standard_errors ) {
              print OUTPUT_FH "\t", $true_profile_stderrs[ $test_descriptor_i ][ $true_profile_i ][ $test_result_i ];
            }
          }
        }
      } # End foreach $test_result_i
    } # End foreach $starting_profile_i
    if( $keep_true_profiles_separated ) {
      print OUTPUT_FH "\n";
    }
    for( my $test_result_i = 0; $test_result_i < scalar( @{ $true_profile_averages[ $test_descriptor_i ][ $true_profile_i ] } ); $test_result_i++ ) {
      if( $true_profile_i == $#{ $data[ $test_descriptor_i ] } ) {
        # .. and print.
        unless( $keep_true_profiles_separated ) {
          print OUTPUT_FH "\t", $test_descriptor_averages[ $test_descriptor_i ]->[ $test_result_i ];
          if( $show_standard_deviations ) {
            print OUTPUT_FH "\t", $test_descriptor_stdevs[ $test_descriptor_i ]->[ $test_result_i ];
          }
          if( $show_standard_errors ) {
            print OUTPUT_FH "\t", $test_descriptor_stderrs[ $test_descriptor_i ]->[ $test_result_i ];
          }
        }
      }
    } # End foreach $test_result_i
  } # End foreach true_profile_i
  unless( $keep_true_profiles_separated ) {
    print OUTPUT_FH "\n";
  }
} # End foreach test descriptor

if( $VERBOSE ) { print ".done.\n"; }

if( $VERBOSE && $output_file ) { print "Closing file \"$output_file\".."; }
close OUTPUT_FH;
if( $VERBOSE && $output_file ) { print ".done.\n"; }

##---------------------------------------------------------------------------##

sub usage {
  print "\tProfillicSimulationTabProcessor.pl [-DV] [-(o|O) <Output File>] <.tab File>\n";
  exit;
}

1;

