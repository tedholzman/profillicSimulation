# DNA
DNA.alltogether.results <- as.matrix( read.table( "results/results_revised/DNA.processed.alltogether.out", header = T ) );
DNA.alltogether.truesseparated.results <- as.matrix( read.table( "results/results_revised/DNA.processed.alltogether.truesseparated.out", header = T ) );
DNA.evenstart.results <- as.matrix( read.table( "results/results_revised/DNA.processed.evenstart.out", header = T ) );
DNA.evenstart.truesseparated.results <- as.matrix( read.table( "results/results_revised/DNA.processed.evenstart.truesseparated.out", header = T ) );
DNA.priorstart.results <- as.matrix( read.table( "results/results_revised/DNA.processed.priorstart.out", header = T ) );
DNA.priorstart.truesseparated.results <- as.matrix( read.table( "results/results_revised/DNA.processed.priorstart.truesseparated.out", header = T ) );
DNA.uniformstart.results <- as.matrix( read.table( "results/results_revised/DNA.processed.uniformstart.out", header = T ) );
DNA.uniformstart.truesseparated.results <- as.matrix( read.table( "results/results_revised/DNA.processed.uniformstart.truesseparated.out", header = T ) );

# AA
AA.alltogether.results <- as.matrix( read.table( "results/results_revised/AA.processed.alltogether.out", header = T ) );
AA.alltogether.truesseparated.results <- as.matrix( read.table( "results/results_revised/AA.processed.alltogether.truesseparated.out", header = T ) );
AA.evenstart.results <- as.matrix( read.table( "results/results_revised/AA.processed.evenstart.out", header = T ) );
AA.evenstart.truesseparated.results <- as.matrix( read.table( "results/results_revised/AA.processed.evenstart.truesseparated.out", header = T ) );
AA.priorstart.results <- as.matrix( read.table( "results/results_revised/AA.processed.priorstart.out", header = T ) );
AA.priorstart.truesseparated.results <- as.matrix( read.table( "results/results_revised/AA.processed.priorstart.truesseparated.out", header = T ) );
AA.uniformstart.results <- as.matrix( read.table( "results/results_revised/AA.processed.uniformstart.out", header = T ) );
AA.uniformstart.truesseparated.results <- as.matrix( read.table( "results/results_revised/AA.processed.uniformstart.truesseparated.out", header = T ) );

library( "plotrix" ) # for "plotCI"
plotIt <- function( .the.results, .profile.length, .conservation.rate, .training.seqs, plot.x.CI = FALSE, plot.y.CI = FALSE, include.legend = TRUE ) {
  .title <- paste( c( paste( "Conservation rate", .conservation.rate, "with", .training.seqs, "training sequences" ), paste( "For profile length", .profile.length ), collapse = "\n" ) );
  .xlim <- c( min( .the.results[ ( .the.results[ , "profile_length" ] == .profile.length ) & ( .the.results[ , "num_training_sequences_per_profile" ] == .training.seqs ) & ( .the.results[ , "conservation_rate" ] == .conservation.rate ), c( "training_conditional_forward", "training_unconditional_forward", "training_conditionalBaldiSiegel_forward", "training_unconditionalBaldiSiegel_forward" ) ] ), max( .the.results[ ( .the.results[ , "profile_length" ] == .profile.length ) & ( .the.results[ , "num_training_sequences_per_profile" ] == .training.seqs ) & ( .the.results[ , "conservation_rate" ] == .conservation.rate ), c( "training_conditional_forward", "training_unconditional_forward", "training_conditionalBaldiSiegel_forward", "training_unconditionalBaldiSiegel_forward" ) ] ) );
  .ylim <- c( min( .the.results[ ( .the.results[ , "profile_length" ] == .profile.length ) & ( .the.results[ , "num_training_sequences_per_profile" ] == .training.seqs ) & ( .the.results[ , "conservation_rate" ] == .conservation.rate ), c( "test_conditional_forward", "test_unconditional_forward", "test_conditionalBaldiSiegel_forward", "test_unconditionalBaldiSiegel_forward" ) ] ), max( .the.results[ ( .the.results[ , "profile_length" ] == .profile.length ) & ( .the.results[ , "num_training_sequences_per_profile" ] == .training.seqs ) & ( .the.results[ , "conservation_rate" ] == .conservation.rate ), c( "test_conditional_forward", "test_unconditional_forward", "test_conditionalBaldiSiegel_forward", "test_unconditionalBaldiSiegel_forward" ) ] ) );
  # Make sure to include (0,0)
  if( .xlim[ 2 ] < 0 ) {
    .xlim[ 2 ] <- 0;
  }
  if( .ylim[ 2 ] < 0 ) {
    .ylim[ 2 ] <- 0;
  }
  
  # Initialize the plot
  plot( NA, xlim = .xlim, ylim = .ylim, main = .title, xlab = "Difference in training set log-liklihood", ylab = "Difference in test set log-liklihood" );

  # Plot those lines at 0
  lines( .xlim, c( 0, 0 ), col = "orange" );
  lines( c( 0, 0 ), .ylim, col = "orange" );
  
  the.methods <- c( "conditional", "unconditional", "conditionalBaldiSiegel", "unconditionalBaldiSiegel" );
  the.methods.colors <- rainbow( length( the.methods ) );
  the.methods.characters <- 20 + 1:length( the.methods );
  the.methods.legends <- c( "CBW", "UBW", "CQA", "UQA" );
  for( .method.i in 1:length( the.methods ) ) {
    .the.method <- the.methods[ .method.i ];
    .training.method.forward <- paste( "training", .the.method, "forward", sep = "_" );
    .test.method.forward <- paste( "test", .the.method, "forward", sep = "_" );
    points( .the.results[ ( .the.results[ , "profile_length" ] == .profile.length ) & ( .the.results[ , "num_training_sequences_per_profile" ] == .training.seqs ) & ( .the.results[ , "conservation_rate" ] == .conservation.rate ), c( .training.method.forward, .test.method.forward ) ], pch = the.methods.characters[ .method.i ], col = the.methods.colors[ .method.i ] );
    if( plot.x.CI ) {
      # x bars
      plotCI( .the.results[ ( .the.results[ , "profile_length" ] == .profile.length ) & ( .the.results[ , "num_training_sequences_per_profile" ] == .training.seqs ) & ( .the.results[ , "conservation_rate" ] == .conservation.rate ), .training.method.forward ], .the.results[ ( .the.results[ , "profile_length" ] == .profile.length ) & ( .the.results[ , "num_training_sequences_per_profile" ] == .training.seqs ) & ( .the.results[ , "conservation_rate" ] == .conservation.rate ), .training.method.forward ], err="x", uiw = ( 1.96 * ( .the.results[ ( .the.results[ , "profile_length" ] == .profile.length ) & ( .the.results[ , "num_training_sequences_per_profile" ] == .training.seqs ) & ( .the.results[ , "conservation_rate" ] == .conservation.rate ), paste( "sd", .training.method.forward, sep = "." ) ] / sqrt( se.n ) ) ), pch = the.methods.characters[ .method.i ], col = the.methods.colors[ .method.i ], add = TRUE )
    } # End if plot.x.CI
    if( plot.y.CI ) {
      # y bars
      plotCI( .the.results[ ( .the.results[ , "profile_length" ] == .profile.length ) & ( .the.results[ , "num_training_sequences_per_profile" ] == .training.seqs ) & ( .the.results[ , "conservation_rate" ] == .conservation.rate ), .training.method.forward ], .the.results[ ( .the.results[ , "profile_length" ] == .profile.length ) & ( .the.results[ , "num_training_sequences_per_profile" ] == .training.seqs ) & ( .the.results[ , "conservation_rate" ] == .conservation.rate ), .test.method.forward ], err="y", uiw = ( 1.96 * ( .the.results[ ( .the.results[ , "profile_length" ] == .profile.length ) & ( .the.results[ , "num_training_sequences_per_profile" ] == .training.seqs ) & ( .the.results[ , "conservation_rate" ] == .conservation.rate ), paste( "sd", .test.method.forward, sep = "." ) ] / sqrt( se.n ) ) ), pch = the.methods.characters[ .method.i ], col = the.methods.colors[ .method.i ], add = TRUE )
    } # End if plot.y.CI
  } # End foreach .method.i
  if( include.legend ) {
    legend( .xlim[ 1 ] + .05*( .xlim[ 2 ] - .xlim[ 1 ] ), .ylim[ 2 ] - .05*( .ylim[ 2 ] - .ylim[ 1 ] ), the.methods.legends, pch = the.methods.characters, col = the.methods.colors, text.col = the.methods.colors );
  } # End if include.legend
} # plotIt (..)
