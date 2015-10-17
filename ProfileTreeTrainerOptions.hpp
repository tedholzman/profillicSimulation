/**
 * \file ProfileTrainerOptions.hpp
 * \author Ted Holzman 
 * \par Library:
 * Galosh profillic
 * \brief Options for ProfileTrainer
 * \copyright &copy; 2008, 2011, 2012, 2013 by Paul T. Edlefsen, Fred Hutchinson Cancer
 *    Research Center.
 *  All rights reserved.
 *****************************************************************************/

///Useful for modifying GALOSH_DEF_OPT behavior.  See vector-valued options below
#ifndef TMP_EXTRA_STUFF
#define TMP_EXTRA_STUFF BOOST_PP_EMPTY()
#endif

/**
 * \note Note the lack of "protective" \#define symbols in this file.  We may
 * want to include it several times, redefining GALOSH_DEF_OPT as needed.
 */

      /**
       * After training a subfamily profile with all positions differing from
       * its parent profile, should we try to retrain after identifying the
       * most closely similar profile positions and forcing the child profile
       * to use the parent profile at those positions?
       *
       * @see shareProfilePositions_percentChangeTo0Threshold
       */
GALOSH_DEF_OPT(shareProfilePositions,bool,false,"After training a subfamily profile with all positions differing from its parent profile, should we try to retrain after identifying the most closely similar profile positions and forcing the child profile to use the parent profile at those positions??" );

      /**
       * After training a subfamily profile with all positions differing from
       * its parent profile, we try to retrain after identifying the most
       * closely similar profile positions and forcing the child profile to use
       * the parent profile at those positions.  We do this until the score
       * after removing a position drops below some threshold.
       *
       * This is the threshold for the percent by which the score has changed
       * from the original (all-positions-different) profile.
       *
       * @see shareProfilePositions
       */
GALOSH_DEF_OPT(shareProfilePositions_percentChangeTo0Threshold,double,-50.0,"After training a subfamily profile with all positions differing from its parent profile, we try to retrain after identifying the most closely similar profile positions and forcing the child profile to use the parent profile at those positions.  We do this until the score after removing a position drops below some threshold. This is the threshold for the percent by which the score has changed from the original (all-positions-different) profile.");

      /**
       * We have to make a decision about whether a sequence belongs to the
       * child profile or the parent profile before recursively breaking the
       * profile into more subfamilies.  If the sequence subfamily mixture
       * parameter exceeds this threshold, we say that the sequence belongs to
       * the child profile.
       */
GALOSH_DEF_OPT(childSequenceMixtureThreshold,double,0.5,"We have to make a decision about whether a sequence belongs to the child profile or the parent profile before recursively breaking the profile into more subfamilies.  If the sequence subfamily mixture parameter exceeds this threshold, we say that the sequence belongs to the child profile.");


/// This is how we're handling vectors.  It is a work-around because vectors are handled specially
/// by boost::program_options.  It allows the command line to look something like
///
/// --profileLengths 10 20 30
///
/// The TMP_EXTRA_STUFF must be set to include (at least) the ->multitoken() thing.
/// It should also be unset at the bottom of the vector initializations.


#undef TMP_EXTRA_STUFF
#define TMP_EXTRA_STUFF ->multitoken()
///Vector definitions go here.

/** do this after the vector definition section */
#undef TMP_EXTRA_STUFF
#define TMP_EXTRA_STUFF BOOST_PP_EMPTY()
