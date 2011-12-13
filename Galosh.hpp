/*---------------------------------------------------------------------------##
##  File:
##      @(#) Galosh.hpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      Global #defines for Galosh.
##
#******************************************************************************
#*  Copyright (C) 2008    Paul Edlefsen                                       *
#*  All rights reserved.                                                      *
#*****************************************************************************/

#if     _MSC_VER > 1000
#pragma once
#endif

#ifndef __GALOSH_GALOSH_HPP__
#define __GALOSH_GALOSH_HPP__

#if defined __MSC_VER
#include <boost/math/special_functions/log1p.hpp>
using boost::math::log1p;
#endif

// OLD
#define INIT_PROBABILITY(t) 

////// What shall we inline?
#define GALOSH_INLINE_INIT inline
#define GALOSH_INLINE_REINITIALIZE inline
#define GALOSH_INLINE_COPY inline
#define GALOSH_INLINE_ACCESSOR inline
#define GALOSH_INLINE_OSTREAM

// "trivial" means it is just a renaming of another fn (so all it does is call that other fn), or perhaps returning simple function of that call.
#define GALOSH_INLINE_TRIVIAL inline

// Probabilities
#ifdef GALOSH_USE_ICSILOG
  #define GALOSH_LOG(x) icsilogv2(x)
#else
  #define GALOSH_LOG(x) ::log(x)
#endif
#if defined __MSC_VER
#define GALOSH_LOG1P(x) boost::math::log1p(x)
#else
#define GALOSH_LOG1P(x) ::log1p(x)
#endif
#define GALOSH_INLINE_PROBABILITIES_ARITHMETIC inline
#define GALOSH_INLINE_PROBABILITIES_CAST inline
#define GALOSH_INLINE_PROBABILITIES_CAST_COMPARE inline
#define GALOSH_INLINE_PROBABILITIES_CAST_ARITHMETIC inline

// Algebra
#define GALOSH_INLINE_ALGEBRA_ARITHMETIC inline
#define GALOSH_INLINE_ALGEBRA_CAST inline
#define GALOSH_INLINE_ALGEBRA_CAST_COMPARE inline
#define GALOSH_INLINE_ALGEBRA_CAST_ARITHMETIC inline

// MultinomialDistribution
#define GALOSH_INLINE_MULTINOMIALDISTRIBUTION_ARITHMETIC inline
#define GALOSH_INLINE_MULTINOMIALDISTRIBUTION_COMPARE inline
#define GALOSH_INLINE_MULTINOMIALDISTRIBUTION_NORMALIZE_ENSUREMIN inline
#define GALOSH_INLINE_MULTINOMIALDISTRIBUTION_COMPLEX_ACCESSOR inline
#define GALOSH_INLINE_MULTINOMIALDISTRIBUTION_EUCLIDEANDISTSQ inline
#define GALOSH_INLINE_MULTINOMIALDISTRIBUTION_DRAW inline
#define GALOSH_INLINE_MULTINOMIALDISTRIBUTION_FROMDIRICHLET inline
#define GALOSH_INLINE_MULTINOMIALDISTRIBUTION_TO_BOLTZMANNGIBBS inline
#define GALOSH_INLINE_MULTINOMIALDISTRIBUTION_FROM_BOLTZMANNGIBBS inline

// Parameters, DynamicProgramming, ProfileTrainer, ProfileTreeTrainer, GaloshTest
#define GALOSH_INLINE_PARAMETERSMODIFIER_APPLY_MODIFICATIONS inline

// Profile
#define GALOSH_INLINE_PROFILE_ARITHMETIC inline
#define GALOSH_INLINE_PROFILE_EUCLIDEANDISTSQ inline
#define GALOSH_INLINE_PROFILE_CROSSENTROPY inline
#define GALOSH_INLINE_PROFILE_FREEPARAMCOUNT inline
#define GALOSH_INLINE_PROFILE_COMPLEX_ACCESSOR inline

// ProfileTrainer
#define GALOSH_INLINE_TRAIN
#define GALOSH_INLINE_TRAIN_PERCENT_CHANGE inline

// ProfileGibbs
#define GALOSH_INLINE_GIBBS

// DynamicProgramming
#define GALOSH_INLINE_ALGORITHM inline
#define GALOSH_INLINE_ALGORITHM_INNERLOOP inline

// ProfuseTest
#define GALOSH_INLINE_PROFUSETEST_START

#endif // __GALOSH_GALOSH_HPP__
