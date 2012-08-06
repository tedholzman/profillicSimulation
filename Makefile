#========================================
# USER RESERVED -- The following are reserved for users to set on the
# command line.  Makefiles should not set these.  These variables are
# for C/C++ compilation, and linking.
## -pg is for use with gprof.  Use on CFLAGS and on LDFLAGS
## note: -pg produces a gmon.out file, which is medium-sized and confusing.
## It is just as well not to use -pg if you're not using gmon.
#CFLAGS		= -Wall -pg
#CFLAGS         = -O3 -funroll-loops -Winline -DNDEBUG=1 \
#                 --param max-inline-insns-single=10000 --param inline-unit-growth=500 --param large-function-growth=1000
CFLAGS	        = -g3 -gdwarf-2
#JFLAGS		=
#LDFLAGS	= -pg

# OPTIMIZE with the -O option.  Override from the command line for
# building debug versions.
#
#OPTFLAGS	= -O

#========================================
### Paul added this for muscle support
# You must uncomment these for muscle support, and note that for now it changes the behavior of ProfileTreeTrainer.hpp from being a non-tree-trainer (just like ProfileTrainer.hpp) to doing "prolific" multiple alignment; this affects both seqantest and profusetest, so watch out!
#MUSCLE_CPPSRC = $(sort $(wildcard muscle/*.cpp))
#MUSCLE_CPPOBJ_TMP = $(subst .cpp,.o,$(MUSCLE_CPPSRC))
#MUSCLE_CPPOBJ = $(subst muscle/main.o,,$(MUSCLE_CPPOBJ_TMP)) ;
#MUSCLE_CFLAGS = -D__HAVE_MUSCLE
#MUSCLE_LDFLAGS =
#MUSCLE_LIBS =

## TODO: Why doesn't this work?
#
#if 1; then \
#  MUSCLE_CPPSRC = $(sort $(wildcard muscle/*.cpp)) \
#  MUSCLE_CPPOBJ_TMP = $(subst .cpp,.o,$(MUSCLE_CPPSRC)) \
#  MUSCLE_CPPOBJ = $(subst muscle/main.o,,$(MUSCLE_CPPOBJ_TMP)) ; \
#  MUSCLE_CFLAGS = -D__HAVE_MUSCLE \
#  MUSCLE_LDFLAGS = \
#  MUSCLE_LIBS	= \
#else \
#  MUSCLE_CPPSRC = \
#  MUSCLE_CPPOBJ = ; \
#  MUSCLE_CFLAGS = \
#  MUSCLE_LDFLAGS = \
#  MUSCLE_LIBS	= \
#fi;


#========================================
### Paul added these for HMMer (and squid) support
HMMER_CFLAGS 	= -I./hmmer/src -I./hmmer/easel
HMMER_LDFLAGS 	= -L./hmmer/src -L./hmmer/easel
HMMER_LIBS	= $(if $(wildcard hmmer/libsquid.*),-lsquid,)  -lhmmer  # TAH temporary: why isn't squid built on hmmer on linux?
#HMMER_LIBS	=  -lhmmer 
#HMMER_CFLAGS 	=
#HMMER_LDFLAGS 	=
#HMMER_LIBS	=

#========================================
### Paul added these for ALGLIB support
ALGLIB_CPPSRC = $(sort $(wildcard alglib/*.cpp))
ALGLIB_CPPOBJ_TMP = $(subst .cpp,.o,$(ALGLIB_CPPSRC))
ALGLIB_CPPOBJ = $(subst alglib/main.o,,$(ALGLIB_CPPOBJ_TMP)) ;
ALGLIB_CFLAGS =
ALGLIB_LDFLAGS =
ALGLIB_LIBS =

#========================================
### Paul added these for Seqan support
SEQAN_CFLAGS 	= -I./seqan
SEQAN_LDFLAGS 	=
SEQAN_LIBS	=
#SEQAN_CFLAGS 	=
#SEQAN_LDFLAGS =
#SEQAN_LIBS	=

#========================================
### Paul added these for Boost support
### For Boost, eg "g++ -I/sw/include/boost-1_44_0/ profile_xml.cpp DNAResidue.o AminoAcidResidue.o DSSPResidue.o Profile.o -o profile_xml -L/sw/lib/ -lboost_serialization":
BOOST_CFLAGS 	= -I./boost-include
BOOST_LDFLAGS 	= -L./boost-lib
BOOST_LIBS	= -lboost_serialization -lboost_graph -lboost_filesystem -lboost_system -lboost_regex \
                   -lboost_program_options \
#                   -lboost_unit_test_framework -lboost_test_exec_monitor

###==============================================
### TAH 11/3/2011
### Added for modularity.  HMMOC, PROLIFIC, etc. can stay in their respective directories
HMMOC_LIB       = ./HMMoC-BFloat-Algebra
PROLIFIC_LIB    = ./prolific
PROFILLIC_LIB   = ./profillic
PROFUSE_LIB     = ./profuse

HMMOC_CFLAGS  = -I$(HMMOC_LIB)
PROLIFIC_CFLAGS = -I$(PROLIFIC_LIB)
PROFILLIC_CFLAGS = -I$(PROFILLIC_LIB)
PROFUSE_CFLAGS = -I$(PROFUSE_LIB)
.PHONY : $(HMMOC_LIB)/Algebra.o

all: tests progs converters

all-with-hmmer: all converters-with-hmmer tests-with-hmmer


###==============================================

INCS = Galosh.hpp \
$(PROLIFIC_LIB)/Ambiguous.hpp \
$(PROLIFIC_LIB)/Sequence.hpp \
$(PROLIFIC_LIB)/MultinomialDistribution.hpp \
$(PROLIFIC_LIB)/ProfileHMM.hpp \
$(PROLIFIC_LIB)/Profile.hpp \
$(PROLIFIC_LIB)/ProfileTree.hpp \
$(PROLIFIC_LIB)/Fasta.hpp \
$(HMMOC_LIB)/Algebra.hpp \
$(PROLIFIC_LIB)/Random.hpp \
$(PROLIFIC_LIB)/DynamicProgramming.hpp \
$(PROFILLIC_LIB)/ProfileTrainer.hpp \
ProfileTreeTrainer.hpp \
ProfileGibbs.hpp \
ProfuseTest.hpp \
ProfuseTest2.hpp \
Parameters.hpp \
CommandlineParameters.hpp \
ProfuseTestOptions.hpp \
$(PROFUSE_LIB)/Profuse.hpp \
$(PROLIFIC_LIB)/AminoAcid20.hpp

TRAIN_INCS = $(INCS) \
$(PROFUSE_LIB)/Profuse.hpp

ALIGN_INCS = $(INCS) \
$(PROFUSE_LIB)/Profuse.hpp

SCORE_INCS = $(INCS) \
$(PROFUSE_LIB)/Profuse.hpp

CREATERANDOMSEQUENCE_INCS = $(PROLIFIC_LIB)/Fasta.hpp \
$(PROLIFIC_LIB)/Random.hpp \
$(PROLIFIC_LIB)/MultinomialDistribution.hpp \
$(HMMOC_LIB)/Algebra.hpp \
$(PROLIFIC_LIB)/Sequence.hpp

DRAWSEQUENCES_INCS = $(HMMOC_LIB)/Algebra.hpp \
$(PROLIFIC_LIB)/Random.hpp \
$(PROLIFIC_LIB)/MultinomialDistribution.hpp \
$(PROLIFIC_LIB)/DynamicProgramming.hpp \
$(PROLIFIC_LIB)/Sequence.hpp \
$(PROLIFIC_LIB)/Fasta.hpp \
$(PROLIFIC_LIB)/Profile.hpp

EVALUATEENTENTES_INCS = $(HMMOC_LIB)/Algebra.hpp \
$(PROLIFIC_LIB)/Random.hpp \
$(PROLIFIC_LIB)/MultinomialDistribution.hpp \
$(PROLIFIC_LIB)/DynamicProgramming.hpp \
$(PROLIFIC_LIB)/Sequence.hpp \
$(PROLIFIC_LIB)/Fasta.hpp \
$(PROLIFIC_LIB)/Profile.hpp \
$(PROLIFIC_LIB)/ProfileTree.hpp \
ProfileTreeTrainer.hpp

PROFILE2HMMER_INCS = $(HMMOC_LIB)/Algebra.hpp \
$(PROLIFIC_LIB)/Random.hpp \
$(PROLIFIC_LIB)/Ambiguous.hpp \
$(PROLIFIC_LIB)/MultinomialDistribution.hpp \
$(PROLIFIC_LIB)/ProfileHMM.hpp \
$(PROLIFIC_LIB)/Sequence.hpp \
$(PROLIFIC_LIB)/Profile.hpp

PROFILE2FASTA_INCS = $(HMMOC_LIB)/Algebra.hpp \
$(PROLIFIC_LIB)/Ambiguous.hpp \
$(PROLIFIC_LIB)/Random.hpp \
$(PROLIFIC_LIB)/MultinomialDistribution.hpp \
$(PROLIFIC_LIB)/ProfileHMM.hpp \
$(PROLIFIC_LIB)/Sequence.hpp \
$(PROLIFIC_LIB)/Fasta.hpp \
$(PROLIFIC_LIB)/Profile.hpp

SEQUENCE2PROFILE_INCS = $(PROFUSETEST_INCS)

ALIGNEDFASTA2PROFILE_INCS = $(PROFUSETEST_INCS)

XML2PROFILE_INCS = $(PROFUSETEST_INCS)

QUICKTEST_INCS = Galosh.hpp \
$(HMMOC_LIB)/Algebra.hpp \
Parameters.hpp

SEQANTEST_INCS =  $(INCS)

PROFUSETEST_INCS = $(INCS) \
ProfuseTest.hpp ProfuseTest2.hpp

MUSCLEHMMERTEST_INCS = $(INCS)

OBJS = $(HMMOC_LIB)/Algebra.o

ALIGN_OBJS = $(OBJS) \
$(PROFUSE_LIB)/Align.o

SCORE_OBJS = $(OBJS) \
$(PROFUSE_LIB)/Score.o

CREATERANDOMSEQUENCE_OBJS = $(OBJS) \
$(PROFUSE_LIB)/CreateRandomSequence.o

DRAWSEQUENCES_OBJS = $(OBJS) \
$(PROFUSE_LIB)/DrawSequences.o

EVALUATEENTENTES_OBJS = $(OBJS) \
EvaluateEntentes.o

PROFILE2HMMER_OBJS = $(OBJS) \
Profile2HMMer.o

PROFILE2FASTA_OBJS = $(OBJS) \
Profile2Fasta.o

SEQUENCE2PROFILE_OBJS = $(OBJS) \
Sequence2Profile.o

ALIGNEDFASTA2PROFILE_OBJS = $(OBJS) \
AlignedFasta2Profile.o

XML2PROFILE_OBJS = $(OBJS) \
XML2Profile.o

PROFUSETEST_OBJS = $(OBJS) \
ProfuseTest.o

PROFUSETEST2_OBJS = $(OBJS) \
ProfuseTest2.o

MUSCLEHMMERTEST_OBJS = $(OBJS) \
MuscleHMMerTests.o

QUICKTEST_OBJS = $(HMMOC_LIB)/Algebra.o \
QuickTests.o

SEQANTEST_OBJS = $(OBJS) \
SeqanTests.o

NEWSEQANTEST_OBJS = $(OBJS) \
NewSeqanTests.o

SOURCES = $(HMMOC_LIB)/Algebra.cpp

ALIGN_SOURCES = $(SOURCES) \
$(PROFUSE_LIB)/Align.cpp

SCORE_SOURCES = $(SOURCES) \
$(PROFUSE_LIB)/Score.cpp

PROFUSETEST_SOURCES = $(SOURCES) \
ProfuseTest.cpp

PROFUSETEST2_SOURCES = $(SOURCES) \
ProfuseTest2.cpp

PROFILE2HMMER_SOURCES = $(SOURCES) \
Profile2HMMer.cpp

CREATERANDOMSEQUENCE_SOURCES = $(SOURCES) \
$(PROFUSE_LIB)/CreateRandomSequence.cpp

DRAWSEQUENCES_SOURCES = $(SOURCES) \
$(PROFUSE_LIB)/DrawSequences.cpp

EVALUATEENTENTES_SOURCES = $(SOURCES) \
EvaluateEntentes.cpp

PROFILE2FASTA_SOURCES = $(SOURCES) \
Profile2Fasta.cpp

SEQUENCE2PROFILE_SOURCES = $(SOURCES) \
Sequence2Profile.cpp

ALIGNEDFASTA2PROFILE_SOURCES = $(SOURCES) \
AlignedFasta2Profile.cpp

XML2PROFILE_SOURCES = $(SOURCES) \
XML2Profile.cpp

MUSCLEHMMERTEST_SOURCES = $(SOURCES) \
MuscleHMMerTests.cpp

QUICKTEST_SOURCES = $(HMMOC_LIB)/Algebra.cpp \
QuickTests.cpp

SEQANTEST_SOURCES = $(HMMOC_LIB)/Algebra.cpp \
SeqanTests.cpp

NEWSEQANTEST_SOURCES = $(HMMOC_LIB)/Algebra.cpp \
NewSeqanTests.cpp

default: quicktest

align: $(ALIGN_SOURCES) $(ALIGN_INCS) $(ALIGN_OBJS) $(MUSCLE_CPPOBJ)
	     $(CXX_LINK) -o align $(ALIGN_OBJS) $(MUSCLE_CPPOBJ)

score: $(SCORE_SOURCES) $(SCORE_INCS) $(SCORE_OBJS) $(MUSCLE_CPPOBJ)
	     $(CXX_LINK) -o score $(SCORE_OBJS) $(MUSCLE_CPPOBJ)

createRandomSequence: $(CREATERANDOMSEQUENCE_SOURCES) $(CREATERANDOMSEQUENCE_INCS) $(CREATERANDOMSEQUENCE_OBJS)
	     $(CXX_LINK) -o createRandomSequence $(CREATERANDOMSEQUENCE_OBJS)

drawSequences: $(DRAWSEQUENCES_SOURCES) $(DRAWSEQUENCES_INCS) $(DRAWSEQUENCES_OBJS)
	     $(CXX_LINK) -o drawSequences $(DRAWSEQUENCES_OBJS)

evaluateEntentes: $(EVALUATEENTENTES_SOURCES) $(EVALUATEENTENTES_INCS) $(EVALUATEENTENTES_OBJS) $(ALGLIB_CPPOBJ)
	     $(CXX_LINK) -o evaluateEntentes $(EVALUATEENTENTES_OBJS) $(ALGLIB_CPPOBJ)

profile2hmmer: $(PROFILE2HMMER_SOURCES) $(PROFILE2HMMER_INCS) $(PROFILE2HMMER_OBJS) $(MUSCLE_CPPOBJ)
	     $(CXX_LINK) -o profile2hmmer $(PROFILE2HMMER_OBJS) $(MUSCLE_CPPOBJ)

profile2fasta: $(PROFILE2FASTA_SOURCES) $(PROFILE2FASTA_INCS) $(PROFILE2FASTA_OBJS) $(MUSCLE_CPPOBJ)
	     $(CXX_LINK) -o profile2fasta $(PROFILE2FASTA_OBJS) $(MUSCLE_CPPOBJ)

sequence2profile: $(SEQUENCE2PROFILE_SOURCES) $(SEQUENCE2PROFILE_INCS) $(SEQUENCE2PROFILE_OBJS) $(MUSCLE_CPPOBJ)
	     $(CXX_LINK) -o sequence2profile $(SEQUENCE2PROFILE_OBJS) $(MUSCLE_CPPOBJ)

alignedFasta2profile: $(ALIGNEDFASTA2PROFILE_SOURCES) $(ALIGNEDFASTA2PROFILE_INCS) $(ALIGNEDFASTA2PROFILE_OBJS) $(MUSCLE_CPPOBJ)
	     $(CXX_LINK) -o alignedFasta2profile $(ALIGNEDFASTA2PROFILE_OBJS) $(MUSCLE_CPPOBJ)

xml2profile: $(XML2PROFILE_SOURCES) $(XML2PROFILE_INCS) $(XML2PROFILE_OBJS) $(MUSCLE_CPPOBJ)
	     $(CXX_LINK) -o xml2profile $(XML2PROFILE_OBJS)

quicktest: $(QUICKTEST_SOURCES) $(QUICKTEST_INCS) $(QUICKTEST_OBJS)
	   $(CXX_LINK) -o quicktest $(QUICKTEST_OBJS)

profusetest: $(PROFUSETEST_SOURCES) $(PROFUSETEST_INCS) $(PROFUSETEST_OBJS) $(MUSCLE_CPPOBJ)
	     $(CXX_LINK) -o profusetest $(PROFUSETEST_OBJS) $(MUSCLE_CPPOBJ)

profusetest2: $(PROFUSETEST2_SOURCES) $(PROFUSETEST_INCS) $(PROFUSETEST2_OBJS) $(MUSCLE_CPPOBJ)
	     $(CXX_LINK) -o profusetest2 $(PROFUSETEST2_OBJS) $(MUSCLE_CPPOBJ)

seqantest: $(SEQANTEST_SOURCES) $(SEQANTEST_INCS) $(SEQANTEST_OBJS) $(MUSCLE_CPPOBJ)
	   $(CXX_LINK) -o seqantest $(SEQANTEST_OBJS) $(MUSCLE_CPPOBJ)

newseqantest: $(NEWSEQANTEST_SOURCES) $(SEQANTEST_INCS) $(NEWSEQANTEST_OBJS) $(MUSCLE_CPPOBJ) 	   
	   $(CXX_LINK) -o newseqantest $(NEWSEQANTEST_OBJS) $(MUSCLE_CPPOBJ)
	   ##$(DOXY_CMD)        ## don't do this with recursive doxygen options turned on, at least until you become immortal and very patient  TAH

musclehmmertest: $(MUSCLEHMMERTEST_SOURCES) $(MUSCLEHMMERTEST_INCS) $(MUSCLEHMMERTEST_OBJS) $(MUSCLE_CPPOBJ)
	     $(CXX_LINK) -o musclehmmertest $(MUSCLEHMMERTEST_OBJS) $(MUSCLE_CPPOBJ)

tests: quicktest seqantest profusetest newseqantest profusetest2

tests-with-hmmer: tests musclehmmertest

converters: sequence2profile alignedFasta2profile xml2profile profile2fasta

converters-with-hmmer: converters profile2hmmer 

progs:  align score createRandomSequence drawSequences evaluateEntentes

## Recompile if the includes are modified ...
$(PROFUSETEST_OBJS): $(PROFUSETEST_SOURCES) $(PROFUSETEST_INCS)
$(PROFUSETEST2_OBJS): $(PROFUSETEST2_SOURCES) $(PROFUSETEST_INCS)
$(PROFILE2HMMER_OBJS): $(PROFILE2HMMER_SOURCES) $(PROFILE2HMMER_INCS)
$(CREATERANDOMSEQUENCE_OBJS): $(CREATERANDOMSEQUENCE_SOURCES) $(CREATERANDOMSEQUENCE_INCS)
$(DRAWSEQUENCES_OBJS): $(DRAWSEQUENCES_SOURCES) $(DRAWSEQUENCES_INCS)
$(EVALUATEENTENTES_OBJS): $(EVALUATEENTENTES_SOURCES) $(EVALUATEENTENTES_INCS)
$(PROFILE2FASTA_OBJS): $(PROFILE2FASTA_SOURCES) $(PROFILE2FASTA_INCS)
$(SEQUENCE2PROFILE_OBJS): $(SEQUENCE2PROFILE_SOURCES) $(SEQUENCE2PROFILE_INCS)
$(ALIGNEDFASTA2PROFILE_OBJS): $(ALIGNEDFASTA2PROFILE_SOURCES) $(ALIGNEDFASTA2PROFILE_INCS)
$(XML2PROFILE_OBJS): $(XML2PROFILE_SOURCES) $(XML2PROFILE_INCS)
$(ALIGN_OBJS): $(ALIGN_SOURCES) $(ALIGN_INCS)
$(SCORE_OBJS): $(SCORE_SOURCES) $(SCORE_INCS)
$(MUSCLEHMMERTEST_OBJS): $(MUSCLEHMMERTEST_SOURCES) $(MUSCLEHMMERTEST_INCS)
$(QUICKTEST_OBJS): $(QUICKTEST_SOURCES) $(QUICKTEST_INCS)
$(SEQANTEST_OBJS): $(SEQANTEST_SOURCES) $(SEQANTEST_INCS)
$(NEWSEQANTEST_OBJS): $(NEWSEQANTEST_SOURCES) $(NEWSEQANTEST_INCS)
.PHONY: clean
clean:
	rm -f $(notdir align score createRandomSequence drawSequences evaluateEntentes profile2hmmer profile2fasta sequence2profile alignedFasta2profile xml2profile quicktest seqantest profusetest musclehmmertest newseqantest profusetest2 $(TRAIN_OBJS) $(ALIGN_OBJS) $(SCORE_OBJS) $(PROFILE2HMMER_OBJS) $(CREATERANDOMSEQUENCE_OBJS) $(DRAWSEQUENCES_OBJS) $(EVALUATEENTENTES_OBJS) $(PROFILE2HMMER_OBJS) $(PROFILE2FASTA_OBJS) $(SEQUENCE2PROFILE_OBJS) $(XML2PROFILE_OBJS) $(ALIGNEDFASTA2PROFILE_OBJS) $(QUICKTEST_OBJS) $(SEQANTEST_OBJS) $(NEWSEQANTEST_OBJS) $(PROFUSETEST_OBJS) $(MUSCLEHMMERTEST_OBJS) $(PROFUSETEST2_OBJS))

#========================================
# FILE EXTENSIONS.  Extensions and prefixes for different types of
# files change from platform to platform.  Hide these in macros so
# that we can more easily cut and paste between makefiles.
o		= .o
EXE_SFX		= 
SCRIPT_SFX 	= 
LIB_PFX		= lib
LIB_SFX		= .a
LIB_SHARED_SFX	= .so
TMPLIB		= libtemp.a

# FILE TOOLS
AR 	= ar qv
CHMOD 	= chmod
CP	= cp
GREP	= grep
MKDIR 	= mkdir
MUNCH 	= stepmunch
MV	= mv
NM 	= nm
RANLIB	= ranlib
RM 	= rm -f
RMDIR 	= rm -rf
STRIP	= strip
UNZIP 	= unzip
ZIP 	= zip


#========================================
# ANSI C Compile and Link
#
CC		= gcc
CC_COMPILE	= $(CC) -c $(OPTFLAGS) $(CFLAGS) $(CC_CFLAGS) $(CC_SYSCFLAGS)
CC_LINK		= $(CC) $(LDFLAGS) $(CC_LDFLAGS) $(CC_SYSLDFLAGS) $(CC_LIBS)
CC_CFLAGS 	= $(BOOST_CFLAGS) $(SEQAN_CFLAGS) $(ALGLIB_CFLAGS) $(HMMER_CFLAGS) $(MUSCLE_CFLAGS)
CC_LDFLAGS	= $(BOOST_LDFLAGS) $(SEQAN_LDFLAGS) $(ALGLIB_LDFLAGS) $(HMMER_LDFLAGS) $(MUSCLE_LDFLAGS)
CC_LIBS		= $(BOOST_LIBS) $(SEQAN_CLIBS) $(ALGLIB_CLIBS) $(HMMER_CLIBS) $(MUSCLE_CLIBS)

# Global system things used for compilation, static linking, etc.
CC_SYSCFLAGS 	= -I.
CC_SYSLDFLAGS 	=
CC_SYSLIBS	=

#========================================
# C++ Compile and Link
#
CXX		= g++
CXX_COMPILE	= $(CXX) -c  $(OPTFLAGS) $(CFLAGS) $(CXX_CFLAGS) $(CXX_SYSCFLAGS)
CXX_LINK	= $(CXX) $(LDFLAGS) $(CXX_LDFLAGS) $(CXX_SYSLDFLAGS) $(CXX_LIBS)
CXX_CFLAGS 	= $(BOOST_CFLAGS) $(SEQAN_CFLAGS) $(ALGLIB_CFLAGS) $(HMMER_CFLAGS) $(MUSCLE_CFLAGS) \
                  $(HMMOC_CFLAGS) $(PROLIFIC_CFLAGS) $(PROFILLIC_CFLAGS) $(PROFUSE_CFLAGS)
CXX_LDFLAGS	= $(BOOST_LDFLAGS) $(SEQAN_LDFLAGS) $(ALGLIB_LDFLAGS) $(HMMER_LDFLAGS) $(MUSCLE_LDFLAGS)
CXX_LIBS	= $(BOOST_LIBS) $(SEQAN_LIBS) $(ALGLIB_LIBS) $(HMMER_LIBS) $(MUSCLE_LIBS)

#========================================
# Doxygen build
#
DOXY_LOC	= /Users/Ted/Development/doxygen-svn/bin/doxygen 
DOXY_CMD	= $(DOXY_LOC) 

# The force flags are used for C/C++ compilers that select the
# language based on the file naming conventions.  Some C++ source
# may be in files with C naming conventions.
CXX_FORCE	= 

# System Flags -- Things for static linking or making sure that the
# compiler understands that a file is a C++ file or whatever.  These
# usually change from platform to platform.
CXX_SYSCFLAGS 	= -I. 
CXX_SYSLDFLAGS 	= 
CXX_SYSLIBS	= 

# Compilation Rules -- Repeat the rules for all of the different
# naming conventions.
#
.cxx.o:	; $(CXX_COMPILE) $<
.cpp.o:	; $(CXX_COMPILE) $<
.cc.o:	; $(CXX_COMPILE) $<
.C.o:	; $(CXX_COMPILE) $<

.cxx:	
	$(CXX_COMPILE) $<
	$(CXX_LINK) -o $@ $*.o $(LIBRARIES)
.cpp:	
	$(CXX_COMPILE) $<
	$(CXX_LINK) -o $@ $*.o $(LIBRARIES)
.cc:	
	$(CXX_COMPILE) $<
	$(CXX_LINK) -o $@ $*.o $(LIBRARIES)
.C:	
	$(CXX_COMPILE) $<
	$(CXX_LINK) -o $@ $*.o $(LIBRARIES)

# for legacy reasons also compile .c as c++
.c.o:	; $(CXX_COMPILE) $(CXX_FORCE) $<
.c:	
	$(CXX_COMPILE) $(CXX_FORCE) $<

