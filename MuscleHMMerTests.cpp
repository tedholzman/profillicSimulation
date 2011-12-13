#include "ProfuseTest.hpp"

#ifdef __HAVE_MUSCLE


// For muscle support
#ifdef	WIN32
#include <windows.h>	// for SetPriorityClass()
#include <io.h>			// for isatty()
#else
#include <unistd.h>		// for isatty()
#endif
#include "muscle/muscle.h"
#include "muscle/seq.h"
#include "muscle/msa.h"
#include "muscle/textfile.h"
#include "muscle/alpha.h"
#include "muscle/profile.h"
#include "muscle/tree.h"

// For HMMer support
namespace hmmer {
extern "C" {
#include "hmmer/src/config.h"		/* compile-time configuration constants */
#include "hmmer/squid/squidconf.h"
#include "hmmer/src/structs.h"		/* data structures, macros, #define's   */
#include "hmmer/src/funcs.h"		/* function declarations                */
#include "hmmer/src/globals.h"		/* alphabet global variables            */
#include "hmmer/squid/squid.h"		/* general sequence analysis library    */
#include "hmmer/squid/msa.h"            /* squid's multiple alignment i/o       */
} // End extern "C"
} // End namespace hmmer

#include <sstream> // For stringstream.
#include <ctime> // for std::time

#include <iostream> // for ofstream stuff
#include <fstream> // for ofstream stuff

/////////////
/// Copied from http://snipplr.com/view.php?codeview&id=1790
#include <errno.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>

namespace galosh {
pid_t
spawn ( const char * const file, const char * const argv[] )
{
    pid_t pid;
    int forkpipe[2];
    int fd_flags, err, ret;

    // set up a pipe for communication between forks
    ret = pipe(forkpipe);
    if (ret == -1)
        return -1;

    pid = fork();
    if (pid == -1) {
        err = errno;

        // close the write pipe
        close(forkpipe[1]);

        goto parent_exit;
    }

    if (pid == 0) {
        // forked child

        // close the read pipe
        ret = close(forkpipe[0]);
        if (ret == -1)
            goto child_err;

        // make the write end close-on-exec
        ret = fcntl(forkpipe[1], F_GETFD);
        if (ret == -1)
            goto child_err;
        fd_flags = ret | FD_CLOEXEC;
        ret = fcntl(forkpipe[1], F_SETFD, fd_flags);
        if (ret == -1)
            goto child_err;

        // try to run the app
        // this const_cast ought to be safe.
        execvp(file, const_cast<char* const*>( argv ) );
        //execvp(file, argv);
        // if we continue executing here, an error occurred

child_err:
        // send the error code through the write pipe
        err = errno;
        write(forkpipe[1], &err, sizeof(err));

        exit(1);
    }

    // parent

    // close the write pipe
    close(forkpipe[1]);

    // get the error code from the read pipe
    ret = read(forkpipe[0], &err, sizeof(err));
    if (ret == -1) {
        err = errno;
        pid = -1;
    } else if (ret == sizeof(err))
        pid = -1;
    else
        err = 0;

parent_exit:
    // close the read pipe
    close(forkpipe[0]);

    if (err)
        errno = err;

    return pid;
} // spawn( const char * const, const char * const * )

}
/////////////

// TODO: REMOVE.  TESTING.
//profile_dna_distribution_tag const Profile::Emission::DNA();


/////////////
/// Copied from muscle's writescorefile.cpp and modified
extern float VTML_SP[32][32];
extern float NUC_SP[32][32];

/**
 * Detect insertion columns by detecting lowercase letters.
 */
double getColScoreAndDetectInsertions ( const ::MSA &msa, unsigned uCol, bool & is_insertion )
{
        is_insertion = false;
	const unsigned uSeqCount = msa.GetSeqCount();
	unsigned uPairCount = 0;
	double dSum = 0.0;
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		if (msa.IsGap(uSeq1, uCol))
			continue;
		unsigned uLetter1 = msa.GetLetterEx(uSeq1, uCol);
		if (uLetter1 >= g_AlphaSize)
			continue;
                if( islower( msa.GetChar(uSeq1, uCol) ) ) {
                  is_insertion = true;
                  // All letters in the column will be lowercase, so this
                  // suffices.
                }
		for (unsigned uSeq2 = uSeq1 + 1; uSeq2 < uSeqCount; ++uSeq2)
			{
			if (msa.IsGap(uSeq2, uCol))
				continue;
			unsigned uLetter2 = msa.GetLetterEx(uSeq2, uCol);
			if (uLetter2 >= g_AlphaSize)
				continue;
			double Score;
			switch (g_Alpha)
				{
			case ALPHA_Amino:
				Score = VTML_SP[uLetter1][uLetter2];
				break;
			case ALPHA_DNA:
			case ALPHA_RNA:
				Score = NUC_SP[uLetter1][uLetter2];
				break;
			default:
				Quit("GetColScore: invalid alpha=%d", g_Alpha);
				}
			dSum += Score;
			++uPairCount;
			}
		}
	if (0 == uPairCount)
		return 0;
	return dSum / uPairCount;
	}
// Note there is no sequence weighting.
double
getNonInsertionSPScore ( const MSA &msa )
{
	const unsigned uColCount = msa.GetColCount();
	const unsigned uSeqCount = msa.GetSeqCount();
        bool is_insertion = false;
        double total_score = 0;
	for (unsigned uCol = 0; uCol < uColCount; ++uCol)
		{
                  double Score = getColScoreAndDetectInsertions(msa, uCol, is_insertion);
                  if( !is_insertion ) {
                    total_score += Score;
                  }
		}
        return total_score;
} // getNonInsertionSPScore( const & MSA )
/////////////

/////////////
/// Copied from (prefab-4's) qscore's qscore.cpp
static void ToUpper ( MSA &msa )
	{
	const int SeqCount = msa.GetSeqCount();
	const int ColCount = msa.GetColCount();

	for (int SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		for (int ColIndex = 0; ColIndex < ColCount; ++ColIndex)
			{
			char c = msa.GetChar(SeqIndex, ColIndex);
			if (isalpha(c))
				{
				c = toupper(c);
				msa.SetChar(SeqIndex, ColIndex, c);
				}
			}
	}
/////////////

/////////////
/// Copied from (prefab-4's) qscore's qscore.h
inline int iabs (int i)
	{
	return i >= 0 ? i : -i;
	}

/////////////
/// Copied from (prefab-4's) qscore's sumpairs.cpp
// Compute the sum of pairs score from two pair maps.
// It is the simplicity of this code that motivates the use of pair maps.
double SumPairs (const int iMapRef[], const int iMapTest[], unsigned uLength)
	{
	unsigned uPairCount = 0;
	unsigned uCorrectPairCount = 0;
	for (unsigned uPos = 0; uPos < uLength; ++uPos)
		{
		int iPosRef = iMapRef[uPos];
		if (-1 == iPosRef)
			continue;

		++uPairCount;

		int iPosTest = iMapTest[uPos];
		if (-1 == iPosTest)
			continue;

		if (iPosRef == iPosTest)
			++uCorrectPairCount;
		}

	if (0 == uPairCount)
		return 0.0;

	return (double) uCorrectPairCount / (double) uPairCount;
	}
/////////////

/////////////
/// Copied from (prefab-4's) qscore's clineshift.cpp
// Compute the Cline shift score [1] from two pair maps.
// It is the relative simplicity of this code that motivates the use of pair maps.
// [1] Cline, Hughey & Karplus (2002), Bioinformatics 18(2) p.306.
double ClineShift (const int iTestMapA[], const int iRefMapA[], unsigned uLengthA,
  const int iTestMapB[], const int iRefMapB[], unsigned uLengthB, double dEpsilon = 0.2)
	{
	unsigned uRefPairCount = 0;
	unsigned uTestPairCount = 0;
	double dTotal = 0.0;

	for (unsigned uPosA = 0; uPosA < uLengthA; ++uPosA)
		{
		int iRefPosB = iRefMapA[uPosA];
		if (-1 == iRefPosB)
			continue;

		++uRefPairCount;

		int iTestPosB = iTestMapA[uPosA];
		if (-1 == iTestPosB)
			continue;

		int iShift = iabs(iRefPosB - iTestPosB);
		double dScore = (1 + dEpsilon)/(1 + iShift) - dEpsilon;
		assert(dScore >= -dEpsilon);
		assert(dScore <= 1.0);

		dTotal += dScore;
		}

	for (unsigned uPosB = 0; uPosB < uLengthB; ++uPosB)
		{
		int iTestPosA = iTestMapB[uPosB];
		if (-1 == iTestPosA)
			continue;

		++uTestPairCount;

		int iRefPosA = iRefMapB[uPosB];
		if (-1 == iRefPosA)
			continue;

		int iShift = iabs(iRefPosA - iTestPosA);
		double dScore = (1 + dEpsilon)/(1 + iShift) - dEpsilon;
		assert(dScore >= -dEpsilon);
		assert(dScore <= 1.0);

		dTotal += dScore;
		}

	if (0 == uRefPairCount)
		Quit("ClineShift: No aligned pair in ref alignment");

	return dTotal / (double) (uTestPairCount + uRefPairCount);
	}
/////////////

/////////////
/// Copied from (prefab-4's) qscore's comparemap.cpp
static bool LocalEq (double x, double y)
	{
	double diff = fabs(x-y);
	return diff < 1e-3;
	}

double ComparePairMapSP (const int iTestMapA[], const int iTestMapB[],
  const int iRefMapA[], const int iRefMapB[], int iLengthA, int iLengthB)
	{
	double dSP = SumPairs(iRefMapA, iTestMapA, iLengthA);

#if	_DEBUG
	{
// Verify symmetry as a correctness check
	double dSP2 = SumPairs(iRefMapB, iTestMapB, iLengthB);
	if (!LocalEq(dSP, dSP2))
		Quit("ComparePairMapSP: dSP=%g dSP2=%g diff=%g", dSP, dSP2, fabs(dSP-dSP2));
	}
#endif

	return dSP;
	}

void ComparePairMap (const int iTestMapA[], const int iTestMapB[],
  const int iRefMapA[], const int iRefMapB[], int iLengthA, int iLengthB,
  double *ptrdSP, double *ptrdPS, double *ptrdCS)
	{
	double dSP = SumPairs(iRefMapA, iTestMapA, iLengthA);
	double dPS = SumPairs(iTestMapA, iRefMapA, iLengthA);
	double dCS = ClineShift(iTestMapA, iRefMapA, iLengthA, iTestMapB,
	  iRefMapB, iLengthB);

#if	_DEBUG
	{
// Verify symmetries as a correctness check
	double dSP2 = SumPairs(iRefMapB, iTestMapB, iLengthB);
	double dPS2 = SumPairs(iTestMapB, iRefMapB, iLengthB);

	double dCS2 = ClineShift(iTestMapB, iRefMapB, iLengthB, iTestMapA,
	  iRefMapA, iLengthA);
	double dCS3 = ClineShift(iRefMapA, iTestMapA, iLengthA, iRefMapB,
	  iTestMapB, iLengthB);
	double dCS4 = ClineShift(iRefMapB, iTestMapB, iLengthB, iRefMapA,
	  iTestMapA, iLengthA);

	if (!LocalEq(dSP, dSP2))
		Quit("CompareSeqs: dSP=%g dSP2=%g diff=%g", dSP, dSP2, fabs(dSP-dSP2));
	if (!LocalEq(dPS, dPS2))
		Quit("CompareSeqs: dPS=%g dPS2=%g diff=%g", dPS, dPS2, fabs(dPS-dPS2));
	if (!LocalEq(dCS, dCS2))
		Quit("CompareSeqs: dCS=%g dCS2=%g diff=%g", dCS, dCS2, fabs(dCS-dCS2));
	if (!LocalEq(dCS, dCS3))
		Quit("CompareSeqs: dCS=%g dCS3=%g diff=%g", dCS, dCS3, fabs(dCS-dCS3));
	if (!LocalEq(dCS, dCS4))
		Quit("CompareSeqs: dCS=%g dCS4=%g diff=%g", dCS, dCS4, fabs(dCS-dCS4));
	}
#endif

	*ptrdSP = dSP;
	*ptrdPS = dPS;
	*ptrdCS = dCS;
	}
/////////////

/////////////
/// Copied from (prefab-4's) qscore's msa.cpp
bool IsGap (char c)
	{
	return '-' == c || '~' == c || '.' == c;
	}
/***
It is sometimes very convenient to represent a pairwise alignment
as a "pair map", which works as follows.

Let iPos1 be the index into ungapped sequence 1, similarly for iPos2.

Then if a pair of letters (iPos1, iPos2) is aligned:

	iMap1[iPos1] = iPos2 and iMap2[iPos2] = iPos1.

If iPos1 is not in an aligned column, or is aligned to a gap, then
iMap1[iPos1] = -1, and similarly for iMap2. This overloads the meaning
of the integer value, so is questionable software engineering practice;
however it's a simple and convenient solution for small applications.
***/
void MSA::GetPairMap (unsigned uSeqIndex1, unsigned uSeqIndex2, int iMap1[],
  int iMap2[]) const
	{
	assert(uSeqIndex1 < GetSeqCount());
	assert(uSeqIndex2 < GetSeqCount());

	int iPos1 = 0;
	int iPos2 = 0;
	const unsigned uColCount = GetColCount();
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		char c1 = GetChar(uSeqIndex1, uColIndex);
		char c2 = GetChar(uSeqIndex2, uColIndex);
		bool bIsGap1 = ::IsGap(c1);
		bool bIsGap2 = ::IsGap(c2);
		if (!bIsGap1 && !bIsGap2)
			{
			if (isupper(c1))
				{
				if (!isupper(c2))
					Quit("Both upper and lower case letters (%c,%c) in ref alignment column %d",
					  c1, c2, uColIndex);
				iMap1[iPos1] = iPos2;
				iMap2[iPos2] = iPos1;
				}
			else
				{
				iMap1[iPos1] = -1;
				iMap2[iPos2] = -1;
				}
			++iPos1;
			++iPos2;
			}
		else if (!bIsGap1 && bIsGap2)
			{
			iMap1[iPos1] = -1;
			++iPos1;
			}
		else if (bIsGap1 && !bIsGap2)
			{
			iMap2[iPos2] = -1;
			++iPos2;
			}
		}

#if	_DEBUG
	{
	int iLength1 = iPos1;
	int iLength2 = iPos2;

	for (int iPos1 = 0; iPos1 < iLength1; ++iPos1)
		{
		int iPos2 = iMap1[iPos1];
		if (-1 == iPos2)
			continue;
		assert(iMap2[iPos2] == iPos1);
		}

	for (int iPos2 = 0; iPos2 < iLength2; ++iPos2)
		{
		int iPos1 = iMap2[iPos2];
		if (-1 == iPos1)
			continue;
		assert(iMap1[iPos1] == iPos2);
		}
	}
#endif
	}
/////////////

/////////////
/// Copied from (prefab-4's) qscore's comparepair.cpp
void ComparePair (const MSA &msaTest, unsigned uTestSeqIndexA,
  unsigned uTestSeqIndexB, const MSA &msaRef, unsigned uRefSeqIndexA,
  unsigned uRefSeqIndexB, double *ptrdSP, double *ptrdPS, double *ptrdCS)
	{
	const int iLengthA = (int) msaTest.GetSeqLength(uTestSeqIndexA);
	const int iLengthB = (int) msaTest.GetSeqLength(uTestSeqIndexB);
	const int iLengthAr = (int) msaRef.GetSeqLength(uRefSeqIndexA);
	const int iLengthBr = (int) msaRef.GetSeqLength(uRefSeqIndexB);

	if (iLengthA != iLengthAr || iLengthB != iLengthBr)
		Quit("Lengths differ");

	int *iRefMapA = new int[iLengthA];
	int *iRefMapB = new int[iLengthB];
	int *iTestMapA = new int[iLengthA];
	int *iTestMapB = new int [iLengthB];

	msaTest.GetPairMap(uTestSeqIndexA, uTestSeqIndexB, iTestMapA, iTestMapB);
	msaRef.GetPairMap(uRefSeqIndexA, uRefSeqIndexB, iRefMapA, iRefMapB);

	ComparePairMap(iTestMapA, iTestMapB, iRefMapA, iRefMapB, iLengthA, iLengthB,
	  ptrdSP, ptrdPS, ptrdCS);

	delete[] iRefMapA;
	delete[] iRefMapB;
	delete[] iTestMapA;
	delete[] iTestMapB;
	}
/////////////

/////////////
/// Copied from (prefab-4's) qscore's comparemsa.cpp
unsigned MSA::GetUngappedColIndex (unsigned uSeqIndex, unsigned uColIndex)
	{
	const unsigned uSeqCount = GetSeqCount();
	const unsigned uColCount = GetColCount();

	assert(uSeqIndex < uSeqCount);
	assert(uColIndex < uColCount);

	if (0 == m_UngapMap)
		{
		m_UngapMap = new unsigned *[uSeqCount];
		memset(m_UngapMap, 0, uSeqCount*sizeof(unsigned *));
		}

	unsigned *ptrMap = m_UngapMap[uSeqIndex];
	if (0 == ptrMap)
		{
		ptrMap = new unsigned[uColCount];
		memset(ptrMap, 0, uColCount*sizeof(unsigned));
		unsigned uUngappedColIndex = 0;
		for (unsigned uGappedColIndex = 0; uGappedColIndex < uColCount;
		  ++uGappedColIndex)
			if (IsGap(uSeqIndex, uGappedColIndex))
				ptrMap[uGappedColIndex] = uInsane;
			else
				{
				ptrMap[uGappedColIndex] = uUngappedColIndex;
				++uUngappedColIndex;
				}
		m_UngapMap[uSeqIndex] = ptrMap;
		}
	unsigned uUngappedColIndex = ptrMap[uColIndex];
	if (uInsane == uUngappedColIndex)
		Quit("GetUngappedColIndex(%u,%u)", uSeqIndex, uUngappedColIndex);
	return uUngappedColIndex;
	}

unsigned MSA::GetGappedColIndex (unsigned uSeqIndex, unsigned uUngappedColIndex)
	{
	unsigned n = 0;
	for (unsigned uColIndex = 0; uColIndex < GetColCount(); ++uColIndex)
		{
		if (!IsGap(uSeqIndex, uColIndex))
			{
			if (n == uUngappedColIndex)
				{
				assert(GetUngappedColIndex(uSeqIndex, uColIndex) == uUngappedColIndex);
				return uColIndex;
				}
			++n;
			}
		}
	assert(false);
	return 0;
	}

void CompareMSA (const MSA &msaTest, const MSA &msaRef, double *ptrdSP,
  double *ptrdPS, double *ptrdCS)
	{
	const unsigned uRefSeqCount = msaRef.GetSeqCount();

	double dTotalSP = 0.0;
	double dTotalPS = 0.0;
	double dTotalCS = 0.0;
	unsigned uPairCount = 0;

	for (unsigned uRefSeqIndexA = 0; uRefSeqIndexA < uRefSeqCount; ++uRefSeqIndexA)
		{
		const char *pstrSeqNameA = msaRef.GetSeqName(uRefSeqIndexA);
		unsigned uTestSeqIndexA;
		bool bFound = msaTest.GetSeqIndex(pstrSeqNameA, &uTestSeqIndexA);
		if (!bFound)
			{
			Quit("Sequence '%s' not found in test alignment", pstrSeqNameA);
			continue;
			}

		for (unsigned uRefSeqIndexB = uRefSeqIndexA + 1; uRefSeqIndexB < uRefSeqCount;
		  ++uRefSeqIndexB)
			{
			const char *pstrSeqNameB = msaRef.GetSeqName(uRefSeqIndexB);
			unsigned uTestSeqIndexB;
			bool bFound = msaTest.GetSeqIndex(pstrSeqNameB, &uTestSeqIndexB);
			if (!bFound)
				{
				Quit("Sequence '%s' not found in test alignment", pstrSeqNameA);
				continue;
				}

			double dSP = dInsane;
			double dPS = dInsane;
			double dCS = dInsane;
			ComparePair(msaTest, uTestSeqIndexA, uTestSeqIndexB, msaRef, uRefSeqIndexA,
			  uRefSeqIndexB, &dSP, &dPS, &dCS);

			dTotalSP += dSP;
			dTotalPS += dPS;
			dTotalCS += dCS;
			++uPairCount;
			}
		}
	if (0 == uPairCount)
		{
		Quit("No sequence pairs in common between test and reference alignment");
		*ptrdSP = 0;
		*ptrdPS = 0;
		*ptrdCS = 0;
		return;
		}

	*ptrdSP = dTotalSP / uPairCount;
	*ptrdPS = dTotalPS / uPairCount;
	*ptrdCS = dTotalCS / uPairCount;
	}
/////////////

/////////////
/// Copied from (prefab-4's) qscore's tc.cpp
static void SeqsDiffer (MSA &msaTest, unsigned uTestSeqIndex1, MSA &msaRef, 
  unsigned uRefSeqIndex1, unsigned uRefColIndex)
	{
	Quit("Test & ref sequences differ, label=%s", msaTest.GetSeqName(uTestSeqIndex1));
	}

static int IsAlignedCol (const MSA &msa, unsigned uColIndex)
	{
	const unsigned uSeqCount = msa.GetSeqCount();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		if (msa.IsGap(uSeqIndex, uColIndex))
			continue;
		return isupper(msa.GetChar(uSeqIndex, uColIndex));
		}
	return false;
	}

double TC (MSA &msaTest, MSA &msaRef)
	{
	const unsigned uRefSeqCount = msaRef.GetSeqCount();
	if (0 == uRefSeqCount)
		Quit("No sequences in ref alignment");

	unsigned *RefSeqIndexToTestSeqIndex = new unsigned[uRefSeqCount];
	for (unsigned uRefSeqIndex = 0; uRefSeqIndex < uRefSeqCount; ++uRefSeqIndex)
		{
		const char *ptrName = msaRef.GetSeqName(uRefSeqIndex);
		unsigned uTestSeqIndex;
		bool bFound = msaTest.GetSeqIndex(ptrName, &uTestSeqIndex);
		if (bFound)
			RefSeqIndexToTestSeqIndex[uRefSeqIndex] = uTestSeqIndex;
		else
			RefSeqIndexToTestSeqIndex[uRefSeqIndex] = uInsane;
		}

	unsigned uRefAlignedColCount = 0;
	unsigned uCorrectlyAlignedColCount = 0;
	unsigned uRefColCount = msaRef.GetColCount();
	for (unsigned uRefColIndex = 0; uRefColIndex < uRefColCount; ++uRefColIndex)
		{
		if (!IsAlignedCol(msaRef, uRefColIndex))
			continue;

		bool bAllAlignedCorrectly = true;
		bool bAllGaps = true;
	// Iterate over all pairs
		for (unsigned uRefSeqIndex1 = 0; uRefSeqIndex1 < uRefSeqCount;
		  ++uRefSeqIndex1)
			{
			unsigned uTestSeqIndex1 = RefSeqIndexToTestSeqIndex[uRefSeqIndex1];
			if (uTestSeqIndex1 == uInsane)
				continue;

			char cRef1 = msaRef.GetChar(uRefSeqIndex1, uRefColIndex);
			if (IsGap(cRef1))
				continue;
			if (!isupper(cRef1))
				Quit("Ref alignment col %d has both upper and lower-case letters", uRefColIndex);
			unsigned uRefUngappedColIndex1 =
			  msaRef.GetUngappedColIndex(uRefSeqIndex1, uRefColIndex);
			unsigned uTestGappedColIndex1 = msaTest.
			  GetGappedColIndex(uTestSeqIndex1, uRefUngappedColIndex1);
			char cTest1 = msaTest.GetChar(uTestSeqIndex1, uTestGappedColIndex1);
			if (cRef1 != toupper(cTest1))
				SeqsDiffer(msaTest, uTestSeqIndex1, msaRef, uRefSeqIndex1, uRefColIndex);
			for (unsigned uRefSeqIndex2 = uRefSeqIndex1 + 1; uRefSeqIndex2 < uRefSeqCount;
			  ++uRefSeqIndex2)
				{
				unsigned uTestSeqIndex2 = RefSeqIndexToTestSeqIndex[uRefSeqIndex2];
				if (uTestSeqIndex2 == uInsane)
					continue;

				char cRef2 = msaRef.GetChar(uRefSeqIndex2, uRefColIndex);
				if (IsGap(cRef2))
					continue;
				bAllGaps = false;
				assert(isupper(cRef2));

				unsigned uRefUngappedColIndex2 =
				  msaRef.GetUngappedColIndex(uRefSeqIndex2, uRefColIndex);
				unsigned uTestGappedColIndex2 = msaTest.
				  GetGappedColIndex(uTestSeqIndex2, uRefUngappedColIndex2);

				char cTest2 = msaTest.GetChar(uTestSeqIndex2, uTestGappedColIndex2);
				if (!(isupper(cTest1) && isupper(cTest2) &&
				  uTestGappedColIndex1 == uTestGappedColIndex2))
					{
					bAllAlignedCorrectly = false;
					goto NextCol;
					}
				}
			}

	NextCol:
		if (!bAllGaps)
			{
			++uRefAlignedColCount;
			if (bAllAlignedCorrectly)
				++uCorrectlyAlignedColCount;
			}
		}
	delete[] RefSeqIndexToTestSeqIndex;

	if (0 == uRefAlignedColCount)
		Quit("No aligned columns (upper case) in ref alignment");

	return (double) uCorrectlyAlignedColCount / (double) uRefAlignedColCount;
	}
/////////////

/////////////
// Copied from hmmer's hmmbuild.c and modified
namespace hmmer {
extern "C" {
/* Function: position_average_score()
 * Date:     Wed Dec 31 09:36:35 1997 [StL]
 * 
 * Purpose:  Calculate scores from tracebacks, keeping them
 *           in a position specific array. The final array
 *           is normalized position-specifically too, according
 *           to how many sequences contributed data to this
 *           position. Used for compensating for sequence 
 *           fragments in ME and MD score optimization. 
 *           Very much ad hoc.
 *           
 *           Code related to (derived from) TraceScore().
 *           
 * Args:     hmm       - HMM structure, scores valid
 *           dsq       - digitized unaligned sequences
 *           wgt       - weights on the sequences
 *           nseq      - number of sequences
 *           tr        - array of nseq tracebacks that aligns each dsq to hmm
 *           pernode   - RETURN: [0]1..M array of position-specific avg scores
 *           ret_avg   - RETURN: overall average full-length, one-domain score
 *           
 * Return:   1 on success, 0 on failure.          
 *           pernode is malloc'ed [0]1..M by CALLER and filled here.
 */
static void
position_average_score (struct plan7_s    *hmm, 
		       unsigned char    **dsq, 
		       float             *wgt,
		       int                nseq,
		       struct p7trace_s **tr,
		       float             *pernode,
		       float             *ret_avg)
{
  unsigned char sym;
  int    pos;                   /* position in seq */
  int    tpos;                  /* position in trace/state sequence */
  float *counts;                /* counts at each position */
  float  avg;                   /* RETURN: average overall */
  int    k;                     /* counter for model position */
  int    idx;                   /* counter for sequence number */

  /* Allocations
   */
  counts = ( float * )MallocOrDie ((hmm->M+1) * sizeof(float));
  FSet(pernode, hmm->M+1, 0.);
  FSet(counts,  hmm->M+1, 0.);

  /* Loop over traces, accumulate weighted scores per position
   */
  for (idx = 0; idx < nseq; idx++)
    for (tpos = 0; tpos < tr[idx]->tlen; tpos++)
      {
	pos = tr[idx]->pos[tpos];
	sym = dsq[idx][tr[idx]->pos[tpos]];
	k   = tr[idx]->nodeidx[tpos];

	/* Counts: how many times did we use this model position 1..M?
         * (weighted)
	 */
	if (tr[idx]->statetype[tpos] == STM || tr[idx]->statetype[tpos] == STD)
	  counts[k] += wgt[idx];

	/* Emission scores.
	 */
	if (tr[idx]->statetype[tpos] == STM) 
	  pernode[k] += wgt[idx] * Scorify(hmm->msc[sym][k]);
	else if (tr[idx]->statetype[tpos] == STI) 
	  pernode[k] += wgt[idx] * Scorify(hmm->isc[sym][k]);
	
	/* Transition scores.
	 */
	if (tr[idx]->statetype[tpos] == STM ||
	    tr[idx]->statetype[tpos] == STD ||
	    tr[idx]->statetype[tpos] == STI)
	  pernode[k] += wgt[idx] * 
	    Scorify(TransitionScoreLookup(hmm, tr[idx]->statetype[tpos], tr[idx]->nodeidx[tpos],
					  tr[idx]->statetype[tpos+1],tr[idx]->nodeidx[tpos+1]));
      }

  /* Divide accumulated scores by accumulated weighted counts
   */
  avg = 0.;
  for (k = 1; k <= hmm->M; k++)
    {
      pernode[k] /= counts[k];
      avg += pernode[k];
    }
  
  free(counts);
  *ret_avg = avg;
  return;
}

/* Function: frag_trace_score()
 * Date:     SRE, Wed Dec 31 10:03:47 1997 [StL]
 * 
 * Purpose:  Allow MD/ME optimization to be used for alignments
 *           that include fragments and multihits -- estimate a full-length
 *           per-domain score.
 *           
 *
 *           
 * Return:   "corrected" score.
 */
static float
frag_trace_score (struct plan7_s *hmm, unsigned char *dsq, struct p7trace_s *tr, 
                 float *pernode, float expected)
{
  float sc;			/* corrected score  */
  float fragexp;		/* expected score for a trace like this */
  int   tpos;			/* position in trace */

                                /* get uncorrected score */
  sc = P7TraceScore(hmm, dsq, tr);

                               /* calc expected score for trace like this */
  fragexp = 0.;
  for (tpos = 0; tpos < tr->tlen; tpos++)
    if (tr->statetype[tpos] == STM || tr->statetype[tpos] == STD)
      fragexp += pernode[tr->nodeidx[tpos]];

				/* correct for multihits */
  fragexp /= (float) TraceDomainNumber(tr);

                                /* extrapolate to full-length, one-hit score */
  sc = sc * expected / fragexp;
  return sc;
}

/* Function: maximum_entropy()
 * Date:     SRE, Fri Jan  2 10:56:00 1998 [StL]
 * 
 * Purpose:  Optimizes a model according to maximum entropy weighting.
 *           See Krogh and Mitchison (1995).
 *
 *           [Actually, we do minimum relative entropy, rather than
 *           maximum entropy. Same thing, though we refer to "ME"
 *           weights and models. The optimization is a steepest
 *           descents minimization of the relative entropy.]
 *           
 *           Expects to be called shortly after a Maxmodelmaker()
 *           or Handmodelmaker(), so that both a new model architecture
 *           (with MAP parameters) and fake tracebacks are available.
 *           
 *           Prints a summary of optimization progress to stdout.
 *           
 * Args:     hmm     - model. allocated, set with initial MAP parameters.
 *           dsq     - dealigned digitized seqs the model is based on
 *           ainfo   - extra info for aseqs
 *           nseq    - number of aseqs
 *           eff_nseq- effective sequence number; weights normalize up to this.
 *           prior   - prior distributions for parameterizing model
 *           tr      - array of fake traces for each sequence        
 *           
 * Return:   (void)
 *           hmm changed to an ME HMM
 *           ainfo changed, contains ME weights          
 */
static void
maximum_entropy (struct plan7_s *hmm, unsigned char **dsq, MSA *msa,
		float eff_nseq, struct p7prior_s *prior, struct p7trace_s **tr)
{
  float *wgt;                  /* current best set of ME weights   */
  float *new_wgt;              /* new set of ME weights to try     */
  float *sc;                    /* log-odds score of each sequence */
  float *grad;                  /* gradient                        */
  float  epsilon;               /* steepness of descent            */
  float  relative_entropy;      /* current best relative entropy   */
  float  new_entropy;           /* relative entropy at new weights */
  float  last_new_entropy;      /* last new_entropy we calc'ed     */
  float  use_epsilon;           /* current epsilon value in use    */
  int    idx;                   /* counter over sequences          */
  int    i1, i2;		/* counters for iterations         */

  float  converge_criterion;
  float  minw, maxw;            /* min, max weight                 */
  int    posw, highw;           /* number of positive weights      */
  float  mins, maxs, avgs;      /* min, max, avg score             */
  float *pernode;               /* expected score per node of HMM  */
  float  expscore;              /* expected score of complete HMM  */
  int    max_iter;		/* bulletproof against infinite loop bugs */

  epsilon  = 0.2;                /* works fine */
  max_iter = 666;
  
  /* Allocations
   */
  sc      = ( float * )MallocOrDie (sizeof(float) * msa->nseq);
  wgt     = ( float * )MallocOrDie (sizeof(float) * msa->nseq);
  new_wgt = ( float * )MallocOrDie (sizeof(float) * msa->nseq);
  grad    = ( float * )MallocOrDie (sizeof(float) * msa->nseq);
  pernode = ( float * )MallocOrDie (sizeof(float) * (hmm->M+1));

  /* Initialization. Start with all weights == 1.0.
   * Find relative entropy and gradient.
   */
  Plan7SWConfig(hmm, 0.5, 0.5);
  P7Logoddsify(hmm, TRUE);

  FSet(wgt, msa->nseq, 1.0);
  position_average_score(hmm, dsq, wgt, msa->nseq, tr, pernode,&expscore);
  for (idx = 0; idx < msa->nseq; idx++) 
    sc[idx] = frag_trace_score(hmm, dsq[idx], tr[idx], pernode, expscore);
  relative_entropy = FSum(sc, msa->nseq) / (float) msa->nseq;
  for (idx = 0; idx < msa->nseq; idx++)
    grad[idx] = relative_entropy - sc[idx];

  
  printf("iter avg-sc min-sc max-sc min-wgt max-wgt +wgt ++wgt rel.ent convergence\n");
  printf("---- ------ ------ ------ ------- ------- ---- ----- ------- -----------\n");
  mins = maxs = avgs = sc[0];
  for (idx = 1; idx < msa->nseq; idx++)
    {
      if (sc[idx] < mins) mins = sc[idx];
      if (sc[idx] > maxs) maxs = sc[idx];
      avgs += sc[idx];
    }
  avgs /= (float) msa->nseq;
  printf("%4d %6.1f %6.1f %6.1f %7.2f %7.2f %4d %5d %7.2f %8s\n",
         0, avgs, mins, maxs, 1.0, 1.0, msa->nseq, 0, relative_entropy, "-");

  
  /* Steepest descents optimization;
   * iterate until relative entropy converges.
   */
  i1 = 0;
  while (++i1 < max_iter)
    {
      /* Gradient gives us a line of steepest descents.
       * (Roughly speaking, anyway. We actually have a constraint
       * that weights are nonnegative and normalized, and the
       * gradient doesn't take these into account.)
       * Look along this line, a distance of epsilon * gradient:
       * if new point is better, accept; if new point is worse,
       * move back along the line by half the distance and re-evaluate.
       */
      use_epsilon = epsilon;
      new_entropy = relative_entropy + 1.0;    /* just ensure new > old */

      i2 = 0; 
      while (new_entropy > relative_entropy && ++i2 < max_iter)
        {
          last_new_entropy = new_entropy;

                                /* find a new point in weight space */
          for (idx = 0; idx < msa->nseq; idx++)
            {
              new_wgt[idx] = wgt[idx] + use_epsilon * grad[idx];
              if (new_wgt[idx] < 0.) new_wgt[idx] = 0.0;
            }
          FNorm(new_wgt, msa->nseq);
          FScale(new_wgt, msa->nseq, (float) msa->nseq);

                                /* Make new HMM using these weights */
          ZeroPlan7(hmm);
          for (idx = 0; idx < msa->nseq; idx++)
            P7TraceCount(hmm, dsq[idx], new_wgt[idx], tr[idx]);
          P7PriorifyHMM(hmm, prior);

  
                                /* Evaluate new point */
	  Plan7SWConfig(hmm, 0.5, 0.5);
	  P7Logoddsify(hmm, TRUE);
	  position_average_score(hmm, dsq, new_wgt, msa->nseq, tr, pernode, &expscore);
          for (idx = 0; idx < msa->nseq; idx++) 
	    sc[idx]      = frag_trace_score(hmm, dsq[idx], tr[idx], pernode, expscore);
	  new_entropy = FDot(sc, new_wgt, msa->nseq) / (float) msa->nseq;

          use_epsilon /= 2.0;
	  /* Failsafe: we're not converging. Set epsilon to zero,
	   * do one more round.
	   */
	  if (use_epsilon < 1e-6) use_epsilon = 0.0; 
	  if (use_epsilon == 0.0) break;
          
          /* Failsafe: avoid infinite loops. Sometimes the
             new entropy converges without ever being better 
             than the previous point, probably as a result
             of minor roundoff error. */
          if (last_new_entropy == new_entropy) break;
        }
      if (i2 == max_iter) printf("   -- exceeded maximum iterations; giving up --\n");

      /* Evaluate convergence before accepting the new weights;
       * then, accept the new point and evaluate the gradient there.
       */
      converge_criterion = fabs((relative_entropy-new_entropy)/relative_entropy);
      relative_entropy = new_entropy;
      FCopy(wgt, new_wgt, msa->nseq);
      for (idx = 0; idx < msa->nseq; idx++)
	grad[idx] = relative_entropy - sc[idx];

      /* Print some statistics about this iteration
       */
      mins = maxs = avgs = sc[0];
      minw = maxw = wgt[0];
      posw = (wgt[0] > 0.0) ? 1 : 0;
      highw = (wgt[0] > 1.0) ? 1 : 0;
      for (idx = 1; idx < msa->nseq; idx++)
        {
          if (sc[idx] < mins) mins = sc[idx];
          if (sc[idx] > maxs) maxs = sc[idx];
          if (wgt[idx] < minw) minw = wgt[idx];
          if (wgt[idx] > maxw) maxw = wgt[idx];
          if (wgt[idx] > 0.0)  posw++;
          if (wgt[idx] > 1.0)  highw++;
          avgs += sc[idx];
        }
      avgs /= (float) msa->nseq;
      printf("%4d %6.1f %6.1f %6.1f %7.2f %7.2f %4d %5d %7.2f %8.5f\n",
             i1, 
             avgs, mins, maxs, 
             minw, maxw, posw, highw,
             relative_entropy, converge_criterion);
      
      if (converge_criterion < 1e-5) break;
    }
  if (i1 == max_iter) printf("   -- exceeded maximum iterations; giving up --\n");

  /* Renormalize weights to sum to eff_nseq, and save.
   */
  FNorm(wgt, msa->nseq);
  FScale(wgt, msa->nseq, (float) eff_nseq);
  FCopy(msa->wgt, wgt, msa->nseq);
			/* Make final HMM using these adjusted weights */
  ZeroPlan7(hmm);
  for (int idx = 0; idx < msa->nseq; idx++)
    P7TraceCount(hmm, dsq[idx], wgt[idx], tr[idx]);
  P7PriorifyHMM(hmm, prior);
                                
  /* Cleanup and return
   */
  free(pernode);
  free(new_wgt);
  free(grad);
  free(wgt);
  free(sc);
  return;
}

enum p7_weight {		/* weighting strategy */
  WGT_NONE, WGT_GSC, WGT_BLOSUM, WGT_PB, WGT_VORONOI, WGT_ME};
/*
 * Run hmmbuild as if you were running it from the command line (with the
 * default parameters).  Allocates hmm; to free, do FreePlan7( hmm ).
 */
struct plan7_s *
hmmbuild ( MSA const * const_msa, char const * name = "" )
{
  struct plan7_s  *hmm;         /* constructed HMM  */

  // We don't actually modify the msa, but the hmmer code doesn't use the const keyword.
  MSA * msa = const_cast<MSA *>( const_msa );

  // Params to hmmbuild procedure
  char *rndfile;		/* random sequence model file to read    */
  char *prifile;		/* Dirichlet prior file to read          */
  char *pamfile;		/* PAM matrix file for heuristic prior   */
  int   do_eff;			/* TRUE to set an effective seq number   */
  float blosumlevel;		/* BLOSUM frac id filtering level [0.62] */
  enum p7_weight w_strategy;	/* weighting strategy */
  int   pbswitch;		/* nseq >= this, switchover to PB weights*/
  enum p7_construction c_strategy; /* construction strategy choice        */
  float gapmax;			/* max frac gaps in mat col for -k       */
  float archpri;		/* "architecture" prior on model size    */
  float swentry;		/* S/W aggregate entry probability       */
  float swexit;			/* S/W aggregate exit probability        */
  enum p7_config {              /* algorithm configuration strategy      */
    P7_BASE_CONFIG, P7_LS_CONFIG, P7_FS_CONFIG, P7_SW_CONFIG } cfg_strategy;

  // Parameter defaults
  rndfile           = NULL;
  prifile           = NULL;
  pamfile           = NULL;
  do_eff            = TRUE;
  blosumlevel       = 0.62;
  w_strategy        = WGT_GSC;
  pbswitch          = 1000;
  c_strategy        = P7_MAP_CONSTRUCTION;
  gapmax            = 0.5;
  archpri           = 0.85;
  swentry           = 0.5;
  swexit            = 0.5;
  cfg_strategy      = P7_LS_CONFIG;

  /* Make alignment upper case, because some symbol counting
   * things are case-sensitive.
   */
  for(int idx = 0; idx < msa->nseq; idx++) {
    s2upper( msa->aseq[ idx ] );
  }
  // Set the alphabet type
  SetAlphabet( hmmAMINO );

  // Local vars for hmmbuild procedure
  float pamwgt;			/* weight on PAM for heuristic prior     */
  struct p7prior_s *pri; /* Dirichlet priors to use                 */
  float  randomseq[MAXABET];	/* null sequence model                     */
  float            p1;		/* null sequence model p1 transition       */
  unsigned char  **dsq;         /* digitized unaligned aseq's              */ 
  float eff_nseq;		/* effective sequence number             */
  int              checksum;	/* checksum of the alignment               */
  struct p7trace_s **tr; /* fake tracebacks for aseq's              */ 

  /* Set up Dirichlet priors */
  if (prifile == NULL)  pri = P7DefaultPrior();
  else                  pri = P7ReadPrior(prifile);
  if (pamfile != NULL)  PAMPrior(pamfile, pri, pamwgt);
  /* Set up the null/random seq model */
  if (rndfile == NULL)  P7DefaultNullModel(randomseq, &p1);
  else                  P7ReadNullModel(rndfile, randomseq, &p1);

  /* Prepare unaligned digitized sequences for internal use 
   */
  DigitizeAlignment(msa, &dsq);
  
  /* In some respects we treat DNA more crudely for now;
   * for example, we can't do eff seq #, because it's
   * calibrated for protein.
   */
  if( Alphabet_type == hmmNUCLEIC ) 
    do_eff = FALSE;

  /* Determine "effective sequence number".
   * The BlosumWeights() routine is now an efficient O(N)
   * memory clustering algorithm that doesn't blow up on,
   * say, Pfam's GP120 alignment (13000+ sequences)
   */
  eff_nseq = (float) msa->nseq;
  if( do_eff ) {
    float *wgt;
    printf("%-40s ... ", "Determining effective sequence number");
    fflush(stdout);
    /* dummy weights array to feed BlosumWeights*/
    wgt = ( float * )MallocOrDie(sizeof(float) * msa->nseq);
    BlosumWeights(msa->aseq, msa->nseq, msa->alen, blosumlevel, wgt);
    eff_nseq = FSum(wgt, msa->nseq);

    free(wgt);
    printf("done. [%.0f]\n", eff_nseq);
  } // End if do_eff

  /* Weight the sequences (optional),
   */
  if( w_strategy == WGT_GSC     || 
      w_strategy == WGT_BLOSUM  || 
      w_strategy == WGT_VORONOI ||
      w_strategy == WGT_PB ) {
    printf("%-40s ... ", "Weighting sequences heuristically");
    fflush(stdout);

    if( w_strategy != WGT_PB && msa->nseq >= pbswitch ) {
      printf("[big alignment! doing PB]... ");
      PositionBasedWeights(msa->aseq, msa->nseq, msa->alen, msa->wgt);
    } else if( w_strategy == WGT_GSC ) {
      GSCWeights(msa->aseq, msa->nseq, msa->alen, msa->wgt);
    } else if( w_strategy == WGT_BLOSUM ) {
      BlosumWeights(msa->aseq, msa->nseq, msa->alen, blosumlevel, msa->wgt);
    } else if( w_strategy == WGT_PB ) {
      PositionBasedWeights(msa->aseq, msa->nseq, msa->alen, msa->wgt);
    } else if( w_strategy ==  WGT_VORONOI ) {
      VoronoiWeights(msa->aseq, msa->nseq, msa->alen, msa->wgt);
    }
    printf("done.\n");
  } // End if w_strategy == WGT_GSC || WGT_BLOSUM  || WGT_VORONOI || WGT_PB

  /* Set the effective sequence number (if do_eff is FALSE, eff_nseq 
   * was set to nseq).
   */
  FNorm(msa->wgt,  msa->nseq);
  FScale(msa->wgt, msa->nseq, eff_nseq);

  /* Build a model architecture.
   * If we're not doing MD or ME, that's all we need to do.
   * We get an allocated, counts-based HMM back.
   * 
   * Because the architecture algorithms are allowed to change
   * gap characters in the alignment, we have to calculate the
   * alignment checksum before we enter the algorithms.
   */
  printf("%-40s ... ", "Constructing model architecture");
  fflush(stdout);
  checksum = GCGMultchecksum(msa->aseq, msa->nseq);
  if( c_strategy == P7_FAST_CONSTRUCTION ) {
    P7Fastmodelmaker(msa, dsq, gapmax, &hmm, &tr);
  } else if( c_strategy == P7_HAND_CONSTRUCTION ) {
    P7Handmodelmaker(msa, dsq, &hmm, &tr);
  } else {
    P7Maxmodelmaker(msa, dsq, gapmax, pri, randomseq, p1, archpri, &hmm, &tr);
  }
  hmm->checksum = checksum;
  printf("done.\n");

  /* Record the null model in the HMM;
   * add prior contributions in pseudocounts and renormalize.
   */
  printf("%-40s ... ", "Converting counts to probabilities");
  fflush(stdout);
  Plan7SetNullModel(hmm, randomseq, p1);
  P7PriorifyHMM(hmm, pri);
  printf("done.\n");

  /* Model configuration, temporary.
   * hmmbuild assumes that it's given an alignment of single domains,
   * and the alignment may contain fragments. So, for the purpose of
   * scoring the sequences (or, optionally, MD/ME weighting),
   * configure the model into hmmsw mode. Later we'll
   * configure the model according to how the user wants to
   * use it.
   */
  Plan7SWConfig(hmm, 0.5, 0.5);

  /* Do model-dependent "weighting" strategies.
   */
  if( w_strategy == WGT_ME ) {
    printf("\n%-40s ...\n", "Maximum entropy weighting, iterative");
    maximum_entropy(hmm, dsq, msa, eff_nseq, pri, tr);
    printf("----------------------------------------------\n\n");
  } // End if w_strategy == WGT_ME

  /* Give the model a name.
   * We deal with this differently depending on whether
   * we're in an alignment database or a single alignment.
   * 
   * If a single alignment, priority is:
   *      1. Use -n <name> if set.
   *      2. Use msa->name (avail in Stockholm or SELEX formats only)
   *      3. If all else fails, use alignment file name without
   *         filename extension (e.g. "globins.slx" gets named "globins"
   *         
   * If a multiple MSA database (e.g. Stockholm/Pfam), 
   * only msa->name is applied. -n is not allowed.
   * if msa->name is unavailable, or -n was used,
   * a fatal error is thrown.
   * 
   * Because we can't tell whether we've got more than one
   * alignment 'til we're on the second one, these fatal errors
   * only happen after the first HMM has already been built.
   * Oh well.
   */
  printf("%-40s ... ", "Setting model name, etc.");
  fflush(stdout);
  Plan7SetName(hmm, const_cast<char *>( name ) ); // const_cast should be safe.

  /* Transfer other information from the alignment to
   * the HMM. This typically only works for Stockholm or SELEX format
   * alignments, so these things are conditional/optional.
   */
  if( msa->acc  != NULL) Plan7SetAccession(hmm,   msa->acc);
  if( msa->desc != NULL) Plan7SetDescription(hmm, msa->desc);

  if( msa->cutoff_is_set[MSA_CUTOFF_GA1] && msa->cutoff_is_set[MSA_CUTOFF_GA2] ) {
    hmm->flags |= PLAN7_GA;
    hmm->ga1 = msa->cutoff[MSA_CUTOFF_GA1];
    hmm->ga2 = msa->cutoff[MSA_CUTOFF_GA2];
  }
  if( msa->cutoff_is_set[MSA_CUTOFF_TC1] && msa->cutoff_is_set[MSA_CUTOFF_TC2] ) {
    hmm->flags |= PLAN7_TC;
    hmm->tc1 = msa->cutoff[MSA_CUTOFF_TC1];
    hmm->tc2 = msa->cutoff[MSA_CUTOFF_TC2];
  }
  if( msa->cutoff_is_set[MSA_CUTOFF_NC1] && msa->cutoff_is_set[MSA_CUTOFF_NC2] ) {
    hmm->flags |= PLAN7_NC;
    hmm->nc1 = msa->cutoff[MSA_CUTOFF_NC1];
    hmm->nc2 = msa->cutoff[MSA_CUTOFF_NC2];
  }

  /* Record some other miscellaneous information in the HMM,
   * like how/when we built it.
   */
  //Plan7ComlogAppend(hmm, argc, argv);
  Plan7SetCtime(hmm);
  hmm->nseq = msa->nseq;
  printf("done. [%s]\n", hmm->name); 
   
  /* Print information for the user
   */
  printf("\nConstructed a profile HMM (length %d)\n", hmm->M);
  PrintPlan7Stats(stdout, hmm, dsq, msa->nseq, tr); 
  printf("\n");

  /* Configure the model for chosen algorithm
   */
  printf("%-40s ... ", "Finalizing model configuration");
  fflush(stdout);
  switch( cfg_strategy ) {
  case P7_BASE_CONFIG:  Plan7GlobalConfig(hmm);              break;
  case P7_SW_CONFIG:    Plan7SWConfig(hmm, swentry, swexit); break;
  case P7_LS_CONFIG:    Plan7LSConfig(hmm);                  break;
  case P7_FS_CONFIG:    Plan7FSConfig(hmm, swentry, swexit); break;
  default:              Die("bogus configuration choice");
  } // End switch( cfg_strategy )
  printf("done.\n");

  P7FreePrior(pri);
  SqdClean();

  return hmm;
} // hmmbuild( MSA const *, char const * )
} // End extern "C"
} // End namespace hmmer

// In C++:
namespace hmmer {
  /*
   * Allocate and return a new hmmer::MSA with the data from the given Muscle
   * MSA object.  Caller must free the returned MSA when done,
   * using MSAFree( MSA * ).
   */
  hmmer::MSA *
  createMSAFromMuscleMSA ( ::MSA const & muscle_msa )
  {
    int muscle_msa_seq_count = muscle_msa.GetSeqCount();
    int muscle_msa_col_count = muscle_msa.GetColCount();
    
    hmmer::MSA * msa =
      hmmer::MSAAlloc( muscle_msa_seq_count, muscle_msa_col_count );
    // TODO: REMOVE
    cout << "Allocated hmmer MSA for " << muscle_msa_seq_count << " gapped sequences, of length " << muscle_msa_col_count << " each." << endl;
    int idx;
    for( uint32_t seq_i = 0; seq_i < muscle_msa_seq_count; seq_i++ ) {
      // TODO: REMOVE
      //cout << "Attempting to convert gapped sequence " << seq_i << " to the hmmer MSA format." << endl;
      idx = // const_cast should be safe here.  hmmer code should declare as const.
        hmmer::GKIStoreKey( msa->index, const_cast<char *>( muscle_msa.GetSeqName( seq_i ) ) );
      if( idx != seq_i ) {
        // TODO: ?
        cout << "ERROR: Expected hmmer::GKIStoreKey( msa->index, name ) to return " << seq_i << " but it instead returned " << idx << endl;
      }
    
      msa->sqname[ seq_i ] = // const_cast should be safe here.
        hmmer::sre_strdup( const_cast<char *>( muscle_msa.GetSeqName( seq_i ) ), -1 );
      msa->aseq[ seq_i ] = // const_cast should be safe here.
        hmmer::sre_strdup( const_cast<char *>( muscle_msa.GetSeqBuffer( seq_i ) ), muscle_msa_col_count );
      msa->sqlen[ seq_i ] = muscle_msa_col_count;
    } // End foreach seq_i
    // Set nseq.
    msa->nseq = muscle_msa_seq_count;
    // TODO: REMOVE
    //cout << "Done converting MSA.  Verifying." << endl;
    hmmer::MSAVerifyParse( msa );
    cout << "Muscle MSA successfully converted to HMMer MSA." << endl;
    cout << "Average unaligned sequence length is " << hmmer::MSAAverageSequenceLength( msa ) << endl;

    return msa;
  } // createMSAFromMuscleMSA( std::MSA const & )
} // End namespace hmmer
/////////////

/////////////
// Copied from muscle's main.cpp and modified
int g_argc;
char **g_argv;
// Just like calling muscle from the command line.
// As a side effect, sets many useful global vars, including:
// ALPHA g_Alpha; // ALPHA_Amino, _RNA, or _DNA
////// WARNING: Right now I'm not using this because it results in an unstable state of the program.  Presumably muscle is doing something it should not be doing.  I have looked for the problem and found others, but nothing that fixes this.  Probably there is a pointer to a temporary variable that is causing the problem.  Anyway, muscle was not designed to be used after Run() finishes, it seems.  So I'm switching to a fork..exec solution.
int
run_muscle ( int argc, char **argv )
{
  #if	WIN32
  // Multi-tasking does not work well in CPU-bound
  // console apps running under Win32.
  // Reducing the process priority allows GUI apps
  // to run responsively in parallel.
    SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
  #endif
    g_argc = argc;
    g_argv = argv;
  
    SetNewHandler();
    SetStartTime();
    ProcessArgVect(argc - 1, argv + 1);
    SetParams();
    SetLogFile();

  //extern void TestSubFams(const char *);
  //TestSubFams(g_pstrInFileName);
  //return 0;

  if ( g_bVersion ) {
    printf(MUSCLE_LONG_VERSION "\n");
    exit(EXIT_SUCCESS);
  }
  
  if( !g_bQuiet ) {
    Credits();
  }
  
  if( MissingCommand() && isatty( 0 ) ) {
    // TODO: Make our own usage(), etc.?
    Usage();
    exit( EXIT_SUCCESS );
  }
  
  
  if( g_bCatchExceptions ) {
    try {
      Run();
    }
    catch (...) {
      OnException();
      exit( EXIT_Except );
    }
  } else {
    Run();
  }
  
  //exit( EXIT_Success );

  return 0;
} // run_muscle( int, char ** )
/////////////

namespace galosh {

/////////////
// Copy rom the given (probs-based) hmmer hmm into the given galosh Profile.
// If the hmm is counts-based, caller should first call hmmer::Plan7Renormalize( hmm )
template <class ProfileType>
void
copyFromHMMerProfile (
  ProfileType & profile,
  struct hmmer::plan7_s const * hmm
)
{
  assert( hmm != 0 );

  typedef typename ProfileType::ProfileResidueType ResidueType;

  profile.reinitialize( static_cast<uint32_t>( hmm->M ) );
    
  cout << "Converting hmmer hmm '" << hmm->name << "' to galosh Profile..";
  cout.flush();
    
  profile.zeroExceptPositions();
    
  profile[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] =
    ( hmm->xt[XTN][LOOP] );
  profile[ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ] =
    ( hmm->xt[XTN][MOVE] );
  profile[ Transition::fromPreAlign ].normalize( 1E-5 );
  
  profile[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ] =
    ( hmm->tbd1 );
  profile[ Transition::fromBegin ][ TransitionFromBegin::toMatch ] =
    ( 1 ) -
    profile[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ];
  
  // TODO: Use null model?  Or overall Insertion dist'n?
  profile[ Emission::PreAlignInsertion ].even();
  
  for( uint32_t pos_i = 0; pos_i < profile.length(); pos_i++ ) {
    cout << '.';
    cout.flush();

    profile[ pos_i ][ Emission::Match ].zero();
    for( uint32_t res_i = 0; res_i < seqan::ValueSize<ResidueType>::VALUE; res_i++ ) {
      // TODO: REMOVE
      //if( pos_i == 309 ) {
      //  ResidueType foo = ResidueType( res_i );
      //  cout << "Residue: " << foo << endl;
      //  char bar = (char)foo;
      //  cout << "Residue as char: " << bar << endl;
      //  unsigned char baz = hmmer::SymbolIndex( bar );
      //  cout << "SymbolIndex( Residue as char ): " << (int)baz << endl;
      //  cout << "profile[ Emission::Insertion ]: " <<
      //    profile[ Emission::Insertion] << endl;
      //  cout << "profile[ Emission::Insertion ][ ResidueType( res_i ) ]: " <<
      //    profile[ Emission::Insertion][ foo ] << endl;
      //  //cout << "profile[ Emission::Insertion ][ Residue as char ]: " <<
      //  //  profile[ Emission::Insertion][ bar ] << endl;
      //  cout << "profile[ Emission::Insertion ][ " << res_i << " ]: " <<
      //    profile[ Emission::Insertion][ res_i ] << endl;
      //  cout << "hmm->ins[ pos_i + 1 ][ hmmer::SymbolIndex( ( char )ResidueType( res_i ) ) ]: " <<
      //    hmm->ins[ pos_i + 1 ][ baz ];
      //}
      profile[ pos_i ][ Emission::Match ][ res_i ] =
        ( hmm->mat[ pos_i + 1 ][ hmmer::SymbolIndex( ( char )ResidueType( res_i ) ) ] );
      // NOTE: In hmmer models there's no insert distribution for the last position.  TODO: Do they use the null distribution for pre-align and post-align insertions?  Do they train that distribution?  Should we glean anything from it, eg include it in our sum here?
      if( pos_i < ( profile.length() - 1 ) ) {
        // Since these are global, we add them.  TODO: Weigh by 1-P(del)
        //profile[ pos_i ][ Emission::Insertion ][ res_i ] =
        profile[ Emission::Insertion ][ res_i ] +=
          //( hmm->ins[ pos_i + 1 ][ hmmer::SymbolIndex( ( char )ResidueType( res_i ) ) ] );
          hmm->ins[ pos_i + 1 ][ hmmer::SymbolIndex( ( char )ResidueType( res_i ) ) ];
      }
    } // End foreach res_i
    profile[ pos_i ][ Emission::Match ].normalize( 1E-5 ); 
    //profile[ pos_i ][ Emission::Insertion ].normalize( 1E-5 ); 
  
    if( pos_i != ( hmm->M - 1 ) ) { // Last pos has no transitions
      //profile[ pos_i ][ Transition::fromMatch ][ TransitionFromMatch::toMatch ] =
      profile[ Transition::fromMatch ][ TransitionFromMatch::toMatch ] +=
        //( hmm->t[ pos_i + 1 ][ 0 ] );
        hmm->t[ pos_i + 1 ][ 0 ];
      //profile[ pos_i ][ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] =
      profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] +=
        //( hmm->t[ pos_i + 1 ][ 1 ] );
        hmm->t[ pos_i + 1 ][ 1 ];
      //profile[ pos_i ][ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] =
      profile[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] +=
        //( hmm->t[ pos_i + 1 ][ 2 ] );
        hmm->t[ pos_i + 1 ][ 2 ];
      //profile[ pos_i ][ Transition::fromMatch ].normalize( 1E-5 );
      //profile[ pos_i ][ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ] =
      profile[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ] +=
        //( hmm->t[ pos_i + 1 ][ 3 ] );
        hmm->t[ pos_i + 1 ][ 3 ];
      //profile[ pos_i ][ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] =
      profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] +=
        //( hmm->t[ pos_i + 1 ][ 4 ] );
        hmm->t[ pos_i + 1 ][ 4 ];
      profile[ pos_i ][ Transition::fromInsertion ].normalize( 1E-5 );
      // TODO: REMOVE
      //cout << "pos_i is " << pos_i << "; D->M is " << hmm->t[ pos_i + 1 ][ 5 ] << "; D->D is " << hmm->t[ pos_i + 1 ][ 6 ] << endl;
      //profile[ pos_i ][ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ] =
      profile[ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ] +=
        //( hmm->t[ pos_i + 1 ][ 5 ] );
        hmm->t[ pos_i + 1 ][ 5 ];
      //profile[ pos_i ][ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ] =
      profile[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ] +=
        //( hmm->t[ pos_i + 1 ][ 6 ] );
        hmm->t[ pos_i + 1 ][ 6 ];
      profile[ pos_i ][ Transition::fromDeletion ].normalize( 1E-5 );
    } // End if this is not the last position  
  } // End foreach pos_i
  
  // Since thse are global, we sum and normalize afterwards
  profile[ Emission::Insertion ].normalize( 1E-5 ); 
  profile[ Transition::fromMatch ].normalize( 1E-5 );
  profile[ Transition::fromInsertion ].normalize( 1E-5 );
  profile[ Transition::fromDeletion ].normalize( 1E-5 );
  
  // TODO: Read this from the hmm? (probly we shouldn't for now)
  //profile[ Transition::fromEnd ][ TransitionFromEnd::toPostAlign ] = ( 1.0 );
  //profile[ Transition::fromEnd ][ TransitionFromEnd::toLoop ] = ( 0.0 );
  
  profile[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] =
    ( hmm->xt[XTC][LOOP] );
  
  profile[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] =
    ( hmm->xt[XTC][MOVE] );
  profile[ Transition::fromPostAlign ].normalize();
  
  // TODO: Use null model?  Or overall Insertion dist'n?
  profile[ Emission::PostAlignInsertion ].even();

  // Use the insertion dist'n for the pre- and post- align insertions too.
  profile[ Emission::PreAlignInsertion ] =
    profile[ Emission::Insertion ];
  profile[ Emission::PostAlignInsertion ] =
    profile[ Emission::Insertion ];
  
  cout << ".done." << endl;
} // copyFromHMMerProfile( ProfileType &, hmmer::hmm const * )

int
prolific ( int argc, char **argv )
{
  // 0th arg is progname
  // 1st arg is infile name
  // 2nd arg is run name
  if( argc != 3 ) {
    cout << "\tUsage: prolific fasta_file run_name" << endl;
    exit( 1 );
  }

#ifdef __PROFUSE_USE_AMINOS
  typedef seqan::AminoAcid20 ResidueType;
  typedef seqan::AminoAcid SequenceResidueType;
#else // __PROFUSE_USE_AMINOS .. else
  typedef seqan::Dna ResidueType;
  typedef seqan::Iupac SequenceResidueType;
#endif // __PROFUSE_USE_AMINOS .. else ..

  //typedef bfloat ProbabilityType;
  //typedef logspace ProbabilityType;
  //typedef floatrealspace ProbabilityType; // For speed (vs doublerealspace)
  typedef doublerealspace ProbabilityType;
  
  typedef bfloat ScoreType; // Preferred
  //typedef logspace ScoreType; // SLOWer than bfloat
  //typedef realspace ScoreType; // Only for very few & small sequences
  
  typedef bfloat MatrixValueType;
  //typedef logspace MatrixValueType;
  //typedef doublerealspace MatrixValueType; // For speed (vs bfloat)
  //typedef floatrealspace MatrixValueType;

  typedef ProfileTreeRoot<ResidueType, ProbabilityType> ProfileType;
  typedef ProfileTreeRoot<ResidueType, ProbabilityType> InternalNodeType;

  string fasta_in_str = argv[ 1 ];
  string run_str = argv[ 2 ];

  string ref_dir_str = "prefab4-ref/";
  static const uint32_t max_fasta_in_len = 256;
  char fasta_in_buf[ max_fasta_in_len ];
  NameFromPath( fasta_in_str.c_str(), fasta_in_buf, max_fasta_in_len );
  string fasta_str = fasta_in_buf;

  string fasta_ref_str = ref_dir_str + fasta_str;
  // TODO: autodetect existence of ref file
  bool have_ref = true; //false;

  string filename_base =
    ( fasta_in_str + ".prolific." + run_str );

  cout << "Run name is \"" << run_str << "\"." << endl;
  cout << "Input fasta filename is \"" << fasta_in_str << "\"." << endl;
  //cout << "Base fasta filename is \"" << fasta_str << "\"." << endl;
  cout << "Reference fasta filename is \"" << fasta_ref_str << "\"." << endl;
  cout << "Output filenames begin with \"" << filename_base << "\"." << endl;

  string filename;

  // TODO: PUT BACK
  //run_muscle( argc, argv ); // Doesn't work.  Leaves state unstable.  Argh.

  filename = fasta_in_str;
  filename += ".muscle";
  const char * const argv_muscle[] =
    { "muscle", "-maxiters", "2", "-stable", "-noanchors", "-in", fasta_in_str.c_str(), "-out", filename.c_str(), 0 };
  const int argc_muscle = 9;
  cout << "Running muscle/muscle";
  for( int arg_i = 1; arg_i < argc_muscle; arg_i++ ) {
    cout << ' ' << argv_muscle[ arg_i ];
  }
  cout << endl;
  pid_t muscle_pid = spawn( "muscle/muscle", argv_muscle );
  int muscle_exit_status;
  wait( &muscle_exit_status );
  cout << "Ran muscle.  Exit status was " << muscle_exit_status << endl;

  // TODO: REMOVE
  //g_pstrInFileName = "prefab4-in/1atiA_1qf6A";
  //g_pstrOutFileName = "prefab4-in/1atiA_1qf6A.muscle";

  // Set up muscle stuff
  SetStartTime();
  ProcessArgVect(argc_muscle - 1, const_cast<char **>( argv_muscle ) + 1 );
  SetParams();
  SetLogFile();
  SetOutputFileName(g_pstrOutFileName);
  SetInputFileName(g_pstrInFileName);
  SetMaxIters(g_uMaxIters);
  SetSeqWeightMethod(g_SeqWeight2);
  SetAlpha( ALPHA_Amino );
  SetPPScore( true );

  bool do_asis = false;
  bool do_split = true;
  bool do_muscle_after_split = false;//true;

  bool use_baldi_siegel = true;
  bool do_Cbw_Ubw = false;
  bool do_Ubw = false;//true;

  bool use_priors = false;
  bool train_globals = true;//false;
  bool use_lengthadjust = true;
  double lengthadjust_threshold = .5;//.15;//.5;//.025;
  double lengthadjust_threshold_increment = .0005;

  bool save_scorefiles = true;//false;


  ProbabilityType very_small( 1E-5 ); // profileValueMinimum

  // First score the original muscled msa
  // TODO: Why do I have to read it in?  Why can't I just use ptrBestMSA?  Crashes.
  MSA msa_muscle;
  filename = g_pstrOutFileName; // fasta_in_str + muscle_suffix_str;
  TextFile msa_muscle_infile =
    TextFile( filename.c_str() );
  msa_muscle.FromFile( msa_muscle_infile );
  msa_muscle_infile.Close();
  cout << "READ in MSA from " << msa_muscle_infile.GetFileName() << endl;

  // Set up
  const unsigned uSeqCount = msa_muscle.GetSeqCount();
  msa_muscle.FixAlpha();
  MSA::SetIdCount(uSeqCount);
  for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
    msa_muscle.SetSeqId(uSeqIndex, uSeqIndex);

  // Compute guide tree (necessary to set muscle weights)
  cout << "Calculating guide tree." << endl;
  Tree tree;
  //TreeFromMSA( msa_muscle, tree, g_Cluster2, g_Distance2, g_Root2, g_pstrDistMxFileName2 );
  //TreeFromMSA( msa_muscle, tree, g_Cluster2, g_Distance2, ROOT_MinAvgLeafDist, g_pstrDistMxFileName2 );
  TreeFromMSA( msa_muscle, tree, CLUSTER_NeighborJoining, g_Distance2, ROOT_MinAvgLeafDist, g_pstrDistMxFileName2 );
  //TreeFromMSA( msa_muscle, tree, CLUSTER_UPGMA, g_Distance2, ROOT_MinAvgLeafDist, g_pstrDistMxFileName2 );
  cout << "Guide tree calculated." << endl;
  SetMuscleTree( tree );

  // Make sure to apply the sequence weights.
  SetMSAWeightsMuscle( msa_muscle );

  cout << "Total SP score is " << getNonInsertionSPScore( msa_muscle ) << endl;
  if( save_scorefiles ) {
    filename += ".spscores";
    g_pstrScoreFileName = filename.c_str();
    WriteScoreFile( msa_muscle );
    cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
    g_pstrScoreFileName = NULL;
  } // End if save_scorefiles

  // Now read the ref seq and qscore the orginial muscle one
  MSA msa_ref;
  double dSP, dPS, dCS, dTC;
  if( have_ref ) {
    filename = fasta_ref_str;
    TextFile msa_ref_infile =
      TextFile( filename.c_str() );
    msa_ref.FromFile( msa_ref_infile );
    msa_ref_infile.Close();
    cout << "READ in MSA from " << msa_ref_infile.GetFileName() << endl;
    dSP = dPS = dCS = dInsane;
    CompareMSA( msa_muscle, msa_ref, &dSP, &dPS, &dCS);
    dTC = TC( msa_muscle, msa_ref );
    printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
           dSP, dPS, dCS, dTC);
  } // End if have_ref

  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType> dp;
  ofstream outfile;

  ScoreType score;
  if( do_asis ) {
    // hmmer
    // Need to convert from a muscle MSA to a hmmer MSA.
    hmmer::MSA * hmmer_msa_muscle =
      hmmer::createMSAFromMuscleMSA( msa_muscle );
    
    // Build HMM using hmmer's hmmbuild
    struct hmmer::plan7_s *hmm = hmmer::hmmbuild( hmmer_msa_muscle );
    if (hmm == NULL) hmmer::Die( "Null HMM" );          // NULL on file parse failure
    // Convert counts to probs
    if( !( hmm->flags | PLAN7_HASPROB ) ) {
      cout << "Calling Plan7Renormalize( hmm )." << endl;
      hmmer::Plan7Renormalize( hmm );
    }
    
    cout << "Done building hmm." << endl;
    
    //exit( 0 );
    
    ProfileType profile;
    copyFromHMMerProfile( profile, hmm );
    //cout << "The profile is:" << endl;
    //cout << profile << endl;
    
    hmmer::FreePlan7(hmm);

    // Load the fasta seqs
    Fasta<SequenceResidueType> fasta_in;
    fasta_in.fromFile( fasta_in_str );
    
    uint32_t num_sequences_to_use = fasta_in.size();
    //cout << "FASTA from file:" << endl;
    //if( num_sequences_to_use ) {
    //  fasta_in.writeFasta( cout, num_sequences_to_use );
    //  cout << endl;
    //} else {
    //  cout << fasta_in << endl;
    //}

    // using the profile we read from hmmer
    // First fix the preAlign and postAlign, which for some reason are excessively high when we read from the hmmer profile.
    //profile[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] =
    //  .1;
    //profile[ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ] =
    //  1 -
    //  profile[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ];
    //profile[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] =
    //  .1;
    //profile[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] =
    //  1 -
    //  profile[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ];
    
    // TODO: REMOVE?  Testing making the indel probs very low before running viterbi.
    profile[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ] =
      very_small;
    profile[ Transition::fromBegin ][ TransitionFromBegin::toMatch ] =
    ( 1 ) -
      profile[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ];
    profile[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] =
      ( .5 ); // very_small;
    profile[ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ] =
      ( 1 ) -
      profile[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ];
    profile[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] =
      ( .5 ); // very_small;
    profile[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] =
      ( 1 ) -
      profile[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ];
    profile[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] =
      very_small;
    profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] =
      very_small;
    profile[ Transition::fromMatch ][ TransitionFromMatch::toMatch ] =
      ( 1 ) -
      profile[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] -
      profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ];  

    ProfileType viterbi_profile = profile;
    ProfileType notraining_profile = profile;

    ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType> trainer_Cbw( &profile, fasta_in, num_sequences_to_use );

    //cout << "The profile (before) is:" << endl;
    //cout << *trainer_Cbw.m_profile << endl;
    
    // TODO: REMOVE:
    //cout << "fasta_in.size() is " << fasta_in.size() << endl;
    //cout << "m_sequence_count is " << trainer_Cbw.m_sequence_count << endl;

    // TODO: Try with the lengthadjust: maxBaumWelchInverseScalar > 0
    trainer_Cbw.m_parameters.minIterations = 1;
    trainer_Cbw.m_parameters.maxIterations = 1000;
    // Having maxPositionCycles_globals > 1 seems ok; takes about the same
    // number of iterations, converges to roughly the same place; takes
    // longer by virtue of having more pos cycles per iteration of course.
    trainer_Cbw.m_parameters.maxPositionCycles = 1;
    // Having maxPositionCycles_globals > 1 seems to make convergence way
    // slower when lengthadjust is on.  Length keeps adjusting..
    trainer_Cbw.m_parameters.maxPositionCycles_globals = 1;
    trainer_Cbw.m_parameters.minBaumWelchInverseScalar = 0; // Straight-up bw.
    trainer_Cbw.m_parameters.maxBaumWelchInverseScalar = 0; // Straight-up bw.
    trainer_Cbw.m_parameters.minBaumWelchInverseScalar_globals = 0; // Straight-up bw.
    trainer_Cbw.m_parameters.maxBaumWelchInverseScalar_globals = 0; // Straight-up bw.
    trainer_Cbw.m_parameters.scorePercentChangeMinimum_position_cycle = 1;
    trainer_Cbw.m_parameters.scorePercentChangeMinimum_iteration = .01;
    // TODO: Debug alwaysAccept..?  It doesn't *always*.. ok?
    trainer_Cbw.m_parameters.alwaysAccept = false;//true;
  
    trainer_Cbw.m_parameters.proposeProfileLengthChanges = use_lengthadjust;
    trainer_Cbw.m_parameters.useAlignmentProfiles = true;
    trainer_Cbw.m_parameters.numIterationsBetweenLengthChanges = 0;
    trainer_Cbw.m_parameters.proposeDeletingThreshold =
      lengthadjust_threshold; //.01;//.025;//.1;
    trainer_Cbw.m_parameters.proposeDeletingThreshold_increment =
      lengthadjust_threshold_increment; //.0005;//.00005;//.0005;//5E-5;//.0003125;//.00625;//.025;
    trainer_Cbw.m_parameters.proposeInsertingThreshold =
      trainer_Cbw.m_parameters.proposeDeletingThreshold;// / 4;//seqan::ValueSize<ResidueType>::VALUE; // TODO: Figure this out...
    trainer_Cbw.m_parameters.proposeInsertingPreAlignThreshold = //.35; //.5;
      trainer_Cbw.m_parameters.proposeInsertingThreshold;
    trainer_Cbw.m_parameters.proposeInsertingPostAlignThreshold = //.35;//.5;
      trainer_Cbw.m_parameters.proposeInsertingThreshold;
    trainer_Cbw.m_parameters.proposeInsertingThreshold_increment =
      trainer_Cbw.m_parameters.proposeDeletingThreshold_increment;// TODO: ERE I AM! PUT BACK? / seqan::ValueSize<ResidueType>::VALUE;
  
    //if( have_trained_profile && start_with_trained_profile ) {
    //  // When we start with the trained profile, we need to get past the
    //  // length modification wait time (which is just one iteration).
    //  trainer_Cbw.m_parameters.minIterations =
    //    max( trainer_Cbw.m_parameters.minIterations, ( uint32_t )2 );
    //}
  
    // Use rabiner scaling? (default true)
    // NOTE: You must change the MatrixValueType to logspace or bfloat iff this is false!\
    // TODO: I don't think I've tested this in a long long time.  Probably safest to disable it for now.  At least make the default false.
    trainer_Cbw.m_parameters.useRabinerScaling = false;

    // Train globals first?
    //trainer_Cbw.m_parameters.trainGlobalsFirst = true; // note: breaks UBW.

    // Use Ubw?
    trainer_Cbw.m_parameters.useUnconditionalBaumWelch = false;//true;
    trainer_Cbw.m_parameters.unconditionalIsolatesGlobals = false;

    trainer_Cbw.m_parameters.trainProfileGlobals = train_globals;
    //trainer_Cbw.m_parameters.maxPositionCycles = 3;
    //trainer_Cbw.m_parameters.useUnconditionalBaumWelch = true;
    trainer_Cbw.m_parameters.usePriors = use_priors;
    //trainer_Cbw.m_parameters.debug = DEBUG_All;
    //trainer_Cbw.m_parameters.verbosity = VERBOSITY_All;
    trainer_Cbw.m_parameters.verbosity = VERBOSITY_Low;
    
    // Use Baldi?
    // For testing Baldi-style gradient ascent
#ifdef ALLOW_BOLTZMANN_GIBBS
    if( use_baldi_siegel ) {
      // NOTE about priors:  since globals are presently not updated using Baldi, you can still usePriors and they will affect the globals *but not the positions*.
      trainer_Cbw.m_parameters.baldiLearningRate = 1; // 0 means noBaldi!
      trainer_Cbw.m_parameters.baldiTemperature = 1;
      trainer_Cbw.m_parameters.baldiHybrid = false;
      trainer_Cbw.m_parameters.siegelMaxFindingThePeakAttempts_positions = 1000; // 0 means Baldi not Baldi / Siegel !!!
      trainer_Cbw.m_parameters.siegelEpsilonScaleFactor = 1.5;
      trainer_Cbw.m_parameters.siegelMaxRefiningThePeakSteps_positions = 1;//1000;
      trainer_Cbw.m_parameters.siegelRefiningThePeakStepsConvergenceThreshold = 1E-5;
      trainer_Cbw.m_parameters.minBaumWelchInverseScalar = 0;
      trainer_Cbw.m_parameters.maxBaumWelchInverseScalar = 0;
      //trainer_Cbw.m_parameters.maxPositionCycles = 10;
    } // use_baldi_siegel
#endif //ALLOW_BOLTZMANN_GIBBS

    DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer viterbi_matrices(
                                                                                                                                                profile,
        fasta_in,
        num_sequences_to_use
      );

    // TODO: REMOVE
    cout << "Calculating viterbi alignments and scores for the 'notraining' profile." << endl;
    viterbi_profile.copyPositions( *trainer_Cbw.m_profile );
    score =
      dp.forward_score_viterbi(
          trainer_Cbw.m_parameters,
          viterbi_profile,
          fasta_in,
          num_sequences_to_use,
          viterbi_matrices
        );
    cout << "The total score for all sequences, using viterbi, is: " << score << endl;
    DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::MultipleAlignment<ProfileType, SequenceResidueType> ma_notraining(
      &viterbi_profile, // trainer_Cbw.m_profile,
      &fasta_in,
      num_sequences_to_use
    );
    dp.forward_viterbiAlign(
      trainer_Cbw.m_parameters,
      viterbi_matrices,
      ma_notraining
    );
    cout << "Done." << endl;

    filename = filename_base + ".notraining.blast";
    outfile.open( filename.c_str() );
    //cout << "Multiple Alignment is:" << endl;
    //ma_notraining.toPairwiseStream( cout, &fasta_in.m_descriptions );
    ma_notraining.toPairwiseStream( outfile, &fasta_in.m_descriptions );
    //ma_notraining.toPileupStream( outfile, &fasta_in.m_descriptions );
    cout << "Wrote out Multiple Alignment to " << filename << endl;
    outfile.close();

    // Now try making an MSA out of it.
    ::MSA msa_notraining;
    ma_notraining.appendToMSA( &msa_notraining, &fasta_in.m_descriptions );
    filename = filename_base + ".notraining.blast.a2m";
    TextFile msa_notraining_outfile =
      TextFile( filename.c_str(), true );
    msa_notraining.ToFile( msa_notraining_outfile );
    msa_notraining_outfile.Close();
    cout << "Wrote MSA to file " << filename << endl;
    
    // Trying score-getting
    cout << "Total SP score is " << getNonInsertionSPScore( msa_notraining ) << endl;
    if( save_scorefiles ) {
      filename += ".spscores";
      g_pstrScoreFileName = filename.c_str();
      WriteScoreFile( msa_notraining );
      cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
      g_pstrScoreFileName = NULL;
    } // End if save_scorefiles
    // Qscore
    if( have_ref ) {
      dSP = dPS = dCS = dInsane;
      CompareMSA( msa_notraining, msa_ref, &dSP, &dPS, &dCS);
      dTC = TC( msa_notraining, msa_ref );
      printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
             dSP, dPS, dCS, dTC);
    } // End if have_ref
    
    // TODO: REMOVE!
    //exit( 0 );
    
    cout << "Training, using conditional Baum-Welch until convergence." << endl;//until convergence." << endl;
    score = trainer_Cbw.train();
    cout << "Now (after training), the score is " << score << endl; //", and the profile is:" << endl;
    //cout << *trainer_Cbw.m_profile << endl;
    //ViterbiDPType::Matrix::SequentialAccessContainer viterbi_matrices =
    //  ViterbiDPType::Matrix::SequentialAccessContainer(
    //    *trainer_Cbw.m_profile,
    //    trainer_Cbw.m_poly_sequences,
    //    trainer_Cbw.m_sequence_count
    //  );
    // TODO: REMOVE
    cout << "Calculating viterbi alignments and scores for the 'Cbw' profile." << endl;
    viterbi_profile.copyPositions( *trainer_Cbw.m_profile );
    score =
      dp.forward_score_viterbi(
          trainer_Cbw.m_parameters,
          viterbi_profile,
          fasta_in,
          num_sequences_to_use,
          viterbi_matrices
        );
    cout << "The total score for all sequences, using viterbi, is: " << score << endl;
    DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::MultipleAlignment<ProfileType, SequenceResidueType> ma_Cbw(
      &viterbi_profile, // *trainer_Cbw.m_profile,
      &fasta_in,
      num_sequences_to_use
    );
    dp.forward_viterbiAlign(
      trainer_Cbw.m_parameters,
      viterbi_matrices,
      ma_Cbw
    );
    cout << "The total score for all sequences, using viterbi, is: " << ma_Cbw.calculateScore() << endl;

    filename = filename_base + ".Cbw.blast";
    outfile.open( filename.c_str() );
    //cout << "Multiple Alignment is:" << endl;
    //ma_Cbw.toPairwiseStream( cout, &fasta_in.m_descriptions );
    ma_Cbw.toPairwiseStream( outfile, &fasta_in.m_descriptions );
    //ma_Cbw.toPileupStream( outfile, &fasta_in.m_descriptions );
    cout << "Wrote out Multiple Alignment to " << filename << endl;
    outfile.close();

    // Now try making an MSA out of it.
    ::MSA msa_Cbw;
    ma_Cbw.appendToMSA( &msa_Cbw, &fasta_in.m_descriptions );
    filename = filename_base + ".Cbw.blast.a2m";
    TextFile msa_Cbw_outfile =
      TextFile( filename.c_str(), true );
    msa_Cbw.ToFile( msa_Cbw_outfile );
    msa_Cbw_outfile.Close();
    cout << "Wrote MSA to file " << filename << endl;
    
    // Trying score-getting
    cout << "Total SP score is " << getNonInsertionSPScore( msa_Cbw ) << endl;
    if( save_scorefiles ) {
      filename += ".spscores";
      g_pstrScoreFileName = filename.c_str();
      WriteScoreFile( msa_Cbw );
      cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
      g_pstrScoreFileName = NULL;
    } // End if save_scorefiles
    // Qscore
    if( have_ref ) {
      dSP = dPS = dCS = dInsane;
      CompareMSA( msa_Cbw, msa_ref, &dSP, &dPS, &dCS);
      dTC = TC( msa_Cbw, msa_ref );
      printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
             dSP, dPS, dCS, dTC);
    } // End if have_ref
    
//    if( do_Cbw_Ubw ) {
//      // Now train using Ubw
//      cout << "Training again using Unconditional Baum-Welch for 2 more iterations" << endl;
//      trainer_Cbw.m_parameters.maxIterations = 2;
//      trainer_Cbw.m_parameters.scorePercentChangeMinimum_position_cycle =
//        DEFAULT_scorePercentChangeMinimum_iteration;
//      trainer_Cbw.m_parameters.scorePercentChangeMinimum_iteration =
//        DEFAULT_scorePercentChangeMinimum_position_cycle;
//      trainer_Cbw.m_parameters.useUnconditionalBaumWelch = true;
//      score = trainer_Cbw.train();
//      cout << "Now (after training again using Ubw for 2 iterations), the score is " << score << endl;//", and the profile is:" << endl;
//      //cout << *trainer_Cbw.m_profile << endl;
//      // TODO: REMOVE
//      cout << "Calculating viterbi scores for the 'Cbw_Ubw2' profile." << endl;
//      viterbi_profile.copyPositions( *trainer_Cbw.m_profile );
//      score =
//        viterbi_dp.forward_score_viterbi(
//          trainer_Cbw.m_parameters,
//          viterbi_profile, // *trainer_Cbw.m_profile,
//          trainer_Cbw.m_poly_sequences,
//          trainer_Cbw.m_sequence_count,
//          viterbi_matrices
//        );
//      cout << "The total score for all sequences, using viterbi, is: " << score << endl;
//      ViterbiDPType::MultipleAlignment<ProfileType> ma_Cbw_Ubw2 =
//        viterbi_dp.forward_viterbiAlign(
//          trainer_Cbw.m_parameters,
//          viterbi_profile, // *trainer_Cbw.m_profile,
//          trainer_Cbw.m_poly_sequences,
//          trainer_Cbw.m_sequence_count,
//          viterbi_matrices // TODO: ARE THESE MATRICES SIZED CORRECTLY?
//        );
//      filename = filename_base + ".Cbw_Ubw2.blast";
//      outfile.open( filename.c_str() );
//      //cout << "Multiple Alignment is:" << endl;
//      //ma_Cbw_Ubw2.toPairwiseStream( cout, &fasta_in.m_descriptions );
//      ma_Cbw_Ubw2.toPairwiseStream( outfile, &fasta_in.m_descriptions );
//      cout << "Wrote out Multiple Alignment to " << filename << endl;
//      outfile.close();
//      
//      // Now try making an MSA out of it.
//      ::MSA msa_Cbw_Ubw2;
//      ma_Cbw_Ubw2.appendToMSA( &msa_Cbw_Ubw2, &fasta_in.m_descriptions );
//      filename = filename_base + ".Cbw_Ubw2.blast.a2m";
//      TextFile msa_Cbw_Ubw2_outfile =
//        TextFile( filename.c_str(), true );
//      msa_Cbw_Ubw2.ToFile( msa_Cbw_Ubw2_outfile );
//      msa_Cbw_Ubw2_outfile.Close();
//      cout << "Wrote MSA to file " << filename << endl;
//      
//      // Trying score-getting
//      cout << "Total SP score is " << getNonInsertionSPScore( msa_Cbw_Ubw2 ) << endl;
//      if( save_scorefiles ) {
//        filename += ".spscores";
//        g_pstrScoreFileName = filename.c_str();
//        WriteScoreFile( msa_Cbw_Ubw2 );
//        cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
//        g_pstrScoreFileName = NULL;
//      } // End if save_scorefiles
//      // Qscore
//      if( have_ref ) {
//        dSP = dPS = dCS = dInsane;
//        CompareMSA( msa_Cbw_Ubw2, msa_ref, &dSP, &dPS, &dCS);
//        dTC = TC( msa_Cbw_Ubw2, msa_ref );
//        printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
//               dSP, dPS, dCS, dTC);
//      } // End if have_ref
//      
//      cout << "Training yet again using Unconditional Baum-Welch for 2 more iterations" << endl;
//      trainer_Cbw.m_parameters.maxIterations = 2;
//      trainer_Cbw.m_parameters.useUnconditionalBaumWelch = true;
//      score = trainer_Cbw.train();
//      cout << "Now (after training again using Ubw for a total of 4 iterations max), the score is " << score << endl;//", and the profile is:" << endl;
//      //cout << *trainer_Cbw.m_profile << endl;
//      // TODO: REMOVE
//      cout << "Calculating viterbi scores for the 'Cbw_Ubw4' profile." << endl;
//      viterbi_profile.copyPositions( *trainer_Cbw.m_profile );
//      score =
//        viterbi_dp.forward_score_viterbi(
//          trainer_Cbw.m_parameters,
//          viterbi_profile, // *trainer_Cbw.m_profile,
//          trainer_Cbw.m_poly_sequences,
//          trainer_Cbw.m_sequence_count,
//          viterbi_matrices
//        );
//      cout << "The total score for all sequences, using viterbi, is: " << score << endl;
//      ViterbiDPType::MultipleAlignment<ProfileType> ma_Cbw_Ubw4 =
//        viterbi_dp.forward_viterbiAlign(
//          trainer_Cbw.m_parameters,
//          viterbi_profile, // *trainer_Cbw.m_profile,
//          trainer_Cbw.m_poly_sequences,
//          trainer_Cbw.m_sequence_count,
//          viterbi_matrices // TODO: ARE THESE MATRICES SIZED CORRECTLY?
//        );
//      filename = filename_base + ".Cbw_Ubw4.blast";
//      outfile.open( filename.c_str() );
//      //cout << "Multiple Alignment is:" << endl;
//      //ma_Cbw_Ubw4.toPairwiseStream( cout, &fasta_in.m_descriptions );
//      ma_Cbw_Ubw4.toPairwiseStream( outfile, &fasta_in.m_descriptions );
//      cout << "Wrote out Multiple Alignment to " << filename << endl;
//      outfile.close();
//      
//      // Now try making an MSA out of it.
//      ::MSA msa_Cbw_Ubw4;
//      ma_Cbw_Ubw4.appendToMSA( &msa_Cbw_Ubw4, &fasta_in.m_descriptions );
//      filename = filename_base + ".Cbw_Ubw4.blast.a2m";
//      TextFile msa_Cbw_Ubw4_outfile =
//        TextFile( filename.c_str(), true );
//      msa_Cbw_Ubw4.ToFile( msa_Cbw_Ubw4_outfile );
//      msa_Cbw_Ubw4_outfile.Close();
//      cout << "Wrote MSA to file " << filename << endl;
//      
//      // Trying score-getting
//      cout << "Total SP score is " << getNonInsertionSPScore( msa_Cbw_Ubw4 ) << endl;
//      if( save_scorefiles ) {
//        filename += ".spscores";
//        g_pstrScoreFileName = filename.c_str();
//        WriteScoreFile( msa_Cbw_Ubw4 );
//        cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
//        g_pstrScoreFileName = NULL;
//      } // End if save_scorefiles
//      // Qscore
//      if( have_ref ) {
//        dSP = dPS = dCS = dInsane;
//        CompareMSA( msa_Cbw_Ubw4, msa_ref, &dSP, &dPS, &dCS );
//        dTC = TC( msa_Cbw_Ubw4, msa_ref );
//        printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
//               dSP, dPS, dCS, dTC);
//      } // End if have_ref
//
//      cout << "Training yet again using Unconditional Baum-Welch until convergence." << endl;
//      trainer_Cbw.m_parameters.maxIterations = 100;
//      trainer_Cbw.m_parameters.scorePercentChangeMinimum_position_cycle = 100;
//      trainer_Cbw.m_parameters.scorePercentChangeMinimum_iteration = 100;
//      trainer_Cbw.m_parameters.useUnconditionalBaumWelch = true;
//      score = trainer_Cbw.train();
//      cout << "Now (after training again using Ubw until convergence), the score is " << score << endl;//", and the profile is:" << endl;
//      //cout << *trainer_Cbw.m_profile << endl;
//      // TODO: REMOVE
//      cout << "Calculating viterbi scores for the 'Cbw_Ubw' profile." << endl;
//      viterbi_profile.copyPositions( *trainer_Cbw.m_profile );
//      score =
//        viterbi_dp.forward_score_viterbi(
//          trainer_Cbw.m_parameters,
//          viterbi_profile, // *trainer_Cbw.m_profile,
//          trainer_Cbw.m_poly_sequences,
//          trainer_Cbw.m_sequence_count,
//          viterbi_matrices
//        );
//      cout << "The total score for all sequences, using viterbi, is: " << score << endl;
//      ViterbiDPType::MultipleAlignment<ProfileType> ma_Cbw_Ubw =
//        viterbi_dp.forward_viterbiAlign(
//          trainer_Cbw.m_parameters,
//          viterbi_profile, // *trainer_Cbw.m_profile,
//          trainer_Cbw.m_poly_sequences,
//          trainer_Cbw.m_sequence_count,
//          viterbi_matrices // TODO: ARE THESE MATRICES SIZED CORRECTLY?
//        );
//      filename = filename_base + ".Cbw_Ubw.blast";
//      outfile.open( filename.c_str() );
//      //cout << "Multiple Alignment is:" << endl;
//      //ma_Cbw_Ubw.toPairwiseStream( cout, &fasta_in.m_descriptions );
//      ma_Cbw_Ubw.toPairwiseStream( outfile, &fasta_in.m_descriptions );
//      cout << "Wrote out Multiple Alignment to " << filename << endl;
//      outfile.close();
//      
//      // Now try making an MSA out of it.
//      ::MSA msa_Cbw_Ubw;
//      ma_Cbw_Ubw.appendToMSA( &msa_Cbw_Ubw, &fasta_in.m_descriptions );
//      filename = filename_base + ".Cbw_Ubw.blast.a2m";
//      TextFile msa_Cbw_Ubw_outfile =
//        TextFile( filename.c_str(), true );
//      msa_Cbw_Ubw.ToFile( msa_Cbw_Ubw_outfile );
//      msa_Cbw_Ubw_outfile.Close();
//      cout << "Wrote MSA to file " << filename << endl;
//      
//      // Trying score-getting
//      cout << "Total SP score is " << getNonInsertionSPScore( msa_Cbw_Ubw ) << endl;
//      if( save_scorefiles ) {
//        filename += ".spscores";
//        g_pstrScoreFileName = filename.c_str();
//        WriteScoreFile( msa_Cbw_Ubw );
//        cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
//        g_pstrScoreFileName = NULL;
//      } // End if save_scorefiles
//      // Qscore
//      if( have_ref ) {
//        dSP = dPS = dCS = dInsane;
//        CompareMSA( msa_Cbw_Ubw, msa_ref, &dSP, &dPS, &dCS);
//        dTC = TC( msa_Cbw_Ubw, msa_ref );
//        printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
//               dSP, dPS, dCS, dTC);
//      } // End if have_ref
//      
//    } // End if do_Cbw_Ubw
//
//    if( do_Ubw ) {
//      cout << "Training again (starting over), using Unconditional Baum-Welch until convergence." << endl;
//
//      profile.copyFrom( notraining_profile );
//
//      ProfileTrainer<ProfileType, ScoreDifferenceType, ScoreType, MatrixValueType> trainer_Ubw =
//        ProfileTrainer<ProfileType, ScoreDifferenceType, ScoreType, MatrixValueType>( &profile, fasta_in, num_sequences_to_use );
//      // TODO: REMOVE:
//      //cout << "fasta_in.size() is " << fasta_in.size() << endl;
//      //cout << "m_sequence_count is " << trainer_Ubw.m_sequence_count << endl;
//      // TODO: REMOVE:
//      //trainer_Ubw.m_parameters.profileValueMinimum = 1e-2; // ?
//      trainer_Ubw.m_parameters.scorePercentChangeMinimum_position_cycle = 100;
//      trainer_Ubw.m_parameters.scorePercentChangeMinimum_iteration = 100;
//      //trainer_Ubw.m_parameters.useRabinerScaling = true;//false;
//      //trainer_Ubw.m_parameters.minIterations = 4;
//      //trainer_Ubw.m_parameters.maxIterations = 2;
//      //trainer_Ubw.m_parameters.maxPositionCycles = 2;
//      trainer_Ubw.m_parameters.maxPositionCycles_globals = 2;
//      //trainer_Ubw.m_parameters.trainGlobalsFirst = true;
//      trainer_Ubw.m_parameters.trainProfileGlobals = train_globals;
//      //trainer_Ubw.m_parameters.maxPositionCycles = 3;
//      trainer_Ubw.m_parameters.useUnconditionalBaumWelch = true;
//      //trainer_Ubw.m_parameters.debug = DEBUG_All;
//      //trainer_Ubw.m_parameters.verbosity = VERBOSITY_All;
//      trainer_Ubw.m_parameters.verbosity = VERBOSITY_Low;
//      
//      //cout << "The profile (before) is:" << endl;
//      //cout << *trainer_Ubw.m_profile << endl;
//      
//      score = trainer_Ubw.train();
//      cout << "Now (after training), the score is " << score << endl; //", and the profile is:" << endl;
//      //cout << *trainer_Ubw.m_profile << endl;
//      //ViterbiDPType::Matrix::SequentialAccessContainer viterbi_matrices =
//      //  ViterbiDPType::Matrix::SequentialAccessContainer(
//      //    *trainer_Ubw.m_profile,
//      //    trainer_Ubw.m_poly_sequences,
//      //    trainer_Ubw.m_sequence_count
//      //  );
//      // TODO: REMOVE
//      cout << "Calculating viterbi scores for the 'Ubw' profile." << endl;
//      viterbi_profile.copyPositions( *trainer_Ubw.m_profile );
//      score =
//        viterbi_dp.forward_score_viterbi(
//          trainer_Ubw.m_parameters,
//          viterbi_profile, // *trainer_Ubw.m_profile,
//          trainer_Ubw.m_poly_sequences,
//          trainer_Ubw.m_sequence_count,
//          viterbi_matrices
//        );
//      cout << "The total score for all sequences, using viterbi, is: " << score << endl;
//      ViterbiDPType::MultipleAlignment<ProfileType> ma_Ubw =
//        viterbi_dp.forward_viterbiAlign(
//          trainer_Ubw.m_parameters,
//          viterbi_profile, // *trainer_Ubw.m_profile,
//          trainer_Ubw.m_poly_sequences,
//          trainer_Ubw.m_sequence_count,
//          viterbi_matrices // TODO: ARE THESE MATRICES SIZED CORRECTLY?
//        );
//      filename = filename_base + ".Ubw.blast";
//      outfile.open( filename.c_str() );
//      //cout << "Multiple Alignment is:" << endl;
//      //ma_Ubw.toPairwiseStream( cout, &fasta_in.m_descriptions );
//      ma_Ubw.toPairwiseStream( outfile, &fasta_in.m_descriptions );
//      cout << "Wrote out Multiple Alignment to " << filename << endl;
//      outfile.close();
//      
//      // Now try making an MSA out of it.
//      ::MSA msa_Ubw;
//      ma_Ubw.appendToMSA( &msa_Ubw, &fasta_in.m_descriptions );
//      filename = filename_base + ".Ubw.blast.a2m";
//      TextFile msa_Ubw_outfile =
//        TextFile( filename.c_str(), true );
//      msa_Ubw.ToFile( msa_Ubw_outfile );
//      msa_Ubw_outfile.Close();
//      cout << "Wrote MSA to file " << filename << endl;
//      
//      // Trying score-getting
//      cout << "Total SP score is " << getNonInsertionSPScore( msa_Ubw ) << endl;
//      if( save_scorefiles ) {
//        filename += ".spscores";
//        g_pstrScoreFileName = filename.c_str();
//        WriteScoreFile( msa_Ubw );
//        cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
//        g_pstrScoreFileName = NULL;
//      } // End if save_scorefiles
//      // Qscore
//      if( have_ref ) {
//        dSP = dPS = dCS = dInsane;
//        CompareMSA( msa_Ubw, msa_ref, &dSP, &dPS, &dCS);
//        dTC = TC( msa_Ubw, msa_ref );
//        printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
//               dSP, dPS, dCS, dTC);
//      } // End if have_ref
//    } // End if do_Ubw

  } // End if do_asis

  // MARK
 if( do_split ) {
   cout << "Splitting based on the guide tree." << endl;
   cout << "g_Cluster2 is " << g_Cluster2 << "; g_Distance2 is " << g_Distance2 << "; g_Root2 is " << g_Root2 << endl;

   const unsigned uSeqCount = msa_muscle.GetSeqCount();

   // Split based on the tree
   const unsigned uInternalNodeCount = uSeqCount - 1;

   unsigned *Leaves1 = new unsigned[uSeqCount];
   unsigned *Leaves2 = new unsigned[uSeqCount];
   const unsigned uRootNodeIndex = tree.GetRootNodeIndex();
   g_uTreeSplitNode1 = uRootNodeIndex;
   g_uTreeSplitNode2 = tree.GetRight( uRootNodeIndex );

   unsigned uCount1;
   unsigned uCount2;
   GetLeaves(tree, g_uTreeSplitNode2, Leaves1, &uCount1);
   GetLeavesExcluding(tree, uRootNodeIndex, g_uTreeSplitNode2,
                      Leaves2, &uCount2);
   cout << endl << "Split at root (node " << uRootNodeIndex << ')' << endl;
   cout << "Group1 = {" << endl;
   for (unsigned n = 0; n < uCount1; ++n) {
     cout << '\t' << Leaves1[n] << '\t' << tree.GetName(Leaves1[n]) << endl;
   }
   cout << " }" << endl;
   cout << "Group2 = {" << endl;
   for (unsigned n = 0; n < uCount2; ++n) {
     cout << '\t' << Leaves2[n] << '\t' << tree.GetName(Leaves2[n]) << endl;
   }
   cout << endl;

   unsigned *Ids1 = new unsigned[uSeqCount];
   unsigned *Ids2 = new unsigned[uSeqCount];
   LeafIndexesToIds(tree, Leaves1, uCount1, Ids1);
   LeafIndexesToIds(tree, Leaves2, uCount2, Ids2);

   cout << "Getting MSA subsets msa1 and msa2" << endl;
   ::MSA msa1;
   ::MSA msa2;
   MSASubsetByIds(msa_muscle, Ids1, uCount1, msa1);
   MSASubsetByIds(msa_muscle, Ids2, uCount2, msa2);

   // ERE I AM
   filename = filename_base + ".1.a2m";
   TextFile msa1_outfile =
     TextFile( filename.c_str(), true );
   msa1.ToFile( msa1_outfile );
   msa1_outfile.Close();
   cout << "Wrote msa1 to file " << filename << endl;
   filename = filename_base + ".2.a2m";
   TextFile msa2_outfile =
     TextFile( filename.c_str(), true );
   msa2.ToFile( msa2_outfile );
   msa2_outfile.Close();
   cout << "Wrote msa2 to file " << filename << endl;

   // Set the ids.
   for (unsigned uSeqIndex = 0; uSeqIndex < uCount1; ++uSeqIndex)
     msa1.SetSeqId(uSeqIndex, Ids1[ uSeqIndex ] );
   for (unsigned uSeqIndex = 0; uSeqIndex < uCount2; ++uSeqIndex)
     msa2.SetSeqId(uSeqIndex, Ids2[ uSeqIndex ] );

   DeleteGappedCols(msa1);
   DeleteGappedCols(msa2);

   cout << "Getting SeqVect subsets v1 and v2" << endl;
   SeqVect v1;
   SeqVect v2;
   SeqVectFromMSA( msa1, v1 );
   SeqVectFromMSA( msa2, v2 );
   cout << "Got SeqVect subsets v1 and v2" << endl;

   filename = filename_base + ".1";
   TextFile v1_outfile =
     TextFile( filename.c_str(), true );
   v1.ToFile( v1_outfile );
   cout << "Wrote v1 to file " << filename << endl;
   v1_outfile.Close();
   filename = filename_base + ".2";
   TextFile v2_outfile =
     TextFile( filename.c_str(), true );
   v2.ToFile( v2_outfile );
   v2_outfile.Close();
   cout << "Wrote v2 to file " << filename << endl;

   cout << "Creating Fasta objects fasta1 and fasta2" << endl;
   Fasta<SequenceResidueType> fasta1;
   Fasta<SequenceResidueType> fasta2;
   fasta1 = v1;
   fasta2 = v2;

   if( do_muscle_after_split ) {
     cout << "Running muscle separately on v1 and v2." << endl;

     // For now, don't propagate argv, just run muscle with default options.
     filename = v1_outfile.GetFileName();
     filename += ".muscle";
     const char * const argv1[] =
       { "muscle", "-maxiters", "2", "-stable", "-noanchors", "-in", v1_outfile.GetFileName(), "-out", filename.c_str(), 0 };
     const int argc1 = 7;
     cout << "Running muscle/muscle";
     for( int arg_i = 1; arg_i < argc1; arg_i++ ) {
       cout << ' ' << argv1[ arg_i ];
     }
     cout << endl;
     muscle_pid =
       spawn( "muscle/muscle", argv1 );
     wait( &muscle_exit_status );
     cout << "Ran muscle on v1.  Exit status was " << muscle_exit_status << endl;
     
     ::MSA msa1_muscle;
     TextFile msa1_muscle_infile =
       TextFile( filename.c_str() );
     msa1_muscle.FromFile( msa1_muscle_infile );
     msa1_muscle_infile.Close();
     cout << "READ in MSA from " << msa1_muscle_infile.GetFileName() << endl;

     // For now, don't propagate argv, just run muscle with default options.
     filename = v2_outfile.GetFileName();
     filename += ".muscle";
     const char * const argv2[] =
       { "muscle", "-maxiters", "2", "-stable", "-noanchors", "-in", v2_outfile.GetFileName(), "-out", filename.c_str(), 0 };
     const int argc2 = 7;
     cout << "Running muscle/muscle";
     for( int arg_i = 1; arg_i < argc2; arg_i++ ) {
       cout << ' ' << argv2[ arg_i ];
     }
     cout << endl;
     muscle_pid =
       spawn( "muscle/muscle", argv2 );
     wait( &muscle_exit_status );
     cout << "Ran muscle on v2.  Exit status was " << muscle_exit_status << endl;
     
     ::MSA msa2_muscle;
     TextFile msa2_muscle_infile =
       TextFile( filename.c_str() );
     msa2_muscle.FromFile( msa2_muscle_infile );
     msa2_muscle_infile.Close();
     cout << "READ in MSA from " << msa2_muscle_infile.GetFileName() << endl;

     msa1_muscle.FixAlpha();
     msa2_muscle.FixAlpha();
     // Set the ids.
     for (unsigned uSeqIndex = 0; uSeqIndex < uCount1; ++uSeqIndex)
       msa1_muscle.SetSeqId(uSeqIndex, Ids1[ uSeqIndex ] );
     for (unsigned uSeqIndex = 0; uSeqIndex < uCount2; ++uSeqIndex)
       msa2_muscle.SetSeqId(uSeqIndex, Ids2[ uSeqIndex ] );

     msa1.Copy( msa1_muscle );
     msa2.Copy( msa2_muscle );
   } // End if do_muscle_after_split

   cout << "========   msa1 =======" << endl;
   cout << "Total SP score is " << getNonInsertionSPScore( msa1 ) << endl;
   if( save_scorefiles ) {
     filename = g_pstrOutFileName;
     filename += ".1.spscores";
     g_pstrScoreFileName = filename.c_str();
     WriteScoreFile( msa1 );
     cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
     g_pstrScoreFileName = NULL;
   } // End if save_scorefiles

   // hmmer
   // Need to convert from a muscle MSA to a hmmer MSA.
   hmmer::MSA * hmmer_msa1 =
     hmmer::createMSAFromMuscleMSA( msa1 );
   
   // Build HMM using hmmer's hmmbuild
   struct hmmer::plan7_s *hmm1 = hmmer::hmmbuild( hmmer_msa1 );
   if (hmm1 == NULL) hmmer::Die( "Null HMM" );          // NULL on file parse failure
   // Convert counts to probs
   if( !( hmm1->flags | PLAN7_HASPROB ) ) {
     cout << "Calling Plan7Renormalize( hmm1 )." << endl;
     hmmer::Plan7Renormalize( hmm1 );
   }
   
   cout << "Done building hmm1." << endl;
   
   //exit( 0 );
   
   ProfileType profile1;
   copyFromHMMerProfile( profile1, hmm1 );
   cout << "The profile (1) as read in from HMMer is:" << endl;
   cout << profile1 << endl;
   
   hmmer::FreePlan7(hmm1);
   
   // Profile training
   //ViterbiDPType dp =
   //  ViterbiDPType();
   //ScoreType score;
   
   uint32_t num_sequences_to_use1 = fasta1.size();
   //cout << "fasta1:" << endl;
   //if( num_sequences_to_use1 ) {
   //  fasta1.writeFasta( cout, num_sequences_to_use1 );
   //  cout << endl;
   //} else {
   //  cout << fasta1 << endl;
   //}

   // using the profile we read from hmmer
   // First fix the preAlign and postAlign, which for some reason are excessively high when we read from the hmmer profile.
   //profile1[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] =
   //  .1;
   //profile1[ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ] =
   //  1 -
   //  profile1[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ];
   //profile1[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] =
   //  .1;
   //profile1[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] =
   //  1 -
   //  profile1[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ];
   
   // TODO: REMOVE?  Testing making the indel probs very low before running viterbi.
   //double very_small = 1E-5; // profileValueMinimum
   //profile1[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ] =
   //  very_small;
   //profile1[ Transition::fromBegin ][ TransitionFromBegin::toMatch ] =
   //  ( 1 ) -
   //  profile1[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ];
   //profile1[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] =
   //  ( .5 );
   //profile1[ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ] =
   //  ( 1 ) -
   //  profile1[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ];
   //profile1[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] =
   //  ( .5 );
   //profile1[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] =
   //  ( 1 ) -
   //  profile1[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ];
   //profile1[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] =
   //  very_small;
   //profile1[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] =
   //  very_small;
   //profile1[ Transition::fromMatch ][ TransitionFromMatch::toMatch ] =
   //  ( 1 ) -
   //  profile1[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] -
   //  profile1[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ];  

   // TODO: REMOVE.  TESTING.  These are default globals from SeqanTests.cpp:
      profile1.normalize( 1E-5 ); // Ensure values aren't too tiny.

#ifndef DISALLOW_FLANKING_TRANSITIONS
      profile1[ galosh::Transition::fromPreAlign ][ galosh::TransitionFromPreAlign::toPreAlign ] =
        .01;
      profile1[ galosh::Transition::fromPreAlign ][ galosh::TransitionFromPreAlign::toBegin ] =
        1.0 -
        profile1[ galosh::Transition::fromPreAlign ][ galosh::TransitionFromPreAlign::toPreAlign ];
#endif // !DISALLOW_FLANKING_TRANSITIONS
      profile1[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletion ] =
        .01;
#ifdef USE_DEL_IN_DEL_OUT
      profile1[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletionIn ] =
        .5;
      profile1[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toMatch ] =
        1.0 -
        (
          profile1[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletion ] +
          profile1[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletionIn ]
        );
      profile1[ galosh::Transition::fromDeletionIn ][ galosh::TransitionFromDeletionIn::toDeletionIn ] =
        .99;//.5;
      profile1[ galosh::Transition::fromDeletionIn ][ galosh::TransitionFromDeletionIn::toMatch ] =
        1.0 -
        profile1[ galosh::Transition::fromDeletionIn ][ galosh::TransitionFromDeletionIn::toDeletionIn ];
#else
      profile1[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toMatch ] =
        1.0 -
        profile1[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletion ];
#endif // USE_DEL_IN_DEL_OUT .. else ..
      
      profile1[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toInsertion ] =
        .01;
      profile1[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletion ] =
        .01;
  #ifdef USE_DEL_IN_DEL_OUT
      profile1[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletionOut ] =
        .5;
      profile1[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toMatch ] =
        1.0 -
        (
          profile1[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toInsertion ] +
          profile1[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletion ] +
          profile1[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletionOut ]
        );
      profile1[ galosh::Transition::fromDeletionOut ][ galosh::TransitionFromDeletionOut::toDeletionOut ] =
        .99;//.5;
      profile1[ galosh::Transition::fromDeletionOut ][ galosh::TransitionFromDeletionOut::toEnd ] =
        1.0 -
        profile1[ galosh::Transition::fromDeletionOut ][ galosh::TransitionFromDeletionOut::toDeletionOut ];
#else
      profile1[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toMatch ] =
        1.0 -
        (
          profile1[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toInsertion ] +
          profile1[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletion ]
        );
#endif // USE_DEL_IN_DEL_OUT .. else ..
      profile1[ 0 ][ galosh::Transition::fromInsertion ][ galosh::TransitionFromInsertion::toInsertion ] =
        .5;
      profile1[ 0 ][ galosh::Transition::fromInsertion ][ galosh::TransitionFromInsertion::toMatch ] =
        1.0 -
        profile1[ 0 ][ galosh::Transition::fromInsertion ][ galosh::TransitionFromInsertion::toInsertion ];
      profile1[ 0 ][ galosh::Transition::fromDeletion ][ galosh::TransitionFromDeletion::toDeletion ] =
        .5;
      profile1[ 0 ][ galosh::Transition::fromDeletion ][ galosh::TransitionFromDeletion::toMatch ] =
        1.0 -
        profile1[ 0 ][ galosh::Transition::fromDeletion ][ galosh::TransitionFromDeletion::toDeletion ];
#ifndef DISALLOW_FLANKING_TRANSITIONS
      profile1[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toPostAlign ] =
        .01;
      profile1[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toTerminal ] =
        1.0 -
        profile1[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toPostAlign ];
#endif // !DISALLOW_FLANKING_TRANSITIONS

    if( use_lengthadjust ) {
      cout << "Initial profile length: " << profile1.length() << endl;
      cout << "Lengthadjust threshold: " << lengthadjust_threshold << endl;
      cout << "Lengthadjust threshold increment: " << lengthadjust_threshold_increment << endl;
    } else {
      cout << "Profile length: " << profile1.length() << endl;
    }

   ProfileType viterbi_profile1 = profile1;
   ProfileType notraining_profile1 = profile1;

    ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType> trainer_Cbw1( &profile1, fasta1, num_sequences_to_use1 );

   // TODO: REMOVE:
   //cout << "fasta1.size() is " << fasta1.size() << endl;
   //cout << "m_sequence_count is " << trainer_Cbw1.m_sequence_count << endl;

    //cout << "The profile (before) is:" << endl;
    //cout << *trainer_Cbw1.m_profile << endl;
    
    // TODO: REMOVE:
    //cout << "fasta_in.size() is " << fasta_in.size() << endl;
    //cout << "m_sequence_count is " << trainer_Cbw1.m_sequence_count << endl;

    // TODO: Try with the lengthadjust: maxBaumWelchInverseScalar > 0
    trainer_Cbw1.m_parameters.minIterations = 1;
    trainer_Cbw1.m_parameters.maxIterations = 1000;
    // Having maxPositionCycles_globals > 1 seems ok; takes about the same
    // number of iterations, converges to roughly the same place; takes
    // longer by virtue of having more pos cycles per iteration of course.
    trainer_Cbw1.m_parameters.maxPositionCycles = 1;
    // Having maxPositionCycles_globals > 1 seems to make convergence way
    // slower when lengthadjust is on.  Length keeps adjusting..
    trainer_Cbw1.m_parameters.maxPositionCycles_globals = 1;
    trainer_Cbw1.m_parameters.minBaumWelchInverseScalar = 0; // Straight-up bw.
    trainer_Cbw1.m_parameters.maxBaumWelchInverseScalar = 0; // Straight-up bw.
    trainer_Cbw1.m_parameters.minBaumWelchInverseScalar_globals = 0; // Straight-up bw.
    trainer_Cbw1.m_parameters.maxBaumWelchInverseScalar_globals = 0; // Straight-up bw.
    trainer_Cbw1.m_parameters.scorePercentChangeMinimum_position_cycle = 1;
    trainer_Cbw1.m_parameters.scorePercentChangeMinimum_iteration = .01;
    // TODO: Debug alwaysAccept..?  It doesn't *always*.. ok?
    trainer_Cbw1.m_parameters.alwaysAccept = false;//true;
  
    trainer_Cbw1.m_parameters.proposeProfileLengthChanges = use_lengthadjust;
    trainer_Cbw1.m_parameters.useAlignmentProfiles = true;
    trainer_Cbw1.m_parameters.numIterationsBetweenLengthChanges = 0;
    trainer_Cbw1.m_parameters.proposeDeletingThreshold =
      lengthadjust_threshold; //.01;//.025;//.1;
    trainer_Cbw1.m_parameters.proposeDeletingThreshold_increment =
      lengthadjust_threshold_increment; //.0005;//.00005;//.0005;//5E-5;//.0003125;//.00625;//.025;
    trainer_Cbw1.m_parameters.proposeInsertingThreshold =
      trainer_Cbw1.m_parameters.proposeDeletingThreshold;// / 4;//seqan::ValueSize<ResidueType>::VALUE; // TODO: Figure this out...
    trainer_Cbw1.m_parameters.proposeInsertingPreAlignThreshold = //.35; //.5;
      trainer_Cbw1.m_parameters.proposeInsertingThreshold;
    trainer_Cbw1.m_parameters.proposeInsertingPostAlignThreshold = //.35;//.5;
      trainer_Cbw1.m_parameters.proposeInsertingThreshold;
    trainer_Cbw1.m_parameters.proposeInsertingThreshold_increment =
      trainer_Cbw1.m_parameters.proposeDeletingThreshold_increment;// TODO: ERE I AM! PUT BACK? / seqan::ValueSize<ResidueType>::VALUE;
  
    //if( have_trained_profile && start_with_trained_profile ) {
    //  // When we start with the trained profile, we need to get past the
    //  // length modification wait time (which is just one iteration).
    //  trainer_Cbw1.m_parameters.minIterations =
    //    max( trainer_Cbw1.m_parameters.minIterations, ( uint32_t )2 );
    //}
  
    // Use rabiner scaling? (default true)
    // NOTE: You must change the MatrixValueType to logspace or bfloat iff this is false!\
    // TODO: I don't think I've tested this in a long long time.  Probably safest to disable it for now.  At least make the default false.
    trainer_Cbw1.m_parameters.useRabinerScaling = false;

    // Train globals first?
    //trainer_Cbw1.m_parameters.trainGlobalsFirst = true; // note: breaks UBW.

    // Use Ubw?
    trainer_Cbw1.m_parameters.useUnconditionalBaumWelch = false;//true;
    trainer_Cbw1.m_parameters.unconditionalIsolatesGlobals = false;

    trainer_Cbw1.m_parameters.trainProfileGlobals = train_globals;
    //trainer_Cbw1.m_parameters.maxPositionCycles = 3;
    //trainer_Cbw1.m_parameters.useUnconditionalBaumWelch = true;
    trainer_Cbw1.m_parameters.usePriors = use_priors;
    //trainer_Cbw1.m_parameters.debug = DEBUG_All;
    //trainer_Cbw1.m_parameters.verbosity = VERBOSITY_All;
    trainer_Cbw1.m_parameters.verbosity = VERBOSITY_Low;
    
    // Use Baldi?
    // For testing Baldi-style gradient ascent
#ifdef ALLOW_BOLTZMANN_GIBBS
    if( use_baldi_siegel ) {
      // NOTE about priors:  since globals are presently not updated using Baldi, you can still usePriors and they will affect the globals *but not the positions*.
      trainer_Cbw1.m_parameters.baldiLearningRate = 1; // 0 means noBaldi!
      trainer_Cbw1.m_parameters.baldiTemperature = 1;
      trainer_Cbw1.m_parameters.baldiHybrid = false;
      trainer_Cbw1.m_parameters.siegelMaxFindingThePeakAttempts_positions = 1000; // 0 means Baldi not Baldi / Siegel !!!
      trainer_Cbw1.m_parameters.siegelEpsilonScaleFactor = 1.5;
      trainer_Cbw1.m_parameters.siegelMaxRefiningThePeakSteps_positions = 1;//1000;
      trainer_Cbw1.m_parameters.siegelRefiningThePeakStepsConvergenceThreshold = 1E-5;
      trainer_Cbw1.m_parameters.minBaumWelchInverseScalar = 0;
      trainer_Cbw1.m_parameters.maxBaumWelchInverseScalar = 0;
      //trainer_Cbw1.m_parameters.maxPositionCycles = 10;
    } // use_baldi_siegel
#endif //ALLOW_BOLTZMANN_GIBBS

    DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer viterbi_matrices1(
                                                                                                                                                profile1,
        fasta1,
        num_sequences_to_use1
      );

    // TODO: REMOVE
    cout << "Calculating viterbi alignments and scores for the 'notraining1' profile." << endl;
    viterbi_profile1.copyPositions( *trainer_Cbw1.m_profile );
    score =
      dp.forward_score_viterbi(
          trainer_Cbw1.m_parameters,
          viterbi_profile1,
          fasta1,
          num_sequences_to_use1,
          viterbi_matrices1
        );
    cout << "The total score for all sequences, using viterbi, is: " << score << endl;
    DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::MultipleAlignment<ProfileType, SequenceResidueType> ma_notraining1(
      &viterbi_profile1, // trainer_Cbw1.m_profile,
      &fasta1,
      num_sequences_to_use1
    );
    dp.forward_viterbiAlign(
      trainer_Cbw1.m_parameters,
      viterbi_matrices1,
      ma_notraining1
    );
    cout << "Done." << endl;

    filename = filename_base + ".1.notraining.blast";
    outfile.open( filename.c_str() );
    //cout << "Multiple Alignment is:" << endl;
    //ma_notraining1.toPairwiseStream( cout, &fasta1.m_descriptions );
    ma_notraining1.toPairwiseStream( outfile, &fasta1.m_descriptions );
    //ma_notraining1.toPileupStream( outfile, &fasta1.m_descriptions );
    cout << "Wrote Multiple Alignment 1 to file " << filename << endl;
    outfile.close();

    // Now try making an MSA out of it.
    ::MSA msa_notraining1;
    ma_notraining1.appendToMSA( &msa_notraining1, &fasta1.m_descriptions );
    filename = filename_base + ".1.notraining.blast.a2m";
    TextFile msa_notraining_outfile1 =
      TextFile( filename.c_str(), true );
    msa_notraining1.ToFile( msa_notraining_outfile1 );
    msa_notraining_outfile1.Close();
    cout << "Wrote MSA 1 to file " << filename << endl;
    
    // Trying score-getting
    cout << "Total SP score (1) is " << getNonInsertionSPScore( msa_notraining1 ) << endl;
    if( save_scorefiles ) {
      filename += ".spscores";
      g_pstrScoreFileName = filename.c_str();
      WriteScoreFile( msa_notraining1 );
      cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
      g_pstrScoreFileName = NULL;
    } // End if save_scorefiles
    // Qscore
    // Commented out because the ref seqs might not be in msa_notraining1!
    //if( have_ref ) {
    //  dSP = dPS = dCS = dInsane;
    //  CompareMSA( msa_notraining1, msa_ref, &dSP, &dPS, &dCS);
    //  dTC = TC( msa_notraining1, msa_ref );
    //  printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
    //         dSP, dPS, dCS, dTC);
    //} // End if have_ref
    
    // TODO: REMOVE!
    //exit( 0 );

    cout << "Training 1, using conditional Baum-Welch until convergence." << endl;//until convergence." << endl;
    score = trainer_Cbw1.train();
    cout << "Now (after training), the score is " << score << endl; //", and the profile is:" << endl;
    //cout << *trainer_Cbw1.m_profile << endl;

    // TODO: REMOVE
    cout << "Calculating viterbi alignments and scores for the 'Cbw' profile (1)." << endl;
    viterbi_profile1.copyPositions( *trainer_Cbw1.m_profile );
    score =
      dp.forward_score_viterbi(
          trainer_Cbw1.m_parameters,
          viterbi_profile1,
          fasta1,
          num_sequences_to_use1,
          viterbi_matrices1
        );
    cout << "The total score for all sequences, using viterbi, is: " << score << endl;
    DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::MultipleAlignment<ProfileType, SequenceResidueType> ma_Cbw1(
      &viterbi_profile1, // *trainer_Cbw1.m_profile,
      &fasta1,
      num_sequences_to_use1
    );
    dp.forward_viterbiAlign(
      trainer_Cbw1.m_parameters,
      viterbi_matrices1,
      ma_Cbw1
    );
    cout << "Done." << endl;
    //cout << "The total score for all sequences, using viterbi, is: " << ma_Cbw1.calculateScore() << endl;

    filename = filename_base + ".1.Cbw.blast";
    outfile.open( filename.c_str() );
    //cout << "Multiple Alignment is:" << endl;
    //ma_Cbw1.toPairwiseStream( cout, &fasta1.m_descriptions );
    ma_Cbw1.toPairwiseStream( outfile, &fasta1.m_descriptions );
    //ma_Cbw1.toPileupStream( outfile, &fasta1.m_descriptions );
    cout << "Wrote out Multiple Alignment Cbw 1 to " << filename << endl;
    outfile.close();

    // Now try making an MSA out of it.
    ::MSA msa_Cbw1;
    ma_Cbw1.appendToMSA( &msa_Cbw1, &fasta1.m_descriptions );
    filename = filename_base + ".1.Cbw.blast.a2m";
    TextFile msa_Cbw_outfile1 =
      TextFile( filename.c_str(), true );
    msa_Cbw1.ToFile( msa_Cbw_outfile1 );
    msa_Cbw_outfile1.Close();
    cout << "Wrote MSA Cbw 1 to file " << filename << endl;
    
    // Trying score-getting
    cout << "Total SP score for Cbw 1 is " << getNonInsertionSPScore( msa_Cbw1 ) << endl;
    if( save_scorefiles ) {
      filename += ".1.Cbw.spscores";
      g_pstrScoreFileName = filename.c_str();
      WriteScoreFile( msa_Cbw1 );
      cout << "Wrote out SP score for Cbw 1 to " << g_pstrScoreFileName << endl;
      g_pstrScoreFileName = NULL;
    } // End if save_scorefiles
    // Qscore
    // Commented out because the ref seqs might not be in msa_Cbw1!
    //if( have_ref ) {
    //  dSP = dPS = dCS = dInsane;
    //  CompareMSA( msa_Cbw1, msa_ref, &dSP, &dPS, &dCS);
    //  dTC = TC( msa_Cbw1, msa_ref );
    //  printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
    //         dSP, dPS, dCS, dTC);
    //} // End if have_ref
   
//   if( do_Cbw_Ubw ) {
//     // Now train using Ubw
//     cout << "Training again using Unconditional Baum-Welch for 2 more iterations" << endl;
//     trainer_Cbw1.m_parameters.maxIterations = 2;
//     trainer_Cbw1.m_parameters.scorePercentChangeMinimum_position_cycle =
//       DEFAULT_scorePercentChangeMinimum_iteration;
//     trainer_Cbw1.m_parameters.scorePercentChangeMinimum_iteration =
//       DEFAULT_scorePercentChangeMinimum_position_cycle;
//     trainer_Cbw1.m_parameters.useUnconditionalBaumWelch = true;
//     score = trainer_Cbw1.train();
//     cout << "Now (after training again using Ubw for 2 iterations), the score is " << score << endl;//", and the profile is:" << endl;
//     //cout << *trainer_Cbw1.m_profile << endl;
//     // TODO: REMOVE
//     cout << "Calculating viterbi scores for the 'Cbw_Ubw2' profile1." << endl;
//     viterbi_profile1.copyPositions( *trainer_Cbw1.m_profile );
//
//
//     score =
//       viterbi_dp.forward_score_viterbi(
//         trainer_Cbw1.m_parameters,
//         viterbi_profile1, // *trainer_Cbw1.m_profile,
//         trainer_Cbw1.m_poly_sequences,
//         trainer_Cbw1.m_sequence_count,
//         viterbi_matrices1
//       );
//     cout << "The total score for all sequences, using viterbi, is: " << score << endl;
//     ViterbiDPType::MultipleAlignment<ProfileType> ma1_Cbw_Ubw2 =
//       viterbi_dp.forward_viterbiAlign(
//         trainer_Cbw1.m_parameters,
//         viterbi_profile1, // *trainer_Cbw1.m_profile,
//         trainer_Cbw1.m_poly_sequences,
//         trainer_Cbw1.m_sequence_count,
//         viterbi_matrices // TODO: ARE THESE MATRICES SIZED CORRECTLY?
//       );
//     filename = filename_base + ".1.Cbw_Ubw2.blast";
//     outfile.open( filename.c_str() );
//     //cout << "Multiple Alignment is:" << endl;
//     //ma1_Cbw_Ubw2.toPairwiseStream( cout, &fasta1.m_descriptions );
//     ma1_Cbw_Ubw2.toPairwiseStream( outfile, &fasta1.m_descriptions );
//     cout << "Wrote out Multiple Alignment to " << filename << endl;
//     outfile.close();
//     
//     // Now try making an MSA out of it.
//     //MSA msa1_Cbw_Ubw2;
//     ma1_Cbw_Ubw2.appendToMSA( &msa1_Cbw_Ubw2, &fasta1.m_descriptions );
//     filename = filename_base + ".1.Cbw_Ubw2.blast.a2m";
//     TextFile msa1_Cbw_Ubw2_outfile =
//       TextFile( filename.c_str(), true );
//     msa1_Cbw_Ubw2.ToFile( msa1_Cbw_Ubw2_outfile );
//     msa1_Cbw_Ubw2_outfile.Close();
//     cout << "Wrote MSA to file " << filename << endl;
//     
//     // Trying score-getting
//     cout << "Total SP score is " << getNonInsertionSPScore( msa1_Cbw_Ubw2 ) << endl;
//     if( save_scorefiles ) {
//       filename += ".spscores";
//       g_pstrScoreFileName = filename.c_str();
//       WriteScoreFile( msa1_Cbw_Ubw2 );
//       cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
//       g_pstrScoreFileName = NULL;
//     } // End if save_scorefiles
//     
//     cout << "Training yet again using Unconditional Baum-Welch for 2 more iterations" << endl;
//     trainer_Cbw1.m_parameters.maxIterations = 2;
//     trainer_Cbw1.m_parameters.useUnconditionalBaumWelch = true;
//     score = trainer_Cbw1.train();
//     cout << "Now (after training again using Ubw for a total of 4 iterations max), the score is " << score << endl;//", and the profile is:" << endl;
//     //cout << *trainer_Cbw1.m_profile << endl;
//     // TODO: REMOVE
//     cout << "Calculating viterbi scores for the 'Cbw_Ubw4' profile1." << endl;
//     viterbi_profile1.copyPositions( *trainer_Cbw1.m_profile );
//     score =
//       viterbi_dp.forward_score_viterbi(
//         trainer_Cbw1.m_parameters,
//         viterbi_profile1, // *trainer_Cbw1.m_profile,
//         trainer_Cbw1.m_poly_sequences,
//         trainer_Cbw1.m_sequence_count,
//         viterbi_matrices
//       );
//     cout << "The total score for all sequences, using viterbi, is: " << score << endl;
//     ViterbiDPType::MultipleAlignment<ProfileType> ma1_Cbw_Ubw4 =
//       viterbi_dp.forward_viterbiAlign(
//         trainer_Cbw1.m_parameters,
//         viterbi_profile1, // *trainer_Cbw1.m_profile,
//         trainer_Cbw1.m_poly_sequences,
//         trainer_Cbw1.m_sequence_count,
//         viterbi_matrices // TODO: ARE THESE MATRICES SIZED CORRECTLY?
//       );
//     filename = filename_base + ".1.Cbw_Ubw4.blast";
//     outfile.open( filename.c_str() );
//     //cout << "Multiple Alignment is:" << endl;
//     //ma1_Cbw_Ubw4.toPairwiseStream( cout, &fasta1.m_descriptions );
//     ma1_Cbw_Ubw4.toPairwiseStream( outfile, &fasta1.m_descriptions );
//     cout << "Wrote out Multiple Alignment to " << filename << endl;
//     outfile.close();
//     
//     // Now try making an MSA out of it.
//     //MSA msa1_Cbw_Ubw4;
//     ma1_Cbw_Ubw4.appendToMSA( &msa1_Cbw_Ubw4, &fasta1.m_descriptions );
//     filename = filename_base + ".1.Cbw_Ubw4.blast.a2m";
//     TextFile msa1_Cbw_Ubw4_outfile =
//       TextFile( filename.c_str(), true );
//     msa1_Cbw_Ubw4.ToFile( msa1_Cbw_Ubw4_outfile );
//     msa1_Cbw_Ubw4_outfile.Close();
//     cout << "Wrote MSA to file " << filename << endl;
//     
//     // Trying score-getting
//     cout << "Total SP score is " << getNonInsertionSPScore( msa1_Cbw_Ubw4 ) << endl;
//     if( save_scorefiles ) {
//       filename += ".spscores";
//       g_pstrScoreFileName = filename.c_str();
//       WriteScoreFile( msa1_Cbw_Ubw4 );
//       cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
//       g_pstrScoreFileName = NULL;
//     } // End if save_scorefiles
//
//     cout << "Training yet again using Unconditional Baum-Welch until convergence." << endl;
//     trainer_Cbw1.m_parameters.maxIterations = 100;
//     trainer_Cbw1.m_parameters.scorePercentChangeMinimum_position_cycle = 100;
//     trainer_Cbw1.m_parameters.scorePercentChangeMinimum_iteration = 100;
//     trainer_Cbw1.m_parameters.useUnconditionalBaumWelch = true;
//     score = trainer_Cbw1.train();
//     cout << "Now (after training again using Ubw until convergence), the score is " << score << endl;//", and the profile is:" << endl;
//     //cout << *trainer_Cbw1.m_profile << endl;
//     // TODO: REMOVE
//     cout << "Calculating viterbi scores for the 'Cbw_Ubw' profile1." << endl;
//     viterbi_profile1.copyPositions( *trainer_Cbw1.m_profile );
//     score =
//       viterbi_dp.forward_score_viterbi(
//         trainer_Cbw1.m_parameters,
//         viterbi_profile1, // *trainer_Cbw1.m_profile,
//         trainer_Cbw1.m_poly_sequences,
//         trainer_Cbw1.m_sequence_count,
//         viterbi_matrices
//       );
//     cout << "The total score for all sequences, using viterbi, is: " << score << endl;
//     ViterbiDPType::MultipleAlignment<ProfileType> ma1_Cbw_Ubw =
//       viterbi_dp.forward_viterbiAlign(
//         trainer_Cbw1.m_parameters,
//         viterbi_profile1, // *trainer_Cbw1.m_profile,
//         trainer_Cbw1.m_poly_sequences,
//         trainer_Cbw1.m_sequence_count,
//         viterbi_matrices // TODO: ARE THESE MATRICES SIZED CORRECTLY?
//       );
//     filename = filename_base + ".1.Cbw_Ubw.blast";
//     outfile.open( filename.c_str() );
//     //cout << "Multiple Alignment is:" << endl;
//     //ma1_Cbw_Ubw.toPairwiseStream( cout, &fasta1.m_descriptions );
//     ma1_Cbw_Ubw.toPairwiseStream( outfile, &fasta1.m_descriptions );
//     cout << "Wrote out Multiple Alignment to " << filename << endl;
//     outfile.close();
//     
//     // Now try making an MSA out of it.
//     //MSA msa1_Cbw_Ubw;
//     ma1_Cbw_Ubw.appendToMSA( &msa1_Cbw_Ubw, &fasta1.m_descriptions );
//     filename = filename_base + ".1.Cbw_Ubw.blast.a2m";
//     TextFile msa1_Cbw_Ubw_outfile =
//       TextFile( filename.c_str(), true );
//     msa1_Cbw_Ubw.ToFile( msa1_Cbw_Ubw_outfile );
//     msa1_Cbw_Ubw_outfile.Close();
//     cout << "Wrote MSA to file " << filename << endl;
//     
//     // Trying score-getting
//     cout << "Total SP score is " << getNonInsertionSPScore( msa1_Cbw_Ubw ) << endl;
//     if( save_scorefiles ) {
//       filename += ".spscores";
//       g_pstrScoreFileName = filename.c_str();
//       WriteScoreFile( msa1_Cbw_Ubw );
//       cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
//       g_pstrScoreFileName = NULL;
//     } // End if save_scorefiles
//     
//   } // End if do_Cbw_Ubw

//   if( do_Ubw ) {
//     cout << "Training again (starting over), using Unconditional Baum-Welch until convergence." << endl;
//
//     profile1.copyFrom( notraining_profile1 );
//
//     ProfileTrainer<ProfileType, ScoreDifferenceType, ScoreType, MatrixValueType> trainer1_Ubw =
//       ProfileTrainer<ProfileType, ScoreDifferenceType, ScoreType, MatrixValueType>( &profile1, fasta1, num_sequences_to_use1 );
//     // TODO: REMOVE:
//     //cout << "fasta1.size() is " << fasta1.size() << endl;
//     //cout << "m_sequence_count is " << trainer1_Ubw.m_sequence_count << endl;
//     // TODO: REMOVE:
//     //trainer1_Ubw.m_parameters.profileValueMinimum = 1e-2; // ?
//     trainer1_Ubw.m_parameters.scorePercentChangeMinimum_position_cycle = 100;
//     trainer1_Ubw.m_parameters.scorePercentChangeMinimum_iteration = 100;
//     //trainer1_Ubw.m_parameters.useRabinerScaling = false;
//     //trainer1_Ubw.m_parameters.minIterations = 4;
//     //trainer1_Ubw.m_parameters.maxIterations = 2;
//     //trainer1_Ubw.m_parameters.maxPositionCycles = 2;
//     //trainer1_Ubw.m_parameters.maxPositionCycles_globals = 2;
//     //trainer1_Ubw.m_parameters.trainGlobalsFirst = true;
//     //trainer1_Ubw.m_parameters.useRabinerScaling = false;
//     trainer1_Ubw.m_parameters.usePriors = use_priors;
//     trainer1_Ubw.m_parameters.trainProfileGlobals = train_globals;
//     //trainer1_Ubw.m_parameters.maxPositionCycles = 3;
//     trainer1_Ubw.m_parameters.useUnconditionalBaumWelch = true;
//     //trainer1_Ubw.m_parameters.debug = DEBUG_All;
//     //trainer1_Ubw.m_parameters.verbosity = VERBOSITY_All;
//     trainer1_Ubw.m_parameters.verbosity = VERBOSITY_Low;
//     
//     //cout << "The profile (before) is:" << endl;
//     //cout << *trainer1_Ubw.m_profile << endl;
//     
//     score = trainer1_Ubw.train();
//     cout << "Now (after training), the score is " << score << endl; //", and the profile is:" << endl;
//     //cout << *trainer1_Ubw.m_profile << endl;
//     //ViterbiDPType::Matrix::SequentialAccessContainer viterbi_matrices =
//     //  ViterbiDPType::Matrix::SequentialAccessContainer(
//     //    *trainer1_Ubw.m_profile,
//     //    trainer1_Ubw.m_poly_sequences,
//     //    trainer1_Ubw.m_sequence_count
//     //  );
//     // TODO: REMOVE
//     cout << "Calculating viterbi scores for the 'Ubw' profile1." << endl;
//     viterbi_profile1.copyPositions( *trainer1_Ubw.m_profile );
//     score =
//       viterbi_dp.forward_score_viterbi(
//         trainer1_Ubw.m_parameters,
//         viterbi_profile1, // *trainer1_Ubw.m_profile,
//         trainer1_Ubw.m_poly_sequences,
//         trainer1_Ubw.m_sequence_count,
//         viterbi_matrices
//       );
//     cout << "The total score for all sequences, using viterbi, is: " << score << endl;
//     ViterbiDPType::MultipleAlignment<ProfileType> ma1_Ubw =
//       viterbi_dp.forward_viterbiAlign(
//         trainer1_Ubw.m_parameters,
//         viterbi_profile1, // *trainer1_Ubw.m_profile,
//         trainer1_Ubw.m_poly_sequences,
//         trainer1_Ubw.m_sequence_count,
//         viterbi_matrices // TODO: ARE THESE MATRICES SIZED CORRECTLY?
//       );
//     filename = filename_base + ".1.Ubw.blast";
//     outfile.open( filename.c_str() );
//     //cout << "Multiple Alignment is:" << endl;
//     //ma1_Ubw.toPairwiseStream( cout, &fasta1.m_descriptions );
//     ma1_Ubw.toPairwiseStream( outfile, &fasta1.m_descriptions );
//     cout << "Wrote out Multiple Alignment to " << filename << endl;
//     outfile.close();
//     
//     // Now try making an MSA out of it.
//     //MSA msa1_Ubw;
//     ma1_Ubw.appendToMSA( &msa1_Ubw, &fasta1.m_descriptions );
//     filename = filename_base + ".1.Ubw.blast.a2m";
//     TextFile msa1_Ubw_outfile =
//       TextFile( filename.c_str(), true );
//     msa1_Ubw.ToFile( msa1_Ubw_outfile );
//     msa1_Ubw_outfile.Close();
//     cout << "Wrote MSA to file " << filename << endl;
//     
//     // Trying score-getting
//     cout << "Total SP score is " << getNonInsertionSPScore( msa1_Ubw ) << endl;
//     if( save_scorefiles ) {
//       filename += ".spscores";
//       g_pstrScoreFileName = filename.c_str();
//       WriteScoreFile( msa1_Ubw );
//       cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
//       g_pstrScoreFileName = NULL;
//     } // End if save_scorefiles
//   } // End if do_Ubw

// mark
   cout << "========   msa2 =======" << endl;
   cout << "Total SP score is " << getNonInsertionSPScore( msa2 ) << endl;
   if( save_scorefiles ) {
     filename = g_pstrOutFileName;
     filename += ".2.spscores";
     g_pstrScoreFileName = filename.c_str();
     WriteScoreFile( msa2 );
     cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
     g_pstrScoreFileName = NULL;
   } // End if save_scorefiles

   // hmmer
   // Need to convert from a muscle MSA to a hmmer MSA.
   hmmer::MSA * hmmer_msa2 =
     hmmer::createMSAFromMuscleMSA( msa2 );
   
   // Build HMM using hmmer's hmmbuild
   struct hmmer::plan7_s *hmm2 = hmmer::hmmbuild( hmmer_msa2 );
   if (hmm2 == NULL) hmmer::Die( "Null HMM" );          // NULL on file parse failure
   // Convert counts to probs
   if( !( hmm2->flags | PLAN7_HASPROB ) ) {
     cout << "Calling Plan7Renormalize( hmm2 )." << endl;
     hmmer::Plan7Renormalize( hmm2 );
   }
   
   cout << "Done building hmm2." << endl;
   
   //exit( 0 );
   
   ProfileType profile2;
   copyFromHMMerProfile( profile2, hmm2 );
   cout << "The profile (2) as read in from HMMer is:" << endl;
   cout << profile2 << endl;
   
   hmmer::FreePlan7(hmm2);
   
   // Profile training
   //ViterbiDPType dp =
   //  ViterbiDPType();
   //ScoreType score;
   
   uint32_t num_sequences_to_use2 = fasta2.size();
   //cout << "fasta2:" << endl;
   //if( num_sequences_to_use2 ) {
   //  fasta2.writeFasta( cout, num_sequences_to_use2 );
   //  cout << endl;
   //} else {
   //  cout << fasta2 << endl;
   //}

   // using the profile we read from hmmer
   // First fix the preAlign and postAlign, which for some reason are excessively high when we read from the hmmer profile.
   //profile2[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] =
   //  .1;
   //profile2[ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ] =
   //  1 -
   //  profile2[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ];
   //profile2[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] =
   //  .1;
   //profile2[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] =
   //  1 -
   //  profile2[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ];
   
   //// TODO: REMOVE?  Testing making the indel probs very low before running viterbi.
   ////double very_small = 1E-5; // profileValueMinimum
   //profile2[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ] =
   //  very_small;
   //profile2[ Transition::fromBegin ][ TransitionFromBegin::toMatch ] =
   //  ( 1 ) -
   //  profile2[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ];
   //profile2[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] =
   //  ( .5 );
   //profile2[ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ] =
   //  ( 1 ) -
   //  profile2[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ];
   //profile2[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] =
   //  ( .5 );
   //profile2[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] =
   //  ( 1 ) -
   //  profile2[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ];
   //profile2[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] =
   //  very_small;
   //profile2[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] =
   //  very_small;
   //profile2[ Transition::fromMatch ][ TransitionFromMatch::toMatch ] =
   //  ( 1 ) -
   //  profile2[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] -
   //  profile2[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ];  


   // TODO: REMOVE.  TESTING.  These are default globals from SeqanTests.cpp:
      profile2.normalize( 1E-5 ); // Ensure values aren't too tiny.

#ifndef DISALLOW_FLANKING_TRANSITIONS
      profile2[ galosh::Transition::fromPreAlign ][ galosh::TransitionFromPreAlign::toPreAlign ] =
        .01;
      profile2[ galosh::Transition::fromPreAlign ][ galosh::TransitionFromPreAlign::toBegin ] =
        1.0 -
        profile2[ galosh::Transition::fromPreAlign ][ galosh::TransitionFromPreAlign::toPreAlign ];
#endif // !DISALLOW_FLANKING_TRANSITIONS
      profile2[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletion ] =
        .01;
#ifdef USE_DEL_IN_DEL_OUT
      profile2[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletionIn ] =
        .5;
      profile2[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toMatch ] =
        1.0 -
        (
          profile2[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletion ] +
          profile2[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletionIn ]
        );
      profile2[ galosh::Transition::fromDeletionIn ][ galosh::TransitionFromDeletionIn::toDeletionIn ] =
        .99;//.5;
      profile2[ galosh::Transition::fromDeletionIn ][ galosh::TransitionFromDeletionIn::toMatch ] =
        1.0 -
        profile2[ galosh::Transition::fromDeletionIn ][ galosh::TransitionFromDeletionIn::toDeletionIn ];
#else
      profile2[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toMatch ] =
        1.0 -
        profile2[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletion ];
#endif // USE_DEL_IN_DEL_OUT .. else ..
      
      profile2[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toInsertion ] =
        .01;
      profile2[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletion ] =
        .01;
  #ifdef USE_DEL_IN_DEL_OUT
      profile2[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletionOut ] =
        .5;
      profile2[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toMatch ] =
        1.0 -
        (
          profile2[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toInsertion ] +
          profile2[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletion ] +
          profile2[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletionOut ]
        );
      profile2[ galosh::Transition::fromDeletionOut ][ galosh::TransitionFromDeletionOut::toDeletionOut ] =
        .99;//.5;
      profile2[ galosh::Transition::fromDeletionOut ][ galosh::TransitionFromDeletionOut::toEnd ] =
        1.0 -
        profile2[ galosh::Transition::fromDeletionOut ][ galosh::TransitionFromDeletionOut::toDeletionOut ];
#else
      profile2[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toMatch ] =
        1.0 -
        (
          profile2[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toInsertion ] +
          profile2[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletion ]
        );
#endif // USE_DEL_IN_DEL_OUT .. else ..
      profile2[ 0 ][ galosh::Transition::fromInsertion ][ galosh::TransitionFromInsertion::toInsertion ] =
        .5;
      profile2[ 0 ][ galosh::Transition::fromInsertion ][ galosh::TransitionFromInsertion::toMatch ] =
        1.0 -
        profile2[ 0 ][ galosh::Transition::fromInsertion ][ galosh::TransitionFromInsertion::toInsertion ];
      profile2[ 0 ][ galosh::Transition::fromDeletion ][ galosh::TransitionFromDeletion::toDeletion ] =
        .5;
      profile2[ 0 ][ galosh::Transition::fromDeletion ][ galosh::TransitionFromDeletion::toMatch ] =
        1.0 -
        profile2[ 0 ][ galosh::Transition::fromDeletion ][ galosh::TransitionFromDeletion::toDeletion ];
#ifndef DISALLOW_FLANKING_TRANSITIONS
      profile2[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toPostAlign ] =
        .01;
      profile2[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toTerminal ] =
        1.0 -
        profile2[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toPostAlign ];
#endif // !DISALLOW_FLANKING_TRANSITIONS

    if( use_lengthadjust ) {
      cout << "Initial profile length: " << profile2.length() << endl;
      cout << "Lengthadjust threshold: " << lengthadjust_threshold << endl;
      cout << "Lengthadjust threshold increment: " << lengthadjust_threshold_increment << endl;
    } else {
      cout << "Profile length: " << profile2.length() << endl;
    }

   ProfileType viterbi_profile2 = profile2;
   ProfileType notraining_profile2 = profile2;

    ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType> trainer_Cbw2( &profile2, fasta2, num_sequences_to_use2 );

   // TODO: REMOVE:
   //cout << "fasta2.size() is " << fasta2.size() << endl;
   //cout << "m_sequence_count is " << trainer_Cbw2.m_sequence_count << endl;

    //cout << "The profile (before) is:" << endl;
    //cout << *trainer_Cbw2.m_profile << endl;
    
    // TODO: REMOVE:
    //cout << "fasta_in.size() is " << fasta_in.size() << endl;
    //cout << "m_sequence_count is " << trainer_Cbw2.m_sequence_count << endl;

    // TODO: Try with the lengthadjust: maxBaumWelchInverseScalar > 0
    trainer_Cbw2.m_parameters.minIterations = 1;
    trainer_Cbw2.m_parameters.maxIterations = 1000;
    // Having maxPositionCycles_globals > 1 seems ok; takes about the same
    // number of iterations, converges to roughly the same place; takes
    // longer by virtue of having more pos cycles per iteration of course.
    trainer_Cbw2.m_parameters.maxPositionCycles = 1;
    // Having maxPositionCycles_globals > 1 seems to make convergence way
    // slower when lengthadjust is on.  Length keeps adjusting..
    trainer_Cbw2.m_parameters.maxPositionCycles_globals = 1;
    trainer_Cbw2.m_parameters.minBaumWelchInverseScalar = 0; // Straight-up bw.
    trainer_Cbw2.m_parameters.maxBaumWelchInverseScalar = 0; // Straight-up bw.
    trainer_Cbw2.m_parameters.minBaumWelchInverseScalar_globals = 0; // Straight-up bw.
    trainer_Cbw2.m_parameters.maxBaumWelchInverseScalar_globals = 0; // Straight-up bw.
    trainer_Cbw2.m_parameters.scorePercentChangeMinimum_position_cycle = 1;
    trainer_Cbw2.m_parameters.scorePercentChangeMinimum_iteration = .01;
    // TODO: Debug alwaysAccept..?  It doesn't *always*.. ok?
    trainer_Cbw2.m_parameters.alwaysAccept = false;//true;
  
    trainer_Cbw2.m_parameters.proposeProfileLengthChanges = use_lengthadjust;
    trainer_Cbw2.m_parameters.useAlignmentProfiles = true;
    trainer_Cbw2.m_parameters.numIterationsBetweenLengthChanges = 0;
    trainer_Cbw2.m_parameters.proposeDeletingThreshold =
      lengthadjust_threshold; //.01;//.025;//.1;
    trainer_Cbw2.m_parameters.proposeDeletingThreshold_increment =
      lengthadjust_threshold_increment; //.0005;//.00005;//.0005;//5E-5;//.0003125;//.00625;//.025;
    trainer_Cbw2.m_parameters.proposeInsertingThreshold =
      trainer_Cbw2.m_parameters.proposeDeletingThreshold;// / 4;//seqan::ValueSize<ResidueType>::VALUE; // TODO: Figure this out...
    trainer_Cbw2.m_parameters.proposeInsertingPreAlignThreshold = //.35; //.5;
      trainer_Cbw2.m_parameters.proposeInsertingThreshold;
    trainer_Cbw2.m_parameters.proposeInsertingPostAlignThreshold = //.35;//.5;
      trainer_Cbw2.m_parameters.proposeInsertingThreshold;
    trainer_Cbw2.m_parameters.proposeInsertingThreshold_increment =
      trainer_Cbw2.m_parameters.proposeDeletingThreshold_increment;// TODO: ERE I AM! PUT BACK? / seqan::ValueSize<ResidueType>::VALUE;
  
    //if( have_trained_profile && start_with_trained_profile ) {
    //  // When we start with the trained profile, we need to get past the
    //  // length modification wait time (which is just one iteration).
    //  trainer_Cbw2.m_parameters.minIterations =
    //    max( trainer_Cbw2.m_parameters.minIterations, ( uint32_t )2 );
    //}
  
    // Use rabiner scaling? (default true)
    // NOTE: You must change the MatrixValueType to logspace or bfloat iff this is false!\
    // TODO: I don't think I've tested this in a long long time.  Probably safest to disable it for now.  At least make the default false.
    trainer_Cbw2.m_parameters.useRabinerScaling = false;

    // Train globals first?
    //trainer_Cbw2.m_parameters.trainGlobalsFirst = true; // note: breaks UBW.

    // Use Ubw?
    trainer_Cbw2.m_parameters.useUnconditionalBaumWelch = false;//true;
    trainer_Cbw2.m_parameters.unconditionalIsolatesGlobals = false;

    trainer_Cbw2.m_parameters.trainProfileGlobals = train_globals;
    //trainer_Cbw2.m_parameters.maxPositionCycles = 3;
    //trainer_Cbw2.m_parameters.useUnconditionalBaumWelch = true;
    trainer_Cbw2.m_parameters.usePriors = use_priors;
    //trainer_Cbw2.m_parameters.debug = DEBUG_All;
    //trainer_Cbw2.m_parameters.verbosity = VERBOSITY_All;
    trainer_Cbw2.m_parameters.verbosity = VERBOSITY_Low;
    
    // Use Baldi?
    // For testing Baldi-style gradient ascent
#ifdef ALLOW_BOLTZMANN_GIBBS
    if( use_baldi_siegel ) {
      // NOTE about priors:  since globals are presently not updated using Baldi, you can still usePriors and they will affect the globals *but not the positions*.
      trainer_Cbw2.m_parameters.baldiLearningRate = 1; // 0 means noBaldi!
      trainer_Cbw2.m_parameters.baldiTemperature = 1;
      trainer_Cbw2.m_parameters.baldiHybrid = false;
      trainer_Cbw2.m_parameters.siegelMaxFindingThePeakAttempts_positions = 1000; // 0 means Baldi not Baldi / Siegel !!!
      trainer_Cbw2.m_parameters.siegelEpsilonScaleFactor = 1.5;
      trainer_Cbw2.m_parameters.siegelMaxRefiningThePeakSteps_positions = 1;//1000;
      trainer_Cbw2.m_parameters.siegelRefiningThePeakStepsConvergenceThreshold = 1E-5;
      trainer_Cbw2.m_parameters.minBaumWelchInverseScalar = 0;
      trainer_Cbw2.m_parameters.maxBaumWelchInverseScalar = 0;
      //trainer_Cbw2.m_parameters.maxPositionCycles = 10;
    } // use_baldi_siegel
#endif //ALLOW_BOLTZMANN_GIBBS

    DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer viterbi_matrices2(
                                                                                                                                                profile2,
        fasta2,
        num_sequences_to_use2
      );

    // TODO: REMOVE
    cout << "Calculating viterbi alignments and scores for the 'notraining2' profile." << endl;
    viterbi_profile2.copyPositions( *trainer_Cbw2.m_profile );
    score =
      dp.forward_score_viterbi(
          trainer_Cbw2.m_parameters,
          viterbi_profile2,
          fasta2,
          num_sequences_to_use2,
          viterbi_matrices2
        );
    cout << "The total score for all sequences, using viterbi, is: " << score << endl;
    DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::MultipleAlignment<ProfileType, SequenceResidueType> ma_notraining2(
      &viterbi_profile2, // trainer_Cbw2.m_profile,
      &fasta2,
      num_sequences_to_use2
    );
    dp.forward_viterbiAlign(
      trainer_Cbw2.m_parameters,
      viterbi_matrices2,
      ma_notraining2
    );
    cout << "Done." << endl;

    filename = filename_base + ".2.notraining.blast";
    outfile.open( filename.c_str() );
    //cout << "Multiple Alignment is:" << endl;
    //ma_notraining2.toPairwiseStream( cout, &fasta2.m_descriptions );
    ma_notraining2.toPairwiseStream( outfile, &fasta2.m_descriptions );
    //ma_notraining2.toPileupStream( outfile, &fasta2.m_descriptions );
    cout << "Wrote Multiple Alignment 2 to file " << filename << endl;
    outfile.close();

    // Now try making an MSA out of it.
    ::MSA msa_notraining2;
    ma_notraining2.appendToMSA( &msa_notraining2, &fasta2.m_descriptions );
    filename = filename_base + ".2.notraining.blast.a2m";
    TextFile msa_notraining_outfile2 =
      TextFile( filename.c_str(), true );
    msa_notraining2.ToFile( msa_notraining_outfile2 );
    msa_notraining_outfile2.Close();
    cout << "Wrote MSA 2 to file " << filename << endl;
    
    // Trying score-getting
    cout << "Total SP score (2) is " << getNonInsertionSPScore( msa_notraining2 ) << endl;
    if( save_scorefiles ) {
      filename += ".2.spscores";
      g_pstrScoreFileName = filename.c_str();
      WriteScoreFile( msa_notraining2 );
      cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
      g_pstrScoreFileName = NULL;
    } // End if save_scorefiles
    // Qscore
    // Commented out because the ref seqs might not be in msa_notraining2!
    //if( have_ref ) {
    //  dSP = dPS = dCS = dInsane;
    //  CompareMSA( msa_notraining2, msa_ref, &dSP, &dPS, &dCS);
    //  dTC = TC( msa_notraining2, msa_ref );
    //  printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
    //         dSP, dPS, dCS, dTC);
    //} // End if have_ref
    
    // TODO: REMOVE!
    //exit( 0 );

    cout << "Training 2, using conditional Baum-Welch until convergence." << endl;//until convergence." << endl;
    score = trainer_Cbw2.train();
    cout << "Now (after training), the score is " << score << endl; //", and the profile is:" << endl;
    //cout << *trainer_Cbw2.m_profile << endl;

    // TODO: REMOVE
    cout << "Calculating viterbi alignments and scores for the 'Cbw' profile (2)." << endl;
    viterbi_profile2.copyPositions( *trainer_Cbw2.m_profile );
    score =
      dp.forward_score_viterbi(
          trainer_Cbw2.m_parameters,
          viterbi_profile2,
          fasta2,
          num_sequences_to_use2,
          viterbi_matrices2
        );
    cout << "The total score for all sequences, using viterbi, is: " << score << endl;
    DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::MultipleAlignment<ProfileType, SequenceResidueType> ma_Cbw2(
      &viterbi_profile2, // *trainer_Cbw2.m_profile,
      &fasta2,
      num_sequences_to_use2
    );
    dp.forward_viterbiAlign(
      trainer_Cbw2.m_parameters,
      viterbi_matrices2,
      ma_Cbw2
    );
    cout << "Done." << endl;
    //cout << "The total score for all sequences, using viterbi, is: " << ma_Cbw2.calculateScore() << endl;

    filename = filename_base + ".2.Cbw.blast";
    outfile.open( filename.c_str() );
    //cout << "Multiple Alignment is:" << endl;
    //ma_Cbw2.toPairwiseStream( cout, &fasta2.m_descriptions );
    ma_Cbw2.toPairwiseStream( outfile, &fasta2.m_descriptions );
    //ma_Cbw2.toPileupStream( outfile, &fasta2.m_descriptions );
    cout << "Wrote out Multiple Alignment Cbw 2 to " << filename << endl;
    outfile.close();

    // Now try making an MSA out of it.
    ::MSA msa_Cbw2;
    ma_Cbw2.appendToMSA( &msa_Cbw2, &fasta2.m_descriptions );
    filename = filename_base + ".2.Cbw.blast.a2m";
    TextFile msa_Cbw_outfile2 =
      TextFile( filename.c_str(), true );
    msa_Cbw2.ToFile( msa_Cbw_outfile2 );
    msa_Cbw_outfile2.Close();
    cout << "Wrote MSA Cbw 2 to file " << filename << endl;
    
    // Trying score-getting
    cout << "Total SP score for Cbw 2 is " << getNonInsertionSPScore( msa_Cbw2 ) << endl;
    if( save_scorefiles ) {
      filename += ".2.Cbw.spscores";
      g_pstrScoreFileName = filename.c_str();
      WriteScoreFile( msa_Cbw2 );
      cout << "Wrote out SP score for Cbw 2 to " << g_pstrScoreFileName << endl;
      g_pstrScoreFileName = NULL;
    } // End if save_scorefiles
    // Qscore
    // Commented out because the ref seqs might not be in msa_Cbw2!
    //if( have_ref ) {
    //  dSP = dPS = dCS = dInsane;
    //  CompareMSA( msa_Cbw2, msa_ref, &dSP, &dPS, &dCS);
    //  dTC = TC( msa_Cbw2, msa_ref );
    //  printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
    //         dSP, dPS, dCS, dTC);
    //} // End if have_ref
   
//   if( do_Cbw_Ubw ) {
//     // Now train using Ubw
//     cout << "Training again using Unconditional Baum-Welch for 2 more iterations" << endl;
//     trainer_Cbw2.m_parameters.maxIterations = 2;
//     trainer_Cbw2.m_parameters.scorePercentChangeMinimum_position_cycle =
//       DEFAULT_scorePercentChangeMinimum_iteration;
//     trainer_Cbw2.m_parameters.scorePercentChangeMinimum_iteration =
//       DEFAULT_scorePercentChangeMinimum_position_cycle;
//     trainer_Cbw2.m_parameters.useUnconditionalBaumWelch = true;
//     score = trainer_Cbw2.train();
//     cout << "Now (after training again using Ubw for 2 iterations), the score is " << score << endl;//", and the profile is:" << endl;
//     //cout << *trainer_Cbw2.m_profile << endl;
//     // TODO: REMOVE
//     cout << "Calculating viterbi scores for the 'Cbw_Ubw2' profile2." << endl;
//     viterbi_profile2.copyPositions( *trainer_Cbw2.m_profile );
//
//
//     score =
//       viterbi_dp.forward_score_viterbi(
//         trainer_Cbw2.m_parameters,
//         viterbi_profile2, // *trainer_Cbw2.m_profile,
//         trainer_Cbw2.m_poly_sequences,
//         trainer_Cbw2.m_sequence_count,
//         viterbi_matrices2
//       );
//     cout << "The total score for all sequences, using viterbi, is: " << score << endl;
//     ViterbiDPType::MultipleAlignment<ProfileType> ma2_Cbw_Ubw2 =
//       viterbi_dp.forward_viterbiAlign(
//         trainer_Cbw2.m_parameters,
//         viterbi_profile2, // *trainer_Cbw2.m_profile,
//         trainer_Cbw2.m_poly_sequences,
//         trainer_Cbw2.m_sequence_count,
//         viterbi_matrices // TODO: ARE THESE MATRICES SIZED CORRECTLY?
//       );
//     filename = filename_base + ".2.Cbw_Ubw2.blast";
//     outfile.open( filename.c_str() );
//     //cout << "Multiple Alignment is:" << endl;
//     //ma2_Cbw_Ubw2.toPairwiseStream( cout, &fasta2.m_descriptions );
//     ma2_Cbw_Ubw2.toPairwiseStream( outfile, &fasta2.m_descriptions );
//     cout << "Wrote out Multiple Alignment to " << filename << endl;
//     outfile.close();
//     
//     // Now try making an MSA out of it.
//     //MSA msa2_Cbw_Ubw2;
//     ma2_Cbw_Ubw2.appendToMSA( &msa2_Cbw_Ubw2, &fasta2.m_descriptions );
//     filename = filename_base + ".2.Cbw_Ubw2.blast.a2m";
//     TextFile msa2_Cbw_Ubw2_outfile =
//       TextFile( filename.c_str(), true );
//     msa2_Cbw_Ubw2.ToFile( msa2_Cbw_Ubw2_outfile );
//     msa2_Cbw_Ubw2_outfile.Close();
//     cout << "Wrote MSA to file " << filename << endl;
//     
//     // Trying score-getting
//     cout << "Total SP score is " << getNonInsertionSPScore( msa2_Cbw_Ubw2 ) << endl;
//     if( save_scorefiles ) {
//       filename += ".spscores";
//       g_pstrScoreFileName = filename.c_str();
//       WriteScoreFile( msa2_Cbw_Ubw2 );
//       cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
//       g_pstrScoreFileName = NULL;
//     } // End if save_scorefiles
//     
//     cout << "Training yet again using Unconditional Baum-Welch for 2 more iterations" << endl;
//     trainer_Cbw2.m_parameters.maxIterations = 2;
//     trainer_Cbw2.m_parameters.useUnconditionalBaumWelch = true;
//     score = trainer_Cbw2.train();
//     cout << "Now (after training again using Ubw for a total of 4 iterations max), the score is " << score << endl;//", and the profile is:" << endl;
//     //cout << *trainer_Cbw2.m_profile << endl;
//     // TODO: REMOVE
//     cout << "Calculating viterbi scores for the 'Cbw_Ubw4' profile2." << endl;
//     viterbi_profile2.copyPositions( *trainer_Cbw2.m_profile );
//     score =
//       viterbi_dp.forward_score_viterbi(
//         trainer_Cbw2.m_parameters,
//         viterbi_profile2, // *trainer_Cbw2.m_profile,
//         trainer_Cbw2.m_poly_sequences,
//         trainer_Cbw2.m_sequence_count,
//         viterbi_matrices
//       );
//     cout << "The total score for all sequences, using viterbi, is: " << score << endl;
//     ViterbiDPType::MultipleAlignment<ProfileType> ma2_Cbw_Ubw4 =
//       viterbi_dp.forward_viterbiAlign(
//         trainer_Cbw2.m_parameters,
//         viterbi_profile2, // *trainer_Cbw2.m_profile,
//         trainer_Cbw2.m_poly_sequences,
//         trainer_Cbw2.m_sequence_count,
//         viterbi_matrices // TODO: ARE THESE MATRICES SIZED CORRECTLY?
//       );
//     filename = filename_base + ".2.Cbw_Ubw4.blast";
//     outfile.open( filename.c_str() );
//     //cout << "Multiple Alignment is:" << endl;
//     //ma2_Cbw_Ubw4.toPairwiseStream( cout, &fasta2.m_descriptions );
//     ma2_Cbw_Ubw4.toPairwiseStream( outfile, &fasta2.m_descriptions );
//     cout << "Wrote out Multiple Alignment to " << filename << endl;
//     outfile.close();
//     
//     // Now try making an MSA out of it.
//     //MSA msa2_Cbw_Ubw4;
//     ma2_Cbw_Ubw4.appendToMSA( &msa2_Cbw_Ubw4, &fasta2.m_descriptions );
//     filename = filename_base + ".2.Cbw_Ubw4.blast.a2m";
//     TextFile msa2_Cbw_Ubw4_outfile =
//       TextFile( filename.c_str(), true );
//     msa2_Cbw_Ubw4.ToFile( msa2_Cbw_Ubw4_outfile );
//     msa2_Cbw_Ubw4_outfile.Close();
//     cout << "Wrote MSA to file " << filename << endl;
//     
//     // Trying score-getting
//     cout << "Total SP score is " << getNonInsertionSPScore( msa2_Cbw_Ubw4 ) << endl;
//     if( save_scorefiles ) {
//       filename += ".spscores";
//       g_pstrScoreFileName = filename.c_str();
//       WriteScoreFile( msa2_Cbw_Ubw4 );
//       cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
//       g_pstrScoreFileName = NULL;
//     } // End if save_scorefiles
//
//     cout << "Training yet again using Unconditional Baum-Welch until convergence." << endl;
//     trainer_Cbw2.m_parameters.maxIterations = 100;
//     trainer_Cbw2.m_parameters.scorePercentChangeMinimum_position_cycle = 100;
//     trainer_Cbw2.m_parameters.scorePercentChangeMinimum_iteration = 100;
//     trainer_Cbw2.m_parameters.useUnconditionalBaumWelch = true;
//     score = trainer_Cbw2.train();
//     cout << "Now (after training again using Ubw until convergence), the score is " << score << endl;//", and the profile is:" << endl;
//     //cout << *trainer_Cbw2.m_profile << endl;
//     // TODO: REMOVE
//     cout << "Calculating viterbi scores for the 'Cbw_Ubw' profile2." << endl;
//     viterbi_profile2.copyPositions( *trainer_Cbw2.m_profile );
//     score =
//       viterbi_dp.forward_score_viterbi(
//         trainer_Cbw2.m_parameters,
//         viterbi_profile2, // *trainer_Cbw2.m_profile,
//         trainer_Cbw2.m_poly_sequences,
//         trainer_Cbw2.m_sequence_count,
//         viterbi_matrices
//       );
//     cout << "The total score for all sequences, using viterbi, is: " << score << endl;
//     ViterbiDPType::MultipleAlignment<ProfileType> ma2_Cbw_Ubw =
//       viterbi_dp.forward_viterbiAlign(
//         trainer_Cbw2.m_parameters,
//         viterbi_profile2, // *trainer_Cbw2.m_profile,
//         trainer_Cbw2.m_poly_sequences,
//         trainer_Cbw2.m_sequence_count,
//         viterbi_matrices // TODO: ARE THESE MATRICES SIZED CORRECTLY?
//       );
//     filename = filename_base + ".2.Cbw_Ubw.blast";
//     outfile.open( filename.c_str() );
//     //cout << "Multiple Alignment is:" << endl;
//     //ma2_Cbw_Ubw.toPairwiseStream( cout, &fasta2.m_descriptions );
//     ma2_Cbw_Ubw.toPairwiseStream( outfile, &fasta2.m_descriptions );
//     cout << "Wrote out Multiple Alignment to " << filename << endl;
//     outfile.close();
//     
//     // Now try making an MSA out of it.
//     //MSA msa2_Cbw_Ubw;
//     ma2_Cbw_Ubw.appendToMSA( &msa2_Cbw_Ubw, &fasta2.m_descriptions );
//     filename = filename_base + ".2.Cbw_Ubw.blast.a2m";
//     TextFile msa2_Cbw_Ubw_outfile =
//       TextFile( filename.c_str(), true );
//     msa2_Cbw_Ubw.ToFile( msa2_Cbw_Ubw_outfile );
//     msa2_Cbw_Ubw_outfile.Close();
//     cout << "Wrote MSA to file " << filename << endl;
//     
//     // Trying score-getting
//     cout << "Total SP score is " << getNonInsertionSPScore( msa2_Cbw_Ubw ) << endl;
//     if( save_scorefiles ) {
//       filename += ".spscores";
//       g_pstrScoreFileName = filename.c_str();
//       WriteScoreFile( msa2_Cbw_Ubw );
//       cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
//       g_pstrScoreFileName = NULL;
//     } // End if save_scorefiles
//     
//   } // End if do_Cbw_Ubw

//   if( do_Ubw ) {
//     cout << "Training again (starting over), using Unconditional Baum-Welch until convergence." << endl;
//
//     profile2.copyFrom( notraining_profile2 );
//
//     ProfileTrainer<ProfileType, ScoreDifferenceType, ScoreType, MatrixValueType> trainer2_Ubw =
//       ProfileTrainer<ProfileType, ScoreDifferenceType, ScoreType, MatrixValueType>( &profile2, fasta2, num_sequences_to_use2 );
//     // TODO: REMOVE:
//     //cout << "fasta2.size() is " << fasta2.size() << endl;
//     //cout << "m_sequence_count is " << trainer2_Ubw.m_sequence_count << endl;
//     // TODO: REMOVE:
//     //trainer2_Ubw.m_parameters.profileValueMinimum = 1e-2; // ?
//     trainer2_Ubw.m_parameters.scorePercentChangeMinimum_position_cycle = 100;
//     trainer2_Ubw.m_parameters.scorePercentChangeMinimum_iteration = 100;
//     //trainer2_Ubw.m_parameters.useRabinerScaling = false;
//     //trainer2_Ubw.m_parameters.minIterations = 4;
//     //trainer2_Ubw.m_parameters.maxIterations = 2;
//     //trainer2_Ubw.m_parameters.maxPositionCycles = 2;
//     //trainer2_Ubw.m_parameters.maxPositionCycles_globals = 2;
//     //trainer2_Ubw.m_parameters.trainGlobalsFirst = true;
//     //trainer2_Ubw.m_parameters.useRabinerScaling = false;
//     trainer2_Ubw.m_parameters.usePriors = use_priors;
//     trainer2_Ubw.m_parameters.trainProfileGlobals = train_globals;
//     //trainer2_Ubw.m_parameters.maxPositionCycles = 3;
//     trainer2_Ubw.m_parameters.useUnconditionalBaumWelch = true;
//     //trainer2_Ubw.m_parameters.debug = DEBUG_All;
//     //trainer2_Ubw.m_parameters.verbosity = VERBOSITY_All;
//     trainer2_Ubw.m_parameters.verbosity = VERBOSITY_Low;
//     
//     //cout << "The profile (before) is:" << endl;
//     //cout << *trainer2_Ubw.m_profile << endl;
//     
//     score = trainer2_Ubw.train();
//     cout << "Now (after training), the score is " << score << endl; //", and the profile is:" << endl;
//     //cout << *trainer2_Ubw.m_profile << endl;
//     //ViterbiDPType::Matrix::SequentialAccessContainer viterbi_matrices =
//     //  ViterbiDPType::Matrix::SequentialAccessContainer(
//     //    *trainer2_Ubw.m_profile,
//     //    trainer2_Ubw.m_poly_sequences,
//     //    trainer2_Ubw.m_sequence_count
//     //  );
//     // TODO: REMOVE
//     cout << "Calculating viterbi scores for the 'Ubw' profile2." << endl;
//     viterbi_profile2.copyPositions( *trainer2_Ubw.m_profile );
//     score =
//       viterbi_dp.forward_score_viterbi(
//         trainer2_Ubw.m_parameters,
//         viterbi_profile2, // *trainer2_Ubw.m_profile,
//         trainer2_Ubw.m_poly_sequences,
//         trainer2_Ubw.m_sequence_count,
//         viterbi_matrices
//       );
//     cout << "The total score for all sequences, using viterbi, is: " << score << endl;
//     ViterbiDPType::MultipleAlignment<ProfileType> ma2_Ubw =
//       viterbi_dp.forward_viterbiAlign(
//         trainer2_Ubw.m_parameters,
//         viterbi_profile2, // *trainer2_Ubw.m_profile,
//         trainer2_Ubw.m_poly_sequences,
//         trainer2_Ubw.m_sequence_count,
//         viterbi_matrices // TODO: ARE THESE MATRICES SIZED CORRECTLY?
//       );
//     filename = filename_base + ".2.Ubw.blast";
//     outfile.open( filename.c_str() );
//     //cout << "Multiple Alignment is:" << endl;
//     //ma2_Ubw.toPairwiseStream( cout, &fasta2.m_descriptions );
//     ma2_Ubw.toPairwiseStream( outfile, &fasta2.m_descriptions );
//     cout << "Wrote out Multiple Alignment to " << filename << endl;
//     outfile.close();
//     
//     // Now try making an MSA out of it.
//     //MSA msa2_Ubw;
//     ma2_Ubw.appendToMSA( &msa2_Ubw, &fasta2.m_descriptions );
//     filename = filename_base + ".2.Ubw.blast.a2m";
//     TextFile msa2_Ubw_outfile =
//       TextFile( filename.c_str(), true );
//     msa2_Ubw.ToFile( msa2_Ubw_outfile );
//     msa2_Ubw_outfile.Close();
//     cout << "Wrote MSA to file " << filename << endl;
//     
//     // Trying score-getting
//     cout << "Total SP score is " << getNonInsertionSPScore( msa2_Ubw ) << endl;
//     if( save_scorefiles ) {
//       filename += ".spscores";
//       g_pstrScoreFileName = filename.c_str();
//       WriteScoreFile( msa2_Ubw );
//       cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
//       g_pstrScoreFileName = NULL;
//     } // End if save_scorefiles
//   } // End if do_Ubw
// endmark


   // Combine two separately-trained muscle files.
   cout << "Combining the two separately trained 'notraining' MSAs." << endl;
   ::MSA msa_combined_notraining;
   ToUpper( msa_notraining1 ); // Make it uppercase, for qscoring
   ToUpper( msa_notraining2 ); // Make it uppercase, for qscoring
   // Set the ids.
   for (unsigned uSeqIndex = 0; uSeqIndex < uCount1; ++uSeqIndex)
     msa_notraining1.SetSeqId(uSeqIndex, Ids1[ uSeqIndex ] );
   for (unsigned uSeqIndex = 0; uSeqIndex < uCount2; ++uSeqIndex)
     msa_notraining2.SetSeqId(uSeqIndex, Ids2[ uSeqIndex ] );
   PWPath path_after;
   AlignTwoMSAs( msa_notraining1, msa_notraining2, msa_combined_notraining, path_after );
   filename = filename_base + ".combined_notraining.muscle";    
   TextFile msa_combined_notraining_outfile =
     TextFile( filename.c_str(), true );
   msa_combined_notraining.ToFile( msa_combined_notraining_outfile );
   msa_combined_notraining_outfile.Close();
   cout << "Wrote out combined_notraining MSA to " << msa_combined_notraining_outfile.GetFileName() << endl;

   // Make sure to apply the sequence weights.
   // TODO.  To do this, we need first to run setup, including SetMuscleTree...?
   //SetMSAWeightsMuscle( msa_combined_notraining );

   cout << "Total SP score is " << getNonInsertionSPScore( msa_combined_notraining ) << endl;
   
   if( save_scorefiles ) {
     filename += ".spscores";
     g_pstrScoreFileName = filename.c_str();
     WriteScoreFile( msa_combined_notraining );
     cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
     g_pstrScoreFileName = NULL;
   } // End if save_scorefiles
   
   // Qscore
   if( have_ref ) {
     dSP = dPS = dCS = dInsane;
     //cout << "About to CompareMSA" << endl;
     CompareMSA( msa_combined_notraining, msa_ref, &dSP, &dPS, &dCS);
     //cout << "About to TC" << endl;
     dTC = TC( msa_combined_notraining, msa_ref );
     //cout << "Just did TC" << endl;
     cout << "Combined_notraining:" << endl;
     printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
            dSP, dPS, dCS, dTC);
     cout << "Muscle alone:" << endl;
     dSP = dPS = dCS = dInsane;
     CompareMSA( msa_muscle, msa_ref, &dSP, &dPS, &dCS);
     dTC = TC( msa_muscle, msa_ref );
     printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
            dSP, dPS, dCS, dTC);
   } // End if have_ref

   cout << "Combining the two separately trained 'Cbw' MSAs." << endl;
   ::MSA msa_combined_Cbw;
   ToUpper( msa_Cbw1 ); // Make it uppercase, for qscoring
   ToUpper( msa_Cbw2 ); // Make it uppercase, for qscoring
   // Set the ids.
   for (unsigned uSeqIndex = 0; uSeqIndex < uCount1; ++uSeqIndex)
     msa_Cbw1.SetSeqId(uSeqIndex, Ids1[ uSeqIndex ] );
   for (unsigned uSeqIndex = 0; uSeqIndex < uCount2; ++uSeqIndex)
     msa_Cbw2.SetSeqId(uSeqIndex, Ids2[ uSeqIndex ] );
   //PWPath path_after;
   AlignTwoMSAs( msa_Cbw1, msa_Cbw2, msa_combined_Cbw, path_after );
   filename = filename_base + ".combined_Cbw.muscle";
   TextFile msa_combined_Cbw_outfile =
     TextFile( filename.c_str(), true );
   msa_combined_Cbw.ToFile( msa_combined_Cbw_outfile );
   msa_combined_Cbw_outfile.Close();
   cout << "Wrote out combined_Cbw MSA to " << msa_combined_Cbw_outfile.GetFileName() << endl;

   // Make sure to apply the sequence weights.
   // TODO.  To do this, we need first to run setup, including SetMuscleTree...?
   //SetMSAWeightsMuscle( msa_combined_Cbw );

   cout << "Total SP score is " << getNonInsertionSPScore( msa_combined_Cbw ) << endl;
   
   if( save_scorefiles ) {
     filename += ".spscores";
     g_pstrScoreFileName = filename.c_str();
     WriteScoreFile( msa_combined_Cbw );
     cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
     g_pstrScoreFileName = NULL;
   } // End if save_scorefiles
   
   // Qscore
   if( have_ref ) {
     dSP = dPS = dCS = dInsane;
     //cout << "About to CompareMSA" << endl;
     CompareMSA( msa_combined_Cbw, msa_ref, &dSP, &dPS, &dCS);
     //cout << "About to TC" << endl;
     dTC = TC( msa_combined_Cbw, msa_ref );
     //cout << "Just did TC" << endl;
     cout << "Combined_Cbw:" << endl;
     printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
            dSP, dPS, dCS, dTC);
     cout << "Muscle alone:" << endl;
     dSP = dPS = dCS = dInsane;
     CompareMSA( msa_muscle, msa_ref, &dSP, &dPS, &dCS);
     dTC = TC( msa_muscle, msa_ref );
     printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
            dSP, dPS, dCS, dTC);
   } // End if have_ref

//   if( do_Cbw_Ubw ) {
//     cout << "Combining the two separately trained 'Cbw_Ubw2' MSAs." << endl;
//     ::MSA msa_combined_Cbw_Ubw2;
//     ToUpper( msa1_Cbw_Ubw2 ); // Make it uppercase, for qscoring
//     ToUpper( msa2_Cbw_Ubw2 ); // Make it uppercase, for qscoring
//     // Set the ids.
//     for (unsigned uSeqIndex = 0; uSeqIndex < uCount1; ++uSeqIndex)
//       msa1_Cbw_Ubw2.SetSeqId(uSeqIndex, Ids1[ uSeqIndex ] );
//     for (unsigned uSeqIndex = 0; uSeqIndex < uCount2; ++uSeqIndex)
//       msa2_Cbw_Ubw2.SetSeqId(uSeqIndex, Ids2[ uSeqIndex ] );
//     //PWPath path_after;
//     AlignTwoMSAs( msa1_Cbw_Ubw2, msa2_Cbw_Ubw2, msa_combined_Cbw_Ubw2, path_after );
//     filename = filename_base + ".combined_Cbw_Ubw2.muscle";    
//     TextFile msa_combined_Cbw_Ubw2_outfile =
//       TextFile( filename.c_str(), true );
//     msa_combined_Cbw_Ubw2.ToFile( msa_combined_Cbw_Ubw2_outfile );
//     msa_combined_Cbw_Ubw2_outfile.Close();
//     cout << "Wrote out combined_Cbw_Ubw2 MSA to " << msa_combined_Cbw_Ubw2_outfile.GetFileName() << endl;
//     
//     // Make sure to apply the sequence weights.
//     // TODO.  To do this, we need first to run setup, including SetMuscleTree...?
//     //SetMSAWeightsMuscle( msa_combined_Cbw_Ubw2 );
//     
//     cout << "Total SP score is " << getNonInsertionSPScore( msa_combined_Cbw_Ubw2 ) << endl;
//     
//     if( save_scorefiles ) {
//       filename += ".spscores";
//       g_pstrScoreFileName = filename.c_str();
//       WriteScoreFile( msa_combined_Cbw_Ubw2 );
//       cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
//       g_pstrScoreFileName = NULL;
//     } // End if save_scorefiles
//     
//     // Qscore
//     if( have_ref ) {
//       dSP = dPS = dCS = dInsane;
//       //cout << "About to CompareMSA" << endl;
//       CompareMSA( msa_combined_Cbw_Ubw2, msa_ref, &dSP, &dPS, &dCS);
//       //cout << "About to TC" << endl;
//       dTC = TC( msa_combined_Cbw_Ubw2, msa_ref );
//       //cout << "Just did TC" << endl;
//       cout << "Combined_Cbw_Ubw2:" << endl;
//       printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
//              dSP, dPS, dCS, dTC);
//       cout << "Muscle alone:" << endl;
//       dSP = dPS = dCS = dInsane;
//       CompareMSA( msa_muscle, msa_ref, &dSP, &dPS, &dCS);
//       dTC = TC( msa_muscle, msa_ref );
//       printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
//              dSP, dPS, dCS, dTC);
//     } // End if have_ref
//
//     cout << "Combining the two separately trained 'Cbw_Ubw4' MSAs." << endl;
//     ::MSA msa_combined_Cbw_Ubw4;
//     ToUpper( msa1_Cbw_Ubw4 ); // Make it uppercase, for qscoring
//     ToUpper( msa2_Cbw_Ubw4 ); // Make it uppercase, for qscoring
//     // Set the ids.
//     for (unsigned uSeqIndex = 0; uSeqIndex < uCount1; ++uSeqIndex)
//       msa1_Cbw_Ubw4.SetSeqId(uSeqIndex, Ids1[ uSeqIndex ] );
//     for (unsigned uSeqIndex = 0; uSeqIndex < uCount2; ++uSeqIndex)
//       msa2_Cbw_Ubw4.SetSeqId(uSeqIndex, Ids2[ uSeqIndex ] );
//     //PWPath path_after;
//     AlignTwoMSAs( msa1_Cbw_Ubw4, msa2_Cbw_Ubw4, msa_combined_Cbw_Ubw4, path_after );
//     filename = filename_base + ".combined_Cbw_Ubw4.muscle";    
//     TextFile msa_combined_Cbw_Ubw4_outfile =
//       TextFile( filename.c_str(), true );
//     msa_combined_Cbw_Ubw4.ToFile( msa_combined_Cbw_Ubw4_outfile );
//     msa_combined_Cbw_Ubw4_outfile.Close();
//     cout << "Wrote out combined_Cbw_Ubw4 MSA to " << msa_combined_Cbw_Ubw4_outfile.GetFileName() << endl;
//     
//     // Make sure to apply the sequence weights.
//     // TODO.  To do this, we need first to run setup, including SetMuscleTree...?
//     //SetMSAWeightsMuscle( msa_combined_Cbw_Ubw4 );
//     
//     cout << "Total SP score is " << getNonInsertionSPScore( msa_combined_Cbw_Ubw4 ) << endl;
//     
//     if( save_scorefiles ) {
//       filename += ".spscores";
//       g_pstrScoreFileName = filename.c_str();
//       WriteScoreFile( msa_combined_Cbw_Ubw4 );
//       cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
//       g_pstrScoreFileName = NULL;
//     } // End if save_scorefiles
//     
//     // Qscore
//     if( have_ref ) {
//       dSP = dPS = dCS = dInsane;
//       //cout << "About to CompareMSA" << endl;
//       CompareMSA( msa_combined_Cbw_Ubw4, msa_ref, &dSP, &dPS, &dCS);
//       //cout << "About to TC" << endl;
//       dTC = TC( msa_combined_Cbw_Ubw4, msa_ref );
//       //cout << "Just did TC" << endl;
//       cout << "Combined_Cbw_Ubw4:" << endl;
//       printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
//              dSP, dPS, dCS, dTC);
//       cout << "Muscle alone:" << endl;
//       dSP = dPS = dCS = dInsane;
//       CompareMSA( msa_muscle, msa_ref, &dSP, &dPS, &dCS);
//       dTC = TC( msa_muscle, msa_ref );
//       printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
//              dSP, dPS, dCS, dTC);
//     } // End if have_ref
//
//     cout << "Combining the two separately trained 'Cbw_Ubw' MSAs." << endl;
//     ::MSA msa_combined_Cbw_Ubw;
//     ToUpper( msa1_Cbw_Ubw ); // Make it uppercase, for qscoring
//     ToUpper( msa2_Cbw_Ubw ); // Make it uppercase, for qscoring
//     // Set the ids.
//     for (unsigned uSeqIndex = 0; uSeqIndex < uCount1; ++uSeqIndex)
//       msa1_Cbw_Ubw.SetSeqId(uSeqIndex, Ids1[ uSeqIndex ] );
//     for (unsigned uSeqIndex = 0; uSeqIndex < uCount2; ++uSeqIndex)
//       msa2_Cbw_Ubw.SetSeqId(uSeqIndex, Ids2[ uSeqIndex ] );
//     //PWPath path_after;
//     AlignTwoMSAs( msa1_Cbw_Ubw, msa2_Cbw_Ubw, msa_combined_Cbw_Ubw, path_after );
//     filename = filename_base + ".combined_Cbw_Ubw.muscle";    
//     TextFile msa_combined_Cbw_Ubw_outfile =
//       TextFile( filename.c_str(), true );
//     msa_combined_Cbw_Ubw.ToFile( msa_combined_Cbw_Ubw_outfile );
//     msa_combined_Cbw_Ubw_outfile.Close();
//     cout << "Wrote out combined_Cbw_Ubw MSA to " << msa_combined_Cbw_Ubw_outfile.GetFileName() << endl;
//     
//     // Make sure to apply the sequence weights.
//     // TODO.  To do this, we need first to run setup, including SetMuscleTree...?
//     //SetMSAWeightsMuscle( msa_combined_Cbw_Ubw );
//     
//     cout << "Total SP score is " << getNonInsertionSPScore( msa_combined_Cbw_Ubw ) << endl;
//     
//     if( save_scorefiles ) {
//       filename += ".spscores";
//       g_pstrScoreFileName = filename.c_str();
//       WriteScoreFile( msa_combined_Cbw_Ubw );
//       cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
//       g_pstrScoreFileName = NULL;
//     } // End if save_scorefiles
//     
//     // Qscore
//     if( have_ref ) {
//       dSP = dPS = dCS = dInsane;
//       //cout << "About to CompareMSA" << endl;
//       CompareMSA( msa_combined_Cbw_Ubw, msa_ref, &dSP, &dPS, &dCS);
//       //cout << "About to TC" << endl;
//       dTC = TC( msa_combined_Cbw_Ubw, msa_ref );
//       //cout << "Just did TC" << endl;
//       cout << "Combined_Cbw_Ubw:" << endl;
//       printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
//              dSP, dPS, dCS, dTC);
//       cout << "Muscle alone:" << endl;
//       dSP = dPS = dCS = dInsane;
//       CompareMSA( msa_muscle, msa_ref, &dSP, &dPS, &dCS);
//       dTC = TC( msa_muscle, msa_ref );
//       printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
//              dSP, dPS, dCS, dTC);
//     } // End if have_ref
//
//   } // End if do_Cbw_Ubw
//
//   if( do_Ubw ) {
//     cout << "Combining the two separately trained 'Ubw' MSAs." << endl;
//     ::MSA msa_combined_Ubw;
//     ToUpper( msa1_Ubw ); // Make it uppercase, for qscoring
//     ToUpper( msa2_Ubw ); // Make it uppercase, for qscoring
//     // Set the ids.
//     for (unsigned uSeqIndex = 0; uSeqIndex < uCount1; ++uSeqIndex)
//       msa1_Ubw.SetSeqId(uSeqIndex, Ids1[ uSeqIndex ] );
//     for (unsigned uSeqIndex = 0; uSeqIndex < uCount2; ++uSeqIndex)
//       msa2_Ubw.SetSeqId(uSeqIndex, Ids2[ uSeqIndex ] );
//     //PWPath path_after;
//     AlignTwoMSAs( msa1_Ubw, msa2_Ubw, msa_combined_Ubw, path_after );
//     filename = filename_base + ".combined_Ubw.muscle";    
//     TextFile msa_combined_Ubw_outfile =
//       TextFile( filename.c_str(), true );
//     msa_combined_Ubw.ToFile( msa_combined_Ubw_outfile );
//     msa_combined_Ubw_outfile.Close();
//     cout << "Wrote out combined_Ubw MSA to " << msa_combined_Ubw_outfile.GetFileName() << endl;
//     
//     // Make sure to apply the sequence weights.
//     // TODO.  To do this, we need first to run setup, including SetMuscleTree...?
//     //SetMSAWeightsMuscle( msa_combined_Ubw );
//     
//     cout << "Total SP score is " << getNonInsertionSPScore( msa_combined_Ubw ) << endl;
//     
//     if( save_scorefiles ) {
//       filename += ".spscores";
//       g_pstrScoreFileName = filename.c_str();
//       WriteScoreFile( msa_combined_Ubw );
//       cout << "Wrote out SP scores to " << g_pstrScoreFileName << endl;
//       g_pstrScoreFileName = NULL;
//     } // End if save_scorefiles
//     
//     // Qscore
//     if( have_ref ) {
//       dSP = dPS = dCS = dInsane;
//       //cout << "About to CompareMSA" << endl;
//       CompareMSA( msa_combined_Ubw, msa_ref, &dSP, &dPS, &dCS);
//       //cout << "About to TC" << endl;
//       dTC = TC( msa_combined_Ubw, msa_ref );
//       //cout << "Just did TC" << endl;
//       cout << "Combined_Ubw:" << endl;
//       printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
//              dSP, dPS, dCS, dTC);
//       cout << "Muscle alone:" << endl;
//       dSP = dPS = dCS = dInsane;
//       CompareMSA( msa_muscle, msa_ref, &dSP, &dPS, &dCS);
//       dTC = TC( msa_muscle, msa_ref );
//       printf("Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
//              dSP, dPS, dCS, dTC);
//     } // End if have_ref
//   } // End if do_Ubw
   
   delete[] Ids1;
   delete[] Ids2;
 } // End if do_split

  // ENDMARK

  return 0;
} // prolific( int, char ** )

} // End namespace galosh

int
main ( int argc, char **argv )
{
  return galosh::prolific( argc, argv );
} // main( int, char** )

    /*
    hmmer::SetAlphabet( hmmAMINO ); // First we must call this to finish setting up.
    
    hmmer::HMMFILE        *hmmfp;
    //char           *hmmfile;
    filename = ( fasta_in_str + muscle_suffix_str + hmm_suffix_str );
    char * hmmfile = // hmmer FileOpen does not have a const modifier for the filename argument, so we need to do a const_cast here.  TODO: Fix hmmer code?
      const_cast<char *>( filename.c_str() );
    //char hmmfile[] =
      //"prefab4-in/1atiA_1g5hA.muscle.hmmbuild";
      //"prefab4-in/1g5hA.muscle.hmmbuild";
      //"prefab4-in/1atiA.from_1atiA_1g5hA.muscle.hmmbuild";
    //struct hmmer::plan7_s *hmm;
    char            env[] = "HMMERDB";  // (a la BLASTDB) 
                                        // 
    cout << "Reading HMMer hmm from \"" << hmmfile << "\"" << endl;
    hmmfp = hmmer::HMMFileOpen(hmmfile, env);   // NULL on failure
    */
    //while (hmmer::HMMFileRead(hmmfp, &hmm)) {    // 0 if no more HMMs
    //}
    //hmmer::HMMFileClose(hmmfp);

    //MSA msa1;
    //TextFile msa1_infile =
    //  TextFile( "prefab4-in/1atiA.from_1atiA_1g5hA.pure.muscle.hmmbuild.profuse.TESTING.Cbw.blast.a2m" );
    //msa1.FromFile( msa1_infile );
    //cout << "READ in MSA from " << msa1_infile.GetFileName() << endl;
    //MSA msa2;
    //TextFile msa2_infile =
    //  TextFile( "prefab4-in/1g5hA.pure.muscle.hmmbuild.profuse.TESTING.Cbw.blast.a2m" );
    //msa2.FromFile( msa2_infile );
    //cout << "READ in MSA from " << msa2_infile.GetFileName() << endl;
    
    // Now make a profile from the msa
    //msa1.FixAlpha(); // Is this necessary?
    //msa2.FixAlpha(); // Is this necessary?
    //SetPPScore(); // Or this?
    //ProfPos * msa1_profile = ProfileFromMSA( msa );
    //cout << "Converted msa to muscle-profile." << endl;
    //// ? Trying this
    //ListProfile( msa1_profile, msa.GetColCount(), &msa );
    
    // Ok, align the two branches:
    // ProfileProfile( msa1, msa2, msa_combined );
    // see also Profile(), which opens the files and writes to a file...
    
    //MSA::SetIdCount( msa1.GetSeqCount() + msa2.GetSeqCount() ); // Is this necessary?
#endif // __HAVE_MUSCLE
