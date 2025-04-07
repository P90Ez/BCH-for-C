//based on https://www.eccpage.com/bch3.c
//If you value you sanity, don't scroll further. I tried to optimize the memory allocations as good as possible, but did not dare to touch any of the algorithms.
//Some of the comments do not make sense anymore, as I tried to rename the variables to something more useful than "bb" or "jj". Please refer to the original implementation for to make use of the comments.

#include "BCH.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define MAX_GaloisFieldOrder 20

static uint8_t galoisFieldOrder = 0;
static uint32_t n = 0;
static uint32_t k = 0;
static uint8_t t = 0;

//GF look up tables
static int* logTable = 0; 
static int* antiLogTable = 0;

//generator polynomial
static int* generatorPolynomial = 0;

static bool BCH_GenerateGaloisField();
static void BCH_GenerateGeneratorPolynomial();

uint8_t GetBitAtIndex(uint8_t const*const Buffer, uint32_t const Index)
{
	if (Buffer == 0) return 0;

	uint32_t const ArrayIndex = Index / 8;
	uint8_t const BitIndex = Index % 8;

	uint8_t byte = *(Buffer + ArrayIndex);
	return (byte >> BitIndex) & 0x01;
}

void SetBitAtIndex(uint8_t * const Buffer, uint32_t const Index, uint8_t const BitValue)
{
	if (Buffer == 0) return;

	uint32_t const ArrayIndex = Index / 8;
	uint8_t const BitIndex = Index % 8;

	*(Buffer + ArrayIndex) &= ~(0x01 << BitIndex);
	*(Buffer + ArrayIndex) |= (BitValue & 0x01) << BitIndex;
}


bool BCH_Init(uint8_t const m, uint8_t const _t)
{
	if (m < 2 || m > MAX_GaloisFieldOrder) return false;

	galoisFieldOrder = m;
	t = _t;
	n = pow(2, galoisFieldOrder) - 1;

	if (!BCH_GenerateGaloisField()) return false;
	BCH_GenerateGeneratorPolynomial();

	return true;
}

void BCH_Cleanup()
{
	free(logTable); logTable = 0;
	free(antiLogTable); antiLogTable = 0;
	free(generatorPolynomial); generatorPolynomial = 0;
}

int BCH_GetCodeLength() { return n; }
int BCH_GetWordLength() { return k; }


static bool BCH_GenerateGaloisField()
/*
 * Generate field GF(2**m) from the irreducible polynomial p(X) with
 * coefficients in p[0]..p[m].
 *
 * Lookup tables:
 *   index->polynomial form: alpha_to[] contains j=alpha^i;
 *   polynomial form -> index form:	index_of[j=alpha^i] = i
 *
 * alpha=2 is the primitive element of GF(2**m)
 */
{
	int i, mask;
	logTable = calloc(n, sizeof(int));
	antiLogTable = calloc(n, sizeof(int));

	if (logTable == 0 || antiLogTable == 0) return false;

	uint8_t primitivePoly[MAX_GaloisFieldOrder];
	memset(primitivePoly, 0, MAX_GaloisFieldOrder);

	switch (galoisFieldOrder)
	{
		case 2: primitivePoly[1] = 1; break;
		case 3: primitivePoly[1] = 1; break;
		case 4: primitivePoly[1] = 1; break;
		case 5: primitivePoly[2] = 1; break;
		case 6: primitivePoly[1] = 1; break;
		case 7: primitivePoly[1] = 1; break;
		case 8: primitivePoly[4] = primitivePoly[5] = primitivePoly[6] = 1; break;
		case 9: primitivePoly[4] = 1; break;
		case 10: primitivePoly[3] = 1; break;
		case 11: primitivePoly[2] = 1; break;
		case 12: primitivePoly[3] = primitivePoly[4] = primitivePoly[7] = 1; break;
		case 13: primitivePoly[1] = primitivePoly[3] = primitivePoly[4] = 1; break;
		case 14: primitivePoly[1] = primitivePoly[11] = primitivePoly[12] = 1; break;
		case 15: primitivePoly[1] = 1; break;
		case 16: primitivePoly[2] = primitivePoly[3] = primitivePoly[5] = 1; break;
		case 17: primitivePoly[3] = 1; break;
		case 18: primitivePoly[7] = 1; break;
		case 19: primitivePoly[1] = primitivePoly[5] = primitivePoly[6] = 1; break;
		case 20: primitivePoly[3] = 1; break;
	}

	mask = 1;
	logTable[galoisFieldOrder] = 0;
	for (i = 0; i < galoisFieldOrder; i++) {
		logTable[i] = mask;
		antiLogTable[logTable[i]] = i;
		if (primitivePoly[i] != 0)
			logTable[galoisFieldOrder] ^= mask;
		mask <<= 1;
	}
	antiLogTable[logTable[galoisFieldOrder]] = galoisFieldOrder;
	mask >>= 1;
	for (i = galoisFieldOrder + 1; i < n; i++) {
		if (logTable[i - 1] >= mask)
			logTable[i] = logTable[galoisFieldOrder] ^ ((logTable[i - 1] ^ mask) << 1);
		else
			logTable[i] = logTable[i - 1] << 1;
		antiLogTable[logTable[i]] = i;
	}
	antiLogTable[0] = -1;

	return true;
}


static void BCH_GenerateGeneratorPolynomial()
/*
 * Compute the generator polynomial of a binary BCH code. Fist generate the
 * cycle sets modulo 2**m - 1, cycle[][] =  (i, 2*i, 4*i, ..., 2^l*i). Then
 * determine those cycle sets that contain integers in the set of (d-1)
 * consecutive integers {1..(d-1)}. The generator polynomial is calculated
 * as the product of linear factors of the form (x+alpha^i), for every i in
 * the above cycle sets.
 */
{
	int	cycleIndex, cycleSetIndex, nextCycleRep, numberOfSelectedSets;
	int	found, nextValue, numberOfCycles, currentRoot, numberOfTerms, redundancyBits;

	int** cycleSet = malloc(sizeof(int*) * n);
	for (int i = 0; i < n; i++)
	{
		*(cycleSet + i) = calloc(MAX_GaloisFieldOrder, sizeof(int));
	}

	int* setSize = calloc(n, sizeof(int));

	/* Generate cycle sets modulo n, n = 2**m - 1 */
	cycleSet[0][0] = 0;
	setSize[0] = 1;
	cycleSet[1][0] = 1;
	setSize[1] = 1;
	cycleSetIndex = 1;			/* cycle set index */

	do {
		/* Generate the jj-th cycle set */
		cycleIndex = 0;
		do {
			cycleIndex++;
			cycleSet[cycleSetIndex][cycleIndex] = (cycleSet[cycleSetIndex][cycleIndex - 1] * 2) % n;
			setSize[cycleSetIndex]++;
			nextValue = (cycleSet[cycleSetIndex][cycleIndex] * 2) % n;
		} while (nextValue != cycleSet[cycleSetIndex][0]);
		/* Next cycle set representative */
		nextCycleRep = 0;
		do {
			nextCycleRep++;
			found = 0;
			for (cycleIndex = 1; ((cycleIndex <= cycleSetIndex) && (!found)); cycleIndex++)
				/* Examine previous cycle sets */
				for (numberOfSelectedSets = 0; ((numberOfSelectedSets < setSize[cycleIndex]) && (!found)); numberOfSelectedSets++)
					if (nextCycleRep == cycleSet[cycleIndex][numberOfSelectedSets])
						found = 1;
		} while ((found) && (nextCycleRep < (n - 1)));
		if (!(found)) {
			cycleSetIndex++;	/* next cycle set index */
			cycleSet[cycleSetIndex][0] = nextCycleRep;
			setSize[cycleSetIndex] = 1;
		}
	} while (nextCycleRep < (n - 1));
	numberOfCycles = cycleSetIndex;		/* number of cycle sets modulo n */

	int d = 2 * t + 1;

	int* selectedCycles = calloc(numberOfCycles, sizeof(int));
	int* rootElements = calloc(n, sizeof(int));

	/* Search for roots 1, 2, ..., d-1 in cycle sets */
	numberOfSelectedSets = 0;
	redundancyBits = 0;
	for (cycleIndex = 1; cycleIndex <= numberOfCycles; cycleIndex++) {
		selectedCycles[numberOfSelectedSets] = 0;
		found = 0;
		for (cycleSetIndex = 0; ((cycleSetIndex < setSize[cycleIndex]) && (!found)); cycleSetIndex++)
			for (currentRoot = 1; ((currentRoot < d) && (!found)); currentRoot++)
				if (currentRoot == cycleSet[cycleIndex][cycleSetIndex]) {
					found = 1;
					selectedCycles[numberOfSelectedSets] = cycleIndex;
				}
		if (selectedCycles[numberOfSelectedSets]) {
			redundancyBits += setSize[selectedCycles[numberOfSelectedSets]];
			numberOfSelectedSets++;
		}
	}
	numberOfTerms = numberOfSelectedSets;
	numberOfSelectedSets = 1;
	for (cycleIndex = 0; cycleIndex < numberOfTerms; cycleIndex++)
		for (cycleSetIndex = 0; cycleSetIndex < setSize[selectedCycles[cycleIndex]]; cycleSetIndex++) {
			rootElements[numberOfSelectedSets] = cycleSet[selectedCycles[cycleIndex]][cycleSetIndex];
			numberOfSelectedSets++;
		}

	k = n - redundancyBits;

	for (int i = 0; i < n; i++)
	{
		free(*(cycleSet + i));
	}
	free(cycleSet); cycleSet = 0;
	free(setSize); setSize = 0;
	free(selectedCycles); selectedCycles = 0;

	generatorPolynomial = calloc(redundancyBits + 1, sizeof(int));

	/* Compute the generator polynomial */
	generatorPolynomial[0] = logTable[rootElements[1]];
	generatorPolynomial[1] = 1;		/* g(x) = (X + zeros[1]) initially */
	for (cycleIndex = 2; cycleIndex <= redundancyBits; cycleIndex++) {
		generatorPolynomial[cycleIndex] = 1;
		for (cycleSetIndex = cycleIndex - 1; cycleSetIndex > 0; cycleSetIndex--)
			if (generatorPolynomial[cycleSetIndex] != 0)
				generatorPolynomial[cycleSetIndex] = generatorPolynomial[cycleSetIndex - 1] ^ logTable[(antiLogTable[generatorPolynomial[cycleSetIndex]] + rootElements[cycleIndex]) % n];
			else
				generatorPolynomial[cycleSetIndex] = generatorPolynomial[cycleSetIndex - 1];
		generatorPolynomial[0] = logTable[(antiLogTable[generatorPolynomial[0]] + rootElements[cycleIndex]) % n];
	}

	free(rootElements); rootElements = 0;
}

uint8_t* BCH_Encode(uint8_t const*const word)
{
	if (word == 0 || generatorPolynomial == 0 || n == 0 || k == 0 || t == 0) return 0;

	int    dataIndex, redundancyIndex;
	int    errorFeedback;

	uint32_t const NBytes = ceil(n / 8.0);
	uint8_t* codeWord = calloc(NBytes, sizeof(uint8_t));

	for (int i = 0; i < k; i++)
	{
		SetBitAtIndex(codeWord, n - k + i, GetBitAtIndex(word, i));
	}

	for (dataIndex = k - 1; dataIndex >= 0; dataIndex--) {
		//errorFeedback = word[dataIndex] ^ codeWord[n - k - 1];
		errorFeedback = GetBitAtIndex(word, dataIndex) ^ GetBitAtIndex(codeWord, n - k - 1);
		if (errorFeedback != 0) {
			for (redundancyIndex = n - k - 1; redundancyIndex > 0; redundancyIndex--)
				if (generatorPolynomial[redundancyIndex] != 0)
				{
					//codeWord[redundancyIndex] = codeWord[redundancyIndex - 1] ^ errorFeedback;
					SetBitAtIndex(codeWord, redundancyIndex, GetBitAtIndex(codeWord, redundancyIndex - 1) ^ errorFeedback);
				}
				else
				{
					//codeWord[redundancyIndex] = codeWord[redundancyIndex - 1];
					SetBitAtIndex(codeWord, redundancyIndex, GetBitAtIndex(codeWord, redundancyIndex - 1));
				}
			//codeWord[0] = generatorPolynomial[0] && errorFeedback;
			SetBitAtIndex(codeWord, 0, generatorPolynomial[0] && errorFeedback);
		}
		else {
			for (redundancyIndex = n - k - 1; redundancyIndex > 0; redundancyIndex--)
			{
				//codeWord[redundancyIndex] = codeWord[redundancyIndex - 1];
				SetBitAtIndex(codeWord, redundancyIndex, GetBitAtIndex(codeWord, redundancyIndex - 1));
			}
			//codeWord[0] = 0;
			SetBitAtIndex(codeWord, 0, 0);
		}
	}

	return codeWord;
}

uint8_t* BCH_Decode(uint8_t const*const codeword)
/*
 * Simon Rockliff's implementation of Berlekamp's algorithm.
 *
 * Assume we have received bits in recd[i], i=0..(n-1).
 *
 * Compute the 2*t syndromes by substituting alpha^i into rec(X) and
 * evaluating, storing the syndromes in s[i], i=1..2t (leave s[0] zero) .
 * Then we use the Berlekamp algorithm to find the error location polynomial
 * elp[i].
 *
 * If the degree of the elp is >t, then we cannot correct all the errors, and
 * we have detected an uncorrectable error pattern. We output the information
 * bits uncorrected.
 *
 * If the degree of elp is <=t, we substitute alpha^i , i=1..n into the elp
 * to get the roots, hence the inverse roots, the error location numbers.
 * This step is usually called "Chien's search".
 *
 * If the number of errors located is not equal the degree of the elp, then
 * the decoder assumes that there are more than t errors and cannot correct
 * them, only detect them. We output the information bits uncorrected.
 */
{
	if (codeword == 0 || antiLogTable == 0 || logTable == 0 || n == 0 || k == 0 || t == 0) return 0;

	int	syndromeIndex, bitIndex, elpStep, tempIndex, errorCount = 0, syndromeError = 0;
	int t2 = 2 * t;

	uint32_t const NBytes = ceil(n / 8.0);
	uint8_t* recd = malloc(sizeof(uint8_t) * NBytes);
	memcpy(recd, codeword, NBytes);

	int* syndromes = calloc(t2 + 1, sizeof(int));

	/* first form the syndromes */
	for (syndromeIndex = 1; syndromeIndex <= t2; syndromeIndex++) {
		int val = 0;
		for (bitIndex = 0; bitIndex < n; bitIndex++)
		{
			//if (recd[bitIndex] != 0)
			if (GetBitAtIndex(recd, bitIndex) != 0)
				val ^= logTable[(syndromeIndex * bitIndex) % n];
		}
		if (val != 0)
			syndromeError = 1; /* set error flag if non-zero syndrome */
		/*
		 * Note:    If the code is used only for ERROR DETECTION, then
		 *          exit program here indicating the presence of errors.
		 */
		 /* convert syndrome from polynomial form to index form  */

		if(val < n) syndromes[syndromeIndex] = antiLogTable[val];
	}

	if (syndromeError > 0) {	/* if there are errors, try to correct them */
		/*
		 * Compute the error location polynomial via the Berlekamp
		 * iterative algorithm. Following the terminology of Lin and
		 * Costello's book :   d[u] is the 'mu'th discrepancy, where
		 * u='mu'+1 and 'mu' (the Greek letter!) is the step number
		 * ranging from -1 to 2*t (see L&C),  l[u] is the degree of
		 * the elp at that step, and u_l[u] is the difference between
		 * the step number and the degree of the elp.
		 */
		 /* initialise table entries */

		int16_t** errorLocatorPolynomial = malloc(sizeof(int16_t*) * (t2+2));
		for (int i = 0; i < t2 + 2; i++)
		{
			errorLocatorPolynomial[i] = calloc(t2, sizeof(int16_t));
		}
		int* discrepancy = calloc(t2 + 1, sizeof(int));
		int* elpDegree = calloc(t2 + 2, sizeof(int));
		int* elpStepDiff = calloc(t2 + 2, sizeof(int));
		int* regPolynomial = calloc(t2, sizeof(int));

		discrepancy[0] = 0;			/* index form */
		discrepancy[1] = syndromes[1];		/* index form */
		errorLocatorPolynomial[0][0] = 0;		/* index form */
		errorLocatorPolynomial[1][0] = 1;		/* polynomial form */
		for (syndromeIndex = 1; syndromeIndex < t2; syndromeIndex++) {
			errorLocatorPolynomial[0][syndromeIndex] = -1;	/* index form */
			errorLocatorPolynomial[1][syndromeIndex] = 0;	/* polynomial form */
		}
		elpDegree[0] = 0;
		elpDegree[1] = 0;
		elpStepDiff[0] = -1;
		elpStepDiff[1] = 0;
		elpStep = 0;

		do {
			elpStep++;
			if (discrepancy[elpStep] == -1) {
				elpDegree[elpStep + 1] = elpDegree[elpStep];
				for (syndromeIndex = 0; syndromeIndex <= elpDegree[elpStep]; syndromeIndex++) {
					errorLocatorPolynomial[elpStep + 1][syndromeIndex] = errorLocatorPolynomial[elpStep][syndromeIndex];
					errorLocatorPolynomial[elpStep][syndromeIndex] = antiLogTable[errorLocatorPolynomial[elpStep][syndromeIndex]];
				}
			}
			else
				/*
				 * search for words with greatest u_lu[q] for
				 * which d[q]!=0
				 */
			{
				tempIndex = elpStep - 1;
				while ((discrepancy[tempIndex] == -1) && (tempIndex > 0))
					tempIndex--;
				/* have found first non-zero d[q]  */
				if (tempIndex > 0) {
					bitIndex = tempIndex;
					do {
						bitIndex--;
						if ((discrepancy[bitIndex] != -1) && (elpStepDiff[tempIndex] < elpStepDiff[bitIndex]))
							tempIndex = bitIndex;
					} while (bitIndex > 0);
				}

				/*
				 * have now found q such that d[u]!=0 and
				 * u_lu[q] is maximum
				 */
				 /* store degree of new elp polynomial */
				if (elpDegree[elpStep] > elpDegree[tempIndex] + elpStep - tempIndex)
					elpDegree[elpStep + 1] = elpDegree[elpStep];
				else
					elpDegree[elpStep + 1] = elpDegree[tempIndex] + elpStep - tempIndex;

				/* form new elp(x) */
				for (syndromeIndex = 0; syndromeIndex < t2; syndromeIndex++)
				{
					errorLocatorPolynomial[elpStep + 1][syndromeIndex] = 0;
				}
				for (syndromeIndex = 0; syndromeIndex <= elpDegree[tempIndex]; syndromeIndex++)
					if (errorLocatorPolynomial[tempIndex][syndromeIndex] != -1)
						errorLocatorPolynomial[elpStep + 1][syndromeIndex + elpStep - tempIndex] =
						logTable[(discrepancy[elpStep] + n - discrepancy[tempIndex] + errorLocatorPolynomial[tempIndex][syndromeIndex]) % n];
				for (syndromeIndex = 0; syndromeIndex <= elpDegree[elpStep]; syndromeIndex++) {
					errorLocatorPolynomial[elpStep + 1][syndromeIndex] ^= errorLocatorPolynomial[elpStep][syndromeIndex];
					errorLocatorPolynomial[elpStep][syndromeIndex] = antiLogTable[errorLocatorPolynomial[elpStep][syndromeIndex]];
				}
			}
			elpStepDiff[elpStep + 1] = elpStep - elpDegree[elpStep + 1];

			/* form (u+1)th discrepancy */
			if (elpStep < t2) {
				/* no discrepancy computed on last iteration */
				if (syndromes[elpStep + 1] != -1)
					discrepancy[elpStep + 1] = logTable[syndromes[elpStep + 1]];
				else
					discrepancy[elpStep + 1] = 0;
				for (syndromeIndex = 1; syndromeIndex <= elpDegree[elpStep + 1]; syndromeIndex++)
				{
					if ((syndromes[elpStep + 1 - syndromeIndex] != -1) && (errorLocatorPolynomial[elpStep + 1][syndromeIndex] != 0))
						discrepancy[elpStep + 1] ^= logTable[(syndromes[elpStep + 1 - syndromeIndex]
							+ antiLogTable[errorLocatorPolynomial[elpStep + 1][syndromeIndex]]) % n];
				}
				/* put d[u+1] into index form */
				discrepancy[elpStep + 1] = antiLogTable[discrepancy[elpStep + 1]];
			}
		} while ((elpStep < t2) && (elpDegree[elpStep + 1] <= t));

		elpStep++;
		if (elpDegree[elpStep] <= t) {/* Can correct errors */
			/* put elp into index form */
			for (syndromeIndex = 0; syndromeIndex <= elpDegree[elpStep]; syndromeIndex++)
			{
				errorLocatorPolynomial[elpStep][syndromeIndex] = antiLogTable[errorLocatorPolynomial[elpStep][syndromeIndex]];
			}

			/* Chien search: find roots of the error location polynomial */
			for (syndromeIndex = 1; syndromeIndex <= elpDegree[elpStep]; syndromeIndex++)
			{
				regPolynomial[syndromeIndex] = errorLocatorPolynomial[elpStep][syndromeIndex];
			}

			int* errorPositions = calloc(t2, sizeof(int));

			errorCount = 0;
			for (syndromeIndex = 1; syndromeIndex <= n; syndromeIndex++) {
				tempIndex = 1;
				for (bitIndex = 1; bitIndex <= elpDegree[elpStep]; bitIndex++)
					if (regPolynomial[bitIndex] != -1) {
						regPolynomial[bitIndex] = (regPolynomial[bitIndex] + bitIndex) % n;
						tempIndex ^= logTable[regPolynomial[bitIndex]];
					}
				if (!tempIndex) {	/* store root and error
						 * location number indices */
					errorPositions[errorCount] = n - syndromeIndex;
					errorCount++;
				}
			}
			//if (errorCount == elpDegree[elpStep])
				/* no. roots = degree of elp hence <= t errors */
				for (syndromeIndex = 0; syndromeIndex < elpDegree[elpStep]; syndromeIndex++)
				{
					//recd[errorPositions[syndromeIndex]] ^= 1;
					SetBitAtIndex(recd, errorPositions[syndromeIndex], GetBitAtIndex(recd, errorPositions[syndromeIndex]) ^ 1);
				}
			

			free(errorPositions);
		}
		else	/* elp has degree >t hence cannot solve */
		{
			//free(recd); recd = 0;
			//Do nothing -> returns noisy code as is
		}

		for (int i = 0; i < t2+2; i++) {
			free(*(errorLocatorPolynomial + i));
		}
		free(errorLocatorPolynomial);
		free(discrepancy);
		free(elpDegree);
		free(elpStepDiff);
		free(regPolynomial);
	}

	free(syndromes);

	//if (recd == 0) return 0;

	uint8_t const KBytes = ceil(k / 8.0);
	uint8_t* Word = malloc(KBytes);

	for (int i = 0; i < k; i++)
	{
		SetBitAtIndex(Word, i, GetBitAtIndex(recd, n - k + i));
	}

	free(recd);
	return Word;
}