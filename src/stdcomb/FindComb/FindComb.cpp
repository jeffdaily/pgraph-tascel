// Note not to be included in the article!


#include <iostream>
#include <string>

#include "FindCombByIdx.h"
#include "FindTotalComb.h"
#include "IndexCombination.h"
#include "BigInteger/BigIntegerLibrary.h"

using namespace std;
using namespace stdcomb;

#define USE_BIGINTEGER 0

int main(int argc, char* argv[])
{
#if USE_BIGINTEGER
	CFindCombByIdx< CFindTotalComb<BigInteger>, BigInteger > findcomb;
#else
	CFindCombByIdx< CFindTotalComb<unsigned long>, unsigned long > findcomb;
#endif

	const unsigned int nComb = 2;
	const unsigned int nSet = 640000;
	
	// Intialize the vector with size nComb
	vector<unsigned int> vec(nComb);
	
	// vector returned by FindCombByIdx
	findcomb.FindCombByIdx( nSet, nComb, 91, vec );

	cout<< vec.at(0) << ",";
	cout<< vec.at(1) << ",";
	cout<< vec.at(2) << ",";
    cout<< endl;

	return 0;
}

