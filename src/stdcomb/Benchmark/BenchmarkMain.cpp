// For benchmarking CIdxComb and next_combination()
// This version uses vectors
#include <iostream>
#include <iomanip>
#include <vector>
#include "combination.h"
#include "IndexCombination.h"
#include "FindTotalComb.h"
#include "FindCombByIdx.h"
#include "timer.h"

#include "BigInteger/BigIntegerLibrary.h"

using namespace std;
using namespace stdcomb;

//#define XINT64 __int64
#define XINT64 BigInteger

// Function Prototypes
int CalPercent( double Curr, double Total );

const int SET = 20;
const int COMB = 6;
const int TOTALCOMB = 38760;

long long Beg;  
long long End;  

int main()
{

    timer_init();

//  Benchmark for template function, next_combination()
	vector<unsigned int> ca;
	int a=0;
	for( a=0; a<SET; ++a )
		ca.push_back (a);
	
	vector<unsigned int> cb;
	for( int b=0; b<COMB; ++b )
		cb.push_back (b);
   
	bool bEnd = false;

	int CurrTotal = 1;
	int CurrPercent = 0;
	
    Beg = timer_start();

	const int ITERATION = 10;
	for( a = 0; a<ITERATION; ++a )
	{
		while( !bEnd )
		{
			
		
			if( !next_combination( ca.begin(), ca.end(), cb.begin(),cb.end() ) )
			{
				bEnd  = true;
			}
	/*
			++CurrTotal;

			int Percent = CalPercent( CurrTotal, TOTALCOMB );
			if( CurrPercent != Percent )
			{
				CurrPercent = Percent;
				cout<< "Percentage completed : " << CurrPercent << endl;
			}
	*/
		}
	}

	End = timer_end(Beg);

	cout << setw(40) << "Inteval for next_combination : " << End << endl;

//  Benchmark for CIdxComb class
	CIdxComb<unsigned int> cic(SET,COMB);

	std::vector<unsigned int> vi(COMB);

	int c=0;
	for( c=0; c<COMB; ++c )
		vi[c] = c;

	bEnd = false;

	CurrTotal = 1;
	CurrPercent = 0;
	
	Beg = timer_start();
	for( a=0; a<ITERATION; ++a )
	{
		while( !bEnd )
		{
		
			if( !cic.GetNextComb( vi ) )
			{
				bEnd  = true;
			}
	/*
			++CurrTotal;

			int Percent = CalPercent( CurrTotal, TOTALCOMB );
			if( CurrPercent != Percent )
			{
				CurrPercent = Percent;
				cout<< "Percentage completed : " << CurrPercent << endl;
			}
	*/
		}
	}
	End = timer_end(Beg);

	cout << setw(40) << "Inteval for CIdxComb : " << End << endl;


//  Benchmark for FindCombByIdx function
////////////////////////////////////////////
	std::vector<unsigned int> vi2(COMB);

	for( c=0; c<COMB; ++c )
		vi2[c] = c;

	bEnd = false;

	CurrTotal = 1;
	CurrPercent = 0;

	unsigned int CIndex = 0;
	
	CFindCombByIdx< CFindTotalComb<XINT64>, XINT64 > fcbi;
		
	Beg = timer_start();
	for( a=0; a<ITERATION; ++a )
	{
		while( !bEnd )
		{

			fcbi.FindCombByIdx( SET, COMB, CIndex, vi2 );
		
			++CIndex;
			
			if( CIndex == TOTALCOMB )
			{
				bEnd  = true;
			}
	
			//++CurrTotal;

			//int Percent = CalPercent( CurrTotal, TOTALCOMB );
			//if( CurrPercent != Percent )
			//{
			//	CurrPercent = Percent;
			//	cout<< "Percentage completed : " << CurrPercent << endl;
			//}
		}
	}
	End = timer_end(Beg);

	cout << setw(40) << "Inteval for CFindComb : " << End << endl;
	
	return 0;
}


int CalPercent( double Curr, double Total )
{
	if( Total == 0 )
		throw string("Division by Zero in CalPercent()");

	int Ans = (int) (( Curr/Total ) * 100.0 + 0.51);

	return Ans;
}


