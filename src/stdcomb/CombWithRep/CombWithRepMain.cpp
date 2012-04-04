#include <iostream>
#include <vector>
#include <string>
#include "CombWithRep.h"


using namespace std;
using namespace stdcomb;

int main()
{
	const int SET = 5;
	const int COMB = 3;
	
	// Initialize the first combination vector to zeros
	std::vector<unsigned int> vi( COMB, 0 );

	// Display the first combination
	int Cnt=0;
	{
		cout<<Cnt<<")";
		++Cnt;
		for( int j=0; j<COMB; ++j )
		{
			cout<< vi[j] << ",";
		}
		cout<<endl;
	}
	// Find and display the subsequent combinations
	while( CombWithRep( SET, COMB, vi ) )
	{
		cout<<Cnt<<")";
		for( int j=0; j<COMB; ++j )
		{
			cout<< vi[j] << ",";
		}
		cout<<endl;
		++Cnt;
	}

	cout<<endl;

	return 0;
}

