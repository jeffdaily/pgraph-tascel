#include <vector>
#include <iostream>
#include "CombFromRepSet.h"

using namespace std;
using namespace stdcomb;


int main()
{
	// Initialize the set
	vector<unsigned int> svi;
	svi.push_back(0); // 0
	svi.push_back(1); // 1
	svi.push_back(2); // 2
	svi.push_back(3); // 3
	svi.push_back(3); // 4
	svi.push_back(3); // 5
	svi.push_back(4); // 6
	svi.push_back(5); // 7

	// Object to find the combinations from set with repeated elements
	CCombFromRepSet cfrs;

	cfrs.SetRepeatSetInfo( svi );
	
	// Set the size of Set and number elements of the combination
	const unsigned int SET = 8;
	const unsigned int COMB = 3;
	cfrs.SetSizes( SET, COMB );

	// Initialize the first combination vector
	vector<unsigned int> vi;
	for( unsigned int j=0; j<COMB; ++j )
		vi.push_back( j );

	// Set the first combination
	cfrs.SetFirstComb( vi );

	// Display the first combination
	int Cnt=0;
	{
		cout<<Cnt<<")";
		++Cnt;
		for( unsigned int i=0; i<vi.size(); ++i)
		{
			cout<<svi[vi[i]]<<",";
		}
		cout<<endl;
	}

	// Find and display the subsequent combinations
	while( cfrs.GetNextComb( vi ) )
	{
		cout<<Cnt<<")";
		++Cnt;
		for( unsigned int i=0; i<vi.size(); ++i)
		{
			cout<<svi[vi[i]]<<",";
		}
		cout<<endl;
	}

	return 0;
}

