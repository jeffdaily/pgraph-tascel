// IntComb.cpp : Defines the entry point for the console application.
//

#include "IndexCombination.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace stdcomb;


int main(int argc, char* argv[])
{
	CIdxComb<unsigned int> cb;

	cb.SetSizes(5,3);


	vector<string> vsAnimal;
	vsAnimal.push_back( "Elephant" );
	vsAnimal.push_back( "Deer" );
	vsAnimal.push_back( "Cat" );
	vsAnimal.push_back( "Bear" );
	vsAnimal.push_back( "Ape" );

	vector<unsigned int> vi(3);

	vi[0] = 0;
	vi[1] = 1;
	vi[2] = 2;

	cout<< vsAnimal[ vi[0] ] << " " 
		<< vsAnimal[ vi[1] ] << " " 
		<< vsAnimal[ vi[2] ] << "\n";

	int Total = 1;
	while ( cb.GetNextComb( vi ) )
	{
		// Do whatever processing you want
		
		
		cout<< vsAnimal[ vi[0] ] << " " 
			<< vsAnimal[ vi[1] ] << " " 
			<< vsAnimal[ vi[2] ] << endl; 
		
		++Total;
	}

	cout<< "\nTotal : " << Total << endl;
	return 0;
}
