#include <iostream>
#include <fstream>
#include "json.hpp"
#include "cbraid.h"
#include "braiding.h"

using json = nlohmann::json;
using namespace CBraid;
using CBraid::sint16;  // Avoid ambiguity with CLN.
using namespace Braiding;
using namespace std;

int main(int argc, char* argv[])
{
    json word, wordlcf;
    ArtinBraid B=ArtinBraid(1);
    
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <filename>" << endl;
        return 1;
    }
    ifstream infile(argv[1]);
    if (!infile.is_open()) {
        cerr << "Could not open file: " << argv[1] << endl;
        return 1;
    }

    // Parse JSON
    infile >> word;     // this automatically parses the entire JSON
    infile.close();

    // Check that it is a list (array)
    if (!word.is_array()) {
        cerr << "JSON is not a list.\n";
        return 1;
    }
    B = ArtinBraid( 4 );
	B = WordToBraid( word, 4);
    B.MakeLCF();
    
    sint16 i, j, k, n=B.Index();
    ArtinFactor F=ArtinFactor(n);
    list<sint16> wlcf;
    list<ArtinFactor>::iterator it;

    for(it=B.FactorList.begin(); it!=B.FactorList.end(); it++)
    {
      if(it!=B.FactorList.begin())
	cout << ". ";

      F=*it;
      for(i=2; i<=n; i++)
	{
	  for(j=i; j>1 && F[j]<F[j-1]; j--)
	    {
	      wlcf.push_back( j-1 );
	      k=F[j];
	      F[j]=F[j-1];
	      F[j-1]=k;
	    }
	}
    }

    for (int x : wlcf) {
        wordlcf.push_back(x);
    }

    ofstream outfile(argv[1]);
    if (!outfile.is_open()) {
        cerr << "Could not open file: " << argv[1] << endl;
        return 1;
    }
    // Write in the json
    outfile << wordlcf;
    outfile.close();

    return 0;
}