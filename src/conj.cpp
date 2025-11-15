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
    json words, wordlcf;
    ArtinBraid B1 = ArtinBraid(1), B2 = ArtinBraid(1), C = ArtinBraid(1);
    sint16 i,j;
    
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
    infile >> words;     // this automatically parses the entire JSON
    infile.close();

    // Get the two lists
    std::list<int> word1 = words["word1"].get<std::list<int>>();
    std::list<int> word2 = words["word2"].get<std::list<int>>();

    // Start the braid
    B1 = ArtinBraid( 4 );
    B2 = ArtinBraid( 4 );
    C  = ArtinBraid( 4 );
	B1 = WordToBraid( word1, 4 );
	B2 = WordToBraid( word2, 4 );
    
    if( AreConjugate( B1, B2, C ) ){
        // Read all the factors and the deltas
        sint16 i, j, k, n=C.Index();
        ArtinFactor F=ArtinFactor(n);
        list<sint16> wlcf;
        list<ArtinFactor>::iterator it;
        
        if(C.LeftDelta!=0){
            wlcf.push_back( 4 );
        }
        
        for(it=C.FactorList.begin(); it!=C.FactorList.end(); it++)
        {
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
    }
    else{
        wordlcf.push_back( false );
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