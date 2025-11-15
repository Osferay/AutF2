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
    ArtinBraid B = ArtinBraid(1), C = ArtinBraid(1);
    list<ArtinBraid> Cent;
    list<ArtinBraid>::iterator it;
    
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

    // Start the braid
    B = ArtinBraid( 4 );
	B = WordToBraid( word, 4);
    Cent = Centralizer( B );

    for(it=Cent.begin(); it!=Cent.end(); it++){
        C = *it;
        sint16 i, j, k, n=C.Index();
        ArtinFactor F=ArtinFactor(n);
        list<sint16> wlcf;
        list<ArtinFactor>::iterator it2;

        if(C.LeftDelta!=0){
            wlcf.push_back( 4 );
        }

        for(it2=C.FactorList.begin(); it2!=C.FactorList.end(); it2++)
        {
            F=*it2;
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

        wordlcf.push_back( wlcf );
        wlcf.clear();
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