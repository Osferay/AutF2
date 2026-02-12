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
    ArtinBraid B1 = ArtinBraid(4), B2 = ArtinBraid(4), C = ArtinBraid(4);
    sint16 i, j, k, n;
    list<sint16> wlcf;
    list<ArtinFactor>::iterator it;
    list<ArtinBraid>::iterator cit;
    list<ArtinBraid> Cent;
    
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
    infile >> word;
    infile.close();

    string func = argv[2];

    // Test cases to see what we have to do
    if( func == "lcf" ){
        
        // In this case we compute the left canonical form
        B1 = WordToBraid( word, 4 );
        B1.MakeLCF();
        n = B1.Index();
        ArtinFactor F = ArtinFactor(n);

        // We translate the left canonical form to a readeable list.
        // 4 represents the element delta, since we are in B_4/Z(B_4) this can only be one.
        if(B1.LeftDelta % 2 != 0){
            wlcf.push_back( 4 );
        }
    
        for(it=B1.FactorList.begin(); it!=B1.FactorList.end(); it++)
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
        if( wlcf.size() == 0 ){
            wordlcf.push_back(0);
        }

        // We put all the elements into a json list.
        for (int x : wlcf) {
            wordlcf.push_back(x);
        }
    }
    else if( func == "conj" ){
        std::list<int> word1 = word["word1"].get<std::list<int>>();
        std::list<int> word2 = word["word2"].get<std::list<int>>();

    	B1 = WordToBraid( word1, 4 );
	    B2 = WordToBraid( word2, 4 );

        if( AreConjugate( B1, B2, C ) ){
        // Read all the factors and the deltas
        n=C.Index();
        ArtinFactor F=ArtinFactor(n);
        
            if(C.LeftDelta % 2 != 0){
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
            if( wlcf.size() == 0 ){
                wordlcf.push_back(0);
            }
            
            for (int x : wlcf) {
                    wordlcf.push_back(x);
            }
        }
        else{
            wordlcf.push_back( false );
        }
    }
    else if( func == "cent" ){
        
        B1 = WordToBraid( word, 4);
        Cent = Centralizer( B1 );

        for(cit=Cent.begin(); cit!=Cent.end(); cit++){
            C = *cit;
            n = C.Index();
            ArtinFactor F=ArtinFactor(n);
            list<ArtinFactor>::iterator it2;

            if(C.LeftDelta % 2 != 0){
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
    }
    else{
        cerr << "Not implemented routine: " << argv[2] << endl;
    }
    
    ofstream outfile(argv[1]);
    if (!outfile.is_open()) {
        cerr << "Could not open file: " << argv[1] << endl;
        return 1;
    }
    
    outfile << wordlcf;
    outfile.close();

    return 0;
}