//#define THOROUGH_CHECK

#include <fstream>
#include <string>
#include <set>
#include <cctype>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <numeric>

#include <cstdio>
#include <array>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <map>
#include <limits>
#include <cassert>

//#include "functions2.hpp"
#include "functions2.cpp"

typedef std::numeric_limits< double > dbl;

using std::vector;
using std::cout;
using std::cerr;
using std::fill;
using std::map;
using std::swap;


// Main function
int main(int argc, char *argv[]) {
    // Initialization of variables...
    char *p, *pp, *ppp, *pppp, *ppppp, *pppppp, *ppppppp, *pppppppp, *ppppppppp,*pppppppppp;
    // Processing the input arguments
    // taking the filename string ...
    char const *tableofranksname (argv[1]);

    // reading how many starts per run without randomizing
    const int Nstarts = (int) strtol(argv[2], &pp, 10);

    // reading the desired number of genes in the resulting solutions
    const int DesiredNumberOfGenes = (int) strtol(argv[3], &ppp, 10);
    
    const long randomSeed = strtol(argv[4], &pppppppppp, 10);


    std::random_device rd; // will be used to obtain a seed for the random number engine
    //static const long long randseed(rd());
    static const long long randseed(randomSeed);
    //static const long long randseed(1876164852);
    std::mt19937_64 rng(randseed); // standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> runif(0.0, 1.0); // random draw from a uniform distribution.

    
    //Turning the inputs into tractable data
    vector<vector<double> > r;
    vector<vector<int> > StartingSignatures;
    vector<std::string> gene_names;
    vector<std::string> cellaffiliations_raw;
    readIn2dData(tableofranksname, r, gene_names, cellaffiliations_raw); // ... turning rank matrix into a data array
    //vector<int> cellaffiliations (readInAffiliationsVector(vectoraffiliationsname));// .. turning cell affiliations to tumors vector into a vector
    //cellaffiliations.push_back(cellaffiliations[cellaffiliations.size()-1]+1); // adding sentinel (final sample marker) at the end
    cellaffiliations_raw.push_back("END_POINTER"); // adding sentinel (final sample marker) at the end
    vector<int> cellaffiliations = ConvertAffiliationsToInts(cellaffiliations_raw);
    
    if(gene_names.size()==0) std::iota(gene_names.begin(),gene_names.end(), 0);

    int k = (int) r.size(); // total number of genes (rows)
    if(k < 10){
        throw std::invalid_argument( "First argument should be a data table of at least 10 rows" ); // checking for the input table to be of a reasonably large size
    }
    
    /*
     ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     Initializing all tumors with random gene indices
     ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     */
    /***** COMPUTE TIES AND SUM OF SQUARED RANKS FOR EACH GENE FOR EACH SAMPLE *****/
    // before permuting the input matrix we can precompute our matrix of tie-adjustments, as it won't change even after a permutation
    // T_ij is the matrix of tie adjustments
    vector<vector <double> > SSR, T_ij; // outputs
    ComputeSquaredRanksAndTieAdjustments(SSR, T_ij, r, cellaffiliations);

    /***** COMPUTE THE SAMPLES' BOUNDARIES *****/
    vector<int> firstCellInSample = SampleBoundaries(cellaffiliations);

    /***** IF RANDOMIZE THE INPUT TABLE *****/
//    if (permute) PermuteRanksWithinEachGeneInEachSample(r, firstCellInSample, runif, rng);
    const int factor=0;

    vector <int> BannedGenes;
    //vector <vector <int>> AlreadySeenSolutions;

    for (int startcounter = 0; startcounter < Nstarts; startcounter++){ // loop for many starts
        // Repeat for the whole range from minExp to (including) maxExp
            vector <double> logpvals, logpvals_0; // current one and previous one.
            vector<int> X, X_0, // States = selected + unselected genes + CurrentNumberOfGenes ("ordered" = selected genes first)
                        G, // signature from which to start the search
                        m_i; // sizes of samples
            std::uniform_real_distribution<> runif(0.0, 1.0); // random draw from a uniform distribution.
            
            vector<int> allgenes(r.size()), startingseed;
            std::iota(allgenes.begin(),allgenes.end(), 0); // integer vector of 0, 1, 2, 3 ... , 1751(-ish). These correspond to gene indices
            
            for (int ll = 0; ll < BannedGenes.size(); ++ll)
                RemoveGeneFromVector(BannedGenes[ll], allgenes);
            //k = allgenes.size();
            
            for (int kek = 0; kek < DesiredNumberOfGenes; kek++){
                int newdraw = kek + floor(runif(rng) * (r.size()-kek));
                //int newdraw = Last50Genes[kek];
                startingseed.push_back(allgenes[newdraw]); // startingseed is an integer vector of length howmanystartinggenes and contains randomly drawn numbers from the 1:1752(-ish) range without replacement
                swap(allgenes[newdraw],allgenes[kek]);
            }
            
            allgenes.push_back(DesiredNumberOfGenes); // adding a partitioning index at the end of the state vector
            //G = allgenes;
            
            
            X_0 = ConstructStateVector(startingseed, r.size(), BannedGenes);
            G = X_0;
            k = X_0.size()-1;
            
            /***** COMPUTE INITIAL INTERMEDIATE FRIEDMAN STATISTIC MEMBERS OF THE EQUATION *****/
            
            vector<double> A, B, C, D, Rc, GG, F, E, S_hats;
            
            double g = InitializeFriedmanMembers(G, A, B, C, D, Rc, GG, F, E, S_hats, logpvals, m_i, r, cellaffiliations, T_ij, factor);

            /***** SIMULATED ANNEALING RUN *****/
            double startingloss = g;
            //The main simulated annealing (SA) part:
            // Variables pertaining to the main SA loop:
            double g_opt = 0.0, g_0 = 0, Tmin = 1e-319;
            int nIter = 0;
            vector <double> A_0 = A, B_0 = B, C_0 = C, D_0 = D, logpvals_opt=logpvals, Rc_0 = Rc, S_hats_0 = S_hats, S_opt = S_hats, A_opt = A, B_opt = B, C_opt = C, D_opt = D, Rc_opt = Rc;
            X = X_0; // currently manipulated solution (and previous one)
            vector <int> X_opt = X;
            double max_loss = 0; // current global best solution value
            int rejectioncounter=0;
            
            //std::ofstream myfile;
            //myfile.open ("/Users/emilezakiev/TRACEofSINGLErun.txt");
            
            //std::ofstream myfile;
            //myfile.open ("/Users/emilezakiev/rngs.txt");
            
            //std::ofstream myfile;
            //myfile.open ("/Users/emilezakiev/CurrentSignatureX.txt");
            //for (int uu = 0; uu < X[k]; ++uu) myfile << X[uu] << "\t";
            //myfile.close();
            
            
            //myfile.close();
            

            /***** PRINTING THE FINAL RESULTS *****/
            //double TotalSum = 0;
            logpvals_opt = logpvals; X_opt = X; S_opt = S_hats; g_opt = g;
            
            // we need to do the final neighbourhood-enumerating search after even the greedy search
            
            
            
            vector <double> A_opt_save = A_opt,
            B_opt_save = B_opt,
            C_opt_save = C_opt,
            D_opt_save = D_opt,
            Rc_opt_save = Rc_opt,
            logpvals_opt_save = logpvals_opt,
            S_opt_save = S_opt;
            vector <int> X_opt_save = X_opt;
            
            double biggest_neighbour = g;
            
            int iter_local_neighbourhood = 0;
            
            while(true){
                iter_local_neighbourhood++;
                vector <vector <double> > neighbourhood;
                for (int genecounter=0; genecounter < X_opt.back(); ++genecounter){
                    int GeneToSwapOut = X_opt[genecounter];
                    vector <double> singleGeneNeighbourhood;
                    for (int replacementIter=0; replacementIter < k-X_opt.back()-1; ++replacementIter){
                        int indexSwapIn =  X_opt.back()+replacementIter;
                        int GeneToSwapIn = X_opt[indexSwapIn];
                        swap(X_opt[genecounter], X_opt[indexSwapIn]);
                        // Calculating the loss function at the new state
                        g = objectivefunction_swap(GeneToSwapIn, GeneToSwapOut, A_opt, B_opt, C_opt, D_opt, T_ij, SSR, Rc_opt, GG, F, E, logpvals_opt, cellaffiliations, r, m_i, X_opt.back(), S_opt, factor, X_opt);
                        singleGeneNeighbourhood.push_back(g);
                        A_opt = A_opt_save;
                        B_opt = B_opt_save;
                        C_opt = C_opt_save;
                        D_opt = D_opt_save;
                        Rc_opt = Rc_opt_save;
                        logpvals_opt = logpvals_opt_save;
                        S_opt = S_opt_save;
                        X_opt = X_opt_save;
                    }
                    neighbourhood.push_back(singleGeneNeighbourhood);
                }
                
                int biggest_neighbour_i = -1, biggest_neighbour_j = -1;
                for (int ii=0; ii < neighbourhood.size(); ++ii)
                    for (int jj=0; jj < neighbourhood[0].size(); ++jj)
                        if (neighbourhood[ii][jj] > biggest_neighbour){
                            biggest_neighbour = neighbourhood[ii][jj];
                            biggest_neighbour_i = ii;
                            biggest_neighbour_j = jj;
                        }
                
                if(biggest_neighbour_i == -1) break;
                
                
                g_opt = objectivefunction_swap(X_opt[biggest_neighbour_j+X_opt.back()], X_opt[biggest_neighbour_i], A_opt, B_opt, C_opt, D_opt, T_ij, SSR, Rc_opt, GG, F, E, logpvals_opt, cellaffiliations, r, m_i, X_opt.back(), S_opt, factor, X_opt);
                
                swap(X_opt[biggest_neighbour_i], X_opt[biggest_neighbour_j+X_opt.back()]);
                A_opt_save = A_opt;
                B_opt_save = B_opt;
                C_opt_save = C_opt;
                D_opt_save = D_opt;
                Rc_opt_save = Rc_opt;
                logpvals_opt_save = logpvals_opt;
                S_opt_save = S_opt;
                X_opt_save = X_opt;
                
            }
            
            
            
            /*
            for (auto it = logpvals_opt_save.begin(); it < logpvals_opt_save.end(); it++) {
                TotalSum += *it;
            }
             */
            
            
            cout << X_opt.back() << "\t|\t";; // printing out the exponent and how many genes are in the optimal solution
            //cout << TotalSum << "\t|\t"; // printing out the cumulative (LogPval)^2
            cout << round2dec(g_opt) << "\t|\t"; // printing out the cumulative (LogPval)^2
            
            /*
            for (int nnn=0; nnn < S_opt.size(); nnn++){
                cout << S_opt[nnn] << "\t"; // printing out the S_opt of each tumor separately
            }
            cout << "\t|\t";
             */
            
            for (int nnn=0; nnn < logpvals_opt.size(); nnn++){
                cout << round2dec(logpvals_opt[nnn]) << "\t"; // printing out the (LogPval)^2 of each tumor separately
            }
            cout << "|\t";
            
            vector<int> GenesToPrint = {X_opt.begin(), X_opt.begin() + X_opt.back()};
            sort(GenesToPrint.begin(), GenesToPrint.end());
            
            //AlreadySeenSolutions.push_back(GenesToPrint);
            vector <int> v=GenesToPrint;
            std::shuffle(v.begin(), v.end(), rng);
            //BannedGenes.push_back(v[0]);

            for (int nnn=0; nnn < GenesToPrint.size(); nnn++){
                cout << gene_names[GenesToPrint[nnn]] << "\t"; // printing out the gene indices in sorted order
            }
            
            cout << "\n";

            //myfile.close();
    } // end of the loop for multiple starts of SA

    return 0;
}
