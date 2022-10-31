#include "functions2.hpp"
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <limits>
#include <numeric>
#include <random>
#include "pgamma_from_R.cpp"


using std::vector;
using std::cout;

void ConvertMatrixToRanks(vector<vector<double> >& expressionmatrix, const vector <int> & cellaffs){
    
}


void RemoveGeneFromVector(int GeneToRemove, vector <int> &X){
    std::vector<int>::iterator position = std::find(X.begin(), X.end(), GeneToRemove);
    if (position != X.end()) // == myVector.end() means the element was not found
        X.erase(position);
}

int RandomlySelectGeneToBan(vector <int> v, std::mt19937_64 rng){
    std::shuffle(v.begin(), v.end(), rng);
    return v[0];
}

vector <int> ConvertAffiliationsToInts(const vector<std::string> cellaffiliations_raw){
    int counter = 0;
    vector <int> cellaffiliations;
    int subcounter = 0;
    for (int i=0; i < cellaffiliations_raw.size()-1; ++i){
        cellaffiliations.push_back(counter);
        subcounter++;
        if(cellaffiliations_raw[i]!=cellaffiliations_raw[i+1]){
            counter++;
            assert(subcounter>9);
            subcounter = 0;
        }
    }
    cellaffiliations.push_back(counter+1);
    return cellaffiliations;
}

double sqrt_logchisq_R(const double score, const int m_i){ // wrapper function that divides the number of degrees of freedom by 2 before feeding into R's normalized incomplete gamma function calculation
        double res = pgamma(score, m_i/2.,2., 0, 1);
        //return (-res)*(-res);
        return sqrt(sqrt(-res));
}

// Auxilary function that takes a char* filename argument and returns a 2d dynamic array containing the data
void readIn2dData(const char* filename, vector<vector<double> >& ranktable, vector<std::string> & gene_names, vector<std::string> & cell_affiliations){
    std::fstream ifs;
    /* open file */
    ifs.open(filename, std::fstream::in);
    
    std::string line;
    std::string charbuf;
    getline(ifs, line);
    std::stringstream ss(line, std::ios_base::out|std::ios_base::in|std::ios_base::binary);
    
    while (ss >> charbuf) cell_affiliations.push_back(charbuf);
    
    while (true){
        std::string line;
        double buf;
        std::string charbuf;
        getline(ifs, line);
        std::stringstream ss(line, std::ios_base::out|std::ios_base::in|std::ios_base::binary);
        if (!ifs){
            // mainly catch EOF
            break;
        }
        if (line[0] == '#' || line.empty())
            // catch empty lines or comment lines
            continue;
        if(!isdigit(line[0])){
            ss >> charbuf;
            gene_names.push_back(charbuf);
        }
        
        vector<double> row;
        
        while (ss >> buf)
                row.push_back(buf);
        ranktable.push_back(row);
        assert(row.size()==cell_affiliations.size());
    }
    if(gene_names.size()>0) assert(gene_names.size()==ranktable.size());
    ifs.close();
}

// Read a table of integers from a space-separated file
void readIn2dData_int(const char* filename, vector<vector<int> >& ranktable){
    std::fstream ifs;
    /* open file */
    ifs.open(filename, std::fstream::in);
    while (true){
        std::string line;
        int buf;
        getline(ifs, line);
        std::stringstream ss(line, std::ios_base::out|std::ios_base::in|std::ios_base::binary);
        if (!ifs){
            // mainly catch EOF
            break;
        }
        
        if (line[0] == '#' || line.empty())
            // catch empty lines or comment lines
            continue;
        vector<int> row;
        while (ss >> buf){
            
            row.push_back(buf);
        }
        ranktable.push_back(row);
    }
    ifs.close();
}


// Smart-update objective function.
/**
 * @param newGene index of gene to be added or substracted
 * @param add if true, add the newGene, else remove it
 * @param A parameters of the Friedman eq.
 * @param B parameters of the Friedman eq.
 * @param C parameters of the Friedman eq.
 * @param D parameters of the Friedman eq.
 * @param T_ij parameters of the Friedman eq.
 * @param SSR parameters of the Friedman eq.
 * @param Rc parameters of the Friedman eq.
 * @param GG parameters of the Friedman eq.
 * @param F parameters of the Friedman eq.
 * @param E parameters of the Friedman eq.
 * @param logpvals log p-values for each sample
 * @param cellaffiliations which cell came from which sample
 * @param r input rank table
 * @param m_i number of cells (for each samples)
 * @param n number of genes in the current solution
 * @param M number of cells in all samples (cumulatively)
 * @param[out] S_hats calculated Frideman statistic values
 * @param factor exponent in the size-dependent penalization factor for the statistic
 * @param PriorState keep track of the previous solution
 */
double objectivefunction(
                    int newGene,
                    bool add,
                    vector<double> &A,
                    vector<double> &B,
                    vector<double> &C,
                    vector<double> &D,
                    const vector<vector<double> > &T_ij,
                    const vector<vector<double> > &SSR,
                    vector<double> &Rc,
                    const vector<double> &GG,
                    const vector<double> &F,
                    const vector<double> &E,
                    vector<double> &logpvals,
                    const vector<int> &cellaffiliations,
                    const vector<vector<double> > &r,
                    const vector<int> m_i,
                    const int n,
                    vector<double> &S_hats,
                    const double factor,
                    const vector<int> &PriorState // State: < selected_gene_0, ..., selected_gene_k, unselected_gene_k+1, ..., unselected_gene_M, CurrentNumberOfGenes>
                ){
    double g = 0.0;
    int whichTumor = 0; // Current sample
    int M = (int) cellaffiliations.size();
    if(add){
        // adding newGene (using += and +1.0)
        // goes through all the cells in all samples
        for (int ll = 0; ll < M-1; ll++){
            A[whichTumor] += 24.0 * Rc[ll] * r[newGene][ll];
            Rc[ll] += r[newGene][ll];

            // If we're going to another sample.
            if(cellaffiliations[ll] != cellaffiliations[ll+1]){
                A[whichTumor]+= 12.0 * SSR[newGene][whichTumor];
                B[whichTumor]+= (E[whichTumor])*(2.0*n + 1.0);
                C[whichTumor]+= F[whichTumor];
                D[whichTumor]+= (T_ij[newGene][whichTumor]/(m_i[whichTumor]-1.0) - GG[whichTumor]);
                const double denominator = (C[whichTumor]-D[whichTumor]);

                // Failsafe in case the denominator goes to zero.
                // Happens when all the genes have the same rank in all the cells.
                // FIXME log when it happens.
                // FIXME check if `long double` would be better.
                double s = (denominator==0) ? 0 : (A[whichTumor]-B[whichTumor])/denominator;
                double s_hat = s/pow(n+1.0, factor);
                if ((isnan(s_hat)) || (s_hat < 0.0) || (isinf(s_hat))){ // There is a computation error, should not happen.
                    // FIXME use assert
                    cout.precision(std::numeric_limits<double>::digits10 + 1);
                    //cout << "Error: randseed = " << randseed << ": ";
                    cout << "S' wasn't good| No. of genes is " << n << " | the genes are: ";
                    for (int k = 0; k < n; k++) cout << PriorState[k] <<"\t";
                    cout << "| attempted to add gene " << newGene << " to tumor " << whichTumor;
                    cout << "| A: " << A[whichTumor] << " | B: " << B[whichTumor] << " | C: " << C[whichTumor] << " | D: " << D[whichTumor];
                    cout <<"\n";
                    s_hat = 0.0;
                }
                
                S_hats[whichTumor] = s_hat;
                logpvals[whichTumor] = sqrt_logchisq_R(s_hat, m_i[whichTumor]);
                g += logpvals[whichTumor];
                whichTumor++;
            } // if cellaffiliations
        } // for ll in M-1

    } else {
        // subtracting newGene (using -= and -1.0)
        for (int ll = 0; ll < M-1; ll++){
            A[whichTumor] -= 24.0*Rc[ll]*r[newGene][ll];
            Rc[ll] -= r[newGene][ll];
            if(cellaffiliations[ll] != cellaffiliations[ll+1]){
                A[whichTumor]+= 12.0*SSR[newGene][whichTumor];
                B[whichTumor]-= (E[whichTumor])*(2.0*n - 1.0);
                C[whichTumor]-= F[whichTumor];
                D[whichTumor]-= (T_ij[newGene][whichTumor]/(m_i[whichTumor]-1.0) - GG[whichTumor]);
                const double denominator = (C[whichTumor]-D[whichTumor]);
                double s = (denominator==0) ? 0 : (A[whichTumor]-B[whichTumor])/denominator;
                double s_hat = s/pow(n-1.0, factor);
                
                if ((isnan(s_hat)) || (s_hat < 0.0) || (isinf(s_hat))) { // Error that should never happen.
                    // FIXME use assert
                    cout.precision(std::numeric_limits<double>::digits10 + 1);
                    //cout << "Error: randseed = " << randseed << ": ";
                    cout << "No. of genes is " << n << " | the genes are: ";
                    for (int k = 0; k < n; k++) cout << PriorState[k] <<"\t";
                    cout << "| attempted to subtract gene " << newGene << " from tumor " << whichTumor;
                    cout << "| A: " << A[whichTumor] << " | B: " << B[whichTumor] << " | C: " << C[whichTumor] << " | D: " << D[whichTumor];
                    cout <<"\n";
                    s_hat = 0.0;
                }
                
                S_hats[whichTumor] = s_hat;
                logpvals[whichTumor] = sqrt_logchisq_R(s_hat, m_i[whichTumor]);
                g += logpvals[whichTumor];
                whichTumor++;
            }
        }
    }
    return g;
}



/** Auxilary function that takes mostly references to existing parts of the Friedman statistic and returns a sum of squares of log(pvalues)
 * Same as the previous objectivefunction, but for swapping two genes instead of adding/substracting a gene.
 *
 * @param geneToSwapIn
 * @param geneToSwapOut
 */
double objectivefunction_swap(int geneToSwapIn,
                         int geneToSwapOut,
                         vector<double> &A,
                         vector<double> &B,
                         vector<double> &C,
                         vector<double> &D,
                         const vector<vector<double> > &T_ij,
                         const vector<vector<double> > &SSR,
                         vector<double> &Rc,
                         const vector<double> &GG,
                         const vector<double> &F,
                         const vector<double> &E,
                         vector<double> &logpvals,
                         const vector<int> &cellaffiliations,
                         const vector<vector<double> > &r,
                         const vector<int> m_i,
                         const int n,
                         vector<double> &S_hats,
                         const double factor,
                         const vector<int> &PriorState){
    double g = 0.0;
    int whichTumor = 0;
    int M = (int) cellaffiliations.size();
    
        // adding geneToSwapIn and immediately subtracting geneToSwapOut
        for (int ll = 0; ll < M-1; ll++){
            A[whichTumor] -= 24.0*Rc[ll]*r[geneToSwapOut][ll];
            Rc[ll] -= r[geneToSwapOut][ll];
            
            A[whichTumor] += 24.0*Rc[ll]*r[geneToSwapIn][ll];
            Rc[ll] += r[geneToSwapIn][ll];
            
            if(cellaffiliations[ll] != cellaffiliations[ll+1]){
                A[whichTumor]+= 12.0*SSR[geneToSwapIn][whichTumor];
                A[whichTumor]+= 12.0*SSR[geneToSwapOut][whichTumor];
                // Updates that are not necessary:
                //B[whichTumor]+= (E[whichTumor])*(2.0*n + 1.0);
                //B[whichTumor]-= (E[whichTumor])*(2.0*n - 1.0);
                //C[whichTumor]+= F[whichTumor];
                //C[whichTumor]-= F[whichTumor];
                D[whichTumor]+= (T_ij[geneToSwapIn][whichTumor]/(m_i[whichTumor]-1.0));
                D[whichTumor]-= (T_ij[geneToSwapOut][whichTumor]/(m_i[whichTumor]-1.0));
                const double denominator = (C[whichTumor]-D[whichTumor]);
                double s = (denominator==0) ? 0 : (A[whichTumor]-B[whichTumor])/denominator;
                double s_hat = s/pow(n, factor);
                
                if ((isnan(s_hat)) || (s_hat < 0.0) || (isinf(s_hat))){
                    // FIXME use assert
                    cout.precision(std::numeric_limits<double>::digits10 + 1);
                    //cout << "Error: randseed = " << randseed << ": ";
                    cout << "S' wasn't at all good| No. of genes is " << n << " | the genes are: ";
                    for (int k = 0; k < n; k++) cout << PriorState[k] <<"\t";
                    cout << "| attempted to swap gene " << geneToSwapOut << " for " << geneToSwapIn << " to tumor " << whichTumor;
                    cout << "| A: " << A[whichTumor] << " | B: " << B[whichTumor] << " | C: " << C[whichTumor] << " | D: " << D[whichTumor];
                    cout <<"\n";
                    s_hat = 0.0;
                }
                
                S_hats[whichTumor] = s_hat;
                logpvals[whichTumor] = sqrt_logchisq_R(s_hat, m_i[whichTumor]-1);
                g += logpvals[whichTumor];
                whichTumor++;
            }
        }
    return g;
}

vector <int> SampleBoundaries(const vector <int> cellaffiliations){
    vector <int> firstCellInSample;
    firstCellInSample.push_back(0);
    for (int tt1 = 0; tt1 < cellaffiliations.size() - 1; tt1++){// compute firstCellInSample vector
        if(cellaffiliations[tt1] != cellaffiliations[tt1+1]) firstCellInSample.push_back(tt1+1);
    }
    return firstCellInSample;
}

void ComputeSquaredRanksAndTieAdjustments(vector<vector <double> > &SSR,
                                          vector<vector <double> > &T_ij,
                                          const vector<vector <double> > &input_table,
                                          const vector <int> &cellaffiliations
                                          ){
    for (int line = 0; line < input_table.size(); line++) { // all genes
        vector <double> TempTn, temprjsq, buf;
        for (int i = 0; i < cellaffiliations.size()-1; i++) { // all cells
            buf.push_back(input_table[line][i]); // running sum
            if (cellaffiliations[i] != cellaffiliations[i+1]){// After this, the sample is done
                std::sort(buf.begin(), buf.end());
                int tieCounter = 1;
                double accumulatedSum = 0;
                double rjsqsums = 0;
                buf.push_back(0); // add sentinel at the end
                double sampleRankSum=0; // counting sum of all ranks
                for (int j = 0; j < (buf.size()-1); j++){ // all cells of this sample
                    sampleRankSum+=buf[j];
                    rjsqsums += buf[j]*buf[j];
                    if(buf[j] != buf[j+1]){ // not a new tie ?
                        accumulatedSum += tieCounter*tieCounter*tieCounter; // dumping accumulated tiecounter to the cumulative sum
                        tieCounter = 1; //resetting the tie counter to 1
                    } else { // extend current tie
                        tieCounter++; // incrementing tie counter
                    }
                } // all sample cells done
                if(sampleRankSum != (buf.size()*(buf.size()-1)/2)){
                    for (int j = 0; j < (buf.size()-1); j++){ // all cells of this sample
                        cout << buf[j] << ",";
                    }
                    cout << "\n";
                    throw "Houston we have a problem in the ranks.\n";
                    //return 0;
                }
                TempTn.push_back(accumulatedSum);
                temprjsq.push_back(rjsqsums);
                buf.clear();
            } // new sample done
        } // new cell done
        T_ij.push_back(TempTn);
        SSR.push_back(temprjsq);
    } // gene done
}

void PermuteRanksWithinEachGeneInEachSample(
                                            vector<vector<double> > &r,
                                            vector<int> &firstCellInSample,
                                            std::uniform_real_distribution <> &runif,
                                            std::mt19937_64 rng){
    for (int iter = 0; iter < firstCellInSample.size()-1; iter++){
        for (int line = 0; line < r.size(); line ++){
            double sampleRankSum=0; int numCellsInSample=0;
            for (int count = firstCellInSample[iter]; count < (firstCellInSample[iter+1]); count ++){
                int swapIndex = count + floor(runif(rng) * (firstCellInSample[iter+1]-count));
                std::swap(r[line][swapIndex], r[line][count]);
                sampleRankSum+=r[line][count]; numCellsInSample++;
            }
            // Sanity check
            if(sampleRankSum != (numCellsInSample * (numCellsInSample+1) / 2)){
                cout << iter << "\n";
                throw "Houston we have a problem in the ranks.\n";
                //return 0;
                
            }
        }
    }

}

double KnownStartingTemperature(const int x){
    // Max temperature depends on the penalization factor
    int joop =(int)(x*100);
                // Max temperature depends on the penalization factor
                switch (joop) {
                    case    10  :   return 4.677445e+07;
                    case    15  :     return 1.818510e+07;
                    case    20    :    return 6.97E+06;
                    case    25    :    return 2.66E+06;
                    case    30    :    return 9.83E+05;
                    case    35    :    return 3.52E+05;
                    case    40    :    return 1.16E+05;
                    case    45    :    return 3.42E+04;
                    case    50    :    return 8.82E+03;
                    case    55    :    return 3.21E+03;
                    case    60    :    return 3.52E+03;
                    case    65    :    return 9.21E+02;
                    case    70    :    return 1.68E+02;
                    case    75    :    return 5.41E+01;
                    case    80    :    return 6.83E+00;
                    case    85    :    return 7.79E-01;
                    case    90    :    return 7.15E-02;
                    case    95    :    return 1.44E-04;
                    case    100    :    return 2.88E-07;
                    case    105    :    return 4.54E-10;
                    case    110    :    return 5.74E-13;
                    case    115    :    return 7.79E-16;
                    case    120    :    return 1.34E-18;
                    case    125    :    return 2.12E-21;
                    case    130    :    return 3.97E-24;
                    case    135    :    return 5.83E-27;
                    case    140    :    return 7.83E-30;
                    case    145    :    return 1.57E-32;
                    case    150    :    return 5.22E-37;
                    default:
                        return 5e+07;
                }
}

vector <int> ConstructStateVector(vector <int> Geneset,
                                int k, vector <int> banned_genes){
    std::sort(Geneset.begin(), Geneset.end());
    vector<int> allgenes(k), X_0;
    std::iota(allgenes.begin(),allgenes.end(), 0); // integer vector of 0, 1, 2, 3 ... , 1751(-ish). These correspond to gene indices

    // Build up the state corresponding to the asked signature
    for (int ll = 0; ll < banned_genes.size(); ++ll)
        RemoveGeneFromVector(banned_genes[ll], allgenes);
    int howmanystartinggenes = (int) Geneset.size();
    
    for (int kek = 0; kek < Geneset.size(); ++kek){
        
        std::swap(allgenes[kek],allgenes[Geneset[kek]]);
    }
    allgenes.push_back(howmanystartinggenes); // adding a partitioning index at the end of the state vector
    X_0 = allgenes;
    return X_0;
}

double InitializeFriedmanMembers(const vector<int> &G,
                                 vector<double> &A,
                                 vector<double> &B,
                                 vector<double> &C,
                                 vector<double> &D,
                                 vector<double> &Rc,
                                 vector<double> &GG,
                                 vector<double> &F,
                                 vector<double> &E,
                                 vector<double> &S_hats,
                                 vector<double> &logpvals,
                                 vector<int> &m_i,
                                 const vector<vector<double> > &r,
                                 const vector <int> &cellaffiliations,
                                 const vector<vector<double> > &T_ij,
                                 const double factor){
    int m = 0, whichTumor = 0; // counts column within sample
    double g = 0.0, sum = 0; // sum: A
    int howmanystartinggenes = G[G.size()-1];
    vector<int> startinggenes = {G.begin(), G.begin() + howmanystartinggenes};
    // asserts for the inputs
    for (int i = 0; i < r[0].size(); i++){ // compute column sums for each column (cell) i
        double colsum = 0.0; // this is going to be the sum
        m++;
        for (auto it=startinggenes.begin(); it<startinggenes.end(); it++) // compute sum of ranks for each cell for the current state (= selection of genes)
            colsum += r[*it][i];

        Rc.push_back(colsum); // column sums Rc
        sum += 12*colsum*colsum; // running sum squared (see `A` when reaching the end of the sample)
                                  
        if (cellaffiliations[i] != cellaffiliations[i+1]){ // last column in sample?
            // Remember statistics for this sample
            A.push_back(sum);
            GG.push_back(m/(m-1.0));
            F.push_back(m*(m+1.0));
            E.push_back(3.0*m*(m+1.0)*(m+1.0)); // instead of pow(,2)
            B.push_back(E[whichTumor]*howmanystartinggenes*howmanystartinggenes);
            C.push_back(howmanystartinggenes*F[whichTumor]);
            double tn = 0;
                                                          
            for (int jej = 0; jej < (int) startinggenes.size(); jej++)
                tn += T_ij[G[jej]][whichTumor];
                                                          
            D.push_back((1.0/(m-1.0))*tn - GG[whichTumor]*howmanystartinggenes);
            const double denominator = (C[whichTumor]-D[whichTumor]);
            double s = (denominator==0) ? 0 : (A[whichTumor]-B[whichTumor])/denominator;
            
            double s_hat = s/pow(startinggenes.size(), factor);

#ifndef NDEBUG
            if ((isnan(s_hat)) || (s_hat < 0.0) || (isinf(s_hat))) { // Sanity check
                // FIXME use assert
                s = 0.0;
                cout.precision(std::numeric_limits<double>::digits10 + 1);
                //cout << "Error: randseed = " << randseed << ": S' wasn't good.\n";
            } // making sure our program doens't crash out of some stupid error
#endif //NDEBUG
            assert(not isnan(s_hat) and s_hat > 0.0 and not isinf(s_hat));
            S_hats.push_back(s_hat);
            logpvals.push_back(sqrt_logchisq_R(s_hat, m-1));
            g += logpvals[whichTumor];
            whichTumor++;
            m_i.push_back(m);
            m = 0;
            sum = 0.0;
        } // if cell affiliations
    } // for i in M-1
    
    //assert that the output data has the desired properties
    
    return g;
}

