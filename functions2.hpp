#ifndef FUNCTIONS_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define FUNCTIONS_H

#include <vector>
#include <random>

using std::vector;

double logchisq_R(const double score, const int m_i);

void RemoveGeneFromVector(int GeneToRemove, vector <int> &X);

int RandomlySelectGeneToBan(const vector <int> vect, std::mt19937_64 rng);

vector <int> ConvertAffiliationsToInts(const vector<std::string> cellaffiliations_raw);

void readIn2dData(const char* filename, vector<vector<double> >& ranktable, vector<std::string> & gene_names,  vector<std::string> & cell_affiliations);

void readIn2dData_int(const char* filename, vector<vector<int> >& ranktable);

vector<int> readInAffiliationsVector(const char* filename);

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
);


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
                              const vector<int> &PriorState);

vector <int> SampleBoundaries(const vector <int> cellaffiliations);

void ComputeSquaredRanksAndTieAdjustments(vector<vector <double> > &SSR,
                                          vector<vector <double> > &T_ij,
                                          const vector<vector <double> > &input_table,
                                          const vector <int> &cellaffiliations
                                          );

void PermuteRanksWithinEachGeneInEachSample(
                                            vector<vector<double> > &r,
                                            vector<int> &firstCellInSample,
                                            std::uniform_real_distribution <> &runif,
                                            std::mt19937_64 rng);

double KnownStartingTemperature(const int x);

vector <int> ConstructStateVector(vector <int> Geneset,
                                  const int k, vector <int> banned_genes);

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
                                 const double factor);

#endif
