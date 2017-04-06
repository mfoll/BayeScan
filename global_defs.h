//    This program, BayeScan, aims at detecting genetics markers under selection,
//	  based on allele frequency differences between population. 
//    Copyright (C) 2010  Matthieu Foll
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

/* this is the header file that contains constants and global variable
definitions for all the functions used in the mixture program.

Each file containing functions should have a
#include "global_defs.h"
statement.

Conditional Preprocessing. In all files except the one in which the
function main() resides, the preprocessor replaces GLOB with extern.
The file containing main should start with
#define UN_EXTERN
*/

#ifndef UN_EXTERN
#define GLOB extern
#else
#define GLOB
#endif

#include "structure.h"

#include "MersenneTwister.h"

#include <fstream>
using namespace std;

//#define sd_prior_alpha 1 // sd of the normal prior for alpha

//#define SCHED_G schedule(dynamic,1) // used for loops over groups (g)
#define SCHED_I schedule(guided,10) // used for loops over loci (i)
#define SCHED_J schedule(guided,1)  // used for loops over populations (j)

//#define SCHED_G schedule(static) // used for loops over groups (g)
//#define SCHED_I schedule(static) // used for loops over loci (i)
//#define SCHED_J schedule(static)  // used for loops over populations (j)


//////////////////////////////////////////
//  global variables declarations
//////////////////////////////////////////
GLOB int I; // nb of loci
GLOB int J; // nb of populations

GLOB int num_threads;

GLOB float codominant; // selection2 indicates the type of markers 0 for dominant, 1 for codominant, 0.5 for AFLP band intensity

GLOB bool SNP_genotypes; // boolean indicating that input file is SNP genotypes (0=aa, 1=Aa, 2=AA) as for AFLPs

GLOB bool fstat; // boolean indicating if we only calculate F-stat (no selection)
GLOB bool all_trace; // boolean indicating if we output the alpha parameters in the MCMC trace

GLOB double *alpha; // vector of locus effect
GLOB bool *alpha_included; // indicates if alpha is not null
//GLOB int nb_alpha_included; // number of non null alpha values

GLOB double prior_odds; // prior odds for neutral model 

GLOB double *mu; // vector of mean intensity for Aa  selection2
GLOB double *delta; // mu+delta = vector of mean intensity for AA selection2
GLOB double *sigma1; // vector of s.d. for intensity for Aa
GLOB double *sigma2; // vector of s.d. for intensity for AA

GLOB float abscence_pc; // fraction of the maximum aflp intensity value to use as a threshold for aa genotype

GLOB double *mean_alpha; // mean and variance of alpha proposal (estimated with pilot run)
GLOB double *var_alpha;

GLOB double *post_alpha; // posterior mean of alpha
GLOB int    *nb_alpha;   // nb of time alpha was included in the model
//GLOB double *p_value; // empirical p-value like in Beaumont and Balding 2004

//GLOB int alpha_updates;

GLOB double *beta; // vector of population effect

// prior for Fis
GLOB bool prior_fis_unif; // true if uniform prior, false if beta
GLOB double prior_fis_lb;  // lower bound for uniform prior
GLOB double prior_fis_hb;  // higher bound
GLOB double prior_fis_a; // a for beta
GLOB double prior_fis_b;  // b for beta

GLOB int nr_out;
GLOB int interval;
GLOB int discard;
GLOB int tot_nr_of_iter;

GLOB int nb_pilot;
GLOB int pilot_length;

GLOB double *cur_fst; // corresponding fst estimation
GLOB double *post_fst; // corresponding fst estimation

//GLOB double **fst; // matrix of estimated mean fst values useless ? selection2

GLOB double *f; // within population inbreeding coefficient

GLOB double *freq_ancestral; // vector of allele frequencies of the ancestral population

GLOB bool *discarded_loci; // vector of boolean indicating if loci must be discarded (=true) or not (=false)

GLOB double *max_intenstiy; // maximum intensity for each locus for normalisation selection2

GLOB double a_p; // parameters of beta prior for allele frequencies of the ancestral population
//GLOB double b_p;

GLOB pop_data *pop; // populations allele count
GLOB allele_freq *freq_locus; // allele frequences

GLOB double log_likelihood; // loglikelihood

GLOB double *locus_likelihood; // partial part of locus likelihood often used

//acceptance rates selection2
GLOB double *acc_alpha;
GLOB double *acc_beta;
GLOB double *acc_freq_ancestral;
GLOB double *acc_f;
GLOB double **acc_freq;
GLOB double acc_a_p;
GLOB double *acc_mu;
GLOB double *acc_delta;
GLOB double *acc_sigma1;
GLOB double *acc_sigma2;

GLOB MTRand randgen; // random generator
GLOB MTRand *randgen_parallel; // random generator

GLOB int iter; // iteration index


GLOB double  *e_ancestral; // value p' is proposed between p-e and p+e (and stay in (0,1))
GLOB double  **e_freq; // same for allele freq

GLOB double  *e_f; // same for f selection2
GLOB double  *e_mu; // same for mu selection2
GLOB double  *e_delta; // same for delta selection2

GLOB double  *var_prop_alpha; // sd of the normal proposal for alpha
GLOB double  *var_prop_beta; // sd of the normal proposal for beta

GLOB double  *var_prop_sigma1; // sd of the normal proposal for sigma selection2
GLOB double  *var_prop_sigma2; // sd of the normal proposal for sigma selection2

GLOB double  var_prop_a_p; // sd of normal proposal for a_p

GLOB double  m1_prior_alpha; // prior for alpha
GLOB double  m2_prior_alpha;
GLOB double  sd_prior_alpha;

//////////////////////////////////////////
//  global functions declarations
//////////////////////////////////////////
GLOB int read_input(std::ifstream& infile,string outcheck);
GLOB int read_input_intensity(string outname,string outcheck); // selection2
GLOB int read_discarded(std::ifstream& infile);
GLOB void write_output(std::ofstream& outfile);
GLOB void write_output_intensity(std::ofstream& outfile); // selection2
GLOB void write_anc_freq(std::ofstream& outfile);
GLOB void write_freq(std::ofstream outfiles[]);
GLOB void write_output_genotype(std::ofstream outfile[]);
GLOB void write_freq_matrix(std::ofstream& freq_pop,int cur_out);

GLOB double genbet(double aa,double bb);

GLOB double gammaln(double xx);
GLOB double factln(int n);

GLOB void dirichlet_dev(double vectX[],double vectA[],int arraySize);
GLOB double rstgam(double alpha1);

GLOB double allelecount_loglikelihood();
GLOB double log_prior_alpha(double alpha);
GLOB double intensity_loglikelihood(float y,int j,int i); // selection2

GLOB void update_ancestral_freq();
GLOB void update_a_p();
GLOB void update_freq();
GLOB void update_beta();
GLOB void update_alpha();

GLOB void update_f_random();

GLOB void update_freq_codominant();
GLOB void update_alpha_codominant();
GLOB void update_beta_codominant();

GLOB void update_f_intensity(); // selection2
GLOB void update_freq_intensity();
GLOB void update_mu_intensity();
GLOB void update_delta_intensity();
GLOB void update_sigma1_intensity();
GLOB void update_sigma2_intensity();

GLOB void jump_model();
GLOB void jump_model_codominant();
