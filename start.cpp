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

#include <time.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <math.h>

#include <omp.h>

#include "anyoption.h"

#define UN_EXTERN
#include "global_defs.h"

#include "errors.cpp"

using namespace std;

/* 1. CREATE AN OBJECT */
AnyOption *opt = new AnyOption();

//////////////////////////////////////
//  main program
//////////////////////////////////////
int main(int argc, char **argv)
{

    int last_printf=0;
    int acc_rate;
    int cur_out=0;

    nr_out=5000;
    interval=10;
    discard=50000;

    nb_pilot=20;
    pilot_length=5000;

    prior_odds=10;

    string file_name;
    string odir="";
    string file_prefix;

    /* 2. SET PREFERENCES  */
    opt->noPOSIX(); /* do not check for POSIX style character options */

    /* 3. SET THE USAGE/HELP   */
	opt->addUsage( " --------------------------- " );
	opt->addUsage( " | BayeScan 2.1 usage:     | " );
	opt->addUsage( " --------------------------- " );
    opt->addUsage( " -help        Prints this help " );	
	opt->addUsage( " --------------------------- " );
	opt->addUsage( " | Input                   | " );
	opt->addUsage( " --------------------------- " );		
    opt->addUsage( " alleles.txt  Name of the genotypes data input file " );
    opt->addUsage( " -d discarded Optional input file containing list of loci to discard" );
    opt->addUsage( " -snp         Use SNP genotypes matrix" );
	opt->addUsage( " --------------------------- " );
	opt->addUsage( " | Output                  | " );
	opt->addUsage( " --------------------------- " );	
    opt->addUsage( " -od .        Output file directory, default is the same as program file" );
    opt->addUsage( " -o alleles   Output file prefix, default is input file without the extension" );
    opt->addUsage( " -fstat       Only estimate F-stats (no selection)" );
    opt->addUsage( " -all_trace   Write out MCMC trace also for alpha paremeters (can be a very large file)" );	
	opt->addUsage( " --------------------------- " );
	opt->addUsage( " | Parameters of the chain | " );
	opt->addUsage( " --------------------------- " );	
	opt->addUsage( " -threads n   Number of threads used, default is number of cpu available " );
    opt->addUsage( " -n 5000      Number of outputted iterations, default is 5000 " );
    opt->addUsage( " -thin 10     Thinning interval size, default is 10 " );
    opt->addUsage( " -nbp 20      Number of pilot runs, default is 20 " );
    opt->addUsage( " -pilot 5000  Length of pilot runs, default is 5000 " );
    opt->addUsage( " -burn 50000  Burn-in length, default is 50000 " );
	opt->addUsage( " --------------------------- " );
	opt->addUsage( " | Parameters of the model | " );
	opt->addUsage( " --------------------------- " );	
    opt->addUsage( " -pr_odds 10  Prior odds for the neutral model, default is 10 " );
    opt->addUsage( " -lb_fis 0    Lower bound for uniform prior on Fis (dominant data), default is 0" );
    opt->addUsage( " -hb_fis 1    Higher bound for uniform prior on Fis (dominant data), default is 1" );
    opt->addUsage( " -beta_fis    Optional beta prior for Fis (dominant data, m_fis and sd_fis need to be set)");
    opt->addUsage( " -m_fis 0.05  Optional mean for beta prior on Fis (dominant data with -beta_fis)" );
    opt->addUsage( " -sd_fis 0.01 Optional std. deviation for beta prior on Fis (dominant data with -beta_fis)" );
    opt->addUsage( " -aflp_pc 0.1 Threshold for the recessive genotype as a fraction of maximum band intensity, default is 0.1" );
	opt->addUsage( " --------------------------- " );
	opt->addUsage( " | Output files            | " );
	opt->addUsage( " --------------------------- " );	
    opt->addUsage( " -out_pilot   Optional output file for pilot runs" );
    opt->addUsage( " -out_freq    Optional output file for allele frequencies" );

    /* 4. SET THE OPTION STRINGS/CHARACTERS */
    opt->setFlag(  "help");
    opt->setOption(  "n" );
    opt->setOption(  "thin" );
    opt->setOption(  "burn" );
    opt->setOption(  "pr_odds" );
    opt->setOption(  "o" );
    opt->setOption(  "d" );
    opt->setOption(  "threads" );
    opt->setFlag(  "out_pilot");
    opt->setFlag(  "out_freq");
    opt->setFlag(  "beta_fis");
    opt->setFlag(  "fstat");
	opt->setFlag(  "all_trace");	
    opt->setFlag(  "snp");
    opt->setOption(  "lb_fis" );
    opt->setOption(  "hb_fis" );
    opt->setOption(  "m_fis" );
    opt->setOption(  "sd_fis" );
    opt->setOption(  "aflp_pc" );
    opt->setOption(  "od" );
    opt->setOption(  "nbp" );
    opt->setOption(  "pilot" );

    /* 5. PROCESS THE COMMANDLINE AND RESOURCE FILE */
    /* go through the command line and get the options  */
    opt->processCommandArgs( argc, argv );

    if ( ! opt->hasOptions() || opt->getArgc()==0)  /* print usage if no options */
    {
        opt->printUsage();
        delete opt;
        exit(0);
    }

    /* 6. GET THE VALUES */

    // create and seed parallels randgens
    num_threads=omp_get_max_threads();
    if ( opt->getValue( "threads" ) != NULL )
        num_threads=max(1,atoi(opt->getValue( "threads" )));
    cout << "Using " << num_threads << " threads (" << omp_get_max_threads() << " cpu detected on this machine)" << endl;
    omp_set_num_threads(num_threads);
    randgen_parallel=new MTRand[num_threads];
    for (int k=0;k<num_threads;k++)
        randgen_parallel[k].seed(randgen.randInt());

    if ( opt->getFlag( "help" ) || opt->getFlag( 'h' ) )
        opt->printUsage();

    if ( opt->getValue( "n" ) != NULL )
        nr_out=atoi(opt->getValue( "n" ));

    if ( opt->getValue( "thin" ) != NULL )
        interval=atoi(opt->getValue( "thin" ));

    if ( opt->getValue( "burn" ) != NULL )
        discard=atoi(opt->getValue( "burn" ));

    if ( opt->getValue( "nbp" ) != NULL )
        nb_pilot=atoi(opt->getValue( "nbp" ));

    if ( opt->getValue( "pilot" ) != NULL )
        pilot_length=atoi(opt->getValue( "pilot" ));

    if ( opt->getValue( "pr_odds" ) != NULL )
        prior_odds=atof(opt->getValue( "pr_odds" ));

    file_name=opt->getArgv(0);

    if (opt->getFlag( "snp" ) )
        SNP_genotypes=true;
    else SNP_genotypes=false;

    if (opt->getFlag( "fstat" ) )
        fstat=true;
    else fstat=false;
	
	if (opt->getFlag( "all_trace" ) )
        all_trace=true;
    else all_trace=false;

    if ( opt->getValue( "od" ) != NULL ) {
        odir=opt->getValue( "od" );
        odir=odir+"/";
    }

    if ( opt->getValue( "o" ) != NULL )
        file_prefix=opt->getValue( "o" );
    else
    {
        file_prefix=file_name;
        if (file_prefix.find(".",1)!= string::npos)
            file_prefix=file_prefix.substr(0,file_prefix.length()-4);
        int position;
        if (position=file_prefix.find_last_of('/'))
            file_prefix=file_prefix.substr(position+1,file_prefix.length());
        //cout << file_prefix << endl;

    }

    // prior for Fis
    double prior_fis_mean;
    double prior_fis_sd;
    prior_fis_unif=true; // true if uniform prior, false if beta
    if (opt->getFlag( "beta_fis" ) )
    {
        prior_fis_unif=false;
        if ( opt->getValue( "m_fis" ) != NULL && opt->getValue( "sd_fis" ) != NULL)
        {
            prior_fis_mean=atof(opt->getValue( "m_fis" )); // mean for beta
            prior_fis_sd=atof(opt->getValue( "sd_fis" )); // sd for beta
            prior_fis_a=prior_fis_mean*(prior_fis_mean*(1-prior_fis_mean)/(prior_fis_sd*prior_fis_sd)-1);
            prior_fis_b=(1-prior_fis_mean)*(prior_fis_mean*(1-prior_fis_mean)/(prior_fis_sd*prior_fis_sd)-1);
            if (prior_fis_a<=0 || prior_fis_b<=0)
            {
                cout << "Error in prior definition for Fis";
                return 1;
            }
        }
        else
        {
            cout << "Error in prior definition for Fis";
            return 1;
        }

    }
    else
    {
        if ( opt->getValue( "lb_fis" ) != NULL && opt->getValue( "hb_fis" ) != NULL)
        {
            prior_fis_lb=atof(opt->getValue( "lb_fis" ));  // lower bound for uniform prior
            prior_fis_hb=atof(opt->getValue( "hb_fis" ));  // higher bound

            if (prior_fis_unif && (prior_fis_lb<0 || prior_fis_hb>1 || prior_fis_lb>=prior_fis_hb))
            {
                cout << "Error in prior definition for Fis";
                return 1;
            }
        }
        else
        {
            prior_fis_lb=0;
            prior_fis_hb=1;
        }

    }

    abscence_pc=0.1;
    if ( opt->getValue( "aflp_pc" ) != NULL )
    {
        abscence_pc=atof(opt->getValue( "aflp_pc" ));
    }

// Calculate length of run based on the Parameters of Markov Chain
    tot_nr_of_iter = nr_out*interval + discard; //number of iterations of the algorithm

// create main output file
    ofstream mainout((odir+file_prefix+".sel").c_str());
    assure(mainout);

// read input file
    ifstream infile(file_name.c_str());
    assure(infile);

// determine if data are dominant or codominant markers
    string line;
    int k=0;
    while (getline(infile,line,'=') )
    {
        if (line.length()>=5 && line.substr(line.length()-5,line.length())=="[pop]") //read allele count
        {
            getline(infile,line);
            getline(infile,line); // read the line of allele count at locus i
            istringstream read_line(line);
            int trash ;
            /*while(read_line.good())
                    {
                    read_line >> trash;
                    k++;
            		}*/
            while (read_line >> trash)
                k++;
            break;
        }
    }

// selection2
    if (k==0)
        codominant=0.5;
    else
        codominant=(float)(k>=4);

    infile.close();
    infile.clear();

    if (codominant==0.5)   // selection2
    {
        string verif_name=odir+file_prefix+"_Verif.txt";
        if (read_input_intensity(file_name,verif_name))
        {
            return 1;
        }
    }
    else
    {
        infile.open(file_name.c_str());
        string verif_name=odir+file_prefix+"_Verif.txt";
        if (read_input(infile,verif_name))
        {
            return 1;
        }
        infile.close();
        infile.clear();
    }

//read discarded input file
    for (int i=0; i < I; i++)
    {
        discarded_loci[i]=false;
    }

    if ( opt->getValue( "d" ) != NULL )
    {
        // read input file
        ifstream infile_discarded(opt->getValue( "d" ));
        assure(infile_discarded);
        read_discarded(infile_discarded);
        infile_discarded.close();
    }

   // public version: remove frequency output, keep only posterior mean matrix
   /* ofstream ancfreq;
// create ouput file for ancestral freq
        if (opt->getFlag( "out_freq" ) )
        {
            ancfreq.open((odir+file_prefix+"_ancestral.txt").c_str());
            assure(ancfreq);
            // write header
            if (codominant<1)
            {
                for (int i=0;i<I;i++)
                {
                    if (!discarded_loci[i])
                        ancfreq << "locus" << (i+1) << " " ;
                }
                ancfreq << endl;
            }
    }*/

    ofstream fst_file;
// create ouput file for fst mean values
    if (!opt->getFlag( "fstat" ))
    {
        fst_file.open((odir+file_prefix+"_fst.txt").c_str());
        assure(fst_file);
    }
// public version: remove frequency output, keep only posterior mean matrix
// create output files for allele frequencies in each population
  //  ofstream *freq_pop2;
    ofstream freq_pop;
    if (opt->getFlag( "out_freq" ) && codominant<1) // selection2
    {
        freq_pop.open((odir+file_prefix+"_freq.txt").c_str());
        /*try
        {
            freq_pop2=new ofstream[J];
        }
        catch (const std::exception & Exp )
        {
            cout << "Not enough memory for freq_pop2" << endl;
            return 1;
        }
        for (int j=0;j<J;j++) //cycle over populations
        {
            ostringstream oss;
            oss << odir << file_prefix << "_pop" << (j+1) << ".txt";
            freq_pop2[j].open(oss.str().c_str());

            if (freq_pop2[j].fail())
            {
                cout << "Cannot not open file " << oss.str().c_str();
                return 1;
            }
            else
            {
                // write header
                if (codominant<1)
                {
                    for (int i=0;i<I;i++)
                    {
                        if (!discarded_loci[i])
                            freq_pop2[j] << "locus" << (i+1) << " " ;
                    }
                    freq_pop2[j] << endl;
                }
            }
        }*/
    }

// create ouput file for acceptance rate
    ofstream outAccRte;
    outAccRte.open((odir+file_prefix+"_AccRte.txt").c_str());
    assure(outAccRte);

    acc_a_p=0;

    if (codominant==0.5)
    {
        for (int i=0; i < I; i++)   // selection2
        {
            acc_mu[i]=0;
            acc_delta[i]=0;
            acc_sigma1[i]=0;
            acc_sigma2[i]=0;
        }
    }

    for (int i=0; i < I; i++)   // selection2
    {
        acc_alpha[i]=0;
        acc_freq_ancestral[i]=0;
    }
    for (int j=0; j < J; j++)   // selection2
    {
        acc_beta[j]=0;
        if (codominant<1)
            acc_f[j]=0;
    }
    for (int i=0; i < I; i++)   // selection2
    {
        for (int j=0; j < J; j++)   // selection2
        {
            acc_freq[i][j]=0;
        }
    }

    if (codominant==1)
        outAccRte <</* "alpha " <<*/ "beta  " << "ances "/* << "f      "* << "a_p  " << "b_p  "*/ <<endl ;
    else if (codominant==0)
        outAccRte <</* "alpha " <<*/ "beta  " << "ances "/* << "f      "*/ << "freq   " /* << "a_p   " << "b_p  "*/ <<endl ;
    else // selection2
        outAccRte <</* "alpha " <<*/ "beta  " << "ances "/* << "f      "*/ << "freq   " /* << "a_p   " << "b_p  "*/ << "f     " << "mu    " << "delta " << "sigma1 " << "sigma2 " <<endl ;

// create output file for proposal variance
    ofstream outprop;
    if (opt->getFlag( "out_pilot" ))
    {
        outprop.open((odir+file_prefix+"_prop"+".txt").c_str());
        assure(outprop);
    }

// initialize f from uniform distribution  selection2
    if (codominant<1)
    {
        if (prior_fis_unif)
        {
            for (int j=0;j<J;j++)
                f[j]=(prior_fis_lb+prior_fis_hb)/2;
        }
        else
        {
            for (int j=0;j<J;j++)
                f[j]=prior_fis_mean;
        }
    }

// initialize alpha with normal(0,1) values
    for (int i=0;i<I;i++)
    {
        alpha[i]=0;//randgen.randNorm(0,1);
        alpha_included[i]=true;
        post_alpha[i]=0;
        post_fst[i]=0;
        nb_alpha[i]=0;
    }
    //alpha_updates=0;
//    nb_alpha_included=0;

// initialize beta with fst=0.1 values
    for (int j=0;j<J;j++)
        beta[j]=log(0.1/(1-0.1));

// initialize beta prior for ancestral allele frequencies
    a_p=1;

    if (codominant==1) // selection2
    {
// initialize allele frequences from non informative dirichlet distribution with parameter 1
        double *alphaP;
        for (int i=0;i<I;i++) //cycle over loci within populations
        {
            alphaP = new double[freq_locus[i].ar];
            for (int k=0;k<freq_locus[i].ar;k++)
            {
                alphaP[k]=1.0;
                for (int j=0;j<J;j++)
                    alphaP[k]+=pop[j].locus[i].data_allele_count[k];
            }
            dirichlet_dev(freq_locus[i].allele,alphaP,freq_locus[i].ar);
            delete[] alphaP;
        }
    }
    else
    {
// initialize ancestral allele frequencies from beta distribution
        for (int i=0;i<I;i++)
            freq_ancestral[i]=randgen.randDblExc();

// initialize allele frequencies from beta(theta*pi,theta*(1-pi))
        double theta;
        for (int i=0;i<I;i++)
        {
            for (int j=0;j<J;j++)
            {
                theta=exp(-(alpha[i]+beta[j]));
                pop[j].locus[i].p=genbet(theta*freq_ancestral[i],theta*(1-freq_ancestral[i])); // selection2 bug : removed a_p*...
                if (pop[j].locus[i].p<=0.0001)
                    pop[j].locus[i].p=0.0001;
                if (pop[j].locus[i].p>=0.9999)
                    pop[j].locus[i].p=0.9999;
                pop[j].locus[i].mean_p=0;
            }
        }
    }

// initialize mu and sigma (to do) selection2
    if (codominant==0.5)
    {
        for (int i = 0; i < I; i++)
        {
            mu[i]=0.2;
            delta[i]=0.2;
            sigma1[i]=0.05;
            sigma2[i]=0.05;
        }
    }

    log_likelihood=allelecount_loglikelihood();

//interval between output of acceptance rates
    acc_rate = tot_nr_of_iter/50;

////////////////////////////////////////////////////////
// Pilot runs for within model updates acceptance rates
// and reversible jump efficient proposal
////////////////////////////////////////////////////////
//Form1->Time->Caption="";

    cout << "Pilot runs..." << endl;

    m1_prior_alpha=0;
    m2_prior_alpha=0;
    sd_prior_alpha=1;

    for (int i=0; i < I; i++)
    {
        e_ancestral[i]=0.2;
        var_prop_alpha[i]=1.0;
    }
    if (codominant==0.5)
    {
        for (int i=0; i < I; i++)
        {
            e_mu[i]=0.01; // selection2 -> to do : add as parameters?
            e_delta[i]=0.01; // selection2 -> to do : add as parameters?
            var_prop_sigma1[i]=0.3; // selection2 -> to do : add as parameters?
            var_prop_sigma2[i]=0.3; // selection2 -> to do : add as parameters?
        }
    }

    for (int j=0; j < J; j++)
    {
        if (codominant<1)
            e_f[j]=0.05;  // selection2 -> to do : add as parameters?
        var_prop_beta[j]=0.7;
    }

    for (int i=0; i < I; i++)
    {
        for (int j=0; j < J; j++)
        {
            e_freq[i][j]=0.2;
        }
    }

    var_prop_a_p=0.2;

    for (int i=0;i<I;i++)
    {
        mean_alpha[i]=0.0;
        var_alpha[i]=5.0;
    }

    int nb_pilot_alpha=0;


    double *m2=new double[I];
    for (int i=0;i<I;i++)
        m2[i]=0;


    for (int k=0;k<nb_pilot;k++)
    {
// pilot runs
        for (iter=0;iter<pilot_length;++iter)
        {
            if  ((100*(pilot_length*k+iter))/(pilot_length*nb_pilot)>last_printf)
            {
                last_printf=(100*(pilot_length*k+iter))/(pilot_length*nb_pilot);
                cout <<  last_printf << "% ";
                cout.flush();
            }
            // do not block application during the long process
            if (codominant<1) //selection2
                update_beta();
            else
                update_beta_codominant();

            if (!fstat)   // update selection related parameters
            {
                /*for (int i=0;i<I;i++)
                {
                    if (!discarded_loci[i])
                    {
                        if (codominant<1) //selection2
                            update_alpha_i(i);
                        else
                            update_alpha_i_codominant(i);
                    }
                }*/
                if (codominant<1) //selection2
                    update_alpha();
                else
                    update_alpha_codominant();
            }

            // update ancestral allele frequences
            if (codominant<1) // selection2
            {
                update_ancestral_freq();
                // update beta prior for ancestral allele frequences
                //update_a_p(Application);
                a_p=1;
            }
            // update f
            //update_f();
            if (codominant==0) // selection2
                update_f_random();
            else if (codominant==0.5) update_f_intensity();

            if (codominant==0.5 && !SNP_genotypes) // selection2
            {
                update_mu_intensity();
                update_delta_intensity();
                update_sigma1_intensity();
                update_sigma2_intensity();
            }

            // update allele frequences
            if (codominant==0) // selection2
                update_freq();
            else if (codominant==1)
                update_freq_codominant();
            else
                update_freq_intensity();


            //calculate mean and variance for alpha
            if (2*(k+1)>=nb_pilot)
            {
                for (int i=0;i<I;i++)
                {
                    if (!discarded_loci[i])
                    {
                        mean_alpha[i]+=alpha[i];
                        m2[i]+=alpha[i]*alpha[i];
                    }
                }
            }

        } // end of one pilot run

        if (2*(k+1)>=nb_pilot)
            nb_pilot_alpha++;

// check acc_rates
        for (int i=0; i < I; i++)
        {
            if (!discarded_loci[i])
            {
                // for alpha
                if (acc_alpha[i]/(pilot_length)>0.4)
                    var_prop_alpha[i]*=1.2;
                if (acc_alpha[i]/(pilot_length)<0.2)
                    var_prop_alpha[i]/=1.2;

                // for ancestral allele frequences
                if (acc_freq_ancestral[i]/(pilot_length)>0.4)
                    e_ancestral[i]*=1.1;
                if (acc_freq_ancestral[i]/(pilot_length)<0.2)
                    e_ancestral[i]/=1.1;

                if (codominant==0.5) // selection2
                {
                    // for mu
                    if (acc_mu[i]/(pilot_length)>0.2)
                        e_mu[i]*=1.1;
                    if (acc_mu[i]/(pilot_length)<0.1)
                        e_mu[i]/=1.1;
                    // for delta
                    if (acc_delta[i]/(pilot_length)>0.2)
                        e_delta[i]*=1.1;
                    if (acc_delta[i]/(pilot_length)<0.1)
                        e_delta[i]/=1.1;
                    // for sigma1
                    if (acc_sigma1[i]/(pilot_length)>0.4)
                        var_prop_sigma1[i]*=1.2;
                    if (acc_sigma1[i]/(pilot_length)<0.2)
                        var_prop_sigma1[i]/=1.2;
                    // for sigma2
                    if (acc_sigma2[i]/(pilot_length)>0.4)
                        var_prop_sigma2[i]*=1.2;
                    if (acc_sigma2[i]/(pilot_length)<0.2)
                        var_prop_sigma2[i]/=1.2;
                }
            }
        }

        for (int j=0; j < J; j++)
        {
            // for beta
            if (acc_beta[j]/(pilot_length)>0.4)
                var_prop_beta[j]*=1.2;
            if (acc_beta[j]/(pilot_length)<0.2)
                var_prop_beta[j]/=1.2;

            if (codominant<1) // selection2
            {
                // for f
                if (acc_f[j]/(pilot_length)>0.4)
                    e_f[j]*=1.1;
                if (acc_f[j]/(pilot_length)<0.2)
                    e_f[j]/=1.1;
            }

        }

        if (codominant<1) // selection2
        {
            for (int j=0; j < J; j++)
            {
                for (int i=0; i < I; i++)
                {
                    if (!discarded_loci[i])
                    {
                        // for allele frequences
                        if (acc_freq[i][j]/(pilot_length)>0.5)
                            e_freq[i][j]*=1.2;
                        if (acc_freq[i][j]/(pilot_length)<0.3)
                            e_freq[i][j]/=1.2;
                    }
                }
            }
            // for a
            if (acc_a_p/pilot_length>0.4)
                var_prop_a_p*=1.2;
            if (acc_a_p/pilot_length<0.2)
                var_prop_a_p/=1.2;
        }



// re-initialize acc rates

        acc_a_p=0;

        if (codominant==0.5)
        {
            for (int i=0; i < I; i++)   // selection2
            {
                acc_mu[i]=0;
                acc_delta[i]=0;
                acc_sigma1[i]=0;
                acc_sigma2[i]=0;
            }
        }

        for (int i=0; i < I; i++)   // selection2
        {
            acc_alpha[i]=0;
            acc_freq_ancestral[i]=0;
        }
        for (int j=0; j < J; j++)   // selection2
        {
            acc_beta[j]=0;
            if (codominant<1)
                acc_f[j]=0;
        }
        for (int i=0; i < I; i++)   // selection2
        {
            for (int j=0; j < J; j++)   // selection2
            {
                acc_freq[i][j]=0;
            }
        }

    } // end of pilot runs

// calculate final mean and variance for alpha proposal in RJ
    if (nb_pilot>0)
    {
        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
            {
                mean_alpha[i]/=(nb_pilot_alpha*pilot_length);
                var_alpha[i]=m2[i]/(nb_pilot_alpha*pilot_length)-mean_alpha[i]*mean_alpha[i];
            }
        }
    }

    if (opt->getFlag( "out_pilot" ))
    {
        outprop << endl;
// write the results of the pilot runs in outprop file
        outprop << setw(3) << setprecision(3) << setiosflags(ios::showpoint) << "Variance of alpha proposal: " << var_prop_alpha[0] << "\n";
        outprop << setw(3) << setprecision(3) << setiosflags(ios::showpoint) << "Variance of beta proposal: " << var_prop_beta[0] << "\n";
// selection2
        if (codominant<1) outprop << setw(3) << setprecision(3) << setiosflags(ios::showpoint) << "Variance of a proposal: " << var_prop_a_p << "\n";
        if (codominant<1) outprop << setw(3) << setprecision(3) << setiosflags(ios::showpoint) << "Variance of allele freq proposal: " << e_freq[0][0] << "\n";
        if (codominant==0.5) outprop << setw(3) << setprecision(3) << setiosflags(ios::showpoint) << "Variance of mu proposal: " << e_mu[0] << "\n";
        if (codominant==0.5) outprop << setw(3) << setprecision(3) << setiosflags(ios::showpoint) << "Variance of delta proposal: " << e_delta[0] << "\n";
        if (codominant==0.5) outprop << setw(3) << setprecision(3) << setiosflags(ios::showpoint) << "Variance of sigma1 proposal: " << var_prop_sigma1[0] << "\n";
        if (codominant==0.5) outprop << setw(3) << setprecision(3) << setiosflags(ios::showpoint) << "Variance of sigma2 proposal: " << var_prop_sigma2[0] << "\n";
        outprop << setw(3) << setprecision(3) << setiosflags(ios::showpoint) << "Variance of ancestral allele freq proposal: " << e_ancestral[0] << "\n";

        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
                outprop << setw(3) << setprecision(3) << setiosflags(ios::showpoint) << "Mean and variance alpha_" << (i+1) << " RJ proposal: " << mean_alpha[i] << " , " << var_alpha[i] << "\n";
        }
        outprop.close();
    }

    // write header
    mainout << "logL ";
    if (codominant<1)
    {
        //mainout << "a  ";
        for (int j=0;j<J;j++)
        {
            mainout << "Fis" << (j+1) << " " ;
        }
    }
    for (int j=0;j<J;j++)
    {
        mainout << "Fst" << (j+1) << " " ;
    }

    if (!fstat && all_trace)
    {
        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
                mainout << "alpha" << (i+1) << " " ;
        }
    }
    mainout << endl;

    /*    ofstream out_intensity; // selection2
    if (codominant==0.5 && !SNP_genotypes)
    {
        // create ouput file
        out_intensity.open((odir+file_prefix+"_intensity.txt").c_str());
        assure(out_intensity);
        // write header
        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
                out_intensity << "muAa_" << (i+1) << " "  << "muAA_" << (i+1) << " ";
        }
        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
                out_intensity << "sigmaAa_" << (i+1) << " "  << "sigmaAA_" << (i+1) << " ";
        }
        out_intensity << endl;
    }*/

    /*ofstream *out_genotype; // selection2
    if (codominant==0.5 && opt->getFlag( "out_freq" ))
    {
        try
        {
            out_genotype=new ofstream[J];
        }
        catch (const std::exception & Exp )
        {
            cout << "Not enough memory for out_genotype" << endl;
            return 1;
        }
        for (int j=0;j<J;j++) //cycle over populations
        {
            ostringstream oss;
            oss << odir << file_prefix << "_genotypes" << (j+1) << ".txt";
            out_genotype[j].open(oss.str().c_str());
            if (out_genotype[j].fail())
            {
                cout << "Cannot not open genotype file for pop"<<(j+1)<< endl;
                return 1;
            }
            else
            {
                // write header
                for (int i=0;i<I;i++)
                {
                    if (!discarded_loci[i])
                        out_genotype[j] << "Aa_" << (i+1) << " " << "AA_" << (i+1) << " " ;
                }
                out_genotype[j] << endl;
            }
        }
    }*/

// total nb of alpha updates to estimate acc rate with RJ
    //alpha_updates=0;

// initialize alpha
    for (int i=0;i<I;i++)
    {
        alpha[i]=0;  // commented this line selection2
        alpha_included[i]=false; // selection2 uncommented this line -> makes more sense, no ?
    }

/////////////////////////////////
// MH algorithm loop
////////////////////////////////

    last_printf=0;
    cout << endl << "Calculation..." << endl;

    for (iter=0;iter<tot_nr_of_iter;++iter)
    {

// print the gauge
        if  ((100*iter/tot_nr_of_iter)>last_printf)
        {
            last_printf=100*iter/tot_nr_of_iter;
            cout << 100*iter/tot_nr_of_iter << "% " ;
            cout.flush();
        }

// do not block application during the long process

         if (codominant<1)  // selection2
             update_beta();
         else
             update_beta_codominant();

        if (!fstat)   // update selection related parameters
        {
            /*for (int i=0;i<I;i++)
            {
                if (!discarded_loci[i])
                {
                    if (alpha_included[i])
                    {
                        if (codominant<1) // selection2
                            update_alpha_i(i);
                        else
                            update_alpha_i_codominant(i);
                        alpha_updates++;
                    }
                }
            }*/
             if (codominant<1) //selection2
                    update_alpha();
                else
                    update_alpha_codominant();
        }

        if (codominant<1) // selection2
        {
            // update ancestral allele frequences
            update_ancestral_freq();
            // update beta prior for ancestral allele frequences
            a_p=1;
            //update_a_p(Application);
        }

// update f selection2
        if (codominant==0) // selection2
            update_f_random();
        else if (codominant==0.5) update_f_intensity();

        if (codominant==0.5 && !SNP_genotypes) // selection2
        {
            update_mu_intensity();
            update_delta_intensity();
            update_sigma1_intensity();
            update_sigma2_intensity();
        }

//update_f(Application);
//f=randgen.randDblExc();
//log_likelihood=allelecount_loglikelihood(Application);
// update allele frequences
        if (codominant==0) // selection2
            update_freq();
        else if (codominant==1)
            update_freq_codominant();
        else
            update_freq_intensity();

        if (!fstat)   // update selection related parameters
        {
            /*for (int i=0;i<I;i++)
            {
                if (!discarded_loci[i])
                {*/
                    if (codominant<1) // selection2
                        jump_model();
                    else
                        jump_model_codominant();
             /*   }
            }*/

        }

        if (!(iter%1000))
            log_likelihood=allelecount_loglikelihood();

// write out parameters for this step if necessary
        if (((iter/interval)*interval == iter)&&(iter>discard))
        {
            cur_out++;
// calculate corresponding fst value
            for (int i=0;i<I;i++)
            {
                if (!discarded_loci[i])
                {
                    cur_fst[i]=0;  // selection2
                    for (int j=0;j<J;j++)
                        cur_fst[i]+=1/(1+exp(-(alpha[i]+beta[j])));
                    cur_fst[i]/=J;
                    post_alpha[i]+=alpha[i];
                    post_fst[i]+=cur_fst[i];
                    if (alpha_included[i])
                    {
                        nb_alpha[i]++;
                    }
                    if (codominant<1)
                    {
                        for (int j=0;j<J;j++)
                            pop[j].locus[i].mean_p+=pop[j].locus[i].p;
                    }
                }
            }

            write_output(mainout);
            /*if (codominant==0.5)
            {
                if (!SNP_genotypes)
                    write_output_intensity(out_intensity); // selection2
                if (opt->getFlag( "out_freq" ))
                    write_output_genotype(out_genotype); // selection2
            }
            */
            // public version: remove freq output
            /*if (opt->getFlag( "out_freq" ))
            {
                write_anc_freq(ancfreq);
                if (codominant<1)  // selection2
                    write_freq(freq_pop2);
            }*/

        }

//print acceptance rate
        if (((iter/acc_rate)*acc_rate == iter) && (iter > 0))
        {
            //double alpha_up=alpha_updates/I;
            //if (alpha_up==0) alpha_up=1;

            if (codominant==0) // selection2
            {
                outAccRte << setw(3) << setprecision(3) << setiosflags(ios::showpoint)
                /*<< acc_alpha[0]/(alpha_up) << " "*/ << acc_beta[0]/(iter)   << " " << acc_freq_ancestral[0]/(iter)
                << " "  << acc_freq[0][0]/(iter) << " " /*<< acc_a_p/iter << " "*/ <<  endl;
            }
            else  if (codominant==1)
            {
                outAccRte << setw(3) << setprecision(3) << setiosflags(ios::showpoint)
                /*<< acc_alpha[0]/(alpha_up) << " "*/ << acc_beta[0]/(iter)   << " " << acc_freq_ancestral[0]/(iter)
                << endl;
            }
            else
            {
                outAccRte << setw(3) << setprecision(3) << setiosflags(ios::showpoint)
               /* << acc_alpha[0]/(alpha_up) << " "*/ << acc_beta[0]/(iter)   << " " << acc_freq_ancestral[0]/(iter)
                << " "  << acc_freq[0][0]/(iter) << " " /*<< acc_a_p/iter << " " */<<  acc_f[0]/(iter) << " " <<  acc_mu[0]/(iter) << " " << acc_delta[0]/(iter) << " "  << acc_sigma1[0]/(iter) << " "  << acc_sigma2[0]/(iter)<< endl;
            }
        }

    } // end of MCMC loop



    cout << endl;

    mainout.close();

//    if (codominant==0.5 && !SNP_genotypes)
//        out_intensity.close();

    if (opt->getFlag( "out_freq" ) && codominant<1)
    {
        // public version: remove freq output
/*        ancfreq.close();
        for (int j=0;j<J;j++)
        {
            freq_pop2[j].close();
            //if (codominant==0.5)
            //    out_genotype[j].close();
        }*/
        if (cur_out>0 /*&& codominant<1*/)
            write_freq_matrix(freq_pop,cur_out);
        freq_pop.close();

    }

// write mean fst values
    if (!opt->getFlag( "fstat" ))
    {
        fst_file << "    " << "prob" << "  log10(PO)" << "   qval" << "  alpha" << "  fst" << endl;
	}

    /*if (cur_out>0)
    {
        double p_rj;
        double bf;

        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
            {
                p_rj=(double)nb_alpha[i]/(double)cur_out;
                if (p_rj==0)
                    bf=-1000;
                else if (p_rj==1)
                    bf=1000;
                else bf=log10l(p_rj/(1-p_rj));
                if (!opt->getFlag( "fstat" ))
                    fst_file << (i+1) << "  " << p_rj << " " << bf << " " <<  setw(5) << setprecision(5) << setiosflags(ios::showpoint) << post_alpha[i]/(double)cur_out << " " <<  post_fst[i]/(double)cur_out << endl;
            }
        }
    }*/

   if (cur_out>0)
    {
        double *p_rj=new double[I];
        double *po=new double[I];
        double *q_val=new double[I];

        // first calculate PO and posterior prob for all loci 
        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
            {
                // between groups
                p_rj[i]=(double)nb_alpha[i]/(double)cur_out;
                if (p_rj[i]==0)
                    po[i]=-1000;
                else if (p_rj[i]==1)
                    po[i]=1000;
                else po[i]=log10l(p_rj[i]/(1-p_rj[i]));
            }
        }
        // now from posterior proba. calculate q-values

        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
            {
                double threshold=p_rj[i];
                int significants=0;
                q_val[i]=0;
                for (int i2=0;i2<I;i2++)
                {
                    if (!discarded_loci[i2] && p_rj[i2]>=threshold) {
                        q_val[i]+=(1-p_rj[i2]);
                        significants++;
                    }
                }
                q_val[i]/=significants;
            }
        }
        
        // now write everything in the output file
        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
            {
                if (!opt->getFlag( "fstat" ))
                    fst_file << (i+1) << "  " << p_rj[i] << " " << po[i] << " " << q_val[i] << " " <<  setw(5) << setprecision(5) << setiosflags(ios::showpoint) << post_alpha[i]/(double)cur_out << " " <<  post_fst[i]/(double)cur_out << " " << endl;
            }
        }
    }

    if (!opt->getFlag( "fstat" ))
        fst_file.close();

// end of process

}

