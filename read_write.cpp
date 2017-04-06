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

#include "global_defs.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include "errors.cpp"

#include <vector>



using namespace std;

int read_input(std::ifstream& infile,string outcheck)
{

    ofstream outfile;


    outfile.open(outcheck.c_str());
    assure(outfile);

    outfile << "Summary of parameters and input files." << endl;
    outfile << "Please check that all is correct while calculation is starting..." << endl << endl;

    string line;

    while (getline(infile,line,'=') )
    {
        if (line.length()>=6 && line.substr(line.length()-6,line.length())=="[loci]") //read nb of loci
        {
            getline(infile,line);
            istringstream read_line(line);
            read_line >> I;
            outfile << "There are " << I << " loci." << endl << endl;
            try
            {
                alpha=new double[I];
            }
            catch ( const std::exception & Exp )
            {
                cout << "Not enough memory for alpha" << endl;
                return 1;
            }
            try
            {
                alpha_included=new bool[I];
            }
            catch ( const std::exception & Exp )
            {
                cout << "Not enough memory for alpha_included" << endl;
                return 1;
            }
            try
            {
                discarded_loci=new bool[I];
            }
            catch ( const std::exception & Exp )
            {
                cout << "Not enough memory for discarded_loci"<< endl;
                return 1;
            }
            try
            {
                mean_alpha=new double[I];
            }
            catch ( const std::exception & Exp )
            {
                cout << "Not enough memory for mean_alpha"<< endl;
                return 1;
            }
            try
            {
                var_alpha=new double[I];
            }
            catch ( const std::exception & Exp )
            {
                cout << "Not enough memory for var_alpha"<< endl;
                return 1;
            }
            try
            {
                post_alpha=new double[I];
            }
            catch ( const std::exception & Exp )
            {
                cout << "Not enough memory for post_alpha"<< endl;
                return 1;
            }
            /*try {
            	  p_value=new double[I];
            }
            catch ( const std::exception & Exp ) {
            	cout << "Not enough memory for p_value"<< endl;
            	return 1;
            } */
            try
            {
                nb_alpha=new int[I];
            }
            catch ( const std::exception & Exp )
            {
                cout << "Not enough memory for nb_alpha"<< endl;
                return 1;
            }
            try
            {
                cur_fst=new double[I];
            }
            catch ( const std::exception & Exp )
            {
                cout << "Not enough memory for cur_fst"<< endl;
                return 1;
            }
            try
            {
                post_fst=new double[I];
            }
            catch ( const std::exception & Exp )
            {
                cout << "Not enough memory for post_fst"<< endl;
                return 1;
            }
            try
            {
                freq_ancestral=new double[I];
            }
            catch ( const std::exception & Exp )
            {
                cout << "Not enough memory for freq_ancestral"<< endl;
                return 1;
            }
            if (codominant==1)   // selection2 : I forgot this before (waste of unused memory) ?
            {
                try
                {
                    freq_locus=new allele_freq[I];
                }
                catch ( const std::exception & Exp )
                {
                    cout << "Not enough memory for allele_freq "<< endl;
                    return 1;
                }
            }

//selection2
            try
            {
                acc_alpha = new double [I];
                var_prop_alpha = new double [I];
                acc_freq_ancestral = new double [I];
                e_ancestral = new double [I];
            }
            catch ( const std::exception & Exp )
            {
                cout << "Not enough memory for acc_rte alpha or freq_ancestral "<< endl;
                return 1;
            }

        }

        else if (line.length()>=13 && line.substr(line.length()-13,line.length())=="[populations]") //read nb of populations
        {
            getline(infile,line);
            istringstream read_line(line);
            read_line >> J;
            outfile << "There are " << J << " populations." << endl << endl;
            try
            {
                pop = new pop_data[J];
            }
            catch ( const std::exception & Exp )
            {
                cout << "Not enough memory for pop_data "<< endl;
                return 1;
            }

            try
            {
                beta=new double[J];
            }
            catch ( const std::exception & Exp )
            {
                cout << "Not enough memory for beta"<< endl;
                return 1;
            }

            if (codominant<1)
            {
                try
                {
                    f=new double[J];
                }
                catch ( const std::exception & Exp )
                {
                    cout << "Not enough memory for f"<< endl;
                    return 1;
                }
            }

            try    // selection2
            {
                acc_beta=new double[J];
                var_prop_beta=new double[J];
            }
            catch ( const std::exception & Exp )
            {
                cout << "Not enough memory for acc_rte beta"<< endl;
                return 1;
            }

            if (codominant<1)   // selection2
            {
                try
                {
                    acc_f=new double[J];
                    e_f=new double[J];
                }
                catch ( const std::exception & Exp )
                {
                    cout << "Not enough memory for acc_rte f"<< endl;
                    return 1;
                }
            }

        }

        else if (line.length()>=5 && line.substr(line.length()-5,line.length())=="[pop]") //read allele count
        {
// initialize fst table
            /*fst=(double**)malloc(I*sizeof(double *));
            for(int i=0;i<I;i++)
            fst[i]=(double*)malloc(J*sizeof(double));*/

            try
            {
                acc_freq = new double*[I];
                e_freq = new double*[I];
                for (int i=0; i < I; i++)
                {
                    acc_freq[i]= new double[J];
                    e_freq[i]= new double[J];
                }
            }
            catch ( const std::exception & Exp )
            {
                cout << "Not enough memory for acc_freq"<< endl;
                return 1;
            }

            getline(infile,line);
            istringstream read_line(line);
            int j;
            read_line >> j; // population index
            j=j-1;

            try
            {
                pop[j].locus = new locus_data[I];
            }
            catch ( const std::exception & Exp )
            {
                cout << "Not enough memory for pop[" << (j+1) << "].locus "<< endl;
                return 1;
            }

            for (int i=0;i<I;i++) // locus index
            {
                getline(infile,line); // read the line of allele count at locus i
                istringstream read_line(line);

                int locus_index;  // check the locus index is correct
                read_line >> locus_index;
                if (locus_index!=i+1)
                {
                    cout << "Could not find allele count for population " << (j+1) << " at locus "<< (i+1) << endl;
                    return 1;
                }

                if (codominant==1) // selection2
                {
                    read_line >> pop[j].locus[i].alleleCount;  // total number of alleles in the sample at locus i

                    read_line >> pop[j].locus[i].ar;   // number of allelic classes at locus i
                    
                    bool ar1=false;
                    if (pop[j].locus[i].ar==1) {
						ar1=true;
						pop[j].locus[i].ar=2;
                    }
                    
                    freq_locus[i].ar=pop[j].locus[i].ar; // the number is the same for all populations
                    try
                    {
                        pop[j].locus[i].data_allele_count = new int[pop[j].locus[i].ar];
                    }
                    catch ( const std::exception & Exp )
                    {
                        cout << "Not enough memory for pop["<< (j+1) <<"locus["<<(i+1)<<"].data_allele_count"<<endl;
                        return 1;
                    }
                    if (j==0)
                    {
                        try
                        {
                            freq_locus[i].allele=new double[freq_locus[i].ar];
                        }
                        catch ( const std::exception & Exp )
                        {
                            cout << "Not enough memory for freq_locus["<< (i+1) <<"].allele"<< endl;
                            return 1;
                        }
                    }
                    if (!ar1)
                    {
                        for (int k=0;k<pop[j].locus[i].ar;k++)
                                read_line >> pop[j].locus[i].data_allele_count[k];
                    }
                    else
                    {
						read_line >> pop[j].locus[i].data_allele_count[0];
						pop[j].locus[i].data_allele_count[1]=0;
		    		}
                }
                else
                {
                    read_line >> pop[j].locus[i].n;
                    read_line >> pop[j].locus[i].nA1;
                }

            }

        }

    }

// write check file

    outfile << "Burn in: " << discard << endl;

    outfile << "Thining interval: " << interval << endl;

    outfile << "Sample size: " << nr_out << endl;

    outfile << "Resulting total number of iterations: " << tot_nr_of_iter << endl ;

    outfile << "Nb of pilot runs: " << nb_pilot << endl;

    outfile << "Length of each pilot run: " << pilot_length << endl << endl;


// write allele counts
    outfile << "Allele counts:" << endl;
    for (int j=0;j<J;j++)
    {
        for (int i=0;i<I;i++)
        {
            outfile << "Pop. " << j+1 << " locus " << i+1 << " : " ;
            if (codominant==1) // selection2
            {
                for (int k=0;k<pop[j].locus[i].ar;k++)
                    outfile << setw(3) << pop[j].locus[i].data_allele_count[k] << " " ;
            }
            else
                outfile << setw(3) << pop[j].locus[i].n << " " << pop[j].locus[i].nA1;
            outfile << endl;
        }
        outfile << endl;
    }
    outfile.close();



    infile.close();

    return 0;
}

int read_input_intensity(string outname,string outcheck)
{
    ofstream outfile;

    outfile.open(outcheck.c_str());
    assure(outfile);

    outfile << "Summary of parameters and input files." << endl;
    outfile << "Please check that all is correct while calculation is starting..." << endl << endl;


    vector<int> individuals; // size should be always the nb of pop , contains the number of indiv in each pop
    double trash;
    int cur_pop;
    int cur_indiv;

// first pass to read number of pop, loci, and indivs
    ifstream infile(outname.c_str());
    assure(infile);

    I=0;
    J=0;

    string line; // contains each line

// use the first line for number of loci
    getline(infile,line);
    istringstream read_line(line);
    read_line >> cur_indiv;
    read_line >> cur_pop;
    if (cur_pop>J)
    {
        J=cur_pop;
        individuals.resize(J,0);
    }
    if (cur_indiv>individuals[cur_pop-1])
        individuals[cur_pop-1]=cur_indiv;

//while(read_line.good())
    while (read_line >> trash)
    {
//	read_line >> trash;
        I++;
    }

//I--;

    while (getline(infile,line))
    {
        istringstream read_line(line);
        read_line >> cur_indiv;
        read_line >> cur_pop;
        if (cur_pop>J)
        {
            J=cur_pop;
            individuals.resize(J,0);
        }
        if (cur_indiv>individuals[cur_pop-1])
            individuals[cur_pop-1]=cur_indiv;
    }


    infile.close();
    infile.clear();
// second pass to fill the database
    infile.open(outname.c_str());

// allocate memory
    try
    {
        locus_likelihood=new double[I];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for locus_likelihood"<< endl;
        return 1;
    }
    try
    {
        alpha=new double[I];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for alpha"<< endl;
        return 1;
    }
    try
    {
        discarded_loci=new bool[I];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for discarded_loci"<< endl;
        return 1;
    }
    try
    {
        alpha_included=new bool[I];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for alpha_included"<< endl;
        return 1;
    }
    try
    {
        mean_alpha=new double[I];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for mean_alpha"<< endl;
        return 1;
    }
    try
    {
        var_alpha=new double[I];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for var_alpha"<< endl;
        return 1;
    }
    try
    {
        post_alpha=new double[I];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for post_alpha"<< endl;
        return 1;
    }
    try
    {
        nb_alpha=new int[I];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for nb_alpha"<< endl;
        return 1;
    }
    try
    {
        cur_fst=new double[I];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for cur_fst"<< endl;
        return 1;
    }
    try
    {
        post_fst=new double[I];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for post_fst"<< endl;
        return 1;
    }
    try
    {
        freq_ancestral=new double[I];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for freq_ancestral"<< endl;
        return 1;
    }
    try
    {
        mu=new double[I];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for mu"<< endl;
        return 1;
    }
    try
    {
        delta=new double[I];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for delta"<< endl;
        return 1;
    }
    try
    {
        sigma1=new double[I];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for sigma1"<< endl;
        return 1;
    }
    try
    {
        sigma2=new double[I];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for sigma2"<< endl;
        return 1;
    }
    try
    {
        max_intenstiy=new double[I];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for max_intenstiy"<< endl;
        return 1;
    }

    try
    {
        acc_alpha = new double [I];
        var_prop_alpha = new double [I];
        acc_freq_ancestral = new double [I];
        e_ancestral = new double [I];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for acc_rte alpha or freq_ancestral "<< endl;
        return 1;
    }

    try
    {
        acc_mu = new double [I];
        acc_delta = new double [I];
        acc_sigma1 = new double [I];
        acc_sigma2 = new double [I];
        e_mu = new double [I];
        e_delta = new double [I];
        var_prop_sigma1 = new double [I];
        var_prop_sigma2 = new double [I];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for acc_rte mu, delta or sigma "<< endl;
        return 1;
    }

    try
    {
        pop = new pop_data[J];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for pop_data "<< endl;
        return 1;
    }
    try
    {
        f=new double[J];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for f"<< endl;
        return 1;
    }
    try
    {
        beta=new double[J];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for beta"<< endl;
        return 1;
    }
    try
    {
        acc_beta=new double[J];
        var_prop_beta=new double[J];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for acc_rte beta"<< endl;
        return 1;
    }

    try
    {
        acc_f=new double[J];
        e_f=new double[J];
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for acc_rte f"<< endl;
        return 1;
    }

    for (int j = 0; j < J; j++)
    {
        try
        {
            pop[j].locus = new locus_data[I];
        }
        catch ( const std::exception & Exp )
        {
            cout << "Not enough memory for pop[" << (j+1) << "].locus ";
            return 1;
        }
    }

    try
    {
        acc_freq = new double*[I];
        e_freq = new double*[I];
        for (int i=0; i < I; i++)
        {
            acc_freq[i]= new double[J];
            e_freq[i]= new double[J];
        }
    }
    catch ( const std::exception & Exp )
    {
        cout << "Not enough memory for acc_freq"<< endl;
        return 1;
    }

    for (int j = 0; j < J; j++)
    {
        for (int i = 0; i < I; i++)
        {
            try
            {
                pop[j].locus[i].indiv = new indiv_data[individuals[j]];
            }
            catch ( const std::exception & Exp )
            {
                cout << "Not enough memory for pop[" << (j+1) << "].locus ";
                return 1;
            }
            pop[j].locus[i].n=individuals[j];
        }
    }

    for (int i = 0; i < I; i++)
        max_intenstiy[i]=0;

    while (getline(infile,line))
    {
        istringstream read_line(line);
        read_line >> cur_indiv;
        read_line >> cur_pop;

        // read intensities
        for (int i = 0; i < I; i++)
        {
            read_line >> pop[cur_pop-1].locus[i].indiv[cur_indiv-1].intensity;
            if (pop[cur_pop-1].locus[i].indiv[cur_indiv-1].intensity>max_intenstiy[i])
                max_intenstiy[i]=pop[cur_pop-1].locus[i].indiv[cur_indiv-1].intensity;
        }
    }

    if (!SNP_genotypes)
    {
        for (int i=0;i<I;i++)
        {
            for (int j=0;j<J;j++)
            {
                for (int k = 0; k < pop[j].locus[i].n; k++)
                {
                    if (max_intenstiy[i]>0)
                        pop[j].locus[i].indiv[k].intensity/=max_intenstiy[i];
                }
            }
        }
    }



    infile.close();
    infile.clear();

// write check file
    outfile << "There are " << I << " loci." << endl << endl;

    outfile << "There are " << J << " populations." << endl << endl;

    for (int j=0; j < J; j++)
    {
        outfile << "There are " << individuals[j] << " individuals in population " << j+1 << "." << endl;
    }
    outfile << endl;

    outfile << "Burn in: " << discard << endl;

    outfile << "Thining interval: " << interval << endl;

    outfile << "Sample size: " << nr_out << endl;

    outfile << "Resulting total number of iterations: " << tot_nr_of_iter << endl ;

    outfile << "Nb of pilot runs: " << nb_pilot << endl;

    outfile << "Length of each pilot run: " << pilot_length << endl << endl;

    outfile.close();


    return 0;
}

//////////////////////////////
// read discarded loci file
//////////////////////////////
int read_discarded(std::ifstream& infile)
{
    int locus;
    while (infile >> locus)
    {
        if (locus>=1 && locus <=I)
            discarded_loci[locus-1]=true;
    }
}

///////////////////////////////////
// write output file
///////////////////////////////////
void write_output(std::ofstream& outfile)
{
// write iteration index
    outfile << iter << "  ";
// write loglikelihood
    outfile << setprecision(8) << " " << log_likelihood;
    if (codominant<1)
    {
// write a
//outfile << setprecision(8) << " " << a_p;
// write fis
        for (int j=0;j<J;j++)
            outfile << setprecision(8) << " " << f[j];
    }
// write fst
    for (int j=0;j<J;j++)
        outfile << setprecision(8) << " " << 1/(1+exp(-beta[j]));
// write prop of non zero alpha
//outfile << setprecision(8) << " " << ((double)nb_alpha_included)/((double)I);
// write alpha
    if (!fstat && all_trace)
    {
        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
                outfile << setprecision(8) << " " << alpha[i];
        }
    }
// write alpha_included
//for (int i=0;i<I;i++)
//outfile << setprecision(8) << " " << (int)alpha_included[i];
// write fst
//for (int i=0;i<I;i++)
//outfile << setprecision(8) << " " << cur_fst[i];

    outfile << endl;
}

///////////////////////////////////
// write output file for mixture model of AFLP band intensity
///////////////////////////////////
void write_output_intensity(std::ofstream& outfile) // selection2
{
// write iteration index
    outfile << iter << "  ";
// write mu and mu+delta
    for (int i=0;i<I;i++)
    {
        if (!discarded_loci[i])
        {
            if (max_intenstiy[i]!=0)
                outfile << setprecision(8) << " " << mu[i]*max_intenstiy[i];
            else
                outfile << setprecision(8) << " " << mu[i];

            if (max_intenstiy[i]!=0)
                outfile << setprecision(8) << " " << (delta[i]+mu[i])*max_intenstiy[i];
            else
                outfile << setprecision(8) << " " << delta[i]+mu[i];
        }
    }
// write sigma1 and sigma2
    for (int i=0;i<I;i++)
    {
        if (!discarded_loci[i])
        {
            if (max_intenstiy[i]!=0)
                outfile << setprecision(8) << " " << sigma1[i]*max_intenstiy[i];
            else
                outfile << setprecision(8) << " " << sigma1[i];
            if (max_intenstiy[i]!=0)
                outfile << setprecision(8) << " " << sigma2[i]*max_intenstiy[i];
            else
                outfile << setprecision(8) << " " << sigma2[i];
        }
    }
    outfile << endl;
}


///////////////////////////////////
// write genotype freq Aa / AA for mixture model of AFLP band intensity
///////////////////////////////////
void write_output_genotype(std::ofstream outfile[]) // selection2
{
    double p,pAa,pAA;
    for (int j=0;j<J;j++)
    {
        // write iteration index
        outfile[j] << iter << "  ";
        // write pAa and pAA
        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
            {
                p=pop[j].locus[i].p;
                pAA=p*p+f[j]*p*(1-p);
                pAa=2*p*(1-p)*(1-f[j]);
                outfile[j] << setprecision(8) << " " << pAa << " " << pAA;
            }
        }
        outfile[j] << endl;
    }
}

/////////////////////////////////////////////
// write output file for allele frequencies
/////////////////////////////////////////////
void write_freq(std::ofstream outfiles[])
{
    for (int j=0;j<J;j++)
    {
        outfiles[j] << iter;
        for (int i=0;i<I;i++) //cycle over loci
        {
            if (!discarded_loci[i])
                outfiles[j] << " " << pop[j].locus[i].p;
        }
        outfiles[j] << endl;
    }
}

/////////////////////////////////////////////
// write single output file for posterior allele frequencies
/////////////////////////////////////////////
// public version: added "popxx" in front of each line
void write_freq_matrix(std::ofstream& freq_pop,int cur_out)
{
    // write header
    freq_pop << "    ";
    for (int i=0;i<I;i++)
    {
        if (!discarded_loci[i])
            freq_pop << "locus" << (i+1) << " ";
    }
    freq_pop << endl;
    for (int j=0;j<J;j++)
    {
        freq_pop << "pop" << (j+1) << " ";
        for (int i=0;i<I;i++)
        {
            if (!discarded_loci[i])
                freq_pop << pop[j].locus[i].mean_p/((double)cur_out) << " ";
        }
        freq_pop << endl;
    }
}

/////////////////////////////////////////////
// write output file for ancestral allele frequencies
/////////////////////////////////////////////
void write_anc_freq(std::ofstream& outfile)
{
    outfile << iter  ;

    for (int i=0;i<I;i++) //cycle over loci
    {
        if (!discarded_loci[i])
        {
            if (codominant==1) // selection2
            {
                outfile << " " << freq_locus[i].ar;
                for (int k=0;k<freq_locus[i].ar;k++)
                    outfile << " " << freq_locus[i].allele[k];
            }
            else
                outfile << " " << freq_ancestral[i];
        }
    }
    outfile << endl;
}
