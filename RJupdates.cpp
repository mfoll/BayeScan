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

#include <algorithm>
#include <math.h>

#include <omp.h>

//#define prop_alpha_mean 0
//#define prop_alpha_var 5

using namespace std;

void jump_model()
{
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double r; // random value to accept/reject the move
    double old_alpha;
    double new_theta,old_theta; // new and old values of theta

    #pragma omp parallel for SCHED_I /*reduction(+:nb_alpha_included)*/ private(r, A, old_alpha, new_theta,old_theta)
    for (int i=0;i<I;i++) // cycle over loci
    {
        if (!discarded_loci[i])
        {
    old_alpha=alpha[i];

// propose new alpha value
    if (!alpha_included[i])
                alpha[i]=randgen_parallel[omp_get_thread_num()].randNorm(mean_alpha[i],sqrt(var_alpha[i]));
    else
        alpha[i]=0;

// change the state of alpha
    alpha_included[i]=!alpha_included[i];

// calculate A
    A=0;

    for (int j=0;j<J;j++)
    {
        // calculate old and new value of theta
        new_theta=exp(-(alpha[i]+beta[j]));
        old_theta=exp(-(old_alpha+beta[j]));

        A+= gammaln(new_theta)-gammaln(old_theta)
            - gammaln(new_theta*freq_ancestral[i]) + gammaln(old_theta*freq_ancestral[i])
            - gammaln(new_theta*(1-freq_ancestral[i])) + gammaln(old_theta*(1-freq_ancestral[i]))
            + freq_ancestral[i]*(new_theta-old_theta)*log(pop[j].locus[i].p)
            + (1-freq_ancestral[i])*(new_theta-old_theta)*log(1-pop[j].locus[i].p);
    }

    if (alpha_included[i]) //if we add parameter
        A+= log_prior_alpha(alpha[i])//-0.5*log(2*M_PI*sd_prior_alpha*sd_prior_alpha)-(alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha)
            -(-0.5*log(2*M_PI*var_alpha[i])-((alpha[i]-mean_alpha[i])*(alpha[i]-mean_alpha[i]))/(2*var_alpha[i]))-log(prior_odds);
// inverse if we remove
    else
        A+= (-0.5*log(2*M_PI*var_alpha[i])-((alpha[i]-mean_alpha[i])*(alpha[i]-mean_alpha[i]))/(2*var_alpha[i]))
            -log_prior_alpha(alpha[i])+log(prior_odds);//-(-0.5*log(2*M_PI*sd_prior_alpha*sd_prior_alpha)-(alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha));

            r=randgen_parallel[omp_get_thread_num()].randDblExc();

// reject proposed value
    if (log(r)>A)
    {
        alpha[i]=old_alpha;
        alpha_included[i]=!alpha_included[i];
    }
   /* else
    {
        if (alpha_included[i])
            nb_alpha_included++;
        else
            nb_alpha_included--;
    }*/
}
    }
}


void jump_model_codominant()
{
    double A; // log of A in MH algorithm : accept move with probability min(1,A)
    double r; // random value to accept/reject the move
    double old_alpha;
    double new_theta,old_theta; // new and old values of theta

    //double old_log_likelihood; // old loglikelihood
    double diff_log_likelihood; // old loglikelihood
    #pragma omp parallel for SCHED_I reduction(+:/*nb_alpha_included,*/log_likelihood) private(r, A, old_alpha, new_theta,old_theta,diff_log_likelihood)
    for (int i=0;i<I;i++) // cycle over loci
    {
        if (!discarded_loci[i])
        {
    old_alpha=alpha[i];

// propose new alpha value
    if (!alpha_included[i])
                alpha[i]=randgen_parallel[omp_get_thread_num()].randNorm(mean_alpha[i],sqrt(var_alpha[i]));
    else
        alpha[i]=0;

// change the state of alpha
    alpha_included[i]=!alpha_included[i];

// calculate A
    A=0;

    double old_l=0;
    for (int j=0;j<J;j++)
    {
        // calculate old and new value of theta
        old_theta=exp(-(old_alpha+beta[j]));

        old_l+=gammaln(old_theta)-gammaln(pop[j].locus[i].alleleCount+old_theta);
        for (int k=0;k<pop[j].locus[i].ar;k++)
            old_l+=gammaln(pop[j].locus[i].data_allele_count[k]+old_theta*freq_locus[i].allele[k])
                   -gammaln(old_theta*freq_locus[i].allele[k]);
    }

    double new_l=0;
    for (int j=0;j<J;j++)
    {
        // calculate old and new value of theta
        new_theta=exp(-(alpha[i]+beta[j]));

        new_l+=gammaln(new_theta)-gammaln(pop[j].locus[i].alleleCount+new_theta);
        for (int k=0;k<pop[j].locus[i].ar;k++)
            new_l+=gammaln(pop[j].locus[i].data_allele_count[k]+new_theta*freq_locus[i].allele[k])
                   -gammaln(new_theta*freq_locus[i].allele[k]);
    }

// store the old loglikelihood and calculate the new loglikelihood
    //old_log_likelihood=log_likelihood;

    //log_likelihood=old_log_likelihood-old_l+new_l;
    diff_log_likelihood=-old_l+new_l;
    A=diff_log_likelihood;

    if (alpha_included[i]) //if we add parameter
        A+= log_prior_alpha(alpha[i])//-0.5*log(2*M_PI*sd_prior_alpha*sd_prior_alpha)-(alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha)
            -(-0.5*log(2*M_PI*var_alpha[i])-((alpha[i]-mean_alpha[i])*(alpha[i]-mean_alpha[i]))/(2*var_alpha[i]))-log(prior_odds);
// inverse if we remove
    else
        A+= (-0.5*log(2*M_PI*var_alpha[i])-((alpha[i]-mean_alpha[i])*(alpha[i]-mean_alpha[i]))/(2*var_alpha[i]))
            -log_prior_alpha(alpha[i])+log(prior_odds);//-(-0.5*log(2*M_PI*sd_prior_alpha*sd_prior_alpha)-(alpha[i]*alpha[i])/(2*sd_prior_alpha*sd_prior_alpha));

            r=randgen_parallel[omp_get_thread_num()].randDblExc();

// reject proposed value
    if (log(r)>A)
    {
        alpha[i]=old_alpha;
        alpha_included[i]=!alpha_included[i];
        //log_likelihood=old_log_likelihood;
    }
    else
    {
        log_likelihood=log_likelihood+diff_log_likelihood;
       /* if (alpha_included[i])
            nb_alpha_included++;
        else
            nb_alpha_included--;*/
    }

    }
}

}
