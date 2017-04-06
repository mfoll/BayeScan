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
#include <stdlib.h>
#include <math.h>
#include <iostream>

#define inv_sqrt_2_PI  0.3989422804014327
#define inv_PI   0.318309886183791

#define halfpi 1.57080
#define sixthpi  0.523598
#define tansixthpi 0.577350
#define tantwelfthpi 0.267949

#define epsilon 1e-6 // limit for allele frequencies (epsilon,1-epsilon)
//
// This is the main arctangent approximation "driver"
// It reduces the input argument's range to [0, pi/12],
// and then calls the approximator.
//
//

inline float atan_66(float xx)
{
    float x=xx;
    float y; // return from atan__s function
    int complement= false; // true if arg was >1
    int region= false; // true depending on region arg is in
    int sign= false; // true if arg was < 0
    if (x <0 )
    {
        x=-x;
        sign=true; // arctan(-x)=-arctan(x)
    }
    if (x > 1.0)
    {
        x=1.0/x; // keep arg between 0 and 1
        complement=true;
    }
    if (x > tantwelfthpi)
    {
        x = (x-tansixthpi)/(1+tansixthpi*x); // reduce arg to under tan(pi/12)
        region=true;
    }
    y=x-x*x*x/3; // run the approximation
    if (region) y+=sixthpi; // correct for region we're in
    if (complement)y=halfpi-y; // correct for 1/x if we did that
    if (sign)y=-y; // correct for negative arg
    return (y);
}

///////////////////////////////////////
//  Returns the value ln[gamma(xx)] for xx > 0.
///////////////////////////////////////
double gammaln(double xx)

{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
                          24.01409824083091,-1.231739572450155,
                          0.1208650973866179e-2,-0.5395239384953e-5
                         };
    int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

///////////////////////////////////
//Returns ln(n!)
/////////////////////////////////
double factln(int n)

{
    double gammaln(double xx);
    static double a[101];
    if (n < 0) fprintf(stderr,"Negative factorial in routine factln");
    if (n <= 1) return 0.0;
    if (n <= 100) return a[n] ? a[n] : (a[n]=gammaln(n+1.0)); // In range of table.
    else return gammaln(n+1.0); // Out of range of table.
}



///////////////////////////////////////
// Returns Log(pi(alpha))
//////////////////////////////////////
double log_prior_alpha(double alpha)
{
//return logl( (1/(sd_prior_alpha*sqrt(2*M_PI)))*exp(-(alpha)*(alpha)/(2*sd_prior_alpha*sd_prior_alpha)) );
//return -0.5*logl(2*M_PI*sd_prior_alpha*sd_prior_alpha)-(alpha*alpha)/(2*sd_prior_alpha*sd_prior_alpha);
    return log( 0.5*(1/(sd_prior_alpha*sqrt(2*M_PI)))*exp(-(alpha-m1_prior_alpha)*(alpha-m1_prior_alpha)/(2*sd_prior_alpha*sd_prior_alpha))
                + 0.5*(1/(sd_prior_alpha*sqrt(2*M_PI)))*exp(-(alpha-m2_prior_alpha)*(alpha-m2_prior_alpha)/(2*sd_prior_alpha*sd_prior_alpha)) );
}

// return density of normal distribution evaluated in y with mean m and sd s
inline float pi_norm(float y,float m,float s) // selection2
{
// normal distribution
//	return	max(inv_sqrt_2_PI * (1/s) * exp(-(y-m)*(y-m)/(2*s*s)),0.00001);
// Cauchy distribution
    float tmp=((y-m)/s);
    return inv_PI*1/(s*(1+tmp*tmp));
}

inline float pi_norm_cumulative(float m,float s) // selection2
{
// Cauchy distribution
    return min(0.5+inv_PI*atan_66((abscence_pc-m)/s),1-epsilon);
}

// Return ln(L(y_{i,j})) i<-locus, j<-locus
double intensity_loglikelihood(float y,int i,int j) // selection2
{
    double loglikelihood=0;
    float paa,pAa,pAA;
    float p=pop[j].locus[i].p;

    pAA=p*p+f[j]*p*(1-p);
    pAa=2*p*(1-p)*(1-f[j]);
    paa=(1-p)*(1-p)+f[j]*p*(1-p);

    if (y>=0)
    {
        if (SNP_genotypes)
        {
            if (y==0) loglikelihood=log(paa);
            else if (y==1) loglikelihood=log(pAa);
            else if (y==2) loglikelihood=log(pAA);
        }
        else
        {
            if (y<=abscence_pc)
                loglikelihood=log(paa);
            else
                loglikelihood=log(pAa*pi_norm(y,mu[i],sigma1[i])/(1-pi_norm_cumulative(mu[i],sigma1[i])) + pAA*pi_norm(y,mu[i]+delta[i],sigma2[i])/(1-pi_norm_cumulative(mu[i]+delta[i],sigma2[i])));
        }
    }

    return loglikelihood;
}

////////////////////////////////////////////
// Returns ln(L)
////////////////////////////////////////////
double allelecount_loglikelihood()
{
    double loglikelihood=0;
    double g;
    double theta;

    #pragma omp parallel for SCHED_I reduction(+:loglikelihood) private(g, theta)
    for (int i=0;i<I;i++)
    {
        if (!discarded_loci[i])
        {
            if (codominant==0.5 && !SNP_genotypes) locus_likelihood[i]=0;
            for (int j=0;j<J;j++)
            {
                if (codominant==0)   // selection2
                {
                    g= pop[j].locus[i].p*pop[j].locus[i].p
                       +2*pop[j].locus[i].p*(1-pop[j].locus[i].p)*(1-f[j])
                       +f[j]*pop[j].locus[i].p*(1-pop[j].locus[i].p);

                    loglikelihood += factln(pop[j].locus[i].n)
                                     -factln(pop[j].locus[i].nA1)
                                     -factln(pop[j].locus[i].n-pop[j].locus[i].nA1)
                                     +pop[j].locus[i].nA1*log(g)
                                     +(pop[j].locus[i].n-pop[j].locus[i].nA1)*log(1-g);
                }
                else if (codominant==1)
                {
                    theta=exp(-(alpha[i]+beta[j]));
                    loglikelihood+=factln(pop[j].locus[i].alleleCount)+gammaln(theta)
                                   -gammaln(pop[j].locus[i].alleleCount+theta);
                    for (int k=0;k<pop[j].locus[i].ar;k++)
                        loglikelihood+=gammaln(pop[j].locus[i].data_allele_count[k]+theta*freq_locus[i].allele[k])
                                       -factln(pop[j].locus[i].data_allele_count[k])-gammaln(theta*freq_locus[i].allele[k]);
                }
                else // intensity
                {
                    for (int k = 0; k < pop[j].locus[i].n; k++)
                    {
                        double tmp=intensity_loglikelihood(pop[j].locus[i].indiv[k].intensity,i,j);
                        if (codominant==0.5 && !SNP_genotypes) locus_likelihood[i]+=tmp;
                        loglikelihood+=tmp;
                    }
                }
            }
        }
    }

    return loglikelihood;
}
