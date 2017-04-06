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
#include <math.h>

double dirichlet(double vectX[],double vectA[],int arraySize)
/*This function returns the log of the dirichlet density with parameters vectA
for a vectorX */
{
    int i;
    double Sumlogxi;
    double Sumloggamma,SumDirAlphas;
    double density; //this is the density that will be returned by the function
//   float gammaln(float xx);

    Sumloggamma=0;
    Sumlogxi=0;
    SumDirAlphas=0;
    for (i=0;i<arraySize;i++)
    {
        Sumloggamma += gammaln(vectA[i]);
        Sumlogxi += (vectA[i]-1.)*logl(vectX[i]);
        SumDirAlphas += vectA[i];
    }

    density = gammaln(SumDirAlphas) - Sumloggamma + Sumlogxi;

    return (density);
}



float DumRan;

//----------------------------------
/* This function generates a vector (of dimension arraySize) of Dirichlet
random variables with parameters vectA and stores them in vect. */

void dirichlet_dev(double vectX[],double vectA[],int arraySize)
{
    int i;
    double sum=0.0;
    double rstgam(double alpha1);


// now we get the deviates

    for (i=0;i<arraySize;i++)
    {
        if (vectA[i]>1.0)
            vectX[i]=rstgam(vectA[i]);
        else if (vectA[i]<1.0)
        {
            /*do {
                DumRan = ranf();
                }
                while (DumRan >= 1);*/
            DumRan=randgen.randDblExc();

            vectX[i]=rstgam(vectA[i]+1.0)*pow((double) DumRan,1.0/vectA[i]);
        }
        else
        {
            /* do {
                 DumRan = ranf();
                 }
                 while (DumRan >= 1);*/
            DumRan=randgen.randDblExc();

            vectX[i]=(-logl((double) DumRan));
        }

        sum += vectX[i];
    }

    for (i=0;i<arraySize;i++)
        vectX[i]=vectX[i]/sum;

    return;
}



/* Chen and Feast algorithm (algorithm 3.20 in Ripley 1988) for generating a
standard Gamma random variable with shape parameter alpha>1) */

double rstgam(double alpha1)
{
    int acpt;
    double c1,c2,c3,c4,c5,u1,u2,w;

    c1=alpha1-1.0;
    c2=(alpha1-(1.0/(6.0*alpha1)))/c1;
    c3=2.0/c1;
    c4=c3+2.0;
    c5=1.0/sqrt(alpha1);

    do
    {
        do
        {
            /*do {
                 DumRan = ranf();
                 }
                 while (DumRan >= 1);
            u1=DumRan;*/
            u1=randgen.randDblExc();

            /*do {
                 DumRan = ranf();
                 }
                 while (DumRan >= 1);
            u2=DumRan;*/
            u2=randgen.randDblExc();

            if (alpha1>2.5)
                u1=u2+c5*(1.0-(1.86*u1));
        }
        while (u1>=1.0 || u1<=0);

        w=(c2*u2)/u1;
        if ((c3*u1+w+(1.0/w))<=c4)
            acpt=1;
        else if ((c3*logl(u1)-logl(w)+w)>=1.0)  //repeat outer do loop
            acpt=0;
        else
            acpt=1;
    }
    while (acpt!=1);

    return c1*w;
}



