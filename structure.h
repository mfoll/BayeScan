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

/*
This header file contains the definitions of data structures.
1) A structure that stores the characteristics of each locus, i.e.
number of alleles and allele frequency spectrum.
2) A structure that stores the characteristics of each population, i.e.
sample size and allele frequency spectrum for each locus.
*/

struct indiv_data  // selection2
{
    float intensity;
};

struct locus_data
{
    int n; // total number of alleles in the sample = n_ij
    int nA1;  /* observed allele frequency spectrum = nA1_ij*/
    float p; // estimated allele frequency = p_ij
    double mean_p; // mean value of p so far (after pilots and burn-in)
    struct indiv_data *indiv; // selection2

    int ar; /*number of allelic classes = K_i*/
    int alleleCount; // total number of alleles in the sample = n_ij
    int *data_allele_count;//[max_nr_alleles];  /* observed allele frequency spectrum = a_ij*/
};

struct pop_data
{
    locus_data *locus;//[max_nr_loci]; // alleles count
};


struct allele_freq
{
    int ar; /*number of allelic classes = K_i*/
    double *allele;//[max_nr_alleles]; /* allele frequency */
};
