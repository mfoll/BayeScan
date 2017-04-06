# 	 This file is used to work with AFLP input files for the software Bayescan in R.

#    This program, BayeScan, aims at detecting genetics markers under selection,
#	 based on allele frequency differences between population. 
#    Copyright (C) 2010  Matthieu Foll
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


# Typical usage: 
# - load this file into R (file/source R code)
# - in R, go to the directory where your input file is (file/change current dir)
# - at the R prompt, type somthing like:
# > pop_subset("input_file.txt","input_file_subset_1-2-5.txt",c(1,2,5))
#
# 	this will read "input_file.txt", only retain populations 1,2,5 and write the result
#	in the new BayeScan input file called "input_file_subset_1-2-5.txt"
#
# > clean_intensity_matrix("input_file_subset_1-2-5.txt","input_file_subset_1-2-5_clean.txt")
#
#	this will remove extreme values from your "input_file_subset_1-2-5.txt" input file
# 	and create a new BayeScan input file called "input_file_subset_1-2-5_clean.txt"
#
# > discard_monomorphic("input_file_subset_1-2-5_clean.txt","discarded.txt")
#
#	this will return you a list of monomorphic markers, ready to use in BayeScan
# 	in a file called "discarded.txt" from your input file "input_file_subset_1-2-5_clean.txt"
#
#

# The function removes extreme values of band intensity (replaced by "-1")
# For a given marker i, we calculate the q-quantile (95% by default) of the 
# distribution of non-zero intensity. Then a value is considered as extreme 
# if it is higher than "distance" (3 by default) times the quantile.
# Extreme values are also printed on the command line of R.
#
# input: the name of a file in the current directory in BayeScan format for
# AFLPs band intensity data (matrix format)
#
# output: the output name for the file modified by the function, will also be
# in BayeScan format for AFLPs band intensity data (matrix format)
#
clean_intensity_matrix<-function(input,output,q=0.95,distance=3)
{
# read the input file
data=read.table(input,colClasses="numeric")
nloci=ncol(data)-2

outloc=TRUE # flag internally used to output on the R prompt
# check each locus at a time
for (i in 1:nloci)
	{
	# find the maximum value as the quantile * distance for this locus
	max=quantile(data[data[,i+2]>0,i+2],q)*distance	
	# check for each individual if the intensity is above the maximum
	for (k in 1:nrow(data))
		{			
			if (!is.na(max) && data[k,i+2]>max) 
				{
				# next 4 lines are used for the output on the R prompt
				if (outloc)
					cat("Max locus",i,"=",max,"\n")
				outloc=FALSE
				cat(" pop",data[k,2],"indiv",data[k,1],"=",data[k,i+2],"removed","\n")
				# put a -1 if the value is extreme
				data[k,i+2]=-1					
				}				
		}
	outloc=TRUE
	}
write.table(data,output,row.names = F, col.names=F)	
sum(data==-1)
}

# The function converts band intensity matrix to binary data in BayeScan format
# negative values (such as "-1") are considered as missing values
# The threshold between presence (1) and absence (0) is calculated independently
# for each markers as a fraction "pc_absence" (10% by default) of the maximum 
# intensity for this marker. It is a good idea to first use the "clean_intensity"
# function before applying this one to discard extreme values.
#
# input: the name of a file in the current directory in BayeScan format for
# AFLPs with band intensity
#
# output: the output file name 
#
convert_to_binary_bayescan<-function(input,outfile,pc_absence=0.1)
{
# read the input file
data=read.table(input,colClasses="numeric")
nb_loci=ncol(data)-2
nb_pops=max(data[,2])

max_intensity=rep(0,nb_loci)
for (i in 1:nb_loci) { 
  max_intensity[i]=max(data[,i+2]) 
}

cat("[loci]=",nb_loci,"\n\n",file=outfile)

cat("[populations]=",nb_pops,"\n\n",file=outfile,append=T)

for (pop in 1:nb_pops)
{
  cur_pop=data[data[,2]==pop,]
  nb_bands=rep(0,nb_loci)
  nb_individuals=rep(0,nb_loci)
  cat("[pop]=",pop,"\n",file=outfile,append=T)
  for (i in 1:nb_loci)
  {
    nb_bands[i]=sum(cur_pop[,i+2]>=pc_absence*max_intensity[i])
    nb_individuals[i]=sum(cur_pop[,i+2]>=0)
    cat(i,nb_individuals[i],nb_bands[i],"\n",file=outfile,append=T)
  }
	cat("\n",file=outfile,append=T)
}
}


# The function converts intensity matrix to binary data in matrix format
# by default "-1" are considered as missing values
# The threshold between presence (1) and absence (0) is calculated independently
# for each markers as a fraction "pc_absence" (10% by default) of the maximum 
# intensity for this marker. It is a good idea to first use the "clean_intensity"
# function before applying this one to discard extreme values.
#
# input: the name of a file in the current directory in BayeScan format for
# AFLPs with band intensity
#
# output: the output name for the file created by the function 
#
#
convert_to_binary_matrix<-function(input,outfile,pc_absence=0.1,miss_char="-1")
{
# read the input file
data=read.table(input,colClasses="numeric")
nb_loci=ncol(data)-2
nb_pops=max(data[,2])

max_intensity=rep(0,nb_loci)
for (i in 1:nb_loci) { 
  max_intensity[i]=max(data[,i+2]) 
}

cat("",file=outfile)

for (pop in 1:nb_pops)
{
  cur_pop=data[data[,2]==pop,]
  nb_individuals=nrow(cur_pop)
  for (i in 1:nb_loci)
  {
    missings=cur_pop[,i+2]==-1
    presence=cur_pop[,i+2]>=pc_absence*max_intensity[i]
    abscence=cur_pop[,i+2]<pc_absence*max_intensity[i]
    cur_pop[presence,i+2]=1
    cur_pop[abscence,i+2]=0
    cur_pop[missings,i+2]=miss_char
  }
  for (k in 1:nb_individuals) {
    cat(as.vector(t(cur_pop[k,])),"\n",file=outfile,append=T)
  }
}
}

# The function converts binary matrix to binary data in BayeScan format
# any value other than 0/1 (such as "-1") is considered as missing value
#
# input: the name of a file in the current directory in binary matrix 
# AFLPDAT format (Ehrich 2006 Molecular Ecology Notes 6: 603-604).
#
# output: the output name for the BayeScan file created by the function 
#
#
convert_from_binary_matrix<-function(input,output)
{
data=read.table(input,header=T)

outfile=output

nb_loci=ncol(data)-2

if (is.integer(data[,2]))
  pops_names=1:max(data[,2])
else
  pops_names=levels(data[,2])
  
nb_pops=length(pops_names)

cat("[loci]=",nb_loci,"\n\n",file=outfile)
cat("[populations]=",nb_pops,"\n\n",file=outfile,append=T)

for (pop in 1:nb_pops)
{
  cur_pop=data[data[,2]==pops_names[pop],]
  nb_bands=rep(0,nb_loci)
  nb_individuals=rep(0,nb_loci)
  cat("[pop]=",pop,"\n",file=outfile,append=T)
  for (i in 1:nb_loci)
  {
    nb_bands[i]=sum(cur_pop[,i+2]==1)
    nb_individuals[i]=sum(cur_pop[,i+2]==0)+sum(cur_pop[,i+2]==1)
    cat(i,nb_individuals[i],nb_bands[i],"\n",file=outfile,append=T)
  }
	cat("\n",file=outfile,append=T)
}
}                          

# This function allows to easily extract a subset of populations ("subpop") 
# from a BayeScan band intensity matrix input file.
# BE CAREFUL, in that case population indices are modified, the order is the
# same but they are shifted to be from 1 to length(subpop) as required by 
# BayeScan
#
# input: the name of a file in the current directory in BayeScan format for
# AFLPs with band intensity
#
# output: the output name for the file modified by the function containing 
# only the subset of populations given in the vector subpop 

pop_subset<-function(input,output,subpop)
{	
# read the input file
data=read.table(input,colClasses="numeric")

# sort the populations we will retain
subpop=sort(subpop)

# output the file but only for populations in subpop 
# first create the data frame structure
data_sub=data.frame(data[1,])
# index of current individual outputed.
cur_ind=1;
# loop over all individuals
for (k in 1:nrow(data))
{
  # check if the pop of the current individual is in subpop
  if (is.element(data[k,2],subpop))
  {
    # copy the individual
    data_sub[cur_ind,]=data[k,]
    # change the population labeling so it will be from 1 to length(subpop)
    data_sub[cur_ind,2]=which(subpop==data[k,2])
    # incerment the counter
    cur_ind=cur_ind+1;
  }
}
write.table(data_sub,output,row.names = F, col.names=F)
}   

# In the case of a BayeScan band intensity matrix input file, this
# function return a list of monomorphic markers, ready to be used as a
# "discarded list" in BayeScan. For band intensity data, the definition of
# polymorphic is more complicated, and you may need this function.
# 1) A locus can be easily considered as polymorphic if the frequence of 
# "absence of band" # is higher than "pc_mono" AND lower than 1-"pc_mono" (5% by default).
# We consider "absence of band" when the intensity is lower than "pc_absence"
# (5% by default) times the maximum intensity value at this locus. It is a good 
# idea to first use the "clean_intensity" function before applying this one to 
# discard extreme values.
# 2) The more complicated: a locus is still considered as polymorphic if the 
# frequency of "absence of band" is lower than pc_mono AND the distribution 
# of "presence of band" is bimodal (contains potentially both heterozygotes homozygotes).
# To test for bimodality we first discard a fraction pc_mode (5% by default) at the two tails 
# of the distribution to avoid extreme values creating signal for multi-modality. 
# Unimodality is tested using a DIP test and the p-value to reject is "pval" (default is 0.05)
#   
# WARNING: this function needs the package "diptest" to be installed
#
# input: the name of a file in the current directory in BayeScan format for
# AFLPs with band intensity
#
# discarded: the output name for the file containing the monomorphic loci 
# 

discard_monomorphic<-function(input,discarded,pc_mono=0.05,pc_absence=0.05,pc_mode=0.05,pval=0.05)
{
data_clean=read.table(input,colClasses="numeric")

# calculate the p-value for the dip test
# need the package "diptest" to be installed
# maximum number of individuals is 5000
# p-value close to zero means unimodality
library(diptest)
data(qDiptab)
dnqd <- dimnames(qDiptab)
nn <- as.integer(dnqd $n)
dip_pval<-function(x)
{
	n1=length(x)
	dip1=dip(x)
	# find the size just above n in the table
	n2=as.integer(dnqd$n[min(which(as.integer(dnqd$n)>=n1))])
  # adjust the dip value
	dip2=dip1*sqrt(n1)/sqrt(n2)
	infvals=which(qDiptab[nn==n2,]<dip2)
	if (length(infvals)==0)
		return(0)
	else
		return(as.numeric(dnqd$Pr[max(infvals)]))
}

# now the part to detect monoporhic loci on the cleaned data 
# vector of booleans indicating for each marker if we should keep it or not
keep=NULL
# loop over each locus
for (i in 1:(ncol(data_clean)-2)) 
{
	# we extract "band presence" values (>pc_absence times max band intensity)
	non_zero=data_clean[,i+2]>(pc_absence*max(data_clean[,i+2]))
	non_zero_data=data_clean[non_zero,i+2]
	# count missing data (-1)
	missing=sum(data_clean[,i+2]==-1)
	# count absence data
	absence=nrow(data_clean)-missing-length(non_zero_data)
	# if band absence >pc_mono and <1-pc_mono we keep the band
	if (absence/(nrow(data_clean)-missing)>pc_mono && absence/(nrow(data_clean)-missing)<1-pc_mono)
		keep[i]=TRUE
	# else if freq absence<pc_mono, perform the dip test on presence data to 
  # detect bimodality 
	else if (absence/(nrow(data_clean)-missing)<pc_mono)
	{
	  # remove the two pc_mode tails of the distribution of non_zero values
		lb=quantile(non_zero_data,pc_mode)
		hb=quantile(non_zero_data,1-pc_mode)
		non_zero_data95=non_zero_data[non_zero_data>lb & non_zero_data<hb]
		keep[i]=dip_pval(non_zero_data95)>=1-pval
		#if (dip_pval(non_zero_data95)>=1-pval) cat(i," ",file="test.txt",append=TRUE)
	}
	else keep[i]=FALSE		
}

# write discarded loci
write.table(which(!keep),discarded,row.names = F, col.names=F)
# output discarded loci as well
length(which(!keep))
}

