TABLE OF CONTENTS
1. INTRODUCTION
2. BACKGROUND
3. INPUT FILE FORMATS
4. PIPELINE ORDER
4a. PEDIGREE TRACING
4b. PEDIGREE IMPUTATION
4c. MULTI-GENERATION TRACING
4d. RECOMBINATION IDENTIFICATION
4e. GENOMIC PREDICTION
5. CHANGELOG

1. INTRODUCTION
PATRIOT, or Parental Allele Tracing, Recombination Identification, and Optimal predicTion and selection, is an R-language pipeline designed to track the inheritance of genome-wide markers from parents to progeny in inbred species. 

2. BACKGROUND
PATRIOT was developed with the initial target of tracking inheritance of chromosomal segments in soybean (Glycine max) based on publicly available genotyping data from two fixed-marker genotyping platforms (SoySNP50k and SoySNP6k chips). The initial codebase was written in R in four steps, in order to allow partial usage based on user needs. The four steps consist of identifying markers which could only have been inherited from a single parent, followed by imputation based on flanking markers of known inheritance, an extension to allow multiple generations of tracing, and genomic prediction based on shared inheritance of markers from an ancestral source.

3. INPUT FILE FORMATS
PATRIOT utilizes two input files to identify inheritance of markers. 
The first is a simplified pedigree file in .csv format, with column headers "Geno", "F", and "M". Each row should contain the simplified pedigree for a single progeny, with the progeny name in column "Geno", the female parent in column "F", and the male parent in column "M". The order of parents does not matter, and any backcross-derived material should have parents presented as though the originated from a single cross. A short example is provided immediately below, where Pro# corresponds to the name of the nth progeny and Par# corresponds to the name of the nth parental line.
#Pedigree file format#
Geno	F	M
Pro1	Par1	Par2
Pro2	Par3	Par4
Pro3	Par3	Par4
#End Pedigree file format#
The second input file for marker inheritance identification is the genotypic data for parents and progeny. This can be stored as .txt or .csv and should adhere to the following general format:
#Genotypic data format#
SNP	Chr	Pos	Par1	Par2	Pro 1	...
SNP1	1	1235	A	T	A	...
SNP2	1	1246	C	G	C	...
SNP3	1	1256	T	C	H	...
...	...	...	...	...	...	...
#End Genotypic data format#
The first three columns should not differ from the presented format, with the exception that marker position may be expressed in map distance (cM) or genetic distance (bp/kbp). The input file should be sorted such that the markers are sorted by chromosome and increasing position. Failure to pre-sort may result in incorrect imputation calls. For the marker data, no specific order of parents or progeny is needed, and marker data may be in A/T/C/G/H format or integer (0/1/2) format. As of present, heterozygotes are not specifically handled, and will be imputed to a single parent if possible. Please check here and/or the changelog at the end of the file for updates if this functionality is needed.

The final file utilized by PATRIOT is the phenotypic data file for genomic prediction/selection. This file should follow the general format outlined below.
#Phenotypic data format#
Env	Geno	Trait1	Trait2	...
Env1	Pro1	45.7	6.2	...
Env1	Pro2	43.8	6.4	...
Env2	Pro1	48.2	6.5	...
...	...	...	...	...
#End Phenotypic data format#
The order of environments or progeny does not matter. Genomic prediction should be possible for any numeric trait.

4. PIPELINE ORDER
The following order represents the "happy path" for use of the full PATRIOT pipeline. However, individual codes may be used on their own or out of order based on user needs.

4a. PEDIGREE TRACING
This code requires the input of a pedigree file and marker file containing genotypic data for both the parents and the progeny. For guidance on file formating, see section 3, INPUT FILE FORMATS. The output file will identify any markers which can be traced to a specific parent by replacing that marker with the parent's name, while any markers which could have come from either parent (or does not directly match either parent) will be replaced with "Par1/Par2", where Par1 and Par2 are the names of the parental lines identified in the pedigree file. For the sake of this README, we will call the former (known parent string) as "classified" markers, and the latter ("Par1/Par2" string) as "unclassified" markers from this point on. 

4b. PEDIGREE IMPUTATION
The pedigree imputation code requires only the file from the previous step as input. "Unclassified" markers (see section 4a) at the ends of a chromosome are imputed to equal the first "classified" marker encountered from that direction. In addition, "unclassified" markers not on the ends of the chromosome are imputed if they are surrounded by the same "classified" marker class (i.e. the nearest "classified" marker class on either side of the "unclassified" marker correspond to the same parent). "Unclassified" markers which are surrounded by "classified" markers of different class (i.e., Par1 on one side and Par2 on the other side) will not be imputed.

4c. MULTI-GENERATION TRACING
For users with multiple generations of genotyped material, this code will replace "classified" markers (see section 4a) with the value from the corresponding position in that parent's record. This allows programs with multiple generations of genotyped lines to trace back to previous generations, which should improve the results of genomic prediction due to the increased representation of ancestral chromosomal segments vs. limited representation of direct parents' chromosomal segments when applied to a breeding program. This information can also be useful in increasing the accuracy of a genomic relationship matrix.

4d. RECOMBINATION IDENTIFICATION
While not required for genomic prediction, a short code is provided to utilize the output from either PEDIGREE IMPUTATION or MULTI-GENERATION TRACING to identify regions of recombination.

4e. GENOMIC PREDICTION
The genomic prediction code utilizes the output file from PEDIGREE IMPUTATION or MULTI-GENERATION TRACING, as well as the Phenotypic data file, in order to generate predictions for untested (or held-out, as in the case of the original paper) lines. The difference between the phenotypic value of the population as a whole and the mean phenotype of those lines which inherited the nth marker from a specific ancestor is calculated and used to replace the raw marker data file as input into rrBLUP.

5. CHANGELOG
Version 0.1 
Initial upload of GP.R, Genomic Prediction Looped.R, GenomicSelectionEffectiveness.R, MGT.R, NAMFullPheno.csv, Pedigree Imputation.R, Pedigree Tracing.R, Recombination identification.R, SoyNAMpedigrees.csv, and datatableGPL.R.

Version 0.2 (current) 
Added code comments for Genomic Prediction Looped.R (now Code 4), Pedigree Imputation.R (now Code 2), Pedigree Tracing.R (now Code 1), MGT.R (now Code 3), and Recombination identification.R (now Code 5). 
Modified Genomic Prediction Looped.R (now Code 4) to automatically handle parent names for calculation of Allele Effect Estimate matrix (aeematrix) rather than hard-coding specific to the SoyNAM example in order to increase viability of use in other crops/programs.