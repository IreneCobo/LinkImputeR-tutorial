# LinkImputeR

## Table of Contents
1. [Introduction](#intro)
2. [SNP Quality Filtering](#filt)
3. [Running LinkImputeR](#run)
    - [Control file in ini format](#control)
    - [Run LinkImputeR in accuracy mode](#accuracy)
    - [Look at the accuracy results and decide which filters you wish to use for the final imputation.](#look)
    - [Perform the final imputation](#imputation)



<a name="intro"></a>
## 1. Introduction  

**SNP quality control filtering** is a key step previous to further data analysis, such as **Genome-Wide Association Analysis (GWAS)**, since it **removes** markers/individuals likely containing **genotyping errors**, ensuring the accuracy of the results. A standard SNP quality filtering for GWAS usually consists on the following steps:

- **Missigness per individual:** Removing individuals having more than a given percentage of missing data.
- **Minor allele count per marker:** Removing SNPs with a minor allele count lower than a given value.
- **Minimum quality score:** Removing SNPs with a quality score lower than a given value.
- **Minimum reads per marker (Depth):** Removing SNPs with less than a given number of reads.
- **Missigness per marker:** Removing SNPs with more than a given percentage of missing individuals per genotype.
- **Minor allele frequency:** Removing SNPs with a minor allele frequency lower than a given value.
- **Mendelian errors (for family-based data only):** Discarding families with more than a given frequency of Mendel errors (considering all SNPs) and also excluding SNPs with more than a given frequency of Mendelian error rate.

However, **deciding the most suitable quality control filtering thresholds for your dataset can be tricky**, since they depend on multiple factors. On the one hand, an **excesively lax quality filtering threshold** can **reduce the quality of the SNP dataset**, and thus, the accuracy of the results, as well as innecessarily increase the computational time. On the other hand, a **too strict threshold** can lead to a **loss of important information**, since it can lead to the removal of good quality SNP data. 

As introduced in the previous Journal Club, **imputation of lacking SNPs is another important step in GWAS**. It consists on the **inference of the lacking SNPs** on the dataset to increase the association power. The imputed genotypes expand the set of SNPs that can be tested for association, and this more comprehensive view of the genetic variation in a study can **enhance true association signals and facilitate meta-analysis**. 

**The selection of suitable quality filtering thresholds is also benefitial for this imputation step**, since imputation methods most often make use only of genotypes that are successfully inferred after having passed these quality filtering threshold.

**Most existing genotype imputation methods use patterns from known genotypes to impute missing genotypes**. For model species, this imputation step is usually performed by combining a reference panel of individuals genotyped at a dense set of polymorphic sites (usually single-nucleotide polymorphisms, or “SNPs”) with a study sample collected from a genetically similar population and genotyped at a subset of these sites. This imputation method predicts unobserved genotypes in the study sample by using a population genetic model to extrapolate allelic correlations measured in the reference panel. To this end, markers from both reference and study panel must be ordered or phased. Consequently, **this imputation method can be tricky to implement in non-model species, since requires high-quality reference genomes**. 

**LinkImputeR** [(Money *et al*. 2017)](https://doi.org/10.1186/s12864-017-3873-5) **is a program that tests for different filter parameters to perform data quality filtering, in order to maximize the quantity and quality of the resulting SNPs, while maintaining accuracy/correlation values**. As a result, LinkImputeR provides a series of combinations of thresholds or "Cases", so that users can decide on thresholds that are most suitable for their purposes. As a consequence, it improves the power of downstream analysis. Once the best SNP quality filtering threshold (aka Case) is selected, LinkImpute uses it to perform imputation of lacking SNPs. This imputation step is specifically designed for non-model organisms since it requires neither ordered markers nor a reference panel of genotypes. 

To this end, **LinkImputeR** uses an imputation method called **LD-kNN Imputation**, which is based on looking for SNPs in **Linkage Disequilibrium (LD)** with the lacking SNP to be imputed. Thus, to impute a genotype at SNP a in sample *b*, LD-kNNi first uses the *l* SNPs most in LD with the SNP to be imputed in order to calculate a distance from sample *b* to every other sample for SNP *a*. The algorithm proceeds by picking the *k* nearest neighbours to *b* that have an inferred genotype at SNP *a* and then scoring each of the possible genotypes, *c g* , as a weighted count of these genotypes.

The data quality filters tested and implemented in LinkImputerR are **depth, minor allele frequency, and missingness by both SNP and sample**. 

Here, we are using a dataset containing **SNP data for 100 Green ash (*Fraxinus pennsylvanica*) individuals belonging to 8 families**. Some of these families are the result of the crossing of two parentals showing a "resistant" or **"lingering" phenotype to the invasive pest Emerald Ash Borer**, whereas others result from crossings between two **sensitive** parentals. Individuals were sequenced using ddRAD and SNPs were called using FreeBayes. A total of **3,565,012 SNPs** were obtained after performing SNP calling.

<a name="filt"></a>
## 2. SNP Quality Filtering

Since LinkImputeR only tests and implements Minor allele frequency, depth and missigness by SNP and sample, it is necessary to perform a SNP quality filtering step for the rest of the quality filtering parameters before running LinkImputeR. 

1. Minor allele count per marker and minimum quality score using *vcftools*

```
vcftools --vcf SNP.vcf --mac 2 --minQ 30 --recode --recode-INFO-all > SNP.mac20.minQ30.vcf
```
- This code discards SNPs with less than 2 counts and a quality score lower than 30
- **Output: 1,325,113 SNPs** after filtering

2. Custom Mendelian Filter: Since it is a family-based study, Mendelian filtering was performed using a python script. However, this filtering can be also performed using *plink* program:

```
plink --vcf SNP.mac20.minQ30.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --me 0.05 0.1 --recode vcf --out SNP.mac20.minQ30.mendelian.vcf
``` 
This code discards families with more than 5 % of Mendelian errors (considering all SNPs) and excludes SNPs with more than a 10 % of Mendelian errors. 

<a name="run"></a>
## 3. Running LinkImputeR (Simple)
In a simple case running LinkImputeR will consist of four steps:
- 1. Create a control file in ini format. The control file includes information on input and output files and the filters to be tested.
- 2. Run LinkImputeR in accuracy mode. This produces data set size and accuracy statistics for each set of filters to be tested. 
- 3. Look at the accuracy results and decide which filters you wish to use for the final imputation. 
- 4. Perform the final imputation. 

<a name="control"></a>
### 1. Control file in *ini* format

LinkImpute makes use of a control file in *.ini* format. This control file includes information on input and output files and the filters to be tested. 

accuracy.ini

```
[Input]
filename = SNP.mac20.minQ30.mendelian.vcf # LinkImputeR takes a standard VCF file as input. LinkImputeR filters input data for biallelic SNPs before performing any other steps. The input file can be gzipped and this should be detected automatically.
save = SNP.filtered.vcf

[Global] #Here are included global parameters that can effect multiple filters
depth = 2,4 #The minimum depth used to call a genotype for use in the filters.

[InputFilters] #These filters are applied before any accuracy run is performed and hence are only performed once no matter how many other filter cases are being tested. In the case of filters also included under CaseFilters it is worthwhile including the least stringent of the cases here (in the case of positionmissing the least stringent of the missing cases in CaseFilters) as this will improve performance.
maf=0.01 #Minor Allele Frequency filter
positionmissing = 0.9 #Position Missing filter. Value is the maximum missingness allowed per SNP.

[CaseFilters] #This section lists the different filters to be tested, for example you may wish to test multiple MAF filters with differing values for the allowed MAF. For each filter the different values to be tested for that filter are separated by commas (this is best illustrated by looking at the example ini file). The test cases then consist of every possible combination of the different filters.
missing = 0.7,0.8,0.9 #Missing filter. Value is the maximum missingness allowed per SNP and sample.
maf=0.01,0.02 #Minor Allele Frequency filter. Value is the minimum MAF for the SNP to be included.

[Stats] #This section defines how the accuracy results should be outputted.
root = ./ #The root directory to which accuracy statistics should be saved.
level = sum #How much output should be written. Options are sum, pretty and table. sum will just produce a summary file, pretty will output extra files including information in a easy human readable way and table will also output the same information in a tab-delimited format.
eachmasked = no

[Output] #Where should the output be written
control = ./impute.xml #The file name of the output control file which will be used in the imputation stage of LinkImputeR.

[Log] #This section controls logging
file = log.txt #The file name of the log file.
level = brief #How much should be logged. Options are, in order of increasing output, critical (default), brief, detail and debug.

[Accuracy] #To estimate accuracy read counts from ‘known’ genotypes are masked at random from across the dataset without replacement. We consider a genotype to be known if it has a read depth ≥30. Accuracy is then defined as the proportion of masked genotypes where the ‘known’ and called genotypes are the same.
numbermasked=10000 #The number of genotypes to mask. Default is 10000
```

<a name="accuracy"></a>
### 2. Run LinkImputeR in accuracy mode

```
java -jar LinkImputeR.jar -s accuracy.ini
```

<a name="look"></a>
### 3. Look at the accuracy results and decide which filters you wish to use for the final imputation

```
Name	Samples	Positions	Accuracy	Correlation	Filters	Additional
Case 1	100	343652	0.9684	0.9428	PositionMiss(0.7),SampleMiss(0.7),MAF(0.01)	Depth(2)
Case 2	100	422531	0.9669	0.9425	PositionMiss(0.8),SampleMiss(0.8),MAF(0.01)	Depth(2)
Case 3	100	930660	0.9558	0.9189	PositionMiss(0.9),SampleMiss(0.9),MAF(0.01)	Depth(2)
Case 4	100	325977	0.9639	0.9405	PositionMiss(0.7),SampleMiss(0.7),MAF(0.02)	Depth(2)
Case 5	100	404026	0.9618	0.9360	PositionMiss(0.8),SampleMiss(0.8),MAF(0.02)	Depth(2)
Case 6	100	907122	0.9494	0.9095	PositionMiss(0.9),SampleMiss(0.9),MAF(0.02)	Depth(2)
Case 7	100	252294	0.9788	0.9578	PositionMiss(0.7),SampleMiss(0.7),MAF(0.01)	Depth(4)
Case 8	100	316696	0.9750	0.9464	PositionMiss(0.8),SampleMiss(0.8),MAF(0.01)	Depth(4)
Case 9	100	808533	0.9616	0.9244	PositionMiss(0.9),SampleMiss(0.9),MAF(0.01)	Depth(4)
Case 10	100	216859	0.9740	0.9527	PositionMiss(0.7),SampleMiss(0.7),MAF(0.02)	Depth(4)
Case 11	100	277145	0.9736	0.9484	PositionMiss(0.8),SampleMiss(0.8),MAF(0.02)	Depth(4)
Case 12	100	761913	0.9571	0.9218	PositionMiss(0.9),SampleMiss(0.9),MAF(0.02)	Depth(4)
```
Case 3 yields the highest number of SNPs (930,660) while maintaining a high accuracy (0.9189), so it was selected to perform the final imputation step. By performing the standard data quality filtering for GWAS  without using LinkImputeR, we obtained 29,833 SNPs out of a original dataset of 3,565,012 SNPs. 

<a name="imputation"></a>
### 4. Perform the final imputation
```
java -jar LinkImputeR.jar impute.xml 'Case 3' imputed.vcf
```
Where *impute.xml* is the name of the control file asked to be created in the ini file, *Case 3* is the name of the case you want to do imputation for and *imputed.vcf* is the name of the imputed vcf file. CASE will have to be in quotes since the default case names include spaces. If the OUTPUT ends in .gz then the file will be gzipped. 

