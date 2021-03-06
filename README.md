# cegwas2-nf
GWA mapping with C. elegans


## Overview of the workflow

![alt text](https://github.com/AndersenLab/cegwas2-nf/blob/master/images/Cegwas2_flow_v2.png)

## Required software packages that should be in users PATH

1. [R-v3.6.0](https://www.r-project.org/)
1. [nextflow-v19.07.0](https://www.nextflow.io/docs/latest/getstarted.html)
1. [BCFtools-v1.9](https://samtools.github.io/bcftools/bcftools.html)
1. [plink-v1.9](https://www.cog-genomics.org/plink2)
1. [R-cegwas2](https://github.com/AndersenLab/cegwas2)
1. [R-tidyverse-v1.2.1](https://www.tidyverse.org/)
1. [R-correlateR](https://github.com/AEBilgrau/correlateR)
1. [R-rrBLUP-v4.6](https://cran.r-project.org/web/packages/rrBLUP/rrBLUP.pdf)
1. [R-RSpectra-v0.13-1](https://github.com/yixuan/RSpectra)
1. [R-ggbeeswarm-v0.6](https://github.com/eclarke/ggbeeswarm)
1. [R-qtl-v1.46-2](https://github.com/kbroman/qtl)
1. [R-sommer-v4.1.1](https://cran.r-project.org/web/packages/sommer/index.html)
1. [R-genetics-v1.3.8.1.2](https://cran.r-project.org/web/packages/genetics/index.html)

## Execution of pipeline using Nextflow
```
git clone https://github.com/AndersenLab/cegwas2-nf.git
cd cegwas2-nf
nextflow main.nf --traitfile=test_traits.tsv --vcf=bin/WI.20180527.impute.vcf.gz --p3d=TRUE --sthresh=BF
```
### Parameters

* `nextflow main.nf --help` - will display the help message

#### One of the two trait parameters is required:

* `--traitfile` - is a tab-delimited formatted (.tsv) file that contains trait information.  Each phenotype file should be in the following format (replace trait_name with the phenotype of interest):

| strain | trait_name_1 | trait_name_2 |
| --- | --- | --- |
| JU258 | 32.73 | 19.34 |
| ECA640 | 34.065378 | 12.32 |
| ... | ... | ... | 124.33 |
| ECA250 | 34.096 | 23.1 |

#### Required mapping parameters
* `--p3d` - This determines what type of kinship correction to perform prior to mapping. `TRUE` corresponds to the EMMAx method and `FALSE` corresponds to the slower EMMA method. We recommend running with `--p3d=TRUE` to make sure all files of the required files are present and in the proper format, then run with `--p3d=FALSE` for a more exact mapping. Default: FALSE.

* `--sthresh` - This determines the signficance threshold required for performing post-mapping analysis of a QTL. `BF` corresponds to Bonferroni correction, `EIGEN` corresponds to correcting for the number of independent markers in your data set, and `user-specified` corresponds to a user-defined threshold, where you replace user-specified with a number. For example `--sthresh=4` will set the threshold to a `-log10(p)` value of 4. We recommend using the strict `BF` correction as a first pass to see what the resulting data looks like. If the pipeline stops at the `summarize_maps` process, no significant QTL were discovered with the input threshold. You might want to consider lowering the threshold if this occurs. 

#### Optional parameters
* `--vcf` - is a VCF file with variant data. All strains with phenotypes should be represented in the VCF used for mapping. There should also abe a tabix-generated index file (.tbi) in the same folder as the specified VCF file that has the same name as the VCF except for the addition of the `.tbi` extension. (generated using `tabix -p vcf vcfname.vcf.gz`). If this flag is not used a VCF for the C. elegans species will be downloaded from [CeNDR](https://elegansvariation.org/data/release/latest)

* `--freqUpper` - Upper bound for variant allele frequency for a variant to be considered for burden mapping. Default = 0.5

* `--minburden` - The number of strains that must share a variant for that variant to be considered for burden mapping. Default = 2

* `--refflat` - Genomic locations for genes used for burden mapping. A default generated from WS245 is provided in the repositories bin. 

* `--genes` - Genomic locations for genes formatted for plotting purposes. A default generated from WS245 is provided in the repositories bin.

* `--fix_names` - This will query the CeNDR strain set an resolve any discrepancies between your strain set and isotype names on CeNDR. This is important if you are not providing your own VCF, however if you provide your own VCF that contains the strains you phenotyped, you do not need to fix strain names (Default = "fix", change to anything but fix to skip).

### R scripts

* `Get_GenoMatrix_Eigen.R` - Takes a genotype matrix and chromosome name as input and identifies the number significant eigenvalues.
* `Fix_Isotype_names.R` - Take sample names present in phenotype data and changes them to isotype names found on [CeNDR](elegansvariation.org) when the `--traitdir` flag is used.
* `Run_Mappings.R` - Performs GWA mapping using the rrBLUP R package and the EMMA or EMMAx algorithm for kinship correction. Generates manhattan plot and phenotype by genotype plot for peak positions.
* `Summarize_Mappings.R` - Generates plot of all QTL identified in nextflow pipeline.
* `Finemap_QTL_Intervals.R` - Run EMMA/EMMAx on QTL region of interest. Generates fine map plot, colored by LD with peak QTL SNV found from genome-wide scan
* `plot_genes.R` - Runs SnpEff and generates gene plot. 
* `makeped.R` - Converts trait `.tsv` files to `.ped` format for burden mapping.
* `rvtest` - Executable to run burden mapping, can be found at the [RVtests homepage](https://github.com/zhanxw/rvtests)
* `plot_burden.R` - Plots the results from burden mapping.
* `Fix_Isotype_names_bulk.R` - Take sample names present in phenotype data and changes them to isotype names found on [CeNDR](elegansvariation.org) when the `--traitfile` flag is used.

### Output Folder Structure

```
Genotype_Matrix
  ├── Genotype_Matrix.tsv
  ├── total_independent_tests.txt
 Mappings
  ├── Data             
      ├── traitname_processed_mapping.tsv
      ├── QTL_peaks.tsv
  ├── Plots   
      ├── traitname_manplot.pdf
      ├── traitname_pxgplot.pdf
      ├── Summarized_mappings.pdf
 Fine_Mappings
  ├── Data             
      ├── traitname_snpeff_genes.tsv
  ├── Plots   
      ├── traitname_qtlinterval_finemap_plot.pdf
      ├── traitname_qtlinterval_gene_plot.pdf
 BURDEN
  ├── VT             
      ├── Data             
          ├── traitname.VariableThresholdPrice.assoc
      ├── Plots   
          ├── traitname_VTprice.pdf
  ├── SKAT   
      ├── Data             
          ├── traitname.Skat.assoc
      ├── Plots   
          ├── traitname_SKAT.pdf
```

#### Genotype_Matrix folder
* `Genotype_Matrix.tsv` - pruned LD-pruned genotype matrix used for GWAS and construction of kinship matrix
* `total_independent_tests.txt` - number of independent tests determined through spectral decomposition of the genotype matrix

#### Mappings folder

##### Data
* `traitname_processed_mapping.tsv` - Processed mapping data frame for each trait mapped
* `QTL_peaks.tsv` - List of signifcant QTL identified across all traits

##### Plots
* `traitname_manplot.pdf` - Manhattan plot for each trait that was analyzed. Two significance threshold lines are present, one for the Bonferronit corrected threshold, and another for the spectral decomposition threshold.
* `traitname_pxgplot.pdf` - Phenotype by genotype split at peak QTL positions for every significant QTL identified
* `Summarized_mappings.pdf` - A summary plot of all QTL identified

#### Fine_Mappings folder

##### Data
* `traitname_snpeff_genes.tsv` - Fine-mapping data frame for all significant QTL

##### Plots
* `traitname_qtlinterval_finemap_plot.pdf` - Fine map plot of QTL interval, colored by marker LD with the peak QTL identified from the genome-wide scan
* `traitname_qtlinterval_gene_plot.pdf` - variant annotation plot overlaid with gene CDS for QTL interval


#### BURDEN folder (Contains two subfolders VT/SKAT with the same structure)

##### Data
* `traitname.VariableThresholdPrice.assoc` - Genome-wide burden mapping result using VT price, see [RVtests homepage](https://github.com/zhanxw/rvtests)
* `traitname.Skat.assoc` - Genome-wide burden mapping result using Skat, see [RVtests homepage](https://github.com/zhanxw/rvtests)

##### Plots
* `traitname_VTprice.pdf` - Genome-wide burden mapping manhattan plot for VTprice
* `traitname_SKAT.pdf` - Genome-wide burden mapping manhattan plot for Skat


# Mediation analysis

## Required software packages that should be in users PATH
1. [R-v3.6.0](https://www.r-project.org/)
1. [R-MultiMed-v3.12](https://www.bioconductor.org/packages/release/bioc/html/MultiMed.html)
1. [R-tidyverse-v1.2.1](https://www.tidyverse.org/)

## different nf
* `main_mediation.nf` - Mapping and Mediation. 
* `mediation.nf` - Mediation on main.nf results.
* `mediation_STDcegwas2.nf` - Mediation on standard cegwas2-nf results.

#### Required mapping parameters for all nf 
* `--transcripteQTL` - eQTL peak file
* `--transcript_exp` - expression file, input of eQTL calling

#### Required mapping parameters for mediation.nf and mediation_STDcegwas2.nf
* `--traitfile` - 
* `--cegwas2dir` -

 
# Permutated_FDR.nf
1. Randomly sample 1 trait from trait file 
1. Permutate the trait for 200 times
1. Run EMMA mapping with BF and EIGEN, respectively
1. Among all the -log10P that passed the threshold, get the 5% FDR

# main_eQTL.nf
* Call expression QTL with EMMA and draw eQTL map
* `--pos` - transcript or gene positions in the genome
