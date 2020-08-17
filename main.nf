#! usr/bin/env nextflow


/*
~ ~ ~ > * USER INPUT PARAMETERS 
*/
date = new Date().format( 'yyyyMMdd' )

params.traitdir  = null
params.traitfile = null
params.vcf     = null
params.p3d     = null
params.sthresh   = null
params.freqUpper = 0.05
params.minburden = 2
params.refflat   = "bin/refFlat.ws245.txt"
params.genes     = "bin/gene_ref_flat.Rda"
params.cendr_v   = "20180527"
params.e_mem   = "32"
params.eigen_mem = params.e_mem + " GB"
params.group_qtl = 1000
params.ci_size   = 150
params.fix_names = "fix"
params.help    = null
params.max_QTL_number = 5
params.max_QTL_size = 2000000


println()

/*
~ ~ ~ > * OUTPUT DIRECTORY 
*/

params.out = "Analysis_Results-${date}"

/*
~ ~ ~ > * INITIATE GENE LOCATION FILE 
*/

genes = Channel.fromPath("${params.genes}")

log.info ""
log.info "------------------------------------------"
log.info "        C. elegans GWAS pipeline "
log.info "------------------------------------------"
log.info ""


if (params.help) {
    log.info "----------------------------------------------------------------"
    log.info "                      USAGE                                     "
    log.info "----------------------------------------------------------------"
    log.info ""
    log.info "nextflow main.nf --traitfile=test_bulk --vcf=bin/WI.20180527.impute.vcf.gz --p3d=TRUE --sthresh=EIGEN # run all traits from a single file"
    log.info "nextflow main.nf --traitdir=test_bulk --p3d=TRUE --sthresh=BF # download VCF from CeNDR"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--traitfile              String                Name of file that contains phenotypes. File should be tab-delimited with the columns: strain trait1 trait2 ..."
    log.info "--vcf                    String                Name of VCF to extract variants from. There should also be a tabix-generated index file with the same name in the directory that contains the VCF. If none is provided, the pipeline will download the latest VCF from CeNDR"
    log.info "--p3d                    BOOLEAN               Set to FALSE for EMMA algortith, TRUE for EMMAx"
    log.info "----------------------------------------------------------------"
    log.info "----------------------------------------------------------------"
    log.info "Optional arguments (General):"
    log.info "--out                    String                Name of folder that will contain the results"
    log.info "--e_mem                  String                Value that corresponds to the amount of memory to allocate for eigen decomposition of chromosomes (DEFAULT = 100)"
    log.info "--cendr_v                String                CeNDR release (DEFAULT = 20180527)"
    log.info "Optional arguments (Marker):"
    log.info "--sthresh                String                Significance threshold for QTL - Options: BF - for bonferroni correction, EIGEN - for SNV eigen value correction, or another number e.g. 4"
    log.info "--group_qtl              Integer               If two QTL are less than this distance from each other, combine the QTL into one, (DEFAULT = 1000)"
    log.info "--ci_size                Integer               Number of SNVs to the left and right of the peak marker used to define the QTL confidence interval, (DEFAULT = 150)"
    log.info "Optional arguments (Burden):"
    log.info "--freqUpper              Float                 Maximum allele frequency for a variant to be considered for burden mapping, (DEFAULT = 0.05)"
    log.info "--minburden              Interger              Minimum number of strains to have a variant for the variant to be considered for burden mapping, (DEFAULT = 2)"
    log.info "--genes                  String                refFlat file format that contains start and stop genomic coordinates for genes of interest, (DEFAULT = bin/gene_ref_flat.Rda)"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info "Information describing the stucture of the input files can be located in input_files/README.txt"
    log.info ""
    log.info ""
    log.info "Flags:"
    log.info "--help                                      Display this message"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info ""
    log.info " Required software packages to be in users path"
    log.info "BCFtools               v1.9"
    log.info "plink                  v1.9"
    log.info "R-cegwas2              Found on GitHub"
    log.info "R-tidyverse            v1.2.1"
    log.info "R-correlateR           Found on GitHub"
    log.info "R-rrBLUP               v4.6"
    log.info "R-sommer               v3.5"
    log.info "R-RSpectra             v0.13-1"
    log.info "R-ggbeeswarm           v0.6.0"
    log.info "--------------------------------------------------------"    
    exit 1
} else {

log.info ""
log.info "Phenotype Directory                     = ${params.traitdir}"
log.info "VCF                                     = ${params.vcf}"
log.info "CeNDR Release                           = ${params.cendr_v}"
log.info "P3D                                     = ${params.p3d}"
log.info "Significance Threshold                  = ${params.sthresh}"
log.info "Max AF for Burden Mapping               = ${params.freqUpper}"
log.info "Min Strains with Variant for Burden     = ${params.minburden}"
log.info "Significance Threshold                  = ${params.sthresh}"
log.info "Gene File                               = ${params.genes}"
log.info "Result Directory                        = ${params.out}"
log.info "Eigen Memory allocation                 = ${params.eigen_mem}"
log.info ""
}

/*
~ ~ ~ > * COMBINE VCF AND VCF INDEX INTO A CHANNEL
*/

if (params.vcf) {
  
  vcf = Channel.fromPath("${params.vcf}")

  vcf_index = Channel.fromPath("${params.vcf}" + ".tbi")

  vcf
    .spread(vcf_index)
    .into{vcf_to_whole_genome;
        vcf_to_fine_map;
        vcf_to_burden;
        vcf_to_query_vcf}

} else {

  process pull_vcf {

    tag {"PULLING VCF FROM CeNDR"}
    executor 'local'

    output:
      file("*.vcf.gz") into dl_vcf
      file("*.vcf.gz.tbi") into dl_vcf_index

    """
      wget https://storage.googleapis.com/elegansvariation.org/releases/${params.cendr_v}/variation/WI.${params.cendr_v}.impute.vcf.gz
      tabix -p vcf WI.${params.cendr_v}.impute.vcf.gz
    """
  }

  dl_vcf
    .spread(dl_vcf_index)
    .into{vcf_to_whole_genome;
        vcf_to_fine_map;
        vcf_to_burden;
        vcf_to_query_vcf}

}

/*
~ ~ ~ > * INITIATE MAPPING QTL GROUPING PARAMETER
*/

Channel
  .from("${params.ci_size}")
  .set{qtl_ci_size}

/*
~ ~ ~ > * INITIATE MAPPING QTL CONFIDENCE INTERVAL SIZE PARAMETER
*/

Channel
  .from("${params.group_qtl}")
  .set{qtl_snv_groupinng}

/*
~ ~ ~ > * INITIATE MAPPING METHOD CHANNEL
*/

Channel
  .from("${params.p3d}")
  .into{p3d_full;
      p3d_fine}

/*
~ ~ ~ > * INITIATE THRESHOLD CHANNEL
*/

Channel
  .from("${params.sthresh}")
  .into{sig_threshold_full;
      sig_threshold_fine}

/*
~ ~ ~ > * INITIATE PHENOTYPE CHANNEL - GENERATES A [trait_name, trait_file] TUPLE
*/



Channel
  .fromPath("${params.traitfile}")
  .set{ traits_to_strainlist }

process fix_strain_names_bulk {

  executor 'local'

  tag {"BULK TRAIT"}

  input:
    file(phenotypes) from traits_to_strainlist

  output:
    file("pr_*.tsv") into fixed_strain_phenotypes
    file("Phenotyped_Strains.txt") into phenotyped_strains_to_analyze

  """
    Rscript --vanilla `which Fix_Isotype_names_bulk.R` ${phenotypes} ${params.fix_names}
  """

}

fixed_strain_phenotypes
    .flatten()
    .map { file -> tuple(file.baseName.replaceAll(/pr_/,""), file) }
  .into{ traits_to_map;
      traits_to_burden }


phenotyped_strains_to_analyze
  .into{strain_list_genome;
      strain_list_finemap}

process vcf_to_geno_matrix {

  executor 'local'

  publishDir "${params.out}/Genotype_Matrix", mode: 'copy'

  cpus 1

  input:
    set file(vcf), file(index) from vcf_to_whole_genome
    file(strains) from strain_list_genome

  output:
    file("Genotype_Matrix.tsv") into geno_matrix

  """

    bcftools view -S ${strains} ${vcf} |\\
    bcftools filter -i N_MISSING=0 -Oz -o Phenotyped_Strain_VCF.vcf.gz

    tabix -p vcf Phenotyped_Strain_VCF.vcf.gz

    plink --vcf Phenotyped_Strain_VCF.vcf.gz \\
      --snps-only \\
      --biallelic-only \\
      --maf 0.05 \\
      --set-missing-var-ids @:# \\
      --indep-pairwise 50 10 0.8 \\
      --geno \\
      --allow-extra-chr

    awk -F":" '\$1=\$1' OFS="\\t" plink.prune.in | \\
    sort -k1,1d -k2,2n > markers.txt

    bcftools query -l Phenotyped_Strain_VCF.vcf.gz |\\
    sort > sorted_samples.txt 

    bcftools view -v snps \\
    -S sorted_samples.txt \\
    -R markers.txt \\
    Phenotyped_Strain_VCF.vcf.gz |\\
    bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' |\\
    sed 's/[[# 0-9]*]//g' |\\
    sed 's/:GT//g' |\\
    sed 's/0|0/-1/g' |\\
    sed 's/1|1/1/g' |\\
    sed 's/0|1/NA/g' |\\
    sed 's/1|0/NA/g' |\\
    sed 's/.|./NA/g'  |\\
    sed 's/0\\/0/-1/g' |\\
    sed 's/1\\/1/1/g'  |\\
    sed 's/0\\/1/NA/g' |\\
    sed 's/1\\/0/NA/g' |\\
    sed 's/.\\/./NA/g' > Genotype_Matrix.tsv

  """

}

geno_matrix
  .into{eigen_gm;
      mapping_gm}


/*
============================================================
~ > *                                                  * < ~
~ ~ > *                                              * < ~ ~
~ ~ ~ > *  EIGEN DECOMPOSITION OF GENOTYPE MATRIX  * < ~ ~ ~
~ ~ > *                                              * < ~ ~
~ > *                                                  * < ~
============================================================
*/

CONTIG_LIST = ["I", "II", "III", "IV", "V", "X"]
contigs = Channel.from(CONTIG_LIST)

/*
------------ Decomposition per chromosome
*/

process chrom_eigen_variants {

  tag { CHROM }

  cpus 1
  memory params.eigen_mem

  input:
    file(genotypes) from eigen_gm
    each CHROM from contigs

  output:
    file("${CHROM}_independent_snvs.csv") into sig_snps_geno_matrix
    file(genotypes) into concat_geno


  """
    cat Genotype_Matrix.tsv |\\
    awk -v chrom="${CHROM}" '{if(\$1 == "CHROM" || \$1 == chrom) print}' > ${CHROM}_gm.tsv
    Rscript --vanilla `which Get_GenoMatrix_Eigen.R` ${CHROM}_gm.tsv ${CHROM}
  """

}

/*
------------ Sum independent tests for all chromosomes
*/

process collect_eigen_variants {

  executor 'local'

  publishDir "${params.out}/Genotype_Matrix", mode: 'copy'

  cpus 1

  input:
    file(chrom_tests) from sig_snps_geno_matrix.collect()

  output:
    file("total_independent_tests.txt") into independent_tests

  """
    cat *independent_snvs.csv |\\
    grep -v inde |\\
    awk '{s+=\$1}END{print s}' > total_independent_tests.txt
  """

}

independent_tests
  .spread(mapping_gm)   // Genotype_Matrix.tsv
  .spread(traits_to_map)   // set val(TRAIT), file(trait_PC1.tsv) single trait file
  .spread(p3d_full)  // p3d value
  .spread(sig_threshold_full)
  .spread(qtl_snv_groupinng)
  .spread(qtl_ci_size)
  .set{mapping_data}

/*
======================================
~ > *                            * < ~
~ ~ > *                        * < ~ ~
~ ~ ~ > *  RUN GWAS MAPPING  * < ~ ~ ~
~ ~ > *                        * < ~ ~
~ > *                            * < ~
======================================
*/

/*
------------ Genome-wide scan
*/

process rrblup_maps {

  cpus 1

  tag { TRAIT }

  publishDir "${params.out}/Mappings/Data", mode: 'copy', pattern: "*processed_mapping.tsv"
//  publishDir "${params.out}/Mappings/Plots", mode: 'copy', pattern: "*.pdf"
  

  input:
  set file("independent_snvs.csv"), file(geno), val(TRAIT), file(pheno), val(P3D), val(sig_thresh), val(qtl_grouping_size), val(qtl_ci_size) from mapping_data

  output:
  file("*raw_mapping.tsv") into raw_map
  set val(TRAIT), file(geno), file(pheno) into processed_map_to_ld
  file("*processed_mapping.tsv") into processed_map_to_summary_plot
  set val(TRAIT), file("*processed_mapping.tsv") into pr_maps_trait
//  file("*.pdf") into gwas_plots

  """
    tests=`cat independent_snvs.csv | grep -v inde`

    Rscript --vanilla `which Run_Mappings.R` ${geno} ${pheno} ${task.cpus} ${P3D} \$tests ${sig_thresh} ${qtl_grouping_size} ${qtl_ci_size}

    if [ -e Rplots.pdf ]; then
        rm Rplots.pdf
    fi
  """
}

/*
------------ Generate GWAS QTL summary plot
*/

// need to run this first to find significant traits, 
// but then you lose track of it when joining channels below this process and therefore need to re run same step to find all peaks
process summarize_maps {

  publishDir "${params.out}/Mappings/Data", mode: 'copy', pattern: "*_peaks.tsv"

  input:
  file(maps) from processed_map_to_summary_plot.collect()

  output:
  file("QTL_peaks.tsv") into qtl_peaks


  """
    cat *processed_mapping.tsv |\\
    awk '\$0 !~ "\\tNA\\t" {print}' |\\
    awk '!seen[\$2,\$5,\$12,\$13,\$14]++' |\\
    awk 'NR>1{print \$5, \$2, \$12, \$13, \$14}' OFS="\\t" > QTL_peaks.tsv

    sig_maps=`wc -l QTL_peaks.tsv | cut -f1 -d' '`

    if [ \$sig_maps = 0 ]; then   # this is for all traits combined
      max_log10=`cat *processed_mapping.tsv | awk 'BEGIN {max = 0} {if (\$4>max && \$4!= "log10p") max=\$4} END {print max}'`
      echo "NO TRAITS HAD SIGNIFICANT MAPPINGS - MAXIMUM -log10p IS \$max_log10 - CONSIDER SETTING BF THRESHOLD BELOW THIS VALUE"
      exit
    fi
  """
}


qtl_peaks  // QTL_peaks.tsv contains all trait
   .splitCsv(sep: '\t')
   .into{peaks;printpeaks}

// peaks    // QTL_peaks.tsv split into individual rows of TRAIT, chrom, start, peak, end. All below join by TRAIT
//   .join(processed_map_to_ld)   // set val(TRAIT), file(TRAIT_genotype_matrix.tsv) file(TRAIT_trait.tsv)
processed_map_to_ld      // set val(TRAIT), file(TRAIT_genotype_matrix.tsv) file(TRAIT_trait.tsv)
   .join(pr_maps_trait)   // TRAIT_processed_mapping.tsv
   .spread(vcf_to_fine_map)   // set cendr.vcf.gz & index
   .spread(strain_list_finemap)   // file("Phenotyped_Strains.txt")
   .into{QTL_peaks; QTL_peaks_print}

/*
======================================
~ > *                            * < ~
~ ~ > *                        * < ~ ~
~ ~ ~ > *  RUN FINE MAPPING  * < ~ ~ ~
~ ~ > *                        * < ~ ~
~ > *                            * < ~
======================================
*/

/*
------------ Extract QTL interval genotype matrix
*/

process prep_ld_files {

  tag {TRAIT}

//  errorStrategy 'ignore'

  input:
    set val(TRAIT), file(complete_geno), file(phenotype), file(pr_map), file(vcf), file(index), file(strains) from QTL_peaks

  output:
    set val(TRAIT), file(complete_geno), file(phenotype), file(pr_map), file(vcf), file(index), file("*ROI_Genotype_Matrix.tsv"), file("*LD.tsv") optional true into LD_files_to_plot

  """
    echo "HELLO"
    cat ${pr_map} |\\
    awk '\$0 !~ "\\tNA\\t" {print}' |\\
    awk '!seen[\$2,\$5,\$12,\$13,\$14]++' |\\
    awk 'NR>1{print \$2, \$5, \$12, \$13, \$14}' OFS="\\t" > ${TRAIT}_QTL_peaks.tsv

    if [[ \$(wc -l <${TRAIT}_QTL_peaks.tsv) -le ${params.max_QTL_number} ]]
    then 

      filename='${TRAIT}_QTL_peaks.tsv'
      echo Start
      while read p; do 
          chromosome=`echo \$p | cut -f1 -d' '`
          trait=`echo \$p | cut -f2 -d' '`
          start_pos=`echo \$p | cut -f3 -d' '`
          peak_pos=`echo \$p | cut -f4 -d' '`
          end_pos=`echo \$p | cut -f5 -d' '`

      if [[ \$[\$end_pos - \$start_pos] -lt ${params.max_QTL_size} ]]
      then
        cat ${phenotype} | awk '\$0 !~ "strain" {print}' | cut -f1 > phenotyped_samples.txt

        bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${vcf} \
        -S phenotyped_samples.txt |\\
        bcftools filter -i N_MISSING=0 |\\
        awk '\$0 !~ "#" {print \$1":"\$2}' > \$trait.\$chromosome.\$start_pos.\$end_pos.txt


        bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${vcf} \
        -S phenotyped_samples.txt |\\
        bcftools filter -i N_MISSING=0 -Oz -o finemap.vcf.gz

        plink --vcf finemap.vcf.gz \\
          --snps-only \\
          --maf 0.05 \\
          --biallelic-only \\
          --allow-extra-chr \\
          --set-missing-var-ids @:# \\
          --geno \\
          --make-bed \\
          --recode vcf-iid bgz \\
          --extract \$trait.\$chromosome.\$start_pos.\$end_pos.txt \\
          --out \$trait.\$chromosome.\$start_pos.\$end_pos

        nsnps=`wc -l \$trait.\$chromosome.\$start_pos.\$end_pos.txt | cut -f1 -d' '`

        plink --r2 with-freqs \\
          --allow-extra-chr \\
          --snps-only \\
          --ld-window-r2 0 \\
          --ld-snp \$chromosome:\$peak_pos \\
          --ld-window \$nsnps \\
          --ld-window-kb 6000 \\
          --chr \$chromosome \\
          --out \$trait.\$chromosome:\$start_pos-\$end_pos.QTL \\
          --set-missing-var-ids @:# \\
          --vcf \$trait.\$chromosome.\$start_pos.\$end_pos.vcf.gz

        sed 's/  */\\t/g' \$trait.\$chromosome:\$start_pos-\$end_pos.QTL.ld |\\
        cut -f2-10 |\\
        sed 's/^23/X/g' | sed 's/\\t23\\t/\\tX\\t/g' > \$trait.\$chromosome.\$start_pos.\$end_pos.LD.tsv

        bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' finemap.vcf.gz |\\
          sed 's/[[# 0-9]*]//g' |\\
          sed 's/:GT//g' |\\
          sed 's/0|0/-1/g' |\\
          sed 's/1|1/1/g' |\\
          sed 's/0|1/NA/g' |\\
          sed 's/1|0/NA/g' |\\
          sed 's/.|./NA/g'  |\\
          sed 's/0\\/0/-1/g' |\\
          sed 's/1\\/1/1/g'  |\\
          sed 's/0\\/1/NA/g' |\\
          sed 's/1\\/0/NA/g' |\\
          sed 's/.\\/./NA/g' |\\
          sed 's/^23/X/g' > \$trait.\$chromosome:\$start_pos-\$end_pos.ROI_Genotype_Matrix.tsv
      fi

      done < \$filename
    fi

  """
}

//p3d_fine
//  .spread(LD_files_to_plot)
//  .set{LD_files_to_finemap}

/*
------------ Run fine mapping
*/

process rrblup_fine_maps {


  tag {TRAIT}


  input:
    set val(TRAIT), file(complete_geno), file(phenotype), file(pr_map), file(vcf), file(index), file(roi_geno_matrix), file(roi_ld) from LD_files_to_plot


  output:
    set val(TRAIT), file(phenotype), file(roi_geno_matrix), file("*prLD_df.tsv") into concat_ld_out

  """
        for i in *ROI_Genotype_Matrix.tsv;
        do
          start_pos=`echo \$i | cut -f2 -d':' | cut -f1 -d'.' | cut -f1 -d'-'`
          end_pos=`echo \$i | cut -f2 -d':' | cut -f1 -d'.' | cut -f2 -d'-'`

          ld_file=`ls *LD.tsv | grep "\$start_pos" | grep "\$end_pos" | tr -d '\\n'`
          echo "\$ld_file"
            Rscript --vanilla `which Finemap_QTL_Intervals.R` ${complete_geno} \$i ${pr_map} \$ld_file ${task.cpus} ${params.p3d}
        done   

  """
}

/*
------------ Process fine mapping : concatenate QTL finemapping files for each each trait
*/

process concatenate_LD_per_trait {

  tag {TRAIT}

  executor 'local'

  input:
    set val(TRAIT), file(phenotype), file(roi_geno_matrix), file(processed_ld) from concat_ld_out

  output:
    set file("${TRAIT}.combined_peak_LD.tsv"), file(phenotype) into combined_ld_data

  """
    for i in *prLD_df.tsv;
    do
    chrom=`head -2 \$i | tail -1 | cut -f2`
    start_pos=`echo \$i | rev | cut -f3 -d"_" | cut -f2 -d'.' | rev`
    end_pos=`echo \$i | rev | cut -f3 -d"_" | cut -f1 -d'.' | rev`

    echo "${TRAIT} \$chrom \$start_pos \$end_pos"

    awk 'BEGIN{OFS="\\t"};FNR == 1 {next};\$1=\$1' OFS="\\t" \$i |\
    awk -v trait="${TRAIT}" \
      -v chrom="\$chrom" \
      -v start_pos="\$start_pos" \
      -v end_pos="\$end_pos" \
      '{print \$0, trait, start_pos, end_pos}' OFS="\t" > ${TRAIT}.\$chrom.\$start_pos.\$end_pos.pr_ld.tsv
    done

    cat *pr_ld.tsv |\
    awk 'BEGIN{OFS="\\t"; print "marker", "CHROM", "POS", "log10p", "peak_marker", "peak_maf", "maf_marker_b", "ld_r2", "trait", "start_pos", "end_pos"}; 
      {print \$0}' > ${TRAIT}.combined_peak_LD.tsv
  """

}

genes
  .spread(combined_ld_data)
  .spread(vcf_to_query_vcf)
  .set{plot_fine_map_genes}

/*
------------ Plot Fine mapping results overlaid on gene locations
*/

process plot_genes {

  cpus 1
  memory '16 GB'

  tag {phenotype}

  publishDir "${params.out}/Fine_Mappings/Data", mode: 'copy', pattern: "*snpeff_genes.tsv"

  input:
    set file(genes), file(ld), file(phenotype), file(vcf), file(vcfindex) from plot_fine_map_genes

  output:
    file("*snpeff_genes.tsv") into gene_plts

  """
    Rscript --vanilla `which plot_genes.R` ${ld} ${phenotype} ${genes} ${vcf} ${params.cendr_v}
  """
}

