#! usr/bin/env nextflow


/*
~ ~ ~ > * USER INPUT PARAMETERS 
*/
date = new Date().format( 'yyyyMMdd' )

params.traitfile = null
params.vcf 		 = null
params.R_libpath = "/projects/b1059/software/R_lib_3.6.0"
params.help 	 = null
params.cegwas2dir = null
params.transcripteQTL = null
params.transcript_exp = null

println()

/*
~ ~ ~ > * OUTPUT DIRECTORY 
*/

params.out = "${params.cegwas2dir}/mediation-${date}"




log.info ""
log.info "------------------------------------------"
log.info "  C. elegans mediation of GWAS pipeline   "
log.info "------------------------------------------"
log.info ""

log.info ""
log.info "Result Directory                        = ${params.out}"
log.info "cegwas2-nf result directory             = ${params.cegwas2dir}"
log.info "eQTL with var_exp ≥ 0.05                = ${params.transcripteQTL}"
log.info "Input expression data of eQTL calling   = ${params.transcript_exp}"
log.info ""




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
		Rscript --vanilla `which Fix_Isotype_names_bulk_new.R` ${phenotypes} ${params.vcf}
	"""

}

fixed_strain_phenotypes
    .flatten()
    .map { file -> tuple(file.baseName.replaceAll(/pr_/,""), file) }
	.into{ traits_to_map;
		  traits_to_burden;
		  traits_to_mediate }






/*
=====================================
~ > *                           * < ~
~ ~ > *                       * < ~ ~
~ ~ ~ > *      Mediation    * < ~ ~ ~
~ ~ > *                       * < ~ ~
~ > *                           * < ~
=====================================
*/


params.qpeak = "${params.cegwas2dir}/Mappings/Data/QTL_peaks.tsv"

Channel
	.fromPath("${params.qpeak}")
	.splitCsv(sep: '\t')
   	.into{peaks;printpeaks; mediate_peaks}



mediate_peaks
.combine(traits_to_mediate, by: 0)
.into{medQTL_peaks; medQTL_peaks_print}


/*
~ ~ ~ > * INITIATE eQTL CHANNEL 
*/

File transcripteqtl_all = new File("${params.transcripteQTL}")
transcripteqtl_all_handle = transcripteqtl_all.getAbsolutePath()




process mediate_data {

	executor 'local'

	tag {TRAIT}


	input:
		set val(TRAIT),val(tch),val(tstart),val(tpeak),val(tend),file(t_file) from medQTL_peaks

	output:
		set val(TRAIT),val(tch),val(tpeak),val(tstart),val(tend), file("${TRAIT}_scaled_mapping.tsv"),file("${TRAIT}_${tch}_${tpeak}_eqtl.tsv") into QTL_phe

		file("${TRAIT}_${tch}_${tpeak}_elist.tsv") 

	"""


    Rscript --vanilla `which mediate.R` ${TRAIT} ${t_file} ${tch} ${tstart} ${tend} ${tpeak} ${transcripteqtl_all_handle}

	"""
}








/*
~ ~ ~ > * INITIATE expression data CHANNEL 
*/

Channel
	.fromPath("${params.transcript_exp}")
	.set{ transcript_exp_file }




params.genoMatri = "${params.cegwas2dir}/Genotype_Matrix/Genotype_Matrix.tsv"

Channel
	.fromPath("${params.genoMatri}")
	.set{ med_gm }


QTL_phe
.spread(med_gm)
.spread(transcript_exp_file)
.into{Gpeak_egene; Gpeak_egene_print}





process multi_mediate {

	
	module = 'R/3.6.0'
	executor 'local'

	cpus 1
	memory '2 GB'

    tag {"${TRAIT}_${tch}_${tpeak}"}



	errorStrategy 'retry'

	input:
		set val(TRAIT),val(tch),val(tpeak), val(tstart),val(tend), file(pheno), file(tr_eqtl), file(geno), file(texpression) from Gpeak_egene


	output:
		file("${TRAIT}_${tch}_${tpeak}_med.tsv") into result_mediate

		

	"""
    Rscript --vanilla `which multi_mediation.R` ${geno} ${texpression} ${pheno} ${tch} ${tpeak} ${TRAIT} ${tr_eqtl}

	"""
}










params.WS_t2g = "/projects/b1059/projects/Gaotian/resource/WS276/WS276_t2g_all.tsv"

File t2g = new File("${params.WS_t2g}")
t2g_handle = t2g.getAbsolutePath()




process multi_summarize_mediate {

	executor 'local'

	publishDir "${params.out}/mediate/summary", mode: 'copy', pattern: "*.tsv"
	publishDir "${params.out}/mediate/plot", mode: 'copy', pattern: "*.png"


	input:
	file(medir) from result_mediate.collect()

	output:
	file("*.tsv") into summarized_medi
	file("*.png") into plot_medi

	"""
	cat *med.tsv |\\
		sort |\\
		uniq  > mediation_analysis.tsv


	Rscript --vanilla `which Summarize_multi_Mediation.R` mediation_analysis.tsv ${t2g_handle}
		
		

	"""
}










/*
=====================================
~ > *                           * < ~
~ ~ > *                       * < ~ ~
~ ~ ~ > *  GENERATE REPORT  * < ~ ~ ~
~ ~ > *                       * < ~ ~
~ > *                           * < ~
=====================================
*/

workflow.onComplete {

    summary = """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]

    """

    println summary

    def outlog = new File("${params.out}/log.txt")
    outlog.newWriter().withWriter {
        outlog << param_summary
        outlog << summary
    }

    // mail summary
    if (params.email) {
        ['mail', '-s', 'cegwas2-nf', params.email].execute() << summary
    }


}

