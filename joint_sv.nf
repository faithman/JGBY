#!/usr/bin/env nextflow 

date = new Date().format( 'yyyyMMdd' )

/*
~ ~ ~ > * Global parameters 
*/

params.out_base               = null
params.bamdir                 = params.bamdir
params.cores                  = params.cores
params.strains                = params.strains
params.cpu                    = params.cores
params.annotation_reference   = params.annotation_reference
params.reference   = params.reference

params.out      = "${date}-${params.out_base}"

if (params.help) {
    log.info '''
_________ _______ _________ _       _________   _______             _        _______ 
\\__    _/(  ___  )\\__   __/( (    /|\\__   __/  (  ____ \\|\\     /|  ( (    /|(  ____ \\
   )  (  | (   ) |   ) (   |  \\  ( |   ) (     | (    \\/| )   ( |  |  \\  ( || (    \\/
   |  |  | |   | |   | |   |   \\ | |   | |     | (_____ | |   | |  |   \\ | || (__    
   |  |  | |   | |   | |   | (\\ \\) |   | |     (_____  )( (   ) )  | (\\ \\) ||  __)   
   |  |  | |   | |   | |   | | \\   |   | |           ) | \\ \\_/ /   | | \\   || (      
|\\_)  )  | (___) |___) (___| )  \\  |   | |     /\\____) |  \\   /    | )  \\  || )      
(____/   (_______)\\_______/|/    )_)   )_(     \\_______)   \\_/     |/    )_)|/       
                                                                                     
'''
    log.info "----------------------------------------------------------------"
    log.info "                      USAGE                                     "
    log.info "----------------------------------------------------------------"
    log.info ""
    log.info "nextflow ID_diverged-regions.nf --out_base Analysis"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--out_base             String                Name of folder to output results"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info "Optional arguments:"
    log.info "Information describing the stucture of the input files can be located in input_files/README.txt"
    log.info ""
    log.info "--cores                INTEGER              Number of cpu to use (default=2)"
    log.info "--email                STRING               email address for job notifications"
    log.info "--bamdir               String                Name of folder where bam files are"
    log.info ""
    log.info "Flags:"
    log.info "--help                                      Display this message"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info ""
    log.info " Required software packages to be in users path"
    log.info "BCFtools               v1.9"
    log.info "BEDtools               v2.27.1"
    log.info "R                      vXX"
    log.info "R-tidyverse            vXX"
    log.info "R-data.table           vXX"
    log.info "--------------------------------------------------------"    
    exit 1
} else {
        log.info '''
_________ _______ _________ _       _________   _______             _        _______ 
\\__    _/(  ___  )\\__   __/( (    /|\\__   __/  (  ____ \\|\\     /|  ( (    /|(  ____ \\
   )  (  | (   ) |   ) (   |  \\  ( |   ) (     | (    \\/| )   ( |  |  \\  ( || (    \\/
   |  |  | |   | |   | |   |   \\ | |   | |     | (_____ | |   | |  |   \\ | || (__    
   |  |  | |   | |   | |   | (\\ \\) |   | |     (_____  )( (   ) )  | (\\ \\) ||  __)   
   |  |  | |   | |   | |   | | \\   |   | |           ) | \\ \\_/ /   | | \\   || (      
|\\_)  )  | (___) |___) (___| )  \\  |   | |     /\\____) |  \\   /    | )  \\  || )      
(____/   (_______)\\_______/|/    )_)   )_(     \\_______)   \\_/     |/    )_)|/       
                                                                                     
'''
    log.info ""
    log.info "BAM Directory                       = ${params.bamdir}"
    log.info "CPUs                                = ${params.cores}"
    log.info "Output folder                       = ${params.out}"
    log.info ""
}

/*
~ ~ ~ > * Reference File 
*/

File reference = new File("${params.reference}")
if (params.reference != "(required)") {
   reference_handle = reference.getAbsolutePath();
   reference_handle_uncompressed = reference_handle.replace(".gz", "")
} else {
   reference_handle = "(required)"
}

/*
~ ~ ~ > * Reference Strain File 
*/

/*
~ ~ ~ > * Initiate BAM file channel 
*/

Channel.fromFilePairs( params.bamdir )
       .map { it -> [it[0], it[1][0], it[1][1]] }
       .into { bam_manta;
        bam_delly;
        bam_delly_recall;
        bam_smoove;
        bam_smoove_recall;
        bam_smoove_svtyper;
        bam_sv_image;
        bam_to_other }


/*
========================================
~ > *                              * < ~
~ ~ > *                          * < ~ ~
~ ~ ~ > * Run Manta SV caller * < ~ ~ ~ 
~ ~ > *                          * < ~ ~ 
~ > *                              * < ~
========================================
*/

/*
~ ~ ~ > * Call Manta
*/

process manta_call {
    
    tag { SM }

    cpus params.cores

    when:
        params.manta_path

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_manta

    output:
        file "*.vcf.gz" into individual_output_vcf_zipped
        file "*.vcf.gz.tbi" into individual_output_index
        set val(SM), file("${SM}_manta.vcf") into manta_to_db


    """
    ${params.manta_path}/configManta.py \\
        --bam ${SM}.bam \\
        --referenceFasta ${reference_handle_uncompressed} \\
        --outputContig \\
        --runDir results/

        python2 results/runWorkflow.py -m local -j ${task.cpus-1}

        cp results/results/variants/diploidSV.vcf.gz .
        cp results/results/variants/diploidSV.vcf.gz.tbi .

        mv diploidSV.vcf.gz ${SM}_manta.vcf.gz 
        mv diploidSV.vcf.gz.tbi ${SM}_manta.vcf.gz.tbi

        bcftools view ${SM}_manta.vcf.gz | \\
        bcftools filter -e 'INFO/SVLEN>100000' | \\
        bcftools filter -e 'INFO/SVLEN<-100000' | \\
        bcftools filter -i 'FILTER="PASS"' -Ov -o ${SM}_manta.vcf

    """

}

individual_output_vcf_zipped
  .toSortedList()
  .set { merged_deletion_vcf }


individual_output_index
  .toSortedList()
  .set { merged_vcf_index }

/*
~ ~ ~ > * Merge Manta
*/

process merge_manta_vcf {

    publishDir "${params.out}/variation", mode: 'copy'

    cpus params.cores

    input:
      file merged_deletion_vcf
      file merged_vcf_index

    output:
      set file("WI.MANTAsv.soft-filter.vcf.gz"), file("WI.MANTAsv.soft-filter.vcf.gz.tbi") into processed_manta_vcf
      file("WI.MANTAsv.soft-filter.stats.txt") into bcf_manta_stats

    """
        bcftools merge -m all \\
            --threads ${task.cpus-1} \\
            -Oz -o WI.MANTAsv.soft-filter.vcf.gz \\
            ${merged_deletion_vcf}
        tabix -p vcf WI.MANTAsv.soft-filter.vcf.gz

        bcftools stats --verbose WI.MANTAsv.soft-filter.vcf.gz > WI.MANTAsv.soft-filter.stats.txt
    """

}


/*
~ ~ ~ > * Process Merged Manta VCF
*/


process prune_manta {

    cpus params.cores
    
    publishDir "${params.out}/variation", mode: 'copy'

    input:
        set file(mantavcf), file(mantaindex) from processed_manta_vcf

    output:
        set file("WI.MANTAsv.LargeRemoved.snpeff.vcf.gz"), file("WI.MANTAsv.LargeRemoved.snpeff.vcf.gz.tbi") into snpeff_manta_vcf
        file("WI.MANTAsv.CONTIGS.tsv.gz") into manta_contigs
        file("MANTAsv_snpeff_out.csv") optional true into manta_snpeff_multiqc

    """
        bcftools plugin setGT -Oz -o manta_gt_filled.vcf.gz -- ${mantavcf} -t . -n 0
        bcftools query -l manta_gt_filled.vcf.gz | sort > sample_names.txt

        bcftools view --samples-file=sample_names.txt -Oz -o manta_gt_filled_sorted.vcf.gz manta_gt_filled.vcf.gz

        bcftools view manta_gt_filled_sorted.vcf.gz | \\
        bcftools filter -e 'INFO/SVLEN>100000' | \\
        bcftools filter -e 'INFO/SVLEN<-100000' | \\
        bcftools view -Oz -o manta_gt_filled_sorted_largeRemoved.vcf.gz

        bcftools view -O v manta_gt_filled_sorted_largeRemoved.vcf.gz | \\
        snpEff eff -v chicken_new \\
                 -csvStats snpeff_out.csv \\
                 -no-downstream \\
                 -no-intergenic \\
                 -no-upstream \\
                 -nodownload \\ | \\
        bcftools view -Oz -o WI.MANTAsv.LargeRemoved.snpeff.vcf.gz

        tabix -p vcf WI.MANTAsv.LargeRemoved.snpeff.vcf.gz

        bcftools query -f '%CHROM\\t%POS\\t%END\\t%SVTYPE\\t%SVLEN\\t%CONTIG[\\t%GT]\\n' WI.MANTAsv.LargeRemoved.snpeff.vcf.gz > WI.MANTAsv.CONTIGS.tsv

        bgzip WI.MANTAsv.CONTIGS.tsv
    """

}

/*
=======================================
~ > *                             * < ~
~ ~ > *                         * < ~ ~
~ ~ ~ > * Run Delly SV caller * < ~ ~ ~ 
~ ~ > *                         * < ~ ~ 
~ > *                             * < ~
=======================================
*/

/*
~ ~ ~ > * First call DELLY
*/




/*
========================================
~ > *                              * < ~
~ ~ > *                          * < ~ ~
~ ~ ~ > * Run smoove SV caller * < ~ ~ ~ 
~ ~ > *                          * < ~ ~ 
~ > *                              * < ~
========================================
*/

/*
~ ~ ~ > * First call smoove
*/

process smoove_call_sv {
    
    cpus params.cores
    
    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_smoove

    output:
        file "${SM}-smoove.genotyped.vcf.gz" into smoove_vcf
        file "${SM}-smoove.genotyped.vcf.gz.csi" into smoove_vcf_index


    """
        smoove call --outdir . --name ${SM} --fasta ${reference_handle_uncompressed} -p 1 --genotype ${SM}.bam
    """

}

smoove_vcf
    .toSortedList()
    .set { smoove_sample_vcfs }

smoove_vcf_index
    .toSortedList()
    .set { smoove_sample_vcf_ind }

/*
~ ~ ~ > * Merge first call smoove
*/

process merge_smoove_vcf {
    
    publishDir "${params.out}/variation", mode: 'copy'

    input:
        file smoove_sample_vcfs
        file smoove_sample_vcf_ind

    output:
        file("merged_smoove.sites.vcf.gz") into smoove_sites

    """
        smoove merge --fasta ${reference_handle_uncompressed} --name merged_smoove --outdir . ${smoove_sample_vcfs}
    """

}

/*
~ ~ ~ > * Recall smoove
*/

process recall_smoove {

    cpus params.cores
    
    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_smoove_recall
        file(merged_smoove) from smoove_sites

    output:
        file("${SM}-joint-smoove.genotyped.vcf.gz") into smoove_vcf_recalled
        file("${SM}-joint-smoove.genotyped.vcf.gz.tbi") into smoove_index_recalled
        set file("${SM}-joint-smoove.genotyped.vcf.gz"), file("${SM}-joint-smoove.genotyped.vcf.gz.tbi") into smoove_2_svtyper
        set val(SM), file("${SM}-joint-smoove.genotyped.vcf") into smoove_to_db


    """
        smoove genotype -x -p 1 --name ${SM}-joint --outdir . --fasta ${reference_handle_uncompressed} --vcf ${merged_smoove} ${SM}.bam

        echo "##INFO=<ID=PREND,Number=.,Type=String,Description="LUMPY probability curve of the END breakend">" > fix_header.txt

        bcftools annotate --header-lines fix_header.txt ${SM}-joint-smoove.genotyped.vcf.gz |\\
        bcftools filter -e 'INFO/SVLEN>100000' | \\
        bcftools filter -e 'INFO/SVLEN<-100000' |\\
        bcftools filter -e 'INFO/END<POS' |\\
        bcftools filter -i 'FILTER="PASS"' -Ov -o ${SM}_temp.vcf

        gatk --java-options "-Xmx4g -Xms4g" \\
             VariantFiltration \\
             -R ${reference_handle_uncompressed} \\
             --variant ${SM}_temp.vcf \\
             --genotype-filter-expression "FILTER != 'PASS'" \\
             --genotype-filter-name "VFILTER" \\
             -O temp_soft_filtered.vcf

        awk 'BEGIN{FS=OFS="\\t"} {gsub("VFILTER",\$7,\$10)} 1' temp_soft_filtered.vcf | \\
        bcftools view -Ov -o ${SM}-joint-smoove.genotyped.vcf

        tabix -p vcf ${SM}-joint-smoove.genotyped.vcf.gz
    """

}



smoove_vcf_recalled
    .toSortedList()
    .set { smoove_sample_recalled_vcfs }

smoove_index_recalled
    .toSortedList()
    .set { smoove_sample_recalled_index }

/*
~ ~ ~ > * Merge re-called smoove
*/

process merge_smoove_cohort {
    
    publishDir "${params.out}/variation", mode: 'copy'

    input:
        file(cohort) from smoove_sample_recalled_vcfs
        file(index) from smoove_sample_recalled_index

    output:
        set file("WI.smoove.square.vcf.gz"), file("WI.smoove.square.vcf.gz.tbi") into smoove_wi_cohort
        set file("WI.LargeRemoved.smoove.square.vcf.gz"), file("WI.LargeRemoved.smoove.square.vcf.gz.tbi") into smoove_largeRem_wi_cohort
        file("WI.smoove.square.stats.txt") into snpeff_smoove_multiqc
        file("WI.LargeRemoved.smoove.square.stats.txt") into snpeff_smoove_large_removed_multiqc

    """
        smoove paste --name WI --outdir . ${cohort}

        tabix -p vcf WI.smoove.square.vcf.gz

        bcftools stats --verbose WI.smoove.square.vcf.gz > WI.smoove.square.stats.txt

        bcftools filter -e 'INFO/SVLEN>100000' WI.smoove.square.vcf.gz | \\
        bcftools filter -e 'INFO/SVLEN<-100000' \\
        -Oz -o WI.LargeRemoved.smoove.square.vcf.gz

        tabix -p vcf WI.LargeRemoved.smoove.square.vcf.gz

        bcftools stats --verbose WI.LargeRemoved.smoove.square.vcf.gz > WI.LargeRemoved.smoove.square.stats.txt
    """

}

/*
~ ~ ~ > * Annotate smoove
*/

process smoove_snpeff {
    
    cpus params.cores

    publishDir "${params.out}/variation", mode: 'copy'

    input:
        set file(smoovevcf), file(smooveindex) from smoove_wi_cohort

    output:
        set file("WI.smoovesv.snpEff.vcf.gz"), file("WI.smoovesv.snpEff.vcf.gz.tbi") into snpeff_smoove_vcf
        file("smoovsv_snpeff_out.csv") optional true into snpeff_smoove_snpeff_multiqc

    script:
      """
        bcftools filter -e 'INFO/END < POS' ${smoovevcf} |\\
        bcftools filter -e 'INFO/AC = 0' |\\
        bcftools view --threads ${task.cpus-1} | \\
        snpEff eff -v chicken_new \\
                 -csvStats snpeff_out.csv \\
                 -no-downstream \\
                 -no-intergenic \\
                 -no-upstream \\
                 -nodownload \\ | \\
        bcftools view -O v | \\
        bcftools view -Oz -o WI.smoovesv.snpEff.vcf.gz

        tabix -p vcf WI.smoovesv.snpEff.vcf.gz
      """
}


/*
==============================================
~ > *                                    * < ~
~ ~ > *                                * < ~ ~
~ ~ ~ > * Generate Multi-caller VCFs * < ~ ~ ~
~ ~ > *                                * < ~ ~
~ > *                                    * < ~
==============================================
*/

manta_to_db
    .join(smoove_to_db)
    .into { variant_db;
            variant_to_survivor;
            variant_to_publish }

/*
~ ~ ~ > * Merge all callers at sample level - SURVIVOR
*/

process survivor_merge {

    tag { SM }

    publishDir "${params.out}/MultiCaller_Isotype_VCFs/SURVIVOR", mode: 'copy'

    input:
        set val(SM), file(mantasv), file(smoovesv) from variant_to_survivor

    output:
        set file("${SM}_survivor_loose.vcf.gz.tbi"),file("${SM}_survivor_loose.vcf.gz") into sm_survivor_merged_vcf
        set val(SM), file("${SM}_survivor_loose_SnpEff.vcf.gz"), file("${SM}_survivor_loose_SnpEff.vcf.gz.tbi") into three_caller_snpeff

    script:
      """
        echo ${SM}_manta > change_manta.txt
        echo ${SM}_delly > change_delly.txt
        echo ${SM}_smoove > change_smoove.txt

        bcftools reheader -s change_manta.txt -o manta_to_merge.vcf ${mantasv}
        bcftools reheader -s change_smoove.txt -o smoove_to_merge.vcf ${smoovesv}

        ls *merge.vcf > sample_files

        SURVIVOR merge sample_files 1000 1 0 0 1 30 ${SM}_survivor_loose.vcf
        bcftools sort -Oz -o ${SM}_survivor_loose.vcf.gz ${SM}_survivor_loose.vcf
        tabix -p vcf ${SM}_survivor_loose.vcf.gz

        bcftools view -O v ${SM}_survivor_loose.vcf.gz | \\
        bcftools filter -e 'POS > INFO/END & INFO/SVTYPE="INV"' |\\
        snpEff eff -v chicken_new \\
                 -csvStats snpeff_out.csv \\
                 -no-downstream \\
                 -no-intergenic \\
                 -no-upstream \\
                 -nodownload \\ | \\
        bcftools view -O v | \\
        bcftools view -Oz -o ${SM}_survivor_loose_SnpEff.vcf.gz

        tabix -p vcf ${SM}_survivor_loose_SnpEff.vcf.gz
      """
}

process survivor2bed {

    tag { SM }

    publishDir "${params.out}/MultiCaller_Isotype_VCFs/SURVIVOR/ThreeCallerBed", mode: 'copy'

    input:
        set val(SM), file(effvcf), file(effvcfindex) from three_caller_snpeff

    output:
        set val(SM), file("${SM}_HQ_SnpEff.bed"), file("${SM}_SnpEff.bed") into three_caller_snpeff_bed

    script:
      """
        bcftools query ${effvcf} -f '[%CHROM\\t%POS\\t%END\\t%INFO/SUPP\\t%INFO/SVTYPE\\t%ST\\t%TY\\t%CO\\t%SAMPLE\\t%GT\\t%INFO/ANN\\n]' |\\
        awk -F"|" '{print \$1, \$2, \$3, \$4, \$5}' OFS="\\t" > ${SM}_SnpEff.bed

        grep "1/1" ${SM}_SnpEff.bed > ${SM}_HQ_SnpEff.bed
      """
}

process plot_SV {

    tag { SM }

    publishDir "${params.out}/MultiCaller_Isotype_VCFs/SURVIVOR/PLOTS", mode: 'copy', pattern: '*.pdf'
    publishDir "${params.out}/MultiCaller_Isotype_VCFs/SURVIVOR/HQ_BED", mode: 'copy', pattern: '*processed_SV.bed'

    input:
        set val(SM), file(HQ_SV), file(allSV) from three_caller_snpeff_bed

    output:
        file("${SM}_SVplot.pdf") into sv_plots
        file("${SM}_processed_SV.bed") into pr_sv_beds

    script:
      """
        Rscript --vanilla `which Plot_SV.R` ${HQ_SV} ${SM}
      """
}


pr_sv_beds
    .toSortedList()
    .set { merged_sample_sv_beds}

process combine_survivor_beds {

    publishDir "${params.out}/MultiCaller_Isotype_VCFs/SURVIVOR/GenotypeMatrix", mode: 'copy'

    input:
        file(beds) from merged_sample_sv_beds

    output:
        set file("Population_SV_GenotypeMatrix.tsv"), file("combinedSV.bed") into sv_survivor_gm

    script:
        """
        cat *.bed > combinedSV.bed

        Rscript --vanilla `which Survivor_SV_to_GM.R` combinedSV.bed
        """

}


