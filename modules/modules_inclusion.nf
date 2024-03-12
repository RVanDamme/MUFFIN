if (params.modular=="full" | params.modular=="assemble" | params.modular=="assem-class" | params.modular=="assem-annot") {
    include {chopper} from './ont_qc' params(output: params.short_qc)
    include {fastp} from './fastp' params(output: params.output)
    include {spades} from './spades' params(output: params.output)
    include {spades_short} from './spades' params(output: params.output)
    include {flye} from './flye' params(output : params.output)
    include {minimap_polish} from'./minimap2'
    include {racon} from './polish'
    include {medaka} from './polish' params(model : params.model)
    include {pilon} from './polish' params(output : params.output)
    include {pilong} from './polish' params(output : params.output)
    include {minimap2} from './minimap2' //mapping for the binning 
    include {extra_minimap2} from './minimap2'
    include {bwa} from './bwa' //mapping for the binning
    include {extra_bwa} from './bwa'
    //include {metabat2_extra} from './metabat2' params(output : params.output)    
    include {metabat2} from './metabat2' params(output : params.output)
    include {semibin2} from './semibin2' params(output : params.output)
    include {comebin} from './comebin' params(output : params.output)
    include {bam_merger} from './samtools_merger' params(output : params.output)
    //include {contig_list} from './list_ids'
    include {cat_all_bins} from './cat_all_bins'
    include {bwa_bin} from './bwa'  
    include {minimap2_bin} from './minimap2'
    include {metaquast} from './quast' params(output : params.output)
}
if (params.modular=="full" | params.modular=="classify" | params.modular=="assem-class" | params.modular=="class-annot"){
    include {sourmash_download_db} from './sourmashgetdatabase'
    include {checkm_download_db} from './checkmgetdatabases'
}
if (params.modular=="full" | params.modular=="classify" | params.modular=="assem-class" | params.modular=="class-annot") {
    include {reads_retrieval} from './seqtk_retrieve_reads' params(output : params.output)
    include {checkm2} from './checkm2' params(output: params.output)
    include {sourmash_bins} from './sourmash' params(output: params.output)
    include {sourmash_checkm_parser} from './checkm_sourmash_parser' params(output: params.output)
    //include {unmapped_retrieve} from './seqtk_retrieve_reads' params(output : params.output)
}
if (params.modular=="full" | params.modular=="annotate" | params.modular=="assem-annot" | params.modular=="class-annot") {
    include {eggnog_download_db} from './eggnog_get_databases'
    include {eggnog_bin} from './eggnog' params(output: params.output)
    include {fastp_rna} from './fastp' params(output: params.output)
    include {de_novo_transcript_and_quant} from './trinity_and_salmon' params(output: params.output)
    include {eggnog_rna} from './eggnog' params(output: params.output)
    include {parser_bin_RNA} from './parser' params(output: params.output)
    include {parser_bin} from './parser' params(output: params.output)
}