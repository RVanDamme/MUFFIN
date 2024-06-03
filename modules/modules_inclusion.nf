if (params.modular=="full" | params.modular=="assemble" | params.modular=="assem-class" | params.modular=="assem-annot") {
    include {chopper} from './modules/ont_qc' params(output: params.short_qc)
    include {fastp} from './modules/fastp' params(output: params.output)
    include {spades} from './modules/spades' params(output: params.output)
    include {spades_short} from './modules/spades' params(output: params.output)
    include {flye} from './modules/flye' params(output : params.output)
    include {minimap_polish} from'./modules/minimap2'
    include {racon} from './modules/polish'
    include {medaka} from './modules/polish' params(model : params.model)
    include {pilon} from './modules/polish' params(output : params.output)
    include {pilong} from './modules/polish' params(output : params.output)
    include {minimap2} from './modules/minimap2' //mapping for the binning 
    include {extra_minimap2} from './modules/minimap2'
    include {bwa} from './modules/bwa' //mapping for the binning
    include {extra_bwa} from './modules/bwa'
    //include {metabat2_extra} from './modules/metabat2' params(output : params.output)    
    include {metabat2} from './modules/metabat2' params(output : params.output)
    include {semibin2} from './modules/semibin2' params(output : params.output)
    include {comebin} from './modules/comebin' params(output : params.output)
    include {bam_merger} from './modules/samtools_merger' params(output : params.output)
    //include {contig_list} from './modules/list_ids'
    include {cat_all_bins} from './modules/cat_all_bins'
    include {bwa_bin} from './modules/bwa'  
    include {minimap2_bin} from './modules/minimap2'
    include {metaquast} from './modules/quast' params(output : params.output)
}
if (params.modular=="full" | params.modular=="classify" | params.modular=="assem-class" | params.modular=="class-annot"){
    include {sourmash_download_db} from './modules/sourmashgetdatabase'
    include {checkm_download_db} from './modules/checkmgetdatabases'
}
if (params.modular=="full" | params.modular=="classify" | params.modular=="assem-class" | params.modular=="class-annot") {
    include {reads_retrieval} from './modules/seqtk_retrieve_reads' params(output : params.output)
    include {checkm2} from './modules/checkm2' params(output: params.output)
    include {sourmash_bins} from './modules/sourmash' params(output: params.output)
    include {sourmash_checkm_parser} from './modules/checkm_sourmash_parser' params(output: params.output)
    //include {unmapped_retrieve} from './modules/seqtk_retrieve_reads' params(output : params.output)
}
if (params.modular=="full" | params.modular=="annotate" | params.modular=="assem-annot" | params.modular=="class-annot") {
    include {eggnog_download_db} from './modules/eggnog_get_databases'
    include {eggnog_bin} from './modules/eggnog' params(output: params.output)
    include {fastp_rna} from './modules/fastp' params(output: params.output)
    include {de_novo_transcript_and_quant} from './modules/trinity_and_salmon' params(output: params.output)
    include {eggnog_rna} from './modules/eggnog' params(output: params.output)
    include {parser_bin_RNA} from './modules/parser' params(output: params.output)
    include {parser_bin} from './modules/parser' params(output: params.output)
}