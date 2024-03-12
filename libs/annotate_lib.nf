lib {
    include {eggnog_download_db} from '../modules/eggnog_get_databases'
    include {eggnog_bin} from '../modules/eggnog' params(output: params.output)
    include {fastp_rna} from '../modules/fastp' params(output: params.output)
    include {de_novo_transcript_and_quant} from '../modules/trinity_and_salmon' params(output: params.output)
    include {eggnog_rna} from '../modules/eggnog' params(output: params.output)
    include {parser_bin_RNA} from '../modules/parser' params(output: params.output)
    include {parser_bin} from '../modules/parser' params(output: params.output)
}