lib {
    include {reads_retrieval} from '../modules/seqtk_retrieve_reads' params(output : params.output)
    include {checkm2} from '../modules/checkm2' params(output: params.output)
    include {sourmash_bins} from '../modules/sourmash' params(output: params.output)
    include {sourmash_checkm_parser} from '../modules/checkm_sourmash_parser' params(output: params.output)
}