#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { blastn } from './modules/blast.nf'
include { filter_best_bitscore } from './modules/blast.nf'

workflow {

  if (params.samplesheet_input != 'NO_FILE') {
    ch_fasta = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['FILE']] }
  } else {
    ch_fasta = Channel.fromPath(params.assembly_search_path)
  }

  ch_db = Channel.fromPath(params.db_dir).combine(Channel.of(params.db_name))

  ch_seqs = ch_fasta.splitFasta(file: true)

  main:
    blastn(ch_seqs.combine(ch_db))
    filter_best_bitscore(blastn.out)
  
}
