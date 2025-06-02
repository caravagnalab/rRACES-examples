#!/usr/bin/env nextflow
// from https://github.com/chelauk/nf-core-mutectplatypus/blob/master/workflows/mutect_platypus.nf
nextflow.enable.dsl=2

include { SEQUENZAUTILS_GCWIGGLE  } from './modules/nf-core/sequenzautils/gcwiggle/main'
include { SEQUENZAUTILS_BAM2SEQZ  } from './modules/nf-core/sequenzautils/bam2seqz/main'
include { SEQUENZAUTILS_MERGESEQZ } from './modules/local/sequenzautils/mergeseqz/main'
include { SEQUENZAUTILS_BINNING   } from './modules/local/sequenzautils/seqzbin/main'
// include { SEQUENZAUTILS_RSEQZ     } from '../modules/nf-core/sequenzautils/seqz_R/main.nf'

include { samplesheetToList } from 'plugin/nf-schema'

workflow SEQUENZA {
    take:
    input_samplesheet
    fasta

    main:
    input_samplesheet
                    .map{ meta, cram, crai ->
                    meta = meta + [id: "${meta.patient}_${meta.sample}"]
                    [meta.patient,meta.sample,meta.status,meta.id,meta.gender,cram,crai] }
                    .branch {
                        tumour: it[2] == 1
                        normal: it[2] == 0
                    }
                    .set{seq_split}

    seq_split.tumour.combine(seq_split.normal, by:0)
                    .set{seq_input_pair}

    seq_input_matched = seq_input_pair
                    .map{ patient, sample1, status1, id1, gender1, files11, files12, sample2, status2, gender2, id2, files21, files22 ->
                    def meta = [patient:patient, sample:sample1, id:id1, gender:gender1]
                    [ meta, files11, files21] }


    if (params.create_wiggle) {
        SEQUENZAUTILS_GCWIGGLE(fasta)
        wiggle = SEQUENZAUTILS_GCWIGGLE.out.wig
    } else {
        wiggle = Channel.fromPath(params.wiggle).collect() 
    }

    SEQUENZAUTILS_BAM2SEQZ(seq_input_matched,
                            fasta,
                            wiggle)

    SEQUENZAUTILS_BAM2SEQZ.out.seqz
                            .groupTuple()
                            .set{merge_seqz_input}

    SEQUENZAUTILS_MERGESEQZ(merge_seqz_input)
    SEQUENZAUTILS_BINNING(SEQUENZAUTILS_MERGESEQZ.out.concat_seqz, params.bin)
}

workflow  {   
    fasta             = params.fasta ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()
    input_samplesheet = params.input ? Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json")) : Channel.empty()
    
    SEQUENZA(input_samplesheet, fasta)

}