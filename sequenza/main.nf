#!/usr/bin/env nextflow
// from https://github.com/chelauk/nf-core-mutectplatypus/blob/master/workflows/mutect_platypus.nf
nextflow.enable.dsl=2

include { SAMTOOLS_CONVERT  as  SAMTOOLS_CONVERT_T    } from './modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT  as SAMTOOLS_CONVERT_N   } from './modules/nf-core/samtools/convert/main'
include { SEQUENZAUTILS_GCWIGGLE  } from './modules/nf-core/sequenzautils/gcwiggle/main'
include { SEQUENZAUTILS_BAM2SEQZ  } from './modules/nf-core/sequenzautils/bam2seqz/main'
include { SEQUENZAUTILS_MERGESEQZ } from './modules/local/sequenzautils/mergeseqz/main'
include { SEQUENZAUTILS_BINNING   } from './modules/local/sequenzautils/seqzbin/main'
include { SEQUENZAUTILS_RSEQZ     } from './modules/local/sequenzautils/seq_R/main'

include { samplesheetToList } from 'plugin/nf-schema'

workflow SEQUENZA {
    take:
    input_samplesheet
    fasta
    fai

    main:
    input_samplesheet.branch { meta, cram, crai ->
                        tumour: meta.status == 1
                        normal: meta.status == 0
                    }.set{seq_split}

    seq_split.tumour.map{ meta, cram, crai -> 
                        meta = meta + [id:meta.patient + meta.sample]
                        [meta, cram, crai]
                    }.set{tumour_cram}
                    
    SAMTOOLS_CONVERT_T(tumour_cram,fasta,fai).bam.map{ meta, bam, bai -> 
                         [meta.patient, [ id:meta.id , status:meta.status], bam, bai]
                    }.set{tumour_bam}
                    
    seq_split.normal.map{ meta, cram, crai -> 
                        meta = meta + [id:meta.patient + meta.sample]
                        [meta, cram, crai]
                     }.set{normal_cram}
                    
    SAMTOOLS_CONVERT_N(normal_cram,fasta,fai).bam.map{ meta, bam, bai -> 
                         [meta.patient, [ id:meta.id , status:meta.status], bam, bai]
                      }.set{normal_bam}
                    

    tumour_bam
            .combine(normal_bam, by:0)
            .map { patient, meta1, files1, meta2, files2 ->
                    [meta1 + [patient:patient], files1, files2]
            }
            .set{seq_input_matched}

    if (params.create_wiggle) {
        SEQUENZAUTILS_GCWIGGLE(fasta)
        wiggle = SEQUENZAUTILS_GCWIGGLE.out.wig // something wrong here
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
    SEQUENZAUTILS_RSEQZ(SEQUENZAUTILS_BINNING.out.seqz_bin)
}

workflow  {   
    fasta           = params.fasta ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()
    fai             = params.fasta ? Channel.fromPath(params.fai).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()

    input_samplesheet = params.input ? Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json")) : Channel.empty()
    
    SEQUENZA(input_samplesheet, fasta, fai)

}
