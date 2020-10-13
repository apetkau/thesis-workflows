#!/usr/bin/env nextflow

nextflow.preview.dsl = 2
 
params.input = "$baseDir/data/*_{1,2}.fastq.gz"
params.threads = 4
params.kmer = 31
params.hashes = 1000

process mash_index {
    input:
        tuple name, file(pair)
 
    output:
        path '*.msh', emit: sketches
 
    """
    mash sketch -r -p ${params.threads} -s ${params.hashes} -k ${params.kmer} -o ${name}.msh -I ${name} ${pair[0]} ${pair[1]}
    """
}

process sketch_paste {
    input:
        path sketch_list

    output:
        path 'sketches.msh'

    """
    mash paste sketches.msh $sketch_list
    """ 
}

workflow {
    infiles = channel.fromFilePairs(params.input)

    mash_index(infiles)
    sketch_paste(mash_index.out.sketches.collect())
}
