#!/usr/bin/env nextflow
 
params.input = "$baseDir/data/*_{1,2}.fastq.gz"
params.threads = 4
params.kmer = 31
params.hashes = 1000

infiles = channel.fromFilePairs(params.input)

/*
 * Index a fastq file with Mash
 */
process mash_index {

    input:
    tuple name, file(pair) from infiles
 
    output:
    path '*.msh' into sketches
 
    """
    mash sketch -r -p ${params.threads} -s ${params.hashes} -k ${params.kmer} -o ${name}.msh -I ${name} ${pair[0]} ${pair[1]}
    """
}
 
/*
 * print the channel content
 */
sketches.subscribe { println it }
