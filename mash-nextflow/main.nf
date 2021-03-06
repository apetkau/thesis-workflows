#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
 
params.input = "$baseDir/data/*_{1,2}.fastq.gz"
params.doassembly = false
params.threads = 4
params.kmer = 31
params.mindepth = 5
params.hashes = 1000

process assemble_reads {
    publishDir 'assemblies'

    input:
        tuple val(name), file(pair)
 
    output:
        path '*.fasta', emit: assemblies

    script:
        """
        skesa --cores ${params.threads} --reads '${pair[0]},${pair[1]}' --contigs_out '${name}.fasta'
        """
}

process mash_index_contigs {
    input:
        path contigs
 
    output:
        path '*.msh', emit: sketches
 
    script:
        file_tokens = contigs.name =~ /^(.*)\./
        name = file_tokens[0][1]

        """
        mash sketch -p ${params.threads} -s ${params.hashes} -k ${params.kmer} -o ${name}.msh -I ${name} ${contigs}
        """
}

process mash_index_reads {
    input:
        tuple val(name), file(pair)
 
    output:
        path '*.msh', emit: sketches
 
    script:
        """
        mash sketch -r -p ${params.threads} -s ${params.hashes} -k ${params.kmer} -m ${params.mindepth} -o ${name}.msh -I ${name} ${pair[0]} ${pair[1]}
        """
}

process sketch_paste {
    publishDir 'sketches'

    input:
        path sketch_list

    output:
        path 'sketches.msh', emit: sketches_file

    script:
        """
        mash paste sketches.msh $sketch_list
        """ 
}

process construct_dist_matrix {
    input:
        path sketches_file

    output:
        path 'matrix.txt', emit: matrix

    script:
        """
        mash triangle ${sketches_file} > matrix.txt
        """
}

process reformat_dist_matrix {
    publishDir 'matrix'

    input:
        path matrix_file

    output:
        path 'matrix.tsv', emit: dist_matrix

    script:
        """
        #!/usr/bin/env python
    
        import numpy as np
        import pandas as pd
        
        def parse_triangle(file):
            with open(file, 'r') as file:
                num_samples = int(file.readline().strip())
        
                rows = []
                names = []
        
                for line in file:
                    line = line.strip()
                    tokens = line.split('\t')
        
                    # Remove sample name
                    name = tokens.pop(0)
        
                    # Pad rest of row with '0'
                    tokens.extend([0] * (num_samples - len(tokens)))
        
                    rows.append(tokens)
                    names.append(name)
        
                matrix = np.matrix(rows).astype('float64')
                
                # Restore to full matrix instead of triangular matrix
                matrix = matrix + matrix.transpose()
                
                return pd.DataFrame(matrix, index=names, columns=names)
            
        matrix = parse_triangle('${matrix_file}')
        matrix.to_csv('matrix.tsv', sep='\t')
        """
}

process construct_tree {
    publishDir 'tree'

    input:
        path matrix_file

    output:
        path 'tree.png', emit: image
        path 'tree.txt', emit: tree

    script:
        """
        #!/usr/bin/env python
    
        import numpy as np
        import pandas as pd
        from skbio import DistanceMatrix
        from skbio.tree import nj
        from ete3 import Tree
        
        matrix = pd.read_csv('${matrix_file}', sep='\t', index_col=0, header=0)
        matrix = DistanceMatrix(matrix.to_numpy(), matrix.index)
    
        # Write tree to file
        nj_tree = nj(matrix, result_constructor=str)
        with open('tree.txt', 'w') as f:
            f.write(nj_tree)
    
        # Draw an image of the tree
        tree = Tree(nj_tree)
        tree.render('tree.png')
        """
}

workflow {
    infiles = channel.fromFilePairs(params.input)

    if (params.doassembly) {
        assemble_reads(infiles)
        mash_index_contigs(assemble_reads.out.assemblies)
        sketch_paste(mash_index_contigs.out.sketches.collect())
    } else {
        mash_index_reads(infiles)
        sketch_paste(mash_index_reads.out.sketches.collect())
    }

    construct_dist_matrix(sketch_paste.out.sketches_file)
    reformat_dist_matrix(construct_dist_matrix.out.matrix)
    construct_tree(reformat_dist_matrix.out.dist_matrix)
}
