#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
 
params.input = "$baseDir/data/*_{1,2}.fastq.gz"
params.threads = 4
params.kmer = 31
params.hashes = 1000

process mash_index {
    input:
        tuple val(name), file(pair)
 
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
        path 'sketches.msh', emit: sketches_file

    """
    mash paste sketches.msh $sketch_list
    """ 
}

process construct_dist_matrix {
    input:
        path sketches_file

    output:
        path 'matrix.txt', emit: matrix

    """
    mash triangle ${sketches_file} > matrix.txt
    """
}

process construct_tree {
    input:
        path matrix

    output:
        path 'tree.png'
        path 'tree.txt'

    """
    #!/usr/bin/env python

    import numpy as np
    from skbio import DistanceMatrix
    from skbio.tree import nj
    from ete3 import Tree
    
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
            
            dist_matrix = DistanceMatrix(matrix, names)
            
            return dist_matrix
        
    matrix = parse_triangle('${matrix}')

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

    mash_index(infiles)
    sketch_paste(mash_index.out.sketches.collect())
    construct_dist_matrix(sketch_paste.out.sketches_file)
    construct_tree(construct_dist_matrix.out.matrix)
}
