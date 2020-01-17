#!/usr/bin/env python3

import os
import sys
import glob

def main ():
    if len( sys.argv ) < 4:
        sys.stderr.write( "Make jellyfish script for all the fasta files in a given directory for parallel processing.\n" )
        sys.stderr.write( "                   -- developmental script by Tao Liu\n\n")
        sys.stderr.write( "need at least 3 parameters: <foldername> <number of scripts> <n mer>\n\n" )
        exit(1)

    dir_name = sys.argv[1]
    if not os.path.isdir( dir_name ):
        sys.stderr.write(f"ERROR: directory {dir_name} can't be found!\n")
        exit(1)

    n_scripts = int( sys.argv[2] )

    n_mer = int( sys.argv[3] )

    # get all fasta file names
    fasta_file_names = glob.glob( f"{dir_name}/*.fasta" )

    n_fasta_files = len(fasta_file_names)

    n_fasta_in_each_script = n_fasta_files // n_scripts

    n_remained = n_fasta_files % n_scripts

    for i in range( n_scripts ):
        with open( f"jellyfish_script{i}.sh", "w" ) as fhd:
            fhd.write("#!/bin/bash\n\n")
            for fn_idx in range( i * n_fasta_in_each_script, (i + 1) * n_fasta_in_each_script ):
                fasta_fn = fasta_file_names[ fn_idx ]
                jf_fn = fasta_fn[:fasta_fn.rfind(".fasta")] + ".jf"
                count_fn = fasta_fn[:fasta_fn.rfind(".fasta")] + ".kmer.txt"
                fhd.write(f"jellyfish count -m {n_mer} -s 100M -C -o {jf_fn} {fasta_fn}; jellyfish dump -c -t -L 10 -o {count_fn} {jf_fn}; rm -f {jf_fn};\n")

    # lastly, remained files go to the last script
    with open( f"jellyfish_script{n_scripts-1}.sh", "a" ) as fhd:
        for fn_idx in range( n_scripts * n_fasta_in_each_script, n_fasta_files ):
            fasta_fn = fasta_file_names[ fn_idx ]
            jf_fn = fasta_fn[:fasta_fn.rfind(".fasta")] + ".jf"
            count_fn = fasta_fn[:fasta_fn.rfind(".fasta")] + ".kmer.txt"
            fhd.write(f"jellyfish count -m {n_mer} -s 100M -C -o {jf_fn} {fasta_fn}; jellyfish dump -c -t -L 10 -o {count_fn} {jf_fn}; rm -f {jf_fn};\n")        

    with open( "run_jellyfish.sh", "w" ) as fhd:
        fhd.write( "#!/bin/bash\n\n" )
        for i in range( n_scripts ):
            fhd.write( f"bash jellyfish_script{i}.sh &\n" )
        fhd.write( "wait;\n" )
            
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted me! ;-) Bye!\n")
        sys.exit(0)
