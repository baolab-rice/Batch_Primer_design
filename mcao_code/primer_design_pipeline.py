#!/usr/bin/env python3
""" OT chromosome index to FASTA sequence for primer design"""

# Pre-install/requirements:
    # Cython for Primer3 operating: pip install Cython
    # Primer3 local version: conda install -c bioconda primer3-py/ pip install primer3-py
    # wget http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/twoBitToFa | chmod 755 twoBitToFa

"""Write a checking function here."""
from subprocess import Popen, PIPE
import argparse
import time 
import primer3

Primer3_global_args = {
    'PRIMER_OPT_SIZE' : 20, 
    'PRIMER_PICK_INTERNAL_OLIGO' : 1,   # As example setting
    'PRIMER_MIN_SIZE' : 18,
    'PRIMER_MAX_SIZE' : 27,
    'PRIMER_OPT_TM' : 60.0,
    'PRIMER_MIN_TM' : 59.0,
    'PRIMER_MAX_TM' : 63.0,
    'PRIMER_PAIR_MAX_DIFF_TM' : 100.0,
    'PRIMER_MIN_GC' : 20.0, 
    'PRIMER_MAX_GC' : 80.0,
    'PRIMER_MAX_POLY_X' : 5,
    'PRIMER_INTERNAL_MAX_POLY_X' : 5,
    'PRIMER_SALT_MONOVALENT' : 50.0,
    'PRIMER_DNA_CONC' : 50.0,
    'PRIMER_MAX_NS_ACCEPTED' : 0, # As example setting
    'PRIMER_MAX_SELF_ANY' : 8.00,
    'PRIMER_MAX_SELF_END' : 3.00,
    'PRIMER_MAX_END_STABILITY' : 9.0,
    'PRIMER_PRODUCT_SIZE_RANGE' : [[150,250],[100,300],[301,400],[401,500],\
        [501,600],[601,700],[701,850],[851,1000]],
    'PRIMER_TM_FORMULA' : 1,    # SantaLucia 1988
    'PRIMER_SALT_CORRECTIONS' : 1,   # SantaLucia 1988
    'PRIMER_NUM_RETURN' : 5,
    'PRIMER_PAIR_MAX_COMPL_ANY' : 12,
    'PRIMER_PAIR_MAX_COMPL_END' : 8,
    # Misprimings, not sure what is MAX Repeat Mispriming
    'PRIMER_MAX_TEMPLATE_MISPRIMING' : 12.00,
    'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING' : 24.00,
    # Primer product size
    'PRIMER_PRODUCT_OPT_SIZE' : 300,
    'PRIMER_PRODUCT_SIZE_RANGE' : [250,350],  # (1)Check format later (2)There's no checkbox(boolean) for Product Size Input and ignore Product Size Range
    # Sequencing
    'PRIMER_SEQUENCING_LEAD' : 50,
    'PRIMER_SEQUENCING_ACCURACY' : 20,
    'PRIMER_SEQUENCING_SPACING' : 500,
    'PRIMER_SEQUENCING_INTERVAL' : 250
    # didn't find fix the __ prime end of the primer
    # didn't find left/right/internal oligo acronym
    # didn't find some other check boxes
}

def readfile(filename):
    input_info = []
    with open(filename, 'r') as f:
        next(f)
        for line in f:
            newline = line.split(',')
            input_info.append(newline)
    f.close()
    return input_info

def retrieving_seq_ucsc(index,ref,amp_max): # (1)How to deal with the minimum amplicon;(2) Need to check with DAS method
    chromosome = index[0]
    ori_start_pos = index[1]
    start_pos = max(1, int(ori_start_pos) - amp_max)
    end_pos = int(ori_start_pos) + 30 + amp_max

    process = Popen(args = ['twoBitToFa',
                            'http://hgdownload.cse.ucsc.edu/gbdb/{}/{}.2bit'.format(ref, ref),
                            'stdout',
                            '-seq={}'.format(chromosome),
                            '-start={}'.format(start_pos),
                            '-end={}'.format(end_pos),],
                    stdout=PIPE)
    stdout, stderr = process.communicate()
    del(stderr)
    show = ''.join(str(stdout).split('\\n')[1:-1])
    return show

def primer_design(ID,Templete):
    Primer3_seq_args = {
        'SEQUENCE_ID':ID,
        'SEQUENCE_TEMPLATE':Templete,
    }
    primers_raw = primer3.bindings.designPrimers(Primer3_seq_args,Primer3_global_args)
    primers = [
                [primers_raw['PRIMER_LEFT_0_SEQUENCE'],primers_raw['PRIMER_RIGHT_0_SEQUENCE']],
                [primers_raw['PRIMER_LEFT_1_SEQUENCE'],primers_raw['PRIMER_RIGHT_1_SEQUENCE']],
                [primers_raw['PRIMER_LEFT_2_SEQUENCE'],primers_raw['PRIMER_RIGHT_2_SEQUENCE']],
            ]
    return primers

def main():
    # Read CSV file    
    index = readfile(args.input)

    __list = []
    for i in range(len(index)):
        # Retrieving genomic sequence data
        genome_seq = retrieving_seq_ucsc(index[i],args.ref,args.max)
        # Primer design
        primers = primer_design("ID:testing",genome_seq)
        
        __list.append([genome_seq,primers])
    
    for element in __list:
        print(element)
    
if __name__ == "__main__":
    start_time = time.time()

    desc = """
    Welcome to PICOTA, this is a beta version.
    """
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--ref', default = "hg38", type=str, help="input your reference genome, default is hg38")
    parser.add_argument('-i' , '--input', required=True, help="input gRNA sequence and chromosome index, CSV format")
    parser.add_argument('-o', '--output', help="output analysis result")
    parser.add_argument('--min', default=280, type=int, help="Minimum optimal amplicon size")
    parser.add_argument('--max', default=320, type=int, help="Maximum optimal amplicon size")
    parser.add_argument('--thread', default=1, help="Run the pipeline with multiple threads")
    args = parser.parse_args()

    main()
    print("This program took %s seconds." % round(time.time() - start_time,2))
    
    
