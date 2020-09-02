'''
Todo:
1 get primer3 settings for NGS primer design
2 run primer3 (3 pairs per position)
3 run in-silico PCR on command-line
4 Check if primers are good.
'''



# - environment requirement
# conda install -c bioconda primer3-py
# conda install -c anaconda cython
# -------------------------------import---------------------------------------

import urllib.request
import primer3
import math
import time
import random
import socket
import re
import sys, csv, operator
import ssl
# ---------------------------------import---------------------------------------


# ---------------------------------functions------------------------------------
def read_by_line(url): # read by line
    lines = url.splitlines()
    return lines


def my_sq_range(start, end, step):
    while start < end:
        yield start
        start += step


# find the flanking sequences for primer design
def flanking_seq(ref, chr, ot_start, amp_max):
    try:
        flank_start = max(1, int(ot_start) - amp_max)  # set to 1 if it is too close to the chr starting position
        flank_end = int(ot_start) + 30 + amp_max  # Assuming the length of gRNA is at most 30
        dna_length = flank_end - flank_start + 1
    except ValueError:
        print("Wrong position information")

    else:
        url = "http://genome.ucsc.edu/cgi-bin/das/"+ref+"/dna?segment=" \
              + chr + "%3A" + str(flank_start) + "," + str(flank_end)

        # print(url)
        # req = urllib.request(url)
        req = urllib.request.Request(url)
        try:
            with urllib.request.urlopen(req) as response:
                resp_html = response.read()
                resp_html_str = str(resp_html)
                flank_seq = resp_html_str[resp_html_str.find('</DNA>') - dna_length-28:resp_html_str.find('</DNA>')]
                seq_array = flank_seq.split("\\n")
        except urllib.error.URLError as e:
            print(e.reason)
    str_empty = ''
    return str_empty.join(seq_array)


# Generate required setting file for primer3
def primer_design_settings(ref_seq, amplicon_min, amplicon_max):
    print("rua")
    # os.system(qry + ' >> ' + output_dir)

# ---------------------------------functions------------------------------------


# testing main
# original csv: chr19,41352729,CCATTAGCACGCGGGTGACCTCCTTGGCGTAGTAGTCGGCCTCAGGCTCG
# output = flanking_seq("hg38", "chr19", 41352729, 320)
# print(output)
output = 'cctctaggggactgcccccacgaccccgcatgtttctgtcgcactctagaagcggtccacttcgctatctcctcctctcca' \
         'agaccagacacctgggtggtagggggctcagtgccatcctctttcggacacccccctcccaccatcacacgttccctttgc' \
         'cccggggtgtcctcttcctccagccagtttcttctgccagtcacttcctacccgtggccccggcactccggcgccccctgg' \
         'gggcccccctcccggctcccctgcccctccgagctcaccgttgtgggtttccaccattagcacgcgggtgacctccttggc' \
         'gtagtagtcggcctcaggctcgggctccggttctgcactctccccggccacccggtcgcgggtgctgttgtacagggcgag' \
         'cacggcctcgggcagcgggccgggcggcacctccccctggctcggggggctggcgagccgcagcttggacaggatctggcc' \
         'gcggatggcctcgatgcgcttccgcttcaccagctccatgtcgatagtcttgcaggtggatagtcccgcggccggccggcc' \
         'aggcgtcagcaccagtagccacagcagcggtagcagcagcggcagcagccgcagcccggagggcggcatgggggaggcggc' \
         'gccccccggcactgccgagagcg'

test = primer3.bindings.designPrimers(
    {
        'SEQUENCE_ID': 'MH1000',
        'SEQUENCE_TEMPLATE': 'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTT' \
                              'AGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCA' \
                              'ACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACG' \
                              'CACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAG' \
                              'TTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAA' \
                              'TGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAA' \
                              'ATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAA' \
                              'TTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCA' \
                              'GATGTTTCCCTCTAGTAG',
        'SEQUENCE_INCLUDED_REGION': [36,342]
    },
    {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_PRODUCT_SIZE_RANGE': [[75,100],[100,125],[125,150],
                                      [150,175],[175,200],[200,225]],
    })
print(test['PRIMER_RIGHT_0_SEQUENCE'], test['PRIMER_LEFT_0_SEQUENCE'])


# g = ['catcaggagagcgcttgctcacccaccccctctctattgatttgttcact', 'gtgttaatggctcccttttgtgtttggctggtcagtgttttgaattggtc', 'tttaattatatctcatgcaaatgactaaggcatatggagaggttgtggag', 'ccaaaagtttagaagcacctattatgtgccactcctactactcttctctg', 'caacataaatagaggaaactgaggctcggaggagttaagtaaaacactgg', 'tcaagttggcccaggggcccacatctcctgaccccagcataggctcttta', 'tgagtttggatggatgagatctttggatgagatatctcagaacaaggagc', 'cagcctgggcaatccgtttgccaagctagagtcactgagcaagccagaca', 'cctacagatgccacaggtaggaacaactgctgtgtactctagctcagctt', 'gagaatagttggtttaatttatgaatttagatggttgaaaaaccaaaagc', 'catagaaacaatcagcacttgtcaagtgcctactatgtgccaggaacttc', 'tgtttcccccttttcatctgtgaaatggggatgacagtaaaataatgcca', 'acctctcctgattgatgaaggatcaagtgagagccatggaatcacttagc', 'acagaggctgtacttaatatg', '']

# str = ''
# print(str.join(g))


'''

f = 'catcaggagagcgcttgctcacccaccccctctctattgatttgttcact\ngtgttaatggctcccttttgtgtttggct' \
    'ggtcagtgttttgaattggtc\ntttaattatatctcatgcaaatgactaaggcatatggagaggttgtggag\nccaaaagtt' \
    'tagaagcacctattatgtgccactcctactactcttctctg\ncaacataaatagaggaaactgaggctcggaggagttaagtaaaacactg' \
    'g\ntcaagttggcccaggggcccacatctcctgaccccagcataggctcttta\ntgagtttggatggatgagatctttggatgagatatctcagaac' \
    'aaggagc\ncagcctgggcaatccgtttgccaagctagagtcactgagcaagccagaca\ncctacagatgccacaggtaggaacaactg' \
    'ctgtgtactctagctcagctt\ngagaatagttggtttaatttatgaatttagatggttgaaaaaccaaaagc\ncatagaaacaatcagca' \
    'cttgtcaagtgcctactatgtgccaggaacttc\ntgtttcccccttttcatctgtgaaatggggatgacagtaaaataatgcca\nacctctcctgat' \
    'tgatgaaggatcaagtgagagccatggaatcacttagc\nacagaggctgtacttaatatg\n'

strings = re.findall(r"\S+", f)
str1 = ''

print(str1.join(strings))


# get parsed information
parser = argparse.ArgumentParser(description='Automatic primer design')
parser.add_argument('--ref', default='hg38', type=str,
                    help='Reference genome. Default: hg38.')
parser.add_argument('--min', default=280, type=int,
                    help='Minimum amplicon size.')
parser.add_argument('--max', default=320, type=int,
                    help='Maximum amplicon size.')
args = parser.parse_args()

'''