# NGS primer design for screened and ranked Off-target sites

## Workflow

**OT chromosome index to Fasta sequence for primer design**

Input from PICOTA: gRNA sequence and chromosome index (in CSV)

User input: Reference genome(default: hg38), optimal amplicon size

Step 1: Chromosome index to Fasta sequence (UCSC) 

Step 2: Label the region that the primers flanking to (based on Primer3 instructions)

**Primer3**  
Local version of gRNA design tool (batch)primer3(plus)

**In-Silico PCR**  
UCSC In-Silico PCR

https://genome.ucsc.edu/cgi-bin/hgPcr

Using BLAST instead will avoid the potential issue of crawling the UCSC server. 
I start with pulling information from the UCSC browser for the initial small-scale usage.
