# BreakPointAssembly
A tool to quickly assembly SV breakpoints using Long Reads

The bp_assemble.py script uses samtools, minimap2 and racon to assemble and polish a list of candidate SV breakpoints. Taking a tsv list of breakpoint positions as input along with the read fastq and bam the script follows 5 steps:

1. Extract reads at the breakpoint positions
2. Find reads that support and span the breakpoint on both chromosome copies
3. Generate scaffold breakpoint sequences using the longest reads that support each arm
4. Align all reads at breakpoint positions to the scaffolds
5. Polish the scaffold sequence using racon

## Dependencies
- pysam & samtools
- mappy & minimap2
- racon

## Setup, Usage
```sh
# Setup:
git clone https://github.com/adcosta17/BreakPointAssembly.git
cd BreakPointAssembly

# Usage: 
python bp_assemble.py --sniffles-input <sniffles_translocation_calls.tsv> \
--input-bam <input.bam> \
--input-fastq <input.fastq> \
--output-folder <path/to/output/folder> \
--reference-genome <reference_genome.fa>

```

## Arguments:
**--sniffles-input** A tsv of SV calls. 6 columns are needed: chromsome_A, start, end, chromsome_B, start, end
**--input-bam** A bam file containing alignments of reads to the reference genome
**--input-fastq** A fastq of the reads
**--reference-genome** A reference genome fasta
**--racon** [Optional] The path to racon. By default assumes racon is in the PATH
