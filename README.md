# 📦 scMOSAIC

## BootLeg Whole Genome Sequencing R package #TODO change this name maybe, something that incorprates ATAC and RNAseq, pacbio, and 10X
*an R package for parsing and combining pacbio long read sequencing with 10X multiomic data*

The function of this package essentially parses 10X multiome ATAC reads to extract the DNA. There are also other basic functions for demultiplexing and counting CAG repeats.

*R*
---

### 📦 R Installation

```R
devtools::install_github("ratoncito/blwgsR", force = TRUE)
library(blwgsR)

#parse the read into inserts barcodes sampleIDs 
read_info <- parse_reads(blwgsR::Fitzwalter2025, out.file = "test.csv")

#count CAG repeats and demultiplex the reads into cells 
cell_info <- demultiplex_reads(reads = read_info, out.file = "None")
```
#construct diagram 
![construct information](./chromium_barcode_info.png)

#  parse reads function

_input_ Path to the input file

_out.file_ the location of the output file, where you want your results to be saved #TODO do we want this or should we just return the table

_plot_ boolean of if you want plots to be created along the way to visualize the process

_P5element_ TODO is this still necessary? #TODO remove

_Spacer_ spacer element from 10X (see construct diagram)

_Read1N_ element from 10X (see construct diagram)

_Read2N_ element from 10X (see construct diagram)

_P7Element_ element from 10X (see construct diagram)

_HTTexon1_p1_ TODO change this name to something more general; desired sequence match start sequence

_HTTexon1_p2_ TODO change this name to something more general; desired sequence match start sequence

_SampleIndexes_ optional used to subset and filter the reads based on if you know what the sample indexes you are expecting

_P5element.max.mm_ the max allowed mismatch for the P5 element TODO still necessary? this should be removed prior to running #TODO remove

_Spacer.max.mm_ the max allowed mismatch for the P5 element TODO still necessary?

_Read1N.max.mm_ the max allowed mismatch for the read1n element

_Read2N.max.mm_ the max allowed mismatch for the read2n element

_P7Element.max.mm_ the max allowed mismatch for the TODO still necessary?

_HTTexon1_p1.max.mm_ the max allowed mismatch for the first part of the desired sequence matched gene #TODO make this not names HTT

_barcode.length_ how long the 10X barcodes are expected to be

_barcode.max.mm_ the max mismatch allowed for barcode from the whitelist to be considered a valid barcode

_sampleid.length_ how long the 10 sample barcode is expected to be

_sampleid.max.mm_ the max mismatch allowed for barcode from the SampleIndexes to be considered a valid sampleID barcode

_file.format_ if the file is in fastq or not, fasta?

