library(scMOSAIC)

#parse the read into inserts barcodes sampleIDs 
read_info <- parse_reads(scMOSAIC::Fitzwalter2025, out.file = "none")

#count CAG repeats and demultiplex the reads into cells 
cell_info <- demultiplex_reads(reads = read_info, out.file = "None")

