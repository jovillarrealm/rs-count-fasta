# count-fasta-rs

This utility gives some stats on assembly reports

    Usage: rs-count-fasta [-c csv_file] <fasta_file>

## Usage 
When a csv file is not specified 

    rs-count-fasta ../cnsg-scripts/GENOMIC/GCA_024699835_Aphelenchoides-besseyi_AORJ.fna 

Output to stdout will look like this

    Total length of sequence:       46759715 bp
    Total number of sequences:      32
    Average contig length is:       1461241 bp
    Largest contig:         18100598 bp
    Shortest contig:                214 bp
    N25 stats:                      25% of total sequence length is contained in the 1 sequences >= 18100598 bp
    N50 stats:                      50% of total sequence length is contained in the 2 sequences >= 16068654 bp
    N75 stats:                      75% of total sequence length is contained in the 3 sequences >= 10965501 bp
    Total GC count:                 19534458 bp
    GC %:                           41.78 %
    Number of Ns:                   2900
    Ns %:                           0.01 %

If the csv file was specified, then the created file will look like this

filename;assembly_length;number_of_sequences;average_length;largest_contig;shortest_contig;N50;GC_percentage;total_N;N_percentage
"GCA_024699835_Aphelenchoides-besseyi_AORJ.fna";46759715;32;1461241.09;18100598;214;16068654;41.78;2900;0.01

### Notes on usage

Concurrency is not recommended when a new file has to be created, otherwise the header line may result out of order.

As a mitigation you could first run it on single file, and then introduce the concurrency on things like xargs or gnu parallel.

If it runs on an empty file, it will not introduce the header line.