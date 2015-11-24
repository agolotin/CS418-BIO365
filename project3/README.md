#Project 3 Mapping Project
The working and the most recent version resides in *dna_mapper_optimal.py_* Works for both .fasta and .fastq files.<br>
The input files can be found in the following locations: <br>
1. Example reference sequence and accompanying reads are locatied in *example_input/*<br>
2. RNASEQ reference sequence and accompanying reads are locatied in *rnaseq_input/*. RNASEQ represents a small example dataset with no errors.<br>
3. Shotgun reference sequence and accompanying reads are located in *shotgun_small_input/*. This is a larger dataset with 100,000 reads and errors.<br>
* Since the working version of the program takes a long time to process all of the 100,000 reads the reads are split into chunks and are located in *shotput_small_input/chunks/*. It is possible to run the program with a chunk number that is being processed.<br>
<br>
Directories with already generated output:<br>
1. *example.fasta.output* contains the output from the example dataset<br>
2. *rnaseq.chr.fasta* contains the output from the RNASEQ dataset<br>
3. *shotgun.small.sam* contains the output from the Shotgun dataset<br>
* There are 3 output files in the *shotgun.small.sam* because of the 3 different kmer lengths tested when checking for errors in reads<br>
<br>
Other directories:<br>
1. *saved_objects/* stores SA, BTW, 1st, LTF in a pickeled form for a specific dataset.<br>
2. *output_sam/* is the default output directory.<br>
3. *versions/* has previous versions of the dna mapper<br>
4. *slurm/* contains the output files from running the large dataset as jobs on the supercomputer<br>
