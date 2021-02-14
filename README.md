*From Assemblies to genes*																											


Detect genes in an assembly, and get their detailed alignments and description using fromAssembly2gene. fromAssembly2gene is a perl program that used several available programs and R packages to identify genes of your interests in an assembled genome and outputs a descriptive table, the alignment of your gene against the gene in the assembly and the predicted peptide.
It runs several steps:

1. Predicts genes using prodigal https://github.com/hyattpd/Prodigal
2. Finds matches of the genes of interest in the predicted genes in the assembly using local blast.
3. Refines the alignments using an R scripts using the "msa", "reshape2", "Biostrings", "seqinr" packages
4. Using the alignments prints out a table that describes the findings. 

The alternative version for Klebsiella:

fromAssembly2gene_Klebsiella: This program works similarly with the exception that before the prediction of genes uses "mlplasmids" to identify plasmids (https://gitlab.com/sirarredondo/mlplasmids)
fromAssembly2gene_Klebsiella separates chromosome and plasmids from assembled genomes and automatically find genes of interest. It generates descriptive tables of presence, absence and truncated genes, curated alignments and peptides predictions.		

 *Requires*
 
The required inputs are assemblies and genes in fasta formats.
The genomes must be located in a folder together with not additional files. 
The genes of interest must be located in a folder together with not additional files.
The genes must have the extension ".fasta" as it is used as a tag for handling file in the program. 
The assemblies can have any kind of extension. 

Please make sure you have the following dependencies in your computer:											
Dependencies: Blast, prodigal and samtools. R packages: "mlplasmids", "msa", "reshape2", "Biostrings", "seqinr"

The program has been parallelised for efficiency.
 
	 							
*Run*
																																	
run like: 
	./fromAssembly2gene.pl -g gene/*fasta -a genome/* -o test	
	./fromAssembly2gene_Klebsiella.pl -g gene/*fasta -a genome/* -o test													

OPTIONS:																																	
	--assemblies  -a Put all your files (assemblies) in a folder and write here the path with * at the end.							
	--genes 	  -g Put all your files (genes) in a folder and write here the path with * at the end. extension ".fasta" required		
	--out 		  -o prefix for output folders																						

OUTPUTS

folder test_output:
This folder contains all the intermedia files used for the program. These files will help you to check in detail where your alignment come from. In case you are puzzle by your final table. Each folder contains:

GenomeName.chr.genes.faa:all predicted genes in AA
GenomeName.chr.genes.fasta: all predicted genes in nt
GenomeName.chr.genes.gff: all predicted genes in gff format with the fasta file at the end
GenomeName.chromosome.GeneName.Blast.txt: blast results of the gene against the genome
GenomeName.chromosome.GeneName.fasta.plusRef.fasta: fasta of gene reference and gene in the genome

folder test_results
This folder contains all final results files and a folder with the predicted peptides that match with the genes of interest
				GeneName.fasta.nt_alignment.fasta: The alignment of each gene of interest for all the genomes analysed
				out.alignments.description.txt: table with the descriptive information of the alignments, stop codons, gaps, insertions, location (chr,plasmid), SNPs, N.copies (numbers of copies)
				Peptides (folder): predicted peptides that match with the genes of interest

Additional files when running fromAssembly2gene_Klebsiella.pl

folder test_output:
Each assembly has two folders one for the chromosome and one for the plasmid. 
GenomeName.chromosome_contigslist.txt:list of contigs in the assembly that are chromosomal
GenomeName.chromosome.fasta: fasta file of chromosomal contigs
GenomeName.fsa_nt.chromosomesummary.txt: summary results from chromosome prediction from mlplasmid
	
