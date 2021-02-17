# From Assemblies to genes


Detect coding genes in an assembly, and get their alignments and description using fromAssembly2gene. fromAssembly2gene is a perl script run in command line that uses several available programs and R packages to identify genes of your interests in an assembled genome and outputs a descriptive table, the alignment of your gene against the gene in the assembly and the predicted peptide.
It runs several steps:

1. Predicts genes using prodigal https://github.com/hyattpd/Prodigal
2. Finds matches of the genes of interest in the predicted genes in the assembly using local blast.
3. Refines the alignments using an R scripts using the "msa", "reshape2", "Biostrings", "seqinr" packages
4. Using the alignments prints out a table that describes the findings. It generates descriptive tables of presence, absence and truncated genes, curated alignments and peptides predictions

The alternative run predicting plasmids:

--Kleb, --Esch or --Ent options. These setting works similarly with the exception that before the prediction of genes, it uses "mlplasmids" to identify plasmids (https://gitlab.com/sirarredondo/mlplasmids)
it separates chromosome and plasmids from assembled genomes and find the genes of your interest.	

## Requires
 
The required inputs are assemblies and genes in fasta formats.
The genomes must be located in a folder together with not additional files. 
The genes of interest must be located in a folder together with not additional files.
The genes must have the extension ".fasta" as it is used as a tag for handling file in the program. 
The assemblies can have any kind of extension. 

Please make sure you have the following dependencies in your computer:											
Dependencies: Blast, prodigal and samtools. R packages: "mlplasmids", "msa", "reshape2", "Biostrings", "seqinr"

The program has been parallelised for efficiency.
 
## INSTALL 
```
git clone https://github.com/LPerlaza/Assembly2Gene.git
cd Assembly2Gene
chmod a+x fromAssembly2gene.pl
```	 							
## RUN
																																	
run like:
```
																																	
 From Assemblies to genes (Version February 2021)																									
 																																	
Detect genes in an assembly, and get their detailed alignments and description using fromAssembly2gene. 		
It generates descriptive tables of presence, absence and truncated  genes, curated alignments and peptides predictions.													
Dependencies: Blast, prodigal and samtools. R packages "msa", "reshape2", "Biostrings", "seqinr".
The R packages "mlplasmids" is needed if --Kleb, --Esch or --Ent are used.
																																	
run like: 

./fromAssembly2gene.pl -g gene/*fasta -a genome/* -o GenesInterest -c 10 -np												


OPTIONS
	--assemblies  -a Put all your files (assemblies) in a folder and write here the path with * at the end.							
	--genes 	  -g Put all your files (genes in nucleotide sequences) in a folder and write here the path with *fasta at the end. extension ".fasta".	
	--out 		  -o prefix for output folders																						
	--cores		  -c number of cores to use (4 default).
	--Kleb        -k When this option is set the program predicts plasmids in Klebsiella pneumonia. 
	--Esch        -e When this option is set the program predicts plasmids in Escherichia coli.
    --Ent         -n When this option is set the program predicts plasmids in  Enterococcus faecium.
    --help        -h print this help
    --version     -v version

NOTE 1: When you write * at the end of the path the program will take every file in the folder. 
NOTE 2: The genes must have the extension ".fasta" as it is used as a tag for handling file in the program. 
NOTE 3: The blast indexes are going to be created inside the genes folder. Have this in mind when runing several times, use *fasta to avoid using the index files as inputs.	
NOTE 4: Make sure your input genes have a starting and stop codon, and are not in reverse. Unexpected mistakes on the identification of the gene can happen when these are not met. We could account for this but the classification of possibly truncated matches is based on the movement of the stop codon in comparison with the reference.
NOTE 5: If you failed to note 4, the program will add a ATG at the start of your sequence and a TAG at the end.
NOTE 6: Note that the prediction of plasmids is only possible for the species that mlplasmids is developed for. If your assemblies are not listed, you don't need to set the option and the program will assumed all your assembly is chromosome. 
NOTE 7: The prediction of Plasmids is done by one species at a time, if you have different species, run the program separately 		

```

## OUTPUTS

folder **test_output**:
	This folder contains all the intermedia files used for the program. These files will help you to check in detail where your alignment come from. In case you are puzzle by your final table. Each folder contains:

  **GenomeName.chr.genes.faa**: all predicted genes in AA
  **GenomeName.chr.genes.fasta**: all predicted genes in nt
  **GenomeName.chr.genes.gff**: all predicted genes in gff format with the fasta file at the end
  **GenomeName.chromosome.GeneName.Blast.txt**: blast results of the gene against the genome
  **GenomeName.chromosome.GeneName.fasta.plusRef.fasta**: fasta of gene reference and gene in the genome

folder *test_results*
	This folder contains all final results files and a folder with the predicted peptides that match with the genes of interest
  **GeneName.fasta.nt_alignment.fasta**: The alignment of each gene of interest for all the genomes analysed
  **test.alignments.description.txt**: table with the descriptive information of the alignments, stop codons, gaps, insertions, SNPs, N.copies (numbers of copies)
  **Peptides** (folder): predicted peptides that match with the genes of interest

Additional files when running fromAssembly2gene_Klebsiella.pl

folder **test_output**:
	Each assembly has two folders one for the chromosome and one for the plasmid. Examples here are about the chromosome, *.plasmid_* are for plasmids  
   **GenomeName.chromosome_contigslist.txt**:list of contigs in the assembly that are chromosomal
   **GenomeName.chromosome.fasta**: fasta file of chromosomal contigs
   **GenomeName.fsa_nt.chromosomesummary.txt**: summary results from chromosome prediction from mlplasmid
	
