# kzfp_annotation
Automatically annotates KRAB zinc-finger protein genes in any input genome. A flexible use allows for the automatic annotation of any zinc-finger protein gene, such as SCAN zinc-finger protein genes.

Performs tailored aligments of protein domains of KZFPs, in particular for the so-called "zinc fingerprints". 

This code was last updated in 2018, and should thus be strictly run with the dependencies below. The best would be to resort to a container to "downgrade" the libraries.

## Publications

Iouranova et al., KRAB zinc finger protein ZNF676 controls the transcriptional influence of LTR12-related endogenous retrovirus sequences, Mobile DNA 2022.

Matsushima et al., in prep.

## Instructions

For the theoretical background, refer to sections 4.1 and 4.1 of KZFPs in ancient humans.pdf in the root folder

Use Python 3.4 (and NOT python 2 as it causes problems with the gestion of integer divisions)

Modules:
- Biopython 1.68
- pybedtools 0.7.9
- pandas 0.19.2
- pyvcf 0.6.8
- numpy 1.11.3


commands:
- HMMER 2 (2.3.2, and NOT HMMER 3)

Procedure:

- Go to UCSC and download 1) the full .fa genome (in a single file) and 2) the .chrom.sizes file that lists chromosomes for your species of interest.

- Note that you can select a subset of chromosomes (or reads) on which you want to run the analysis by editing the .chrom.sizes text file. Do NOT leave additional backspaces at the end of the file.

- Typically, the .fa genome file will be compressed in .2bit, so download and use the command twoBitToFa to uncompress it into a .fa file.

- Open mainMultiprocessing.py in your favorite text edition program (microsoft word <3) 

- Edit line 10: genomePath = "a", a = full path to the uncompressed .fa genome file.

- Edit line 11:  chromListPath = "b", b = full path to the text file listing chromosomes

- Edit line 27: pool = multiprocessing.Pool(). You can add an integer as an argument if you wish to limit the number of parallel processes

- Go on Biomart and download the annotation for your genome of interest, by respecting the following column order and using TAB separators:

	Chromosome/scaffold name	Gene start (bp)	Gene end (bp)	Gene name	Gene % GC content	Strand	Gene stable ID	Gene type

- Place the downloaded file in the root folder (same level as this readme) and rename it to "EnsemblGenesBadFormat.txt"
	Annotation.py tries to map Ensembl chromosome names to UCSC chromosome names in a very naive way: 
	1) adds "chrom" in front of chromosome names
	2) transforms the strand from -1/1 to -/+
	3) fills in empty gene names with the string "no_Ensembl_name"
	4) replaces spaces in gene names with underscores

- If you wish to control chromosome names (handy when working with scaffolds) and gene names yourself, you need to edit the annotation file you have downloaded. Make sure that:
	1) there are no spaces or tabs in the values of the annotation file
	2) tab formatting is still preserved
	3) you replace -1/1 strands with -/+
	4) you rename the annotation file to "EnsemblGenes.txt"
	5) you open "Annotation.py" and set naiveEnsemblToUCSCFormatting to False (line 56)

- Run mainMultiprocessing.py (> python mainMultiprocessing.py)

- Open Annotation.py in a text editor

- If you wish to include out of frame zinc finger arrays in downstream analysis: 
	set includeOutOfFrameZFAs = True
	specificy the minimum number of zinc fingers to be contained in the array for it to be considered for downstream analysis (minOFZFALength = 2 by default)

- Run Annotation.py

- Find results in the Tables folder, and bed files in the Bed folder.
	If you are awesome and are working with python, you may directly load the KZFP table as a pandas dataframe by pickle loading "KZFPTable.dat" from the root folder.
	Note that if want to reannotate a subset of chromosomes/scaffolds (e.g. to drop alternate chromosomes from a future KZFP homology comparison), you don't need to run mainMultiprocessing.py again.
	Just edit the text file listing chromosomes and run Annotation.py again

To compare KZFPs between species:
- Open SpeciesComparison.py
- Edit lines 10 and 11: species1 and species2 values are cosmetics and only influence the name of the final homology file
- Edit lines 12 and 13: path 1 and path2: indicate the position of the KZFPTable.dat file for each of the two species you want to compare.
	Note that if you use the same KZFPTable.dat,the script understands it as the same species and drops duplicates automatically 
	(e.g. does not score human ZNF221 to human ZNF221 homology to 1 but skips it instead)
- Edit line 14: canonicalThreshold: between 0 and 1
	KZFPs with a canonical score (proportion of canonical zinc fingers to all zinc fingers) equal to or lower than the threshold will be ignored. Handy to remove biases due to fully degenerated KZFPs
- Edit line 15: AlignSpacers: True or False
	Computes the alignment between spacers when set to True. Note that aligning spacers equates to alignming very long DNA sequences (>10'000kb), which is too heavy for a personal computer, even
	between a few KZFPs. 
- Run SpeciesComparison.py
- Find results in the Tables folder (csv) or in the SpeciesComparisons folder (pickle compressed .dat file)



