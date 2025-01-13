import subprocess as ss
import pickle
from Functions import separateChromosomes
from Functions import mergeBedFiles
from Functions import mergeDomains
from Functions import getChromosomeNamesList
from Main import mainAnalysisForChromosome
import multiprocessing
#Cosmetics: genome name
genomePath = "/home/bopekno/Documents/Genomes/hg19test/hg19.fa"
chromListPath = "/home/bopekno/Documents/Genomes/hg19test/chromosomeNamesList"
genome = genomePath.split(".")[0].split("/")[-1]
#First argument: whole genome UCSC format fasta file. Second argument: UCSC formatted .chrom.sizes file (chromosomes with their names and sizes)


chromosomeList = separateChromosomes(genomePath, chromListPath)

#use this next line instead of the one above if the chromosomes have already been seperated and you don't need to separate them again
#chromosomeList = getChromosomeNamesList("/home/bopekno/Documents/Genomes/mm9/mm9.chrom.sizes")

#parameters must be zipped, thus there will be as many parameters as rows in KZFPTable1, and the other elements will remain the same.
nParam = len(chromosomeList)
param2 = [genome]*nParam
params = zip(chromosomeList, param2)
if __name__ == '__main__':
    #you can add an integer as an argument to multiprocessing.Pool() to limit the number of parallel processes (if you're using your computer for something else)
    pool = multiprocessing.Pool()

    #results is a list of pd dataframes. It needs to be concatenated
    pool.map(mainAnalysisForChromosome, params)

print("merging bed files and domain files")
#merging bed files
bed_domains = ["KZFPDNAList", "KRABList", "singleZFList", "inFrameZFArrayList", "outOfFrameZFArrayList", "KRABBboxDNAList", "KRABBDivDNAList"]
params2 = zip([chromosomeList]*len(bed_domains), bed_domains)

if __name__ == '__main__':
    #you can add an integer as an argument to multiprocessing.Pool() to limit the number of parallel processes (if you're using your computer for something else)
    pool = multiprocessing.Pool()

    #results is a list of pd dataframes. It needs to be concatenated
    pool.map(mergeBedFiles, params2)

#merging domain files
domains = ["inFrameZFArrayListDNA", "KRABBboxDNAList", "KRABBDivDNAList", "KRABListDNA", "KZFPDNAList", "outOfFrameZFArrayList", "singleZFListDNA"]
params3 = zip([chromosomeList]*len(domains), domains)

if __name__ == '__main__':
    #you can add an integer as an argument to multiprocessing.Pool() to limit the number of parallel processes (if you're using your computer for something else)
    pool = multiprocessing.Pool()

    #results is a list of pd dataframes. It needs to be concatenated
    pool.map(mergeDomains, params3)
