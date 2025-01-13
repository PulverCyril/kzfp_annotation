#!/home/bopekno/anaconda3/envs/python3.4/bin/python
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import translate
from Functions import sixFrameTranslation
from Functions import KZFPProfileSearch
from Functions import assembleZFArrays
from Functions import extractKRABs
from Functions import constructKZFPs
from Functions import KZFPListToBed
from Functions import singleZFListToBed
from Functions import zfArrayListToBed
from Functions import KRABToBed
from Functions import protDomainListToDNADomainList
from Functions import constructDNAKZFPs
from Functions import DNAKZFPListToBed
from Functions import assembleDNAZFArrays
from Functions import mergeDNAKZFPLists
from Functions import DNAzfArrayListToBed
from Functions import mergeDNAZFALists
from Functions import assembleOutOfFrameZFArrays
from Functions import extractKRABBBoxes
from Functions import KRABBboxListToBed
from Functions import addKRABBboxToKZFP
from Functions import extractKRABBDivs
from Functions import addKRABBDivToKZFP
from Functions import KRABBDivListToBed
from Functions import findZFArrayFingerprint
import pickle
import subprocess


def mainAnalysisForChromosome(params):
	chromosomeName, genome = params
	print("singleChromosomeAnalysis Started !")
	N = chromosomeName
	
	#reading the chromosome DNA file
	chromosome = SeqIO.read("./Chromosomes/{}.fa".format(N), "fasta")
	print(chromosome.id)

	#translating the chromosome DNA file the six possible reading frames and storing translations into fasta files
	print("Start translation...")
	proteins = sixFrameTranslation(chromosome, genome)
	for protein in proteins:
		SeqIO.write(protein, "./Translations/{}.fasta".format(protein.id), "fasta")
	print("Translation finished.")
	
	#Reading frames: FF = forward frame, RF = reverse frame, # = shift
	readingFrames = ["FF0", "FF1", "FF2", "RF0", "RF1", "RF2"]
	forwardReadingFrames = ["FF0", "FF1", "FF2"]
	reverseReadingFrames = ["RF0", "RF1", "RF2"]
	KZFPProtList = []
	zfArrayList = []
	singleZFList = []
	KRABList = []
	KRABBBoxList = []
	KRABBDivList = []
	for readingFrame in readingFrames:
		print("Search for domains with hmmpfam...")
		KZFPProfileSearch("./Translations/{}_{}.fasta".format(N, readingFrame))
		print("Domain search finished.")
		#parsing the hmmpfam result files
		print("Assembling KZFPs...")
		#the assembleZFarrays function produces the in frame zinc finger arrays (iZFAs)
		singleZincFingers, zincFingerArrays = assembleZFArrays("./hmmpfam/{}_{}_zf-C2H2.txt".format(N, readingFrame), "./Translations/{}_{}.fasta".format(N, readingFrame))
		KRABdomains = extractKRABs("./hmmpfam/{}_{}_KRAB.txt".format(N, readingFrame), "./Translations/{}_{}.fasta".format(N, readingFrame))
		KRABBBoxes = extractKRABBBoxes("./hmmpfam/{}_{}_KRAB_B.txt".format(N, readingFrame), "./Translations/{}_{}.fasta".format(N, readingFrame))
		KRABBDivs = extractKRABBDivs("./hmmpfam/{}_{}_KRAB_B_div.txt".format(N, readingFrame), "./Translations/{}_{}.fasta".format(N, readingFrame))
		zfArrayList += zincFingerArrays
		singleZFList += singleZincFingers
		KRABList += KRABdomains
		KRABBBoxList += KRABBBoxes
		KRABBDivList += KRABBDivs
		

	singleZFListDNA = protDomainListToDNADomainList (singleZFList, "./Chromosomes/{}.fa".format(N))
	KRABListDNA = protDomainListToDNADomainList(KRABList, "./Chromosomes/{}.fa".format(N))
	zfArrayListDNA = protDomainListToDNADomainList(zfArrayList, "./Chromosomes/{}.fa".format(N))
	KRABBBoxListDNA = protDomainListToDNADomainList(KRABBBoxList, "./Chromosomes/{}.fa".format(N))
	KRABBDivListDNA = protDomainListToDNADomainList(KRABBDivList, "./Chromosomes/{}.fa".format(N))


	#building KZFPs based on the out of frame zinc finger arrays (1000bp distance, not in frame, not overlapping, if available in a 40'000 bp interval: more than one single zinc finger, otherwise 1 single zinc finger is ok)
	KZFPDNAList = []
	zfArrayListDNACopy = list(zfArrayListDNA)
	#DNAzfArrayListToBed(zfArrayListDNA, "swagbefore.bed", "./Chromosomes/{}.fa".format(N))
	OOFZFArrays = assembleOutOfFrameZFArrays(zfArrayListDNACopy, singleZFListDNA, "./Chromosomes/{}.fa".format(N))
	#DNAzfArrayListToBed(zfArrayListDNA, "swag.bed", "./Chromosomes/{}.fa".format(N))
	KZFPDNAList = constructDNAKZFPs(KRABListDNA, OOFZFArrays, "./Chromosomes/{}.fa".format(N))

	#Adding a KRAB B Box (based on human KRAB B Boxes) to already built KZFPs
	addKRABBboxToKZFP(KZFPDNAList, KRABBBoxListDNA, KRABListDNA, OOFZFArrays)

	#Adding a KRAB b divergent box (human based) to already built KZFPs
	addKRABBDivToKZFP(KZFPDNAList, KRABBDivListDNA, KRABListDNA, OOFZFArrays)

	#finding fingerprints:
	findZFArrayFingerprint(zfArrayListDNA, singleZFListDNA)
	findZFArrayFingerprint(OOFZFArrays, singleZFListDNA)

	#Saving SeqRecord objects in binary objects for all the DNA converted sequences (for use in MutationsMain.py)
	#needed to store: single zinc fingers, in frame Zinc finger arrays, ouf of frame zinc finger arrays DNA built zinc fingers, KRABs, KZFPs
	with open("./Domains/{}_singleZFListDNA.dat".format(N), "wb") as f:
		pickle.dump(singleZFListDNA,f)
	with open("./Domains/{}_inFrameZFArrayListDNA.dat".format(N), "wb") as f:
		pickle.dump(zfArrayListDNA,f)
	with open("./Domains/{}_outOfFrameZFArrayList.dat".format(N), "wb") as f:
		pickle.dump(OOFZFArrays,f)
	with open("./Domains/{}_KRABListDNA.dat".format(N), "wb") as f:
		pickle.dump(KRABListDNA,f)
	with open("./Domains/{}_KZFPDNAList.dat".format(N), "wb") as f:
		pickle.dump(KZFPDNAList,f)
	with open("./Domains/{}_KRABBboxDNAList.dat".format(N), "wb") as f:
		pickle.dump(KRABBBoxListDNA,f)
	with open("./Domains/{}_KRABBDivDNAList.dat".format(N), "wb") as f:
		pickle.dump(KRABBDivListDNA,f)
	


	#Saving bed files to visualize results:
	#KZFPListToBed(KZFPProtList, "./Bed/{}_KZFPProtList.bed".format(N), "./Chromosomes/{}.fa".format(N))
	singleZFListToBed(singleZFList, "./Bed/{}_singleZFList.bed".format(N), "./Chromosomes/{}.fa".format(N))
	zfArrayListToBed(zfArrayListDNA, "./Bed/{}_inFrameZFArrayList.bed".format(N), "./Chromosomes/{}.fa".format(N))
	DNAzfArrayListToBed (OOFZFArrays, "./Bed/{}_outOfFrameZFArrayList.bed".format(N), "./Chromosomes/{}.fa".format(N))
	KRABToBed(KRABList, "./Bed/{}_KRABList.bed".format(N), "./Chromosomes/{}.fa".format(N))
	DNAKZFPListToBed(KZFPDNAList, "./Bed/{}_KZFPDNAList.bed".format(N))
	KRABBboxListToBed(KRABBBoxListDNA, "./Bed/{}_KRABBboxDNAList.bed".format(N))
	KRABBDivListToBed(KRABBDivListDNA, "./Bed/{}_KRABBDivDNAList.bed".format(N))

	#deleting temporary files:
	#dna chromosomes:
	subprocess.call("rm ./Chromosomes/{}.fa".format(N), shell=True)

	#translated chromosomes
	subprocess.call("rm ./Translations/{}*".format(N), shell=True)

	#hmmpfam text files:
	subprocess.call("rm ./hmmpfam/{}*".format(N), shell=True)

