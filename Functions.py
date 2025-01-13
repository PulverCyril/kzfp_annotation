import os
import subprocess
import re
import pickle
import os.path
import operator
import multiprocessing
from linecache import getline
from bisect import bisect
from bisect import bisect_left
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import translate
from Bio import SearchIO
from Bio.Alphabet import generic_dna
from Bio.Seq import MutableSeq
from Bio import AlignIO
import vcf
import numpy as np
#import seaborn as sns
import matplotlib.pyplot as plt
import pickle
import pandas as pd
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.pairwise2 import format_alignment

def sixFrameTranslation (sequence, genome):
	"""Translation in all six frames, using the standard table (index = 1).
	Input: SeqRecord DNA object, output: list of SeqRecord protein objects, first foward frame 0, 1, 2 then reverse frame 0, 1, 2
	in each SeqRecord.id: sequence.id (id of the input SeqRecord object) followed by _FF0 (for roward frame 0) (or RF0) or 1, 2"""

	rawTranslation =[]
	"""i = reading frame index, the standard codon table is used and the principle to move the reading frame to the right is to start reading the sequence at nucleotide 0, then 1 then 2
	Reverse reading frames are generated exactly in the same way, with a reverse complement function beforehand."""
	i=0
	while i < 3:
		#this sequence compensates the total length of the translated sequence for the missing nucleotide, somehow having 1 difference in length of nucleotide caused
		#HMMER3SCAN to get infinite bit scores everytime...
		seqToAdd = str()
		if i == 1:
			seqToAdd = "N"
		if i == 2:
			seqToAdd = "NN"
		rawTranslation.append(SeqRecord(translate(sequence.seq[i:]+seqToAdd, table = "Standard", stop_symbol = "*", to_stop = False, gap = None), \
		id="{}_FF{}".format(sequence.id, i), name = "{}_FF{}".format(sequence.id, i), description = genome))

		#print(translate(sequence.seq[i:], table = "Standard", stop_symbol = "*", to_stop = False, gap = None))
		i+=1

	#Three reverse strands:
	reverseSequence = sequence.reverse_complement()
	i=0
	while i < 3:
		seqToAdd = str()
		if i == 1:
			seqToAdd = "N"
		if i == 2:
			seqToAdd = "NN"
		rawTranslation.append(SeqRecord(translate(reverseSequence.seq[i:]+seqToAdd, table = "Standard", stop_symbol = "*", to_stop = False, gap = None), \
		id="{}_RF{}".format(sequence.id, i), name = "{}_RF{}".format(sequence.id, i), description=genome))
		
		i+=1
	return rawTranslation

def KZFPProfileSearch(filePath):
	"""passes through various unix commands to create a new hmmpfam file from C2H2 and KRAB domain seed sequences, 
	which is used to search the query protein sequence. The filepath of this protein is given as a string argument"""
	#Building the hmm profiles is only required if the KZFPfam file does not alreay exist, which is generally the case,
	if not os.path.isfile("KRAB.hmm"):
		subprocess.call("hmmbuild KRAB.hmm ./PFAM/KRAB/KRAB_seed.sto", shell=True)
	if not os.path.isfile("zf-C2H2.hmm"):
		subprocess.call("hmmbuild zf-C2H2.hmm ./PFAM/zf-C2H2/zf-C2H2_seed.sto", shell=True)
	if not os.path.isfile("KRAB_B.hmm"):
		#run MUSCLE on the sequences of KZFPs containing a KRAB B box according to literature:
		subprocess.call("muscle -in ./KRAB_B_box/KRAB_B_box_KZFPs.fasta -out ./KRAB_B_box/KRAB_B_box_KZFPs.clw -clwstrict", shell = True)
		#import the resulting clustalw aligment file into biopython:
		KRAB_B_boxAlignment = AlignIO.read("./KRAB_B_box/KRAB_B_box_KZFPs.clw", "clustal")
		#slice the alignment file at the wanted position (this is actually a parameter which has to be specified, 108 to 124)
		KRAB_B_boxAlignment = KRAB_B_boxAlignment[:, 108:124]
		#write the sliced alignment containing only the KRAB B box to a clustal file so that hmm can build on it
		AlignIO.write(KRAB_B_boxAlignment, "./KRAB_B_box/KRAB_B_box.clw", "clustal")
		#feeding the alignment to hmmbuild
		subprocess.call("hmmbuild KRAB_B.hmm ./KRAB_B_box/KRAB_B_box.clw", shell=True)

	if not os.path.isfile("KRAB_B_div.hmm"):
		#run MUSCLE on the sequences of KZFPs containing a KRAB B box according to literature:
		subprocess.call("muscle -in ./KRAB_B_div/KRAB_B_div_KZFPs.fasta -out ./KRAB_B_div/KRAB_B_div_KZFPs.clw -clwstrict", shell = True)
		#import the resulting clustalw aligment file into biopython:
		KRAB_B_divAlignment = AlignIO.read("./KRAB_B_div/KRAB_B_div_KZFPs.clw", "clustal")
		#slice the alignment file at the wanted position (this is actually a parameter which has to be specified, 108 to 124)
		KRAB_B_divAlignment = KRAB_B_divAlignment[:, 69:102]
		#write the sliced alignment containing only the KRAB B box to a clustal file so that hmm can build on it
		AlignIO.write(KRAB_B_divAlignment, "./KRAB_B_div/KRAB_B_div.clw", "clustal")
		#feeding the alignment to hmmbuild
		subprocess.call("hmmbuild KRAB_B_div.hmm ./KRAB_B_div/KRAB_B_div.clw", shell=True)

	#executing the profile search and storing it in a new file:

	fullFilePath = str(filePath)
	a = './Translations/'
	b = '.fasta'

	ChromosomeName = fullFilePath.split(a)[-1].split(b)[0]
	
	#searching for ZFPs
	subprocess.call("hmmpfam zf-C2H2.hmm {} > ./hmmpfam/{}_zf-C2H2.txt".format(filePath, ChromosomeName), shell = True)
	#searching for KRABs:
	subprocess.call("hmmpfam KRAB.hmm {} > ./hmmpfam/{}_KRAB.txt".format(filePath, ChromosomeName), shell = True)
	#searching for KRAB B Boxes:
	subprocess.call("hmmpfam KRAB_B.hmm {} > ./hmmpfam/{}_KRAB_B.txt".format(filePath, ChromosomeName), shell = True)

	subprocess.call("hmmpfam KRAB_B_div.hmm {} > ./hmmpfam/{}_KRAB_B_div.txt".format(filePath, ChromosomeName), shell = True)

	#subprocess.call("hmmscan KZFPfam {}".format(filePath), shell = True)

def assembleZFArrays(hmmpfamOutput, chromosomeQuery):
	"""
	takes the text file output from a hmmpfam command and returns a list of zinc finger SeqcRecord objects as single zinc fingers and a list of 
	zinc finger arrays. Sequences are retrieved from chromosomeQuery reading frame."""
	#opening file
	chromosome = SeqIO.read(chromosomeQuery, "fasta")
	#storing the chromosome and reading frame as chr#_RF#
	chromRF = chromosome.id
	output = open(hmmpfamOutput, "r")
	#getting the number of hits in the sequence to know how many lines to read: it's on line index 16
	content = output.readlines()
	# you may also want to remove whitespace characters like `\n` at the end of each line
	content = [x.strip() for x in content] 
	#check if there are results here
	#print(content)
	line16 = content[16]
	line16list = line16.split(" ")
	numberOfHits = int(line16list[-1])
	spacer = re.compile("\s+")
	#going through all the hit lines, checking that the bitscore is above threshold, and if so, adding a SeqRecord object to the singleZFList object
	hitNumber = 0
	index = 0
	singleZFList = []
	while hitNumber < numberOfHits:
		currentLine = content[hitNumber+21]
		hitNumber+=1
		lineList = spacer.split(currentLine)
		
		#checking score threshold:
		if float(lineList[8])>=0:
			#hmmerpfam returns indexes in one-based coordinate, so the start coordinate needs to be reduced by one
			start = int(lineList[2])-1
			end = int(lineList[3])
			singleZFList.append(SeqRecord(chromosome.seq[start:end], id="{}_sZF_{}".format(chromRF, index), name = "{}_sZF_{}".format(chromRF, index), \
		description = "single zinc finger", annotations = {"AAStart": start, "AAEnd": end, "chromosome": chromosome.id[:-4], "Genome": chromosome.description.split()[-1]}))
		index += 1
	output.close()
	

	#building new ZFP arrays:
	zfArrayList = []
	skippedIndex = -1
	zfArrayCounter = 0
	for index, singleZF in enumerate(singleZFList):
		#At the start, the array candidate is simply the first zinc finger alone
		#If the zinc finger array is several zinc fingers long, no need to do all the "sub arrays" of the main array -> we must skip some indexes
		
		if index > skippedIndex:
			#This list contains all the ids of the contained single zinc fingers, and will be stored in the annotations attribute
			zfContained = [singleZFList[index].id]
			
			endOfArrayReached = False
			i = 0
			while endOfArrayReached == False:
				#if condition not to get out of range of the list
				
				if index+i+1 <= len(singleZFList)-1:
					
					#we check whether the start of the next single zinc finger starts 4-6 AA from the end of the previous zinc finger
					#in an alternate version, we can assemble single zinc fingers which are within a distance of 1kb from one another in a single array
					#to mitigate the skewing effect of zinc finger array encoded by exons seperated by an intron
					#print("Distance: {}".format(singleZFList[index+i+1].annotations["AAStart"] - singleZFList[index+i].annotations["AAEnd"]))
					if 1 <= singleZFList[index+i+1].annotations["AAStart"] - singleZFList[index+i].annotations["AAEnd"] <= 333:
						#print("Length is ok")
						zfContained.append(singleZFList[index+i+1].id)
						i += 1
						
					else:
						#print("End of array reached")
						endOfArrayReached = True
						
				else:
					endOfArrayReached = True
				#added condition to build a zinc finger array: it must contain at least two zinc fingers:
				if endOfArrayReached == True:
					#the ZF array sequence is the chromosome sequence, start AA = start of the first sequence, stop AA = end of the last sequence 
						newSequence = chromosome.seq[singleZF.annotations["AAStart"]:singleZFList[index+i].annotations["AAEnd"]]
						zfArrayList.append(SeqRecord(newSequence, id="{}_protZFA_{}".format(chromosome.id, zfArrayCounter), name="{}_protZFA_{}".format(chromosome.id, zfArrayCounter), \
						description = "zinc finger array built by the protein approach", annotations = {"AAStart": singleZF.annotations["AAStart"], "AAEnd": singleZFList[index+i].annotations["AAEnd"], "singleZFContained": zfContained, "chromosome": chromosome.id[:-4], "Genome": chromosome.description.split()[-1]}))
						skippedIndex = index+i
						zfArrayCounter +=1

	return singleZFList, zfArrayList

#old version, not used anymore (check assembleOutOfFrameZFArrays)
def assembleDNAZFArrays(singleZFList, chromosomeDNA):
	"""This function is an alternative version of the assembleZFArrays function, in the sense that it builds ZFarrays based on DNA distance irrespective of reading frames, with Conditions
	that they are not further away than 1000 bp and that they are on the same strand. Of course since those ZF arrays are not in frame, they don't have an amino acid sequence."""
	#converting the protein coordinates of the single zinc finger list into DNA coordinates:
	chromosome = SeqIO.read(chromosomeDNA, "fasta")
	chromosomeLengthDNA = len(str(chromosome.seq))
	for singleZF in singleZFList:
		singleZF.annotations["Shift"], singleZF.annotations["Strand"] = getShiftStrand(singleZF)
		singleZF.annotations["DNAStart"], singleZF.annotations["DNAEnd"] = proteinToDNACoordinates(singleZF, chromosomeLengthDNA)

	#seperating forward strand single zinc fingers from reverse strand zinc fingers:
	FFzfDNAList = [singleZF for singleZF in singleZFList if singleZF.annotations["Strand"] == "+"]
	FFzfDNAList = sorted(FFzfDNAList, key=lambda array: array.annotations["DNAStart"])
	print([singleZF.id for singleZF in FFzfDNAList])

	RFzfDNAList = [singleZF for singleZF in singleZFList if singleZF.annotations["Strand"] == "-"]
	RFzfDNAList = sorted(RFzfDNAList, key=lambda array: array.annotations["DNAStart"])
	print([singleZF.id for singleZF in RFzfDNAList])
	#building DNAZFarrays based on bp distance: 
	#building new ZFP arrays:
	#foward strand:
	zfArrayList = []
	skippedIndex = -1
	zfArrayCounter = 0
	for index, singleZF in enumerate(FFzfDNAList):
		#At the start, the array candidate is simply the first zinc finger alone
		#If the zinc finger array is several zinc fingers long, no need to do all the "sub arrays" of the main array -> we must skip some indexes
		
		if index > skippedIndex:
			#This list contains all the ids of the contained single zinc fingers, and will be stored in the annotations attribute
			zfContained = [FFzfDNAList[index].id]
			
			endOfArrayReached = False
			i = 0
			while endOfArrayReached == False:
				#if condition not to get out of range of the list
				
				if index+i+1 <= len(FFzfDNAList)-1:
					
					#the DNA distance must be 1000bp or less, and more than 1 bp to avoid overlap
					print(FFzfDNAList[index+i+1].annotations["DNAStart"] - FFzfDNAList[index+i].annotations["DNAEnd"])
					if 1 <= FFzfDNAList[index+i+1].annotations["DNAStart"] - FFzfDNAList[index+i].annotations["DNAEnd"] <= 1000:
						#print("Length is ok")
						zfContained.append(FFzfDNAList[index+i+1].id)
						print(zfContained)
						i += 1
						
					else:
						#print("End of array reached")
						endOfArrayReached = True
						
				else:
					endOfArrayReached = True
				#added condition to build a zinc finger array: it must contain at least two zinc fingers:
				if endOfArrayReached == True and len(zfContained)>1:
					#the ZF array sequence is the chromosome sequence, start AA = start of the first sequence, stop AA = end of the last sequence 
						newSequence = chromosome.seq[singleZF.annotations["DNAStart"]:FFzfDNAList[index+i].annotations["DNAEnd"]]
						zfArrayList.append(SeqRecord(newSequence, id="{}_FF_DNAZFA_{}".format(chromosome.id, zfArrayCounter), name="{}_FF_DNAZFA_{}".format(chromosome.id, zfArrayCounter), \
						description = "DNA assembled zinc finger array", annotations = {"DNAStart": singleZF.annotations["DNAStart"], "DNAEnd": FFzfDNAList[index+i].annotations["DNAEnd"], \
						"Strand": singleZF.annotations["Strand"], "singleZFContained": zfContained, "DNASequence": newSequence}))
						skippedIndex = index+i
						zfArrayCounter +=1
	#reverse strand:
	skippedIndex = -1
	for index, singleZF in enumerate(RFzfDNAList):
		#At the start, the array candidate is simply the first zinc finger alone
		#If the zinc finger array is several zinc fingers long, no need to do all the "sub arrays" of the main array -> we must skip some indexes
		
		if index > skippedIndex:
			#This list contains all the ids of the contained single zinc fingers, and will be stored in the annotations attribute
			zfContained = [RFzfDNAList[index].id]
			
			endOfArrayReached = False
			i = 0
			while endOfArrayReached == False:
				#if condition not to get out of range of the list
				
				if index+i+1 <= len(RFzfDNAList)-1:
					
					#the DNA distance must be 1000bp or less, and more than 1 bp to avoid overlap
					if 1 <= RFzfDNAList[index+i+1].annotations["DNAStart"] - RFzfDNAList[index+i].annotations["DNAEnd"] <= 1000:
						#print("Length is ok")
						zfContained.append(RFzfDNAList[index+i+1].id)
						i += 1
						
					else:
						#print("End of array reached")
						endOfArrayReached = True
						
				else:
					endOfArrayReached = True
				#added condition to build a zinc finger array: it must contain at least two zinc fingers:
				if endOfArrayReached == True and len(zfContained)>1:
					#the ZF array sequence is the chromosome sequence, start AA = start of the first sequence, stop AA = end of the last sequence 
						newSequence = chromosome.seq[singleZF.annotations["DNAStart"]:RFzfDNAList[index+i].annotations["DNAEnd"]]
						zfArrayList.append(SeqRecord(newSequence, id="{}_RF_DNAZFA_{}".format(chromosome.id, zfArrayCounter), name="{}_RF_DNAZFA_{}".format(chromosome.id, zfArrayCounter), \
						description = "DNA assembled zinc finger array", annotations = {"DNAStart": singleZF.annotations["DNAStart"], "DNAEnd": RFzfDNAList[index+i].annotations["DNAEnd"], \
						"Strand": singleZF.annotations["Strand"], "singleZFContained": zfContained, "DNASequence": newSequence}))
						skippedIndex = index+i
						zfArrayCounter +=1
	return zfArrayList



def assembleOutOfFrameZFArrays(InFrameZFAList, singleZincFingerList, chromosomeDNA):
	"""This function is very similar to assembleDNAZFArrays, but it assembles zinc finger arrays which are in distance of 1kb, and not single zinc fingers. This allows to keep track of in frame zinc finger arrays contained
	within the out of frame zinc finger array, for a more efficient prediction of the outcome of genetic changes."""
	#converting the protein coordinates of the single zinc finger list into DNA coordinates:
	chromosome = SeqIO.read(chromosomeDNA, "fasta")
	chromosomeLengthDNA = len(str(chromosome.seq))
	#inFrameZFAListCopy = list(InFrameZFAList)
	for InFrameZFA in InFrameZFAList:
		InFrameZFA.annotations["Shift"], InFrameZFA.annotations["Strand"] = getShiftStrand(InFrameZFA)
		InFrameZFA.annotations["DNAStart"], InFrameZFA.annotations["DNAEnd"] = proteinToDNACoordinates(InFrameZFA, chromosomeLengthDNA)

	#seperating forward strand in frame zinc finger arrays from reverse strand in frame zinc finger arrays:
	FFZFADNAList = [InFrameZFA for InFrameZFA in InFrameZFAList if InFrameZFA.annotations["Strand"] == "+"]
	FFZFADNAList = sorted(FFZFADNAList, key=lambda array: array.annotations["DNAStart"])

	RFZFADNAList = [InFrameZFA for InFrameZFA in InFrameZFAList if InFrameZFA.annotations["Strand"] == "-"]
	RFZFADNAList = sorted(RFZFADNAList, key=lambda array: array.annotations["DNAStart"])
	#building DNAZFarrays based on bp distance between already detected in frame ZFAs: 
	#building new ZFP arrays:
	#foward strand:
	zfArrayList = []
	skippedIndex = -1
	zfArrayCounter = 0
	for index, inFrameZFA in enumerate(FFZFADNAList):
		#At the start, the array candidate is simply the first in frame zinc finger array alone
		#If the out of frame zinc finger array encompasses several in frame zinc finger arrays, no need to do all the "sub arrays" of the main array -> we must skip some indexes
		
		if index > skippedIndex:
			#This list contains all the ids of the contained single zinc fingers, and will be stored in the annotations attribute
			zfContained = list(FFZFADNAList[index].annotations["singleZFContained"])
			#This list contains all the ids of the contained in frame zinc finger arrays, and will be stored in the annotations attribute
			inFrameZFAContained = [FFZFADNAList[index].id]
			DNAEndList = [FFZFADNAList[index].annotations["DNAEnd"]]
			endOfArrayReached = False
			i = 0
			while endOfArrayReached == False:
				#if condition not to get out of range of the list
				
				if index+i+1 <= len(FFZFADNAList)-1:
					
					#the DNA distance must be 1000bp or less, but on the opposite of building in frame KZFPs, overlap is allowed here. (thus the <= 1 condition is gone)
					#another complication: sometimes (like for ZNF726): the DNAEnd coordinate of the last ZFA added to the out of frame ZFA is not actually the end coordinate of the out of frame ZFA.
					#To solve this problem, DNAEnd coordinates must be stored and sorted at every round.
					#if FFZFADNAList[index+i+1].annotations["DNAStart"] - FFZFADNAList[index+i].annotations["DNAEnd"] <= 1000:
					if FFZFADNAList[index+i+1].annotations["DNAStart"] - max(DNAEndList) <= 1000:
						#print("Length is ok")
						inFrameZFAContained.append(FFZFADNAList[index+i+1].id)
						zfContained += FFZFADNAList[index+i+1].annotations["singleZFContained"]
						DNAEndList.append(FFZFADNAList[index+i+1].annotations["DNAEnd"])
						i += 1
						
					else:
						#print("End of array reached")
						endOfArrayReached = True
						
				else:
					endOfArrayReached = True
				#added condition to build a zinc finger array: it must contain at least two zinc fingers:
				if endOfArrayReached == True and len(zfContained)>0:
					#The out of frame zinc finger arrays does not have an amino acid sequence, as it might contain by definition different reading frames. It is only defines by a DNA sequence.
					#To overcome the problem of downstream mutation fate, the out of frame zinc finger array contains a list of in frame zinc finger arrays, to make mutation predictions possible.
					#Overlapping zinc finger arrays cause an additionnal problem, which is that the DNA start and DNA end coordinates are not always those of the first and last zinc finger array found to belong in the sequence
					#(for example in the case where one zinc finger array is completely contained inside the boundaries of another zinc finger array, which is on the same strand but on a different reading frame)
					DNAStartList = [ZFA.annotations["DNAStart"] for ZFA in FFZFADNAList if ZFA.id in inFrameZFAContained]
					DNAEndList = [ZFA.annotations["DNAEnd"] for ZFA in FFZFADNAList if ZFA.id in inFrameZFAContained]
					DNAStart = min(DNAStartList)
					DNAEnd = max(DNAEndList)
					newSequence = chromosome.seq[DNAStart:DNAEnd+1]

					#lastly, the list of single zinc fingers have to be reordered based on their DNAStart. Ordering by in frame zinc finger array DNAStart is 
					#not sufficient, as sometimes, a smaller in frame zinc finger array starts after but ends before a longer in frame zinc finger array. This causes
					#the single zinc fingers of shorter zinc finger array to be placed after those of the longer zinc finger array.
					#The ordering of the single zinc fingers depends on the reading frame (forward: sort by DNAStart, increasing, reverse: DNAEnd, decreasing)
					

					#forward strand:
					sZFDNAStartList = [(sZF.id, sZF.annotations["DNAStart"]) for sZF in [sZF for sZF in singleZincFingerList if sZF.id in zfContained]]
					#reordering:
					sZFDNAStartList.sort(key = lambda tup: tup[1])

					#creating a new list of id:
					zfContained = [sZF[0] for sZF in sZFDNAStartList]


					zfArrayList.append(SeqRecord(newSequence, id="{}_FF_OOFZFA_{}".format(chromosome.id, zfArrayCounter), name="{}_FF_OOFZFA_{}".format(chromosome.id, zfArrayCounter), \
					description = "Out of frame assembled zinc finger array", annotations = {"DNAStart": DNAStart, "DNAEnd": DNAEnd, \
					"Strand": inFrameZFA.annotations["Strand"], "singleZFContained": zfContained, "inFrameZFAContained": inFrameZFAContained, "DNASequence": newSequence, "chromosome": chromosome.id, "Genome": inFrameZFA.annotations["Genome"]}))
					
					skippedIndex = index+i
					zfArrayCounter +=1
					
	#reverse strand:
	skippedIndex = -1
	for index, inFrameZFA in enumerate(RFZFADNAList):
		
		if index > skippedIndex:

			zfContained = list(RFZFADNAList[index].annotations["singleZFContained"])
			inFrameZFAContained = [RFZFADNAList[index].id]
			DNAEndList = [RFZFADNAList[index].annotations["DNAEnd"]]
			endOfArrayReached = False
			i = 0
			while endOfArrayReached == False:

				
				if index+i+1 <= len(RFZFADNAList)-1:
					
					if RFZFADNAList[index+i+1].annotations["DNAStart"] - max(DNAEndList) <= 1000:
						inFrameZFAContained.append(RFZFADNAList[index+i+1].id)
						zfContained += RFZFADNAList[index+i+1].annotations["singleZFContained"]
						DNAEndList.append(RFZFADNAList[index+i+1].annotations["DNAEnd"])
						i += 1
						
					else:
						endOfArrayReached = True
						
				else:
					endOfArrayReached = True
				if endOfArrayReached == True and len(zfContained)>0:
					DNAStartList = [ZFA.annotations["DNAStart"] for ZFA in RFZFADNAList if ZFA.id in inFrameZFAContained]
					DNAEndList = [ZFA.annotations["DNAEnd"] for ZFA in RFZFADNAList if ZFA.id in inFrameZFAContained]
					DNAStart = min(DNAStartList)
					DNAEnd = max(DNAEndList)
					newSequence = chromosome.seq[DNAStart:DNAEnd+1]

					#reverse strand:
					sZFDNAStartList = [(sZF.id, sZF.annotations["DNAEnd"]) for sZF in [sZF for sZF in singleZincFingerList if sZF.id in zfContained]]
					#reordering:
					sZFDNAStartList.sort(key = lambda tup: tup[1], reverse = True)

					#creating a new list of id:
					zfContained = [sZF[0] for sZF in sZFDNAStartList]
					
					zfArrayList.append(SeqRecord(newSequence, id="{}_RF_OOFZFA_{}".format(chromosome.id, zfArrayCounter), name="{}_RF_OOFZFA_{}".format(chromosome.id, zfArrayCounter), \
					description = "Out of frame assembled zinc finger array", annotations = {"DNAStart": DNAStart, "DNAEnd": DNAEnd, \
					"Strand": inFrameZFA.annotations["Strand"], "singleZFContained": zfContained, "inFrameZFAContained": inFrameZFAContained, "DNASequence": newSequence, "chromosome": chromosome.id, "Genome": inFrameZFA.annotations["Genome"]}))
					
					skippedIndex = index+i
					zfArrayCounter +=1
	return zfArrayList


def extractKRABs(hmmpfamOutput, chromosomeQuery):
	"""takes a hmmpfam output file and a translated chromosome and returns a list containing the KRAB domains found by hmmpfam with a threshold of 13"""
		#opening file
	chromosome = SeqIO.read(chromosomeQuery, "fasta")
	#storing the chromosome and reading frame as chr#_RF#
	chromRF = chromosome.id
	output = open(hmmpfamOutput, "r")
	#getting the number of hits in the sequence to know how many lines to read: it's on line index 16
	content = output.readlines()
	# you may also want to remove whitespace characters like `\n` at the end of each line
	content = [x.strip() for x in content] 
	#print(content)
	line16 = content[16]
	line16list = line16.split(" ")
	numberOfHits = int(line16list[-1])
	spacer = re.compile("\s+")
	#going through all the hit lines, checking that the bitscore is above threshold, and if so, adding a SeqRecord object to the KRABList object
	hitNumber = 0
	index = 0
	KRABList = []
	while hitNumber < numberOfHits:
		currentLine = content[hitNumber+21]
		hitNumber+=1
		lineList = spacer.split(currentLine)
		
		#checking score threshold:
		if float(lineList[8])>=13:
			start = int(lineList[2])-1
			end = int(lineList[3])
			KRABList.append(SeqRecord(chromosome.seq[start:end], id="{}_KRAB_{}".format(chromRF, index), name = "{}_KRAB_{}".format(chromRF, index), \
		description = "KRAB domain", annotations = {"AAStart": start, "AAEnd": end, "chromosome": chromosome.id[:-4], "Genome": chromosome.description.split()[-1]}))
		index += 1
	output.close()
	return KRABList

def extractKRABBBoxes(hmmpfamOutput, chromosomeQuery):
	"""takes a hmmpfam output file and a translated chromosome and returns a list containing the KRAB B box domains found by hmmpfam with a threshold of 13"""
		#opening file
	chromosome = SeqIO.read(chromosomeQuery, "fasta")
	#storing the chromosome and reading frame as chr#_RF#
	chromRF = chromosome.id
	output = open(hmmpfamOutput, "r")
	#getting the number of hits in the sequence to know how many lines to read: it's on line index 16
	content = output.readlines()
	# you may also want to remove whitespace characters like `\n` at the end of each line
	content = [x.strip() for x in content] 
	#print(content)
	line16 = content[16]
	line16list = line16.split(" ")
	numberOfHits = int(line16list[-1])
	spacer = re.compile("\s+")
	#going through all the hit lines, checking that the bitscore is above threshold, and if so, adding a SeqRecord object to the KRABList object
	hitNumber = 0
	index = 0
	KRABBboxList = []
	while hitNumber < numberOfHits:
		currentLine = content[hitNumber+21]
		hitNumber+=1
		lineList = spacer.split(currentLine)
		
		#checking score threshold:
		if float(lineList[8])>=13:
			start = int(lineList[2])-1
			end = int(lineList[3])
			KRABBboxList.append(SeqRecord(chromosome.seq[start:end], id="{}_KRABBbox_{}".format(chromRF, index), name = "{}_KRABBbox_{}".format(chromRF, index), \
		description = "KRABBbox domain", annotations = {"AAStart": start, "AAEnd": end, "chromosome": chromosome.id[:-4], "Genome": chromosome.description.split()[-1]}))
		index += 1
	output.close()
	return KRABBboxList

def extractKRABBDivs(hmmpfamOutput, chromosomeQuery):
	"""takes a hmmpfam output file and a translated chromosome and returns a list containing the KRAB B divergent domains found by hmmpfam with a threshold of 13"""
		#opening file
	chromosome = SeqIO.read(chromosomeQuery, "fasta")
	#storing the chromosome and reading frame as chr#_RF#
	chromRF = chromosome.id
	output = open(hmmpfamOutput, "r")
	#getting the number of hits in the sequence to know how many lines to read: it's on line index 16
	content = output.readlines()
	# you may also want to remove whitespace characters like `\n` at the end of each line
	content = [x.strip() for x in content] 
	#print(content)
	line16 = content[16]
	line16list = line16.split(" ")
	numberOfHits = int(line16list[-1])
	spacer = re.compile("\s+")
	#going through all the hit lines, checking that the bitscore is above threshold, and if so, adding a SeqRecord object to the KRABList object
	hitNumber = 0
	index = 0
	KRABBDivList = []
	while hitNumber < numberOfHits:
		currentLine = content[hitNumber+21]
		hitNumber+=1
		lineList = spacer.split(currentLine)
		
		#checking score threshold:
		if float(lineList[8])>=13:
			start = int(lineList[2])-1
			end = int(lineList[3])
			KRABBDivList.append(SeqRecord(chromosome.seq[start:end], id="{}_KRABBDiv_{}".format(chromRF, index), name = "{}_KRABBDiv_{}".format(chromRF, index), \
		description = "KRABBDiv domain", annotations = {"AAStart": start, "AAEnd": end, "chromosome": chromosome.id[:-4], "Genome": chromosome.description.split()[-1]}))
		index += 1
	output.close()
	return KRABBDivList



def constructKZFPs (KRABdomainList, zfArrayList, chromosomeQuery):
	"""This function constructs putative krab zinc fingers in a naive way by finding the closest zfp array downstream of each KRAB domain and regrouping them in a single sequence 
	by reading in the chromosomeQuery full chromosome translation sequence. Works only on one of the six protein reading frame !"""

	if len(KRABdomainList) == 0 or len(zfArrayList) == 0:
		return []
	#loading the chromosome (the entry query for HMMER3SCAN):
	chromosome = SeqIO.read(chromosomeQuery, "fasta")

	#building a list of the same size and same indexing as zfArrayList, but only containing the AAstart value:
	zfAAStartList =[array.annotations["AAStart"] for array in zfArrayList]

	KZFPList = []
	KZFPcounter = 0
	for KRAB in KRABdomainList:
		#finds the next zinc finger after the end amino acid of the KRAB domain. Bisect returns the next index in the list which value is higher than a limit value.
		#If the limit value is above the highest value of the list, it returns an index = len(zfAAStartList) (thus pointing out of the list). This situation can 
		#be handled with a simple if condition, which if not respected yields no krab zinc finger.
		nextArrayIndex = bisect(zfAAStartList, KRAB.annotations["AAEnd"])

		if nextArrayIndex != len(zfArrayList):
			#constructing the sequnece of the KZFP: going from the first AA of the KRAB domain to the last AA of the zf array
	 		newSequence = chromosome.seq[KRAB.annotations["AAStart"]:zfArrayList[nextArrayIndex].annotations["AAEnd"]]
	 		KZFPList.append(SeqRecord(newSequence, id = "{}_KZFP_{}".format(chromosome.id, KZFPcounter), name = "{}_KZFP_{}".format(chromosome.id, KZFPcounter), description = "KRAB zinc finger",  \
	 		annotations = {"AAStart": KRAB.annotations["AAStart"], "AAEnd": zfArrayList[nextArrayIndex].annotations["AAEnd"], "KRAB": KRAB.id, "zfArrayContained": zfArrayList[nextArrayIndex].id,  \
	 		"zfStart":  zfArrayList[nextArrayIndex].annotations["AAStart"], "singleZFContained": zfArrayList[nextArrayIndex].annotations["singleZFContained"]}))
	 		KZFPcounter += 1
	return KZFPList

def getChromosome(sequence):
	return sequence.annotations["chromosome"]

def getShiftStrand(sequence):
	"""Takes a protein SeqRecord object which .id attribute is of syntax chr#_FF#"and returns the shift corresponding to the reading frame
	and the DNA strand as a tuple"""
	chromosomeName = getChromosome(sequence)
	readingFrames = {"FF0": (0, "+"), "FF1": (1, "+"), "FF2": (2, "+"), "RF0": (0, "-"), "RF1": (1, "-"), "RF2": (2, "-")}
	strIndex = len(chromosomeName)+1
	print(sequence.id)
	print(sequence.id[strIndex:strIndex+3])
	return readingFrames[sequence.id[strIndex:strIndex+3]]

def proteinToDNACoordinates(proteinSequence, chromosomeLengthDNA):
	"""Takes a SeqRecord protein object and returns the corresponding the nucleotide start and end coordinates as a tuple. Of course the chromosome length must correspond to the
	chromsome of origin of the protein sequence. The toBed argument is given as a bool, and if == True adds 1 to the end coordinate for proper display in a .bed file"""
	shift, strand = getShiftStrand(proteinSequence)
	#Because it is needed for calculation of the start and end coordinates for the reverse strand, we recuperate the length 
		#of the appropriate chromosome:
	start = 0
	end = 0
	#if the protein was translated from the forward strand, then the start of the AA corresponds to the start of the DNA
	if strand == "+":
		start = (proteinSequence.annotations["AAStart"])*3+shift
		#needs -1, because if start is 0 on a modulo 3, end should be 2 on a modulo 3, and not 0 again
		end = (proteinSequence.annotations["AAEnd"])*3+shift - 1
	#if on the contrary the protein was translated from the reverse strand, then the start of the AA corresponds to the end of the DNA !
	else:
		start = chromosomeLengthDNA-((proteinSequence.annotations["AAEnd"])*3+shift)
		end = chromosomeLengthDNA -1 - ((proteinSequence.annotations["AAStart"])*3+shift)
		#end += 1
	
	
	return (start, end)


def KZFPListToBed(KZFPList, targetFile, chromosomeDNAQuery):
	"""This function takes a list of krab zinc finger (in the SeqRecord format, with special annotations done as in the constructKZFPs function) and writes the main
	information in a bed file which name is targetFile.bed"""
	bedFile = open(targetFile, "w")
	chromosomeDNA = SeqIO.read(chromosomeDNAQuery, "fasta")
	chromosomeLengthDNA = len(chromosomeDNA.seq)
	for KZFP in KZFPList:
		#getting the shift corresponding to the reading frame to go back to DNA coordinates, as well as the strand (+ means forwads, - means reverse):
		shift, strand = getShiftStrand(KZFP)

		chromosome = getChromosome(KZFP)

		start, end = proteinToDNACoordinates(KZFP, chromosomeLengthDNA)
		#+1 because in a bed file, the end position is not part of the feature...
		end += 1
		nameList = [KZFP.id]+[KZFP.annotations["KRAB"]]+[KZFP.annotations["zfArrayContained"]]+KZFP.annotations["singleZFContained"]
		name = "-".join(nameList)
		#0: chromosome, 1: start, 2: end, 3: name, 4: score, 5: strand
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", str(strand)])+"\n"
		bedFile.write(bedLine)
	bedFile.close()

def singleZFListToBed(singleZFFullList, targetFile, chromosomeDNAQuery):
	bedFile = open(targetFile, "w")
	chromosomeDNA = SeqIO.read(chromosomeDNAQuery, "fasta")
	chromosomeLengthDNA = len(chromosomeDNA.seq)
	for sZF in singleZFFullList:
		shift, strand = getShiftStrand(sZF)

		chromosome = getChromosome(sZF)

		start, end = proteinToDNACoordinates(sZF, chromosomeLengthDNA)
		#+1 because in a bed file, the end position is not part of the feature...
		end += 1
		name = sZF.id
		#0: chromosome, 1: start, 2: end, 3: name, 4: score, 5: strand
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", str(strand)])+"\n"
		bedFile.write(bedLine)
	bedFile.close()

def zfArrayListToBed (zfArrayList, targetFile, chromosomeDNAQuery):
	bedFile = open(targetFile, "w")
	chromosomeDNA = SeqIO.read(chromosomeDNAQuery, "fasta")
	chromosomeLengthDNA = len(chromosomeDNA.seq)
	for zfArray in zfArrayList:
		shift, strand = getShiftStrand(zfArray)

		chromosome = getChromosome(zfArray)

		start, end = proteinToDNACoordinates(zfArray, chromosomeLengthDNA)
		#+1 because in a bed file, the end position is not part of the feature...
		end += 1
		#creating a list for the name, containing the id followed by all the single zinc fingers contained in the array
		nameList = [zfArray.id]+zfArray.annotations["singleZFContained"]
		name = "-".join(nameList)
		#0: chromosome, 1: start, 2: end, 3: name, 4: score, 5: strand
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", str(strand)])+"\n"
		bedFile.write(bedLine)
	bedFile.close()
def DNAzfArrayListToBed(zfArrayList, targetFile, chromosomeDNAQuery):
	bedFile = open(targetFile, "w")
	chromosomeDNA = SeqIO.read(chromosomeDNAQuery, "fasta")
	chromosomeLengthDNA = len(chromosomeDNA.seq)
	for zfArray in zfArrayList:
		#creating a list for the name, containing the id followed by all the single zinc fingers contained in the array
		chromosome = getChromosome(zfArray)
		start = (zfArray.annotations["DNAStart"])
		#+1 because in a bed file, the end position is not part of the feature...
		end = (zfArray.annotations["DNAEnd"])+1
		strand = zfArray.annotations["Strand"]
		nameList = [zfArray.id]+zfArray.annotations["singleZFContained"]
		name = "-".join(nameList)
		#0: chromosome, 1: start, 2: end, 3: name, 4: score, 5: strand
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", str(strand)])+"\n"
		bedFile.write(bedLine)
	bedFile.close()
def KRABToBed (KRABList, targetFile, chromosomeDNAQuery):
	bedFile = open(targetFile, "w")
	chromosomeDNA = SeqIO.read(chromosomeDNAQuery, "fasta")
	chromosomeLengthDNA = len(chromosomeDNA.seq)
	for KRAB in KRABList:
		#getting the shift corresponding to the reading frame to go back to DNA coordinates, as well as the strand (+ means forwads, - means reverse):
		shift, strand = getShiftStrand(KRAB)


		chromosome = getChromosome(KRAB)

		start, end = proteinToDNACoordinates(KRAB, chromosomeLengthDNA)
		#+1 because in a bed file, the end position is not part of the feature...
		end += 1
		name = KRAB.id
		#0: chromosome, 1: start, 2: end, 3: name, 4: score, 5: strand
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", str(strand)])+"\n"
		bedFile.write(bedLine)
	bedFile.close()


def protDomainListToDNADomainList (protDomainList, chromosomeDNAQuery):
	"""This function edits the annotations contained in the SeqRecord protein objects passed to as a list to adding
	DNA related information, i.e. DNAStart, DNAEnd, Shift and Strand"""
	DNADomainList = []
	chromosomeDNA = SeqIO.read(chromosomeDNAQuery, "fasta")
	chromosomeLengthDNA = len(chromosomeDNA.seq)
	for protDomain in protDomainList:

		protDomain.annotations["Shift"], protDomain.annotations["Strand"] = getShiftStrand(protDomain)
		#converting AA coordinates into DNA coordinates:
		
		protDomain.annotations["DNAStart"], protDomain.annotations["DNAEnd"] = proteinToDNACoordinates(protDomain, chromosomeLengthDNA)
		protDomain.annotations["DNASequence"] = chromosomeDNA.seq[protDomain.annotations["DNAStart"]:protDomain.annotations["DNAEnd"]+1]
		DNADomainList.append(protDomain)
	return DNADomainList

def addDNASequenceToKZFPDNAList(KZFPDNAList, chromosomeDNA):
	chromosome = SeqIO.read(chromosomeDNA, "fasta")
	for KZFPDNA in KZFPDNAList:
		KZFPDNA.annotations["DNASequence"] = chromosome.seq[KZFPDNA.annotations["DNAStart"]:KZFPDNA.annotations["DNAEnd"]+1]

def constructDNAKZFPs(KRABdomainList, zfArrayList, chromosomeQuery):
	"""This function constructs KRABs based on DNA information, thus allowing to build KRABs which are not 
	detected as in frame by the protein KZFP assembly approach. It requires lists of KRABs and zinc finger arrays
	which HAVE DNA information (annotations dic entries DNAStart, DNAEnd, Shift and Strand). It builds KZFP in a very
	naive fashion: each KRAB is assigned to the next zinc finger array containing at least 2 single zinc fingers downstream (for + strand) or upstream( for - strand).
	If no zinc finger array containing at least two single zinc fingers is found within a distance threshold, the KZFP constructs itself with the closest zinc finger array
	containing a single zinc finger. The KZFP sequence is built by reading the chromosome fasta DNA file between start and end positions"""
	chromosome = SeqIO.read(chromosomeQuery, "fasta")

	#building a list of the same size and same indexing as zfArrayList, but only containing the DNAStart value for the + strand (or DNAEnd for the - strand):
	FFzfDNAList = [array for array in zfArrayList if array.annotations["Strand"] == "+"]
	FFzfDNAList = sorted(FFzfDNAList, key=lambda array: array.annotations["DNAStart"])

	RFzfDNAList = [array for array in zfArrayList if array.annotations["Strand"] == "-"]
	RFzfDNAList = sorted(RFzfDNAList, key=lambda array: array.annotations["DNAStart"])

	FFzfDNAStartList =[array.annotations["DNAStart"] for array in FFzfDNAList]
	RFzfDNAStartList =[array.annotations["DNAEnd"] for array in RFzfDNAList]
	#print("Forward list: {}".format(FFzfDNAStartList))
	#print("Reverse list: {}".format(RFzfDNAStartList))
	#problem is those lists are not sorted !
	KZFPList = []
	KZFPcounter = 0
	
	for KRAB in KRABdomainList:

		#first possibility: the KRAB is positive strand:
		if KRAB.annotations["Strand"] == "+":
			#finds the next zinc finger after the end DNA base of the KRAB domain. Bisect returns the next index in the list which value is higher than a limit value.
			#If the limit value is above the highest value of the list, it returns an index = len(FFzfDNAStartList) (thus pointing out of the list). This situation can 
			#be handled with a simple if condition, which if not respected yields no krab zinc finger. The zinc finger array found with the bisect function contains at least one single zinc finger.
			#Thus, it is the zinc finger array to fall back to in case no zinc finger array with at least 2 single zinc finger is found in the distance threshold from the KRAB domain.
			nextArrayIndex = bisect(FFzfDNAStartList, KRAB.annotations["DNAEnd"])
			#i will allow indexation of the zinc finger arrays found whithin the distance threshold
			if nextArrayIndex != len(FFzfDNAList):
				i = 0
				iLimitReached = False
				multipleZFAFound = False
				while not iLimitReached:
					#print("Limit distance: {}".format(FFzfDNAList[nextArrayIndex+i].annotations["DNAEnd"] - KRAB.annotations["DNAStart"]))
					if nextArrayIndex+i != len(FFzfDNAList) and (FFzfDNAList[nextArrayIndex+i].annotations["DNAEnd"] - KRAB.annotations["DNAStart"]) < 40000:
						#print("limit ok")
						#checking if the zinc finger array has at least 2 single zinc fingers, otherwise looking for the next one. This is used to prevent KZFP building from stalling on zinc finger arrays containing single zinc fingers.
						if len(FFzfDNAList[nextArrayIndex+i].annotations["singleZFContained"])>1:
							#constructing the sequence of the KZFP: going from the first DNA base of the KRAB domain to the last DNA base of the zf array
							newSequence = chromosome.seq[KRAB.annotations["DNAStart"]:FFzfDNAList[nextArrayIndex+i].annotations["DNAEnd"]+1]
							KZFPList.append(SeqRecord(newSequence, id = "{}_KZFP_{}".format(chromosome.id, KZFPcounter), name = "{}_KZFP_{}".format(chromosome.id, KZFPcounter), description = "KRAB zinc finger (DNA assembly)",  \
							annotations = {"DNAStart": KRAB.annotations["DNAStart"], "DNAEnd": FFzfDNAList[nextArrayIndex+i].annotations["DNAEnd"], "Strand": "+", "KRAB": KRAB.id, "zfArrayContained": FFzfDNAList[nextArrayIndex+i].id,  \
							"zfStart":  FFzfDNAList[nextArrayIndex+i].annotations["DNAStart"], "singleZFContained": FFzfDNAList[nextArrayIndex+i].annotations["singleZFContained"], "DNASequence": newSequence, "chromosome": chromosome.id, "Genome": KRAB.annotations["Genome"]}))

							#print(KZFPList[KZFPcounter].id)
							#print("sequence length: {}".format(len(KZFPList[KZFPcounter].seq)))
							KZFPcounter += 1
							#no need to keep searching if one KZFP has been found
							iLimitReached = True
							multipleZFAFound = True

						else:
							#the previous zinc finger array candidate did not contain at least one zinc finger array: keep on searching
							i += 1
					else: 
						iLimitReached = True
				#if no zinc finger array containing more than 1 zinc finger was found: we append the closest one containing one zinc finger:
				if not multipleZFAFound:
					#control on distance is stil necessary: 
					if (FFzfDNAList[nextArrayIndex].annotations["DNAEnd"] - KRAB.annotations["DNAStart"]) < 40000:
						newSequence = chromosome.seq[KRAB.annotations["DNAStart"]:FFzfDNAList[nextArrayIndex].annotations["DNAEnd"]+1]
						KZFPList.append(SeqRecord(newSequence, id = "{}_KZFP_{}".format(chromosome.id, KZFPcounter), name = "{}_KZFP_{}".format(chromosome.id, KZFPcounter), description = "KRAB zinc finger (DNA assembly)",  \
						annotations = {"DNAStart": KRAB.annotations["DNAStart"], "DNAEnd": FFzfDNAList[nextArrayIndex].annotations["DNAEnd"], "Strand": "+", "KRAB": KRAB.id, "zfArrayContained": FFzfDNAList[nextArrayIndex].id,  \
						"zfStart":  FFzfDNAList[nextArrayIndex].annotations["DNAStart"], "singleZFContained": FFzfDNAList[nextArrayIndex].annotations["singleZFContained"], "DNASequence": newSequence, "chromosome": chromosome.id, "Genome": KRAB.annotations["Genome"]}))
						KZFPcounter += 1

		#second possibility: the KRAB is on the negative strand: 
		if KRAB.annotations["Strand"] == "-":
			print(KRAB.id)
			nextArrayIndex = bisect_left(RFzfDNAStartList, KRAB.annotations["DNAStart"])-1
			#bisect left returns the index at which the query element (here the DNAStart of the KRAB domain) should be inserted for the list to stay ordered. Therefore, if the Start of the KRAB domain is to the left of the 
			#first zinc finger array, it will still return the 0 index, which is not correct. This is why we substract one in the indexing. Also, care must be taken in order not to go out of the list
			if nextArrayIndex != -1:
				i=0
				iLimitReached = False
				multipleZFAFound = False
				while not iLimitReached:
					if nextArrayIndex-i != -1 and (KRAB.annotations["DNAEnd"] - RFzfDNAList[nextArrayIndex-i].annotations["DNAStart"]) < 40000:
						print(RFzfDNAList[nextArrayIndex-i].id)
						print(RFzfDNAList[nextArrayIndex-i].annotations["DNAEnd"])
						print(RFzfDNAStartList[nextArrayIndex-i])
						print(RFzfDNAList[nextArrayIndex-i].annotations["inFrameZFAContained"])
						print(len(RFzfDNAList[nextArrayIndex-i].annotations["singleZFContained"]))
						if len(RFzfDNAList[nextArrayIndex-i].annotations["singleZFContained"])>1:
							#constructing the sequence of the KZFP: going from the first DNA base of the zfArray to the last DNA base of the KRAB domain (because negative strand)
							newSequence = chromosome.seq[RFzfDNAList[nextArrayIndex-i].annotations["DNAStart"]:KRAB.annotations["DNAEnd"]+1]
							KZFPList.append(SeqRecord(newSequence, id = "{}_KZFP_{}".format(chromosome.id, KZFPcounter), name = "{}_KZFP_{}".format(chromosome.id, KZFPcounter), description = "KRAB zinc finger (DNA assembly)",  \
							annotations = {"DNAStart": RFzfDNAList[nextArrayIndex-i].annotations["DNAStart"], "DNAEnd": KRAB.annotations["DNAEnd"], "Strand": "-", "KRAB": KRAB.id, "zfArrayContained": RFzfDNAList[nextArrayIndex-i].id,  \
							"zfStart":  RFzfDNAList[nextArrayIndex-i].annotations["DNAEnd"], "singleZFContained": RFzfDNAList[nextArrayIndex-i].annotations["singleZFContained"], "DNASequence": newSequence, "chromosome": chromosome.id, "Genome": KRAB.annotations["Genome"]}))
							#print(KZFPList[KZFPcounter].id)
							#print("sequence length: {}".format(len(KZFPList[KZFPcounter].seq)))
							#print(nextArrayIndex-i)
							#print(KRAB.annotations["DNAEnd"] - RFzfDNAList[nextArrayIndex-i].annotations["DNAStart"])
							KZFPcounter += 1
							iLimitReached = True
							multipleZFAFound = True
						else: 
							i += 1
					else:
						iLimitReached = True
				if not multipleZFAFound:
					#control on distance is stil necessary: 
					if KRAB.annotations["DNAEnd"] - RFzfDNAList[nextArrayIndex].annotations["DNAStart"] < 40000:
						newSequence = chromosome.seq[RFzfDNAList[nextArrayIndex].annotations["DNAStart"]:KRAB.annotations["DNAEnd"]+1]
						KZFPList.append(SeqRecord(newSequence, id = "{}_KZFP_{}".format(chromosome.id, KZFPcounter), name = "{}_KZFP_{}".format(chromosome.id, KZFPcounter), description = "KRAB zinc finger (DNA assembly)",  \
						annotations = {"DNAStart": RFzfDNAList[nextArrayIndex].annotations["DNAStart"], "DNAEnd": KRAB.annotations["DNAEnd"], "Strand": "-", "KRAB": KRAB.id, "zfArrayContained": RFzfDNAList[nextArrayIndex].id,  \
						"zfStart":  RFzfDNAList[nextArrayIndex].annotations["DNAEnd"], "singleZFContained": RFzfDNAList[nextArrayIndex].annotations["singleZFContained"], "DNASequence": newSequence, "chromosome": chromosome.id, "Genome": KRAB.annotations["Genome"]}))
						KZFPcounter += 1
	return KZFPList

def mergeDNAKZFPLists(KZFPDNAListDNAZFA, KZFPDNAListProtZFA):
	"""This function intends to build a definitive KZFP list out of two methods of finding single ZFs.
	The definitive KZFPDNAList will contain all the KZFPs found by the DNA ZFA approach, as well as non overlapping KZFPs found by the protein ZFA approach (to be less restrictive than Michael).
	The function first generates bed files corresponding to both lists, and finds protein ZFA based KZFPs not overlapping with any DNA ZFA KZFP. It adds those to the DNA ZFA KZPFs to form the 
	definitive KZFP list"""
	#bed file for KZFPDNAListDNAZFA:
	bedFile = open("./Bed/Temporary/KZFPDNAListDNAZFA.bed", "w")
	for KZFP in KZFPDNAListDNAZFA:
		#finding the chromosome
		chromosome = getChromosome(KZFP)
		start = (KZFP.annotations["DNAStart"])
		#+1 because in a bed file, the end position is not part of the feature...
		end = (KZFP.annotations["DNAEnd"])+1
		name = KZFP.id
		strand = KZFP.annotations["Strand"]
		#0: chromosome, 1: start, 2: end, 3: name, 4: score, 5: strand
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", str(strand)])+"\n"
		bedFile.write(bedLine)
	bedFile.close()

	#bed file for KZFPDNAListProtZFA
	bedFile = open("./Bed/Temporary/KZFPDNAListProtZFA.bed", "w")
	for KZFP in KZFPDNAListProtZFA:
		#finding the chromosome
		chromosome = getChromosome(KZFP)
		start = (KZFP.annotations["DNAStart"])
		#+1 because in a bed file, the end position is not part of the feature...
		end = (KZFP.annotations["DNAEnd"])+1
		name = KZFP.id
		strand = KZFP.annotations["Strand"]
		#0: chromosome, 1: start, 2: end, 3: name, 4: score, 5: strand
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", str(strand)])+"\n"
		bedFile.write(bedLine)
	bedFile.close()

	#comparing the bed files:
	subprocess.call("bedtools intersect -v -a ./Bed/Temporary/KZFPDNAListProtZFA.bed -b ./Bed/Temporary/KZFPDNAListDNAZFA.bed > ./Bed/Temporary/NonOverlappingKZFPDNAProtList.bed", shell = True)

	#reading the IDs in the created bed file (it's the 4th column, i.e. column index 3 in python)
	bedFile = open("./Bed/Temporary/NonOverlappingKZFPDNAProtList.bed", 'r')
	#getting lines
	content = bedFile.readlines()
	# you may also want to remove whitespace characters like `\n` at the end of each line
	content = [x.strip() for x in content] 

	#we want to what's in each line before the first spacer or any kind (space or tab delimitation)
	spacer = "\t"
	domainSNPList = []
	singleZFListToAdd = []
	for field in content:
		ZFID = field.split(spacer)[3]
		singleZFListToAdd.append([singleZF for singleZF in KZFPDNAListProtZFA if ZFID == singleZF.id][0])
	bedFile.close()
	#renaming the IDs in the protein approach so that we can differentiate them from the ones bearing the same ID in the DNA ZFA approach.:
	for singleZF in singleZFListToAdd:
		singleZF.id+= "_protZFA"

	#merging the lists
	KZFPDNAList = KZFPDNAListDNAZFA+singleZFListToAdd
	print(singleZFListToAdd)
	return KZFPDNAList


def DNAKZFPListToBed (KZFPList, targetFile):
	bedFile = open(targetFile, "w")
	for KZFP in KZFPList:
		#finding the chromosome
		chromosome = getChromosome(KZFP)
		start = (KZFP.annotations["DNAStart"])
		#+1 because in a bed file, the end position is not part of the feature...
		end = (KZFP.annotations["DNAEnd"])+1
		nameList = [KZFP.id]+[KZFP.annotations["KRAB"]]+[KZFP.annotations["zfArrayContained"]]+KZFP.annotations["singleZFContained"]
		name = "-".join(nameList)
		strand = KZFP.annotations["Strand"]
		#0: chromosome, 1: start, 2: end, 3: name, 4: score, 5: strand
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", str(strand)])+"\n"
		print(bedLine)
		bedFile.write(bedLine)
	bedFile.close()
def KRABBboxListToBed(KRABBBoxListDNA, targetFile):
	bedFile = open(targetFile, "w")
	for KRABBbox in KRABBBoxListDNA:
		#finding the chromosome
		chromosome = getChromosome(KRABBbox)
		start = (KRABBbox.annotations["DNAStart"])
		#+1 because in a bed file, the end position is not part of the feature...
		end = (KRABBbox.annotations["DNAEnd"])+1
		name = KRABBbox.id
		strand = KRABBbox.annotations["Strand"]
		#0: chromosome, 1: start, 2: end, 3: name, 4: score, 5: strand
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", str(strand)])+"\n"
		print(bedLine)
		bedFile.write(bedLine)
	bedFile.close()

def KRABBDivListToBed(KRABBDivListDNA, targetFile):
	bedFile = open(targetFile, "w")
	for KRABBDiv in KRABBDivListDNA:
		#finding the chromosome
		chromosome = getChromosome(KRABBDiv)
		start = (KRABBDiv.annotations["DNAStart"])
		#+1 because in a bed file, the end position is not part of the feature...
		end = (KRABBDiv.annotations["DNAEnd"])+1
		name = KRABBDiv.id
		strand = KRABBDiv.annotations["Strand"]
		#0: chromosome, 1: start, 2: end, 3: name, 4: score, 5: strand
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", str(strand)])+"\n"
		print(bedLine)
		bedFile.write(bedLine)
	bedFile.close()

def DNAKZFPIDListToBed(KZFPList, targetFile):
	bedFile = open(targetFile, "w")
	for KZFP in KZFPList:
		#finding the chromosome
		chromosome = getChromosome(KZFP)
		start = (KZFP.annotations["DNAStart"])
		#+1 because in a bed file, the end position is not part of the feature...
		end = (KZFP.annotations["DNAEnd"])+1
		name = KZFP.id
		strand = KZFP.annotations["Strand"]
		#0: chromosome, 1: start, 2: end, 3: name, 4: score, 5: strand
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", str(strand)])+"\n"
		print(bedLine)
		bedFile.write(bedLine)
	bedFile.close()

def mergeDNAZFALists(DNABuiltZFArrays, protBuiltzfArrayListDNA):
	"""Intends to complement the list of DNA built zinc finger arrays with the zinc finger arrays found by the protein method, in order to be less restrictive
	in downstream analysis. The body of the function is basically the same than for the mergeDNAKZFPLists function above."""
	#bed file for DNABuiltZFArrays:
	bedFile = open("./Bed/Temporary/DNABuiltZFArrays.bed", "w")
	for ZFA in DNABuiltZFArrays:
		#finding the chromosome
		chromosome = getChromosome(ZFA)
		start = (ZFA.annotations["DNAStart"])
		#+1 because in a bed file, the end position is not part of the feature...
		end = (ZFA.annotations["DNAEnd"])+1
		name = ZFA.id
		strand = ZFA.annotations["Strand"]
		#0: chromosome, 1: start, 2: end, 3: name, 4: score, 5: strand
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", str(strand)])+"\n"
		bedFile.write(bedLine)
	bedFile.close()

	#bed file for protBuiltzfArrayListDNA
	bedFile = open("./Bed/Temporary/protBuiltzfArrayListDNA.bed", "w")
	for ZFA in protBuiltzfArrayListDNA:
		#finding the chromosome
		chromosome = getChromosome(ZFA)
		start = (ZFA.annotations["DNAStart"])
		#+1 because in a bed file, the end position is not part of the feature...
		end = (ZFA.annotations["DNAEnd"])+1
		name = ZFA.id
		strand = ZFA.annotations["Strand"]
		#0: chromosome, 1: start, 2: end, 3: name, 4: score, 5: strand
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", str(strand)])+"\n"
		bedFile.write(bedLine)
	bedFile.close()

	#comparing the bed files:
	subprocess.call("bedtools intersect -v -a ./Bed/Temporary/protBuiltzfArrayListDNA.bed -b ./Bed/Temporary/DNABuiltZFArrays.bed > ./Bed/Temporary/NonOverlappingZFAList.bed", shell = True)

	#reading the IDs in the created bed file (it's the 4th column, i.e. column index 3 in python)
	bedFile = open("./Bed/Temporary/NonOverlappingZFAList.bed", 'r')
	#getting lines
	content = bedFile.readlines()
	# you may also want to remove whitespace characters like `\n` at the end of each line
	content = [x.strip() for x in content] 

	#we want to what's in each line before the first spacer or any kind (space or tab delimitation)
	spacer = "\t"
	domainSNPList = []
	ZFAsToAdd = []
	for field in content:
		ZFID = field.split(spacer)[3]
		ZFAsToAdd.append([singleZF for singleZF in protBuiltzfArrayListDNA if ZFID == singleZF.id][0])
	bedFile.close()

	#merging the lists
	ZFAList = DNABuiltZFArrays+ZFAsToAdd
	print(ZFAsToAdd)
	return ZFAList




def getVCFinRegion(regionList, VCF, sortedBedFileName, VCFOutput):
	"""uses the tabix command line to extract VCF information from a big chromosome specific VCF file, restricted to the regions specified in the lisf of SeqRecord objects given"""
	#sorting regionList by DNAStart and then DNAEnd (required for tabix):

	regionList = sorted(regionList, key=lambda region: (region.annotations["DNAStart"], region.annotations["DNAEnd"]))

	bedFile = open(sortedBedFileName, "w")


	#storing the sorted file in a new bed:
	for region in regionList:
		chromosome = getChromosome(region)
		#bed format for VCF must not be "chr19" but only "19" in the first column...
		chromosomeNumber = chromosome[3:]
		start = (region.annotations["DNAStart"])
		#+1 because in a bed file, the end position is not part of the feature...
		end = (region.annotations["DNAEnd"])+1
		name = region.id
		strand = region.annotations["Strand"]
		#0: chromosome, 1: start, 2: end, 3: name, 4: score, 5: strand
		bedLine = "\t".join([chromosomeNumber, str(start), str(end), name, "0", str(strand)])+"\n"
		bedFile.write(bedLine)
	bedFile.close()
	#creating the VCF file
	subprocess.call("tabix -0 -R {} -h {} > {}".format(sortedBedFileName, VCF, VCFOutput), shell=True)
	
	#compressing for later tabix usage
	subprocess.call("bgzip -f {}".format(VCFOutput), shell=True)
	
	#sorting the compressed VCF file but VCF tools is not compatible with the extended VCF format from the publications...
	#subprocess.call("vcf-sort {}.gz > {}".format(VCFOutput, VCFOutput), shell=True)
	
	#compressing the sorted file
	#subprocess.call("bgzip -f {} > {}.gz".format(VCFOutput, VCFOutput), shell=True)
	
	#creating an index:
	subprocess.call("tabix -0f -p vcf {}.gz".format(VCFOutput), shell=True)
	

def extractLinkerDomain(KZFP, KRAB, zfArray):
	"""Given a KZFP SeqRecord object, which has corresponding zfArray and KRAB domains, constructs a SeqRecord object containing the 
	DNA sequence of the linker, as well as the DNAStart and DNAEnd positions stored in the .annotations dictionary"""
	#positive strand:
	#finding mutations in the space (thereafter caled "linker" domain) between the end of the KRAB and the beginning of the zfArray
	if KZFP.annotations["Strand"] == "+":
		linkerStart = KRAB.annotations["DNAEnd"]+1
		linkerEnd = zfArray.annotations["DNAStart"]-1
		newSequence = KZFP.annotations["DNASequence"][(linkerStart-KZFP.annotations["DNAStart"]):(linkerEnd-KZFP.annotations["DNAStart"])]



		linker = SeqRecord(newSequence, id = "{}_Linker".format(KZFP.id), name = "{}_Linker".format(KZFP.id), description = "Linker domain (DNA)",  \
				 		annotations = {"DNAStart": linkerStart, "DNAEnd": linkerEnd, "Strand": KZFP.annotations["Strand"], "BelongsToKZFP": KZFP.id, "chromosome": getChromosome(KZFP), "Genome": KRAB.annotations["Genome"]})
		KZFP.annotations["Linker"] = "{}_Linker".format(KZFP.id)

	#On the negative strand, the linker is found between the DNAEnd coordinate + 1 of the zinc finger array and the DNAStart coordinate of the KRAB domain
	elif KZFP.annotations["Strand"] == "-":
		linkerStart = zfArray.annotations["DNAEnd"]+1
		linkerEnd = KRAB.annotations["DNAStart"]-1
		newSequence = KZFP.annotations["DNASequence"][(linkerStart-KZFP.annotations["DNAStart"]):(linkerEnd-KZFP.annotations["DNAStart"])]



		linker = SeqRecord(newSequence, id = "{}_Linker".format(KZFP.id), name = "{}_Linker".format(KZFP.id), description = "Linker domain (DNA)",  \
				 		annotations = {"DNAStart": linkerStart, "DNAEnd": linkerEnd, "Strand": KZFP.annotations["Strand"], "BelongsToKZFP": KZFP.id, "chromosome": getChromosome(KZFP), "Genome": KRAB.annotations["Genome"]})
		KZFP.annotations["Linker"] = "{}_Linker".format(KZFP.id)

	return linker

def SNPToAAChange(domain, SNP):
	"""Takes a domain and a SNP, and returns a dictionary indicating the consequence of the SNP on the amino acid sequence of the domain ({"codonNumber": codonNumber, "altCodon": altCodon, "mutationType": mutationType})"""

	SNPMutationType = {}
	#DNASequence = Seq()
	#Checking if the strand is negative: the DNA sequence needs to be reversed, and the way the position of the codon is calculated is not the same
	if domain.annotations["Strand"] == "-":
		DNASequence = domain.annotations["DNASequence"].reverse_complement()
		codonNumber = (domain.annotations["DNAEnd"]-SNP[0])//3

		intraCodonPosition = (domain.annotations["DNAEnd"]-SNP[0])%3

		#Getting the reference codon
		refCodon = DNASequence[codonNumber*3:codonNumber*3+3]

		#translating the reference codon:
		refAA = refCodon.translate()

		#building the alternative codon
		refCodonList = list(str(refCodon))
		#if the strand is negative, the alternative DNA base must be changed to its complement		
		altBase = Seq(str(SNP[2]), generic_dna)

		altBase = altBase.complement()

		refCodonList[intraCodonPosition] = str(altBase)

		altCodon= Seq("".join(refCodonList), generic_dna)

		#translating the alternative codon
		altAA = altCodon.translate()

		#finding out which kind of mutation it is:
		mutationType = str()
		if str(altAA) == str(refAA):
			mutationType = "Synonymous"
		elif str(refAA) !="*":
			if str(altAA) == "*":
				mutationType = "Nonsense"
			else:
				mutationType = "Missense"
		elif str(refAA) == "*":
			mutationType = "Nonstop"
		#storing the kind of mutation linked to the SNP in a new dictionary 
		SNPMutationType = {"codonNumber": codonNumber, "refCodon": refCodon, "refAA": refAA, "altCodon": altCodon, "altAA": altAA, "mutationType": mutationType}
	#for a forward strand reding:
	else:
		DNASequence = domain.annotations["DNASequence"]

		codonNumber = (SNP[0]-domain.annotations["DNAStart"])//3

		intraCodonPosition = (SNP[0]-domain.annotations["DNAStart"])%3

		#Getting the reference codon
		refCodon = DNASequence[codonNumber*3:codonNumber*3+3]
		
		#translating the reference codon:
		refAA = refCodon.translate()
		#building the alternative codon
		refCodonList = list(str(refCodon))
		altBase = Seq(str(SNP[2]), generic_dna)

		refCodonList[intraCodonPosition] = str(altBase)

		altCodon= Seq("".join(refCodonList), generic_dna)
		#translating the alternative codon
		altAA = altCodon.translate()

		#finding out which kind of mutation it is:
		mutationType = str()
		if str(altAA) == str(refAA):
			mutationType = "Synonymous"
		elif str(refAA) !="*":
			if str(altAA) == "*":
				mutationType = "Nonsense"
			else:
				mutationType = "Missense"
		elif str(refAA) == "*":
			mutationType = "Nonstop"
		#storing the kind of mutation linked to the SNP in a new dictionary
		SNPMutationType = {"codonNumber": codonNumber, "refCodon": refCodon, "refAA": refAA, "altCodon": altCodon, "altAA": altAA, "mutationType": mutationType}
	return SNPMutationType

def findsingleZFFingerprint(singleZF):
	"""Takes a single zinc finger and returns its fingerprint and its negativeFingerprint as a string tuple (fingerprint, negativeFingerprint) (not a sequence or seqrecord object)"""
	#if str(singleZF.seq).count("H") >= 2 and str(singleZF.seq).count("C") >= 2:
	canonical = re.compile(".{0,2}C.{2,5}C.{12,15}H.{3,5}H")
	#building the same regex, but with a group encompassing the H used for fingerprint indexing
	canonicalEnd = re.compile(".{0,2}C.{2,5}C.{12,15}(H).{3,5}H")
	if re.match(canonical, str(singleZF.seq)) is not None and str(singleZF.seq).find("*") == -1: #question: does Michael allow stop codons in a fingerprint
		print(singleZF.seq)
		#finding the index corresponding to the first H of the regex
		HPenultimate = re.match(canonicalEnd, str(singleZF.seq)).start(1)
		fingerPrintCharList = []
		#fingerprints are amino acids at positions 6, 3, 2 and -1 given that the first H of the regex is index 7, and that index 0 amino acid does not exist (indexing jumps from -1 to 1)
		fingerPrintCharList.append(str(singleZF.seq)[HPenultimate-7])
		fingerPrintCharList.append(str(singleZF.seq)[HPenultimate-5])
		fingerPrintCharList.append(str(singleZF.seq)[HPenultimate-4])
		fingerPrintCharList.append(str(singleZF.seq)[HPenultimate-1])
		fingerPrint = "".join(fingerPrintCharList)

		#building the negativeFingerprint, which is the zinc finger w/o the base contacting residues, i.e. everything except the fingerprint.
		negativeFingerprint = list(str(singleZF.seq))
		Fingerprintindexes = [HPenultimate-1, HPenultimate-4, HPenultimate-5, HPenultimate-7]
		for index in sorted(Fingerprintindexes, reverse=True):
			del negativeFingerprint[index]
		negativeFingerprint = "".join(negativeFingerprint)
	else:
		fingerPrint = "XXXX"
		negativeFingerprint = (len(singleZF.seq)-4)*"X"
	return (fingerPrint, negativeFingerprint)


def findZFArrayFingerprint(zfArrayList, singleZFList):
	"""Takes a zfArray SeqRecord Object (single chromosome restricted) list and for each array: finds the fingerprint for each of the zinc fingers inside the array, and stores it in a list under .annotations["Fingerprint"]
	Conditions are: at least 2 C and 2 H otherwise non-canonical."""
	#for each array
	for array in zfArrayList:
		#find the single zfs corresponding to the ids stored in .annotations["SingleZFContained"]:
		arraySingleZFList = [singleZF for singleZF in singleZFList if singleZF.id in array.annotations["singleZFContained"]]
		#ordering the singleZFs by amino acid position, which solves the issue of strandedness (no need to write two different codes for the plus and minus strands)
		arraySingleZFList = sorted(arraySingleZFList, key = lambda singleZF: singleZF.annotations["AAStart"])
		array.annotations["Fingerprint"] = []
		array.annotations["NegativeFingerprint"] = []
		#for each single Zinc finger:
		for index, singleZF in enumerate(arraySingleZFList):

			fingerPrint, negativeFingerprint = findsingleZFFingerprint(singleZF)

			#appending the single zf fingerprint to the zfArray.annotations["Fingerprint"] entry
			if index == 0:
				array.annotations["Fingerprint"] = [fingerPrint]
				array.annotations["NegativeFingerprint"] = [negativeFingerprint]
			else:
				array.annotations["Fingerprint"].append(fingerPrint)
				array.annotations["NegativeFingerprint"].append(negativeFingerprint)

def findSNPList(domain, chromosomeNumber):
		domain.annotations["Mutated"] = False
		#creating a VCF file corresponding to the domain, i.e. a sub VCF file of the whole chromosome VCF with the range specified:
		#subprocess.call("tabix -0f -h ./VCF/denisova/VCF/hg19_1000g/T_hg19_1000g.{}.mod.vcf.gz {}:{}-{} > ./Mutations/VCF/{}.vcf".format(chromosomeNumber, chromosomeNumber, domain.annotations["DNAStart"]+1, domain.annotations["DNAEnd"], domain.id), shell=True)

		#strategy: we will store the SNPs in a list of tuple in the domain.annotations["SNPList"] dictionary entry
		domainSNPList = []
		vcf_reader = vcf.Reader(open("./Mutations/VCF/{}.vcf".format(domain.id), 'r'))

		for record in vcf_reader:
			#replacing the SNP in the sequence
			if record.is_snp and record.samples[0].data.GT != "./.":
				if record.samples[0].data.GQ >=30 and record.QUAL>=50 and record.FILTER != "LowQual":
					domainSNPList.append((record.POS-1, record.ALT[0]))
		#if the domain was mutated, we add it to the list of mutated domains
		if len(domainSNPList) > 0:
			domain.annotations["Mutated"] = True
			domain.annotations["SNPList"] = domainSNPList

def findsingleZFAltFingeprint(mutatedSingleZF):
	"""Takes a mutated singleZF (i.e. requires SNPMutated as True), builds the alternative fingerprint for each SNP of the single zinc finger. 
	If the alternative fingerprint is different from the original fingerprint, adds an entry to the ["SNPMutationType][SNP] dictionary"""
	if mutatedSingleZF.annotations["SNPMutated"]:
		for SNP in mutatedSingleZF.annotations["SNPMutationType"]:
			#replacing the AA mutation resulting from the SNP (as determined by the previously called SNPToAAChange function)

			altAASequenceList = list(mutatedSingleZF.seq)
			altAASequenceList[mutatedSingleZF.annotations["SNPMutationType"][SNP]["codonNumber"]] = str(mutatedSingleZF.annotations["SNPMutationType"][SNP]["altCodon"].translate())
			altAASequence = "".join(altAASequenceList)

			#building the alternative single zinc finger seqRecord object, which shares all attributes with the original zinc finger seqRecord object except for the amino acid sequence:
			altSingleZF = SeqRecord(Seq(altAASequence), id = "altSingleZF")
			altFingerprint, negativeAltFingerprint = findsingleZFFingerprint(altSingleZF)
			originalFingerprint, negativeFingerprint = findsingleZFFingerprint(mutatedSingleZF)
			changesFingerprint = False
			if originalFingerprint != altFingerprint:
				changesFingerprint = True
				mutatedSingleZF.annotations["SNPMutationType"][SNP]["altFingerprint"] = altFingerprint
			
			mutatedSingleZF.annotations["SNPMutationType"][SNP]["changesFingerprint"] = changesFingerprint

def findAltZFArrayFingerprint(mutatedZFArrayList, mutatedSingleZFList, chromosomeNumber):
	"""Takes a mutated zinc finger array list on which findZFArrayFingerprint has already been called (so that the reference fingerprint is already defined),
	 a mutated single zinc finger list on which findsingleZFAltFingeprint has already been called (so that "changesFingerprint" and "altFingerprint" of
	 	each single zinc finger already exists) finds the corresponding mutated 
		single zinc fingers and builds the alternative fingerprint"""
	#opening the file containing the corresponding single zinc fingers (by chromosome)
	singleZFList = []	
	with open("./Domains/{}_singleZFListDNA.dat".format(chromosomeNumber), "rb") as f:
		singleZFList = pickle.load(f)
	#for each array
	for array in mutatedZFArrayList:
		print(array.id)
		#find the mutated single zfs corresponding to the ids stored in .annotations["SingleZFContained"]. If the single zf is not found in the mutated
		#single zfs, it fetches the corresponding single zf in the nonmutated list.
		arrayMutatedSingleZFList = [singleZF for singleZF in mutatedSingleZFList if singleZF.id in array.annotations["singleZFContained"]]
		#completing the list with the single ZFs not found in the mutated single zf list:
		arrayMergedSingleZFList = [singleZF for singleZF in singleZFList if (singleZF.id in array.annotations["singleZFContained"] and singleZF.id not in [zf.id for zf in arrayMutatedSingleZFList])	] + arrayMutatedSingleZFList
		#ordering the singleZFs by amino acid position:
		arrayMergedSingleZFList = sorted(arrayMergedSingleZFList, key = lambda singleZF: singleZF.annotations["AAStart"])

		array.annotations["altFingerprint"] = array.annotations["Fingerprint"]
		#array.annotations["altNegativeFingerprint"] = array.annotations["NegativeFingerprint"]
		#going throuh the single zinc fingers in the zinc finger array:
		for index, singleZF in enumerate(arrayMergedSingleZFList):
			#if the zf is mutated, for each of its SNP, we check whether it changes the fingerprint
			if singleZF.id in [zf.id for zf in arrayMutatedSingleZFList]:
				#the SNPMutationType annotation is a dict, so we are running through keys here and must use them to explicitely get the information about the mutation
				for SNP in singleZF.annotations["SNPMutationType"]:
					if SNP["changesFingerprint"]:
						#note: if there are several non synonymous SNPs that change the fingerprint in the same out of frame zinc finger array, only the last one will
						#show in the alternative fingerprint. In practice, this is never the case.
						array.annotation["altFingerprint"][index] = SNP["altFingerprint"]

def extractMutationsFromVEPOutput(VEPOutputFile, MutationFolder, Indel = False):
	"""this function takes a VEP output file which first column ("Uploaded variation") is of format Chromosome_Position_Reference/alternative 
	and stores it in a new tab delimited text file with columns chromosome position reference alternative, one file per chromosome to make downstream analysis faster.
	The bool Indel is here to indicate the VEP output file taken as an input reports Indels, and not SNPs. By default: SNPs"""
	mutationType = "SNP"
	if Indel:
		mutationType = "Indel"

	VEP = open(VEPOutputFile, "r")
	#getting lines
	content = VEP.readlines()
	# you may also want to remove whitespace characters like `\n` at the end of each line
	content = [x.strip() for x in content] 
	
	#we want to what's in each line before the first spacer or any kind (space or tab delimitation)
	spacer = re.compile("\s+")
	#Line 0 must be jumped
	previousLine = str()
	previousChromosome = str()
	chromosomeCounter = 0
	for line in content[1:]:
		MutationEntry = spacer.split(line)[0]
		if MutationEntry != previousLine:
			chromosome, position, refAlt = MutationEntry.split("_")
			chromosome = "chr"+chromosome
			if chromosome != previousChromosome:
				if chromosomeCounter > 0:
					print("file closed")
					MutationOutput.close()
				MutationOutput = open("{}/{}_{}".format(MutationFolder, chromosome, mutationType), "w")
				chromosomeCounter += 1
			previousChromosome = chromosome
			reference, alternative = refAlt.split("/")
			
			MutationLine = "\t".join([chromosome, position, reference, alternative])+"\n"
			MutationOutput.write(MutationLine)
		previousLine = MutationEntry
	MutationOutput.close()
	VEP.close()


def findVEPDerivedSNPList(domain, chromosomeNumber, VEPDerivedFolder):
	"""Takes an output file from the function extractSNPFromVEPOutput and annotates the domain entered as an input if some SNPs are found within the domain"""
	domain.annotations["SNPMutated"] = False
	print(domain.id)
	print("Start {}, End {}".format(domain.annotations["DNAStart"], domain.annotations["DNAEnd"]))

	chromosomeSNPs = open("{}/{}_SNP".format(VEPDerivedFolder, chromosomeNumber), "r")

	#getting lines
	content = chromosomeSNPs.readlines()
	# you may also want to remove whitespace characters like `\n` at the end of each line
	content = [x.strip() for x in content] 

	#we want to what's in each line before the first spacer or any kind (space or tab delimitation)
	spacer = "\t"
	domainSNPList = []
	for SNP in content:
		position, ref, alt = SNP.split(spacer)[1:]
		position = int(position)
		position -= 1
		#checking if the SNP is in the limit of the DNA coordinates of the domain:
		if domain.annotations["DNAStart"] <= position < domain.annotations["DNAEnd"]:
			#since the VEP was derived from VCF used for ensembl, it's 1 based coordinate:
			domainSNPList.append((position, ref, alt))
			print("SNP found in domain {}".format(domain.id))
	#if the domain was mutated, we add it to the list of mutated domains
	if len(domainSNPList) > 0:
		domain.annotations["SNPMutated"] = True
		domain.annotations["SNPList"] = domainSNPList
	chromosomeSNPs.close()


def findVEPDerivedIndelList(domain, chromosomeNumber, VEPDerivedFolder):
	"""Takes an output file from the function extractSNPFromVEPOutput and annotates the domain entered as an input if some Indels are found within the domain"""
	domain.annotations["IndelMutated"] = False

	chromosomeIndels = open("{}/{}_Indel".format(VEPDerivedFolder, chromosomeNumber), "r")

	#getting lines
	content = chromosomeIndels.readlines()
	# you may also want to remove whitespace characters like `\n` at the end of each line
	content = [x.strip() for x in content] 

	#we want to what's in each line before the first spacer or any kind (space or tab delimitation)
	spacer = "\t"
	domainIndelList = []
	for Indel in content:
		position, ref, alt = Indel.split(spacer)[1:]
		position = int(position)
		position -= 1
		#checking if the SNP is in the limit of the DNA coordinates of the domain:
		if domain.annotations["DNAStart"] <= position < domain.annotations["DNAEnd"]:
			#since the VEP was derived from VCF used for ensembl, it's 1 based coordinate:
			domainIndelList.append((position, ref, alt))

	#if the domain was mutated, we add it to the list of mutated domains
	if len(domainIndelList) > 0:
		print(domain.id)
		domain.annotations["IndelMutated"] = True
		domain.annotations["IndelList"] = domainIndelList
	chromosomeIndels.close()

def buildVCFFromVEP(MutationFolder, VEPOutputFile, Indel = False):
	"""Takes a VEPOutputFile and reconstructs a basic VCF file for IGV viewing purposes
	Required files: #CHROM  POS ID  REF ALT QUAL    FILTER  INFO"""
	mutationType = "SNP"
	if Indel:
		mutationType = "Indel"
	VEP = open("{}/{}".format(MutationFolder, VEPOutputFile), "r")
	#getting lines
	content = VEP.readlines()
	# you may also want to remove whitespace characters like `\n` at the end of each line
	content = [x.strip() for x in content] 
	
	#we want to what's in each line before the first spacer or any kind (space or tab delimitation)
	spacer = re.compile("\s+")
	#Line 0 must be jumped
	previousLine = str()
	previousChromosome = str()
	chromosomeCounter = 0
	for line in content[1:]:
		MutationEntry = spacer.split(line)[0]
		if MutationEntry != previousLine:
			chromosome, position, refAlt = MutationEntry.split("_")
			chromosome = "chr"+chromosome
			if chromosome != previousChromosome:
				if chromosomeCounter > 0:
					print("file closed")
					MutationOutput.close()
				MutationOutput = open("{}/{}_{}.vcf".format(MutationFolder, chromosome, mutationType), "w")
				MutationOutput.write("##fileformat=VCFv4.1\n")
				MutationOutput.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])+"\n")
				chromosomeCounter += 1
			previousChromosome = chromosome
			reference, alternative = refAlt.split("/")
			
			MutationLine = "\t".join([chromosome, position, "id", reference, alternative, ".", "PASS", "info"])+"\n"
			MutationOutput.write(MutationLine)
		previousLine = MutationEntry
	MutationOutput.close()
	VEP.close()
def getKRABStart(KZFP, KRABList):
	"""Takes a KZFP and a list of KRAB domains (must contain the KRAB present in the KZFP) and returns the DNA start coordinate of the KRAB domain corresponding to the KZFP"""
	KRABid = KZFP.annotations["KRAB"]
	if KRABid is not None:
		KRAB = [KRAB for KRAB in KRABList if KRABid == KRAB.id][0]
		return KRAB.annotations["DNAStart"]
	else:
		return None

def getKRABEnd(KZFP, KRABList):
	KRABid = KZFP.annotations["KRAB"]
	if KRABid is not None:
		KRAB = [KRAB for KRAB in KRABList if KRABid == KRAB.id][0]
		return KRAB.annotations["DNAEnd"]
	else:
		return None



def getKRABLength(KZFP, KRABList):
	"""Takes a KZFP and a list of KRAB domains (which must contain the KRAB belonging to the KZFP) and returns the length of the KRAB domain 
	in DNA coordinates"""
	KRABid = KZFP.annotations["KRAB"]
	if KRABid is not None:
		KRAB = [KRAB for KRAB in KRABList if KRABid == KRAB.id][0]
		return KRAB.annotations["DNAEnd"] - KRAB.annotations["DNAStart"]
	else:
		return None

def getKRABSequence(KZFP, KRABList):
	"""Takes a KZFP and a list of KRAB domains (which must contain the KRAB belonging to the KZFP) and returns amino acid sequence of the KRAB"""
	KRABid = KZFP.annotations["KRAB"]
	if KRABid is not None:
		KRAB = [KRAB for KRAB in KRABList if KRABid == KRAB.id][0]
		return str(KRAB.seq)
	else:
		return None
def getLinkerSequence(KZFP, LinkerList):
	Linkerid = KZFP.annotations["Linker"]
	if Linkerid is not None:
		linker = [linker for linker in LinkerList if Linkerid == linker.id][0]
		return str(linker.seq)
	else:
		return None
def getLinkerStart(KZFP, LinkerList):
	Linkerid = KZFP.annotations["Linker"]
	if Linkerid is not None:
		linker = [linker for linker in LinkerList if Linkerid == linker.id][0]
		return linker.annotations["DNAStart"]
	else:
		return None

def getLinkerEnd(KZFP, LinkerList):
	Linkerid = KZFP.annotations["Linker"]
	if Linkerid is not None:
		linker = [linker for linker in LinkerList if Linkerid == linker.id][0]
		return linker.annotations["DNAEnd"]
	else:
		return None


def getLinkerLength(KZFP, LinkerList):
	"""Takes a KZFP and a list of Linker domains (which must contain the Linker belonging to the KZFP) and returns the length of the Linker domain 
	in DNA coordinates"""
	linkerid = KZFP.annotations["Linker"]
	if linkerid is not None:
		Linker = [Linker for Linker in LinkerList if linkerid == Linker.id][0]
		return Linker.annotations["DNAEnd"] - Linker.annotations["DNAStart"]
	else:
		return None

def getOutOfFrameZFArrayLength(KZFP, outOfFrameZFArrayListDNA):
	"""Takes a KZFP and a list of OutOfFrameZFArray domains (which must contain the OutOfFrameZFArray belonging to the KZFP) and returns the length of the OutOfFrameZFArray domain 
	in DNA coordinates"""
	OutOfFrameZFArray = [OutOfFrameZFArray for OutOfFrameZFArray in outOfFrameZFArrayListDNA if KZFP.annotations["zfArrayContained"] == OutOfFrameZFArray.id][0]
	return OutOfFrameZFArray.annotations["DNAEnd"] - OutOfFrameZFArray.annotations["DNAStart"]

def getZincFingers(KZFP, outOfFrameZFArrayListDNA, singleZFDNAList):
	"""Takes a KZFP and a list of single zinc fingers and returns a list of strings containing all zinc fingers ordered form N to C terminus"""
	ZincFingers = []
	for sZFid in getSingleZFContained(KZFP, outOfFrameZFArrayListDNA):
		sZF = [sFZ for sFZ in singleZFDNAList if sZFid == sFZ.id][0]
		ZincFingers.append(str(sZF.seq))
	return ZincFingers

def getKZFPFingerprint(KZFP, outOfFrameZFArrayListDNA):
	"""Takes a KZFP containing an outOfFrameZFArray and returns the fingerprint as a list of characters string"""
	return [outOfFrameZFArray.annotations["Fingerprint"] for outOfFrameZFArray in outOfFrameZFArrayListDNA if outOfFrameZFArray.id == KZFP.annotations["zfArrayContained"]][0]

def getKZFPNegativeFingerprint(KZFP, outOfFrameZFArrayListDNA):
	"""Takes a KZFP containing an outOfFrameZFArray and returns the negative fingerprint as a list of characters string"""
	return [outOfFrameZFArray.annotations["NegativeFingerprint"] for outOfFrameZFArray in outOfFrameZFArrayListDNA if outOfFrameZFArray.id == KZFP.annotations["zfArrayContained"]][0]

def getInFrameZFArrayNumber(KZFP, OutOfFrameZFArrayListDNA):
	OutOfFrameZFArray = [OutOfFrameZFArray for OutOfFrameZFArray in OutOfFrameZFArrayListDNA if KZFP.annotations["zfArrayContained"] == OutOfFrameZFArray.id][0]
	return len(OutOfFrameZFArray.annotations["inFrameZFAContained"])

def getInFrameZFAContained(KZFP, OutOfFrameZFArrayListDNA):
	OutOfFrameZFArray = [OutOfFrameZFArray for OutOfFrameZFArray in OutOfFrameZFArrayListDNA if KZFP.annotations["zfArrayContained"] == OutOfFrameZFArray.id][0]
	return (OutOfFrameZFArray.annotations["inFrameZFAContained"])

def getOutOfFrameZFArrayStart(KZFP, OutOfFrameZFArrayListDNA):
	OutOfFrameZFArray = [OutOfFrameZFArray for OutOfFrameZFArray in OutOfFrameZFArrayListDNA if KZFP.annotations["zfArrayContained"] == OutOfFrameZFArray.id][0]
	return(OutOfFrameZFArray.annotations["DNAStart"])

def getOutOfFrameZFArrayEnd(KZFP, OutOfFrameZFArrayListDNA):
	OutOfFrameZFArray = [OutOfFrameZFArray for OutOfFrameZFArray in OutOfFrameZFArrayListDNA if KZFP.annotations["zfArrayContained"] == OutOfFrameZFArray.id][0]
	return(OutOfFrameZFArray.annotations["DNAEnd"])

def getSingleZFContained(KZFP, OutOfFrameZFArrayListDNA):
	OutOfFrameZFArray = [OutOfFrameZFArray for OutOfFrameZFArray in OutOfFrameZFArrayListDNA if KZFP.annotations["zfArrayContained"] == OutOfFrameZFArray.id][0]
	return (OutOfFrameZFArray.annotations["singleZFContained"])

def getKRABBbox(KZFP, KRABBboxDNAList):
	"""retuns of tuple: (id, start, end, length) only if the KZFP contains a KRABBbox"""
	if KZFP.annotations["ContainsKRABBbox"]:
		KRABBbox = [KRABBbox for KRABBbox in KRABBboxDNAList if KZFP.annotations["KRABBboxContained"] == KRABBbox.id][0]
		return(KRABBbox.id, KRABBbox.annotations["DNAStart"], KRABBbox.annotations["DNAEnd"], KRABBbox.annotations["DNAEnd"]-KRABBbox.annotations["DNAStart"])
	else:
		return(None, None, None, None)
def getKRABBDiv(KZFP, KRABBDivDNAList):
	"""retuns of tuple: (id, start, end, length) only if the KZFP contains a KRABBbox"""
	if KZFP.annotations["ContainsKRABBDiv"]:
		KRABBDiv = [KRABBDiv for KRABBDiv in KRABBDivDNAList if KZFP.annotations["KRABBDivContained"] == KRABBDiv.id][0]
		return(KRABBDiv.id, KRABBDiv.annotations["DNAStart"], KRABBDiv.annotations["DNAEnd"], KRABBDiv.annotations["DNAEnd"]-KRABBDiv.annotations["DNAStart"])
	else:
		return(None, None, None, None)


def KZFPidZFAToBed(targetFile):
	"""Takes the whole KZFP dataframe and builds a bed file for oofZFAs, identified by the ID of the KZFP"""
	KZFPTable = []
	with open("KZFPTable.dat", "rb") as f:
		KZFPTable = pickle.load(f)

	bedFile = open(targetFile, "w")
	for index, KZFP in KZFPTable.iterrows():
		#finding the chromosome
		chromosome = KZFP["chromosome"]
		start = (KZFP["OutOfFrameZFAStart"])
		#+1 because in a bed file, the end position is not part of the feature...
		end = (KZFP["OutOfFrameZFAEnd"])+1
		name = KZFP["id"]
		strand = KZFP["Strand"]
		#0: chromosome, 1: start, 2: end, 3: name, 4: score, 5: strand
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", str(strand)])+"\n"
		print(bedLine)
		bedFile.write(bedLine)
	bedFile.close()


def KZFPToBed(targetFile):
	"""Takes the whole KZFP dataframe and builds a bed file for it, start and end are the whole KZFP length, identified by the ID of the KZFP"""
	KZFPTable = []
	with open("KZFPTable.dat", "rb") as f:
		KZFPTable = pickle.load(f)

	bedFile = open(targetFile, "w")
	for index, KZFP in KZFPTable.iterrows():
		#finding the chromosome
		chromosome = KZFP["chromosome"]
		start = (KZFP["DNAStart"])
		#+1 because in a bed file, the end position is not part of the feature...
		end = (KZFP["DNAEnd"])+1
		name = KZFP["id"]
		strand = KZFP["Strand"]
		#0: chromosome, 1: start, 2: end, 3: name, 4: score, 5: strand
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", str(strand)])+"\n"
		print(bedLine)
		bedFile.write(bedLine)
	bedFile.close()

def mergeBedFiles(params):
	"""merges bedfiles across all chromosomes for the given domain"""
	chromosomeList, domain = params
	for index, N in enumerate(chromosomeList):
		if index == 0: 
			subprocess.call("cat ./Bed/{}_{}.bed > ./Bed/{}.bed".format(N, domain, domain), shell=True)
		else:
			subprocess.call("cat ./Bed/{}_{}.bed >> ./Bed/{}.bed".format(N, domain, domain), shell=True)
		subprocess.call("rm ./Bed/{}_{}.bed".format(N, domain), shell=True)

def mergeDomains(params):
	"""merges pickled files across all chromosomes for the given domain"""
	chromosomeList, domain = params
	domainList = []

	for index, N in enumerate(chromosomeList):
		with open("./Domains/{}_{}.dat".format(N, domain), "rb") as f:
			domainList += pickle.load(f)
		subprocess.call("rm ./Domains/{}_{}.dat".format(N, domain, domain), shell=True)


	#saving merged dna domain list
	with open("./Domains/{}.dat".format(domain), "wb") as f:
		pickle.dump(domainList,f)
def getGeneFromAnnotation(KZFP, Annotation):
	"""Takes a KZFP seqrecord object and an annotation pandas dataframe (as constructed in the plotting.py script) and returns the 
	GeneName and the geneID in a tuple if a match is found, otherwise returns noname"""
	if KZFP.id in list(Annotation["id"]):
		return (Annotation.ix[KZFP.id]["GeneName"], Annotation.ix[KZFP.id]["GeneID"], Annotation.ix[KZFP.id]["GeneType"])
	else:
		#return("noGene_{}".format(KZFP.id), "noGene_{}".format(KZFP.id), "noGene_{}".format(KZFP.id), "noGene_{}".format(KZFP.id))
		return("notAnnotated", "notAnnotated", "notAnnotated")

def generateSNPInDelBed(variantFolder, targetFile):
	"""Takes a SNP or InDel containing file generated by the extractMutationsFromVEPOutput function and creates a BED6 file for both InDels and SNPs"""
	variantList = ["SNP", "Indel"]
	#remark: no Y chromosome
	chromosomeList = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]
	
	for variantForm in variantList:
		BedFile = open(targetFile+variantForm+".bed", 'w')
		for chromosomeNumber in chromosomeList:
			#opening the variant containing file:
			variantFile = open("{}/{}_{}".format(variantFolder, chromosomeNumber, variantForm), "r")
			#getting lines
			content = variantFile.readlines()
			# you may also want to remove whitespace characters like `\n` at the end of each line
			content = [x.strip() for x in content] 
	
			#we want to what's in each line before the first spacer or any kind (space or tab delimitation)
			spacer = re.compile("\s+")

			for line in content:
				chromosome, start, ref, alt = spacer.split(line)
				#note: no score and no strand, as variants are supposed to be homozygous
				bedLine = "\t".join([chromosome, start, str(int(start)+1), "{}/{}".format(ref, alt)])+"\n"
				BedFile.write(bedLine)
			variantFile.close()
		BedFile.close()
	
def generateKZFPKRABBed(KZFPTable, targetFile):
	"""Generates a BED6 file for all KRABs in KZFPs"""
	BedFile = open(targetFile, 'w')
	for index, KZFP in KZFPTable.iterrows():
		chromosome = KZFP["chromosome"]
		start = KZFP["KRABStart"]
		end = KZFP["KRABEnd"]
		name = KZFP["KRABid"]
		strand = KZFP["Strand"]
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", strand])+"\n"
		BedFile.write(bedLine)
	BedFile.close()

def generateKZFPLinkerBed(KZFPTable, targetFile):
	BedFile = open(targetFile, 'w')
	for index, KZFP in KZFPTable.iterrows():
		chromosome = KZFP["chromosome"]
		start = KZFP["SpacerStart"]
		end = KZFP["SpacerEnd"]
		name = KZFP["Spacerid"]
		strand = KZFP["Strand"]
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", strand])+"\n"
		BedFile.write(bedLine)
	BedFile.close()	

def generateKZFPoofZFABed(KZFPTable, targetFile):
	BedFile = open(targetFile, 'w')
	for index, KZFP in KZFPTable.iterrows():
		chromosome = KZFP["chromosome"]
		start = KZFP["OutOfFrameZFAStart"]
		end = KZFP["OutOfFrameZFAEnd"]
		name = KZFP["OutOfFrameZFAid"]
		strand = KZFP["Strand"]
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", strand])+"\n"
		BedFile.write(bedLine)
	BedFile.close()

def generateKZFPGeneIDoofZFABed(KZFPTable, targetFile):
	BedFile = open(targetFile, 'w')
	for index, KZFP in KZFPTable.iterrows():
		chromosome = KZFP["chromosome"]
		start = KZFP["OutOfFrameZFAStart"]
		end = KZFP["OutOfFrameZFAEnd"]
		name = KZFP["GeneID"]
		strand = KZFP["Strand"]
		additional = KZFP["GeneName"]
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", strand, additional])+"\n"
		BedFile.write(bedLine)
	BedFile.close()

def generateKZFPGeneIDKRABBed(KZFPTable, targetFile):
	BedFile = open(targetFile, 'w')
	for index, KZFP in KZFPTable.iterrows():
		chromosome = KZFP["chromosome"]
		start = KZFP["KRABStart"]
		end = KZFP["KRABEnd"]
		name = KZFP["GeneID"]
		strand = KZFP["Strand"]
		additional = KZFP["GeneName"]
		bedLine = "\t".join([chromosome, str(start), str(end), name, "0", strand, additional])+"\n"
		BedFile.write(bedLine)
	BedFile.close()
def getBedLength(bedFile):
	"""Takes the path to a bedfile as an argument and returns the sum of the lengths of the intervals"""
	totalLength = 0
	with open(bedFile, "r") as f:
			content = f.readlines()
			# you may also want to remove whitespace characters like `\n` at the end of each line
			content = [x.strip() for x in content] 
	
			#we want to what's in each line before the first spacer or any kind (space or tab delimitation)
			spacer = re.compile("\s+")

			for line in content:
				#avoiding commented lines:
				if line[0] != "#":
					#Start and stop are columns indexes 1 and 2 (0 based):
					totalLength += int(spacer.split(line)[2]) - int(spacer.split(line)[1])
	return totalLength

def file_len(fname):
	"""returns the line count of any file"""
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1
def checkForStopCodon(domain):
	"""Takes any domain with an AA amino sequence as a seqrecord .seq attribute and checks for stop codon, that is "*" character.
	Returns True if one is found, false otherwise"""
	if "*" in str(domain.seq):
		return True
	else:
		return False

def separateChromosomes(fullGenomeFile, chromosomeListFile):
	"""Genomes are assembled in a full file. This function takes a list of chromosomenames as a text file, reads it and 
	writes only those listed chromosomes in seperated files."""
	chromosomeNamesList = getChromosomeNamesList(chromosomeListFile)
	
	#limit number of files in which we will be running the analysis, only for tests
	#betaLimit = 4
	genome = open(fullGenomeFile, "r")
	#getting lines
	content = genome.readlines()
	# you may also want to remove whitespace characters like `\n` at the end of each line
	content = [x.strip() for x in content] 
	
	spacer = re.compile("\s+")
	numberOfFiles = 0
	for lineIndex, line in enumerate(content):
		
		#if numberOfFiles > betaLimit:
		#	break
		
		#that's the new fasta sequence indicator. The scaffold name starts right after the > and ends with a space (i.e. spacer)
		if line[0] == ">":
			#recuperating the name of the fasta sequence
			scaffoldName = spacer.split(line)[0]
			#removing the ">" from the name of the fasta sequence: 
			scaffoldName = scaffoldName[1:]
			#we check that the sequence name is only in the list
			if scaffoldName in chromosomeNamesList:
				print(scaffoldName)

				#opening a new fasta file for the fasta sequence. We are using these as "chromosomes" were used to seperate the analysis
				#in the human genome
				newFasta =  open("./Chromosomes/{}.fa".format(scaffoldName), "w")
				#writing the chromosome name
				newFasta.write(line+"\n")
				#reading the fasta file until the next > to catch the DNA sequence
				for localIndex, localLine in enumerate(content[lineIndex+1:]):
					if localLine[0] == ">":
						break
					newFasta.write(localLine+"\n")
				newFasta.close()
				numberOfFiles = numberOfFiles+1
	genome.close()
	
	return chromosomeNamesList

def getChromosomeNamesList(chromosomeListFile):
	"""reads the txt file given and returns the chromosome list corresponding"""
	chromosomeNamesList = []
	#reading the list:
	with open(chromosomeListFile, "r") as f:
		content = f.readlines()
		# you may also want to remove whitespace characters like `\n` at the end of each line
		content = [x.strip() for x in content] 
		#we want to what's in each line before the first spacer or any kind (space or tab delimitation)
		spacer = re.compile("\s+")
		for line in content:
			chromosomeNamesList.append(spacer.split(line)[0])
	return chromosomeNamesList

def addKRABBboxToKZFP(KZFPDNAList, KRABBboxListDNA, KRABDNAList, outOfFrameZFArrayListDNA):
	"""
	for each KZFP in the list, if a KRAB b box is in the same strand, fully contained in the spacer (that is between the end of the KRAB domain
		and the beginning of the oZFA and if the distance between the KRAB and the KRAB b box is less than 2000 bp, the KRAB b box is added).
	Should be done chromosome by chromosome
	"""
	#KRABAtoBMaxDistance = 1500

	#negative strand
	for KZFP in [KZFP for KZFP in KZFPDNAList if KZFP.annotations["Strand"] =="-"]:
		ContainsKRABBbox = False
		for KRABBbox in [KRABBbox for KRABBbox in KRABBboxListDNA if KRABBbox.annotations["Strand"] == "-"]:
			#same strand is already checked before -> only limits of the domain and the right order need to be respected:
			if (getKRABStart(KZFP, KRABDNAList) > KRABBbox.annotations["DNAEnd"]) and (KRABBbox.annotations["DNAStart"] > getOutOfFrameZFArrayEnd(KZFP, outOfFrameZFArrayListDNA)):	
				KZFP.annotations["KRABBboxContained"] = KRABBbox.id
				ContainsKRABBbox = True
		KZFP.annotations["ContainsKRABBbox"] = ContainsKRABBbox
	#positive strand
	for KZFP in [KZFP for KZFP in KZFPDNAList if KZFP.annotations["Strand"] =="+"]:
		ContainsKRABBbox = False
		for KRABBbox in [KRABBbox for KRABBbox in KRABBboxListDNA if KRABBbox.annotations["Strand"] == "+"]:
			if (KRABBbox.annotations["DNAStart"] > getKRABEnd(KZFP, KRABDNAList)) and (getOutOfFrameZFArrayStart(KZFP, outOfFrameZFArrayListDNA) > KRABBbox.annotations["DNAEnd"]):
				KZFP.annotations["KRABBboxContained"] = KRABBbox.id
				ContainsKRABBbox = True
		KZFP.annotations["ContainsKRABBbox"] = ContainsKRABBbox
def addKRABBDivToKZFP(KZFPDNAList, KRABBDivListDNA, KRABDNAList, outOfFrameZFArrayListDNA):
	"""
	for each KZFP in the list, if a KRAB b divergent box is in the same strand, fully contained in the spacer (that is between the end of the KRAB domain
		and the beginning of the oZFA and if the distance between the KRAB and the KRAB b box is less than 2000 bp, the KRAB b box is added).
	Should be done chromosome by chromosome
	"""
	#KRABAtoBMaxDistance = 1500

	#negative strand
	for KZFP in [KZFP for KZFP in KZFPDNAList if KZFP.annotations["Strand"] =="-"]:
		ContainsKRABBDiv = False
		for KRABBDiv in [KRABBDiv for KRABBDiv in KRABBDivListDNA if KRABBDiv.annotations["Strand"] == "-"]:
			#same strand is already checked before -> only limits of the domain and the right order need to be respected:
			if (getKRABStart(KZFP, KRABDNAList) > KRABBDiv.annotations["DNAEnd"]) and (KRABBDiv.annotations["DNAStart"] > getOutOfFrameZFArrayEnd(KZFP, outOfFrameZFArrayListDNA)):	
				KZFP.annotations["KRABBDivContained"] = KRABBDiv.id
				ContainsKRABBDiv = True
		KZFP.annotations["ContainsKRABBDiv"] = ContainsKRABBDiv
	#positive strand
	for KZFP in [KZFP for KZFP in KZFPDNAList if KZFP.annotations["Strand"] =="+"]:
		ContainsKRABBDiv = False
		for KRABBDiv in [KRABBDiv for KRABBDiv in KRABBDivListDNA if KRABBDiv.annotations["Strand"] == "+"]:
			if (KRABBDiv.annotations["DNAStart"] > getKRABEnd(KZFP, KRABDNAList)) and (getOutOfFrameZFArrayStart(KZFP, outOfFrameZFArrayListDNA) > KRABBDiv.annotations["DNAEnd"]):
				KZFP.annotations["KRABBDivContained"] = KRABBDiv.id
				ContainsKRABBDiv = True
		KZFP.annotations["ContainsKRABBDiv"] = ContainsKRABBDiv

def mergeVCFs(MutationFolder):
	#note: no Y chromosome as the sequenced skeletons were females
	chromosomeList = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]
	IndelFileNames = ["{}_Indel.vcf".format(N) for N in chromosomeList]
	SNPFileNames = ["{}_SNP.vcf".format(N) for N in chromosomeList]

	with open(MutationFolder+"/AllIndels.vcf", 'w') as outfile:
	    for fname in IndelFileNames:
	        with open(MutationFolder+"/"+fname) as infile:
	            for line in infile:
	                outfile.write(line)
	
	with open(MutationFolder+"/AllSNPs.vcf", 'w') as outfile:
	    for fname in SNPFileNames:
	        with open(MutationFolder+"/"+fname) as infile:
	            for line in infile:
	                outfile.write(line)

def createBedsFromVEP(MutationFolder):
	#note: no Y chromosome as the sequenced skeletons were females
	chromosomeList = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]
	IndelFileNames = ["{}_Indel".format(N) for N in chromosomeList]
	SNPFileNames = ["{}_SNP".format(N) for N in chromosomeList]
	spacer = re.compile("\s+")
	with open(MutationFolder+"/AllIndels.bed", 'w') as outfile:
		for fname in IndelFileNames:
			with open(MutationFolder+"/"+fname) as infile:
				for line in infile:
					print (spacer.split(line))
					chromosome, DNAStart, ref, alt, space = spacer.split(line)
					outfile.write("\t".join([chromosome, str(DNAStart), str(int(DNAStart)+1), ref+"/"+alt])+"\n")
	
	with open(MutationFolder+"/AllSNPs.bed", 'w') as outfile:
		for fname in SNPFileNames:
			with open(MutationFolder+"/"+fname) as infile:
				for line in infile:
					chromosome, DNAStart, ref, alt, space = spacer.split(line)
					outfile.write("\t".join([chromosome, str(DNAStart), str(int(DNAStart)+1), ref+"/"+alt])+"\n")
#only for fingerprints, however we will just expand it to everything and produce one big table for KRAB homology, fingerprint homology and negative fingerprint homology.
def getHomologousFingerprints(KZFPTable1path, KZFPTable2path, canonicalThreshold, identityThreshold, NegativeFingerprint = False):
	"""Takes two pathes to KZFPTables and a threshold between 0 and 1 and returns, for each KZFP in table 1, the corresponding KZFPs in table 2.
	If table1 == table2: the same KZFP id is avoided for each KZFP. If the tables are different, the tables are assumed to belong to 
	different species and thus no result is filtered based on KZFP id"""
	if NegativeFingerprint:
		toCompare = "NegativeFingerprint"
	else:
		toCompare = "Fingerprint"
	#we need to modify the matrix and add the scores corresponding to stop codons which are not included in the blosum80 default matrix of biopython
	scoreMatrix = matlist.blosum80
	#listing the unique letters of the tuples:
	allResidues = [residue for residuePair in scoreMatrix.keys() for residue in residuePair]
	#removing duplicates
	uniqueResidues = list(set(allResidues))
	#we need to add the score for the stop codons ("*"): they have a -6 alignment score with every residue except themselves: -1.
	#for the fingerprints, we also force an alignment by quadruplets by introducing a character ("_") which marks the beginning and 
	#end of each single zinc finger fingerprint. We assign a score called sameFrameScore to matching "_" and a negative one for all "_" missmatches
	#we also penalize gap opening 
	#we need to construct all possible tuples of combinations:
	sameFrameScore = 100
	for residue in uniqueResidues:
		newKey = (residue, "*")
		scoreMatrix[newKey]= -6
		newKey2 = (residue, "_")
		scoreMatrix[newKey2]=-1*sameFrameScore
	#adding the final touch: the stop vs himself: it has a gap of -1:
	scoreMatrix[("*", "*")]= -1

	#we also add the fictive character "_" to the blosum matrix, and we will assign a high score when it matches itself and a null score when it is matched with any other residue
	scoreMatrix[("_", "_")] = sameFrameScore
	KZFPTable1 = []
	KZFPTable2 = []
	with open(KZFPTable1path, "rb") as f:
		KZFPTable1 = pickle.load(f)
	with open(KZFPTable2path, "rb") as f:
		KZFPTable2 = pickle.load(f)
	homologyTable = pd.DataFrame()

	for index1, KZFP1 in KZFPTable1.iterrows():

		print(KZFP1["id"])
		if KZFP1["CanonicalScore"] > canonicalThreshold:
			sequence1List = [sequence for sequence in KZFP1["{}".format(toCompare)] if "X" not in sequence]
			sequence1 = "_".join(sequence1List)
			
			for index2, KZFP2 in KZFPTable2.iterrows():
				sequence2 = str()
				avoidSameID = False
				if KZFPTable1path == KZFPTable2path:
					if KZFP1["id"] == KZFP2["id"]:
						avoidSameID = True
				if KZFP2["CanonicalScore"] > canonicalThreshold and not avoidSameID: #and not KZFP1["OutOfFrameZFAid"] == KZFP2["OutOfFrameZFAid"]:
					sequence2List = [sequence for sequence in KZFP2["{}".format(toCompare)] if "X" not in sequence]
					sequence2 = "_".join(sequence2List)
					for a in pairwise2.align.globalds(sequence1, sequence2, scoreMatrix, -20, 0, one_alignment_only=True, penalize_end_gaps = (False, False)):
						#calcul du pourcentage de similitude:
						zippedAlignment = zip(a[0], a[1])
						alignmentLength = a[4]
						score = a[2]
						identity = 0
						toSubstract = 0
						for alignPos in zippedAlignment:
							#score correction for "_":
							if alignPos == ("_", "_"):
								score -= 100
							#we correct the score for everything except a match between the "_" and a gap
							elif "_" in alignPos and not "-" in alignPos:
								score += 100
	
							#percentage identity: we check brutally character by character but avoid counting matchint "_" in the identity:
							if alignPos[0] == alignPos[1]:
								if alignPos[0] == "_":
									toSubstract += 1
								else:
									identity += 1
						identity = identity/(alignmentLength-toSubstract)
						score = score/alignmentLength
						if identity >= identityThreshold:
							homologyTable = homologyTable.append(other = pd.DataFrame([KZFP1["id"], KZFP1["GeneName"], KZFP1["{}".format(toCompare)], KZFP1["CanonicalScore"], KZFP1["GeneType"], KZFP1["KRABSequence"], KZFP2["id"], KZFP2["GeneName"], KZFP2["{}".format(toCompare)], KZFP2["CanonicalScore"], KZFP2["GeneType"], KZFP2["KRABSequence"], identity, score]).T)
	homologyTable.columns = ["KZFP1", "GeneName1", "{}1".format(toCompare), "CanonicalScore1", "GeneType1", "KRABSequence1", "KZFP2", "GeneName2", "{}2".format(toCompare), "CanonicalScore2", "GeneType2", "KRABSequence2", "identity", "score"]
	homologyTable = homologyTable.reset_index(drop = True)
	return homologyTable

def getHomology(KZFPTable1path, KZFPTable2path, canonicalThreshold, withLinkers, nProcesses):
	#checking if the tables are the same:
	sameSpecies = False
	if KZFPTable1path == KZFPTable2path:
		sameSpecies = True
	#we need to modify the matrix and add the scores corresponding to stop codons which are not included in the blosum80 default matrix of biopython
	scoreMatrix = matlist.blosum80
	#listing the unique letters of the tuples:
	allResidues = [residue for residuePair in scoreMatrix.keys() for residue in residuePair]
	#removing duplicates
	uniqueResidues = list(set(allResidues))
	#we need to add the score for the stop codons ("*"): they have a -6 alignment score with every residue except themselves: -1.
	#for the fingerprints, we also force an alignment by quadruplets by introducing a character ("_") which marks the beginning and 
	#end of each single zinc finger fingerprint. We assign a score called sameFrameScore to matching "_" and a negative one for all "_" missmatches
	#we also penalize gap opening 
	#we need to construct all possible tuples of combinations:
	sameFrameScore = 100
	for residue in uniqueResidues:
		newKey = (residue, "*")
		scoreMatrix[newKey]= -6
		newKey2 = (residue, "_")
		scoreMatrix[newKey2]=-1*sameFrameScore
	#adding the final touch: the stop vs himself: it has a gap of -1:
	scoreMatrix[("*", "*")]= -1
	#and the stop codon vs the special character:
	scoreMatrix[("_", "*")] = -1*sameFrameScore

	#we also add the fictive character "_" to the blosum matrix, and we will assign a high score when it matches itself and a null score when it is matched with any other residue
	scoreMatrix[("_", "_")] = sameFrameScore
	KZFPTable1 = []
	KZFPTable2 = []
	with open(KZFPTable1path, "rb") as f:
		KZFPTable1 = pickle.load(f)
	with open(KZFPTable2path, "rb") as f:
		KZFPTable2 = pickle.load(f)
	homologyTable = pd.DataFrame()


	#parameters must be zipped, thus there will be as many parameters as rows in KZFPTable1, and the other elements will remain the same.
	param1 = []
	for index1, KZFP1 in KZFPTable1.iterrows():
			param1.append(KZFP1)
	nParam = len(param1)
	param2 = [KZFPTable2]*nParam
	param3 = [sameSpecies]*nParam
	param4 = [scoreMatrix]*nParam
	param5 = [canonicalThreshold]*nParam
	param6 = [withLinkers]*nParam
	params = zip(param1, param2, param3, param4, param5, param6)
	if nProcesses != 0:
		pool = multiprocessing.Pool(nProcesses)
	else:
		pool = multiprocessing.Pool()

	#results is a list of pd dataframes. It needs to be concatenated
	results = pool.map(compareKZFPs, params)

	for dataframe in results:
		homologyTable = homologyTable.append(dataframe)


	homologyTable = homologyTable.reset_index(drop = True)
	return homologyTable


def compareKZFPs(params):
	KZFP1, KZFPTable2, sameSpecies, scoreMatrix, canonicalThreshold, withLinkers = params
	print(KZFP1["id"])
	homologyTable = pd.DataFrame()

	if KZFP1["CanonicalScore"] > canonicalThreshold:
		ZincFingers1List = [sequence for sequence in KZFP1["ZincFingers"]]
		ZincFingers1 = "_".join(ZincFingers1List)

		Fingerprint1List = [sequence for sequence in KZFP1["Fingerprint"] if "X" not in sequence]
		Fingerprint1 = "_".join(Fingerprint1List)

		NegativeFingerprint1 = [sequence for sequence in KZFP1["NegativeFingerprint"] if "X" not in sequence]
		NegativeFingerprint1 = "_".join(NegativeFingerprint1)

		for index2, KZFP2 in KZFPTable2.iterrows():
			Fingerprint2 = str()
			NegativeFingerprint2 = str()
			avoidSameID = False
			if sameSpecies:
				if KZFP1["id"] == KZFP2["id"]:
					avoidSameID = True
			if KZFP2["CanonicalScore"] > canonicalThreshold and not avoidSameID:
				ZincFingers2List = [sequence for sequence in KZFP2["ZincFingers"]]
				ZincFingers2 = "_".join(ZincFingers2List)

				Fingerprint2List = [sequence for sequence in KZFP2["Fingerprint"] if "X" not in sequence]
				Fingerprint2 = "_".join(Fingerprint2List)

				NegativeFingerprint2List = [sequence for sequence in KZFP2["NegativeFingerprint"] if "X" not in sequence]
				NegativeFingerprint2 = "_".join(NegativeFingerprint2List)

				ZincFingerAlignment = zip("", "")
				ZincFingerScore = 0
				ZincFingerIdentity = 0

				FingerprintAlignment = zip("", "")
				FingerprintScore = 0
				FingerprintIdentity = 0

				NegativeFingerprintScore = 0
				NegativeFingerprintIdentity = 0

				KRABAlignment = zip("", "")
				KRABScore = 0
				KRABIdentity = 0

				LinkerAlignment = zip("", "")
				LinkerScore = 0
				LinkerIdentity = 0

				alignmentLength = 0
				toSubstract = 0

				##Full length zinc fingers: (does not care whether zinc fingers are degenerated or not, but aligns analogous to the Fingerprint (see below and
					#methods in my report))
				#note, although the for loop, there is only one alignment (note one_alignment_only = True)
				for a in pairwise2.align.globalds(ZincFingers1, ZincFingers2, scoreMatrix, -20, 0, one_alignment_only=True, penalize_end_gaps = (False, False)):
					#calcul du pourcentage de similitude:
					zippedAlignment = zip(a[0], a[1])
					ZincFingerAlignment = format_alignment(*a)
					alignmentLength = a[4]
					ZincFingerScore = a[2]
					for alignPos in zippedAlignment:
						#score correction for "_":
						if alignPos == ("_", "_"):
							ZincFingerScore -= 100
						#we correct the score for everything except a match between the "_" and a gap
						elif "_" in alignPos and not "-" in alignPos:
							ZincFingerScore += 100

						#percentage identity: we check brutally character by character but avoid counting matchint "_" in the identity:
						if alignPos[0] == alignPos[1]:
							if alignPos[0] == "_":
								toSubstract += 1
							else:
								ZincFingerIdentity += 1
				ZincFingerIdentity = ZincFingerIdentity/(alignmentLength-toSubstract)
				ZincFingerScore = ZincFingerScore/(alignmentLength-toSubstract)



				#fingerprint
				#note, although the for loop, there is only one alignment (note one_alignment_only = True)
				for a in pairwise2.align.globalds(Fingerprint1, Fingerprint2, scoreMatrix, -20, 0, one_alignment_only=True, penalize_end_gaps = (False, False)):
					#calcul du pourcentage de similitude:
					zippedAlignment = zip(a[0], a[1])
					FingerprintAlignment = format_alignment(*a)
					toSubstract = 0
					alignmentLength = a[4]
					FingerprintScore = a[2]
					for alignPos in zippedAlignment:
						#score correction for "_":
						if alignPos == ("_", "_"):
							FingerprintScore -= 100
						#we correct the score for everything except a match between the "_" and a gap
						elif "_" in alignPos and not "-" in alignPos:
							FingerprintScore += 100

						#percentage identity: we check brutally character by character but avoid counting matchint "_" in the identity:
						if alignPos[0] == alignPos[1]:
							if alignPos[0] == "_":
								toSubstract += 1
							else:
								FingerprintIdentity += 1
				FingerprintIdentity = FingerprintIdentity/(alignmentLength-toSubstract)
				FingerprintScore = FingerprintScore/(alignmentLength-toSubstract)

				#negative fingerprint
				for a in pairwise2.align.globalds(NegativeFingerprint1, NegativeFingerprint2, scoreMatrix, -20, 0, one_alignment_only=True, penalize_end_gaps = (False, False)):
					#calcul du pourcentage de similitude:
					zippedAlignment = zip(a[0], a[1])
					alignmentLength = a[4]
					NegativeFingerprintScore = a[2]
					toSubstract = 0
					for alignPos in zippedAlignment:
						#score correction for "_":
						if alignPos == ("_", "_"):
							NegativeFingerprintScore -= 100
						#we correct the score for everything except a match between the "_" and a gap
						elif "_" in alignPos and not "-" in alignPos:
							NegativeFingerprintScore += 100

						#percentage identity: we check brutally character by character but avoid counting matchint "_" in the identity:
						if alignPos[0] == alignPos[1]:
							if alignPos[0] == "_":
								toSubstract += 1
							else:
								NegativeFingerprintIdentity += 1
				NegativeFingerprintIdentity = NegativeFingerprintIdentity/(alignmentLength-toSubstract)
				NegativeFingerprintScore = NegativeFingerprintScore/(alignmentLength-toSubstract)

				#KRAB domain: note that the settings of the pairwise alignment are more classic, and that the "_" character doesnt need to be corrected for
				#as it is simply not present
				if KZFP1["KRABSequence"] is not None and KZFP2 ["KRABSequence"] is not None:
					for a in pairwise2.align.globalds(KZFP1["KRABSequence"], KZFP2["KRABSequence"], scoreMatrix, -1, -1, one_alignment_only=True):
						#calcul du pourcentage de similitude:
						zippedAlignment = zip(a[0], a[1])
						KRABAlignment = format_alignment(*a)
						alignmentLength = a[4]
						KRABScore = a[2]
						for alignPos in zippedAlignment:
							#percentage identity: we check brutally character by character but avoid counting matchint "_" in the identity:
							if alignPos[0] == alignPos[1]:

								KRABIdentity += 1					

					KRABScore = KRABScore/(alignmentLength)
					KRABIdentity = KRABIdentity/(alignmentLength)
				else:
					KRABScore = None
					KRABIdentity = None
					KRABAlignment = None

				if withLinkers is True: 
					#Linker: only a DNA sequence (reminder: linker = spacer, NOT the linker which separates individual zinc fingers)
					#note that out of frame zinc finger arrays do not have a linker.
					if KZFP1["SpacerSequence"] is not None and KZFP2["SpacerSequence"] is not None:
						for a in pairwise2.align.globalxs(KZFP1["SpacerSequence"], KZFP2["SpacerSequence"], -1, -1, one_alignment_only=True):
							#calcul du pourcentage de similitude:
							LinkerAlignment = format_alignment(*a)
							zippedAlignment = zip(a[0], a[1])
							alignmentLength = a[4]
							LinkerScore = a[2]

						#in that case, matching DNA bases = 1, not matching = 0, therefore score and identity are the same
						LinkerScore = LinkerScore/(alignmentLength)
					else: 
						LinkerScore = None
						LinkerAlignment = None


				homologyTable = homologyTable.append(other = pd.DataFrame([KZFP1["Genome"], KZFP1["id"], KZFP1["GeneName"], KZFP1["GeneID"], KZFP1["KRABSequence"], KZFP1["ZincFingers"], KZFP1["Fingerprint"], KZFP1["NegativeFingerprint"], KZFP1["CanonicalScore"], KZFP1["GeneType"], \
					KZFP2["Genome"], KZFP2["id"], KZFP2["GeneName"], KZFP2["GeneID"], KZFP2["KRABSequence"], KZFP2["ZincFingers"], KZFP2["Fingerprint"], KZFP2["NegativeFingerprint"], KZFP2["CanonicalScore"], KZFP2["GeneType"], \
					KRABAlignment, KRABScore, KRABIdentity, ZincFingerAlignment, ZincFingerScore, ZincFingerIdentity, FingerprintAlignment, FingerprintScore, FingerprintIdentity, NegativeFingerprintScore, NegativeFingerprintIdentity]).T)
				if withLinkers is True:
					homologyTable = homologyTable.append(other = pd.DataFrame([LinkerAlignment, LinkerScore]).T)

		homologyTable.columns = ["Genome1", "KZFP1", "GeneName1", "GeneID1", "KRABSequence1", "ZincFingers1", "Fingerprint1", "NegativeFingerprint1", "CanonicalScore1", "GeneType1", \
		"Genome2", "KZFP2", "GeneName2", "GeneID2", "KRABSequence2", "ZincFingers2", "Fingerprint2", "NegativeFingerprint2", "CanonicalScore2", "GeneType2", \
		"KRABAlignment", "KRABScore", "KRABIdentity","ZincFingerAlignment", "ZincFingerScore", "ZincFingerIdentity", "FingerprintAlignment", "FingerprintScore", "FingerprintIdentity", "NegativeFingerprintScore", \
		"NegativeFingerprintIdentity"]

		if withLinkers is True:
			homologyTable.columns.append("LinkerAlignment", "LinkerIdentity")

		homologyTable = homologyTable.reset_index(drop = True)
	return homologyTable