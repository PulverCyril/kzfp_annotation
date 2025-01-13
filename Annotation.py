import numpy as np
import pickle
import pandas as pd
import subprocess
import math
from Functions import findZFArrayFingerprint
from Functions import getChromosome
from Functions import getKRABLength
from Functions import getLinkerLength
from Functions import getOutOfFrameZFArrayLength
from Functions import getKZFPFingerprint
from Functions import getInFrameZFArrayNumber
from Functions import extractLinkerDomain
from Functions import getOutOfFrameZFArrayStart
from Functions import getOutOfFrameZFArrayEnd
from Functions import DNAKZFPIDListToBed
from Functions import getGeneFromAnnotation
from Functions import getOutOfFrameZFArrayEnd
from Functions import getKRABStart
from Functions import getKRABEnd
from Functions import getLinkerStart
from Functions import getLinkerEnd
from Functions import generateSNPInDelBed
from Functions import generateKZFPKRABBed
from Functions import generateKZFPLinkerBed
from Functions import generateKZFPoofZFABed
from pybedtools import BedTool
from Functions import getBedLength
from Functions import file_len
from Functions import generateKZFPGeneIDoofZFABed
from Functions import generateKZFPGeneIDKRABBed
from Functions import getInFrameZFAContained
from Functions import getSingleZFContained
from Functions import checkForStopCodon
from Functions import findsingleZFFingerprint
from Functions import getKRABBbox
from Functions import getKRABBDiv
from Functions import getKZFPNegativeFingerprint
from Functions import getChromosomeNamesList
from Functions import getKRABSequence
from Functions import getZincFingers
from Functions import getLinkerSequence

#Specify the chromosomes on which you want to run the annotation:
chromosomeList = getChromosomeNamesList("/home/bopekno/Documents/Genomes/hg19test/chromosomeNamesList")

#Set to False if you have already created an EnsemblGenes.txt file which is formatted to UCSC standards.
#When set to true, tries to do a naive conversion from Ensembl to UCSC annotation by: adding "chr" in front 
#of all chromosome names, transforming 1 and -1 strands into + and -, filling empty gene names and filling spaces in gene names 
#(otherwise parsing fails)
naiveEnsemblToUCSCFormatting = True

#set to True if you would like to run downstream analysis (such as fingerprint alignments) with out of frame zinc finger arrays not contained in any KZFPs
includeOutOfFrameZFAs = True


singleZFDNAList = []
inFrameZFArrayListDNA = []
outOfFrameZFArrayListDNA = []
#only used if includeOutOfFrameZFAs is True
outOfFrameZFArrayToAnalyze = []
KRABDNAList = []
KRABBboxDNAList = []
KRABBDivDNAList = []
KZFPDNAList = []

#loading domain files:
with open("./Domains/singleZFListDNA.dat", "rb") as f:
        singleZFDNAList=pickle.load(f)
with open("./Domains/inFrameZFArrayListDNA.dat", "rb") as f:
        inFrameZFArrayListDNA=pickle.load(f)
with open("./Domains/outOfFrameZFArrayList.dat", "rb") as f:
        outOfFrameZFArrayListDNA=pickle.load(f)
with open("./Domains/KRABListDNA.dat", "rb") as f:
        KRABDNAList=pickle.load(f)
with open("./Domains/KRABBboxDNAList.dat", "rb") as f:
        KRABBboxDNAList=pickle.load(f)
with open("./Domains/KRABBDivDNAList.dat", "rb") as f:
        KRABBDivDNAList=pickle.load(f)
with open("./Domains/KZFPDNAList.dat", "rb") as f:
        KZFPDNAList=pickle.load(f)

#building the spacer domains:
linkerDNAList = [] 
for KZFP in KZFPDNAList:
    KRAB = [KRAB for KRAB in KRABDNAList if KZFP.annotations["KRAB"] == KRAB.id][0]
    outOfFrameZFArray = [OutOfFrameZFArray for OutOfFrameZFArray in outOfFrameZFArrayListDNA if KZFP.annotations["zfArrayContained"] == OutOfFrameZFArray.id][0]
    linkerDNAList.append(extractLinkerDomain(KZFP, KRAB, outOfFrameZFArray))

#we select the subgroup of out of frame zinc finger arrays which are not already contained in KZFPs
if includeOutOfFrameZFAs is True:
    #setting a threshold on a minimal number of single zinc finger to consider an array a true out of frame ZFA array filters most of them. This parameter may
    #be set to any positive integer including 0.
    minOFZFALength = 2
    outOfFrameZFArrayToAnalyze = [OFZFA for OFZFA in outOfFrameZFArrayListDNA if OFZFA.id not in [KZFP.annotations["zfArrayContained"] for KZFP in KZFPDNAList] \
    and len(OFZFA.annotations["singleZFContained"]) >= minOFZFALength]

    #quickfix to make OFZFA compatible with the annotation process for KZFPs:
    for OFZFA in outOfFrameZFArrayToAnalyze:
        OFZFA.annotations["KRAB"] = None
        OFZFA.annotations["Linker"] = None
        OFZFA.annotations["zfArrayContained"] = OFZFA.id
        OFZFA.annotations["ContainsKRABBbox"] = False
        OFZFA.annotations["ContainsKRABBDiv"] = False

    KZFPDNAList.extend(outOfFrameZFArrayToAnalyze)
    #TODO: merge selected out of frame ZFAs with KZFPs in a way that 1) allows to use the same annotation pipeline and 2) still allows to build the KZFPTable
    #Modify the KZFP compare function to allow for comparisons between out of frame zinc finger arrays and KZFPs, add a comparison between full-length KZFPs, 
    #add a comparison between spacer domains.


#Ensembl annotation:
#building a bed file for all KZFPs:
DNAKZFPIDListToBed(KZFPDNAList, "./Bed/KZFPDNAList.bed")
if naiveEnsemblToUCSCFormatting:
    #formatting the ensembl gene output to match bed specifications:
    EnsemblGenesDF = pd.read_csv("EnsemblGenesBadFormat.txt", sep = "\t")
    #missing gene names
    EnsemblGenesDF["Gene name"] = EnsemblGenesDF["Gene name"].fillna("no_Ensembl_name")
    #gene names with spaces
    EnsemblGenesDF["Gene name"] = EnsemblGenesDF["Gene name"].str.replace("\s+", "_")
    EnsemblGenesDF.to_csv("EnsemblGenes.txt", sep="\t", index = False, index_label = False)

with open("EnsemblGenes.txt", "r") as ensemblOutput:
    content = ensemblOutput.readlines()
    # you may also want to remove whitespace characters like `\n` at the end of each line
    content = [x.strip() for x in content] 

    #we want to what's in each line before the first spacer or any kind (space or tab delimitation) to modify it to chr->content
    spacer = "\t"

    with open("EnsemblGenes.bed", "w") as bedFile:
        bedFile.write("\t".join(["#chromosome", "start", "end", "GeneName", "score", "strand", "GeneID", "GeneType"])+"\n")
        for field in content[1:]:
            chromosome, start, end, GeneName, score, strand, GeneID, GeneType = field.split(spacer)
            if naiveEnsemblToUCSCFormatting:
                #modifiying the chromosome and the strand so that it fits bedfile specifications
                chromosome = "chr"+ chromosome
                realStrand = "+"
                if strand == str(-1):
                    realStrand = "-"
            bedLine = "\t".join([chromosome, str(start), str(end), GeneName, "0", realStrand, GeneID, GeneType])+"\n"
            bedFile.write(bedLine)

#overlapping the ensembl annotation with the KZFPs
subprocess.call("""bedtools intersect -s -a ./Bed/KZFPDNAList.bed -b EnsemblGenes.bed -wo -f 0.5 | awk -v OFS="\t" '{print $4,$10,$13,$15,$8,$9,$14}' > KZFPAnnotation.tsv""", shell = True)
#subprocess.call("""bedtools intersect -s -a ./Bed/KZFPDNAList.bed -b EnsemblProteinCoding.bed -wo -f 0.5 | awk -v OFS="\t" '{print $4,$10,$13,$14,$8,$9}' > KZFPAnnotation.tsv""", shell = True)

#Creating a dataframe for the annotation file:
KZFPAnnotation = pd.read_table("KZFPAnnotation.tsv", sep = "\t", engine = "python", index_col = 0, names = np.array(["GeneName", "GeneID", "Overlap", "EnsemblGeneStart", "EnsemblGeneEnd", "GeneType"]))
#adding the id in a column to use drop_duplicates on it:
KZFPAnnotation["id"] = KZFPAnnotation.index
#adding a column corresponding to the ratio between the overlap and the length of the ensembl gene
KZFPAnnotation["OverlapEnsemblGeneRatio"] = KZFPAnnotation["Overlap"]/(KZFPAnnotation["EnsemblGeneEnd"]- KZFPAnnotation["EnsemblGeneStart"])
#sorting by index (or id) and then overlap, from lowest to highest:
KZFPAnnotation=KZFPAnnotation.sort_values(by = ["id", "OverlapEnsemblGeneRatio"], ascending = [True, False])
KZFPAnnotation.to_csv(path_or_buf = "./Tables/PreTrimAnnotation.csv", sep = "\t")
#dropping duplicates
KZFPAnnotation = KZFPAnnotation.drop_duplicates(subset = ["id"], keep = "first")



KZFPTable = pd.DataFrame.from_records([[KZFP.annotations["Genome"] for KZFP in KZFPDNAList], \
    [getChromosome(KZFP) for KZFP in KZFPDNAList], \
    [KZFP.annotations["DNAStart"] for KZFP in KZFPDNAList], \
    [KZFP.annotations["DNAEnd"] for KZFP in KZFPDNAList], \
    [getGeneFromAnnotation(KZFP, KZFPAnnotation)[0] for KZFP in KZFPDNAList], \
    [KZFP.annotations["Strand"] for KZFP in KZFPDNAList], \
    [KZFP.id for KZFP in KZFPDNAList], \
    [getGeneFromAnnotation(KZFP, KZFPAnnotation)[1] for KZFP in KZFPDNAList], \
    [getGeneFromAnnotation(KZFP, KZFPAnnotation)[2] for KZFP in KZFPDNAList], \
    [KZFP.annotations["KRAB"] for KZFP in KZFPDNAList], \
    [getKRABStart(KZFP, KRABDNAList) for KZFP in KZFPDNAList], \
    [getKRABEnd(KZFP, KRABDNAList) for KZFP in KZFPDNAList], \
    [getKRABLength(KZFP, KRABDNAList) for KZFP in KZFPDNAList], \
    [getKRABSequence(KZFP, KRABDNAList) for KZFP in KZFPDNAList], \
    [KZFP.annotations["ContainsKRABBbox"] for KZFP in KZFPDNAList], \
    [getKRABBbox(KZFP, KRABBboxDNAList)[0] for KZFP in KZFPDNAList], \
    [getKRABBbox(KZFP, KRABBboxDNAList)[1] for KZFP in KZFPDNAList], \
    [getKRABBbox(KZFP, KRABBboxDNAList)[2] for KZFP in KZFPDNAList], \
    [getKRABBbox(KZFP, KRABBboxDNAList)[3] for KZFP in KZFPDNAList], \
    [KZFP.annotations["ContainsKRABBDiv"] for KZFP in KZFPDNAList], \
    [getKRABBDiv(KZFP, KRABBDivDNAList)[0] for KZFP in KZFPDNAList], \
    [getKRABBDiv(KZFP, KRABBDivDNAList)[1] for KZFP in KZFPDNAList], \
    [getKRABBDiv(KZFP, KRABBDivDNAList)[2] for KZFP in KZFPDNAList], \
    [getKRABBDiv(KZFP, KRABBDivDNAList)[3] for KZFP in KZFPDNAList], \
    [KZFP.annotations["Linker"] for KZFP in KZFPDNAList], \
    [getLinkerSequence(KZFP, linkerDNAList) for KZFP in KZFPDNAList], \
    [getLinkerStart(KZFP, linkerDNAList) for KZFP in KZFPDNAList], \
    [getLinkerEnd(KZFP, linkerDNAList) for KZFP in KZFPDNAList], \
    [getLinkerLength(KZFP, linkerDNAList) for KZFP in KZFPDNAList], \
    [KZFP.annotations["zfArrayContained"] for KZFP in KZFPDNAList], \
    [getOutOfFrameZFArrayStart(KZFP, outOfFrameZFArrayListDNA) for KZFP in KZFPDNAList], \
    [getOutOfFrameZFArrayEnd(KZFP, outOfFrameZFArrayListDNA) for KZFP in KZFPDNAList], \
    [getOutOfFrameZFArrayLength(KZFP, outOfFrameZFArrayListDNA) for KZFP in KZFPDNAList], \
    [getZincFingers(KZFP, outOfFrameZFArrayListDNA, singleZFDNAList) for KZFP in KZFPDNAList], \
    [getKZFPFingerprint(KZFP, outOfFrameZFArrayListDNA) for KZFP in KZFPDNAList], \
    [getKZFPNegativeFingerprint(KZFP, outOfFrameZFArrayListDNA) for KZFP in KZFPDNAList], \
    [len([szfFingerprint for szfFingerprint in getKZFPFingerprint(KZFP, outOfFrameZFArrayListDNA) if "X" not in szfFingerprint]) for KZFP in KZFPDNAList], \
    [getInFrameZFArrayNumber(KZFP, outOfFrameZFArrayListDNA) for KZFP in KZFPDNAList], \
    [getInFrameZFAContained(KZFP, outOfFrameZFArrayListDNA) for KZFP in KZFPDNAList], \
    [getSingleZFContained(KZFP, outOfFrameZFArrayListDNA) for KZFP in KZFPDNAList], \
    [len([szfFingerprint for szfFingerprint in getKZFPFingerprint(KZFP, outOfFrameZFArrayListDNA)]) for KZFP in KZFPDNAList], \
    [len([szfFingerprint for szfFingerprint in getKZFPFingerprint(KZFP, outOfFrameZFArrayListDNA) if "X" not in szfFingerprint]) for KZFP in KZFPDNAList], \
    [len([szfFingerprint for szfFingerprint in getKZFPFingerprint(KZFP, outOfFrameZFArrayListDNA) if "X" not in szfFingerprint])/len([szfFingerprint for szfFingerprint in getKZFPFingerprint(KZFP, outOfFrameZFArrayListDNA)]) for KZFP in KZFPDNAList]]).T

KZFPTable.columns = ["Genome", "chromosome", "DNAStart", "DNAEnd", "GeneName", "Strand", "id", "GeneID", "GeneType", \
                    "KRABid", "KRABStart", "KRABEnd", "KRABLength", "KRABSequence", "ContainsKRABBbox", "KRABBboxid", \
                    "KRABBboxStart", "KRABBboxEnd", "KRABBboxLength", "ContainsKRABBDiv", "KRABBDivid", "KRABBDivStart", \
                    "KRABBDivEnd", "KRABBDivLength", "Spacerid", "SpacerSequence", "SpacerStart", "SpacerEnd", "SpacerLength", \
                    "OutOfFrameZFAid", "OutOfFrameZFAStart", "OutOfFrameZFAEnd", "OutOfFrameZFALength", "ZincFingers", "Fingerprint", \
                    "NegativeFingerprint", "CanonicalsZF", "numberOfInFrameZFAs", "inFrameZFAsContained", "singleZFsContained", \
                    "numberOfSingleZFs", "nonCanonicalsZFs", "CanonicalScore"]

KZFPTable.index=[KZFP.id for KZFP in KZFPDNAList]

print(str(KZFPTable["Strand"]))

#saving the table in python binary
with open("KZFPTable.dat", "wb") as f:
    pickle.dump(KZFPTable, f, protocol = 2)
#saving the table in csv
KZFPTable.to_csv(path_or_buf = "./Tables/KZFPTable.csv", sep = "\t")


#Saving useful info about detected domains***************************************************************************************************
#*************************************************************KRABS**************************************************************************
KRABTable = pd.DataFrame.from_records([[KRAB.annotations["Genome"] for KRAB in KRABDNAList], \
    [getChromosome(KRAB) for KRAB in KRABDNAList], \
    [KRAB.annotations["DNAStart"] for KRAB in KRABDNAList], \
    [KRAB.annotations["DNAEnd"] for KRAB in KRABDNAList], \
    [KRAB.id for KRAB in KRABDNAList], \
    [KRAB.annotations["DNAEnd"] - KRAB.annotations["DNAStart"] for KRAB in KRABDNAList], \
    [KRAB.annotations["Strand"] for KRAB in KRABDNAList], \
    [KRAB.annotations["DNASequence"] for KRAB in KRABDNAList], \
    [str(KRAB.seq) for KRAB in KRABDNAList], \
    [checkForStopCodon(KRAB) for KRAB in KRABDNAList]]).T
if not KRABTable.empty:
    KRABTable.columns = ["Genome", "chromosome", "DNAStart", "DNAEnd", "id", "DNALength", "Strand", "DNASequence", "AASequence", "ContainsStopCodon"]
    KRABTable.index = [KRAB.id for KRAB in KRABDNAList]
    #saving the table in csv
    KRABTable.to_csv(path_or_buf = "./Tables/KRABTable.csv", sep = "\t")

#*************************************************************KRABBboxes**************************************************************************
#KRABBDivDNAList
KRABBboxesTable = pd.DataFrame.from_records([[KRABBbox.annotations["Genome"] for KRABBbox in KRABBboxDNAList], \
    [getChromosome(KRABBbox) for KRABBbox in KRABBboxDNAList], \
    [KRABBbox.annotations["DNAStart"] for KRABBbox in KRABBboxDNAList], \
    [KRABBbox.annotations["DNAEnd"] for KRABBbox in KRABBboxDNAList], \
    [KRABBbox.id for KRABBbox in KRABBboxDNAList], \
    [KRABBbox.annotations["DNAEnd"] - KRABBbox.annotations["DNAStart"] for KRABBbox in KRABBboxDNAList], \
    [KRABBbox.annotations["Strand"] for KRABBbox in KRABBboxDNAList], \
    [KRABBbox.annotations["DNASequence"] for KRABBbox in KRABBboxDNAList], \
    [str(KRABBbox.seq) for KRABBbox in KRABBboxDNAList], \
    [checkForStopCodon(KRABBbox) for KRABBbox in KRABBboxDNAList]]).T
if not KRABBboxesTable.empty:
    KRABBboxesTable.columns = ["Genome", "chromosome", "DNAStart", "DNAEnd", "id", "DNALength", "Strand", "DNASequence", "AASequence", "ContainsStopCodon"]
    KRABBboxesTable.index = [KRABBbox.id for KRABBbox in KRABBboxDNAList]
    #saving the table in csv
    KRABBboxesTable.to_csv(path_or_buf = "./Tables/KRABBboxTable.csv", sep = "\t")

#*************************************************************KRABBDiv**************************************************************************
#KRABBDivDNAList
KRABBDivTable = pd.DataFrame.from_records([[KRABBDiv.annotations["Genome"] for KRABBDiv in KRABBDivDNAList], \
    [getChromosome(KRABBDiv) for KRABBDiv in KRABBDivDNAList], \
    [KRABBDiv.annotations["DNAStart"] for KRABBDiv in KRABBDivDNAList], \
    [KRABBDiv.annotations["DNAEnd"] for KRABBDiv in KRABBDivDNAList], \
    [KRABBDiv.id for KRABBDiv in KRABBDivDNAList], \
    [KRABBDiv.annotations["DNAEnd"] - KRABBDiv.annotations["DNAStart"] for KRABBDiv in KRABBDivDNAList], \
    [KRABBDiv.annotations["Strand"] for KRABBDiv in KRABBDivDNAList], \
    [KRABBDiv.annotations["DNASequence"] for KRABBDiv in KRABBDivDNAList], \
    [str(KRABBDiv.seq) for KRABBDiv in KRABBDivDNAList], \
    [checkForStopCodon(KRABBDiv) for KRABBDiv in KRABBDivDNAList]]).T
if not KRABBDivTable.empty:
    KRABBDivTable.columns = ["Genome", "chromosome", "DNAStart", "DNAEnd", "id", "DNALength", "Strand", "DNASequence", "AASequence", "ContainsStopCodon"]
    KRABBDivTable.index = [KRABBDiv.id for KRABBDiv in KRABBDivDNAList]
    #saving the table in csv
    KRABBDivTable.to_csv(path_or_buf = "./Tables/KRABBDivTable.csv", sep = "\t")
#*************************************************************singleZFs**************************************************************************
singleZFsTable = pd.DataFrame.from_records([[sZF.annotations["Genome"] for sZF in singleZFDNAList], \
    [getChromosome(sZF) for sZF in singleZFDNAList], \
    [sZF.annotations["DNAStart"] for sZF in singleZFDNAList], \
    [sZF.annotations["DNAEnd"] for sZF in singleZFDNAList], \
    [sZF.id for sZF in singleZFDNAList], \
    [sZF.annotations["DNAEnd"] - sZF.annotations["DNAStart"] for sZF in singleZFDNAList], \
    [sZF.annotations["Strand"] for sZF in singleZFDNAList], \
    [sZF.annotations["DNASequence"] for sZF in singleZFDNAList], \
    [str(sZF.seq) for sZF in singleZFDNAList], \
    [checkForStopCodon(sZF) for sZF in singleZFDNAList], \
    [findsingleZFFingerprint(sZF)[0] for sZF in singleZFDNAList], \
    [findsingleZFFingerprint(sZF)[1] for sZF in singleZFDNAList]]).T
if not singleZFsTable.empty:
    singleZFsTable.columns = ["Genome", "chromosome", "DNAStart", "DNAEnd", "id", "DNALength", "Strand", "DNASequence", "AASequence", "ContainsStopCodon", "Fingerprint", "NegativeFingerprint"]
    singleZFsTable.index = [sZF.id for sZF in singleZFDNAList]
    #saving the table in csv
    singleZFsTable.to_csv(path_or_buf = "./Tables/singleZFsTable.csv", sep = "\t")
#*************************************************************ifZFAs**************************************************************************
inFrameZFArrayTable = pd.DataFrame.from_records([[ifZFA.annotations["Genome"] for ifZFA in inFrameZFArrayListDNA], \
    [getChromosome(ifZFA) for ifZFA in inFrameZFArrayListDNA], \
    [ifZFA.annotations["DNAStart"] for ifZFA in inFrameZFArrayListDNA], \
    [ifZFA.annotations["DNAEnd"] for ifZFA in inFrameZFArrayListDNA], \
    [ifZFA.id for ifZFA in inFrameZFArrayListDNA], \
    [ifZFA.annotations["DNAEnd"] - ifZFA.annotations["DNAStart"] for ifZFA in inFrameZFArrayListDNA], \
    [ifZFA.annotations["Strand"] for ifZFA in inFrameZFArrayListDNA], \
    [ifZFA.annotations["DNASequence"] for ifZFA in inFrameZFArrayListDNA], \
    [str(ifZFA.seq) for ifZFA in inFrameZFArrayListDNA], \
    [checkForStopCodon(ifZFA) for ifZFA in inFrameZFArrayListDNA], \
    [ifZFA.annotations["Fingerprint"] for ifZFA in inFrameZFArrayListDNA], \
    [ifZFA.annotations["singleZFContained"] for ifZFA in inFrameZFArrayListDNA]]).T
if not inFrameZFArrayTable.empty:
    inFrameZFArrayTable.columns = ["Genome", "chromosome", "DNAStart", "DNAEnd", "id", "DNALength", "Strand", "DNASequence", "AASequence", "ContainsStopCodon", "Fingerprint", "SingleZFContained"]
    inFrameZFArrayTable.index = [ifZFA.id for ifZFA in inFrameZFArrayListDNA]
    #saving the table in csv
    inFrameZFArrayTable.to_csv(path_or_buf = "./Tables/inFrameZFArrayTable.csv", sep = "\t")

#*************************************************************oofZFAs**************************************************************************
outOfFrameZFArrayTable = pd.DataFrame.from_records([[oofZFA.annotations["Genome"] for oofZFA in outOfFrameZFArrayListDNA], \
    [getChromosome(oofZFA) for oofZFA in outOfFrameZFArrayListDNA], \
    [oofZFA.annotations["DNAStart"] for oofZFA in outOfFrameZFArrayListDNA], \
    [oofZFA.annotations["DNAEnd"] for oofZFA in outOfFrameZFArrayListDNA], \
    [oofZFA.id for oofZFA in outOfFrameZFArrayListDNA], \
    [oofZFA.annotations["DNAEnd"] - oofZFA.annotations["DNAStart"] for oofZFA in outOfFrameZFArrayListDNA], \
    [oofZFA.annotations["Strand"] for oofZFA in outOfFrameZFArrayListDNA], \
    [str(oofZFA.seq) for oofZFA in outOfFrameZFArrayListDNA], \
    [oofZFA.annotations["Fingerprint"] for oofZFA in outOfFrameZFArrayListDNA], \
    [oofZFA.annotations["singleZFContained"] for oofZFA in outOfFrameZFArrayListDNA], \
    [oofZFA.annotations["inFrameZFAContained"] for oofZFA in outOfFrameZFArrayListDNA]
    ]).T
if not outOfFrameZFArrayTable.empty:
    outOfFrameZFArrayTable.columns = ["Genome", "chromosome", "DNAStart", "DNAEnd", "id", "DNALength", "Strand", "DNASequence", "Fingerprint", "SingleZFContained", "InFrameZFArrayContained"]
    outOfFrameZFArrayTable.index = [oofZFA.id for oofZFA in outOfFrameZFArrayListDNA]
    #saving the table in csv
    outOfFrameZFArrayTable.to_csv(path_or_buf = "./Tables/outOfFrameZFArrayTable.csv", sep = "\t")
