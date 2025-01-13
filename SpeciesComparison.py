import pickle
import pandas as pd
from Functions import getHomology

#species1 and species2 values are cosmetics and only influence the name of the final homology file
#path 1 and path2: indicate the position of the KZFPTable.dat file for each of the two species you want to compare. Note that if you use the same KZFPTable.dat,
#the script understands it as the same species and drops duplicates automatically (e.g. does not score ZNF221 to ZNF221 homology to 1 but skips it instead)
#canonicalThreshold: (between 0 and 1, KZFPs with a canonical score =[proportion of canonical zinc fingers to all zinc fingers] 
#lower than the threshold will be ignored)
species1 = "hg19"
species2 = "hg19"
path1 = "./KZFPTable.dat"
path2 = "./KZFPTable.dat"
canonicalThreshold = 0

#by default, 0 automatically allocates processes
nProcesses = 0

homologyTable = getHomology(path1, path2, canonicalThreshold, False, nProcesses)

#saving results
homologyTable.to_csv(path_or_buf = "./Tables/{}{}Homology.csv".format(species1, species1), sep = "\t")

with open("./SpeciesComparisons/{}{}HomologyTable.dat".format(species1, species2), "wb") as f:
   pickle.dump(homologyTable, f)
