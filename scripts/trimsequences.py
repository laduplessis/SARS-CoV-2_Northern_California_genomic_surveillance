import os, sys, time, csv
import numpy as np
from optparse import OptionParser
from Bio import SeqIO, AlignIO
from dropcolumns import dropcolumns


# Remove columns from alignment according to two criteria
###################################################################################################
#
# - All positions marked with an "X" in a fasta sequence file 
#	(e.g. sites under selection for resistance mutations)
# - All positions where more than a certain proportion of sequences have an unknown or ambiguous 
#	assignment
# - If an alignment file is provided the columns are directly trimmed from it
#
#
###################################################################################################
# Parameters
###################################################################################################

usage = "usage: %prog [option]"
parser = OptionParser(usage=usage)

parser.add_option("-H","--histfile",
                  dest = "histfile",
                  default = "",
                  metavar = "path",
                  help = "Path to input csv file with histogram of alignment [required]")

parser.add_option("-f","--fastafile",
                  dest = "fastafile",
                  default = "",
                  metavar = "path",
                  help = "Path to template fasta file with columns to delete marked [optional]")

parser.add_option("-a","--alignfile",
                  dest = "alignfile",
                  default = "",
                  metavar = "path",
                  help = "Path to input alignment file to extract sequences from [optional]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "path",
                  help = "Path to save output file in [required]")

parser.add_option("-p","--prefix",
                  dest = "prefix",
                  default = "output",
                  metavar = "string",
                  help = "Prefix for the output file(s) [default: %default] [required]")

parser.add_option("-c","--cutoff",
                  dest = "cutoff",
                  default = "0.75",
                  metavar = "float",
                  help = "Cutoff proportion for ambiguous columns [default: %default] [required]")

parser.add_option("-m","--markedchar",
                  dest = "markedchar",
                  default = "X",
                  metavar = "character",
                  help = "Character to search for in template fasta file [default: %default] [required]")

parser.add_option("-I","--inverse",
                  dest = "inverse",
                  action = "store_true",
                  help = "Do the inverse (remove all but the selected columns) [default: %default] [optional]")


(options,args) = parser.parse_args()

histfile   = os.path.abspath(options.histfile)  if options.histfile  != "" else ""
fastafile  = os.path.abspath(options.fastafile) if options.fastafile != "" else ""
alignfile  = os.path.abspath(options.alignfile) if options.alignfile != "" else ""
outputpath = os.path.abspath(options.outputpath)+"/"
prefix     = options.prefix
markedchar = options.markedchar
cutoff     = float(options.cutoff)
inverse    = options.inverse

###################################################################################################	


def getHistDropCols(msahist, checkchars=["N", "-"], cutoff=0.75):
	"""Read histogram and return all columns in the alignment with a proportion > cutoff of 
	sequences with one of checkchars at the site

	"""

	histfile = open(msahist, "r")
	histraw  = list(csv.reader(histfile, delimiter=","))
	colnames = histraw[0]
	histmat  = np.array(histraw[1:]).astype("float")
	dropcols = []

	(rows, cols) = histmat.shape
	for i in range(0,rows):
		n = sum(histmat[i,])
		count = 0
		for j in range(0, cols):
			if (colnames[j] in checkchars):
				count += histmat[i,j]

		prop = count/float(n)
		if (prop > cutoff):
			dropcols.append(i)

	return(dropcols)
#


def getFastaDropCols(fastafile, checkchar="X"):
	"""Read sequence from fastafile and return all sites marked with checkchar

	"""

	template = SeqIO.read(fastafile, "fasta")
	dropcols = []
	for i in range(0, len(template.seq)):
		if (template.seq[i] == "X"):
			dropcols.append(i)

	return(dropcols)

#

###################################################################################################	

start = time.time()

dropcols = []

# Read alignment histogram and find all positions with unknown/ambiguous > cutoff
if (histfile != ""):
	ambiguous = ["R", "Y", "W", "S", "K", "M", "D", "V", "H", "B", "N", "-"]
	dropcols += getHistDropCols(histfile, checkchars=ambiguous, cutoff=cutoff)

# Read fasta sequence and find positions with "X"
if (fastafile != ""):
	dropcols += getFastaDropCols(fastafile, checkchar=markedchar)
dropcols = sorted(set(dropcols))

if (inverse): 
	sys.stdout.write("%d columns found to keep in the alignment\n" % len(dropcols))	
	fname = "kept"
else: 
	sys.stdout.write("%d columns found to trim from alignment\n" % len(dropcols))
	fname = "dropped"

# Save output
if (not os.path.exists(outputpath)):
	os.mkdir(outputpath)

with open(outputpath+prefix+"."+fname+".csv","w") as f:
    f.write(",".join(map(str, dropcols))+"\n")

# Remove from alignment
if (alignfile != ""):
	msa = AlignIO.read(alignfile, "fasta")

	if (inverse):
		n = len(msa[0])
		savecols = dropcols
		dropcols = list(range(0,n))
		for i in savecols:
			dropcols.remove(i)

	msa = dropcolumns(msa, dropcols, verbose=False)
	AlignIO.write(msa, outputpath+prefix+".fas", "fasta")


end = time.time()
sys.stdout.write("Total time taken: "+str(end-start)+" seconds\n")

