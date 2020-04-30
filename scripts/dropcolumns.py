
import os, sys, time
from optparse import OptionParser
from Bio import AlignIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

# Given alignment and list of columns removes those columns from the alignment
###################################################################################################
#
# As standalone program:
# - List of columns can be provided as a comma-separated file
# - Alignment should be in fasta format
#
###################################################################################################


def dropcolumns(msa, dropcols, verbose=False):
	"""Remove list of columns from an msa and return a new msa object

	   Start numbering columns from 0
	"""

	# Create blank alignment
	blanks = []
	for seq in msa:
		blanks.append(SeqRecord(Seq("", generic_dna), id=seq.id, description=""))
	newmsa = MultipleSeqAlignment(blanks)


	# Build new alignment
	index = 0
	for i in sorted(dropcols):
		if (index != i):
			if (verbose):
				sys.stdout.write("[%d,%d)\n" % (index,i))
			newmsa += msa[:,index:i]
		index = i+1
	if (index < len(msa[0])):
		if (verbose):
			sys.stdout.write("[%d,%d)\n" % (index,len(msa[0])))
		newmsa += msa[:,index:]

	# Remove <unknown description> 
	for seq in newmsa:
		seq.description = ""

	return(newmsa)
#


def testdropcolumns(teststr, dropcols):
	"""Test the function

	"""

	seqs = []
	for i in range(0,10):
	 	seqs.append(SeqRecord(Seq(teststr, generic_dna), id="seq%d" % i))
	msa = MultipleSeqAlignment(seqs)
	print(msa)

	msa = dropcolumns(msa, dropcols, verbose=True)
	print(msa)
#



###################################################################################################

if __name__ == "__main__":


	###################################################################################################
	# Parameters
	###################################################################################################

	usage = "usage: %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--inputfile",
	                  dest = "inputfile",
	                  default = "",
	                  metavar = "path",
	                  help = "Path to input csv file with columns to drop [required]")

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

	parser.add_option("-t","--test",
	                  dest = "test",
	                  action = "store_true",
	                  help = "Run a test [default: %default] [optional]")


	(options,args) = parser.parse_args()

	inputfile  = os.path.abspath(options.inputfile)
	alignfile  = os.path.abspath(options.alignfile) if options.alignfile != "" else ""
	outputpath = os.path.abspath(options.outputpath)+"/"
	prefix     = options.prefix
	test       = options.test

	###################################################################################################	

	start = time.time()

	if (test): 
		# Test cases
		testdropcolumns(("A"*3)+"GCG"+("A"*10)+"GCGA", [17,4])
		testdropcolumns("CG"+("A"*3)+"GCCCCG"+("A"*10)+"GCG", [0,6,7,8,9,22])
		testdropcolumns("CCCG"+("A"*5)+"GC", [0,1,2,10])
		testdropcolumns("CCCG"+("A"*5)+"GCCC", [0,1,2,10,11,12])
		testdropcolumns(("A"*5)+"GCCCG"+("A"*5), [6,7,8])
	else:

		# Read input
		msa = AlignIO.read(alignfile, "fasta")

		with open(inputfile, 'r') as f:
			templist = []
			for line in f:
				templist += line.strip().split(',')
		dropcols = list(map(int, templist))

		print(dropcols)
		print(len(dropcols))
		
		# Drop columns
		msa = dropcolumns(msa, dropcols, verbose=False)

		# Save output
		if (not os.path.exists(outputpath)):
			os.mkdir(outputpath)

		AlignIO.write(msa, outputpath+prefix+".fas", "fasta")
	#

	end = time.time()
	sys.stdout.write("Total time taken: "+str(end-start)+" seconds\n")
