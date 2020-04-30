import os, sys, time
from optparse import OptionParser
from Bio import AlignIO
from Bio.Seq import Seq


# Script to convert alignments between formats
# BioPython doesn't write Phylip files with full and padded filenames.
#
# Can also change all ambiguity codes to N (option -r) and change all 
# sequences to uppercase (option -u)
####################################################################################################################################

if __name__ == "__main__":

	################################################################################################################################
	# Parameters
	################################################################################################################################
	usage = "usage: %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--inputfile",
	                  dest = "inputfile",
	                  default = "",
	                  metavar = "path",
	                  help = "Path to input alignment [required]")

	parser.add_option("-o","--outputfile",
	                  dest = "outputfile",
	                  default = "",
	                  metavar = "path",
	                  help = "Path to output file in [required]")

	parser.add_option("-1","--informat",
	                  dest = "informat",
	                  default = "fasta",
	                  metavar = "string",
	                  help = "File format of input alignment [default: %default] [optional]")

	parser.add_option("-2","--outformat",
	                  dest = "outformat",
	                  default = "fasta",
	                  metavar = "string",
	                  help = "File format of sequence files [default: %default] [optional]")

	parser.add_option("-a","--alphabet",
	                  dest = "alphabet",
	                  default = "A,C,G,T,R,Y,W,S,K,M,D,V,H,B,N,-,_",
	                  metavar = "string",
	                  help = "Comma-separated list of all characters (nucleotides) allowed in an alignment [default: %default] [optional]")

	parser.add_option("-r","--removeAmbiguities",
	                  dest = "removeAmbiguities",
	                  action = "store_true",
	                  help = "Change all ambiguity codes to N [default: %default] [optional]")

	parser.add_option("-u","--uppercase",
	                  dest = "uppercase",
	                  action = "store_true",
	                  help = "Change all nucleotides to upper case [default: %default] [optional]")

	(options,args) = parser.parse_args()

	inputfile     = os.path.abspath(options.inputfile)
	outputfile    = os.path.abspath(options.outputfile)
	informat      = options.informat
	outformat     = options.outformat
	remove        = options.removeAmbiguities
	upper         = options.uppercase
	
	ambiguities = ["R", "Y", "W", "S", "K", "M", "D", "V", "H", "B"]

	start = time.time()

	msa = AlignIO.read(inputfile, informat)


	if (upper):
		for i in range(0, len(msa)):
			msa[i].seq = msa[i].seq.upper()


	if (remove):
		for i in range(0, len(msa)):
			sequence = str(msa[i].seq)
			for c in ambiguities:		
				sequence  = sequence.replace(c, "N")

			msa[i].seq = Seq(sequence)



	AlignIO.write(msa, outputfile, outformat)


	end = time.time()
	sys.stdout.write("Total time taken: "+str(end-start)+" seconds\n")



