import argparse
import mappy as mp



def parseArgs() :
	parser = argparse.ArgumentParser(description='compare assemblies with minimap2')
	parser.add_argument('-r', help='input assembly 1 (ref) fasta', dest='asm1Filename',required=True)
	parser.add_argument('-q', help='input assembly 2 (query) fasta', dest='asm2Filename',required=True)
	parser.add_argument('-m', help='minimum gap size, default=25', dest='minGapSize',required=False, default=25, type=int)
	parser.add_argument('-n', help='minimum sequence length in assembly 2 (query), default=20Mb', dest='minQueryLen',required=False, default=20000000, type=int)
	parser.add_argument('-o', help='output prefix, default=output', dest='outputPrefix', required=False, default='output')

	arguments = parser.parse_args()
	print("arguments:")
	print("-r, input assembly 1 (ref) fasta: %s" % arguments.asm1Filename)
	print("-q, input assembly 2 (query) fasta: %s" % arguments.asm2Filename)
	print("-m, minimum gap size: %i" % arguments.minGapSize)
	print("-n, minimum sequence length in assembly 2: %i" % arguments.minQueryLen)
	print("-o, output prefix: %s" % arguments.outputPrefix)
	
	print("-- \n")
		
	return arguments

def runIndex(refFastaFilname) :
	print("building index for assembly 1 (ref): %s\n" % refFastaFilname)
	a = mp.Aligner(refFastaFilname, preset="asm5")
	return a
	

	
def runMapper(referenceIndex, asm2Filename, minQueryLen) :
	print("running minimap2 and finding top hit per query sequence\n")             
	hits=[]
	for name, seq, qual in mp.fastx_read(asm2Filename):
		print("... query: %s" % name)
		if len(seq) < minQueryLen:
			print("...... Skipping, query too short (seq len of %i is less than minimum: %i)\n" % (len(seq), minQueryLen))
			continue
		                 
		for hit in referenceIndex.map(seq):
			
			hits.append(name+"\t"+str(len(seq))+"\t"+str(hit))

	return (hits)

def writeOutput(outputPrefix, scaffoldMapListOut):
	print("writing results")
	with open(outputPrefix + "_hit_results", 'w') as r:
		for i in scaffoldMapListOut:
			r.write(i+"\n")

def run() :
	"run script"
	# parse arguments
	args = parseArgs()
	
	# generate index for reference fasta
	refIndex = runIndex(refFastaFilname = args.asm1Filename)
	
	
	scaffoldMapList = runMapper(referenceIndex = refIndex, asm2Filename = args.asm2Filename, minQueryLen = args.minQueryLen)

	
	writeOutput(outputPrefix=args.outputPrefix,scaffoldMapListOut=lsscaffoldMapList)

	print("DONE")

if __name__ == '__main__':
    run()
