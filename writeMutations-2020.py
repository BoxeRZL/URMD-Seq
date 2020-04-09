import sys
import os
from Bio import SeqIO, AlignIO
import shutil
from Bio.Align import AlignInfo
import csv
import glob

#***is it >75% or >=75% ?????

def main():		
##USAGE -- python writeCons-6-jan.py ./MID-folders-location grpSize(5) consensus(0.75) cambridgeRef.txt
##identifying mutations post grouping analysis ()
	print(str(sys.argv))
	#midList = ['3','4','5','6','7','10','13','15','16','17', '19', '20','21','22','24','30','33','36','39','40']
	#	midList = [("AGACGCACTC","3"), ("AGCACTGTAG","4"), ("ATCAGACACG","5"), ("ATATCGCGAG","6"), ("CGTGTCTCTA", "7"), ("CATAGTAGTG","13"), ("ATACGACGTA","15"), ("TCACGTACTA","16"), ("CGTCTAGTAC","17"), ("TGTACTACTC","19"), ("ACGACTACAG","20"),("CGTAGACTAG","21"), ("TACGAGTATG","22"), ("TAGAGACGAG","24"), ("AGACTATACT","30"), ("ATAGAGTACT","33"), ("CGACGTGACT","36"), ("TACAGATCGT","39"), ("TACGCTGTCT","40")]
	midList = ['15','1']#,'4','5','6','7','13','15','16','17']

	outPath = str(sys.argv[1])+"/mutationSummary.csv"
	try:
		os.remove(outPath)
	except:
		print("old output file removed")

	outputFile = open(outPath, 'w')
	outputFile.write(','.join(['position', 'mid', 'reference','coverage', 'a', 't', 'c', 'g', '-', 'not indel', 'identified mutation?']))
	outputFile.write("\n")
	outputFile.close()

	for mid in midList:
		midsRoot = sys.argv[1]
		filepath = midsRoot+"/"+mid

		try:
			unzipCall = "unzip -q "+str(filepath)+"/alignedGroups"+str(mid)+".zip -d "+str(filepath)+"/"
			print(str(unzipCall))
			zipReturn = os.system(unzipCall)
			maxGrps = len(glob.glob1(filepath,"aligned*.aln"))		
			mkSubCons(filepath, mid, maxGrps, int(sys.argv[2]), float(sys.argv[3]), str(sys.argv[4]), outPath)
			
		except:
			print("Unable to find zipped file")
		"""int(sys.argv[2]),"""
		#CCTCGAGGTCGACGGTATCGACGAGTGCGTNNNNCANNNNGTNNNNNgcccacacgttccccttaaataaga
		try:
			rezipCall = "zip -mq alignedGroups"+str(mid)+".zip "+str(filepath)+"/aligned*.aln"
			print(rezipCall)
			zipReturn = os.system(rezipCall)
			#print("shutil.rmtree(filepath)")
		except:
			print("could not execute zipping")


def mkSubCons(filepath, mid, maxGrpNo, grpSizeThreshold, consThreshold, refSeqPath, outPath):        ##make subcons, input filepath and the total number of groups
	print("getsubcons max "+str(maxGrpNo)+"  group size: "+str(grpSizeThreshold)+"  consensus threshold: "+str(consThreshold))
	print("analyzing MID "+str(mid))

	consIn = "/Consensus_long_"+str(grpSizeThreshold)+"_"+str(consThreshold)+".txt"
	consOut = "/Consensus_long_"+str(grpSizeThreshold)+"_"+str(consThreshold)+".aln"
	try:
		os.remove(filepath+consIn)
		os.remove(filepath+consOut)
	except:
		print("nothing to remove")

	#posnCount = 16534
	writeFile = open(filepath+consIn, 'a')
	refSeq = SeqIO.read(refSeqPath, 'fasta')

	writeFile.write(">refSeq\n"+str(refSeq.seq)+"\n")
	print(str(maxGrpNo)+"\n")
	for i in range (1, maxGrpNo+1):
		alnRaw = AlignIO.read(filepath+"/aligned"+str(i)+".aln", "fasta")		##
		if (len(alnRaw) >=grpSizeThreshold):										##if groups are over threshold size
			alnObj = AlignInfo.SummaryInfo(alnRaw)
			cons = alnObj.gap_consensus(threshold=consThreshold, ambiguous="-")	##consensus for bases over consensus threshold

			#trimPosn = cons.seq.find("AATAAGA")								##trim degenerate primers
			#trimmedCons = cons[trimPosn:]
			trimmedCons = cons[51:400]		#16559 < x < 280
			#print(str(trimPosn)+"  "+str(i)+"  "+str(trimmedCons))				#print each subcons
			writeFile.write(">group"+str(i)+"\n"+str(trimmedCons).upper()+"\n")
		print("\r writing groups:  {0:.2f} %".format(i/maxGrpNo*100)))

	print("closing aln cons")
	writeFile.close()



	cmdCall = "mafft --retree 2 --quiet --maxiterate 10 --op 0.05 --ep 0.05 "+filepath+consIn+" > "+filepath+consOut
	print(str(cmdCall))
	returnOut = os.system(cmdCall)

	mutPrint(filepath, mid, consOut, outPath)


def mutPrint(filepath, mid, rawFile, outPath):

	posn = 16533				#refseq before start site -1 for zero indexing in python. ***ln 112 check indexing
	chrSize = 16568				#mtDNA size (16568)-1 for zero indexing
	insert = 0

	posnMin = 34
	posnMax = 280				## **** change region of interest start and end position here

	"""try:						#move this to main, plus header write for output
		os.remove(outPath)
		outputFile
	except:
		print("no output to remove")"""	
	refAlignment = AlignIO.read(filepath+rawFile, "fasta")

	outputFile = open(outPath, 'a')
	#outputFile.write(','.join(['position', 'mid', 'reference','coverage', 'a', 't', 'c', 'g', '-', 'mutation?']))
	#outputFile.write("\n")

	for base in range(0,len(refAlignment[0])):
		if (refAlignment[0][base] == "-"):
			insert += 1
		else:
			posn += 1
			insert = 0

		a = refAlignment[1:,base].count("a")
		t = refAlignment[1:,base].count("t")
		c = refAlignment[1:,base].count("c")
		g = refAlignment[1:,base].count("g")
		indel = refAlignment[1:,base].count("-")

		coverage = a+t+c+g+indel
		#print(str(coverage))

		#if ((a != coverage) or (t != coverage) or (c != coverage) or (g != coverage) or (indel != coverage)):
		if (not((a/coverage == 1) or (t/coverage == 1) or (c/coverage == 1) or (g/coverage == 1) or (indel/coverage == 1))):
		
			#outputFile.write(str(posn+(0.001*insert))+","+str(converage)+","+str(a)+",")
			outputFile.write(",".join([str(((posn+(0.001*insert))%chrSize)),str(mid),str(refAlignment[0][base]),str(coverage),str(a),str(t),str(c),str(g),str(indel)]))
			if ((indel == 0)):
				outputFile.write(",1")
			else:
				outputFile.write(",0")

#** jan 6
			posnAgreement = 0
			if (a > 0): posnAgreement += 1
			if (t > 0): posnAgreement += 1
			if (c > 0): posnAgreement += 1
			if (g > 0): posnAgreement += 1

			if ((posnAgreement > 1) and (posn >= posnMin) and (posn <= posnMax)):
				outputFile.write(",1")
			else:
				outputFile.write(",0")

			outputFile.write("\n")

	outputFile.close()

	#call bash with Consensus.txt

	#readfile = open(filepath+"/Consensus.txt", 'r')
	"""consAlnRaw = AlignIO.read(filepath+consOut, "fasta")

	consSummaryNoAlign = AlignInfo.SummaryInfo(consAlnRaw)
	consSummary = consSummaryNoAlign.gap_consensus(ambiguous="-")
	#consensus = consSummary.gap_consensus()
	#Now, we want to make the PSSM, but ignore any N ambiguity residues when calculating this:

	individualPSSM = consSummaryNoAlign.pos_specific_score_matrix(consSummary)

	writePSSM = open(filepath+"PSSM.txt", 'w')
	print("writing PSSM")



	#for i in len(individualPSSM):
	#	if (individualPSSM[0,i-1] == "-"
	writePSSM.write(str(individualPSSM))
	writePSSM.close()

	writeRef(filepath, consSummary)
"""

def writeRef(filepath, consSummary):

	reference = SeqIO.read("./cambridgeRef.txt", "fasta")

	writeFile.open(filepath+"/refAlign.txt", 'w')
	writeFile.write(">"+str(reference.id)+"\n"+str(reference.seq)+"\n")
	writeFile.write(">grpCons\n"+str(consSummary)+"\n")
	writeFile.close()

	cmdCall = "mafft --retree 2 --maxiterate 10 --op 0.05 --ep 0.05 "+filepath+"refAlign.txt > "+filepath+"refAlign.aln"
	returnOut = os.system(cmdCall)	

	readPosn.open(filepath+"refAlign.aln", 'r')
	refAlignment = readPosn.readlines()

	readPSSM.open(filepath+"PSSM.txt", 'r')
	PSSMlines = readPSSM.readlines()

	posn = 16534
	chrSize = 16568
	for base in len(refAlignment[0]):
		if (refAlignment[0][base] != "-"):				#case ref(1) posn++
			posn += 1									#	  cons(1) insert=0
			#insertC = 0
			if (refAlignment[1][base] != "-"):
				insertC = 0

		else:											#case ref(0) posn=posn
			if (refAlignment[1][base] != "-"):			#	  cons(1) insert++
				insertC +=1
		



	"""	writeFile = open(filepath+"/ConsensusSub"+str(subNo)+".txt", 'a')
		for i in range(c,maxGrpNo):
			#grpFile = filepath+"/aligned"+str(i)+"aln"
			alnRaw = AlignIO.read(filepath+"/aligned"+str(i)+".aln", "fasta")

			alnObj = AlignInfo.SummaryInfo(alnRaw)
			cons = alnObj.gap_consensus(threshold=0.51,ambiguous="-")
			print("write "+str(i))
			writeFile.write(">"+str(i)+"\n"+str(cons).upper()+"\n")

		cmdCall = "mafft --retree 2 --maxiterate 10 --op 0.05 --ep 0.05 "+filepath+"/ConsensusSub"+str(subNo)+".txt > "+filepath+"/ConsensusSub"+str(subNo)+".aln"
		print(str(cmdCall))
		returnOut = os.system(cmdCall)

		c += 50
		subNo += 1



	i = 0
	writeFile = open(filepath+"/Consensus.txt", 'a')
	for i in range (1, subNo):
		alnRaw = AlignIO.read(filepath+"/ConsensusSub"+str(i)+".aln", "fasta")
		alnObj = AlignInfo.SummaryInfo(alnRaw)
		cons = alnObj.gap_consensus(threshold=0.51, ambiguous="-")
		print(str(i)+"  "+str(cons))				#print each subcons
		writeFile.write(">SubConsensus"+str(i)+"\n"+str(cons).upper()+"\n")

	writeFile.close()
	cmdCall = "mafft --retree 2 --maxiterate 10 --op 0.05 --ep 0.05 "+filepath+"/Consensus.txt > "+filepath+"/Consensus.aln"
	print("complete consensus calc: "+str(cmdCall))
	returnOut = os.system(cmdCall)

	subNo -= 1

	return subNo


	"""

if __name__=='__main__':
	main()

