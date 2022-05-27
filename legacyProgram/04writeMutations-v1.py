import sys
import os
from Bio import SeqIO, AlignIO
import shutil
from Bio.Align import AlignInfo
import csv
import glob
from datetime import date
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# python 04writeMutations-v1.py $dataLocation/$laneName-$trimQuality-$flashOverlap $groupSize $frequencyThreshold 00cambridgeRefLong.txt $dataLocation $outputFileName $laneName- $WTadjustedFlag

def logDebugging(debugStr):
	debugBool = True
	if (debugBool): print(str(debugStr))

def main():		
	midList = ['1','3','4','5','6','7','10','13','14','15','16','17','18','19','20','21','22','24','30','33','36','39','40']

	midList = [str(sys.argv[7])+mid for mid in midList]
	doAdjustment = True
	analysisName = str(sys.argv[6])
	outPath = analysisName+str(date.today().strftime("%d-%b-%Y"))+".csv"	#str(sys.argv[6])	#str(sys.argv[1])+"/"+str(sys.argv[5])+"mutationSummary-N-"+str(sys.argv[2])+"-"+str(sys.argv[3])+".csv"
	print(str(outPath))
	try:
		os.remove(str(sys.argv[1]+"/../"+str(analysisName)+"-contamination.txt"))
	except:
		print("no prior contamination file to remove")

	outputFile = open(outPath, 'a')
	outputFile.write(','.join(['groupSize','mutationFreq','lane','position', 'mid', 'reference','coverage', 'a', 't', 'c', 'g', 'n', '-', 'nonWTbases']))
	outputFile.write("\n")
	outputFile.close()

	processingVar = str(','.join([str(sys.argv[2]),str(sys.argv[3]),str(sys.argv[5])]))


	for mid in midList:
		midsRoot = sys.argv[1]
		filepath = midsRoot+"/"+mid

		unzipCall = "unzip -q "+str(filepath)+"/WTgroupsAligned-"+str(mid)+".zip -d "+str(filepath)+"/"
		logDebugging(str(unzipCall))
		zipReturn = os.system(unzipCall)
		logDebugging(str(filepath)+"/WTgroup*.aln")
		maxGrps = len(glob.glob1(filepath,"WTgroup*.aln"))		

		try:
			mkConsWTadjusted(filepath, mid, maxGrps, int(sys.argv[2]), float(sys.argv[3]), str(sys.argv[4]), outPath, processingVar, analysisName, doAdjustment)
			rmUnzippedFiles = "rm "+str(filepath)+"/WTgroup*.aln"
			logDebugging(str(rmUnzippedFiles))
			rmReturn = os.system(rmUnzippedFiles)
		except:
			print("Unable to complete analysis")
		
		rmUnzippedFiles = "rm "+str(filepath)+"/WTgroup*.aln"
		logDebugging(str(rmUnzippedFiles))
		rmReturn = os.system(rmUnzippedFiles)
		#CCTCGAGGTCGACGGTATCGACGAGTGCGTNNNNCANNNNGTNNNNNgcccacacgttccccttaaataag

def mkConsWTadjusted(filepath, mid, maxGrpNo, grpSizeThreshold, consThreshold, refSeqPath, outPath, processingVar, analysisName, doAdjustment):
	#make WTadjConsensus.txt file and align it
	writeFile = open(filepath+"/WTadjConsensus.txt","w+")
	refSeq = SeqIO.read(refSeqPath, "fasta")
	writeFile.write(">refSeq\n"+str(refSeq.seq)+"\n")
	writeFile.close()
	logDebugging("opened consensus file")

	for i in range(1,maxGrpNo+1):
		#adjust for Ns if WT present at >=25% plus contamination filter
		adjWT(filepath, i, grpSizeThreshold, mid, str(refSeq.seq), analysisName, doAdjustment)
	
	cmdCall = "mafft --retree 2 --quiet --maxiterate 10 --op 0.05 --ep 0.05 "+filepath+"/WTadjConsensus.txt > "+filepath+"/WTadjConsensus.aln"
	aligning = os.system(cmdCall)
	try:
		mutPrint(filepath, mid, "/WTadjConsensus.aln", outPath,processingVar)
	except:
		print("could not find required alignment files")

def adjWT(filepath, groupNo, grpSizeThreshold, mid, refSeqStr, analysisName, doAdjustment):
	grpAlnRaw = AlignIO.read(filepath+"/WTgroup"+str(groupNo)+".aln", "fasta")
	grpAlnObj = AlignInfo.SummaryInfo(grpAlnRaw)
	WTgroupCons = str(grpAlnRaw[1].seq).upper()

	for c in range(0,len(WTgroupCons)):
		if (doAdjustment) and (WTgroupCons[c] == "N"):
			wtCount = grpAlnRaw[2:,c].count(grpAlnRaw[0,c])
			wtPrevalence = wtCount/len(grpAlnRaw[2:,c])
			if(wtPrevalence > 0.25):
				WTgroupCons = WTgroupCons[:c] + str(grpAlnRaw[0,c]) + WTgroupCons[c+1:]
	mismatches = 0
	mismatches = str(grpAlnRaw[0].seq).count('-') + str(grpAlnRaw[1].seq).count('-')
	
	sequenceIdentity = 0
	referencePairwise = str(grpAlnRaw[0].seq).replace("-","")[:315]
	#print(referencePairwise)
	groupPairwise = str(grpAlnRaw[1].seq).replace("-","")
	#print(groupPairwise)
	#print(format_alignment(*pairwise2.align.globalxx("attgt","attcg")[0]).count("|"))
	sequenceIdentity = int(format_alignment(*pairwise2.align.globalxx(referencePairwise, groupPairwise)[0]).count("|"))
	if (sequenceIdentity >= 290) and (len(grpAlnRaw[2:,0]) >= grpSizeThreshold):
		#print(str(sequenceIdentity)+"  from group: "+str(groupNo))
		writeFile = open(filepath+"/WTadjConsensus.txt", 'a')
		writeFile.write(">"+str(groupNo)+"\n"+str(WTgroupCons).upper()+"\n")
		writeFile.close()
	if (sequenceIdentity < 290):
		print(str(sequenceIdentity)+" from group: "+str(groupNo))
		print(format_alignment(*pairwise2.align.globalxx(referencePairwise, groupPairwise)[0]))
		#logDebugging("opening "+str(filepath)+"/../../contamination.txt")
		writefile = open(filepath+"/../../"+str(analysisName)+"-contamination.txt",'a')
		if (os.stat(filepath+"/../../"+str(analysisName)+"-contamination.txt").st_size == 0): writefile.write(">refSeq\n"+str(refSeqStr)+"\n")
		writefile.write(">"+str(mid)+" group "+str(groupNo)+"\n"+str(WTgroupCons).upper()+"\n")
		writefile.close()

	return 0

def mutPrint(filepath, mid, rawFile, outPath, processingVar):
	posn = 16533				#refseq before start site -1 for zero indexing in python. ***ln 112 check indexing
	chrSize = 16568				#mtDNA size (16568)-1 for zero indexing
	insert = 0

	posnMin = 34
	posnMax = 254				## **** change region of interest start and end position here

	refAlignment = AlignIO.read(filepath+rawFile, "fasta")

	outputFile = open(outPath, 'a')

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
		n = refAlignment[1:,base].count("n")
		indel = refAlignment[1:,base].count("-")
		coverage = a+t+c+g+n+indel

		if (not((a/coverage == 1) or (t/coverage == 1) or (c/coverage == 1) or (g/coverage == 1) or (n/coverage == 1) or (indel/coverage == 1))):
		
			if (posn == 0):
				posn = int(chrSize)+1
			outputFile.write(",".join([str(processingVar),str(((posn+(0.001*insert))%chrSize)),str(mid),str(refAlignment[0][base]),str(coverage),str(a),str(t),str(c),str(g),str(n),str(indel)]))
			posnAgreement = 0
			if (a > 0): posnAgreement += 1
			if (t > 0): posnAgreement += 1
			if (c > 0): posnAgreement += 1
			if (g > 0): posnAgreement += 1
			wtBases = 0

			wtBases=max(a,t,c,g)
			outputFile.write(","+str(coverage-n-indel-wtBases))
			if (insert == 0): outputFile.write(",1")
			else: outputFile.write(",")

			outputFile.write("\n")

	outputFile.close()
	print("# groups:    "+str(len(refAlignment[1:,1])))


if __name__=='__main__':
	main()
