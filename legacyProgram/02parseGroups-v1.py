#!/usr/bin/python

import sys
import os
from Bio import SeqIO, AlignIO
import shutil
from Bio.Align import AlignInfo
import csv
import multiprocessing

# python ./02parseGroups-v1.py $2/$1-$3-$4 $1-$3-$4.extendedFrags.fastq 00cambridgeRefLong.txt m $1
##folder location should only contain raw fasta for analysis

def main():
	filepathSplit = sys.argv[1]								#folderLocation
	filepath = filepathSplit+"/"+str(sys.argv[2])			#raw data
	print(str(filepath))
	refSeqPath = sys.argv[3]
	midPrefix = sys.argv[5]+"-"

	try:													#open raw data file
		lineNo = 0
		fp = open(filepath, "r")
		print("file opened: "+str(filepath))
	except:
		print("ERROR: invalid filepath to raw fasta")

	midList = [("ACGAGTGCGT",midPrefix+"1"), ("AGACGCACTC",midPrefix+"3"), ("AGCACTGTAG",midPrefix+"4"), ("ATCAGACACG",midPrefix+"5"), ("ATATCGCGAG",midPrefix+"6"), ("CGTGTCTCTA",midPrefix+"7"), ("TCTCTATGCG",midPrefix+"10"), ("CATAGTAGTG",midPrefix+"13"), ("CGAGAGATAC",midPrefix+"14"), ("ATACGACGTA",midPrefix+"15"), ("TCACGTACTA",midPrefix+"16"), ("CGTCTAGTAC",midPrefix+"17"),("TCTACGTAGC",midPrefix+"18")]#, ("TGTACTACTC","19"), ("ACGACTACAG","20"),("CGTAGACTAG","21"), ("TACGAGTATG","22"), ("TAGAGACGAG","24"), ("AGACTATACT","30"), ("ATAGAGTACT","33"), ("CGACGTGACT","36"), ("TACAGATCGT","39"), ("TACGCTGTCT","40")]

	try:
		for item in midList:                    			#delete anything in the folder from previous analyses
			shutil.rmtree(filepathSplit+"/"+item[1])
			os.remove(filepathSplit+"/"+item[1]+"groupSize.csv")
		os.remove(filepathSplit+"/mutated.txt")
		print("old analyses data cleaned")
	except:
		print("nothing to clean")

	parseByMid(filepathSplit, filepath, midList)			#for each sequence in lane parse into mid group, everything else is in mutated.txt


	if (sys.argv[4] == 'm'):	
		jobs = []
		for mid in midList:
			p = multiprocessing.Process(target=speedUpWrapper, args=(filepathSplit, filepath, refSeqPath, mid))
			jobs.append(p)
			p.start()
		print(str(jobs))
	else:
		for mid in midList:
			speedUpWrapper(filepathSplit, filepath, refSeqPath, mid)

	print('all fn complete')


def speedUpWrapper(filepathSplit, filepath, refSeqPath, mid):
	try:
		filepathFull = filepathSplit+"/"+mid[1]
		noGrps = mkGrp(filepathFull, mid)      						#from each mid make subgroups and keep track of how many created
		print(noGrps)
		grp = 0

		for grp in range(1, noGrps+1):
			cmdCall = "mafft --retree 2 --quiet --maxiterate 10 --op 0.05 --ep 0.05 "+filepathFull+"/group"+str(grp)+".txt > "+filepathFull+"/aligned"+str(grp)+".aln"
			returnOut = os.system(cmdCall)
			print("\r MID "+mid[1]+" aligning:  {0:.2f} %".format(grp/noGrps*100), end='')#{0:.2f}".format(noSeq/fileLength*100)+" %", end='\r')

		noSubCons = mkCons(filepathFull, noGrps)
		print("SUBCONS FORMED    "+str(noSubCons))

		WTediting = defaultWT(filepathFull, noGrps, 0.75)

	except:
		print("MID not found")

	print("exit: "+str(mid))

def mkGrpCons(noGrps, filepathFull, refSeqPath):

	refSeq = SeqIO.read(refSeqPath, "fasta")
	print("mkGrpCons  " +str(filepathFull)+"  "+str(refSeqPath))

	supAlnRaw = AlignIO.read(filepathFull+"/Consensus.aln", "fasta")
	supAlnObj = AlignInfo.SummaryInfo(supAlnRaw)
	supCons = supAlnObj.gap_consensus(threshold=0.5,ambiguous="-")  
	 
	for grp in range(1,noGrps+1):      #noGrps is list of numbers
						#create file, write refSeq, write super consensus
		print("\r writing ConsGroup:  {0:.2f} %".format(grp/noGrps*100), end='')
		writeFile = open(filepathFull+"/ConsGroup"+str(grp)+".txt", 'w')
		writeFile.write(">refSeq\n"+str(refSeq.seq)+"\n") 
		writeFile.write(">superConsensus\n"+str(supCons)+"\n")

						#get group consensus and write to file
		grpPath = filepathFull+"/aligned"+str(grp)+".aln"
		grpAlnRaw = AlignIO.read(grpPath, "fasta")
		grpAlnObj = AlignInfo.SummaryInfo(grpAlnRaw)
		grpCons = grpAlnObj.gap_consensus(threshold = 0.5, ambiguous="-")
		#print("grp cons check  "+str(grpCons)+"  "+str(grp))
		writeFile.write(">"+str(grp)+"\n"+str(grpCons)+"\n")

		for seqFrag in SeqIO.parse(grpPath, "fasta"):
			writeFile.write(">"+seqFrag.id+"\n"+str(seqFrag.seq)+"\n")
		#print("writing grp:  "+str(grp))

	return noGrps



def parseByMid(filepathSplit, filepath, midList):

	noSeq = 0 										#***may need to return value if used later

	fileLength = len(list(SeqIO.parse(filepath, "fastq")))
	for seqFragment in SeqIO.parse(filepath, "fastq"):				#for each sequence in lane parse into mid group, everything else is in mutated.txt
		noSeq += 1
		found = 0
		for mid in midList:
			posn = seqFragment.seq.find(mid[0])     				#search for mid in seq
			if (posn >= 0):                         				#if mid sequence is found
				if ((str(seqFragment.seq)[posn+14:posn+16] == "CA") and (str(seqFragment.seq)[posn+20:posn+22] == "GT") and (str(seqFragment.seq)[posn+27:posn+29] == "GC")):# and (len(seqFragment.seq) > 300)):		#check for presence of conserved sites in start seq
					seqFragment.description = str(mid[0])
					found = str(mid[1])

		if ((found) and (len(seqFragment.seq) >= 300)):         										#if matching MID is found append seq to 'mid/ungrouped.txt'
			writePath = filepathSplit+"/"+found+"/ungrouped.txt" 
			try:
				writeFile = open(writePath, 'a')
				writeFile.write(">"+seqFragment.id+"\n")
				writeFile.write(str(seqFragment.seq).upper()+"\n")
			except:
				os.mkdir(filepathSplit+"/"+found)
				writeFile = open(writePath, 'w')
				writeFile.write(">"+seqFragment.id+"\n")
				writeFile.write(str(seqFragment.seq).upper()+"\n")
		else:														#if matching MID not found append seq to 'mutated.txt'
			writePath = filepathSplit+"/mutated.txt"
			try:
				writeFile = open(writePath, 'a')
				writeFile.write(">"+seqFragment.id+"\n")
				writeFile.write(str(seqFragment.seq).upper()+"\n")
			except:
				os.mkdir(filepathSplit+"/mutated.txt")
				writeFile = open(writePath, 'w')
				writeFile.write(">"+seqFragment.id+"\n")
				writeFile.write(str(seqFragment.seq).upper()+"\n")

	print("MID searching:  100 %")


def mkGrp(filepath, mid):					#filepath=directory with ungrouped txt file, mid=(midseq, number)
	
	grpList = {}							#"startSeq":("groupNumber","groupSize")}
										#####  "startSeq":["grpNo", "groupSize", seq1, seq2, ... seqN]
	grpNo = 0
	ungrpedFile = filepath+"/ungrouped.txt"
	print(ungrpedFile)
	grpFound = False
	fileLength = len(list(SeqIO.parse(ungrpedFile, "fasta")))
	noSeq = 1

	for seqFragment in SeqIO.parse(ungrpedFile, "fasta"):
		posn = seqFragment.seq.find(mid[0])
		seqStart = str(seqFragment.seq)[posn:posn+29]		#seqStart is the first 29 bases (the PID)
		
		if seqStart in grpList.keys():
			grpFound = True
		else:
			grpFound = False

		if (grpFound):
			grpList[seqStart][1] += 1
			grpList[seqStart].append(seqFragment)
		else:
			grpNo += 1
			grpList[seqStart] = [grpNo, 1, seqFragment]
		noSeq += 1

	print(len(grpList))                         ##### no grps
	sortedK = sorted(grpList.keys())

	grpCounter = 0
	for sortedGrp in sortedK:					##reorder groups alphabetically and keep total
		if (len(grpList[sortedGrp]) > 3):
			grpCounter += 1
			grpList[sortedGrp][0] = grpCounter
		else:
			grpList[sortedGrp][0] = 0

	for sortedGrp in sortedK:
		i = 2
		writeFile = open(filepath+"/group"+str(grpList[sortedGrp][0])+".txt", 'a')

		while i < len(grpList[sortedGrp]):
			writeFile.write(">"+grpList[sortedGrp][i].id+"\n"+str(grpList[sortedGrp][i].seq)+"\n")
			i += 1

	###print group sizes output
	outGrp = {}
	for grp in sortedK:			#categorize all groups by size
		try:
			outGrp[grpList[grp][1]] += 1
		except:
			outGrp[grpList[grp][1]] = 1
	outK = sorted(outGrp.keys())

	writeFile = open(filepath+"groupSize.csv", 'w')
	writeFile.write("Size, No of groups\n")
	for grp in outK:			#print groups by size
		writeFile = open(filepath+"groupSize.csv", 'a')
		writeFile.write(str(grp)+", "+str(outGrp[grp])+"\n")

	return grpCounter

def mkCons(filepath, maxGrpNo):
	print("cons max: "+str(maxGrpNo))

	writeFile = open(filepath+"/Consensus.txt", 'a')
	for i in range(1,maxGrpNo+1):
		alnRaw = AlignIO.read(filepath+"/aligned"+str(i)+".aln", "fasta")
		alnObj = AlignInfo.SummaryInfo(alnRaw)
		cons = alnObj.gap_consensus(threshold=0.51,ambiguous="-")
		writeFile.write(">"+str(i)+"\n"+str(cons[50:]).upper()+"\n")
	
	writeFile.close()
	cmdCall = "mafft --retree 2 --quiet --maxiterate 10 --op 0.05 --ep 0.05 "+filepath+"/Consensus.txt > "+filepath+"/Consensus.aln"
	print("complete consensus calc: "+str(cmdCall))
	returnOut = os.system(cmdCall)

def defaultWT(filepath, maxGrpNo, grpThreshold):
	sampleCons = getCons(filepath+"/Consensus.aln", 0.5)
	print(str(maxGrpNo)+"  sampleCons created:"+str(sampleCons))

	for i in range(1,maxGrpNo+1):
		writeFile = open(filepath+"/WTgroup"+str(i)+".txt", 'w')
		writeFile.write(">sampleCons\n"+sampleCons+"\n")
		writeFile.write(">groupCons\n"+str(getCons(filepath+"/aligned"+str(i)+".aln", grpThreshold)[50:])+"\n")
		groupRaw = open(filepath+"/group"+str(i)+".txt")

		for line in groupRaw:
			if (line[0] == '>'):
				writeFile.write(str(line)+"\n")
			else:
				writeFile.write(str(line[50:])+"\n")
		writeFile.close()
		cmdCall = "mafft --retree 2 --quiet --maxiterate 10 --op 0.05 --ep 0.05 "+filepath+"/WTgroup"+str(i)+".txt > "+filepath+"/WTgroup"+str(i)+".aln"
		aligning = os.system(cmdCall)

	return sampleCons

def getCons(filepath, cutoff):
	consRaw = AlignIO.read(filepath, "fasta")
	consObj = AlignInfo.SummaryInfo(consRaw)
	consFull = str(consObj.gap_consensus(threshold=cutoff, ambiguous="N"))
	print("CONS:  "+str(filepath)+"  "+str(cutoff))
	return consFull


if __name__=='__main__':
	main()

