#!/usr/bin/python

import sys
import os
from Bio import SeqIO, AlignIO
import shutil
from Bio.Align import AlignInfo
import csv
import multiprocessing

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

	midList = [("ACGAGTGCGT",midPrefix+"1"), ("ACGCTCGACA",midPrefix+"2"), ("AGACGCACTC",midPrefix+"3"), ("AGCACTGTAG",midPrefix+"4"), ("ATCAGACACG",midPrefix+"5"), ("ATATCGCGAG",midPrefix+"6"), ("CGTGTCTCTA",midPrefix+"7"), ("CTCGCGTGTC",midPrefix+"8"), ("TCTCTATGCG",midPrefix+"10"), ("TGATACGTCT",midPrefix+"11"), ("CATAGTAGTG",midPrefix+"13"), ("CGAGAGATAC",midPrefix+"14"), ("ATACGACGTA",midPrefix+"15"), ("TCACGTACTA",midPrefix+"16"), ("CGTCTAGTAC",midPrefix+"17"), ("TCTACGTAGC",midPrefix+"18"), ("TGTACTACTC",midPrefix+"19"), ("ACGACTACAG",midPrefix+"20"), ("CGTAGACTAG",midPrefix+"21"), ("TACGAGTATG",midPrefix+"22"), ("TACTCTCGTG",midPrefix+"23"), ("TAGAGACGAG",midPrefix+"24"), ("TCGTCGCTCG",midPrefix+"25"), ("ACATACGCGT",midPrefix+"26"), ("ACGCGAGTAT",midPrefix+"27"), ("ACTACTATGT",midPrefix+"28"), ("ACTGTACAGT",midPrefix+"29"), ("AGACTATACT",midPrefix+"30"), ("AGCGTCGTCT",midPrefix+"31"), ("AGTACGCTAT",midPrefix+"32"), ("ATAGAGTACT",midPrefix+"33"), ("CACGCTACGT",midPrefix+"34"), ("CAGTAGACGT",midPrefix+"35"), ("CGACGTGACT",midPrefix+"36"), ("TACACACACT",midPrefix+"37"), ("TACACGTGAT",midPrefix+"38"), ("TACAGATCGT",midPrefix+"39"), ("TACGCTGTCT",midPrefix+"40"), ("TAGTGTAGAT",midPrefix+"41"), ("TCGATCACGT",midPrefix+"42"), ("TCGCACTAGT",midPrefix+"43"), ("TCTAGCGACT",midPrefix+"44"), ("TCTATACTAT",midPrefix+"45")]

	try:
		for item in midList:                    			#delete anything in the folder from previous analyses
			shutil.rmtree(filepathSplit+"/"+item[1])
			os.remove(filepathSplit+"/"+item[1]+"groupSize.csv")
		os.remove(filepathSplit+"/mutated.txt")
		print("old analyses data cleaned")
	except:
		print("nothing to clean")

	parseByMid(filepathSplit, filepath, midList)			#for each sequence in lane parse into mid group, everything else is in mutated.txt


	jobs = []
	for mid in midList:
		p = multiprocessing.Process(target=speedUpWrapper, args=(filepathSplit, filepath, refSeqPath, mid))
		jobs.append(p)
		p.start()
	print(str(jobs))

	print('all fn complete')


def speedUpWrapper(filepathSplit, filepath, refSeqPath, mid):
	try:
		filepathFull = filepathSplit+"/"+mid[1]
		noGrps = mkGrp(filepathFull, mid)      						#from each mid make subgroups and keep track of how many created
		print(noGrps)
		grp = 0

		for grp in range(1, noGrps+1):
			cmdCall = "mafft --6merpair --retree 2 --quiet --maxiterate 10 --op 0.05 --ep 0.05 "+filepathFull+"/group"+str(grp)+".txt > "+filepathFull+"/aligned"+str(grp)+".aln"
			returnOut = os.system(cmdCall)
			print("\r MID "+mid[1]+" aligning:  {0:.2f} %".format(grp/noGrps*100), end='')#{0:.2f}".format(noSeq/fileLength*100)+" %", end='\r')

		noSubCons = mkCons(filepathFull, noGrps)
		print("SUBCONS FORMED    "+str(noSubCons))

		WTediting = defaultWT(filepathFull, noGrps, 0.75)

	except:
		print("MID not found")

	print("exit: "+str(mid))


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
	cmdCall = "mafft --retree 1 --quiet --maxiterate 0 "+filepath+"/Consensus.txt > "+filepath+"/Consensus.aln"
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
		cmdCall = "mafft --6merpair --retree 2 --quiet --maxiterate 10 --op 0.05 --ep 0.05 "+filepath+"/WTgroup"+str(i)+".txt > "+filepath+"/WTgroup"+str(i)+".aln"
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
