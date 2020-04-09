#!/usr/bin/python


import sys
import os
from Bio import SeqIO, AlignIO
import shutil
from Bio.Align import AlignInfo
import csv
import multiprocessing

####USAGE --  grpParseClean.py /folder/Location rawTextInFolder refSeq multiprocessFlag
##folder location should only contain raw fasta for analysis

def main():
	filepathSplit = sys.argv[1]								#folderLocation
	filepath = filepathSplit+"/"+str(sys.argv[2])			#raw data
	print(str(filepath))
	refSeqPath = sys.argv[3]

	try:													#open raw data file
		lineNo = 0
		fp = open(filepath, "r")
		print("file opened: "+str(filepath))
	except:
		print("ERROR: invalid filepath to raw fasta")

	#midList = [("ACGAGTGCGT","1"), ("AGACGCACTC","3"), ("AGCACTGTAG","4"), ("ATCAGACACG","5"), ("ATATCGCGAG","6"), ("CGTGTCTCTA", "7"), ("TCTCTATGCG","10"), ("CATAGTAGTG","13"), ("CGAGAGATAC", "14"), ("ATACGACGTA","15"), ("TCACGTACTA","16"), ("CGTCTAGTAC","17"), ("TCTACGTAGC","18"), ("TGTACTACTC","19"), ("ACGACTACAG","20"),("CGTAGACTAG","21"), ("TACGAGTATG","22"), ("TAGAGACGAG","24"), ("AGACTATACT","30"), ("ATAGAGTACT","33"), ("CGACGTGACT","36"), ("TACAGATCGT","39"), ("TACGCTGTCT","40")]
	midList = [("ATACGACGTA","15"), ("ACGAGTGCGT","1")]

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
		#speedUpWrapper(filepathSplit, filepath, refSeqPath, mid)
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

		#zipGroups = "zip -mq "+str(filepathFull)+"/rawGroups.zip "+str(filepathFull)+"/group*.txt"
		#print(zipGroups)
		#zipReturn = os.system(zipGroups)
		#zipAligned = "zip -mq "+str(filepathFull)+"/alignedGroups.zip "+str(filepathFull)+"/aligned*.aln"
		#print(zipAligned)
		#zipReturn = os.system(zipAligned)

	except:
		print("MID not found")
#
#		noGrps = mkGrpCons(noGrps, filepathFull, refSeqPath)	
#		print("MKGRPCONS DONE, CONSGROUP WRITTEN   "+str(noGrps))

#		for grp in range(1, noGrps+1):
#			cmdCall = "mafft --retree 2 --quiet --maxiterate 10 --op 0.05 --ep 0.05 "+filepathFull+"/ConsGroup"+str(grp)+".txt > "+filepathFull+"/ConsAligned"+str(grp)+".aln"
#			returnOut = os.system(cmdCall)
#			print("MID "+mid[1]+" consensus aligning:  {0:.2f} %".format(grp/noGrps*100), end='\r')#{0:.2f}".format(noSeq/fileLength*100)+" %", end='\r')


	print("exit: "+str(mid))


def mkGrpCons(noGrps, filepathFull, refSeqPath):

	refSeq = SeqIO.read(refSeqPath, "fasta")
	print("mkGrpCons  " +str(filepathFull)+"  "+str(refSeqPath))

	supAlnRaw = AlignIO.read(filepathFull+"/Consensus.aln", "fasta")
	supAlnObj = AlignInfo.SummaryInfo(supAlnRaw)
	#print(str(supAlnObj))
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
	#print("in parseByMid"+ str(filepath))
	#fileLength = (sum(1 for line in open(filepath)))/2
	fileLength = len(list(SeqIO.parse(filepath, "fastq")))
	for seqFragment in SeqIO.parse(filepath, "fastq"):				#for each sequence in lane parse into mid group, everything else is in mutated.txt
		print("\r searching for MID:  {0:.2f} %".format(noSeq/fileLength*100), end='')
		noSeq += 1
		found = 0
		#print("\rMID searching:  {0:.2f} %".format(noSeq/fileLength*100))#, end='\r')#(noSeq/fileLength*100," percent complete         \r")
		for mid in midList:
			posn = seqFragment.seq.find(mid[0])     				#search for mid in seq
			if (posn >= 0):                         				#if mid sequence is found
				if ((str(seqFragment.seq)[posn+14:posn+16] == "CA") and (str(seqFragment.seq)[posn+20:posn+22] == "GT") and (str(seqFragment.seq)[posn+27:posn+29] == "GC")):# and (len(seqFragment.seq) > 300)):		#check for presence of conserved sites in start seq
					#	***how to print status non-permanently (update as each mid is found)
					#print(str(mid[1])+" FOUND  "+str(posn)+" mid:  "+str(mid[0])+"       next conserved:"+str(seqFragment.seq)[posn+14:posn+16]+str(seqFragment.seq)[posn+20:posn+22]+str(seqFragment.seq)[posn+27:posn+29])
					seqFragment.description = str(mid[0])
					found = str(mid[1])

		#print("id: "+str(found))
		if ((found) and (len(seqFragment.seq) >= 300)):         										#if matching MID is found append seq to 'mid/ungrouped.txt'
			writePath = filepathSplit+"/"+found+"/ungrouped.txt" 
			#print("writing ungrouped")         
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
		print("\r grouping:  {0:.2f} %".format(noSeq/fileLength*100), end='')
		posn = seqFragment.seq.find(mid[0])
		seqStart = str(seqFragment.seq)[posn:posn+29]		#seqStart is the first 29 bases (the PID)
		
		#print("\rsearching:  "+filepath+"/ungrouped.txt   for: "+seqStart)

		if seqStart in grpList.keys():
			grpFound = True
		else:
			grpFound = False

		if (grpFound):
			#grpList.append(seqStart)
			grpList[seqStart][1] += 1
			grpList[seqStart].append(seqFragment)
		else:
			grpNo += 1
			grpList[seqStart] = [grpNo, 1, seqFragment]
		noSeq += 1

	#print("grouping:  {0:.2f} %".format(noSeq/fileLength*100))
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
		#print("group done: "+str(grpList[sortedGrp][0]))


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
	for i in range(1,maxGrpNo):
		alnRaw = AlignIO.read(filepath+"/aligned"+str(i)+".aln", "fasta")
		alnObj = AlignInfo.SummaryInfo(alnRaw)
		cons = alnObj.gap_consensus(threshold=0.51,ambiguous="-")
		writeFile.write(">"+str(i)+"\n"+str(cons).upper()+"\n")
	
	writeFile.close()
	cmdCall = "mafft --retree 2 --quiet --maxiterate 10 --op 0.05 --ep 0.05 "+filepath+"/Consensus.txt > "+filepath+"/Consensus.aln"
	print("complete consensus calc: "+str(cmdCall))
	returnOut = os.system(cmdCall)


if __name__=='__main__':
	main()


##ISSUES TO CHANGE IF NEEDED
##mkGrp returns number of all groups or only groups > 1 in size
##does parseByMid need to reurn a value of how many seq used

