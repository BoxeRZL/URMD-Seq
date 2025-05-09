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
import multiprocessing


def main():
    print(str(sys.argv))
    midList = ['1', '2', '3', '4', '5', '6', '7', '8', '10', '11', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45']

    midList = [str(sys.argv[7]) + mid for mid in midList]
    doAdjustment = True

    analysisName = str(sys.argv[6])
    outPath = analysisName + str(date.today().strftime("%d-%b-%Y")) + ".csv"
    print(str(outPath))

    outputFile = open(outPath, 'a')
    outputFile.write(','.join(['groupSize', 'mutationFreq', 'lane', 'position', 'mid', 'reference', 'coverage', 'a', 't', 'c', 'g', 'n', '-', 'nonWTbases']))
    outputFile.write("\n")
    outputFile.close()

    processingVar = str(','.join([str(sys.argv[2]), str(sys.argv[3]), str(sys.argv[5])]))

    # New lines: assign command-line arguments to descriptive variables
    grpSizeThreshold = int(sys.argv[2])
    consThreshold = float(sys.argv[3])
    refSeqPath = str(sys.argv[4])

    jobs = []
    for mid in midList:
        midsRoot = sys.argv[1]
        filepath = midsRoot + "/" + mid
        # Modified line: use descriptive variable names in function call
        p = multiprocessing.Process(target=multiprocessWrapper, args=(midsRoot, filepath, mid, grpSizeThreshold, consThreshold, refSeqPath, outPath, processingVar, analysisName, doAdjustment))
        jobs.append(p)
        p.start()
    print(str(jobs))


    print("all fn complete")



def multiprocessWrapper(midsRoot, filepath, mid, grpSizeThreshold, consThreshold, refSeqPath, outPath, processingVar, analysisName, doAdjustment):
    try:
        unzipCall = "unzip -q " + str(filepath) + "/WTgroupsAligned-" + str(mid) + ".zip -d " + str(filepath) + "/"
        zipReturn = os.system(unzipCall)
        maxGrps = len(glob.glob1(filepath, "WTgroup*.aln"))

    except:
        print("no groups to align")
    try:
        mkConsWTadjusted(filepath, mid, maxGrps, grpSizeThreshold, consThreshold, refSeqPath, outPath, processingVar, analysisName, doAdjustment)
        rmUnzippedFiles = "rm " + str(filepath) + "/WTgroup*.aln"
        rmReturn = os.system(rmUnzippedFiles)
    except:
        print("Unable to complete analysis")

    try:
        rmUnzippedFiles = "rm " + str(filepath) + "/WTgroup*.aln"
        rmReturn = os.system(rmUnzippedFiles)
    except:
        print("Unable to remove zipped files")
        # CCTCGAGGTCGACGGTATCGACGAGTGCGTNNNNCANNNNGTNNNNNgcccacacgttccccttaaataag


def mkConsWTadjusted(filepath, mid, maxGrpNo, grpSizeThreshold, consThreshold, refSeqPath, outPath, processingVar, analysisName, doAdjustment):
    writeFile = open(filepath + "/WTadjConsensus.txt", "w+")
    refSeq = SeqIO.read(refSeqPath, "fasta")
    writeFile.write(">refSeq\n" + str(refSeq.seq) + "\n")
    writeFile.close()

    for i in range(1, maxGrpNo + 1):
        adjWT(filepath, i, grpSizeThreshold, consThreshold, mid, str(refSeq.seq), analysisName, doAdjustment)

    cmdCall = "mafft --retree 1 --quiet " + filepath + "/WTadjConsensus.txt > " + filepath + "/WTadjConsensus.aln"
    aligning = os.system(cmdCall)
    try:
        mutPrint(filepath, mid, "/WTadjConsensus.aln", outPath, processingVar)
    except:
        print("could not find required alignment files")


def adjWT(filepath, groupNo, grpSizeThreshold, consThreshold, mid, refSeqStr, analysisName, doAdjustment):
    grpAlnRaw = AlignIO.read(filepath + "/WTgroup" + str(groupNo) + ".aln", "fasta")
    grpAlnObj = AlignInfo.SummaryInfo(grpAlnRaw)
    WTgroupCons = str(grpAlnObj.gap_consensus(threshold=consThreshold)).upper()  # Apply consThreshold here

    # Replace ambiguous nucleotides 'X' that doesn't pass the consThreshold with 'N'
    WTgroupCons = WTgroupCons.replace("X", "N")
    
    # Debug print statement to check the generated consensus sequence
    print(f"Generated consensus sequence for group {groupNo}: {WTgroupCons}")

    for c in range(0, len(WTgroupCons)):
        if (doAdjustment) and (WTgroupCons[c] == "N"):
            wtCount = grpAlnRaw[2:, c].count(grpAlnRaw[0, c])
            wtPrevalence = wtCount / len(grpAlnRaw[2:, c])
            WTgroupCons = WTgroupCons[:c] + str(grpAlnRaw[0, c]) + WTgroupCons[c + 1:]

    sequenceMismatch = 0
    referencePairwise = str(grpAlnRaw[0].seq).replace("-", "").upper()
    groupPairwise = str(WTgroupCons).replace("-", "").upper()
    sequenceMismatch = int(format_alignment(*pairwise2.align.globalms(referencePairwise, groupPairwise, 1, 0, -1, 0)[0]).count(" "))
    sequenceMismatch = sequenceMismatch + int(format_alignment(*pairwise2.align.globalms(referencePairwise, groupPairwise, 1, 0, -1, 0)[0]).count("."))

    if (sequenceMismatch < 15) and (len(grpAlnRaw[2:, 0]) >= grpSizeThreshold):
        writeFile = open(filepath + "/WTadjConsensus.txt", 'a')
        writeFile.write(">" + str(groupNo) + "\n" + str(WTgroupCons).upper() + "\n")
        writeFile.close()

    if (sequenceMismatch >= 15):
        print(str(sequenceMismatch) + " from group: " + str(groupNo) + " mismatches:" + str(sequenceMismatch))
        print(format_alignment(*pairwise2.align.globalms(referencePairwise, groupPairwise, 1, 0, -1, 0)[0]))
        writefile = open(filepath + "/../../" + str(analysisName) + "-contamination.txt", 'a')
        if os.stat(filepath + "/../../" + str(analysisName) + "-contamination.txt").st_size == 0:
            writefile.write(">refSeq\n" + str(refSeqStr) + "\n")
        writefile.write(">" + str(mid) + " group " + str(groupNo) + " mismatches:" + str(sequenceMismatch) + "\n" + str(WTgroupCons).upper() + "\n")
        writefile.close()

    return 0



def mutPrint(filepath, mid, rawFile, outPath, processingVar):
    posn = 16533  # refseq before start site -1 for zero indexing in python
    chrSize = 16568  # mtDNA size (16568)-1 for zero indexing
    insert = 0

    refAlignment = AlignIO.read(filepath + rawFile, "fasta")

    # Generate unique filename for multiprocessing-safe writing
    pid = os.getpid() # Process ID: a unique ID that prevent overlapping leading to garbled or misaligned rows
    outFilePath = f"{outPath}.worker_{pid}.csv"
    outputFile = open(outFilePath, 'a')

    for base in range(0, len(refAlignment[0])):
        if refAlignment[0][base] == "-":
            insert += 1
        else:
            posn += 1
            insert = 0

        a = refAlignment[1:, base].count("a")
        t = refAlignment[1:, base].count("t")
        c = refAlignment[1:, base].count("c")
        g = refAlignment[1:, base].count("g")
        n = refAlignment[1:, base].count("n")
        indel = refAlignment[1:, base].count("-")

        coverage = a + t + c + g + n + indel

        if ((a / coverage == 1) or (t / coverage == 1) or (c / coverage == 1) or (g / coverage == 1) or (n / coverage == 1) or (indel / coverage == 1)):
            continue

        if posn == 0:
            posn = int(chrSize) + 1

        wtBases = max(a, t, c, g)

        outputFile.write(",".join([
            str(processingVar),
            str(((posn + (0.001 * insert)) % chrSize)),
            str(mid),
            str(refAlignment[0][base]),
            str(coverage),
            str(a), str(t), str(c), str(g), str(n), str(indel),
            str(coverage - n - indel - wtBases)
        ]))

        if insert == 0:
            outputFile.write(",1")
        else:
            outputFile.write(",")

        outputFile.write("\n")

    outputFile.close()


    print("# groups:    " + str(len(refAlignment[1:, 1])))



if __name__ == '__main__':
    main()
