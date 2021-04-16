import csv
import sys

laneName = sys.argv[1]
dataLocation = sys.argv[2]

outputMatrix = []
midListRaw = ['1','3','4','5','6','7','10','13','14','15','16','17','18']
midList = []
for mid in midListRaw:
	mid = laneName+"-"+mid
	print(mid)
	midList.append(str(mid))
print(midList)
outputList = list(range(101))
outputList[0] = str(laneName)

outputMatrix.append(outputList)

for mid in midList:
	try:
		outputList = [0]*101		
		outputList[0] = str(mid)
		lineList = [line.rstrip('\n') for line in open(str(dataLocation)+'/'+str(mid)+'groupSize.csv')] 
		sumAll = 0
		print(len(lineList))
		for item in lineList[1:]:
			tempList = item.rsplit(',')
			if (int(tempList[0]) < 100):
				outputList[int(tempList[0])] = str(tempList[1])
			else:
				outputList[100] = int(outputList[100])+int(tempList[1])

		outputMatrix.append(outputList)
	except:
		print("mid not found "+str(mid))

with open(dataLocation+"/groupSizesSummary-"+str(laneName)+".csv", "w") as f:
	writer = csv.writer(f,delimiter=',',lineterminator='\n')
	writer.writerows(outputMatrix)
