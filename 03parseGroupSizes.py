import csv
import sys

laneName = sys.argv[1]
dataLocation = sys.argv[2]

outputMatrix = []
midListRaw = ['1','2','3','4','5','6','7','8','10','11','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45']
midList = []
for mid in midListRaw:
	#mid = dataLocation+"/"+laneName+"-"+mid
	mid = laneName+"-"+mid
	print(mid)
	midList.append(str(mid))
print(midList)
outputList = list(range(1001))
outputList[0] = str(laneName)

outputMatrix.append(outputList)

for mid in midList:
	try:
		outputList = [0]*1001		
		outputList[0] = str(mid)
		lineList = [line.rstrip('\n') for line in open(str(dataLocation)+'/'+str(mid)+'groupSize.csv')] 
		sumAll = 0
		print(len(lineList))
		for item in lineList[1:]:
			tempList = item.rsplit(',')
			#print(tempList[0])
			#print(str(tempList))
			if (int(tempList[0]) < 1000):
				outputList[int(tempList[0])] = str(tempList[1])
			else:
				outputList[1000] = int(outputList[1000])+int(tempList[1])
			#sumAll = sumAll+int(tempList[1])
			#print(sumAll)
		#outputList[202] = sumAll-sum(outputList[1:201])
		outputMatrix.append(outputList)
	except:
		print("mid not found "+str(mid))

with open(dataLocation+"/groupSizesSummary-"+str(laneName)+".csv", "w+") as f:
	writer = csv.writer(f,delimiter=',',lineterminator='\n')
	writer.writerows(outputMatrix)
