import sys
import csv

tableName = sys.argv[1]

addressMain = sys.argv[2]
address = {}

address["MiblpXu"] = addressMain + "/MIBLP-XU/"
address["IblpDen"] = addressMain + "/RANDOM/RAND_BILEVEL/"
address["IblpFis"] = addressMain + "/IBLP-FIS/"

if tableName == "table3.csv":
    mainColNum = 5
    dataSet = ["IblpDen", "IblpFis"]
    methodList = ["1", "2", "3", "4", "5"]
elif tableName == "table4.csv":
    mainColNum = 5
    dataSet = ["MiblpXu"]
    methodList = ["1", "2", "3", "4", "5"]
elif tableName == "table5.csv":
    mainColNum = 2
    dataSet = ["IblpDen", "IblpFis"]
    methodList = ["4", "7"]
elif tableName == "table6.csv":
    mainColNum = 2
    dataSet = ["MiblpXu", "IblpDen", "IblpFis"]
    methodList = ["4", "7"]
elif tableName == "table7.csv":
    mainColNum = 4
    dataSet = ["MiblpXu", "IblpDen", "IblpFis"]
    methodList = ["11", "6", "9", "7"]
elif tableName == "table8.csv":
    mainColNum = 4
    dataSet = ["MiblpXu", "IblpDen", "IblpFis"]
    methodList = ["10", "3", "8", "4"]
elif tableName == "table9.csv":
    mainColNum = 4
    dataSet = ["MiblpXu", "IblpDen", "IblpFis"]
    methodList = ["15", "12", "13", "14"]

#open input files
fileInInstance = []
fileInULVar = []
fileInLLVar = []
fileInBestSol = []
fileInTime = []
fileInNode = []
instanceDict = {}
bestSolDict = {}
timeDict = {}
nodeDict = {}
uLVarDict = {}
lLVarDict = {}

count = 0
for i in dataSet:
    name = "instanceList" + i + ".summary"
    fileInInstance.append(open(name, "r"))
    instanceDict[i] = fileInInstance[count].read().split()
    count = count + 1

count = 0
if tableName in ["table5.csv", "table6.csv"]:
    for i in dataSet:
        name = "uLIntVarNum" + i + ".summary"
        fileInULVar.append(open(name, "r"))
        uLVarDict[i] = fileInULVar[count].read().split()

        name = "lLIntVarNum" + i + ".summary"
        fileInLLVar.append(open(name, "r"))
        lLVarDict[i] = fileInLLVar[count].read().split()

        count = count + 1
    


count = 0
methodNum = 1
for i in methodList:
    for j in dataSet:
        name = "cost" + j + "Method" + i + ".summary"
        fileInBestSol.append(open(name, "r"))
        name = j + str(methodNum)
        bestSolDict[name] = fileInBestSol[count].read().split()
        
        name = "time" + j + "Method" + i + ".summary"
        fileInTime.append(open(name, "r"))
        name = j + str(methodNum)
        timeDict[name] = fileInTime[count].read().split()

        name = "node" + j + "Method" + i + ".summary"
        fileInNode.append(open(name, "r"))
        name = j + str(methodNum)
        nodeDict[name] = fileInNode[count].read().split()

        count = count + 1

    methodNum = methodNum + 1

#make field name
if tableName in ["table5.csv", "table6.csv"]:
    fieldNames = ["Instance", "r1" , "r2"]
else:
    fieldNames = ["Instance"]

for i in range(mainColNum):
    newCol = "BestSol" + str(i + 1) 
    fieldNames.append(newCol)

    newCol = "Time" + str(i + 1)
    fieldNames.append(newCol)

    newCol = "Nodes" + str(i + 1)
    fieldNames.append(newCol)

#write output
with open(tableName, 'w') as fileOut:
    writer = csv.DictWriter(fileOut, fieldnames=fieldNames)

    writer.writeheader()
    
    #if tableName not in ["table5.csv", "table6.csv"]:
    for i in dataSet:
        instanceList = instanceDict[i]
        for j in range(len(instanceList)):
            newRow = {}
            instance = instanceList[j].replace(address[i], "")
            newRow['Instance'] = instance
            if tableName in ["table5.csv", "table6.csv"]:
                uLVarList = uLVarDict[i] 
                newRow['r1'] = uLVarList[j]

                lLVarList = lLVarDict[i]
                newRow['r2'] = lLVarList[j]
                
            for k in range(mainColNum):
                name = i + str(k + 1)
                bestSolList = bestSolDict[name]
                timeList = timeDict[name]
                nodeList = nodeDict[name]

                newCol = "BestSol" + str(k + 1)
                newRow[newCol] = bestSolList[j]
                newCol = "Time" + str(k + 1)
                if float(timeList[j]) < 0:
                    time = " > 3600"
                    node = str(abs(int(nodeList[j])))
                else:
                    time = timeList[j]
                    node = nodeList[j]
                newRow[newCol] = time
                newCol = "Nodes" + str(k + 1)
                newRow[newCol] = node

            writer.writerow(newRow)

for i in range(len(dataSet)):
    fileInInstance[i].close()
    fileInBestSol[i].close()
    fileInTime[i].close()
    fileInNode[i].close()
    if tableName in ["table5.csv", "table6.csv"]:
        fileInULVar[i].close()
        fileInLLVar[i].close()
                
            
                
            
    
