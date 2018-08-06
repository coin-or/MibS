import sys
import math

setSize = int(sys.argv[1])

methodType = sys.argv[2]

path = sys.argv[3]

name = path + "/finalInstanceList"
fileAcceptedInst = open(name, "r")

name = path + "/r1GeR2ListMiblpXu"
fileGeXu = open(name, "r")

name = path + "/r1LeqR2ListIblpDen"
fileLeqDen = open(name, "r")

name = path + "/r1GeR2ListIblpDen"
fileGeDen = open(name, "r")

name = path + "/r1LeqR2ListIblpFis"
fileLeqFis = open(name, "r")

name = path + "/r1GeR2ListIblpFis"
fileGeFis = open(name, "r")

fileTime = open("time.summary", "r")
fileInst = open("instanceList.summary", "r")
fileVFNum = open("VFNum.summary", "r")
fileVFTime = open("VFTime.summary", "r")
fileUBNum = open("UBNum.summary", "r")
fileUBTime = open("UBTime.summary", "r")

name = "tmpOutputMethod" + methodType
fileTmpOut = open(name, "w+")

timeList = fileTime.read().split()

instList = fileInst.read().split()

VFNumList = fileVFNum.read().split()

VFTimeList = fileVFTime.read().split()

UBNumList = fileUBNum.read().split()

UBTimeList = fileUBTime.read().split()

acceptedInstList = fileAcceptedInst.read().split()

geXuList = fileGeXu.read().split()

leqDenList = fileLeqDen.read().split()
geDenList = fileGeDen.read().split()

leqFisList = fileLeqFis.read().split()
geFisList = fileGeFis.read().split()

totalTime = {}
VFNum = {}
VFTime = {}
UBNum = {}
UBTime = {}
timePercent = {}

dataSet = {"geXu", "leqDen", "geDen", "leqFis", "geFis"}
for i in dataSet:
    totalTime[i] = 0.0
    VFNum[i] = 0
    VFTime[i] = 0.0
    UBNum[i] = 0
    UBTime[i] = 0.0

for i in range(setSize):
    instance = instList[i]
    if instance in acceptedInstList:
        if instance in geXuList:
            data = "geXu"
        elif instance in leqDenList:
            data = "leqDen"
        elif instance in geDenList:
            data = "geDen"
        elif instance in leqFisList:
            data = "leqFis"
        else:
            data = "geFis"

        totalTime[data] = totalTime[data] + abs(float(timeList[i]))
        VFNum[data] = VFNum[data] + int(VFNumList[i])
        VFTime[data] = VFTime[data] + float(VFTimeList[i])
        UBNum[data] = UBNum[data] +int(UBNumList[i])
        UBTime[data] = UBTime[data] + float(UBTimeList[i])

for data in dataSet:
    timePercent[data] = int(100 * ((VFTime[data] + UBTime[data]) /totalTime[data]))

for data in dataSet:
    fileTmpOut.write(data + '\n')
    fileTmpOut.write(str(VFNum[data]) + '\n')
    fileTmpOut.write(str(UBNum[data]) + '\n')
    fileTmpOut.write(str(timePercent[data]) + '\n')
            
fileTime.close()
fileInst.close()
fileVFNum.close()
fileVFTime.close()
fileUBNum.close()
fileUBTime.close()
fileTmpOut.close()






