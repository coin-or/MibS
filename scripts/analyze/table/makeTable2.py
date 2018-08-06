import math

fileMethod4 = open("tmpOutputMethod4", "r")
fileMethod7 = open("tmpOutputMethod7", "r")
fileMethod8 = open("tmpOutputMethod8", "r")
fileMethod9 = open("tmpOutputMethod9", "r")

#gathering data
#methodSet = {"leqDen4", "geDen4", "leqFis4", "geFis4", "geXu4",
#             "leqDen7", "geDen7", "leqFis7", "geFis7", "geXu7",
#             "leqDen8", "geDen8", "leqFis8", "geFis8", "geXu8",
#             "leqDen9", "geDen9", "leqFis9", "geFis9", "geXu9"}
VFNum = {}
UBNum = {}
timePercent = {}


listMethod4 = fileMethod4.read().split()
listMethod7 = fileMethod7.read().split()
listMethod8 = fileMethod8.read().split()
listMethod9 = fileMethod9.read().split()

for i in {4, 7, 8, 9}:
    if i == 4:
        listTmp = listMethod4
    elif i == 7:
        listTmp = listMethod7
    elif i == 8:
        listTmp = listMethod8
    else:
        listTmp = listMethod9
        
    for j in range(20):
        if "leq" in listTmp[j] or "ge" in listTmp[j]:
            methodName = listTmp[j] + str(i)
            VFNum[methodName] = listTmp[j+1]
            UBNum[methodName] = listTmp[j+2]
            timePercent[methodName] = str(int(listTmp[j+3]))

#making table2

outFile = open("table2", "w+")

outFile.write('\n')

outFile.write("Data Set")
outFile.write("         withPoolWhenXYIntOr")
outFile.write("           withoutPoolWhenXYIntOr")
outFile.write("       withPoolWhenXYIntOr")
outFile.write("             withoutPoolWhenXYIntOr" + '\n')

outFile.write("                 LFixed-LFixed(Linking)")
outFile.write("        LFixed-LFixed(Linking)")
outFile.write("       LFixed-LFixed(Fractional)")
outFile.write("       LFixed-LFixed(Fractional)" + '\n')

outFile.write("              ---------------------------")
outFile.write("   ---------------------------")
outFile.write("   ---------------------------")
outFile.write("     -----------------------------")
outFile.write('\n')

outFile.write("                  SL        UB     Time")
outFile.write("       SL        UB     Time")
outFile.write("         SL        UB     Time")
outFile.write("           SL        UB     Time")
outFile.write('\n')

for dataSet in ["Den", "Fis", "Xu"]:
    for sign in ["leq", "ge"]:
        if dataSet != "Xu" or sign == "ge":
            if dataSet == "Den":
                outFile.write("IBLP-DEN")
            elif dataSet == "Fis":
                outFile.write("IBLP-FIS")
            else:
                outFile.write("IBLP-XU")
            outFile.write("        ")

            count = 0
            for i in [4, 8, 7, 9]:
                methodNum = sign + dataSet + str(i)
                value = VFNum[methodNum]
                length = len(value)
                spaceNum = 3 - int(math.ceil(length/2))
                for j in range(spaceNum):
                    outFile.write(" ")
                outFile.write(value)
                spaceNum = 7 - spaceNum - length
                for j in range(spaceNum):
                    outFile.write(" ")
                outFile.write("   ")

                value = UBNum[methodNum]
                length = len(value)
                spaceNum = 3 - int(math.ceil(length/2))
                for j in range(spaceNum):
                    outFile.write(" ")
                outFile.write(value)
                spaceNum = 7 - spaceNum - length
                for j in range(spaceNum):
                    outFile.write(" ")
                outFile.write("   ")

                value = timePercent[methodNum]
                outFile.write(value)
                length = len(value)
                if length == 1:
                    outFile.write(" ")
                spaceNum = 6 + 2 * count
                count = count + 1
                for j in range(spaceNum):
                    outFile.write(" ")

            outFile.write('\n')

            if dataSet != "Xu":
                if sign == "leq":
                    outFile.write("(r1 <= r2)")
                    outFile.write('\n')
                    outFile.write('\n')
                else:
                    outFile.write("(r1 > r2)")
                    outFile.write('\n')
                    outFile.write('\n')

outFile.close()
                
                
    



