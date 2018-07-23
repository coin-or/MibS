import sys

setSize = int(sys.argv[1])

shouldRefineInstList = 0
if len(sys.argv) == 3:
  shouldRefineInstList = 1
  instanceFile = open(sys.argv[2],"r")
  finalInstanceFile = open("finalInstanceList","w+")
  

methodNum = 11

file1 = open("timeMethod1.summary","r")

file2 =open("timeMethod2.summary","r")

file3 =open("timeMethod3.summary","r")

file4 =open("timeMethod4.summary","r")

file5 =open("timeMethod5.summary","r")

file6 =open("timeMethod6.summary","r")

file7 =open("timeMethod7.summary","r")

file8 =open("timeMethod8.summary","r")

file9 =open("timeMethod9.summary","r")

file10 =open("timeMethod10.summary","r")

file11 =open("timeMethod11.summary","r")


finalFile1 = open("finalTimeMethod1","w+")

finalFile2 =open("finalTimeMethod2","w+")

finalFile3 =open("finalTimeMethod3","w+")

finalFile4 =open("finalTimeMethod4","w+")

finalFile5 =open("finalTimeMethod5","w+")

finalFile6 =open("finalTimeMethod6","w+")

finalFile7 =open("finalTimeMethod7","w+")

finalFile8 =open("finalTimeMethod8","w+")

finalFile9 =open("finalTimeMethod9","w+")

finalFile10 =open("finalTimeMethod10","w+")

finalFile11 =open("finalTimeMethod11","w+")


mat = [[0 for x in range(methodNum)] for y in range(setSize)]

mat[0] = file1.read().split()

mat[1] = file2.read().split()

mat[2] = file3.read().split()

mat[3] = file4.read().split()

mat[4] = file5.read().split()

mat[5] = file6.read().split()

mat[6] = file7.read().split()

mat[7] = file8.read().split()

mat[8] = file9.read().split()

mat[9] = file10.read().split()

mat[10] = file11.read().split()

if shouldRefineInstList == 1:
  indList = instanceFile.read().split()
  


accepted = []
cnt = 0
for j in range(setSize):
  for i in range(methodNum):
    if float(mat[i][j]) >= 5.00:
      accepted.append(j)
      cnt += 1
      break
    elif float(mat[i][j]) < 0:
      for k in range(methodNum):
        if float(mat[k][j]) >= 0:
          accepted.append(j)
          cnt += 1
          break
      break


for i in range(cnt):
  index = accepted[i]
  finalFile1.write(mat[0][index] + '\n')
  finalFile2.write(mat[1][index] + '\n')
  finalFile3.write(mat[2][index] + '\n')
  finalFile4.write(mat[3][index] + '\n')
  finalFile5.write(mat[4][index] + '\n')
  finalFile6.write(mat[5][index] + '\n')
  finalFile7.write(mat[6][index] + '\n')
  finalFile8.write(mat[7][index] + '\n')
  finalFile9.write(mat[8][index] + '\n')
  finalFile10.write(mat[9][index] + '\n')
  finalFile11.write(mat[10][index] + '\n')
  if shouldRefineInstList == 1:
    finalInstanceFile.write(indList[index] + '\n')

file1.close()
file2.close()
file3.close()
file4.close()
file5.close()
file6.close()
file7.close()
file8.close()
file9.close()
file10.close()
file11.close()

finalFile1.close()
finalFile2.close()
finalFile3.close()
finalFile4.close()
finalFile5.close()
finalFile6.close()
finalFile7.close()
finalFile8.close()
finalFile9.close()
finalFile10.close()
finalFile11.close()

if shouldRefineInstList == 1:
  instanceFile.close()
  finalInstanceFile.close()

