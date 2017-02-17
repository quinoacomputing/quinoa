#! /usr/bin/env python


""" Need a doc string """

import sys
import os


nargin = len(sys.argv);


showUsage = 0
if nargin < 4 or nargin > 5:
   showUsage = 1
else:
   argv1 = sys.argv[1].lower()
   if argv1.rfind('help') > -1:
      showUsage = 1

if showUsage == 1:
   a = 'NOX-executable '
   b = 'parameter-file '
   c = 'NOX-output-file '
   d = 'expected-NOX-output-file(optional)'
   os.system('echo ')
   os.system('echo Usage:')
   os.system('echo     ' +sys.argv[0]+' '+a+b+c+d)

if nargin > 1:
   NOXexecutable = sys.argv[1]
if nargin > 2:
   pythonInputFile = sys.argv[2]

os.system("rm -f hostdata;hostname > hostdata;")
hostdata = open("hostdata",'r')
hostname = hostdata.read().rstrip()
hostdata.close()
os.system("rm -f hostdata")

if nargin > 3:
   NOXoutput = sys.argv[3]
if NOXoutput.find(hostname) == -1:
   NOXoutput += '_'+hostname

os.system('rm -f '+NOXoutput)
#os.execv(rm,[-f,NOXoutput])

#
# if NOXexpectedOutput is provided, generate single files from it,
# one for each run
#
if nargin == 5:
   NOXexpectedOutput = sys.argv[4]
   os.system('rm -f diffOutput_'+hostname)

   try:
      NOXexpectedOutputFile = open(NOXexpectedOutput,'r')
   except IOError:
      try:
         NOXexpectedOutputFile = open(NOXexpectedOutput+'_'+hostname,'r')
      except:
         os.system("echo ERROR: Cannot open file "+NOXexpectedOutput)


   NOXexpectedOutputLines = NOXexpectedOutputFile.readlines()
   if NOXexpectedOutputLines[0].rfind('#NOX TEST TAG') == -1:
      os.system('echo ERROR: #NOX TEST TAG 1# not found in first line of '+ NOXexpectedOutput)
   expectedTestNumberList = []
   for line in NOXexpectedOutputLines:
      if line.rfind('#NOX TEST TAG') > -1:
         splitLine = line.split()
         expectedTestNumber = int(splitLine[3][0:-1])
         expectedTestNumberList.append(expectedTestNumber)
         if expectedTestNumber > 1:
            os.fsync(tmpNOXexpectedOutputFile.fileno())
            tmpNOXexpectedOutputFile.close()
         tmpNOXexpectedOutputFileName = 'tmpExpectedOutput'+splitLine[3][0:-1]
         tmpNOXexpectedOutputFile = open(tmpNOXexpectedOutputFileName,'w')
      else:
         tmpNOXexpectedOutputFile.write(line)
   os.fsync(tmpNOXexpectedOutputFile.fileno())
   tmpNOXexpectedOutputFile.close()

# read parameter input file, and place the lines in paramInputLines
try:
   input = open(pythonInputFile,'r')
except IOError:
   os.system('echo ERROR: Cannot open file '+pythonInputFile)

paramInputLines = input.readlines()
input.close

#
# remove blank lines and comments from the paramInputLines data structure
#
phraseID = 0
deleteList = []
for phrase in paramInputLines:
   # tag comment and blank lines... these will be deleted below
   if phrase.isspace() or phrase[0] == '#':
       deleteList.append(phraseID)
   phraseID = phraseID+1

# delete tagged lines
deleteList.reverse()
for deleteItem in deleteList:
   del paramInputLines[deleteItem];

#
# Identify Lines indicating multiple sublist parameters via "VAR x"
#
VARlist = []      # contains all VAR values
VARindexList = [] # contains the line numbers containg respective VAR values
numCombinations = 1
for phrase in paramInputLines:
   # find lines containing VAR'
   truncAt = phrase.rfind("VAR")

   if truncAt > -1:
      phraseLoc = paramInputLines.index(phrase)
      VARindexList.append(phraseLoc)
      splitPhrase = phrase.split()
      indx = splitPhrase.index('VAR')
      VARval = int(splitPhrase[indx+1])
      numCombinations *= VARval
      VARlist.append(VARval)

      paramInputLines[phraseLoc] = phrase[0:truncAt]

if nargin == 5:
   if numCombinations != len(expectedTestNumberList):
      os.system('echoERROR: The number of tagged runs in "'+NOXexpectedOutput+'" does not match the number')
      os.system('echo ERROR: of combinations of parameters in "'+pythonInputFile+'".')
      os.system('echo ERROR: \t\tNumber of tagged runs in "'+NOXexpectedOutput+'" is ',len(expectedTestNumberList))
      os.system('echo ERROR: \t\tNumber of combinations in "'+pythonInputFile+'" is ',numCombinations)
      os.system('echo ""')

if numCombinations > 1:
   os.system('echo Starting suite of tests for '+os.getcwd()+NOXexecutable.replace('.',''))
else:
   os.system('echo Starting test for '+os.getcwd()+NOXexecutable.replace('.',''))
   VARlist = [1]
   VARindexList = [0]

lenVARlist = len(VARlist)

#
# Check if enough sublist parameter values are provided... i.e. do these numbers
# match the corresponding VAR values?
#
for i in range(lenVARlist):
   for j in range(1,VARlist[i]+1):
      indx = VARindexList[i] + j
      if paramInputLines[indx][0:2] == '@@':
         os.system('echo ERROR: Too few sublist parameters placed in '+ pythonInputFile)
         os.system('echo ERROR: See line: '+paramInputLines[VARindexList[i]]+' VAR '+str(VARlist[i]))


##
## Loop through all possible combinations of sublists parameters. x contains
## the given combination. for each combination a temporary nox parameter input
## file is generated and the nox executable is called. output is placed in
## temporary files to be diffed against the expected output. the output is
## also accumulated in NOXoutput for saving.
##
x = []
for q in range(lenVARlist):
   x.append(0)

testNumber = 1
testNumberList = []
statusList = []
entireParamList = []
while x[lenVARlist-1] < VARlist[lenVARlist-1]:

   currParamList = []
   doNotPrint = [];
   for q in range(len(VARlist)):
      for r in range(VARlist[q]):
         if r != x[q]:
            doNotPrint.append(r+VARindexList[q]+1)
         else:
            currParamList.append(paramInputLines[r+VARindexList[q]+1].rstrip())

   entireParamList.append(currParamList)

   #print '\nBeginning test '+str(testNumber)+'...'

   #
   # write temporary nox parameter input files using combinations of information
   # within pythonInputFile.
   #
   testNumberList.append(str(testNumber))
   tmpNOXInput =  open('./tmpNOXInputFile','w');

   numItemsOnLine = 0
   for y in range(len(paramInputLines)):
      if doNotPrint.count(y) == 0:
         if paramInputLines[y][0:2] == '@@':
            if numItemsOnLine == 1:
               tmpNOXInput.write('\n')
            tmpNOXInput.write(paramInputLines[y][0:2].rstrip()+'\n')
            numItemsOnLine = 0
         elif paramInputLines[y][0] == "@" and paramInputLines[y][1] != "@":
            tmpNOXInput.write('\n' + paramInputLines[y].rstrip()+'\n')
            numItemsOnLine = 0
         elif numItemsOnLine == 1:
            tmpNOXInput.write(' '+paramInputLines[y].rstrip()+'\n')
            numItemsOnLine = 0
         else:
            tmpNOXInput.write(paramInputLines[y].rstrip())
            numItemsOnLine = numItemsOnLine + 1

   os.fsync(tmpNOXInput.fileno())
   tmpNOXInput.close()

   # execute nox using the temporary parameter input files
   command = 'echo \"#NOX TEST TAG '+str(testNumber)+'#\" >> ' + NOXoutput
   os.system(command)
   command = NOXexecutable+' -p ./tmpNOXInputFile  >>  tmp'+NOXoutput+str(testNumber)
   # print 'Executing ' + command
   status = os.system(command)
   statusList.append(status)
   # print 'Completed run with status ',status
   os.system('cat tmp'+NOXoutput+str(testNumber)+' >> '+NOXoutput)
   os.system('rm -f tmpNOXInputFile');
   testNumber = testNumber+1

   # update the combination of parameters
   x[0] = x[0]+1
   i = 0
   while x[i] == VARlist[i] and i < lenVARlist-1:
      x[i] = 0
      x[i+1] = x[i+1]+1
      i = i+1

# print '\nDone running tests\n'

##
## diff the temporary nox output files against the temporary expected nox
## output files
##
if nargin == 5:
   os.system("sync")
   diffList = []
   for i in testNumberList:
      old = 'tmpExpectedOutput'+str(i)
      new = 'tmp'+NOXoutput+str(i)

      command = 'echo \"#NOX TEST TAG '+str(i)+'#\" >> diffOutput_'+hostname
      os.system(command)
      command = 'diff '+ new + ' ' + old + '>> diffOutput_'+hostname+' &> /dev/null'
      diffSuccess = os.system(command)
      diffList.append(diffSuccess)
      #os.system('rm -f '+old)

   # remove the temporary expected nox output files
   for i in expectedTestNumberList:
      os.system('rm -f '+ 'tmpExpectedOutput'+str(i))

statusTestSuccessful = 1
# remove the temporary nox output files
for i in testNumberList:
   if statusList[int(i)-1] != 0:
      statusTestSuccessful = 0
      os.system('echo Status test '+i+' failed.')
      os.system("sync")
      os.system('cat tmp'+NOXoutput+str(i))
   os.system('rm -f '+ 'tmp'+NOXoutput+str(i))

if statusTestSuccessful == 1:
   os.system("echo Status tests were successful!")


if nargin == 5:
   # check to see which diff tests were successful
   diffTestSuccessful = 1
   indx = 1
   for i in diffList:
      if i != 0:
         os.system('echo Diff test '+str(indx)+' failed, params: "'+str(entireParamList[indx-1])+'"')
         diffTestSuccessful = 0
         os.system('echo See '+hostname+':'+os.getcwd()+'/diffOutput_'+hostname+' for diff results.')
      indx += 1
   if diffTestSuccessful == 1:
      os.system("echo Diff tests were successful!")





