# This is a module with the intention of having it called up from outside
# This module is copyF.py
def copyFiles(inFileName, outFileName):
	   inFile = open(inFileName, 'r')
	   outFile = open(outFileName, 'w')
	   # Read & Copy each line   
	   for inLine in inFile.readlines():
	       outLine = inLine 
	       #inLinemod = modifyLine( inLine )    
	       outFile.writelines( outLine )
	   inFile.close()
	   outFile.close()
#######	  