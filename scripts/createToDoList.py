#! /usr/bin/python

import os
import re
import sys
from math import sqrt

# Global constants

#__________________ showStates __________________
def showStates(constraints):	
	stateFileName = ''
	defaultBaseDir = "dataBase/data/results/pHMC/"
	stateDir = defaultBaseDir+"states"
	states = []
	if(os.path.exists(stateDir)):
		print "List of available state files matching the constraints\n\n"
		stateFiles = os.listdir(stateDir)		
#		StateDescriptorL8x8x8x16Nf1Kap0.13500Lam0.05000Y0.71100Rho1.000R0.500PolDeg64PolAl0.500_level6.dat		
			
		p = re.compile( "StateDescriptorL(?P<l0>[0-9]*)x(?P<l1>[0-9]*)x(?P<l2>[0-9]*)x(?P<l3>[0-9]*)Nf(?P<nf>[0-9]*)Kap(?P<kappa>[0-9.]*)Lam(?P<lambda>[0-9.]*)Y(?P<y>[0-9.]*)Rho(?P<rho>[0-9.]*)R(?P<r>[0-9.]*)PolDeg(?P<polyDegree>[0-9]*)PolAl(?P<polyAlpha>[0-9.]*)_level(?P<level>[0-9]*)" )
		counter = 0
		for entry in stateFiles:							
			matchObj = p.match(entry)			 
			if (matchObj):
				param = matchObj.groupdict()
				param["confCount"] = getConfigCount(defaultBaseDir, entry)
				 
				matchFound = 1
				for cons in constraints:
					if (cons[1] == "="):
						if (float(param[cons[0]]) != cons[2] ): 
							matchFound = 0;
					elif (cons[1] == ">"):
						if (float(param[cons[0]]) < cons[2] ): 
							matchFound = 0;
					elif (cons[1] == "<"):
						if (float(param[cons[0]]) > cons[2] ): 
							matchFound = 0;
				
				if (matchFound):								
					print repr(counter+1) + ". " + entry
					print "Lattice size:\t" + param['l0'] + " x " + param['l1'] + " x " + param['l2'] + " x " + param['l3']
					print "No. fermion doublets:\t" + param['nf']
					print "Kappa:\t" + param['kappa']
					print "Lambda:\t" + param['lambda']
					print "Y:\t" + param['y']
					print "Rho:\t" + param['rho']
					print "R:\t" + param['r']
					print "Polynom degree:\t" + param['polyDegree']
					print "Alpha:\t" + param['polyAlpha']
					print "Level:\t" + param['level']			
					states.append(entry)	
					counter = counter + 1
					print "NO. CONFIGURATIONS CORRESPONDING TO THIS STATE: ", param["confCount"] 
					print "_______________________________________\n"								
	print("End of List")
	return states
#__________________ getConfigCount __________________
def getConfigCount(path, stateDescriptorFileName):
	fileNameExt = re.match( "StateDescriptor([a-zA-Z._0-9]+).dat", stateDescriptorFileName)
	if (fileNameExt):
		fileNameExt = fileNameExt.group(1)
	else:
		fileNameExt = ""
	confCount = 0
	if (len(fileNameExt) > 0):		
		subFolderDir = "dataBase/data/results/pHMC/configurations"
		subFolders = os.listdir(subFolderDir)
		searchString = "subFolder"+fileNameExt;
		for subFolder in subFolders:
			if (subFolder.find(searchString) != -1):
				confDir = subFolderDir+"/"+subFolder
				confs = os.listdir(confDir)
				for confName in confs:
					if confName.find(fileNameExt):
						confCount = confCount + 1;
	return confCount
#__________________ setContraints __________________
def setConstraints(constraintString):
	commands= []
	commands.append( ["kappa", "float"] ) 
	commands.append( ["lambda", "float"] )
	commands.append( ["polyDegree", "int"] )
	commands.append( ["level", "int"] )
	commands.append( ["confCount", "int"] )	
	commands.append( ["l0", "int"] )
	commands.append( ["l1", "int"] )
	commands.append( ["l2", "int"] )
	commands.append( ["l3", "int"] )
		
	string = raw_input("Enter search string (enter 'help' press <RETURN> to see help):"+constraintString+" ")
	while (string == "help" or string == ""):		
		print("Possible commands: ( cmd[<>=]value i.e. confCount>1)")
		helpMsg = ""
		for cmd in commands:
			helpMsg += cmd[0]+" ("+cmd[1]+"), "
		print helpMsg
		string = raw_input("Enter search string:"+constraintString+" ")
		
	string = string.lower()
	string = constraintString + " " + string	
	
	currConstraints = []	
		
	for cmd in commands:		
		if(string.find(cmd[0].lower()) != -1):
			if (cmd[1] == "int"):
				matchObj = re.match(".*"+cmd[0]+"\s?(?P<rel>[<>=])\s?(?P<value>[0-9]*).*", string, re.I)
			elif (cmd[1] == "float"):
				matchObj = re.match(".*"+cmd[0]+"\s?(?P<rel>[<>=])\s?(?P<value>[0-9.]*).*", string, re.I)
						
			if (matchObj):
				if(cmd[1] == "float"):
					currConstraints.append( [cmd[0], matchObj.group("rel"), float(matchObj.group("value")) ])
				elif (cmd[1] == "int"):
					currConstraints.append( [cmd[0], matchObj.group("rel"), int(matchObj.group("value")) ])
			else:
				print "Error in parsing", cmd[0]
	return [currConstraints, string]

#__________________ getConfigCount __________________
def displaySubMenu():
	res = "no";
	while (not res.isdigit()):
		print "\t 1. Specify further constraints"
		print "\t 2. Create ToDo-List"		
		print "\t 0. Quit"		
		res = raw_input("Enter choice: ")		
	return int(res); 
#______________________________ M a i n ______________________________
os.system("clear")
choice = 1
currSearchString = ""
states = [];

[constraints, currSearchString] = setConstraints(currSearchString)
states = showStates(constraints)

while(not choice == 0):
	choice = displaySubMenu();
	if (choice == 1):	
		[constraints, currSearchString] = setConstraints(currSearchString)
		states = showStates(constraints)
	if (choice == 2):
		if(len(states)>0):
			cmds = ""
			cmds = raw_input("Enter command list for module analyzer (i.e. det weight): ")
			analyzerToDoFile = open("AnalyzerToDoList.dat", "a")
			for state in states:
				analyzerToDoFile.write(state+" "+cmds+"\n")
			analyzerToDoFile.close()
			print("End of program. Output written to AnalyzerToDoList.dat")			

