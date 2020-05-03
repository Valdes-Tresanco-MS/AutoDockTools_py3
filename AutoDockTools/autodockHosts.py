# ##################################################################################################
#  Disclaimer                                                                                      #
#  This file is a python3 translation of AutoDockTools (v.1.5.7)                                   #
#  Modifications made by Valdes-Tresanco MS (https://github.com/Valdes-Tresanco-MS)                #
#  Tested by Valdes-Tresanco-MS and Valdes-Tresanco ME                                             #
#  There is no guarantee that it works like the original distribution,                             #
#  but feel free to tell us if you get any difference to correct the code.                         #
#                                                                                                  #
#  Please use this cite the original reference.                                                    #
#  If you think my work helps you, just keep this note intact on your program.                     #
#                                                                                                  #
#  Modification date: 2/5/20 19:51                                                                 #
#                                                                                                  #
# ##################################################################################################

#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/autodockHosts.py,v 1.2 2012/08/07 17:36:20 rhuey Exp $
#
# $Id: autodockHosts.py,v 1.2 2012/08/07 17:36:20 rhuey Exp $
#
#
#

from collections import UserDict

class AutoDockHosts(UserDict):

	def __init__(self, localHostDict):
		UserDict.__init__(self)
		self.update(localHostDict)

	def buildEntry(self, host=None,agPath=None, adPath=None, vinaPath=None, qType='int', userSpecific=0):
		d={}
		d['host']=host
		d['autogrid']=agPath
		d['autodock']=adPath
		d['vina']=vinaPath
		d['queuetype']=qType
		d['userSpecific']=userSpecific
		return d

	def addHost(self,macroName, hostdict=None, **kw):
		host = kw['host']
		agPath = kw['autogrid']
		adPath = kw['autodock']
		vinaPath = kw['vina']
		qType = kw['queuetype']
		userSpecific = kw['userSpecific']
		if not hostdict:
			hostdict=self.buildEntry(host=host, agPath=agPath, adPath=adPath,
				qType=qType,userSpecific=userSpecific)
		self[macroName]=hostdict
		#nb: preexisting macroName entry is overwritten

	def saveHostFile(self, filename, whichOnes='all'):
		#will be in file called .adthosts.py
		#and consist of python code for dictionary called 'newhosts'
		fptr = open(filename, 'w')
		outstr = 'hostMacros={'
		#outstr = 'hosts={'
		fptr.write(outstr)
		#always write a localHost line
		#get the correct macroList here
		if whichOnes=='all':
			macroList=list(self.keys())
		elif whichOnes=='userSpecific':
			#get the one with userSpecific=1, only
			macroList=[]
			for item in list(self.items()):
				if item[1]['userSpecific']:
					macroList.append(item)
		else:
			macroList=[]
			for item in list(self.items()):
				if not item[1]['userSpecific']:
					macroList.append(item)
			#get the other ones...
		for i in range(len(macroList)):
			h=macroList[i][0]
			self.writeEntry(h,fptr)
			if i<len(macroList)-1:
				fptr.write('\t\t},\n')
			else:
				fptr.write('\t\t}\n')
		#close the whole thing
		outstr = '\t}\n'
		fptr.write(outstr)
		fptr.close()

	def writeEntry(self, macroName,fptr):
		outstr='\t\''+macroName+'\': {\n'
		#outstr='\t\''+hostName+'\': {\n'
		fptr.write(outstr)
		d = self[macroName]
		#d = self[hostName]
		klist = list(d.keys())
		for i in range(len(klist)):
			k=klist[i]
			if k=='userSpecific':
				outstr=	'\t\t\''+k + '\': '+ str(d[k])
			else:
				outstr=	'\t\t\''+k + '\': \''+ str(d[k])+'\''
			fptr.write(outstr)
			if i<len(klist)-1:
				outstr= ',\n'
			else:
				outstr= '\n'
			fptr.write(outstr)


	def loadHostFile(self, filename):
		newStuff=__import__(filename)
		self.update(filename.adhosts)
		#self.update(adhosts)


