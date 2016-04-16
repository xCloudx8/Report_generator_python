#!/usr/bin/python
#Creator: Daniele Castrovilli
#Feb 2015

import urllib2
import subprocess
import datetime
import sys
import re

def getRefPat(pato): #Pato code
	URL = 'http://www.omim.org/search?index=entry&start=1&limit=1&search='+ pato
	
	p = subprocess.Popen(['/usr/bin/lynx', '-dump', URL ], stdout=subprocess.PIPE, stderr=subprocess.PIPE) #Bash from python

	out, err = p.communicate()
	x = out.find('# ') #Truncate text
	if x < 0 : #Check symbols only when looking for phenotype and gene
    		x = out.find('* ') 
	y = out.find('.')  
	out = out[x+2:y] #Overwrite text
	return out

def getDescription(refCode): #Getting desc 
	URL = 'http://www.omim.org/entry/' + refCode
	p = subprocess.Popen(['/usr/bin/lynx', '-dump', URL ], stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
	
	outt, err = p.communicate()
	
	x = outt.find('Description') 
	y = outt.find('REFERENCES')
	outt = outt[x:y] 
	
	outt = re.sub("\[\d+\]","",outt, count=0)
	
	
	return outt
	
