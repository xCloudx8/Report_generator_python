#!/usr/bin/python
#Feb 2015

import subprocess
import datetime
import re

def getRefPat(pato): #Getting RefCode of pathology
	URL = 'http://www.omim.org/search?index=entry&start=1&limit=3&search='+ pato
	
	p = subprocess.Popen(['/usr/bin/lynx', '-dump', URL ], stdout=subprocess.PIPE, stderr=subprocess.PIPE) #Bash from python
	
	out, err = p.communicate()
	x = out.find('# ') #Truncate string
	if x < 0:
		x = out.find('* ')
	elif x < 0:
		x = out.find('% ')
	y = out.find('.')
	x = x+2  
	out = out[x:y] #Overwrite text
	
	return out


def getDescription(refCode): #Getting description
	URL = 'http://www.omim.org/entry/'+refCode+'?search='+refCode
	p = subprocess.Popen(['/usr/bin/lynx', '-dump', URL ], stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
	outt, err = p.communicate()
	
	#outt = re.sub("\[\d+\]","",outt, count=0)
	return outt

def getUrl(refCode): #Getting description
	URL = 'http://www.omim.org/entry/' + refCode
	
	return URL

