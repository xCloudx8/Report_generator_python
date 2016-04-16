#! /usr/bin/python

import subprocess

def getRefSeq( gene ):
	"""Module RefSeq.

	Getting gene description thanks to lynx (must be installed)
	"""
	URL = "http://www.ncbi.nlm.nih.gov/IEB/Research/Acembly/av.cgi?db=human&c=Gene&a=fiche&l=" + gene #Url usage

	p = subprocess.Popen(['/usr/bin/lynx', '-dump', URL ], stdout=subprocess.PIPE, stderr=subprocess.PIPE) #bash command in python

	out, err = p.communicate()

	x = out.find('['+gene+']') #Truncate text
	y = out.find('[provided')  
	out = out[x:y] #Overwrite

	return out #output



