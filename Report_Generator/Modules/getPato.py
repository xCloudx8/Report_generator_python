#!/usr/bin/python

#Feb 2015
import subprocess
import datetime
import re

def getRefPat(pato): #Getting RefCode of pathology
	URL = 'http://www.omim.org/search?index=entry&start=1&limit=3&search='+ pato
	
	p = subprocess.Popen(['/usr/bin/lynx', '-dump', URL ], stdout=subprocess.PIPE, stderr=subprocess.PIPE) #Bash from python

	out, err = p.communicate()
	x = out.find('# ') #Truncate text
	y = out.find('-')  
	out = out[x:y] #Overwrite text
	x = out.find(' ')
	out = out[x:]
	return out

def getDescription(refCode): #Getting description
	URL = 'http://www.omim.org/entry/' + refCode
	p = subprocess.Popen(['/usr/bin/lynx', '-dump', URL ], stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
	
	outt, err = p.communicate()
	
	x = outt.find('Description')
	if x < 0 : #Check symbols only when looking for phenotype and gene
    		x = out.find('Text') 
	y = outt.find('.')
	outt = outt[x:y] 
	outt = re.sub("\[\d+\]","",outt, count=0)
	return outt

def getMed(gene, pato):
	URL = 'http://www.ncbi.nlm.nih.gov/medgen?term='+gene+'%20'+ pato
	p = subprocess.Popen(['/usr/bin/lynx', '-dump', URL ], stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
	outt, err = p.communicate()
	
	return URL

def getClin(gene,pato):
	URL = 'http://www.ncbi.nlm.nih.gov/clinvar?term='+gene+'%20'+ pato
	p = subprocess.Popen(['/usr/bin/lynx', '-dump', URL ], stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
	
	outt, err = p.communicate()
	
	x = outt.find('m') 
	y = outt.find('[')
	outt = outt[x:y] 
	outt = re.sub("\[\d+\]","",outt, count=0)
	return outt

def getClinUrl(url):
	URL = 'http://www.ncbi.nlm.nih.gov/clinvar/variation/' + url
	return URL

def getUrl(pato): #Getting RefCode of pathology
	URL = 'http://www.omim.org/search?index=entry&start=1&limit=1&search='+ pato
	
	return URL

def getDbSnp(gene):
	URL = 'http://www.ncbi.nlm.nih.gov/snp/?term=' + gene
	return URL

def getOTG(gene):
	URL = 'http://browser.1000genomes.org/Homo_sapiens/Search/Results?site=ensembl&q=' +gene
	return URL
