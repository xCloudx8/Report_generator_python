#!/usr/bin/python
import subprocess

def pdf(filename):
        """Generates a .pdf file.
		
		Use the pdf latex instance (must be installed)
        """
        p = subprocess.Popen(['/usr/bin/pdflatex', filename +".tex" ], stdout=subprocess.PIPE, stderr=subprocess.PIPE) #bash command in python
        out, err = p.communicate()  
        return out #output
        
