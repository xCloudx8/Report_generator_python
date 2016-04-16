#!/usr/bin/python
# Uses new software VarSeq
import sys
sys.path.insert(0, '/Private/Development_Branch/Modules')
import gene
import generate
import getPato
import getPathology
import phenogene
import csv
import glob
import os
import os.path
import time
import urllib
from collections import defaultdict

"""Report generator.

   Requirements:
    -Genelist
    -Target Investigaton
"""
class EmptyArgs(Exception): #Class for personal exceptions messages
	def __init__(self,msg):
		"""Creating the error message"""
		self.message = msg

	def __str__(self):
		"""Returning a textual representation of the message"""
		return self.message

def getColumns(inFile, delim="\t", header=True):
    """
    Get columns of data from inFile. The order of the rows is respected

    :param inFile: column file separated by delim
    :param header: if True the first line will be considered a header line
    :returns: a tuple of 2 dicts (cols, indexToName). cols dict has keys that
    are headings in the inFile, and values are a list of all the entries in that
    column. indexToName dict maps column index to names that are used as keys in
    the cols dict. The names are the same as the headings used in inFile. If
    header is False, then column indices (starting from 0) are used for the
    heading names (i.e. the keys in the cols dict)
    """
    cols = {}
    indexToName = {}
    for lineNum, line in enumerate(inFile):
        if lineNum == 0:
            headings = line.split(delim)
            i = 0
            for heading in headings:
                heading = heading.strip()
                if header:
                    cols[heading] = []
                    indexToName[i] = heading
                else:
                    #In this case the heading is actually just a cell
                    cols[i] = [heading]
                    indexToName[i] = i
                i += 1
        else:
            cells = line.split(delim)
            i = 0
            for cell in cells:
                cell = cell.strip()
                cols[indexToName[i]] += [cell]
                i += 1

    return cols, indexToName

#Creating lists
nomeg = []
res = []
i = 0

csvname = raw_input("Insert the .txt filename: ")
naame = csvname #pdf name
data ='/home/daniele/Dropbox/Private/Development_Branch/Data/'

#Main file together with new_report.py
if os.path.exists(csvname+".txt") is True:
	final = file(csvname+".txt", 'r')
	cols, indexToName = getColumns(final)
	final.close()
else:
	raise EmptyArgs("You haven't inserted a valid file: "+csvname+".txt")

#Retrieve from main file
alt = cols['Alternates']
ref = cols['Reference']
Clinvar = cols['Accession']
Phast = cols['PhastCons 46way Primate']
Eff = cols['EFFECT']
hgmd = cols['HGMDID']
rsid = cols['Identifier']
omim = cols['OMIM']
phylop = cols['PhyloP 46way Primate']
polyph = cols['Polyphen2 HDIV Score']
sift = cols['SIFT Score']
genes = cols['Gene Names']
phenotype = cols['PHENOTYPE']
refamino = cols['Ref Amino Acid']
altamino = cols['Alt Amino Acid']
poschr = cols['Chr:Pos']

#Getting columns from every .csv file
if os.path.exists(data+"denovo.txt") is True:
	finalD = file(data+"denovo.txt", 'r')
	colDN, indexToName = getColumns(finalD)
	finalD.close()

	#Denovo
	idDN = colDN['GeneID']
	geneDN = colDN['Gene Names']
	phenotypeDN = colDN['Disorders']

	geneDN = filter(None,geneDN)
	idDN = filter(None,idDN)

else:
	idDN = []
	geneDN = []
	phenotypeDN = []

if os.path.exists(data+"ad.txt") is True:
	finald = file(data+"ad.txt", 'r')
	colAD, indexToName = getColumns(finald)
	finald.close()

	#Autosomical Dominant
	idAD = colAD['GeneID']
	geneAD = colAD['Gene Names']
	phenotypeAD = colAD['Disorders']

	geneAD = filter(None,geneAD)
	idAD = filter(None,idAD)
else:
	idAD = []
	geneAD = []
	phenotypeAD = []

if os.path.exists(data+"ar.txt") is True:
	finalR = file(data+"ar.txt", 'r')
	colAR, indexToName = getColumns(finalR)
	finalR.close()

	#Autosomic Recessive
	idAR = colAR['GeneID']
	geneAR = colAR['Gene Names']
	phenotypeAR = colAR['Disorders']

	geneAR = filter(None,geneAR)
	idAR = filter(None,idAR)
else:
	idAR = []
	geneAR = []
	phenotypeAR = []

if os.path.exists(data+"ch.txt") is True:
	finalC = file(data+"ch.txt", 'r')
	colCH, indexToName = getColumns(finalC)
	finalC.close()

	#Compound heterozygous
	idCH = colCH['GeneID']
	geneCH = colCH['Gene Names']
	phenotypeCH = colCH['Disorders']

	geneCH = filter(None,geneCH)
	idCH = filter(None,idCH)
else:
	idCH = []
	geneCH = []
	phenotypeCH = []

if os.path.exists(data+"xl.txt") is True:
	finalx = file(data+"xl.txt", 'r')
	colXL, indexToName = getColumns(finalx)
	finalx.close()

	#X-Linked
	idXL = colXL['GeneID']
	geneXL = colXL['Gene Names']
	phenotypeXL = colXL['Disorders']

	geneXL = filter(None,geneXL)
	idXL = filter(None,idXL)
else:
	idXL = []
	geneXL = []
	phenotypeXL = []

#Remove any special character that latex does not allow
for i in range(len(phenotype)):
	phenotype[i] = phenotype[i].replace("_"," ")
	phenotype[i] = phenotype[i].replace("/"," ")
	phenotype[i] = phenotype[i].replace("\\x2c"," ")
	phenotype[i] = phenotype[i].replace("open angle"," ")

genes = list(set(genes)) #Prevents duplicates

#Creating variables with the selected columns
otg = cols['All Indiv Freq']
ESP = cols['All MAF']
#phast =cols['ConservationScorePhastCons']

#Public urls but leave the other where they are
omimurl='http://www.omim.org/search?index=entry&start=1&limit=1&search='
clinvarurl='http://www.ncbi.nlm.nih.gov/clinvar/variation/'
dbsnpurl = 'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs='

#Paths to be change in order to be used everywhere
targetpath = '/home/daniele/Dropbox/Private/Development_Branch/Target/'
listpath = '/home/daniele/Dropbox/Private/Development_Branch/GeneList/'

if os.path.exists(targetpath+"target.txt") is True: #Check if the file exist
	with open(targetpath+"target.txt","r") as myFile:
		targetinv = myFile.read().replace('\n', '') #Replace the chr
		print("File target founded!")
else:
	raise ValueError("The file specified does not exist: target.txt")

for i in range (len(genes)):
	nomeg.append(genes[i])
	print genes[i] #printing list of genes analized

	name = nomeg[i]
	nomeg[i] = name.strip().upper() #Upper case

	desc = gene.getRefSeq(nomeg[i]) #Recupero la descrizione
	res.append(desc)

ult = []
for n in range(len(genes)):
	lengen = len(genes[n])+2
	r = "\\begin{huge}\n\color{black}\n"+genes[n]+"\n\end{huge}\n\\\ \\vspace{0.1cm}\\\ \n" +res[n][lengen:]+"\n \\textbf{}\\\ " #Leave the textbf{} blank!
	ult.append(r)

#Change the path to your logo
logopath = '/home/daniele/Dropbox/Private/Development_Branch/Logo'

#Generic Informations
name = raw_input("Insert ordering phisician's name: ")
namepati = raw_input("Insert patient's name: ")
sex = raw_input("Insert patient sex: ")
ancestry = raw_input("Insert patient ancestry: ")
patId = raw_input("Insert patient ID: ")
dob = raw_input("Insert patient date of birth (dd/mm/yyy): ")
infos = raw_input("Insert indications for testing: ")
phenot = raw_input("Insert observed phenotype: ")
inheritance = raw_input("Insert type of inheritance: ")
frequency = raw_input("Insert disease frequency: ")
penetrance = raw_input("Insert disease penetrance: ")
reported = raw_input("Insert date reported (dd/mm/yyy): ")

q = raw_input("Autosomic recessive (y/n) empty for blank: ")
if q == 'y':
	ARding = '\ding{51}'
elif q == 'n':
	ARding = '\ding{56}'
else:
	ARding = '\ding{113}'

q = raw_input("Autosomic dominant (y/n) empty for blank: ")
if q == 'y':
	ADding = '\ding{51}'
elif q == 'n':
	ADding = '\ding{56}'
else:
	ADding = '\ding{113}'

q = raw_input("Compound heterozyous (y/n) empty for blank: ")
if q == 'y':
	CHding = '\ding{51}'
elif q == 'n':
	CHding = '\ding{56}'
else:
	CHding = '\ding{113}'

q = raw_input("DeNovo (y/n) empty for blank: ")
if q == 'y':
	DNding = '\ding{51}'
elif q == 'n':
	DNding = '\ding{56}'
else:
	DNding = '\ding{113}'

q = raw_input("x-Linked (y/n) empty for blank: ")
if q == 'y':
	Xding = '\ding{51}'
elif q == 'n':
	Xding = '\ding{56}'
else:
	Xding = '\ding{113}'

#Writing on file
out_file = open( naame + ".tex","w")
out_file.write(
		"\documentclass[hidelinks, a4paper, 11pt]{article} \n"+
		"\usepackage{helvet} \n"+
		"\usepackage[utf8]{inputenc} \n"+
		"\usepackage{amssymb} \n"+
		"\usepackage{amsfonts} \n"+
		"\usepackage{amsthm} \n"+
		"\usepackage{latexsym} \n"+
		"\usepackage{fancyhdr} \n"+
		"\usepackage{wrapfig,lipsum,booktabs} \n"+
		"\usepackage{mathtools} \n"+
		"\usepackage{python} \n"+
		"\usepackage{hyperref} \n"+
		"\pagestyle{fancy} \n" +
		"\usepackage{sfmath} % sans serif in math \n"+
		"\usepackage{floatflt} % per mettere testo in fianco alle tabelle / figure \n"+
		"\\renewcommand{\headrulewidth}{0pt} \n"+
		"\\renewcommand{\\familydefault}{\sfdefault} \n"+
		"\hyphenpenalty=999999999 \n"+
		"\\tolerance=10000 \n"+
		"\sloppy \n"+
		"\setlength{\parindent}{0pt} \n"+
		"\usepackage{graphicx} \n"+
		"\usepackage[top=2cm, bottom=2cm, left=2cm, right=2cm]{geometry} \n"+
		"\n"+
		"\usepackage{xcolor,colortbl} \n"+
		"\usepackage{pifont} \n"+
		"\definecolor{grigio}{rgb}{0.9,0.9,0.9} \n"+
		"\definecolor{grigiom}{rgb}{0.6,0.6,0.6} \n"+
		"\definecolor{grigiol}{rgb}{0.3,0.3,0.3} \n"+
		"\definecolor{grigioxl}{rgb}{0.2,0.2,0.2} \n"+
		"\n"+
		"\\begin{document} \n"+
		"\let\cleardoublepage\clearpage \n"+
		"\\newgeometry{bottom=7mm} \n"+
		"\\thispagestyle{empty} \n"+
		"\\begin{center} \n"+
		"\graphicspath{{"+logopath+"}} \n"+
		"\includegraphics[width=10cm]{"+logopath+"/pg.pdf} \n"+
		"{\color{grigiom}{--------------------------------------------------------------------------------}\\\}\n"+
		"\huge{Genome Insight-MED report } \n"+
		"\end{center}"+
		"\\normalsize \n"+
		"\\vspace{0.7cm} \n"+
		"\\arrayrulecolor{grigiom} \\begin{table}[!h] \n"+
		"\\begin{tabular}{ | p{5cm} | p{10cm} |}\n"+
		"\hline \n"+
		"\cellcolor{grigio} Ordering phisician &  \cellcolor{grigio}"+name+" \\\ \n"+
		"\hline \n"+
		" Patient name & "+namepati+" \\\ \n"+
		"\hline \n"+
		"\cellcolor{grigio} Patient Sex &  \cellcolor{grigio}"+ sex +" \\\ \n"+
		"\hline \n"+
		" Patient ancestry & "+ancestry+" \\\ \n"+
		"\hline \n"+
		"\cellcolor{grigio} Patient ID &  \cellcolor{grigio}"+  patId+" \\\ \n"+
		"\hline \n"+
		" Patient date of birth & "+ dob +" \\\ \n"+
		"\hline \n"+
		"\cellcolor{grigio} Indication for testing &  \cellcolor{grigio}"+infos+" \\\ \n"+
		"\hline \n"+
		" Observed phenotype &"+ phenot +" \\\ \n"+
		"\hline \n"+
		"\cellcolor{grigio} Mode of inheritance &  \cellcolor{grigio}"+ inheritance +" \\\ \n"+
		"\hline \n"+
		" Disease frequency &"+frequency+" \\\ \n"+
		"\hline \n"+
		"\cellcolor{grigio} Disease penetrance &  \cellcolor{grigio}"+penetrance+" \\\ \n"+
		"\hline \n"
		"Inheritance pattern &"+ ARding+" Autosomal-Recessive "+ ADding+" Autosomal-dominant \\\ "  +
		"&"+CHding+" Compound-heterozygous "+ DNding +" De-Novo "+ Xding+" X-Linked   \\\ " +
		"\hline \n"+
		"\cellcolor{grigio} Date reported &  \cellcolor{grigio}"+ reported +" \\\ \n"+
		"\hline \n"+
		"\end{tabular} \n"+
		"\end{table} \n"+
		"\\section*{Target investigation} \n"+
		targetinv +"\n"+
		"\\vfill{} \n"+
		"\\begin{center} \n"+
		"\\footnotesize{"+"Footnote \\\ \n"+
		"Footnote infos} \n"+
		"\end{center} \n"+
		"\\newpage \n"+
		"\\begin{floatingtable}[l]{ \n"+
		"\\begin{tabular}{r  c} \n"+
		"\cellcolor{grigio}Reference Genome & \cellcolor{grigio}HG19/GRCH37 \\\ \n"+
		"Phastcons & v.46 \\\ \n"+
		"\cellcolor{grigio}Pylop & \cellcolor{grigio}V.46 \\\ \n"+
		"Ensembl/proteins: & v.75 \\\ \n"+
	    "\cellcolor{grigio}Ensembl/genes: & \cellcolor{grigio}V.75 \\\ \n"+
		"dbSNP: & v.141 \\\ \n"+
		"\cellcolor{grigio}HGMD: & \cellcolor{grigio}2014.1 \\\ \n"+
		"HPO: & Build-2014.05 \\\ \n"+
		"\cellcolor{grigio}1000Genomes: & \cellcolor{grigio}v.3-20110521 \\\ \n"+
		"ESP/EVS: & ESP6500S \\\ \n"+
		"\cellcolor{grigio}ClinVar: & \cellcolor{grigio}2014.6 \\\ \n"+
		"COSMIC & v.70 \\\ \n"+
		"\cellcolor{grigio}PGVD & \cellcolor{grigio}V1.0.25 \\\ \n"+
		"\end{tabular}} \n"+
		"\end{floatingtable} \n"+
		"\large{\\textbf{PGVD v.1.0.25}} \n"+
		"\\vspace{1mm} \n"+
		" \n"+
		"\\vspace{3mm} \n"+
		"\large{\\textbf{INCIDENTAL FINDINGS}} \n"+
		"\\vspace{1mm} \n"+
		"56 genes published by the ACMG in 2013 regarding incidental findings in exome and genome \n"+
		"sequencing. The term 'incidental finding' indicate unexpected positive findings as \n"+
		"results of a deliberate search for pathogenic or likely pathogenic alterations in genes \n"+
		"that are not apparently relevant to a diagnostic indication for which the sequencing \n"+
		"test was ordered.\\\ \n"+
		"\large{\\textbf{TARGET GENE LIST PROVIDED BY THE PHYSICIAN}} \\\ \n"+
		"\\vspace{1mm} \n")

#Chek existence of the file
if os.path.exists(listpath+"genelist.txt") is True: #Check if the file exist
	with open(listpath+"genelist.txt","r") as myFile:
		targetinv = myFile.read().replace('\n', '')
		print "Gene list found!"
	out_file.write(targetinv+"\n")
else:
	raise ValueError("The file you have specified does not exist: genelist.txt")

out_file.write(
	"\section*{Clinically significant variants}\n"+
	"Variants already known to be damaging for a particular phenotype.\n"+
	"Selection of variants was based on their presence in clinically related source \n"+
	"like ClinVar and/or HGMD. Only variants with a frequency $<5\%$ in 1000 Genome \n"+
	"project and Exome Sequencing Project were reported.\n"+
	"\\begin{table}[!h] \n"+
	"\\begin{tabular}{|c c c c c c|} \n"+
	"\hline\n"+
	"\cellcolor{grigio}\\textbf{ID} &  \cellcolor{grigio}\\textbf{Gene} & \cellcolor{grigio}\\textbf{Clinvar} & \cellcolor{grigio}\\textbf{HGMD} &\cellcolor{grigio}\\textbf{Phenotype} & \cellcolor{grigio}\\textbf{In provided gene list}\\\ \n")

for i in range(len(rsid)):

	with open(listpath+"genelist.txt") as openfile:
		for line in openfile:
		    for part in line.split():
		        if genes[i] in part:
		            inlist = '\ding{51}'
		        else:
		        	inlist = '\ding{56}'

	if phenotype[i] > 10:
		phenotype2 = phenotype[i][:10]

	out_file.write(
	"\href{"+dbsnpurl+rsid[i][2:]+"}{"+rsid[i]+"}" +" & "+ genes[i] + " & "+ "\href{"+clinvarurl+Clinvar[i]+"}{"+Clinvar[i]+"}"+ " & "+ hgmd[i]+ " & "+ phenotype2 + " & "+ inlist + "\\\ \n")

out_file.write(
	"\\hline \n"+
    "\end{tabular} \n"+
    "\end{table}\n"+
    "\\begin{table}[h]\n"+
	"\\begin{tabular}{p{5cm} | p{10cm}}\n"+
	"\cellcolor{grigio}Pathogenic DM & \cellcolor{grigio} Reported in multiple unrelated cases, with control data. \n"+
	"Functional or expression evidence suggests deleterious effect on gene function. \\\ \n"+
  	"\hline \n"+
 	"\cellcolor{grigio}Likely pathogenic DM? & \cellcolor{grigio} Reported in limited cases, or in single family cohort, with or without control data. \n"+
 	"Limited or no functional evidence available, but overall biological expectations suggestive of deleterious effects.\\\ \n"+
 	"\hline \n"+
 	"\cellcolor{grigio}Risk Factor & \cellcolor{grigio} Something that increases a person's chances of developing a disease.\\\ \n"+
	"\end{tabular} \n"+
	"\end{table}\n"+
	"\section*{Predicted significant variants} \n"+
	"Variants predicted to be damaging.\n"+
	"Evaluation of variants that are not reported in the clinically related sources.\\\ \n"+
	"Variants are filtered according to the following criteria: \n"+
	"\\begin{enumerate} \n"+
	"\item Selection of low frequency variants using the 1000 Genomes project and Exome Sequencing Project (ESP6500). The filtering value was: $ <1\% $ \n"+
	"\item Selection of variants predicted to be damaging to the encoded protein based on evolutionary conservation of the mutated residue and presumed effect on protein structure. \n"+
	"The overall damaging effect was estimated based on three in silico prediction tools: SIFT, Polyphen2 Humavar, PhyloP (SIFT score $< 0.05$ or PolyPhen2 HumVar score $> 0.15$; \n"+
	"PhyloP $> 0.95$). \n"+
	"\item Selection of variants based on Inheritance Pattern: DeNovo, Autosomal recessive, X-linked, Compound heterozygous. Heterozygous variants were excluded if present in the healthy \n"+ 	"parents and if not follow the compound heterozygous mode of inheritance \n"+
	"\end{enumerate} \n"+
	"\\arrayrulecolor{grigiom} \\begin{table}[!h] \n"+
	"\\begin{tabular}{|c c c c|} \n")

out_file.write(
			"\hline \n"+
			"\cellcolor{grigio} \\textbf{ID} & \cellcolor{grigio} \\textbf{GENE} & \cellcolor{grigio}\\textbf{Gene phenotype} & \cellcolor{grigio}\\textbf{In provided gene list}\\\ \n"+
			"\cellcolor{grigio}  & \cellcolor{grigio} \\textbf{Autosomal recessive} & \cellcolor{grigio} & \cellcolor{grigio}\\\ \n")

geneAR = list(set(geneAR))
print "Autosomal recessive gene(s)"
for l in range(len(geneAR)):
	with open(listpath+"genelist.txt") as openfile:
		for line in openfile:
		    for part in line.split():
		        if geneAR[l] in part:
		            inlist = '\ding{51}'
		        else:
		        	inlist = '\ding{56}'

	x = phenotypeAR[l].find('(')
	y = phenotypeAR[l].find(')')
	x = x+1

	w = idAR[l].find(',')

	out_file.write(idAR[l][:w] +" & "+ geneAR[l]+ " & "+phenotypeAR[l][x:y]+ " & "+ inlist + "\\\ " + "\n")

out_file.write("\cellcolor{grigio}  & \cellcolor{grigio} \\textbf{Autosomal Dominant} & \cellcolor{grigio} & \cellcolor{grigio}\\\ \n")
out_file.write("\\hline \n")

print "Autosomal Dominant gene(s)"
geneAD = list(set(geneAD))
for l in range(len(geneAD)):
	with open(listpath+"genelist.txt") as openfile:
		for line in openfile:
		    for part in line.split():
		        if geneAD[l] in part:
		            inlist = '\ding{51}'
		        else:
		        	inlist = '\ding{56}'


	out_file.write(idAD[l]+" & "+ geneAD[l]+ " & "+phenotypeAD[l]+ " & "+ inlist + "\\\ " + "\n")


out_file.write("\cellcolor{grigio}  & \cellcolor{grigio} \\textbf{X-Linked} & \cellcolor{grigio} & \cellcolor{grigio}\\\ \n")
out_file.write("\\hline \n")

print "X-Linked gene(s)"
geneXL = list(set(geneXL))
for l in range(len(geneXL)):
	with open(listpath+"genelist.txt") as openfile:
		for line in openfile:
		    for part in line.split():
		        if geneXL[l] in part:
		            inlist = '\ding{51}'
		        else:
		        	inlist = '\ding{56}'



	out_file.write(idXL[l] +" & "+ geneXL[l]+ " & "+phenotypeXL[l]+ " & "+ inlist + "\\\ " + "\n")

out_file.write("\cellcolor{grigio}  & \cellcolor{grigio} \\textbf{Compound heterozygous} & \cellcolor{grigio} & \cellcolor{grigio}\\\ \n")
out_file.write("\\hline \n")

print "Compound heterozygous gene(s)"
geneCH = list(set(geneCH))
for l in range(len(geneCH)):
	with open(listpath+"genelist.txt") as openfile:
		for line in openfile:
		    for part in line.split():
		        if geneCH[l] in part:
		            inlist = '\ding{51}'
		        else:
		        	inlist = '\ding{56}'



	out_file.write(idCH[l] +" & "+ geneCH[l]+ " & "+phenotypeCH[l]+ " & "+ inlist + "\\\ " + "\n")

out_file.write(
	"\\hline \n"+
	"\end{tabular} \n"+
	"\end{table} \n"+
	"\\newpage \n"+
	"\\vspace{3mm} \n"+
	"\\normalsize \n"+
	"\\newpage \n"+
	"\\begin{huge} \n"+
	"Genetic risk factor \n"+
	"\end{huge} \\\ \n"+
	"{\color{grigiom}{-------------------------------------------------------------}}\\\ \n"+
	"\\vspace{0.1cm} \\\ \n")

for j in range(len(phenotype)):
	nomeg[j]
	out_file.write(ult[j][:500])
	print("\n")
	#moi =
	dbsnp = 'http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs='#Resetting after each iteration

	colorPh = ''
	colorS = ''
	colorP = ''

	effect = phenotype[j]
	ESP1 = ESP[j][:6]
	otg1 = otg[j][:6]

	phastcon = Phast[j]
	Phylop = phylop[j]
	SIFT = sift[j]
	PolyPhen = polyph[j][:7]

	if phylop > 0.95:
		colorPh = '\color{red}'

	if SIFT < 0.05:
		colorS = '\color{red}'

	if PolyPhen >= 0.85:
		colorP = '\color{red}'
	elif PolyPhen > 0.16 and PolyPhen < 0.85:
		colorP = '\color{yellow}'
	else:
		colorP = '\color{green}'

	HGMD = hgmd[j]

	CLINVAR = nomeg[j]
	patoClin = phenotype[j]
	urlclin =  getPato.getClin(CLINVAR, patoClin)
	urlClin = 'http://www.ncbi.nlm.nih.gov/clinvar?term='

	resPat = getPathology.getRefPat(patoClin)
	newsearch = resPat[2:]
	desc = getPathology.getDescription(newsearch)
	x = desc.find('Description') #Truncate string
	y = desc.find('.')
	out = desc[x:y] #Overwrite text
	urlPat = getPathology.getUrl(newsearch)

	Medgene = nomeg[j]
	medpato = patoClin
	urlMed = getPato.getMed(Medgene, medpato)

	Pubmed = 'http://www.ncbi.nlm.nih.gov/pubmed/?term='

	OMIM = getPato.getRefPat(patoClin)
	urlomim = 'http://www.omim.org/search?index=entry&start=1&limit=1&search='+OMIM
	if len(OMIM) < 1:
		OMIM = 'NA'
	phenotype[j].replace(" ","+")

	out_file.write(
		"\\vspace{0.3cm} \\\ \n" +
		"{\color{grigiom}{-------------------------------------------------------------}}\\\ \n"+
		"{\Large{\\textbf{\href{"+dbsnp+rsid[j]+"}{"+rsid[j]+"}}} \hfill \color{red}{"+alt[j] +"} \color{black}{"+ref[j]+"}}\\\ \n"+
		"\\\ \n" +
		"\color{grigiom}{"+ poschr[j] +" \hfill {\color{grigiom}{Reference:}}{\color{black}"+ref[j]+"} \\\ \n"+
		"\\vspace{0.7cm} \\\ \n"+
		#"\\textbf{MOI --}"+ moi+"  \\\ \n"+)

		"\\begin{table}[!h] \n"+
		"{\color{grigiom}{------------------------------------------------------------------}}\\\ \n"+
		"\\raggedright{\color{grigioxl}\Large{Frequency}}\\\ \n"+
		"\\vspace{0.5cm} \n"+

		"\\begin{tabular}{c c} \n"+
		"\\vspace{0.15cm}\n"+
		"\color{grigiol}Effect & " + effect +"\\\ \n"+
		"\\vspace{0.15cm} \n"+
		"\color{grigiol}ESP\% & "+ ESP1 +"\\\ \n"+
		"\\vspace{0.15cm} \n"+
		"\color{grigiol}1kG & "+ otg1 + "\\\ \n"+
	    "\end{tabular}\n \\\ \n"+

	    "\\vspace{0.5cm} \n"+
		"{\color{grigiom}{------------------------------------------------------------------}}\\\ \n"+
		"\\raggedright{\color{grigioxl}\Large{Conservation and damaging}}\\\ \n"+
		"\\vspace{0.5cm} \n"+

		"\\begin{tabular}{c c}\n"+
		"\\vspace{0.15cm}\n"+
		"\color{grigiol}PhastCon & "+ phastcon +"\\\ \n"+
		"\\vspace{0.15cm}\n"+
		"\color{grigiol}Phylop & " +colorPh+ Phylop +"\\\ \n"+
		"\\vspace{0.15cm} \n"+
		"\color{grigiol}SIFT & "+colorS+SIFT+"\\\ \n"+
		"\\vspace{0.15cm} \n"+
		"\color{grigiol}Polyphen & "+colorP+PolyPhen+"\\\ \n"+
		"\end{tabular} \n \\\ \n"+
		"\\vspace{0.5cm} \n"+
		"{\color{grigiom}{------------------------------------------------------------------}}\\\ \n"+
		"\\raggedright{\color{grigioxl}\Large{Resources}}\\\ \n "+
		"\\vspace{0.5cm} \n"+

		"\\begin{tabular}{c c}\n"+
		"\\vspace{0.15cm}\n"+
		"\color{grigiol}HGDM & "+ HGMD+"\\\ \n"+
		"\\vspace{0.15cm}\n"+
		"\color{grigiol}CLINVAR & \href{"+urlClin+genes[j]+"}"+"{"+CLINVAR+"}"+"\\\ \n"+
		"\\vspace{0.15cm} \n"+
		"\color{grigiol}OMIM & \href{"+urlomim+"}{"+ OMIM+"}"+"\\\ \n"+
		"\\vspace{0.15cm}\n"+
		"\color{grigiol}MedGene & \href{"+urlMed+"}"+"{"+Medgene+"}"+"\\\ \n"+
		"\\vspace{0.15cm}\n"+
		"\color{grigiol}PubMed & \href{"+Pubmed+medpato+" "+genes[j]+"}"+"{"+genes[j]+"}"+"\\\ \n"+
		"\end{tabular}\n"+
		"\end{table}\n"+

		"\color{black}{"+
		"\\textbf{"+patoClin+"} }\\\ \n "+
		"\color{grigiol}{"+
		"\\textit{\\textbf{\color{black}{Description}}} \\\ \n"+
		out[12:]+"} \\\ \n" +
		" \\\ \n"+
		"For more information: \href{"+urlomim+"}{\\textbf{"+patoClin+"}} \\\ \n"+
		"\\newpage\n")

	print "Elaborating file please wait..."

out_file.write("\end{document} \n")

print ("Creating PDF...")

out_file.close() #Closing the file
generate.pdf(naame) #Generating pdf file with pdflatex
if os.path.exists(naame +".pdf") is True: #Check if the file is correctly generated
	print(naame+".pdf succesfully created!")
	print("Removing unnecessary files!")
	time.sleep(1)
	#Removing unnecessary files
	os.remove(naame+".log")
	os.remove(naame+".aux")
	os.remove(naame+".out")
	#os.remove(naame+".tex")
	print("\n")
	sys.exit() #Exit the program
else:
	print("File not created check what's wrong") #Check if the file exists
