#!/usr/bin/env python3.5

##__Updated__: 18_08_2017
##__Author__: Xyrus Maurer-Alcala; maurerax@gmail.com
##__Usage__: python 2a_remove_rDNA.py --help

##########################################################################################
## This script is intended to identify and isolate SSU/LSU sequences 					##
## Prior to running this script, ensure the following:									##
##																						##
## 1. You have assembled your transcriptome and COPIED the 'assembly' file 				##
##    (contigs.fasta, or scaffolds.fasta) to the PostAssembly Folder					##
## 2. Removed small sequences (usually sequences < 300bp) with ContigFilterPlusStats.py	##
## 3. Have the Databases set up correctly (e.g. with BLAST or USearch) and in their 	##
##	  respective folders! See the manual if you need help								##
##																						##
## 								COMMAND Example Below									##
##																						##
## 				E-mail Xyrus (author) for help if needed: maurerax@gmail.com			##
##																						##
##							Next Script(s) to Run: 										##
##								2b_removeBact.py										##
##																						##
##########################################################################################


import argparse, os, sys
from argparse import RawTextHelpFormatter,SUPPRESS
from Bio import SeqIO


#------------------------------ Colors For Print Statements ------------------------------#

class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   ORANGE = '\033[38;5;214m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'


#------------------------------- Main Functions of Script --------------------------------#

###########################################################################################
###--------------------- Parses and Checks Command-Line Arguments ----------------------###
###########################################################################################
	
def check_args():

	parser = argparse.ArgumentParser(description=
	color.BOLD+'\nThis script will remove '+color.RED+'rDNA contigs (both SSU and LSU)'+color.END\
	+color.BOLD+'\nfrom your Assembly using a set of '+color.RED+'SSU/LSU rDNAs '+color.END\
	+color.BOLD+'from diverse\n'+color.ORANGE+'Eukaryotes, Bacteria and Archaea'+color.END\
	+color.BOLD+'.'+color.END+usage_msg(), usage=SUPPRESS,formatter_class=RawTextHelpFormatter)
	
	required_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Required Options'+color.END)
 
	required_arg_group.add_argument('--input_file','-in', action='store', required=True,
	help=color.BOLD+color.GREEN+"Fasta file of Nucleotide sequences"+color.END)

	optional_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Options'+color.END)
	optional_arg_group.add_argument('--threads','-t', default='2',
	help=color.BOLD+color.GREEN+' Number of threads to use for BLAST\n (default = 2)\n'+color.END)
	optional_arg_group.add_argument('-author', action='store_true',
	help=color.BOLD+color.GREEN+' Print author contact information\n'+color.END)

	if len(sys.argv[1:]) == 0:
		print (parser.description)
		print ('\n')
		sys.exit()

	elif '-author' in sys.argv[1:]:
		print (parser.description)
		print (color.BOLD+color.ORANGE+'\n\n\tQuestions/Comments? Email Xyrus (author) at maurerax@gmail.com\n\n'+color.END)
		sys.exit()

	elif 0<len(sys.argv[1:])<2:
		if '-h' or '--help' in sys.argv[1:]:
			pass
		else:
			print (color.RED+color.BOLD+'\n\nMissing Required Command Line Arguments!'+color.END)
			print (color.BOLD+'\n\nFor usage and all options: '+color.DARKCYAN+'python '\
			+'python 2a_remove_rDNA.py -h\n\n'+color.END)
			sys.exit()

	args = parser.parse_args()
	
	if args.input_file.endswith('bp.fasta') != True:
		print (color.BOLD + '\n\nCheck that you are giving an appropriately Named/Processed'\
		'Fasta file(s) to this script\n\nNOTE that this script CURRENTLY expects your'\
		' Fasta files to contain '+color.RED+ '"rna"'+color.END+color.BOLD+' in \nthe Fasta File'\
		' Name and must end with ' + color.RED + '"bp.fasta"\n\n' + color.END)
		sys.exit()
		
	return args


###########################################################################################
###------------------------------- Script Usage Message --------------------------------###
###########################################################################################

def usage_msg():
	return color.BOLD+color.RED+'\n\nExample usage:'+color.CYAN+' python 2a_remove_rRNA.py --input_file ../Op_me_Xxma_rna.200bp.fasta'+color.END


###########################################################################################
###--------------------------- Does the Inital Folder Prep -----------------------------###
###########################################################################################

def prep_folders(args):
	rRNA_folder = '../'+args.input_file.split('/')[-1].split('_rna')[0]+'/rRNA_Removal/'
	
	if os.path.isdir(rRNA_folder) != True:
		os.system('mkdir '+rRNA_folder)


###########################################################################################
###---------------------- Uses BLAST to identify SSU/LSU Sequences ---------------------###
###########################################################################################

def remove_rDNA(args): 	

	rRNA_folder = '../'+args.input_file.split('/')[-1].split('_rna')[0]+'/rRNA_Removal/'

	blast_output = rRNA_folder+args.input_file.split('/')[-1].split('_rna')[0]+'_allSSULSUresults.tsv'
	
	BLASTN_cmd = 'blastn -query '+args.input_file+' -evalue 1e-10 -max_target_seqs 1 -outfmt'\
	' 6 -db ../../db_BvsE/SSULSUdb -num_threads '+args.threads+' -out '+blast_output
		
	print (color.BOLD+'\n\nBLASTing '+color.DARKCYAN+args.input_file.split('/')[-1]+color.END\
		+color.BOLD+ ' against the rDNA database\n\n' + color.END)

	os.system(BLASTN_cmd)		

	rDNA_Hits = list(set([i.split('\t')[0] for i in open(blast_output).readlines()]))

	print (color.BOLD+'Binning Sequences from '+color.DARKCYAN+args.input_file.split('/')[-1]\
		+color.END+color.BOLD+' as rDNA OR Potentially Protein-Coding\n\n'+color.END)

	no_SSULSU = 0
	with_SSULSU = 0

	inFasta = [seq_rec for seq_rec in SeqIO.parse(args.input_file,'fasta')]

	with open(rRNA_folder+args.input_file.split('/')[-1].split('_rna')[0]+'_rRNAseqs.fasta','w+') as HasSSU:
		for seq_rec in inFasta:
			if seq_rec.description in rDNA_Hits:
				HasSSU.write('>'+seq_rec.description+'\n'+str(seq_rec.seq)+'\n')
				with_SSULSU += 1

	with open(rRNA_folder+args.input_file.split('/')[-1].split('_rna')[0] + '_NorRNAseqs.fasta','w+') as NoSSU:
		for seq_rec in inFasta:
			if seq_rec.description not in rDNA_Hits:
				NoSSU.write('>'+seq_rec.description+'\n'+str(seq_rec.seq)+'\n')
				no_SSULSU += 1

	return str(with_SSULSU), str(no_SSULSU)
	
	
###########################################################################################
###--------------------------- Updates Log of SSU/LSU Removal --------------------------###
###########################################################################################
	
def update_log(args, with_SSU, no_SSU):

	if os.path.isdir('../PostAssembly_Logs/') != True:
		os.system('mkdir ../PostAssembly_Logs/')
	
	print (color.BOLD+'There are '+color.RED+with_SSU+' rRNA contigs'+color.END+color.BOLD\
	+' and '+color.PURPLE+no_SSU+' Putative Protein-coding contigs'+color.END+color.BOLD\
	+' in '+color.DARKCYAN+args.input_file.split('/')[1]+'\n' + color.END)
	
	with open('../PostAssembly_Logs/'+args.input_file.split('/')[1].split('.fas')[0]+'.Log.txt','a') as LogFile:
		LogFile.write('rDNA Contigs\t'+with_SSU+'\tn/a\tn/a\n')
		LogFile.write('Non-rDNA Contigs\t'+no_SSU+'\tn/a\tn/a\n')


###########################################################################################
###-------------------------------- Next Script Message --------------------------------###
###########################################################################################

def next_script(args):

	print (color.BOLD+'\nLook for '+color.DARKCYAN+args.input_file.split('/')[1].split('_rna')[0]\
		+ '_NorRNAseqs.fasta'+color.END+color.BOLD+' in the '+args.input_file.split('/')[1].split('_rna')[0]\
		+' Folder\n\n' + color.END)
	print (color.BOLD + 'Next Script is: ' + color.GREEN + '2b_remove_Bact.py\n\n'+ color.END)


###########################################################################################
###-------------------------- Cleans Up the PostAssembly Folder ------------------------###
###########################################################################################

def clean_up(args):
	
	home_folder = '../'+args.input_file.split('/')[-1].split('_rna')[0]+'/'
	
	os.system('mv '+args.input_file+' '+home_folder+'SizeFiltered/')

	os.system('cp '+home_folder+'rRNA_Removal/*NorRNA*.fasta '+home_folder)

##########################################################################################
###--------------- Checks Command Line Arguments and Calls on Functions ---------------###
##########################################################################################

def main():

	args = check_args()

	prep_folders(args)

	with_SSULSU, no_SSULSU = remove_rDNA(args)
 	
	update_log(args, with_SSULSU, no_SSULSU)
 	
	clean_up(args)
 
	next_script(args)

main()
