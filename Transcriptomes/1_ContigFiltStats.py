#!/usr/bin/env python3.5

##__Updated__: 25_08_2017
##__Author__: Xyrus Maurer-Alcala; maurerax@gmail.com
##__Usage__: python XPlateContamID.py --help

from __future__ import print_function

import argparse, os, sys, time
from distutils import spawn
from datetime import datetime, date, time
from argparse import RawTextHelpFormatter,SUPPRESS

from Bio import SeqIO


#----------------------------- Colors For Print Statements ------------------------------#
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
###---------------------------- UPDATE USEARCH PATH BELOW! -----------------------------###
###########################################################################################
	## IF Usearch is IN YOUR PATH then no updating is needed...

def check_usearch_path():

	usearch_path = ''

	if usearch_path == '':
		usearch_path = spawn.find_executable("usearch")
	else:
		pass

	if usearch_path == None:
		print (color.BOLD + '\n\nPlease open this script and check that you have included'\
		+' the PATH to the'+color.BLUE+' "usearch" '+color.END+color.BOLD+'executable.\n\n'+color.END)
		print (color.BOLD+color.BLUE+'LOOK FOR:\n\n'+color.RED\
		+'#------------------------------ UPDATE USEARCH PATH BELOW! -------------------------------#'\
		+color.BLUE+'\n\nThis is somewhere around lines 50 - 80...\n\n'+color.END)

		sys.exit()
	else:
		pass
	
	return usearch_path


###########################################################################################
###------------------------- Checks the Command Line Arguments -------------------------###
###########################################################################################

def check_args():
	
	parser = argparse.ArgumentParser(description=
	color.BOLD + '\n\nThis script is intended to help manage '+color.RED+'Cross-Contamination\n'\
	+color.END+color.BOLD+'among samples sequenced from a '+color.CYAN+'shared HTS library'\
	+' \npreparation/sequencing run!\n'+color.END+color.BOLD+'\nRunning this script should'\
	' occur between steps '+color.RED+'7 '+color.END+color.BOLD+'and '+color.RED+'8 '+color.END+\
	color.BOLD+'of the \nPost-Assembly work-flow (from the KatzLab)\n'+color.END+usage_msg(),\
	usage=SUPPRESS,formatter_class=RawTextHelpFormatter)
	
	required_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Required Options'+color.END)
 
	required_arg_group.add_argument('--id', '-id', action='store', default = '0.98',
	help=color.BOLD+color.GREEN+' Identity threshold to bin contaminant sequences\n '\
	'(default = 0.98)'+color.END)

	optional_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Options'+color.END)

	optional_arg_group.add_argument('-author', action='store_true',
	help=color.BOLD+color.GREEN+' Prints author contact information\n'+color.END)

	if len(sys.argv[1:]) == 0:
		print (parser.description)
		print ('\n')
		sys.exit()

	args = parser.parse_args()

	more_info = return_more_info(args)
	if more_info != None:
		print (parser.description)
		print (more_info)
		sys.exit()
		
	if args.id == None:
		print (color.BOLD+'\n\nMissing necessary identity threshold\n')
		print ('\nFor help and list of options:\n'+color.CYAN+'\n\tpython XPlateContamID.py'\
		' --help'+color.END+'\n')
		sys.exit()

	return args


###########################################################################################
###------------------------------- Script Usage Message --------------------------------###
###########################################################################################

def usage_msg():
	return (color.BOLD+color.RED+'\n\nExample usage:'+color.CYAN+' python'\
	' XPlateContamID.py -id 0.98'+color.END)


##########################################################################################
###-------- Storage for LARGE (Annoying) Print Statements for Flagged Options ---------###
##########################################################################################

def return_more_info(args):

	author = (color.BOLD+color.ORANGE+'\n\n\tQuestions/Comments? Email Xyrus (author) at'\
	' maurerax@gmail.com\n\n'+color.END)

	if args.author == True:
		return author


###########################################################################################
###--------------------------- Does the Inital Folder Prep -----------------------------###
###########################################################################################

def prep_folders():
	
	if not os.path.isdir('Nucleotide') or not os.path.isdir('Protein') or \
	not os.path.isdir('SpreadSheet'):
		if os.path.isdir('Nucleotide') != True:
			os.system('mkdir Nucleotide')			
		if os.path.isdir('Protein') != True:
			os.system('mkdir Protein')			
		if os.path.isdir('SpreadSheet') != True:
			os.system('mkdir SpreadSheet')			
		missing = 'y'
	else:
		missing = 'n'

	if os.path.isdir('Post_XContam_Nucleotide/') != True:
		os.makedirs('Post_XContam_Nucleotide/')
	if os.path.isdir('Post_XContam_Protein/') != True:
		os.makedirs('Post_XContam_Protein/')
	if os.path.isdir('Post_XContam_SpreadSheets/') != True:
		os.makedirs('Post_XContam_SpreadSheets/')
	
	if missing == 'y':	

		print (color.BOLD+'\nNucleotide, Protein and TSV Files MUST be in the '+color.ORANGE+\
		'Nucleotide, Protein, SpreadSheet'+color.END+color.BOLD+' folders!\n\n'+color.END)
		sys.exit()
	
	
###########################################################################################
###----------- Creates OG-Based Fasta Files for ALL Taxa in the "Nuc" Folder -----------###
###########################################################################################

def createOGfiles():

#	OGS_of_interest = ['OG5_126558','OG5_126561','OG5_126568','OG5_126569','OG5_126570','OG5_126572','OG5_126573','OG5_126583','OG5_126588','OG5_126595','OG5_126596','OG5_126599','OG5_126600','OG5_126605','OG5_126607','OG5_126610','OG5_126611','OG5_126613','OG5_126616','OG5_126618','OG5_126623','OG5_126628','OG5_126631','OG5_126641','OG5_126647','OG5_126654','OG5_126655','OG5_126659','OG5_126663','OG5_126665','OG5_126668','OG5_126674','OG5_126676','OG5_126678','OG5_126679','OG5_126681','OG5_126686','OG5_126687','OG5_126688','OG5_126691','OG5_126692','OG5_126694','OG5_126696','OG5_126698','OG5_126706','OG5_126708','OG5_126711','OG5_126714','OG5_126722','OG5_126724','OG5_126731','OG5_126734','OG5_126737','OG5_126741','OG5_126748','OG5_126752','OG5_126757','OG5_126758','OG5_126760','OG5_126762','OG5_126766','OG5_126769','OG5_126773','OG5_126774','OG5_126775','OG5_126783','OG5_126792','OG5_126794','OG5_126795','OG5_126799','OG5_126800','OG5_126803','OG5_126812','OG5_126818','OG5_126820','OG5_126821','OG5_126823','OG5_126826','OG5_126828','OG5_126841','OG5_126845','OG5_126846','OG5_126852','OG5_126853','OG5_126857','OG5_126862','OG5_126867','OG5_126872','OG5_126874','OG5_126877','OG5_126878','OG5_126884','OG5_126886','OG5_126888','OG5_126890','OG5_126895','OG5_126896','OG5_126900','OG5_126906','OG5_126907']

	if not os.path.isdir('OGfiles/'):
		os.makedirs('OGfiles/')

	taxon_file_list = []
	OG_file_list = []	
	LKHdict = {}

	seq_count = 0
	file_count = 0

	print ('\n\n')

	for LKH in os.listdir('Nucleotide'):
		file_count += 1
		taxon_file_list += [LKH]
		infasta = [i for i in SeqIO.parse('Nucleotide/'+LKH,'fasta')]

		for seq in infasta:
			if seq.description.split('OG5_')[-1].split('_')[0] not in LKHdict.keys() and seq.description.split('OG5_')[-1].split('_')[0][0:2] < '16':
				LKHdict.setdefault(seq.description.split('OG5_')[-1].split('_')[0],[]).append(seq)
				seq_count += 1

			elif seq.description.split('OG5_')[-1].split('_')[0][0:2] < '16':
				LKHdict[seq.description.split('OG5_')[-1].split('_')[0]].append(seq)
				seq_count += 1

			print (color.BOLD+'Parsed through '+color.ORANGE+str(seq_count)+color.END+color.BOLD\
			+' Nucleotide Sequences from '+color.CYAN+str(file_count)+color.END+color.BOLD\
			+' Fasta files'+color.END, end='\r')
	print (color.END+'\n\n')					

	for k, v in LKHdict.items():
		num_taxa = list(set([i.description.split('_XX')[0] for i in v]))
		if len(num_taxa) >= 3:
			OG_file_list.append('OG5_'+k+'.fasta')
			temp = list(set([i for i in v]))
			temp.sort(key=lambda x: (-len(x.seq),-int(x.description.split('Cov')[-1].split('_')[0])))			

			with open('OGfiles/OG5_'+k+'.fasta','w+') as w:
				for i in temp:
					w.write('>'+i.description+'\n'+str(i.seq)+'\n')
	
	return taxon_file_list, OG_file_list


###########################################################################################
###-------------- Generates Clusters using a User-defined Minimum Identity -------------###
###########################################################################################

def self_blast(args):
	time.sleep(.1)

	if 'log.txt' not in os.listdir(os.curdir):
		with open('log.txt','w+') as w:
			w.write("OGname\tsequences\tLKH\ttaxa\n")

	if not os.path.exists('UBlast_Output/'):
		os.makedirs('UBlast_Output/')

	for file in os.listdir('OGfiles/'):
		if file.endswith('.fasta'):
			with open('log.txt','a') as w:
				temp_infasta = [i for i in open('OGfiles/'+file).read().split('\n') if i.startswith('>')]
				taxa = list(set([i.split('_XX')[0] for i in temp_infasta]))
				LKH_number = list(set([i.split('_XX_')[-1].split('_')[0] for i in temp_infasta]))
				w.write(str(file)+'\t'+str(len(temp_infasta))+'\t'+str(len(LKH_number))+'\t'+str(len(taxa))+'\n')

			ublast_cmd = usearch_path+' -ublast OGFiles/'+file+' -db OGFiles/'+file+' -self -id '+\
			str(args.id)+' -query_cov 0.7 -evalue 1e-5 -strand plus -blast6out UBlast_Output/'+\
			file.split('.fas')[0]+'.Self.98ID.tsv'
			
			os.system(ublast_cmd)

	time.sleep(.1)


###########################################################################################
###---------- Processes the Binning-Files and Attempts to keep the BEST Data -----------###
###########################################################################################

def Post_SelfBlast_Prep(OG_file_list):

	data_by_taxon = {}

	print (color.BOLD+'\n\nParsing through the OG File List and BLAST outputs\n\n'+color.END)
	
	xcontam_seqs = []		
	OGs_Processed = 0
	
	for file in OG_file_list:
		print (color.BOLD+'Parsed through '+color.ORANGE+str(OGs_Processed)+color.END+color.BOLD\
		+' OG Fasta files for '+color.CYAN+'cross-contaminants'+color.END+color.BOLD+' from a'\
		+color.RED+' SINGLE multiplexed Illumina H(M)iSeq Lane'+color.END, end='\r')

		self_blast_dict = {}	
		seen_seqs_queries = []
		best_seqs = []		
		all_data = [i for i in open('UBlast_Output/'+file.split('.fas')[0]+'.Self.98ID.tsv').read().split('\n') if i != '']

		if len(all_data) < 1:
			continue
		else:
			pass

		target_seq = list(set([i.split('\t')[1] for i in all_data]))
		target_seq.sort(key=lambda x: (-int(x.split('Len')[-1].split('_')[0]), -int(x.split('Cov')[-1].split('_')[0])))
		pass_counter = 1

		for target in target_seq:
			if target not in seen_seqs_queries:
				data_to_parse = list(set([query.split('\t')[0] for query in all_data if query.split('\t')[1] == target]))
				self_blast_dict[pass_counter] = [target]
				self_blast_dict[pass_counter] += data_to_parse
				seen_seqs_queries += data_to_parse
				pass_counter += 1

		all_seqs_seen =  list(set([i for j in self_blast_dict.values() for i in j]))

		for k, v in self_blast_dict.items():
			v.sort(key=lambda x: (-int(x.split('Cov')[-1].split('_')[0]),-int(x.split('Len')[-1].split('_')[0])))		

		stringency_check = {}
		compared = []
		pass2 = []

		for k, v in self_blast_dict.items():
			if int(v[0].split('Cov')[-1].split('_')[0]) >= 10 * int(v[1].split('Cov')[-1].split('_')[0]):
				best_seqs.append(v[0])
				compared.append(k)
			elif v[0].split('_XX')[0] == v[1].split('_XX')[0]:
				best_seqs.append(v[0])
				compared.append(k)
			elif len(v) == 2 and int(v[-1].split('Cov')[-1].split('_')[0]) >= 25:
				best_seqs += [i for i in v]
				compared.append(k)
			else:
				stringency_check[k] = v

		for k, v in stringency_check.items():
			vals = [int(i.split('Cov')[-1].split('_')[0]) for i in v]
			diff = [abs(vals[i+1]-vals[i]) for i in range(len(vals)-1)]	
			avg = sum(diff) / len(diff)
			m = [[vals[0]]]

			for x in vals[1:]:
				if abs(x - m[-1][-1]) < avg:
					m[-1].append(x)
				else:
					m.append([x])

			if len(m) == 1 and m[-1] >= 20:
				best_seqs += [i for i in v]
				pass2.append(k)

			elif m[0][-1]/float(m[1][-1]) >=5:
				best_seqs += [i for i in v[:len(m[0])]]
				pass2.append(k)

		xcontam_seqs += [seq_name for seq_name in all_seqs_seen if seq_name not in best_seqs]
		
		for i in best_seqs:
			data_by_taxon.setdefault(i.split('_XX_')[0],[]).append(i)
			
		OGs_Processed += 1
	
	current_date = datetime.now().strftime('%d_%m_%Y')

	with open('XContam_Seqs_Removed.'+current_date+'.txt','w+') as w:
		for i in xcontam_seqs:
			w.write(i+'\n')	
	print ('\n\n')

	return data_by_taxon, xcontam_seqs


###########################################################################################
###-------------------- Produces the Final Set of "Pruned" Files -----------------------###
###########################################################################################
	
def prep_Fasta_Files(data_by_taxon, xcontam_seqs, taxon_file_list):

	Taxa_written_out = 0

	print (color.BOLD+'Writting out data for for all taxa after removing '+color.CYAN\
	+'cross-contaminant sequences (this can take a while!)\n\n'+color.END)
	
	for file in taxon_file_list:
		if file.split('_XX_')[0] in data_by_taxon.keys():
	
			in_NucFasta = ['>'+i.description+'\n'+str(i.seq)+'\n' for i in SeqIO.parse('Nucleotide/'+file,'fasta') if i.description not in xcontam_seqs]
			in_ProtFasta = ['>'+i.description+'\n'+str(i.seq)+'\n' for i in SeqIO.parse('Protein/'+file.replace('NTD','AA'),'fasta') if i.description  not in xcontam_seqs]
			in_TSV = [i for i in open('SpreadSheets/'+file.split('fas')[0].replace('NTD.ORF','allOGCleanresults')+'tsv').read().split('\n') if i != '' and i.split('\t')[0] not in xcontam_seqs]

			with open('Post_XContam_Nucleotide/'+file.split('.fas')[0]+'.XContClean.fasta','w+') as w:
				for i in in_NucFasta:
					w.write(i)
			with open('Post_XContam_Protein/'+file.split('.fas')[0].replace('NTD','AA')+'.XContClean.fasta','w+') as w:
				for i in in_ProtFasta:
					w.write(i)
			with open('Post_XContam_SpreadSheets/'+file.split('fas')[0].replace('NTD.ORF','allOGCleanresults')+'XContClean.tsv','w+') as w:
				for i in in_TSV:
					w.write(i+'\n')
		else:
			pass
		
		Taxa_written_out += 1


###########################################################################################
###----------- Checks How Well or Poorly a Given Taxon Performed (Survived)-------------###
###########################################################################################

def check_survivorship(taxon_file_list):
	
	current_date = datetime.now().strftime('%d_%m_%Y')
	
	survived = [f.split('_XX')[0] for f in os.listdir('Post_XContam_Nucleotide/')]
	Taxa_died = list(set([t.split('_XX')[0] for t in taxon_file_list if t.split('_XX')[0] not in survived]))
	with open('Dead_Taxa_XContam.'+current_date+'.txt','w+') as w:
		for taxon in Taxa_died:
			w.write(taxon+'\n')	


###########################################################################################
###-------------------------------- Next Script Message --------------------------------###
###########################################################################################

def next_script():

	print (color.BOLD+'\n\nCopy your data for each taxon from the '+color.PURPLE+\
	'"Post_XContam"'+color.END+color.BOLD+' folders into their respective folders (with'\
	' a copy of the SpreadSheet in the '+color.PURPLE+'"TSVtoXML"'+color.END+color.BOLD+\
	' folder!'+color.END)
	
	print(color.BOLD+'\n\nNext Script is: '+color.GREEN+'8_TSVtoXML.py\n\n'+color.END)


##########################################################################################
###--------------- Checks Command Line Arguments and Calls on Functions ---------------###
##########################################################################################
							
def main():

	prep_folders()

	args = check_args()
	
	taxon_file_list, OG_file_list = createOGfiles()
	
	self_blast(min_ID)

	data_by_taxon, xcontam_seqs = Post_SelfBlast_Prep(OG_file_list)

	prep_Fasta_Files(data_by_taxon, xcontam_seqs, taxon_file_list)

	check_survivorship(taxon_file_list)

	next_script()
		
main()	
