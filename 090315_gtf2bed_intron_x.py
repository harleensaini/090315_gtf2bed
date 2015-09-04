__author__ = 'proline'

import re
import argparse

parser=argparse.ArgumentParser(description='DEXSeq_intron_annotation GFF to BED for coverage pattern recognition')

parser.add_argument('-i', help='input DEXSeq gff', dest='inputfile', default=None, type=str)
parser.add_argument('-o1', help='output bed filename for gtf2bed_introns', dest='outfile1')
parser.add_argument('-o2', help='output bed filename for intron start, end, center positions', dest='outfile2')
parser.add_argument('-int', help='interval at intron start, end, center position', type=int, dest='interval')

args=parser.parse_args()

inputfile=args.inputfile
outfile1=args.outfile1
outfile2=args.outfile2
interval=args.interval

# Defining a function to make a bed file from gtf file for DEXSeq introns
def gtf2bed_geneintronID():
	intron_bed_output = open(outfile1, 'w')

	with open(inputfile, 'r') as infile:
		for eachline in infile:
			strippedline = eachline.strip('\n')
			eachrowbytab = strippedline.split('\t')
			if eachrowbytab[2] == "intronic_part":
				geneid=re.findall(r'gene_id \"([\w\.]+)\"',eachrowbytab[8])
				intronid=re.findall(r'exonic_part_number \"([\w\.]+)\"',eachrowbytab[8])
				geneintron_id = geneid+intronid

				bed_col1 = eachrowbytab[0]
				bed_col2 = eachrowbytab[3]
				bed_col3 = eachrowbytab[4]
				bed_col4 = "_".join(geneintron_id)
				bed_col5 = '.'
				bed_col6 = eachrowbytab[6]
				bed_file_out = bed_col1, bed_col2, bed_col3, bed_col4, bed_col5, bed_col6
				intron_bed_output.write("\t".join(bed_file_out)+"\n")
		intron_bed_output.close()
		return


gtf2bed_geneintronID()

# Defining a function to make a bed file containing start, center, and end of intron coordinates

#intron_bed_output = open('test_intron.bed', 'r')


def gtf2bed_intron_pattern_ref():
#	interval = int(50)
	intron_startendcenter_output = open(outfile2, 'w')
	with open(outfile1, 'r') as infile4intron:
		for bedline in infile4intron:
			strippedbedline = bedline.strip('\n')
			eachbedlinebytab = strippedbedline.split('\t')

			# add intron size conditional here - chrEnd - chrStart
			if int(eachbedlinebytab[2])-int(eachbedlinebytab[1]) > (interval*3):

				# Conditional NOT NEEDED - if minus strand, then chrStart -50 else +50
				intronchr_col1 = eachbedlinebytab[0]
				chrstart = int(eachbedlinebytab[1])
				chrstartplusint = chrstart+interval
				intronID = eachbedlinebytab[3]
				intron_col5 = eachbedlinebytab[4]
				intronStrand = eachbedlinebytab[5]

				chrend = int(eachbedlinebytab[2])
				chrendminusint = chrend-interval

				chrcenter = int((chrstart+chrend)/2)
				chrcenterplusint = chrcenter+(interval/2)
				chrcenterminusint = chrcenter-(interval/2)

				intronchrstart_A_col2 = str(chrstart)
				intronchrstart_A_col3 = str(chrstartplusint)
				intronID_A_col4 = intronID+"_A"

				intronchrend_B_col2 = str(chrendminusint)
				intronchrend_B_col3 = str(chrend)
				intronID_B_col4 = intronID+"_B"

				intronchrcenter_C_col2 = str(chrcenterminusint)
				intronchrcenter_C_col3 = str(chrcenterplusint)
				intronID_C_col4 = intronID+"_C"

				bedlinePosA = intronchr_col1, intronchrstart_A_col2, intronchrstart_A_col3, intronID_A_col4, intron_col5, intronStrand
				bedlinePosB = intronchr_col1, intronchrend_B_col2, intronchrend_B_col3, intronID_B_col4, intron_col5, intronStrand
				bedlinePosC = intronchr_col1, intronchrcenter_C_col2, intronchrcenter_C_col3, intronID_C_col4, intron_col5, intronStrand

				bedlinePosA_tab = "\t". join(bedlinePosA)
				bedlinePosB_tab = "\t". join(bedlinePosB)
				bedlinePosC_tab = "\t". join(bedlinePosC)

				bedlinesABC = bedlinePosA_tab, bedlinePosC_tab, bedlinePosB_tab

				intron_startendcenter_output.write("\n".join(bedlinesABC)+"\n")
		intron_startendcenter_output.close()
		return

gtf2bed_intron_pattern_ref()