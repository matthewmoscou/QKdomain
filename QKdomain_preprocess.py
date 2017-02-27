#! /usr/bin/python

"""
%prog interproscan_output preprocess_summary
Reads InterProScan output and identifies the breadth of domains present within a set of proteins.

Author: Matthew Moscou <matthew.moscou@tsl.ac.uk> and Helen Brabham <helen.brabham@tsl.ac.uk>
This performs the following:
   1. Reads InterProScan output
   2. Identifies the breadth of domains present within a set of proteins
      This is used for subsequent analyses, allowing Pfam or other domains to be merged

Improvements from existing script set (22 February 2017):

Future improvements to include:
"""

# modules
import optparse
from optparse import OptionParser

import sets

import string


# function
def average(number_list):
	return float(sum(number_list)) / len(number_list)


# import arguments and options
usage = "usage: %prog interproscan_output preprocess_summary [abbreviations]"
parser = OptionParser(usage=usage)
(options, args) = parser.parse_args()

domain_abbreviation = {}

# import domain abbreviations
if len(args) > 2:
	abbreviation_file = open(args[2], 'r')
	
	for line in abbreviation_file.readlines():
		line = string.replace(line, '\n', '')
		sline = string.split(line, '\t')
	
		domain_abbreviation[sline[0]] = sline[1]
	
	abbreviation_file.close()

# import interproscan output
interproscan_file = open(args[0], 'r')

# initialize gene position domain dictionary
# domainID -> [software, short name, long name, [geneID.1, geneID.2, ...]]
domain_software_annotation_frequency_length = {}

for line in interproscan_file.readlines():
	line = string.replace(line, '\n', '')
	sline = string.split(line, '\t')
	
	if sline[4] not in domain_software_annotation_frequency_length.keys():
		if len(sline) < 13:
			domain_software_annotation_frequency_length[sline[4]] = [sline[3], sline[5], '', [], []]
		else:
			domain_software_annotation_frequency_length[sline[4]] = [sline[3], sline[5], sline[12], [], []]
	
	domain_software_annotation_frequency_length[sline[4]][3].append(sline[0])
	domain_software_annotation_frequency_length[sline[4]][4].append(float(int(sline[7]) - int(sline[6]) + 1))

interproscan_file.close()

preprocess_summary = open(args[1], 'w')

preprocess_summary.write('domainID' + '\t' + 'software' + '\t' + 'shortname' + '\t' + 'longname' + '\t' + 'genes_with_domain' + '\t' + 'total_observed_domains' + '\t' + 'average_domain_length')

if len(args) > 2:
	preprocess_summary.write('\t' + 'abbreviations')

preprocess_summary.write('\n')

for domain in domain_software_annotation_frequency_length.keys():
	preprocess_summary.write(domain + '\t' + '\t'.join(domain_software_annotation_frequency_length[domain][:3]) + '\t' + str(len(sets.Set(domain_software_annotation_frequency_length[domain][3]))) + '\t' + str(len(domain_software_annotation_frequency_length[domain][3])) + '\t' + str(average(domain_software_annotation_frequency_length[domain][4])))
	
	if len(args) > 2:
		if domain in domain_abbreviation.keys():
			preprocess_summary.write('\t' + domain_abbreviation[domain])
		else:
			preprocess_summary.write('\t')
	
	preprocess_summary.write('\n')

preprocess_summary.close()
