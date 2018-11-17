#! /usr/bin/python
"""
%prog interproscan_results abbreviations domain output
Reads InterProScan output and user-defined abbreviations, performs the following action:
    1. Extracts all InterProScan output that contains a specific domain of interest

Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
"""

# modules
import optparse
from optparse import OptionParser
import string


# import arguments and options
usage = "usage: %prog interproscan_results abbreviations domain output"
parser = OptionParser(usage=usage)
(options, args) = parser.parse_args()

# import domain abbreviations
domain_abbreviation = {}
domain_group_identifiers = {}
color_scheme = {}

abbreviation_file = open(args[1], 'r')

for line in abbreviation_file.readlines():
    line = string.replace(line, '\n', '')
    sline = string.split(line)

    if len(sline) > 1:
        domain_abbreviation[sline[0]] = sline[1]

        if sline[1] not in domain_group_identifiers.keys():
            domain_group_identifiers[sline[1]] = []

        domain_group_identifiers[sline[1]].append(sline[0])

abbreviation_file.close()

# interproscan analysis
interproscan_file = open(args[0], 'r')

interproscan_data = interproscan_file.readlines()

gene_interproscan = {}

for line in interproscan_data:
    sline = string.split(line, '\t')

    if sline[0] not in gene_interproscan.keys():
        gene_interproscan[sline[0]] = [[], []]

    if sline[4] in domain_abbreviation.keys():
        gene_interproscan[sline[0]][0].append(domain_abbreviation[sline[4]])

    gene_interproscan[sline[0]][1].append(line)

interproscan_file.close()

# export genes with domain
interproscan_output = open(args[3], 'w')

for gene in gene_interproscan.keys():
    if args[2] in gene_interproscan[gene][0]:
        for line in gene_interproscan[gene][1]:
            interproscan_output.write(line)

interproscan_output.close()

