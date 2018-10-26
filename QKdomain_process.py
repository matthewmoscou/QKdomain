#! /usr/bin/python

"""
%prog fasta interproscan_output abbreviations process_summary
Reads InterProScan output and user-defined abbreviations, performs the following analyses:
	1. Defines the non-overlapping domain structure
	2. Exports individual or multiple sequential domains (NB, NB-LRR, CC-NB)
	3. Permits extended domain export (+/- based on length of domain or fixed values)

Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
"""

# modules
import optparse
from optparse import OptionParser
import random
import sets
import string


# import arguments and options
usage = "usage: %prog fasta interproscan_output abbreviations process_summary [selected_domain_output]"
parser = OptionParser(usage=usage)
parser.add_option("-d", "--domain", action="store", type="string", dest="domain", default="", help="Individual domain(s) to export")
parser.add_option("-n", "--nextend", action="store", type="float", dest="nextend", default=-1, help="Extended N-terminal export for selected domain")
parser.add_option("-c", "--cextend", action="store", type="float", dest="cextend", default=-1, help="Extended C-terminal export for selected domain")
parser.add_option("-u", "--undefined", action="store", dest="undefined", default="", help="Export undefined regions (i.e. without annotation)")
parser.add_option("-p", "--plot", action="store_true", dest="plot", default="", help="Plot domain structure of all proteins")
parser.add_option("-i", "--iTOL", action="store", dest="iTOL", default="", help="Generate iTOL domain structure, requires shape and color file")

(options, args) = parser.parse_args()


# random seed
random.seed()


# import protein sequence (FASTA)
fasta_file = open(args[0], 'r')
	
ID_sequence = {}
longest_protein = 0
	
for line in fasta_file.readlines():
	if len(line) > 0:
		if line[0] == '>':
			ID = string.split(line)[0][1:]
			ID_sequence[ID] = ''
		else:
			ID_sequence[ID] += string.split(line)[0]
			
			if len(ID_sequence[ID]) > longest_protein:
				longest_protein = len(ID_sequence[ID])

fasta_file.close()
	

# import domain abbreviations
domain_abbreviation = {}
domain_group_identifiers = {}
color_scheme = {}

abbreviation_file = open(args[2], 'r')

for line in abbreviation_file.readlines():
	line = string.replace(line, '\n', '')
	sline = string.split(line)

	if len(sline) > 1:
		domain_abbreviation[sline[0]] = sline[1]

		if sline[1] not in domain_group_identifiers.keys():
			domain_group_identifiers[sline[1]] = []
	
		domain_group_identifiers[sline[1]].append(sline[0])
	
	if len(sline) > 2:
		color_scheme[sline[1]] = sline[2]

abbreviation_file.close()


# initialize color scheme for total number of domains

hexidecimal = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E']

for index in range(len(domain_group_identifiers.keys())):
	color = '#'

	for cindex in range(6):
		random.shuffle(hexidecimal)
		color += hexidecimal[0]

	if domain_group_identifiers.keys()[index] not in color_scheme.keys():
		color_scheme[domain_group_identifiers.keys()[index]] = color

# initialize color and shape associations, if iTOL command provided
domain_shape_color = {}

if len(options.iTOL) > 0:
	iTOL_associations = open(options.iTOL, 'r')

	for line in iTOL_associations.readlines():
		sline = string.split(line)
		domain_shape_color[sline[0]] = [sline[1], sline[2]]

	iTOL_associations.close()

	iTOL = open(args[3] + '_iTOL.txt', 'w')
	iTOL.write('DATASET_DOMAINS' + '\n')
	iTOL.write('SEPARATOR COMMA' + '\n')
	iTOL.write('DATASET_LABEL,Protein domains' + '\n')
	iTOL.write('COLOR,#ff0000' + '\n')
	iTOL.write('DATASET_SCALE,500,1000,1500' + '\n')
	iTOL.write('DATA' + '\n')


# initialize gene position domain dictionary
gene_position_domain = {}
	
for gene in ID_sequence.keys():
	gene_position_domain[gene] = []

	for index in range(len(ID_sequence[gene])):
		gene_position_domain[gene].append([])


# interproscan analysis
interproscan_file = open(args[1], 'r')
	
interproscan_data = interproscan_file.readlines()

gene_domains = [0]
geneID = string.split(interproscan_data[0])[0]
	
for line in interproscan_data:
	sline = string.split(line, '\t')
	
	for position_index in range(int(sline[6]) - 1, int(sline[7])):
		if sline[4] in domain_abbreviation.keys():
			gene_position_domain[sline[0]][position_index].append(domain_abbreviation[sline[4]])
	
interproscan_file.close()

# domain analysis
process_summary_file = open(args[3], 'w')

if len(args) > 4:
	individual_domain_file = open(args[4], 'w')

if len(options.undefined) > 0:
	undefined_region_file = open(options.undefined, 'w')

if options.plot:
	plot_file = open('gene_structure.R', 'w')
	plot_file.write('plot(c(0,' + str(longest_protein) + '), c(0,-' + str(len(gene_position_domain.keys())) + '), type="n", main="", xlab="", ylab="", axes=F)' + '\n')

gene_index = 1

for gene in gene_position_domain.keys():
	gene_structure = []
	gene_structure_start_stop = []

	positions = range(len(gene_position_domain[gene]))
	local_domains = []
	start = {}

	# if undefined, initialize start position
	if len(options.undefined) > 0:
		if len(gene_position_domain[gene][0]) == 0:
			undefined_start = 0
		else:
			undefined_start = -1

		undefined_domain_index = 1

	for position in positions:
		if len(start.keys()) > 0:
			# if there is an immediate transition from one domain to the next
			if len(sets.Set(local_domains) & sets.Set(gene_position_domain[gene][position])) == 0:
				if len(local_domains) > 0:
					for domain_group in domain_group_identifiers.keys():
						if len(sets.Set(local_domains) & sets.Set([domain_group])) > 0:
							gene_structure.append(domain_group)
							gene_structure_start_stop.append([start[domain_group], position])

				local_domains = []
				start = {}

				for domain in gene_position_domain[gene][position]:
					local_domains.append(domain)
					start[domain] = position

			# else if a new domain is present, add start position
			elif len(sets.Set(gene_position_domain[gene][position]) - sets.Set(local_domains)) > 0:
				for domain_group in domain_group_identifiers.keys():
					if len((sets.Set(gene_position_domain[gene][position]) - sets.Set(local_domains)) & sets.Set([domain_group])) > 0:
						local_domains.append(domain_group)
						start[domain_group] = position

			# else if a domain finished
			elif len(sets.Set(local_domains) - sets.Set(gene_position_domain[gene][position])) > 0:
				for domain_group in domain_group_identifiers.keys():
					if len((sets.Set(local_domains) - sets.Set(gene_position_domain[gene][position])) & sets.Set([domain_group])) > 0:
						gene_structure.append(domain_group)
						gene_structure_start_stop.append([start[domain_group], position])

						local_domains.remove(domain_group)
						del start[domain_group]

		
		# if position contains one or more domains
		elif len(gene_position_domain[gene][position]) > 0:
			# if initializing the start of one or more domains, add to start dictionary
			if len(start.keys()) == 0:
				domain_start = True
			else:
				domain_start = False

			# add domains to local_domains
			for domain in list(sets.Set(gene_position_domain[gene][position])):
				if domain not in local_domains:
					local_domains.append(domain)

					# initialize start of domain if first instance
					if domain_start:
						start[domain] = position
	
			# if undefined was active, export position between domains, reset undefined
			if len(options.undefined) > 0:
				if undefined_start >= 0:
					undefined_region_file.write('>'  + gene + '_' + str(undefined_domain_index) + '_' + str(undefined_start) + '_' + str(position) + '\n')
					undefined_region_file.write(ID_sequence[gene][undefined_start:position] + '\n')

					undefined_start = -1
					undefined_domain_index += 1
		
		# if position does not contain a domain
		elif len(gene_position_domain[gene][position]) == 0:
			# add domain(s) to gene structure that ended, reinitialize
			if len(local_domains) > 0:
				for domain_group in domain_group_identifiers.keys():
					if len(sets.Set(local_domains) & sets.Set([domain_group])) > 0:
						gene_structure.append(domain_group)
						gene_structure_start_stop.append([start[domain_group], position])
				
				start = {}
				local_domains = []

			# if undefined, initialize start position
			if len(options.undefined) > 0:
				if undefined_start < 0:
					undefined_start = position

	# at end of sequence, if domains reach end, add domain(s) to gene structure
	if len(local_domains) > 0:
		for domain_group in domain_group_identifiers.keys():
			if len(sets.Set(local_domains) & sets.Set(domain_group_identifiers[domain_group])) > 0:
				gene_structure.append(domain_group)
				gene_structure_start_stop.append([start[domain_group], position])
	
	if options.plot:
		plot_file.write('lines(c(0,' + str(len(ID_sequence[gene])) + '), c(-' + str(gene_index) + ',-' + str(gene_index) + '), col="black")'+ '\n')

		for domain_index in range(len(gene_structure)):
			plot_file.write('rect(' + str(gene_structure_start_stop[domain_index][0]) + ', -' + str(gene_index + 0.25) + ', ' + str(gene_structure_start_stop[domain_index][1]) + ', -' + str(gene_index - 0.25) + ', col="' + color_scheme[gene_structure[domain_index]] + '")' + '\n')

	if len(options.iTOL) > 0:
		if len(gene_structure) > 0:
			iTOL.write(gene + ',' + str(len(ID_sequence[gene])))

			for domain_index in range(len(gene_structure)):
				if gene_structure[domain_index] in domain_shape_color.keys():
					iTOL.write(',' + domain_shape_color[gene_structure[domain_index]][0] + '|' + str(gene_structure_start_stop[domain_index][0]) + '|' + str(gene_structure_start_stop[domain_index][1]) + '|' + domain_shape_color[gene_structure[domain_index]][1] + '|' + gene_structure[domain_index])
				else:
					iTOL.write(',' + domain_shape_color['default'][0] + '|' + str(gene_structure_start_stop[domain_index][0]) + '|' + str(gene_structure_start_stop[domain_index][1]) + '|' + domain_shape_color['default'][1] + '|' + gene_structure[domain_index])

			iTOL.write('\n')

	# export ordered domain structure
	process_summary_file.write(gene + '\t' + '-'.join(gene_structure) + '\n')

	# if exporting specific domain(s), scan for multiple structures in protein sequence
	if len(options.domain) > 0:
		num_domains = options.domain.count('-') + 1

		for index in range(0, len(gene_structure) - num_domains + 1):
			if options.domain == '-'.join(gene_structure[index:(index+num_domains)]):
				if options.nextend > 0:
					if options.nextend >= 1.0:
						if (gene_structure_start_stop[index][0] - int(options.nextend)) >= 0:
							local_start = gene_structure_start_stop[index][0] - int(options.nextend)
						else:
							local_start = 0
					else:
						if (gene_structure_start_stop[index][0] - (options.nextend * (gene_structure_start_stop[index + num_domains - 1][1] - gene_structure_start_stop[index][0]))) >= 0:
							local_start = gene_structure_start_stop[index][0] - int(options.nextend * (gene_structure_start_stop[index + num_domains - 1][1] - gene_structure_start_stop[index][0]))
						else:
							local_start = 0

				else:
					local_start = gene_structure_start_stop[index][0]

				if options.cextend > 0:
					if options.nextend >= 1.0:
						if (gene_structure_start_stop[index + num_domains - 1][1] + int(options.cextend)) < len(ID_sequence[gene]):
							local_stop = gene_structure_start_stop[index + num_domains - 1][1] + int(options.cextend)
						else:
							local_stop = len(ID_sequence[gene])
					else:
						if (gene_structure_start_stop[index + num_domains - 1][1] + (options.cextend * (gene_structure_start_stop[index + num_domains - 1][1] - gene_structure_start_stop[index][0]))) < len(ID_sequence[gene]):
							local_stop = gene_structure_start_stop[index + num_domains - 1][1] + int(options.cextend * (gene_structure_start_stop[index + num_domains - 1][1] - gene_structure_start_stop[index][0]))
						else:
							local_stop = len(ID_sequence[gene])
				else:
					local_stop = gene_structure_start_stop[index + num_domains - 1][1]

				if len(args) > 4:
					individual_domain_file.write('>' + gene + '_' + str(local_start) + '_' + str(local_stop) + ' ' + '-'.join(gene_structure) + '\n')
					individual_domain_file.write(ID_sequence[gene][local_start:local_stop] + '\n')

	gene_index += 1

process_summary_file.close()

if len(args) > 4:
	individual_domain_file.close()

if len(options.undefined) > 0:
	undefined_region_file.close()

if len(options.iTOL) > 0:
	iTOL.close()
