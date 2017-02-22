# QKdomain
A set of scripts that can be used for domain analysis. The majority of which are developed to process files for input/output to/from InterProScan, MEME Suite, and phylogenetic analysis.

## Scripts
<i>QKdomain_preprocess.py</i> reads InterProScan output and identifies the breadth of domains present within a set of proteins. The user can then specify abbreviations and groupings in subsequent analyses.

<i>QKdomain_extract.py</i> reads protein sequence and InterProScan output, then performs one or more of the following analyses:

1. Annotate domain composition in all proteins
2. Extracts protein sequence not associated with a known domain
3. Extracts protein sequence encompassing specific domains (optional)
4. Extract proteins with a specific domain composition, strict or relaxed (optional)

<i>QKdomain_alignment.py</i> assesses the quality of a multiple sequence alignment in representation and coverage.
