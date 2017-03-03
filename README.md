# QKdomain
A set of scripts that can be used for domain analysis. The majority of which are developed to process files for input/output to/from InterProScan, MEME Suite, and phylogenetic analysis.

## Scripts
<i>QKdomain_preprocess.py</i> reads InterProScan output and identifies the breadth of domains present within a set of proteins. The user can then specify abbreviations and groupings in subsequent analyses.

<i>QKdomain_process.py</i> reads protein sequence and InterProScan output, then performs one or more of the following analyses:

1. Annotate domain composition in all proteins
2. Extracts protein sequence not associated with a known domain (optional)
3. Extracts protein sequence encompassing specific domains (optional)
4. Extract proteins with a specific domain composition, strict or relaxed (optional)

<i>QKdomain_alignment.py</i> assesses the quality of a multiple sequence alignment in representation and coverage.

## Example
### Identify all thioredoxins in <i>Arabidopsis thaliana</i>
The user is required to identify a set of protein sequences to be evaluated. Next, use InterProScan to evaluate these sequences. The resulting output in tab-separated format (tsv) is used as input for domain analysis. `QKdomain_preprocess.py` is used to assess the domain content present in the protein data set.
```bash
python QKdomain_preprocess.py examples/Athaliana_167_TAIR10.protein.TRX.fa.tsv At_TRX_preprocess.txt
```
The output file `At_TRX_preprocess.txt` is used to manually curate identifiers to generate redundant abbreviations for the next component of the analysis. Below is an example of a set of abbreviations generated from `At_TRX_preprocess.txt`.
```
Coil	CC
PF00085	TRX
PF00462	GLRX
PF00515	TPR
PF00789	UBX
...
```

Next, we identify the domain structure of proteins in the data set.
```bash
python QKdomain_process.py examples/Athaliana_167_TAIR10.protein.TRX.fa examples/Athaliana_167_TAIR10.protein.TRX.fa.tsv examples/At_TRX_abbreviations.txt At_TRX_process.txt
```
The output will be structured with the first column the protein identifier and the second column the protein domain structure.
```
AT1G35620.1	TRX-TRX
AT3G56420.1	TRX
AT5G66410.1	TRX
AT1G14570.1	UBA-TRX-UIM-UBX-UBI
AT3G53220.1	TRX
AT3G15360.1	TRX
AT2G32920.1	TRX-TRX-TRX
AT1G08570.1	CC-TRX
...
```

To export all thioredoxin (TRX) domains in the set of proteins, use the following command:
```bash
python QKdomain_process.py -d TRX examples/Athaliana_167_TAIR10.protein.TRX.fa examples/Athaliana_167_TAIR10.protein.TRX.fa.tsv examples/At_TRX_abbreviations.txt At_TRX_process.txt At.TRX.domain.fa
```
The file `At.TRX.domain.fa` will contain individual non-overlapping domains delineated by InterProScan. Overlapping domains would be integrated into a single domain and would require manual curation in a multiple sequence alignment to correct. The protein identifier is modified to include the start and end positions selected and the domain structure of the original protein.

`QKdomain_process.py` has several additional parameters that allow the user to extend the size of the extracted domain. The commands `--nextend` will extend N-terminal from the start of the domain, whereas `--cextend` will extend C-terminal from the end of the domain. If the value provided by the user is between 0.0 and 1.0, the extension will be weighted based on the size of the domain itself. Otherwise for values greater than 1 will be treated as a fixed number of amino acids.

Example of extracting thioredoxin (TRX) domains with an additional 10 amino acids from both N-terminal and C-terminal ends of the domain.
```bash
python QKdomain_process.py -d TRX --nextend 10 --cextend 10 examples/Athaliana_167_TAIR10.protein.TRX.fa examples/Athaliana_167_TAIR10.protein.TRX.fa.tsv examples/At_TRX_abbreviations.txt At_TRX_process.txt At.TRX.domain.fa
```

InterProScan is limited by the set of predicted motifs that have been identified and annotated. When focusing on a set of proteins with a shared domain, their may be novel motifs in the sequence. Using the command `--undefined` will extract all regions in the protein sequence that are not annotated by a InterProScan domain. A filename is required to export FASTA sequence for all extracted regions. An example is shown below for the <i>A. thaliana</i> thioredoxin data set.
```bash
python QKdomain_process.py -d TRX examples/Athaliana_167_TAIR10.protein.TRX.fa examples/Athaliana_167_TAIR10.protein.TRX.fa.tsv examples/At_TRX_abbreviations.txt At_TRX_process.txt At.TRX.domain.fa --undefined At.TRX.undefined.regions.fa
```
