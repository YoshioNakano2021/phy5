# Phy 5, a generator of phylogenetic tree based on genome-wide comparisons of pentamer profiles.

### Phy5

Phy5 is a Shiny web applicaion, a generator of phylogenetic tree based on genome-wide comparisons of pentamer profiles.
It is written in R using the Shiny, Biostrings, ape and pvclust libraries.

Phy5 online version is can be used at http://nakano.no-ip.org:3838/phy5/phy5/ or 
https://phy5.shinyapps.io/Phy5R/.

Upload files of fasta-formated DNA sequences, one file for one strain, and Phy5 
calculates pentamaer profiles.
The pentanucleotide profiles can be downloaded as csv files in raw numbers or proportional numbers,
The resultant phylogenetic tree can be downloaded as a pdf file.

### Phy5cli

Phy5cli is the command line version of Phy5.

Phy5cli is a python script including R function with rpy2 library.  It requires python 3 with the following libraries:
numpy, 
pandas, 
rpy2,
R (version >=4.0) and its packages: importr, 
ocalconverter, 
py2rpy, 
rpy2py, 
pandas2ri, 
numpy2ri, 
Biostrings, 
ape, 
pvclust,
ctc.

Phy5cli is the executable file that works on Linux platform but R and the packages, Biostrings, ape, and pvclust, are required.  The file is saved in the $PATH and Phy5cli executable without needing to specify the full path to the file.

Phy5cli can operate on a single FASTA-formatted file, treating each sequence as a sample.

To run Phy5cli for a single fasta file, run

`$python3 Phy5cli.py -f fasta.fna -s manhattan -g ward`

- `fasta.fna`:Your input single FASTA file.
- `manhattan`:Distance method. It should be chosen one out of 'euclidean','manhattan', 'maximum', 'canberra', and 'binary'.
- `ward`: Agglomeration method. It should be chosen one out of 'ward', 'average', 'single', 'complete', 'mcquitty', 'median', and 'centroid'.

To run the executable version, run

`$Phy5cli -f fasta.fna -s manhattan -g ward`

To run Phy5cli for on a set of FASTA-formatted files, treating each file as a sample, run

`$python3 Phy5cli.py -d fastafolder -s manhattan -g ward`

- `fastafolder`:The directory containing your input fasta files.

Phy5cli produces two csv tabls of frequencies of penta nucleotides in each sequence in sigle fasta file or each file in a directory and proportional values of the frequencies, and the phylogenetic tree whose file name containing distance method and agglomeration method.
