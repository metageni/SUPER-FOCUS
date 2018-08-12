# SUPER-FOCUS
SUPER-FOCUS 0.30: A tool for agile functional analysis of shotgun metagenomic data
(c) Silva, G. G. Z., Green K., B. E. Dutilh, and R. A. Edwards.

If you use SUPER-FOCUS in your research, please cite:

    Silva, G. G. Z., Green K., B. E. Dutilh, and R. A. Edwards: 
    SUPER-FOCUS: A tool for agile functional analysis of shotgun metagenomic data. 
	Bioinformatics. 2015 Oct 9. pii: btv584. Website: https://edwards.sdsu.edu/SUPERFOCUS

# Installation
You can now easily install SUPER-FOCUS using [conda](https://conda.io) via the
[Bioconda](https://bioconda.github.io/) channel. It is as easy as:

    conda install -c bioconda super-focus

Note that SUPER-FOCUS currently runs on Python 3.6 and if you are using a conda
environment based on Python 3+ it might be better to create a separate conda
environment to use with SUPER-FOCUS:

    conda create -n super-focus -c bioconda super-focus
	conda activate super-focus

This will create a conda environment called `super-focus` (as specified by the
`-n` argument), and install SUPER-FOCUS along with all its dependencies. The second
line activates the newly created `super-focus` conda environment.

## Dependencies
- [Python 3.6](http://www.python.org/download)
- [Jellyfish 2.2.6](https://github.com/gmarcais/Jellyfish/releases/tag/v2.2.6). If using macOS, use [bioconda](https://anaconda.org/bioconda/jellyfish)
- [Numpy 1.12.1](https://github.com/numpy/numpy)
- [SciPy 0.19.0](https://github.com/scipy/scipy)

## Sequence aligners
One of the below aligners:
- [DIAMOND 0.9.14](http://ab.inf.uni-tuebingen.de/software/diamond)
- [RAPSearch2 2.24](http://rapsearch2.sourceforge.net). If using macOS, install `brew`, and then `brew install brewsci/science/rapsearch2`



## Download SUPER-FOCUS database
Use `superfocus_downloadDB.py` to download and format the SUPER-FOCUS database
for the available aligners:

```
python superfocus_downloadDB.py <aligner>
```
where `<aligner>` is `rapsearch`, or `diamond` (or all of them). You
may choose as many aligners as you want among the three, as long as they are
installed.

**NOTE**: RAPSearch2 and DIAMOND won't work properly if you are trying to use a
database formatted with an incorrect version of the aligner. Thus, please
re-run `superfocus_downloadDB.py` in case any aligner was updated on your
system.


# Run SUPER-FOCUS
The main SUPER-FOCUS program is `superfocus.py`. Here is a list of the
available command line options:

	-h print help
	
	-q FASTA/FASTQ
		Path to directory with FASTA/FASTQ file(s)

	-dir string
		output directory

	-o string
		output prefix (default 'output_')
	
	-mi float
		minimum identity (default 60 %)

	-ml int
		minimum alignment (amino acids) (default: 15)

	-focus int
		runs FOCUS; 1 does run; 0 does not run: default 0

	-t int
		number of threads (default 4)

	-e float
		e-value (default 0.00001)

	-db string
		database (DB_90, DB_95, DB_98, or DB_100; default DB_98)

	-p int
		amino acid input; 0 nucleotides; 1 amino acids (default 0)

	-a string
		aligner choice (rapsearch (only fasta files) or diamond; default rapsearch)

	-fast int
		runs RAPSearch2 or DIAMOND on fast mode - 0 (False) / 1 (True) (default: 1)	
  
	-n int
		normalizes each query counts based on number of hits; 0 doesn't normalize; 1 normalizes (default: 1)

## Usage example
```
superfocus.py -q query.fasta -dir output_dir
```
This will run SUPER-FOCUS using `query.fasta` as a query, and output results into `myOutputDirectory`.

### General recommendations
- The FOCUS reduction is not necessary if not wanted (it is off by default: set `-focus 1` to run FOCUS reduction).
- Run RAPSearch for short sequences, it is less sensitive for long sequences.
- Primarily use DIAMOND for large datasets only. It is slower than blastx for small datasets.
	 
## Output
SUPER-FOCUS output will be add the folder selected by the `-dir` argument.
