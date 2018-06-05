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

Note that SUPER-FOCUS currently runs on Python 2.7 and if you are using a conda
environment based on Python 3+ it might be better to create a separate conda
environment to use with SUPER-FOCUS:

    conda create -n super-focus -c bioconda super-focus 
	conda activate super-focus

This will create a conda environment called `super-focus` (as specified by the
`-n` argument), and install SUPER-FOCUS along with all its dependencies. The second
line activates the newly created `super-focus` conda environment.

## Dependencies
- Python >= 2.6.X < 3.Y: http://www.python.org/download
- Jellyfish: http://www.cbcb.umd.edu/software/jellyfish
- Numpy: http://sourceforge.net/projects/numpy/files/NumPy
- SciPy: http://sourceforge.net/projects/scipy

## Sequence aligners
One of the below aligners:
- RAPSearch2: http://rapsearch2.sourceforge.net
- DIAMOND: http://ab.inf.uni-tuebingen.de/software/diamond
- BLAST: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST


## Download SUPER-FOCUS database
Use `superfocus_downloadDB.py` to download and format the SUPER-FOCUS database
for the available aligners:

```
python superfocus_downloadDB.py <aligner>
```
where `<aligner>` is `rapsearch`, `blast`, or `diamond` (or all of them). You
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
		query file (FASTA or FASTQ format) or folder with multiple FASTA/FASTQ files when -m 1

	-dir string
		output directory

	-o string
		project name (default 'my_project')
	
	-mi float
		minimum identity (default 60 %)

	-ml int
		minimum alignment (amino acids) (default: 15)

	-focus int
		runs FOCUS; 1 does run; 0 does not run: default 0

	-t int
		number of threads (default 8	)

	-e float
		e-value (default 0.00001)

	-db string
		database (DB_90, DB_95, DB_98, or DB_100; default DB_98)

	-p int
		amino acid input; 0 nucleotides; 1 amino acids (default 0)
		
	-k int
		keep original tabular output. 0 delete it / 1 keep it (default 0)

	-a string
		aligner choice (rapsearch, blast (only fasta files), diamond; default rapsearch)

	-fast int
		runs RAPSearch2 or DIAMOND on fast mode - 0 (False) / 1 (True) (default: 1)	
  
	-n int
		normalizes each query counts based on number of hits; 0 doesn't normalize; 1 normalizes (default: 1)

	-r string
		use only the subsystems in the organisms predicted by "-focus"– ncbi / rast annotation  (default: ncbi)
		
## Usage example
```
superfocus.py -q query.fasta -dir myOutputDirectory
```
This will run SUPER-FOCUS using `query.fasta` as a query, and output results into `myOutputDirectory`.

### General recommendations
- The FOCUS reduction is not necessary if not wanted (it is off by default: set `-focus 1` to run FOCUS reduction).
- Run RAPSearch for short sequences, it is less sensitive for long sequences.
- Use BLAST if you want the results to be as sensitive as possible.
- Primarily use DIAMOND for large datasets only. It is slower than blastx for small datasets.
	 
## Output
SUPER-FOCUS output will be add the folder selected by the `-dir` argument.

# Plotting output
Please read [plotting output](https://github.com/metageni/SUPER-FOCUS/tree/master/plotting_output)
for examples on how to plot your output.


# Copyright and License
Copyright (C) 2014-2017  Genivaldo Gueiros Z. Silva

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.
