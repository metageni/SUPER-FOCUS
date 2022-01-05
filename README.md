![](logo/superfocus_logo_small.png "Logo")

### SUPER-FOCUS: A tool for agile functional analysis of metagenomic data
* [Installation](#installation)
* [Dependencies](#dependencies)
* [Aligners](#aligners)
* [Download SUPER-FOCUS Database](#database)
* [Running SUPER-FOCUS](#run)
* [General Recomendations](#recomendations)
* [Ouput](#output)
* [Citing](#citing)

## Is SUPER-FOCUS right for you?
This [blog post](https://onestopdataanalysis.com/metagenome-functional-profile/) talks about SUPER_FOCUS. Please read it and make sure the tool is right for you.

## Installation
This will give you command line program:

	pip3 install superfocus

or

	# clone super-focus
	git clone https://github.com/metageni/SUPER-FOCUS.git

	# install super-focus
	cd SUPER-FOCUS && python setup.py install

	# if you do not have super user privileges, you can install it like this
	cd SUPER-FOCUS && python setup.py install --user


## Dependencies
- [Python >= 3.6](http://www.python.org/download)
- [Numpy 1.12.1](https://github.com/numpy/numpy)
- [SciPy 0.19.0](https://github.com/scipy/scipy)  

If you have Python 3.6, you can install both dependencies with:  
`pip3 install -r requirements.txt`

## Aligners
One of the below aligners, which can easily be installed with [`conda`](https://conda.io/docs/):
- [DIAMOND 0.9.14](http://ab.inf.uni-tuebingen.de/software/diamond)
- [RAPSearch2 2.24](http://rapsearch2.sourceforge.net)
- [MMSEQS2](https://github.com/soedinglab/MMseqs2)
- [BLAST 2.6.0](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

To install the aligners, having [`conda`](https://conda.io/docs/installation.html) installed, simply run:  
 `conda install -c bioconda <aligner>`

 Note that they are all available from the [`bioconda`](https://bioconda.github.io/) channel.

## Database

If you have the superfocus databases downloaded already, you can set the `SUPERFOCUS_DB` environment variable to point
to that directory. Alternatively, you can provide the `--alternate_directory` flag to point to that location.

### Installing the databases

Some of the steps below could be automatized. However, many users have had problem with the database formatting, and
it was requested for the initial steps to be manual.

### Downloading prebuilt databases

We have prebuilt several of the databases, so if you have made a `conda` install, choose the right version and 
you should be able to download the databases

***Diamond***

Please check your `diamond` version with `diamond --version` and then read the [diamond documentation](https://github.com/bbuchfink/diamond/wiki/5.-Advanced-topics#database-format-versions) to know which version to download. You can also find out the database version you have installed with `diamond dbinfo`.

Cluster Size | diamond version 1 databases | diamond version 2 databases | diamond version 3 databases
--- | --- | --- | ---
90 |  [90 v1](https://edwards.sdsu.edu/SUPERFOCUS/downloads/conda/diamond_v1/90_clusters.db.dmnd.zip)  | [90 v2](https://edwards.sdsu.edu/SUPERFOCUS/downloads/conda/diamond_v2/90_clusters.db.dmnd.zip)  | [90 v3](https://edwards.sdsu.edu/SUPERFOCUS/downloads/conda/diamond_v3/90_clusters.db.dmnd.zip) 
95 |  [95 v1](https://edwards.sdsu.edu/SUPERFOCUS/downloads/conda/diamond_v1/95_clusters.db.dmnd.zip)  | [95 v2](https://edwards.sdsu.edu/SUPERFOCUS/downloads/conda/diamond_v2/95_clusters.db.dmnd.zip)  | [95 v3](https://edwards.sdsu.edu/SUPERFOCUS/downloads/conda/diamond_v3/95_clusters.db.dmnd.zip) 
98 |  [98 v1](https://edwards.sdsu.edu/SUPERFOCUS/downloads/conda/diamond_v1/98_clusters.db.dmnd.zip)  | [98 v2](https://edwards.sdsu.edu/SUPERFOCUS/downloads/conda/diamond_v2/98_clusters.db.dmnd.zip)  | [98 v3](https://edwards.sdsu.edu/SUPERFOCUS/downloads/conda/diamond_v3/98_clusters.db.dmnd.zip) 
100 |  [100 v1](https://edwards.sdsu.edu/SUPERFOCUS/downloads/conda/diamond_v1/100_clusters.db.dmnd.zip)  | [100 v2](https://edwards.sdsu.edu/SUPERFOCUS/downloads/conda/diamond_v2/100_clusters.db.dmnd.zip)  | [100 v3](https://edwards.sdsu.edu/SUPERFOCUS/downloads/conda/diamond_v3/100_clusters.db.dmnd.zip) 

After downloading, you need to copy these to `lib/python3.8/site-packages/superfocus_app/db/static/diamond` in the same location as superfocus:

e.g. for `90_clusters`:
```bash
mkdir -p  $(which superfocus | sed -e 's#bin/superfocus$#lib/python3.8/site-packages/superfocus_app/db/static/diamond#') &&
unzip -d  $(which superfocus | sed -e 's#bin/superfocus$#lib/python3.8/site-packages/superfocus_app/db/static/diamond#') 90_clusters.db.dmnd.zip
```

***MMSEQS2***

There is only one version of the MMSEQS2 databases and so the installation is easier!

Cluster Size | mseqs2 databases
--- | ---
90 |  [mmseqs_90.zip](https://edwards.sdsu.edu/SUPERFOCUS/downloads/conda/mmseqs2/mmseqs_90.zip)
95 |  [mmseqs_95.zip](https://edwards.sdsu.edu/SUPERFOCUS/downloads/conda/mmseqs2/mmseqs_95.zip)
98 |  [mmseqs_98.zip](https://edwards.sdsu.edu/SUPERFOCUS/downloads/conda/mmseqs2/mmseqs_98.zip)


After downloading, you need to copy these to `lib/python3.8/site-packages/superfocus_app/db/static/diamond` in the same location as superfocus:

e.g. for `90_clusters`:
```bash
mkdir -p  $(which superfocus | sed -e 's#bin/superfocus$#lib/python3.8/site-packages/superfocus_app/db/static/mmseqs2#') &&
unzip -d  $(which superfocus | sed -e 's#bin/superfocus$#lib/python3.8/site-packages/superfocus_app/db/static/mmseqs2#') mmseqs_90.zip
```

### Manual installation

If those databases don't work, you might need to try the manual installation:

#### Download and uncompress
First download the database with the steps below or using your favorite method to download and uncompress files:
```
# download
wget edwards.sdsu.edu/superfocus/downloads/db.zip
# uncompress
unzip db.zip
```
**NOTE**: You can also download the small file named `db_small.zip` and test the instalation before downloading the large file.

#### Format
Now that you downloaded the database, please use the instructions below to format it and move into the database folder.
```
superfocus_downloadDB -i <clusters_folder> -a <aligner> -c <clusters>
```
where
- `<clusters_folder>` is the path to the database you downloaded and uncompressed above (folder `clusters/`)
- `<aligner>` is `rapsearch`, `diamond`, or `blast` (or all of them separated by `,`). You
may choose as many aligners as you want among the three, as long as they are
installed.
- `<clusters>` is the cluster of the database you want to format which are `90`, `95`, `98`, and/or `100`. Default: `90`. If more than one, please separe by comma (e.g. 90,95,98,100).

**NOTE**: RAPSearch2 and DIAMOND won't work properly if you are trying to use a
database formatted with an incorrect version of the aligner. Thus, please
re-run `superfocus_downloadDB` in case any aligner was updated on your
system.


## Run
The main SUPER-FOCUS program is `superfocus`. Here is a list of the
available command line options:

    usage: superfocus    [-h] [-v] -q QUERY -dir OUTPUT_DIRECTORY
                         [-o OUTPUT_PREFIX] [-a ALIGNER] [-mi MINIMUM_IDENTITY]
                         [-ml MINIMUM_ALIGNMENT] [-t THREADS] [-e EVALUE]
                         [-db DATABASE] [-p AMINO_ACID] [-f FAST]
                         [-n NORMALISE_OUTPUT] [-m FOCUS] [-b ALTERNATE_DIRECTORY]
                         [-d] [-l LOG]

    SUPER-FOCUS: A tool for agile functional analysis of shotgun metagenomic data.

    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit
      -q QUERY, --query QUERY
                            Path to FAST(A/Q) file or directory with these files.
      -dir OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                            Path to output files
      -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                            Output prefix (Default: output).
      -a ALIGNER, --aligner ALIGNER
                            aligner choice (rapsearch, diamond, mmseqs2, or blast; default
                            rapsearch).
      -mi MINIMUM_IDENTITY, --minimum_identity MINIMUM_IDENTITY
                            minimum identity (default 60 perc).
      -ml MINIMUM_ALIGNMENT, --minimum_alignment MINIMUM_ALIGNMENT
                            minimum alignment (amino acids) (default: 15).
      -t THREADS, --threads THREADS
                            Number Threads used in the k-mer counting (Default:
                            4).
      -e EVALUE, --evalue EVALUE
                            e-value (default 0.00001).
      -db DATABASE, --database DATABASE
                            database (DB_90, DB_95, DB_98, or DB_100; default
                            DB_90)
      -p AMINO_ACID, --amino_acid AMINO_ACID
                            amino acid input; 0 nucleotides; 1 amino acids
                            (default 0).
      -f FAST, --fast FAST  runs RAPSearch2 or DIAMOND on fast mode - 0 (False) /
                            1 (True) (default: 1).
      -n NORMALISE_OUTPUT, --normalise_output NORMALISE_OUTPUT
                            normalises each query counts based on number of hits;
                            0 doesn't normalize; 1 normalizes (default: 1).
      -m FOCUS, --focus FOCUS
                            runs FOCUS; 1 does run; 0 does not run: default 0.
      -b ALTERNATE_DIRECTORY, --alternate_directory ALTERNATE_DIRECTORY
                            Alternate directory for your databases.
      -d, --delete_alignments
                            Delete alignments
      -l LOG, --log LOG     Path to log file (Default: STDOUT).

    superfocus -q input_folder -dir output_dir

## Query

The query can be one or more fasta or fastq files, or a directory containing those files. We filter for
files that end `.fasta`, `.fastq`, or `.fna`, so please ensure any file that you want processed has one
of those file extensions.

You can provide a mixture of input files or directories, and we will filter the files as appropriate.

For example:

```bash
superfocus -q fastq1.fastq -q fastq2.fastq -q directory/ -dir output
```

will process the two fastq files `fastq1.fastq` and `fastq2.fastq` as well as any `fasta` or `fastq` files in `directory`
and put the output in `output`.

We currently do not handle `gzipped` or otherwise compressed input files.

## Recomendations
- The FOCUS reduction is not necessary if not wanted (it is off by default: set `-focus 1` to run FOCUS reduction)
- Run RAPSearch for short sequences, it is less sensitive for long sequences
- Primarily use DIAMOND for large datasets only. It is slower than blastx for small datasets
- Run mmseqs2 if you are running multiple jobs in parallel (e.g. on a cluster).
- BLAST is known for being really slow

## Output
SUPER-FOCUS output will be add the folder selected by the `-dir` argument.

## Citing
SUPER-FOCUS was written by Genivaldo G. Z. Silva. Feel free to create an [issue or ask questions](https://github.com/metageni/SUPER-FOCUS/issues)

If you use SUPER-FOCUS in your research, please cite:

#### Paper

    Silva, G. G. Z., Green K., B. E. Dutilh, and R. A. Edwards:
    SUPER-FOCUS: A tool for agile functional analysis of shotgun metagenomic data.
	Bioinformatics. 2015 Oct 9. pii: btv584. Website: https://edwards.sdsu.edu/SUPERFOCUS

#### Extended tool manual
    Silva, G. G. Z., F. A. Lopes, and R. A. Edwards
    An Agile Functional Analysis of Metagenomic Data Using SUPER-FOCUS.
	Protein Function Prediction: Methods and Protocols, 2017.
