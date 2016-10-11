#SUPER-FOCUS
A tool for agile functional analysis of shotgun metagenomic data || version 0.26
---------------------------------------------------------------------------------------------------------------------------------------
(c)   Silva, G. G. Z., Green K., B. E. Dutilh, and R. A. Edwards: 
		SUPER-FOCUS: A tool for agile functional analysis of shotgun metagenomic data. Bioinformatics. 2015 Oct 9. pii: btv584.
website: 	https://edwards.sdsu.edu/SUPERFOCUS
---------------------------------------------------------------------------------------------------------------------------------------

#############################################################################################################
Program: superfocus__downloadDB.py: Downloads and formats the SUPER_FOCUS database for the available aligners
#############################################################################################################
(1) USAGE
python superfocus__downloadDB.py aligner
Example: python superfocus__downloadDB.py rapsearch blast diamond

You may choose as many aligners as you want among the three, as long as they are installed.

(2) RECOMMENDATIONS
	- RAPSearch2 and DIAMOND don't work properly using a already formatted database with a newer version of the 
	  aligner. Thus, please re-run 'superfocus__downloadDB.py' in the case of any aligner was updated in the 
	  system.

#############################################################################################################
Program: superfocus.py: SUPER-FOCUS main program
#############################################################################################################

(1) USAGE
-----
SUPER-FOCUS: A tool for agile functional analysis of shotgun metagenomic data
--------------------------------------------------------------------------------------------------------------------------------
Options:
-h     ------: print help

-q     query file (FASTA or FASTQ format) or folder with multiple FASTA/FASTQ files when -m 1

-dir   string: output directory

-m     int:    run the program for multiple files – 0 (False) / 1 ( True) (default: 0)	 
 
-o     string: project name (default 'my_project')

-mi    float:  minimum identity (default 60 %)

-ml    int:    minimum alignment (amino acids) (default: 15)

-focus int:    runs FOCUS; 1 does run; 0 does not run: default 0

-t     int:    number of threads (default 8)

-e     float:  e-value (default 0.00001)

-db    string: database (DB_90, DB_95, DB_98, or DB_100; default DB_98)

-p     int:    amino acid input; 0 nucleotides; 1 amino acids (default 0)

-a     string: aligner choice (rapsearch, blast, diamond; default rapsearch)

-fast  int:    runs RAPSearch2 or DIAMOND on fast mode - 0 (False) / 1 (True) (default: 1)	
  
-n     int:    normalizes each query counts based on number of hits; 0 doesn't normalize; 1 normalizes (default: 1)

-r     string: use only the subsystems in the organisms predicted by "-focus"– ncbi / rast annotation  (default: ncbi)

--------------------------------------------------------------------------------------------------------------------------------
example> python superfocus.py -q query.fasta -dir myOutputdirectory
	 
(2) OUTPUT
SUPER-FOCUS output will be add the folder selected in -dir

(3) RECOMMENDATIONS
	- The FOCUS reduction is not necessary if not wanted (set -focus 0)
	- Run RAPSearch for short sequences. it is less sensitive for long sequences
	- How BLAST if you want the result to be the most sensitive as possible
	- Only use DIAMOND for large datasets. It is slower than blastx for small datasets

(4) DEPENDENCIES
------------
Python >= 2.6.X < 3.Y: http://www.python.org/download
Jellyfish: http://www.cbcb.umd.edu/software/jellyfish
Numpy: http://sourceforge.net/projects/numpy/files/NumPy
SciPy: http://sourceforge.net/projects/scipy

One of the below aligners:
RAPSearch2: http://rapsearch2.sourceforge.net
DIAMOND: http://ab.inf.uni-tuebingen.de/software/diamond
BLAST: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST

COPYRIGHT AND LICENSE
---------------------
Copyright (C) 2014-2015  Genivaldo Gueiros Z. Silva

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.