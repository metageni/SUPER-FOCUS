# Scaffold_builder
 Combining de novo and reference-guided assembly with Scaffold_builder
-----

(c)            Silva GG, Dutilh BE, Matthews TD, Elkins K, Schmieder R, Dinsdale EA, Edwards RA
Please cite:   "Combining de novo and reference-guided assembly with Scaffold_builder", Source Code for Biology and Medicine 2013

WEBSERVER
-----
http://edwards.sdsu.edu/scaffold_builder/



USAGE
-----
python scaffold_builder.py -q query_contigs.fna -r reference_genome.fna -p output_prefix [-t] [-i] [-a] [-b]

	-q fasta file of contigs
		Required. Query contigs in Fasta format. These contigs may be the output of a de novo
		assembly program such as Newbler, Velvet or MIRA.

	-r fasta file containing reference genome
		Required. Reference genome in Fasta format. This should preferably be a completed genome
		sequence.

	-p prefix output files
		Required. All the output files have this project name as prefix.

	-t length of terminus that will be aligned (default 300 nt)
		At any break between two contigs, scaffold_builder checks whether the termini
		of the adjacent contigs are homologous by aligning them using Smith-Waterman's Algorithm, and
		combines them if that is the case.

	-i minimum identity for merging contigs (default 80%)
		If the termini are similar, scaffold_builder assumes that the contigs should
		have been combined by the assembly program, but the similarity was probably
		below the assembly thresholds, or the contigs were not merged due to ambiguous
		read mapping. The sequences are combined and in the case that non-identical
		nucleotides are aligned, the IUPAC code of their consensus is placed in the
		resulting sequence.

	-a minimum length for ambiguously mapped contigs (default 95%)
		If a contig maps to more than one location on the reference genome, it will
		not be scaffolded because it's location is ambiguous. This parameter defines
		how much of the length of a contig should be mapped in more than one location
		for it to be considered ambiguously mapped.

	-b 0/1 dictates behavior for rearrangements (default 0)
		0: place end-to-end
		1: create new scaffold sequence
		If the mapping of the contigs onto the reference suggests that they overlap,
		but the contig termini are too dissimilar to join them, this option dictates
		whether scaffold_builder places the contigs end-to-end (default; deletions
		expected) or to start a new scaffold sequence (inversions expected).

	-g maximum gap length allowed (default 5000nt)



DEPENDENCIES
------------
python 2.7.X or 2.9.Y

COPYRIGHT AND LICENSE
---------------------
Copyright (C) 2011-2013  Genivaldo GZ Silva

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.
