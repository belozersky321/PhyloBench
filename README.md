# PhyloBench
Scripts for a benchmark for phylogenetic programs

This repository is a part of Supplementary material to a paper describing PhyloBench, 
a benchmark for testing phylogenetic programs. It includes the following scripts:
 * maketable.py
 * selectmnems.py
 * extractseq.py
 * blast_ort.py
 * muscle.sh
 * makegoodlist.py
 * infoalign.py
 * fastme.py
 * tnt.sh
 * tnt_unput.txt
 * convert.py
 * raxml.sh
 * consense-ncbi.py
 * newick.py
 * selection.py
 * pfs-c.py

-----

Here the procedure of generating files consisting PhyloBench is described. 

Files needed from Internet 
-----
 * `names.dmp` and `nodes.dmp` from the archive https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.tar.gz
 * `swisspfam` from http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/swisspfam.gz
 * `pfamseq` from http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/pfamseq.gz
 * `speclist.txt`: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/speclist.txt
 * `stable.py` from http://www.drive5.com/muscle/stable.tar.gz

Programs and packages needed to be installed
-----
 * Python 2.7
 * EMBOSS
 * BLAST+
 * FastME
 * RAxML
 * TNT
 * Muscle
 * PHYLIP

Preliminary actions
-----
Put `swisspfam`, `nodes.dmp`, `names.dmp` and `speclist.txt` to the working directory.

Run `swisspfam_reduce.py`. The result is two files: `swisspfam_reduced.txt` and `swisspfam_coord.txt`.

Run `maketable.py`. The result is the file `table.txt`.

Index `pfamseq` as EMBOSS sequence database named "pfamseq" (program `dbxfasta` of the EMBOSS package).

Make the directory named "Scripts" and put all files from this repository (scripts and `tnt_input.txt`) there.

To create files for each species set, make a directory for this species set 
and perform the following steps in this directory. 
It is assuming that files `swisspfam`, `nodes.dmp`, `names.dmp`, `speclist.txt`, 
`swisspfam_reduced.txt`, `swisspfam_coord.txt`, `table.txt` 
are in the parent directory and the directory Scripts is a sister directory during performing all these steps.

Step 1: Selection of species
-----
Selection of species can be made with one of three procedures. 
All three procedures require files `names.dmp`, `table.txt` and `swisspfam_reduced.txt` in the parent folder.

To select species set with the Procedure 1, choose an "embracing" taxon (e.g., Metazoa), 
a taxon level to select one species from each taxon of this level (e.g., family), 
and a number of species to select (e.g., 60). 
Create a directory for the new species set and make it active. 
Run `selectmnems.py` specifying the embracing taxon, the taxon level and the number, for example:

`python2.7 ../Scripts/selectmnems.py tax=Metazoa level=family n=60`

The result is three files: `list.txt` with selected organism mnemonics, 
`pfamlist.txt` with accession numbers of selected Pfam families and 
`<level>.txt` with names of selected taxons of the chosen level 
(for the example above the third file name would be `family.txt`). 
If Pfam contains less than *n* taxons of the chosen level in the embracing taxon, 
then the files `list.txt` and `<level>.txt` will contain less than *n* items.

To select species with the Procedure 2, 
choose several embracing taxons and for each a taxon level and a number. 
This procedure is a generalization of the Procedure 1 and consists of repeating Procedure 1, 
at each repeat using a set of Pfam domains restricted by the previous repeat. 
Example: if taxons are Nematoda, Arthropoda and Chordata, levels are genus, family and order 
and numbers are 10, 25, and 25, then the procedure consists of the following chain of commands:

    python2.7 ../Scripts/selectmnems.py tax=Nematoda level=genus n=10
    mv pfamlist.txt pfamlist1.txt  
    mv list.txt list1.txt
    python2.7 ../Scripts/selectmnems.py tax=Arthropoda level=family n=25 pfam=pfamlist1.txt
    mv pfamlist.txt pfamlist2.txt
    mv list.txt list2.txt
    python2.7 ../Scripts/selectmnems.py tax=Chordata level=order n=25 pfam=pfamlist2.txt
    cat list1.txt list2.txt >> list.txt

Procedure 3 differs from the Procedure 2 by excluding some subtaxons from the embracing taxon. Example:

    python2.7 ../Scripts/selectmnems.py tax=Metazoa exc=Chordata,Arthropoda,Nematoda level=genus n=20 pfam=pfamlist3.txt

Step 2: extracting sequences of domains
-----
Run `extractseq.py`, which creates the directory "Domains" with subdirectories named by 
Pfam accessions from the file `pfamlist.txt`. 
Each subdirectory contain files with names consisting of organism mnemonics 
from the file `list.txt` with extension ".fasta". 
Each file contains sequences of all domains of the certain Pfam family from the certain organism, 
in fasta format.

Step 3: creating alignments of orthologous groups
-----
Choose some prefix for the files with sequences and alignments of orthologous groups of this set.

Run `blast_ort.py`. 
This script creates the directory "Sequences" with fasta files containing sequences of found orthologous groups.

Copy the file `stable.py` to the current directory and run `muscle.sh`.
This script creates the directory "Alignments" with fasta files containing
alignments of sound orthologous groups. 
The extension of the names of these files is ".afa".

Step 4: building the reference tree
-----
Run `makegoodlist.py` to create the files `goodlist.txt` and `badlist.txt`.

Run `fastme.py` to infer trees from orthologous groups with the program FastME.

Run `tnt.sh` to infer trees from orthologous groups with the program TNT.
This script requires the files `tnt_input.txt`and `convert.py` in the sister
directory "Scripts".

Run `raxml.sh` to infer trees from orthologous groups with the program RAxML. 
Because RAxML is very slow, the execution may take several hours (run in background with nohup).

...

Step 5: create subalignments
-----
Run `selection.py`.

...

Combined sets
-----
The above steps 1–5 create alignments of one species set. 
To create Combined sets of alignments of 15, 30 and 45 each, the script `pfs-c.py` was used. 
In this script, it a number of parameters of the current version of Phylobench are fixed. 
These parameters are:
 * the lists of archaeal, bacterial and eukaryotic species sets
 * the numbers of "only bacterial", "only eukaryotic" and "universal" Pfam families to select, each is 325
Edit this script to change these settings.