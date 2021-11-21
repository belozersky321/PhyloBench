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
 * fastme.py
 * tnt.sh
 * convert.py
 * raxml.sh
 * consense-ncbi.py
 * selection.py
 * pfs-c.py

and the file "tnt-input.txt" needed for "tnt.sh".

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

Run `maketable.py`. The result is the file `table.txt`.

Index `pfamseq` as EMBOSS sequence database named "pfamseq" (program `dbxfasta` of the EMBOSS package).
In the original file "pfamseq", ID and AC of sequences are in non-standard order (first AC, then ID).
Thus we changed their order with the command:

`sed 's/>\([A-Z0-9\.]*\) \([A-Z0-9]*_[A-Z0-9]*\) />\2 \1 /' pfamseq > pfamseq.fasta`

Now, the indexing command is as follows:

`dbxfasta -dbname pfamseq -dbresource pfamseq -idformat idacc -directory $pfamdir -filenames pfamseq.fasta -fields id,acc -compressed -outfile pfamseq.dbxfasta -release 33.1 -date $today`

(here "$pfamdir" is for the folder with "pfamseq.fasta" and "$today" is for the current date).
The script `extractseq.py` refers to "pfamseq-id" assuming that pfamseq is indexed this way.
The file emboss.default should contain the following lines:

    DB pfamseq [
        type: P
        method: emboss
        format: fasta
        directory: <specify here the directory with "pfamseq.fasta">
        file: pfamseq.fasta
        indexdirectory: <specify here the directory were dbxfasta was run>
        release: "33.1"
        fields: "id acc"
        comment: "Pfam sequences from pfamseq release 33.1"
    ]

Make the directory named "Scripts" and put all files from this repository (scripts and `tnt_input.txt`) there.

To create files for each species set, make a directory for this species set 
and perform the following steps in this directory. 
It is assuming that files `swisspfam`, `nodes.dmp`, `names.dmp`, `speclist.txt` and `table.txt` 
are in the parent directory and the directory "Scripts" is a sister directory 
during performing all these steps.

Step 1: Selection of species
-----
Selection of species can be made with one of three procedures. 
All three procedures require files `names.dmp`, `table.txt` and `swisspfam` 
in the parent folder.

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
Pfam accessions listed in the file `pfamlist.txt`. 
Each subdirectory contains files with names consisting of organism mnemonics 
listed in the file `list.txt`, with extension ".fasta". 
Each file contains sequences of all domains of the certain Pfam family from the certain organism, 
in fasta format.

Step 3: creating alignments of orthologous groups
-----
Choose some prefix for the files with sequences and alignments of orthologous groups of this set.

Run `blast_ort.py`. 
This script creates the directory "Sequences" with fasta files containing sequences of found orthologous groups.
It takes one argument, which is a short name of the taxonomic set. 
For example, for Actinobacteria you may choose the short name "AC" and run:

`python2.7 ../Scripts/blast_ort.py AC > blast_ort.log`

Copy the file `stable.py` to the directory "Scripts". 
If necessary, change the first line of `muscle.sh` (specify the actual path to the Python 2.7 interpreter).
Run `muscle.sh`:

`bash ../Scripts/muscle.sh`

This script creates the directory "Alignments" with fasta files containing
alignments of found orthologous groups. 
The extension of the names of these files is ".afa".

Step 4: building the reference tree
-----
Run `makegoodlist.py` to create the files `goodlist.txt` and `badlist.txt`.

Run `fastme.py` to infer trees from orthologous groups with the program FastME.

Run the bash script  `tnt.sh` to infer trees from orthologous groups with the program TNT.
This script requires the files `tnt_input.txt`and `convert.py` in the sister
directory "Scripts".

Run the bash script `raxml.sh` to infer trees from orthologous groups with the program RAxML. 
Because RAxML is very slow, the execution may take several hours or even days (run in background with nohup).
A good idea would be to parallelize this script by spliting the file "goodlist.txt"
to several sublists and running commands of this script in parallel with different sublists as input.

After all trees are inferred, run the script `consense-ncbi.py` to create 
two files: `ncbi.tre` and `ncbi-consense.tre`. Both contains trees in Newick format,
the first one is the (usually unresolved) tree from the NCBI taxonomy, 
the second one is the (usually fully resolved) tree obtained from the first one
by adding branches form the inferred trees.

Step 5: create subalignments
-----
Run `selection.py` three times to generate subalignments of 15, 30 and 45 sequences each.

    python2.7 ../Scripts selection.py 15
    python2.7 ../Scripts selection.py 30
    python2.7 ../Scripts selection.py 45

Combined sets
-----
The above steps 1–5 create alignments of one species set. 
To create Combined sets of alignments of 15, 30 and 45 each, the script `pfs-c.py` was used. 
In this script, the lists of archaeal, bacterial and eukaryotic species sets and their short names 
are fixed according to the current version of PhyloBench.
Namely,  archaeal sets are "AG" and "AR", bacterial sets are "AC", "FI", "OB" and "PB", and eukaryotic
sets are "AS", "CH", "EB", "FB", "MA" and "ST".
Edit this script to create Combined sets from different taxonimic sets.

The script `pfs-c.py` must be run from the main working directory, which is parent to all 
directories of the species sets.

`python2.7 Scripts/pfs-c.py`

The result is the folder "Combined" with three subfolders, with 15- 30- and 45-sequence alignments.
