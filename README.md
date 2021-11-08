# PhyloBench
Python 2.7 and bash scripts for a benchmark for phylogenetic programs

This repository is a supplement to the paper describing PhyloBench, a benchmark for testing phylogenetic programs. It includes the following scripts:
 * selectmnems.py
 * swisspfam_reduce.py
 * maketable.py
 * extractseq.py
 * blast_ort.py
 * muscle.sh
 * makegoodlist.py
 * infoalign.py
 * fastme.py
 * raxml.sh
 * tnt.sh
 * tnt_unput.txt
 * convert.py
 * consense-ncbi.py
 * newick.py
 * selection.py
 * pfs-c.py

-----
The procedure of generating files consisting PhyloBench is as follows. 

Files from Internet
-----
 * `names.dmp` and `nodes.dmp` from the archive https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.tar.gz
 * `swisspfam` from http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/swisspfam.gz
 * `pfamseq` from http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/pfamseq.gz
 * `speclist.txt`: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/speclist.txt
 * `stable.py` from http://www.drive5.com/muscle/stable.tar.gz
 

Prerequisites
-----
 * Python 2.7
 * EMBOSS
 * BLAST+
 * FastME
 * RAxML
 * TNT
 * Muscle

Preliminary actions
-----
Put `swisspfam`, `nodes.dmp`, `names.dmp` and `speclist.txt` to the working directory.

Run `swisspfam_reduce.py`. The result is two files: `swisspfam_reduced.txt` and `swisspfam_coord.txt`.

Run `makefile.py`. The result is the file `table.txt`.

Index `pfamseq` as EMBOSS sequence database (program `dbxfasta`).

To create files for each species set, make a directory for this species set and perform the following steps in this directory. It is assuming that files `swisspfam`, `nodes.dmp`, `names.dmp`, `speclist.txt`, `swisspfam_reduced.txt`, `swisspfam_coord.txt`, `table.txt` are in parent directory during performing all these steps.

Step 1: Selection of species
-----
Selection of species can be made with one of three procedures. All three procedures requires files `names.dmp`, `table.txt` and `swisspfam_reduced.txt` in the parent folder.

To select species set with the Procedure 1, choose an "embracing" taxon (e.g., Metazoa), a taxon level to select one species from each taxon of this level (e.g., family), and a number of species to select (e.g., 60). Create a directory for the new species set and make it active. Run `selectmnems.py` specifying the embracing taxon, the taxon level and the number, for example:

`python2.7 selectmnems.py tax=Metazoa level=family n=60`

The result is three files: `list.txt` with selected organism mnemonics, `pfamlist.txt` with accession numbers of selected Pfam families and `<level>.txt` with names of selected taxons of the chosen level (for the example above the third file name would be `family.txt`). If Pfam contains less than n taxons of the chosen level in the embracing taxon, then the files `list.txt` and `<level>.txt` will contain less than n items.

To select species with the Procedure 2, choose several embracing taxons and for each a taxon level and a number. This procedure is a generalization of the Procedure 1 and consists of repeating Procedure 1, at each repeat using a set of Pfam domains restricted by the previous repeat. Example: if taxons are Nematoda, Arthropoda and Chordata, levels are genus, family and order and numbers are 10, 25, and 25, then the procedure consists of the following chain of commands:

`python2.7 selectmnems.py tax=Nematoda level=genus n=10`

mv pfamlist.txt pfamlist1.txt

mv list.txt list1.txt

`python2.7 selectmnems.py tax=Arthropoda level=family n=25 pfam=pfamlist1.txt`

mv pfamlist.txt pfamlist2.txt

mv list.txt list2.txt

`python2.7 selectmnems.py tax=Chordata level=order n=25 pfam=pfamlist2.txt`

The result

Step 2: extracting sequences of domains
-----


Run 
