####################################################################################
# This bash script is a part of PhyloBench
# Run it after makegoodlist.py
# The optional parameter is the name of a file with the list of orthologous groups.
# If the parameter is omitted, then "goodlist.txt" is used.
####################################################################################

if [[ ! -d Trees ]]
  then mkdir Trees
fi
if [[ ! -d Trees/RAxML ]]
  then mkdir Trees/RAxML
fi
if [[ ! -d Trees/RAxML/logs ]]
  then mkdir Trees/RAxML/logs
fi
if [[ ! -d Trees/RAxML/parsimony ]]
  then mkdir Trees/RAxML/parsimony
fi
if [[ ! -d Alignments/reduced ]]
  then mkdir Alignments/reduced
fi
if [[ $# -eq 0 ]]
  then inlist=goodlist.txt
  else inlist=$1
fi 
for x in `cat $inlist` # !!! Parallelize if possible!
do
  seqret Alignments/$x.afa phylip::temp.phy -auto
  echo $x `date '+%D %H:%M:%S'`
  raxmlHPC -m PROTGAMMAAUTO -s temp.phy -p 12345 -n tmp > Trees/RAxML/logs/$x.info2
  mv RAxML_bestTree.tmp Trees/RAxML/$x.tre
  mv RAxML_parsimonyTree.tmp Trees/RAxML/parsimony/$x.tre
  rm RAxML_result.tmp
  mv RAxML_log.tmp Trees/RAxML/logs/$x.log
  mv RAxML_info.tmp Trees/RAxML/logs/$x.info
  rm temp.phy
  if [[ -f temp.phy.reduced ]]
    then mv temp.phy.reduced Alignments/reduced/$x.phy
  fi
done
