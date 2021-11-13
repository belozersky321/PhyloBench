####################################################
# This bash script is a part of PhyloBench
# Run it after blast_ort.py
####################################################

python=/usr/bin/python2.7 # !!! Change this to the actual path to Python 2.7 interpreter
if [[ ! -d Alignments ]]
  then mkdir Alignments
fi
for x in `ls Sequences | sed 's/\.fasta$//'`
do
  echo $x
  muscle -in Sequences/$x.fasta -out temp.afa 2>> muscle.log
  $python ../Scripts/stable.py Sequences/$x.fasta temp.afa > Alignments/$x.afa
  rm temp.afa
done

