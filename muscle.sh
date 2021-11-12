if [[ ! -d Alignments ]]
  then mkdir Alignments
fi
for x in `ls Sequences | sed 's/\.fasta$//'`
do
  muscle -in Sequences/$x.fasta -out temp.afa
  python2.7 stable.py Sequences/$x.fasta temp.afa > Alignments/$x.afa
  rm temp.afa
done

