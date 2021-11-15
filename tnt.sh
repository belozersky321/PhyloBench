####################################################
# This bash script is a part of PhyloBench
# Run it after makegoodlist.py
# It requires the file "tnt_input.txt"
# and the Python script "convert.py"
####################################################

if [[ ! -d Trees ]]
  then mkdir Trees
fi
if [[ ! -d Trees/TNT ]]
  then mkdir Trees/TNT
fi
if [[ ! -d Trees/TNT/logs ]]
  then mkdir Trees/TNT/logs
fi
for x in `cat goodlist.txt`
do
  echo $x `date '+%D %H:%M:%S'` 
  seqret Alignments/$x.afa nexus::infile.nex -auto
  tnt < ../Scripts/tnt_input.txt >> tnt.tmp 2> tnt.stderr
  python2.7 ../Scripts/convert.py
  mv newick.tre Trees/TNT/$x.tre
  mv tnt.out Trees/TNT/logs/$x.tnt
  mv tnt.stderr Trees/TNT/logs/$x.log
  rm infile.nex temp.tre
done
