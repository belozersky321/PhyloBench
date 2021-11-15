#############################################################################
#    This script is written on Python 2.7
#    It is a part of PhyloBench. 
#    It requires directory "Alignment" created by "muscle.sh".
#    The command line argument is a threshold for the number
#    of good columns to include the alignment into "goodlist.txt".
# 
#    Copyright (C) 2021 Sergey Spirin
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#    Contact: Sergey sas@belozersky.msu.ru
############################################################################

from os import listdir
from sys import argv
from subprocess import Popen, PIPE

class Sequence:
  def __init__(self, seqtype="U", name="Sequence"):
    self.seqtype = seqtype
    self.name = name
    self.des = ""
    self.seq = ""

  def __str__(self):
    """
    This function determines the format of printing for elements of Sequence class
    """
    _st = list()
    _st.append(">" + self.name + " " + self.des)
    for _i in range(0, len(self.seq)/60): 
      _st.append(self.seq[60 * _i : 60 * (_i + 1)])
    if len(self.seq)%60 != 0:
      _st.append(self.seq[-(len(self.seq)%60):])
    return "\n".join(_st)
# class Sequence

def readfasta(fastafile):
    '''readfasta returns the list of sequences from the open file in fasta format'''
    fastafile.seek(0)
    result = list()
    seqnum = -1
    for line in fastafile:
      if len(line.strip()) > 0:
        if line[0] == ">":
          newseq = Sequence("P", "")
          seqnum += 1
          tmp = line.strip().split(" ")
          if len(tmp[0]) == 1:
            newseq.name = " "
            print "Sequence with empty name!"
          else:
            newseq.name = tmp[0][1:]
          result += [newseq]
          if len(tmp) > 1 : newseq.des = " ".join(tmp[1:])
        if line[0] != ">" and seqnum >= 0:
          tmpline = line.replace(" ", "")
          newseq.seq += tmpline.strip()
          result[seqnum] = newseq
    return result
# readfasta

def infoalign(fastafile):
  '''
  Argument is the name of file with sequence alignment.
  Alignment should be in fasta format.
  Result is three numbers: good columns, conserved columns, gapped columns.
  "Good" means not conserved and without gaps.
  '''

  indata = open(fastafile, "r")
  inali = readfasta(indata)
  indata.close()

  lenali = len(inali[0].seq)
  numseq = len(inali)
  good = 0
  gapped = 0
  conserved = 0
  for i in range(lenali):
    letter = inali[0].seq[i]
    if (letter == '-'):
      isgapped = True
      isconserved = False
    else:
      isgapped = False
      isconserved = True
      j = 1
      while (j < numseq) and (not isgapped):
        if inali[j].seq[i] == "-":
          isgapped = True
        if inali[j].seq[i] != letter:
          isconserved = False
        j += 1
    if isconserved: conserved += 1
    if isgapped: gapped += 1
    if (not isconserved) and (not isgapped): good += 1
  # for
  return (good, conserved, gapped)
# infoalign

#############################
# main
#############################
try:
  threshold = int(argv[1])
except Exception:
  print "Threshold is incorrect or not specifyed, set to 9"
  threshold = 9

outdata = open("goodlist.txt", "w")
trashdata = open("badlist.txt", "w")

gn = 0
bn = 0
for filename in listdir("Alignments"):
  (good, conserved, gapped) = infoalign("Alignments/" + filename)
  if int(good) > threshold: 
    gn += 1
    outdata.write(filename.split(".")[0] + "\n")
  else: 
    bn += 1
    trashdata.write(filename.split(".")[0] + "\n")
trashdata.close()
outdata.close()
print gn, "good alignments", bn, "bad alignments"
# end