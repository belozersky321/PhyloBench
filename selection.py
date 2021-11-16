#############################################################################
#    This script is written on Python 2.7
#    It is a part of PhyloBench. 
#    It makes subalignments of the alignments of benchmark.
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

import random
import sys
import os
import math
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

###################################################
# main
###################################################
python = "python2.7" # Python 2.7 interpreter
threshold = 9        # threshold for number of good columns to regard alignment as good

speciesnumber = int(sys.argv[1])   # Number of sequences in subalignments
outdir = "Selection" + sys.argv[1]
if not os.access(outdir, os.F_OK):
  os.mkdir(outdir)
if not os.access(outdir + "/Sequences", os.F_OK):
  os.mkdir(outdir + "/Sequences")
if not os.access(outdir + "/Alignments", os.F_OK):
  os.mkdir(outdir + "/Alignments")

alidir = "Alignments/"
speciesflow = open("list.txt", "r")
familiesflow = open("goodlist.txt", "r")

specieslist = list()
for x in speciesflow:
  species = x.strip()
  if len(species) > 0:
    specieslist.append(species)
speciesflow.close()

familiestuple = ()
for x in familiesflow:
  family = x.strip()
  if len(family) > 0:
    familiestuple += (family,)
familiesflow.close()

successes = 0
faults = 0
random.seed()
for family in familiestuple:
  ok = True
  if os.access(alidir + "reduced/" + family + ".phy", os.F_OK):
    redali = open(alidir + "reduced/" + family + ".phy", "r")
    n = int(redali.readline().split()[0])
    if n < speciesnumber: 
      ok = False
    else: 
      os.system("seqret " + alidir + "reduced/" + family + ".phy fasta::tmpred.afa -auto")
      inali = "tmpred.afa"
      specieslist = list()
      for k in range(n):
        specieslist.append(redali.readline().split()[0])
    redali.close()
  else: # no reduced alignment
    inali = alidir + family + ".afa"
  if ok: 
    flag = False
    outseq = outdir + "/Sequences/" + family + ".fasta"
    outali = outdir + "/Alignments/" + family + ".afa"
    popytka = 0
    while not flag and popytka < 13:
      popytka += 1
      specsublist = random.sample(specieslist, speciesnumber)
      outflow = open(outseq, "w")
      for x in specsublist:
        command = ["degapseq", alidir + family + ".afa:" + x, "-filter"]
        degapseq = Popen(command, stdout = PIPE, stderr = PIPE)
        (out, err) = degapseq.communicate()
        if err: print err
        outflow.write(out)
      # for x
      outflow.close()
      command = ["muscle", "-in", outseq, "-out", "temp.afa"]
      muscle = Popen(command, stdout = PIPE, stderr = PIPE)
      (out, err) = muscle.communicate()
      print family, out
#      command = [python, "../Scripts/infoalign.py", "temp.afa"]
#      infoalign = Popen(command, stdout = PIPE, stderr = PIPE)
#      (out, err) = infoalign.communicate()
#      if err: print err
#      (good, conserved, gapped) = tuple(out.split('\t'))
      (good, conserved, gapped) = infoalign("temp.afa")
      if int(good) > threshold: 
        flag = True
        os.system("{} ../Scripts/stable.py {} temp.afa > {}".format(python, outseq, outali))
    os.remove("temp.afa")
    if os.access("tmpred.afa", os.F_OK): os.remove("tmpred.afa")
    if flag:
      successes += 1
    else:
      faults += 1
      print "For", family, "cannot select a good subalignment"
  else:
    print "For", family, "there are no", speciesnumber, "different sequences"
# for i
print successes, "successes; ", faults, "faults"
# end
