#############################################################################
#    This script is written on Python 2.7
#    It is a part of PhyloBench. 
#    It requires file "swisspfam" from Pfam
#    and files "list.txt" and "pfamlist.txt" created by the script "selectmnems.py".
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

# pfamseq should be indexed for EMBOSS
from os import system, mkdir, access, F_OK
from time import ctime

def ispf(x):
  if len(x) > 8:
    if len(x.split('.')) == 2 and len(x.split('.')[0]) == 7 and x[:2] == "PF":
      return True
    else:
      return False
  else:
    return False
# ispf

speclist = "list.txt"
pfamlist = "pfamlist.txt"
swisspfam = "../swisspfam"

pfamset = set()
pfamin = open(pfamlist, "r")
for line in pfamin:
  if line.strip():
    pfamset.add(line.strip())
pfamin.close()

specset = set()
specin = open(speclist, "r")
for line in specin:
  if line.strip():
    specset.add(line.strip())
specin.close()

if not access("Domains", F_OK):
  mkdir("Domains")
coord = open(swisspfam, "r")
k = 0
d = 0
p = 0
for line in coord:
  k += 1
  ell = line.strip().split()
  if ell:
    if ell[0][0] == '>':
      seqid = ell[0][1:]
    else:
      if seqid.split('_')[1] in specset:
        for z in ell:
          if ispf(z):
            pfamid = z
            if pfamid in pfamset:
              d += 1
              begin = ell[-1].split('-')[0]
              end = ell[-1].split('-')[1]
              pfamdir = "Domains/" + pfamid
              if not access(pfamdir, F_OK):
                p += 1
                mkdir(pfamdir)
              seqfile = pfamdir + "/" + seqid.split("_")[1] + ".fasta"
              name = seqid + "/" + begin + "-" + end
              if not access(seqfile, F_OK):
                system("descseq pfamseq-id:" + seqid + "[" + begin + ":" + end + "] " + seqfile + " -name " + name + " -auto")
              else:
                system("descseq pfamseq-id:" + seqid + "[" + begin + ":" + end + "] -name " + name + " -filter >> " + seqfile)
  if k % 1000000 == 0:
    print ctime(), "{} lines, {} domains, {} pfam families".format(k, d, p)
coord.close()
# end
