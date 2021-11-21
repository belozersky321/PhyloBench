#############################################################################
#    This script is written on Python 2.7
#    It is a part of PhyloBench. 
#    It requires subdirectories with taxonomic alignment sets
#    including subsubdirectories "Selection15", "Selection30" and "Selection45" with alignments.
#    It is assumed that short names "AR" and "AG" are for archeael alignments,
#    short names "AC", "FI", "OB" and "PB" are for bacterial alignments
#    and short names "AS", "CH", "EB", "FB", "MA", "ST" are for eukaryotic alignments
#    This script creates the subfolder "Combined" with combined alignments sets.
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

from random import sample, choice
import os

arcs = ["AR", "AG"]
euks = ["AS", "CH", "EB", "FB", "MA", "ST"]
bacs = ["AC", "FI", "OB", "PB"]

arc = {"15":set(), "30":set(), "45":set()}
euk = {"15":dict(), "30":dict(), "45":dict()}
bac = {"15":dict(), "30":dict(), "45":dict()}
for num in ["15", "30", "45"]:
  for fl in os.listdir("."):
    if os.path.isdir(fl):
      if "Selection" + num in os.listdir(fl):
        if "Alignments" in os.listdir(fl + "/Selection" + num):
          for x in os.listdir(fl + "/Selection{}/Alignments".format(num)):
            ttt = x.split('_')
            short = ttt[0]
            pf = ttt[1]
            if short in arcs:
              arc[num].add(x.split('.')[0])
            if short in bacs: 
              if pf in bac[num]: bac[num][pf].append(x.split('.')[0])
              else: bac[num][pf] = [x.split('.')[0]]
            if short in euks: 
              if pf in euk[num]: euk[num][pf].append(x.split('.')[0])
              else: euk[num][pf] = [x.split('.')[0]]

arcall = arc["15"] & arc["30"] & arc["45"]
print "{} archaeal alignmens".format(len(arcall))

eukall = dict()
eukset = set(euk["15"].keys()) & set(euk["30"].keys()) & set(euk["45"].keys())
for pf in eukset:
  pfset = set(euk["15"][pf]) & set(euk["30"][pf]) & set(euk["45"][pf])
  if pfset:
    eukall[pf] = list(pfset)
bacall = dict()
bacset = set(bac["15"].keys()) & set(bac["30"].keys()) & set(bac["45"].keys())
for pf in bacset:
  pfset = set(bac["15"][pf]) & set(bac["30"][pf]) & set(bac["45"][pf])
  if pfset:
    bacall[pf] = list(pfset)

eukset = set(eukall.keys())
bacset = set(bacall.keys())
oeukset = eukset - bacset
obacset = bacset - eukset
univset = eukset & bacset
print "{} only eukaryotic Pfam families".format(len(oeukset))
print "{} only bacterial Pfam families".format(len(obacset))
print "{} universal Pfam families".format(len(univset))

n = len(arcall)
if n > 0:
  n1 = min(n/2, len(oeukset), len(obacset))
  n2 = min(n - n1, len(univset))
  combueuk = sample(univset, n2)
  combubac = sample(univset, n2)
else:
  n1 = min(len(oeukset), len(obacset))
  combueuk = univset
  combubac = univset

comboeuk = sample(oeukset, n1)
combobac = sample(obacset, n1)

if not os.path.isdir("Combined"): 
  os.mkdir("Combined")
for num in ["15", "30", "45"]:
  if not os.path.isdir("Combined/" + num): 
    os.mkdir("Combined/" + num)
  if not os.path.isdir("Combined/" + num + "/Alignments"): 
    os.mkdir("Combined/" + num + "/Alignments")

if n > 0:
  if n1 + n2 < len(arcall):
    arclist = sample(arcall, n1 + n2)
  else:
    arclist = list(arcall)
  for y in arclist:
    for num in ["15", "30", "45"]:
      os.system("cp */Selection{}/Alignments/{}.afa Combined/{}/Alignments/".format(num, y, num))

for x in (comboeuk + combueuk):
  taxa = dict()
  for a in eukall[x]:
    ttt = a.split('_')
    s = ttt[0]
    if s in taxa:
      taxa[s].append(a)
    else:
      taxa[s] = [a]
  o = choice(taxa.keys())
  y = choice(taxa[o])
  for num in ["15", "30", "45"]:
    os.system("cp */Selection{}/Alignments/{}.afa Combined/{}/Alignments/".format(num, y, num))
print "{} eukaryotic alignments".format(len(comboeuk) + len(combueuk))

for x in (combobac + combubac):
  taxa = dict()
  for a in bacall[x]:
    ttt = a.split('_')
    s = ttt[0]
    if s in taxa:
      taxa[s].append(a)
    else:
      taxa[s] = [a]
  o = choice(taxa.keys())
  y = choice(taxa[o])
  for num in ["15", "30", "45"]:
    os.system("cp */Selection{}/Alignments/{}.afa Combined/{}/Alignments/".format(num, y, num))
print "{} bacterial alignments".format(len(combobac) + len(combubac))
# end
