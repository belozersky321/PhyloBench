#############################################################################
#    This script is written on Python 2.7
#    It is a part of PhyloBench.
#    It created two trees in Newick format. 
#    The tree in file "ncbi.tre" reflects NCBI taxonomy of species
#    whose mnemoniucs are in "list.txt".
#    The tree in file "ncbi-consense.tre" is obtained from "ncbi.tre"
#    by adding branches from trees inferred from individual orthologous groups. 
#    This script requires the files "nodes.dmp" from NCBI Taxonomy database,
#    "speclist.txt" from Uniprot (these two in the parent folder),
#    "list.txt" created by "selectmnems.py", "goodlist.txt" created by "makegoodlist.txt"
#    (in the current folder) and folders "Trees/FastME", "Trees/TNT" and "Trees/RAxML"
#    with trees inferred from orthologous groups listed in "goodlist.txt".
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

from sys import argv, stderr
from subprocess import Popen, PIPE
import os

def extractnames(newick):
  name = ""
  names = frozenset()
  mode = True
  for i in range(len(newick)):
    if newick[i] == ':': mode = False
    if newick[i] in "(),":
      mode = True
      if name:
        names = names | frozenset([name])
      name = ""
    if mode and (newick[i] not in "(),:;") and not newick[i].isspace():
      name += newick[i]
  if name: names = names | frozenset([name])
  return names
# extractnames

def dotsnstars(branch):
  result = ""
  left = branch[0]
  right = branch[1]
  listnames = sorted(list(left | right))
  for name in listnames:
    if name in left: result += '.'
    if name in right: result += '*'
  return result
# dotsnstars

def reverse(branch):
  return (branch[1], branch[0])

def readnewick(newick):
  result = frozenset()
  allnames = extractnames(newick)
  for name in allnames:
    left = frozenset([name])
    right = allnames - left
    result = result | frozenset([(left, right)]) # add trivial branches
  stack = list()
  for i in range(len(newick)):
    if newick[i] == '(':
      stack.append(i)
    if newick[i] == ')':
      j = stack.pop()
      left = extractnames(newick[j:i]) # all names in the pair of brackets
      right = allnames - left
      result = result | frozenset([(left, right)]) # non-trivial branch
  return result
# readnewick

def writenewick(branches, root):
  # reverse some branches for consistency
  tmplist = list()
  for branch in branches:
    if branch != root:
      if branch[1].issubset(root[0]) or branch[1].issubset(root[1]): tmplist.append(branch)
  for branch in tmplist:
    branches = branches | frozenset([reverse(branch)])
    branches = branches - frozenset([branch])

  # add reverse copy of the root
  root2 = reverse(root)
  branches = branches | frozenset([root2])

  # create tree as a dictionary with parent nodes as values  
  clades = set()
  tree = dict()
  rc = root[0] | root[1]  # the root clade
  clades.add(rc) 
  for (left, right) in branches:
    clade = frozenset(left)
    if clade not in clades: 
      for cl2 in clades:
        if clade.issubset(cl2):
          if clade in tree:
            if tree[clade].issuperset(cl2): tree[clade] = cl2
          else:
            tree[clade] = cl2
        if clade.issuperset(cl2):
          if cl2 in tree:
            if tree[cl2].issuperset(clade): tree[cl2] = clade
          else:
            tree[cl2] = clade
      clades.add(clade)

  # forming Newick string as result
  result = ";"
  k = 0
  stack = list()  # stack of clades
  stack.append((rc, False)) # "False" means no comma after it in Newick
  tree[rc] = None
  while stack:
    (cl, flag) = stack.pop()
    if len(cl) == 1: # species
      species = list(cl)[0]
      result = result[:k] + species + result[k:]
      k += len(species)
      if flag:
        result = result[:k] + ', ' + result[k:]
        k += 2
      else:
        while result[k] == ')': k += 1
        if result[k] == ',': k += 2
    else: # clade
      if flag:
        result = result[:k] + "(), " + result[k:]
      else:
        result = result[:k] + "()" + result[k:]
      k += 1
      tmplist = list()
      s = False      # the first clade into the stack = the least out: no comma after it
      for clade in clades:
        if tree[clade] == cl: 
          stack.append((clade, s))
          s = True   # there should be comma after each clade except of the first 
  return result
# writenewick

def conc(branch1, branch2):
  result = False
  if branch1[0].isdisjoint(branch2[0]): result = True
  if branch1[1].isdisjoint(branch2[0]): result = True
  if branch1[0].isdisjoint(branch2[1]): result = True
  if branch1[1].isdisjoint(branch2[1]): result = True
  return result
# conc

def readcons(filename):
  species = list()
  branchlist = list()
  inf = open(filename, "r")
  line = inf.readline()
  while not line.startswith("Sets"):
    line = inf.readline()
    ttt = line.split()
    if len(ttt) == 2:
      if ttt[0][:-1].isdigit(): species.append(ttt[1].strip())
  while not line.startswith("Extended"):
    line = inf.readline()
    if line.startswith('.') or line.startswith('*'):
      lefts = set()
      rights = set() 
      i = 0
      j = 0
      while(j < len(species)):
        if line[i] == '.':
          lefts.add(species[j])
          j += 1
        if line[i] == '*':
          rights.add(species[j])
          j += 1
        i += 1
      branchlist.append((frozenset(lefts), frozenset(rights)))
  inf.close()
  return branchlist
# readcons

###############################
# main
###############################
treedirs = dict()
treedirs["bionj"] = "Trees/FastME"
treedirs["mp"] = "Trees/TNT"
treedirs["ml"] = "Trees/RAxML"

nodes = open("../nodes.dmp", "r")
parent = dict()
for x in nodes:
  ell = x.split("\t|\t")
  if len(ell) >= 3:
    parent[ell[0]] = ell[1]
nodes.close()
stderr.write("Parents OK, {} nodes\n".format(len(parent.keys())))

specset = set()
inf = open("list.txt", "r")
for line in inf:
  spec = line.strip()
  if spec: specset.add(spec)
inf.close()

taxid = dict()
taxonset = set()
specs = open("../speclist.txt", "r")
x = specs.readline()
while not x.startswith("(1) Real organism codes"):
  x = specs.readline()
while not x.startswith("_____"):
  x = specs.readline()
while not x.startswith("======================================================================="):
  x = specs.readline()
  if "N=" in x:
    y = x.split()
    mnem = y[0]
    if mnem in specset:
      taxid[mnem] = y[2][:-1]
      taxonset.add(y[2][:-1])
specs.close()
stderr.write("{} taxids for the species\n".format(len(taxonset)))

alltaxons = dict()
for mnem in taxid:
  t = taxid[mnem]
  alltaxons[mnem] = {t}
  while t in parent and t != parent[t]:
    t = parent[t]
    alltaxons[mnem].add(t)

tottaxons = set()
for mnem in specset:
  tottaxons |= alltaxons[mnem]
stderr.write("{} taxons\n".format(len(tottaxons)))

branches = frozenset()
root = None
for t in tottaxons:
  left = frozenset()
  for mnem in specset:
    if t in alltaxons[mnem]:
      left = left | frozenset([mnem])
  right = frozenset(tottaxons - left)
  if len(left) < len(specset): 
    right = frozenset(specset - left)
    branches = branches | frozenset([(left, right)])
    if root:
      if len(left) > len(root[0]): root = (left, right) # find root branch as the branch with maximum left side
    else:
      root = (left, right)
stderr.write("{} initial branches\n".format(len(branches)))

# output
outf = open("ncbi.tre", "w")
outf.write(writenewick(branches, root) + '\n')
outf.close()
stderr.write("NCBI tree written to \"ncbi.tre\"\n")

for phylipfile in ["intree", "outtree", "outfile"]:
  if os.access(phylipfile, os.F_OK): os.remove(phylipfile)
inf = open("goodlist.txt")
for line in inf:
  og = line.strip()
  if og:
   for method in treedirs:
     os.system("cat {}/{}.tre >> intree".format(treedirs[method], og)) 
consense = Popen("consense", stdin = PIPE, stdout = PIPE, stderr = PIPE)
(out, err) = consense.communicate("y\n")
if err: stderr.write(err)
os.remove("intree")
os.rename("outtree", "consense.tre")
os.rename("outfile", "consense.txt")

consbranchlist = readcons("consense.txt")
stderr.write("{} branches in consense.txt\n".format(len(consbranchlist)))
n = 0
for br1 in consbranchlist:
  if br1 not in branches and reverse(br1) not in branches:
    flag = True
    for br2 in branches:
      if not conc(br1, br2): flag = False
    if flag: 
      branches = branches | frozenset([br1])
      n += 1
stderr.write("{} branches added to the NCBI tree\n".format(n))

# output
outf = open("ncbi-consense.tre", "w")
outf.write(writenewick(branches, root) + '\n')
outf.close()
stderr.write("Supplemented NCBI tree written to \"ncbi-consense.tre\"\n")
# end
