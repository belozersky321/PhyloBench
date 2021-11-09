#############################################################################
#    This script is written on Python 2.7
#    It is a part of PhyloBench. 
#    It requires files "nodes.dmp" from NCBI Taxonomy Database and "speclist.txt" from Uniprot.
#    It creates the file "table.txt", which is used by the script "selectmnems.py".
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

nodes = open("nodes.dmp", "r")
taxons = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom']

taxlevel = dict()
for taxon in taxons:
  taxlevel[taxon] = dict()
parent = dict()

for x in nodes:
#  ell = x.strip().split('\t')
  ell = x.split("\t|\t")
  if len(ell) >= 3:
    parent[ell[0]] = ell[1]
    for taxon in taxons:
      if ell[2] == taxon:
        taxlevel[taxon][ell[0]] = ell[0]
nodes.close()

print "Parents OK,", len(parent.keys()), "nodes"

n = 0
for s in parent:
  if s not in taxlevel['superkingdom']:
    n += 1

print n, "nodes without superkingdom"
m = n + 1
i = 0
while n < m:
  i += 1
  m = n
  n = 0
  for s in parent:
    t = parent[s]
    for taxon in taxons:
      if t in taxlevel[taxon]:
        taxlevel[taxon][s] = taxlevel[taxon][t]
    if s not in taxlevel['superkingdom']:
      n += 1
  print "Iter", i, n, "of", m, "nodes remains"
# while

print "Taxons OK"

specs = open("speclist.txt", "r")
table = open("table.txt", "w")
outline = ["Mnemonics", "taxid"] + taxons
table.write('\t'.join(outline) + '\n')

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
    taxid = y[2][:-1]
    outline = [mnem, taxid]
    for taxon in taxons:
      if taxid in taxlevel[taxon]:
        outline += [taxlevel[taxon][taxid]]
      else:
        outline += ["NA"]
    table.write('\t'.join(outline) + '\n')
specs.close()
table.close()
# end
 