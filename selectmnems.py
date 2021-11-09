#############################################################################
#    This script is written on Python 2.7
#    It is a part of PhyloBench. 
#    It requires files "names.dmp" from NCBI Taxonomy Database, "pfamseq" from Pfam
#    and "table.txt" made by the script "maketable.py".
#    It creates the files "list.txt" and "pfamlist.txt".
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

fs = frozenset

def create_graph(pfs):
  result = dict()
  result["vertices"] = set(pfs.keys())
  result["power"] = dict()
  for _u in result["vertices"]:
    for _v in result["vertices"]:
      _edge = fs((_u, _v))
      if _edge not in result["power"]:
        result["power"][_edge] = len(pfs[_u] & pfs[_v])
  return result
# create_graph

def remove_least(graph):
  _valence = dict()
  for _v in graph["vertices"]:
    _valence[_v] = 0
    for _u in graph["vertices"]:
      if _v != _u:
         _edge = fs((_u, _v))
         _valence[_v] += graph["power"][_edge]
  minval = -1
  for _v in graph["vertices"]:
    if minval < 0:
      minval = _valence[_v]
      least = _v
    if _valence[_v] < minval:
      minval = _valence[_v]
      least = _v

  result = dict()
  result["vertices"] = graph["vertices"] - {least}
  result["power"] = dict()
  for _u in result["vertices"]:
    for _v in result["vertices"]:
      _edge = fs((_u, _v))
      if _edge not in result["power"]:
        result["power"][_edge] = graph["power"][_edge] 
  return result
# remove_least

def vertices(graph):
  return graph["vertices"]

def select_max(d, tax): 
  '''
  d is dictionary with keys = organism mnemonics and values = sets of Pfam AC's presented in the organism
  tax is dictionary with keys = organism mnemonics and values = taxon id's (of some level) of this organism
  The result is a subdictionary of d, its keys are organism mnemonics having maximum number of Pfam AC's among the same taxon
  '''
  result = dict()
  mn = dict() # keys = taxon id's of some level
  for x in d:
    if tax[x] in mn:
      mn[tax[x]] |= {x}
    else:
      mn[tax[x]] = {x}
  for y in mn:
    m = 0
    for x in mn[y]:
      if len(d[x]) > m:
        m = len(d[x])
        repr = x
    if m:
      result[repr] = d[repr]
  return result

def ne(d):
  '''
  Return the subdictionary of d with non-empty values
  '''
  result = dict()
  for x in d:
    if d[x]:
      result[x] = d[x]
  return result

def utax(d, tax):
  '''
  d is dictionary with keys = organism mnemonics and values = sets of Pfam AC's presented in the organism
  tax is dictionary with keys = organism mnemonics and values = taxon id's (of some level) of this organism
  The result is a dictionary with keys = taxons and values = union of values of d for all members of the taxon
  '''
  result = dict() # keys = taxon id's of some level
  for x in d:
    if tax[x] in result:
      result[tax[x]] |= d[x]
    else:
      result[tax[x]] = set(d[x])
  return result

def ispf(x):
  '''Does the row of "swisspfam" contain a Pfam accesion number?'''
  if len(x) > 8:
    if len(x.split('.')) == 2 and len(x.split('.')[0]) == 7 and x[:2] == "PF":
      return True
    else:
      return False
  else:
    return False

####################################
# main                             #
####################################
'''
Input: files 
	"names.dmp" (from NCBI)
	"table.txt" (made with maketable.py)
	"swisspfam" (from Pfam)
        All files are presumed to be in the parent folder.
Parameters:
	taxon, from which organisms should be
	n, the disired number of organisms
Optional parameters:
	level, the taxonomy level, organisms should be not more than one from each taxon of this level
	exc, the blacklist of subtaxons of the taxon
	pfam, a file with Pfam AC's that should be used (by default, use all Pfam AC's)
Output: files
	"list.txt", organism mnemonics, one in a row, a subset of organisms having maximum common Pfam families
	"pfamlist.txt", Pfam AC's, one in a row, that are presented in all selected organisms
'''
from sys import argv, stderr, exit

taxname = None
n = None
extaxes = []
level = None
inpfam = set()
for arg in argv[1:]:
  if '=' in arg:
    (name, value) = tuple(arg.split('='))
    if name.startswith("tax"):
      taxname = value
    if name == 'n':
      n = int(value)
    if name.startswith("exc"):
      extaxes = value.split(',')
    if name == 'level':
      level = value
    if name.startswith("pfam"):
      inpf = open(value, "r")
      for x in inpf:
        pf = x.strip()
        if pf:
          inpfam |= {pf}
      inpf.close()
if (not taxname) or (not n):
  stderr.write("Usage:\npython2.7 " + argv[0] + " tax=<taxon> n=<number> [ext=<comma-separated list>] [level=<level>] [pfam=<file>]\n")
  stderr.write("Example:\npython2.7 " + argv[0] + " tax=Eukaryota n=40 exc=Metazoa,Fungi level=genus\n")
  exit(1)

# Determining the taxon id (taxid) and id's of blacklist taxons (extaxids)
names = open("../names.dmp", "r")
taxid = ""
extaxids = set()
if level:
  taxdict = dict()
  for x in names:
    ell = x.split('\t|\t')
    if 'scientific name' in ell[3]:
      if ell[1] == taxname:
        taxid = ell[0]
      for ext in extaxes:
        if ell[1] == ext:
          extaxids |= {ell[0]}
      taxdict[ell[0]] = ell[1]
else: # we do not need all taxon id's, only the taxon's and from the blacklist (extaxes)
  x = names.readline()
  while x and ((not taxid) or (len(extaxids) < len(extaxes))):
    ell = x.split('\t|\t')
    if 'scientific name' in ell[3]:
      if  ell[1] == taxname:
        taxid = ell[0]
      for ext in extaxes:
        if ell[1] == ext:
          extaxids |= {ell[0]}
    x = names.readline()
names.close()

if not taxid:
  stderr.write("Cannot find scientific name \"" + taxname + "\"!\n")
  exit(1)
if (len(extaxids) < len(extaxes)):
  stderr.write("Not all scientific names found: " + ", ".join(extaxes) + "\n")

print "Name OK,", taxname, "=", taxid

table = open("../table.txt", "r")
pfs = dict()
ti = 0
if level:
  header = table.readline().strip()
  hh = header.split('\t')
  for i in range(len(hh)):
    if hh[i] == level:
      ti = i
  if ti:
    tax = dict()
  else:
    stderr.write("Unknown level: " + level + ", ignored\n")
for x in table:
  ell = x.strip().split('\t')
  be = False
  if taxid in ell and not (extaxids & set(ell)):
    be = True
  if be:
    pfs[ell[0]] = set()
    if ti:
      tax[ell[0]] = ell[ti]
# for
table.close()
if not pfs:
  stderr.write("Cannot find taxon in table (unsupported level?)\n")
  exit(1)
else:
  print "Speclist OK,", len(pfs.keys()), "species"

indata = open("../swisspfam", "r")
mnem = ""
for x in indata:
  ell = x.strip().split()
  if ell:
    if ell[0][0] == '>':
      mnem = ell[0][1:].split('_')[1]
    else:
      if mnem in pfs:
        for z in ell:
          if ispf(z): pfs[mnem].add(z)
indata.close()

if inpfam: # a list of Pfam families is set
  for mnem in pfs:
    pfs[mnem] &= inpfam

pfs = ne(pfs)
print len(pfs.keys()), "mnemonics in swisspfam"
speclist = sorted(pfs.keys())
if not speclist:
  stderr.write("Cannot find species in swisspfam\n")
  exit(1)

if ti:
  taxpfs = utax(pfs, tax)
else:
  taxpfs = pfs

print len(taxpfs.keys()), "different taxons of level", level

pfset = set()
for taxon in taxpfs:
  pfset |= taxpfs[taxon]

print len(pfset), "total Pfam families"

if len(taxpfs.keys()) <= n:
  for taxon in taxpfs: pfset &= taxpfs[taxon]
  for mnem in pfs: pfs[mnem] &= pfset
  if ti: 
    pfs = select_max(pfs, tax)
    for mnem in pfs: pfset &= pfs[mnem]
  pflist = sorted(list(pfset))
  speclist = sorted(pfs.keys())
else:
  # first variant
  graph = create_graph(taxpfs)
  taxlist1 = list(vertices(graph))
  while len(taxlist1) > n:
    graph = remove_least(graph)
    taxlist1 = list(vertices(graph))

  pfset1 = set(pfset)
  for taxon in taxlist1: pfset1 &= taxpfs[taxon]
  print "Variant 1:", len(pfset1), "Pfam families (preliminary)"
  pfs1 = dict()
  if ti: 
    for mnem in pfs:
      if tax[mnem] in taxlist1: pfs1[mnem] = pfs[mnem] & pfset1
    pfs1 = select_max(pfs1, tax)
  else:
    for mnem in pfs:
      if mnem in taxlist1: pfs1[mnem] = pfs[mnem] & pfset1
  for mnem in pfs1: pfset1 &= pfs1[mnem]
  print "Variant 1:", len(pfset1), "Pfam families"

  popfs = dict()
  for taxon in taxpfs:
    for pf in taxpfs[taxon]:
      if pf in popfs: popfs[pf] += 1
      else: popfs[pf] = 1

  # second variant
  goodpfset = set()
  for pf in popfs:
    if popfs[pf] >= n/2:
      goodpfset.add(pf)
  print len(goodpfset), "more or less represented Pfam families"
  for taxon in taxpfs: taxpfs[taxon] &= goodpfset
  taxpfs = ne(taxpfs)

  graph = create_graph(taxpfs)
  taxlist2 = list(vertices(graph))
  while len(taxlist2) > n:
    graph = remove_least(graph)
    taxlist2 = list(vertices(graph))

  pfset2 = set(pfset)
  for taxon in taxlist2: pfset2 &= taxpfs[taxon]
  print "Variant 2:", len(pfset2), "Pfam families (preliminary)"
  pfs2 = dict()
  if ti: 
    for mnem in pfs:
      if tax[mnem] in taxlist2: pfs2[mnem] = pfs[mnem] & pfset2
    pfs2 = select_max(pfs2, tax)
  else:
    for mnem in pfs:
      if mnem in taxlist2: pfs2[mnem] = pfs[mnem] & pfset2
  for mnem in pfs2: pfset2 &= pfs2[mnem]
  print "Variant 2:", len(pfset2), "Pfam families"

  # third variant
  goodpfset = set()
  for pf in popfs:
    if popfs[pf] >= n:
      goodpfset.add(pf)
  print len(goodpfset), "well represented Pfam families"
  for taxon in taxpfs: taxpfs[taxon] &= goodpfset
  taxpfs = ne(taxpfs)

  graph = create_graph(taxpfs)
  taxlist3 = list(vertices(graph))
  while len(taxlist3) > n:
    graph = remove_least(graph)
    taxlist3 = list(vertices(graph))

  pfset3 = set(pfset)
  for taxon in taxlist3: pfset3 &= taxpfs[taxon]
  print "Variant 3:", len(pfset3), "Pfam families (preliminary)"
  pfs3 = dict()
  if ti: 
    for mnem in pfs:
      if tax[mnem] in taxlist3: pfs3[mnem] = pfs[mnem] & pfset3
    pfs3 = select_max(pfs3, tax)
  else:
    for mnem in pfs:
      if mnem in taxlist3: pfs3[mnem] = pfs[mnem] & pfset3
  for mnem in pfs3: pfset3 &= pfs3[mnem]
  print "Variant 3:", len(pfset3), "Pfam families"

  if len(pfset1) >= len(pfset2) and len(pfset1) >= len(pfset3):
    pflist = sorted(list(pfset1))
    speclist = sorted(pfs1.keys())
  if len(pfset2) >= len(pfset1) and len(pfset2) >= len(pfset3):
    pflist = sorted(list(pfset2))
    speclist = sorted(pfs2.keys())
  if len(pfset3) >= len(pfset2) and len(pfset3) >= len(pfset1):
    pflist = sorted(list(pfset3))
    speclist = sorted(pfs3.keys())

print "Choice OK,", len(speclist), "species,", len(pflist), "Pfam families"

outspecs = open("list.txt", "w")
for species in speclist:
  outspecs.write(species + '\n')
outspecs.close()

outpfam = open("pfamlist.txt", "w")
for pf in pflist:
  outpfam.write(pf + '\n')
outpfam.close()

if ti:
  outtax = open(level + ".txt", "w")
  taxdict['NA'] = 'NA'
  for species in speclist:
    outtax.write(species + '\t' + taxdict[tax[species]] + '\n')
  outtax.close() 

# end
