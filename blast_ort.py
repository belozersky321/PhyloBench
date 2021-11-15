#############################################################################
#    This script is written on Python 2.7
#    It is a part of PhyloBench. 
#    It requires files "list.txt" and "pfamlist.txt" created by "selectmnems.py"
#    and directory "Domains" created by "extractseq.py".
#    The command line argument is (an arbitrary) short prefix for names of created fasta files.
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

from sys import argv, exit
import os
from random import seed, shuffle
from subprocess import Popen, PIPE

pfamlistname = "pfamlist.txt"  # Pfam ID's
mnemlistname = "list.txt"  # Mnemonics of species

if not os.access(mnemlistname, os.F_OK):
  print "The file {} must be in the current folder".format(mnemlistname)
  exit(1)
if not os.access(pfamlistname, os.F_OK):
  print "The file {} must be in the current folder".format(pfamlistname)
  exit(1)
try:
  setid = argv[1] # AC Actinobacteria, AG for Archaea-genus etc.
except Exception:
  print "Short name is not specified, set to \"XX\""
  setid = "XX"

specfile = open(mnemlistname, "r")
speclist = specfile.readlines()
specfile.close()

numbers = dict()

seed()
N = 0   # succesful ortgroups
nbh = 0 # no blast hits
nog = 0 # no ortologous groups
ok = dict()
if not os.access("Sets", os.F_OK):
  os.mkdir("Sets")
if not os.access("Sequences", os.F_OK):
  os.mkdir("Sequences")
pfamfile = file(pfamlistname, "r")
for pfamline in pfamfile:
  pfamid = pfamline.strip()  # ID of Pfam family 
  if len(pfamid) > 0:
    print pfamid
    for s in speclist:
      ss = s.strip() # species mnemonics
      command = ["grep", "-c", ">", "Domains/" + pfamid + "/" + ss + ".fasta"]
      grepc = Popen(command, stdout = PIPE, stderr = PIPE)
      (out, err) = grepc.communicate()
      numbers[ss] = int(out) # How many representatives of the species in the Pfam family
    # for 
    
    # find the minimum number
    minspec = speclist[0].strip()
    minnum = numbers[minspec]
    for species in numbers:
      if numbers[species] < minnum:
        minnum = numbers[species]
        minspec = species
    # for
    if minnum == 0:
      print "No sequences for %s in %s" % (minspec, pfamid)
    else:
      print "%d sequences for %s" % (minnum, minspec)
      query = "Domains/" + pfamid + "/" + minspec + ".fasta"  # sequences of minspec
      db = "Domains/" + pfamid + "/tmpblast"
      dblist = []
      for s in speclist:
        ss = s.strip()
        if ss != minspec:
          dblist += ["Domains/" + pfamid + "/" + ss + ".fasta"]
      command = ["makeblastdb", "-dbtype", "prot", "-in", " ".join(dblist), "-out", db]
      makeblastdb = Popen(command, stdout = PIPE, stderr = PIPE) # blast db from sequences of all species except minspec
      (out, err) = makeblastdb.communicate()
      if err: print err
      command = ["blastp", "-query", query, "-db", db]
      command += ["-word_size", "2"]
      command += ["-window_size", "0"]
      command += ["-threshold", "7"]
      command += ["-comp_based_stats", "0"]
      command += ["-num_alignments", "1000"]
      command += ["-outfmt", "6 qseqid sseqid score"]
      blastp = Popen(command, stdout = PIPE, stderr = PIPE)
      (blastres, err) = blastp.communicate()
      if err:
        print err
      else:
        series = dict() # keys are sequences of minspec, values are OGs as sets of sequences
        scores = dict() # keys are sequences of minspec, values are dictionaries with keys = sequences of other species and values = scores
        for line in blastres.split('\n'):
          if line:
            ell = line.split()
            q = ell[0] # query
            s = ell[1] # subject
            score = int(ell[2])  # blast score
            if q in series:
              if s in scores[q]: scores[q][s] = max(score, scores[q][s]) # s-variant; m-variant is scores[q][s] += score
              else: scores[q][s] = score
            else:
              series[q] = set([q])
              scores[q] = dict()
              scores[q][s] = score
        for q in series: 
          ok[q] = True 
          maxscore = dict()
          maxseq = dict()
          for species in numbers:
            if species != minspec:
              maxscore[species] = 0
              maxseq[species] = ""
          # for
          for s in scores[q]: # sequences of other species
            cursp = s.split("/")[0].split("_")[1] # species mnemonics of s
            if scores[q][s] > maxscore[cursp]:
              maxscore[cursp] = scores[q][s]
              maxseq[cursp] = s
          for species in maxscore:
            if maxscore[species] > 0:
              series[q].add(maxseq[species])
            else:
              print "No hit for %s in %s!" % (q, species)
              ok[q] = False
          if not ok[q]: nbh += 1

      # Check series for intersecting 
      num = 0
      serlist = series.keys() # sequences of the minimal species 
      sernum = len(serlist)
      if sernum > 1:
        for i in range(sernum):
          q1 = serlist[i]
          if ok[q1]:
            for j in range(i + 1, sernum):
              q2 = serlist[j]
              if ok[q2]:
                intersection = series[q1] & series[q2]
                for s in intersection:
                  if scores[q1][s] < scores[q2][s]:
                    ok[q1] = False
                  if scores[q1][s] > scores[q2][s]:
                    ok[q2] = False
                if intersection and ok[q1]:
                    ok[q2] = False
          if ok[q1]: 
            num += 1
      else:
        if sernum == 1: num += 1
      if num == 0:
        nog += 1
      if num > 0:
        # renumeration and output
        if not os.access("Sets/" + pfamid, os.F_OK):
          os.mkdir("Sets/" + pfamid)
        k = 0
        for q in series:
          if ok[q]:
            k += 1
            print ("k = %d; q = %s; orgroup: %d sequences" % (k, q, len(series[q])))
            pfamfasta = "Sequences/" + setid + "_" + pfamid.split(".")[0] + "_" + str(k) + ".fasta"
            outf = open(pfamfasta, "w")
            currorlist = list(series[q])
            shuffle(currorlist)
            for seq in currorlist:
              org = seq.split("/")[0].split("_")[1]
              fasta = "Domains/" + pfamid + "/" + org + ".fasta"
              usa = fasta + ":" + seq.replace('/', '_')
              os.system("seqret " + usa + " -filter >> Sets/" + pfamid + "/set" + str(k) + ".fasta")
              command = ["seqret", usa, "-osformat2", "plain", "-filter"]
              seqret = Popen(command, stdout = PIPE, stderr = PIPE)
              (plainseq, err) = seqret.communicate()
              if err: print "seqret", q, seq, err
              else:
                command = ["descseq", "-name", org, "-sformat1", "plain", "-filter"]
                descseq = Popen(command, stdin = PIPE, stdout = PIPE, stderr = PIPE)
                (fastaseq, err) = descseq.communicate(plainseq)
                if err: print "descseq", q, seq, err
                else:
                  outf.write(fastaseq)
            # for
            outf.close()
          # if good
        # for
        N += k
      # if num > 0:
      print pfamid, "%d orgroups, total %d\n" % (num, N)
    # if minnum > 0
  # if pfamid is not empty
# for pfamid
pfamfile.close()

print nbh, "cases of BLAST failure"
print nog, "families with no orthologous groups"
# end