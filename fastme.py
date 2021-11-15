from subprocess import Popen, PIPE
from time import strftime, gmtime
from sys import stderr, exit
from os import access, F_OK, rename, remove, mkdir

names = open("goodlist.txt", "r")

if not access("Trees", F_OK):
  mkdir("Trees")
if not access("Trees/FastME", F_OK):
  mkdir("Trees/FastME")
if not access("Trees/FastME/logs", F_OK):
  mkdir("Trees/FastME/logs")
if not access("Distances", F_OK):
  mkdir("Distances")  
for line in names:
  x = line.strip()
  print x, strftime("%d%b%Y %H:%M:%S", gmtime())
  command = ["seqret", "Alignments/" + x + ".afa", "phylip::infile", "-auto"]
  seqret = Popen(command, stdout = PIPE, stderr = PIPE)
  (out, err) = seqret.communicate()
  if err:
    stderr.write("{}: seqret error\n{}\n".format(x, err))
    exit(1)
  command = ["fastme", "-p", "-i", "infile"]
  command += ["-o", "Trees/FastME/" + x + ".tre", "-O", "Distances/" + x + ".dist", "-I",\
              "Trees/FastME/logs/" + x + ".info", "-f", "8"]
  i = 0
  while not access("Trees/FastME/" + x + ".tre", F_OK):
    if i > 0:
      stderr.write(x + ": " + str(i) + '\n')
    fastme = Popen(command, stdout = PIPE, stderr = PIPE)
    (out, err) = fastme.communicate() 
    if err: stderr.write("{}: fastme error\n{}\n".format(x, err))
    i += 1
  logf = open("Trees/FastME/logs/" + x + ".log", "w")
  logf.write(out)
  logf.close()
  remove("infile")

names.close()