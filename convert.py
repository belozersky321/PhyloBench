#############################################################################
#    This script is written on Python 2.7
#    It is a part of PhyloBench. 
#    It reads the file "temp.tre", which is an output of TNT
#    and makes the file "newick.tre" with the same tree in Newick format.
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

'''Convert output of TNT to Newick format'''

instream = file("temp.tre", "r")
for line in instream.readlines():
  if line[0] == '(':
    outstream = file("newick.tre", "w")
    tmpline = ''
    flag = 0
    for x in line:
      if flag == 0:
        if x != ' ':
          f = 0
          outstream.write(x)
        else:
          f = 1
      else:
        if x == ')':
          f = 1
          outstream.write(')')
        else:
          f = 0
          if x == ';' or x == '\n':
            outstream.write(x)
          else:
            outstream.write(',' + x)
      flag = f
    # for
    outstream.close()
  # if
# for
instream.close()
#end    
          
