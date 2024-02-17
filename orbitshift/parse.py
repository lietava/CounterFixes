#!/usr/bin/env python3
import sys
dir = ""
def openfile(filename):
  filename = dir+filename
  print("openfile Opening:",filename)
  try:
    with open(filename) as f:
      lines = f.readlines()
  except IOError:
    print("Can not open file, exiting") 
    exit(1)
  print("Opened ok")
  return lines
def parseExcel(periods):
  ptable = {}  
  filltable = {}
  afill = ''
  for line in periods:
    items = line.split("\t")
    if len(items) != 3 :
      print("items != 3:",items)
      continue
    period = items[0]
    fill = items[1]
    run = items[2].strip()
    if period != '':
      print("New period:",period)
      filltable = {}
      if fill != '':
        afill = fill  
        filltable[fill] = [run]  
        ptable[period] = filltable
      else:
        print("Error: new period without fill")
        exit(1)
    elif fill != '':
      afill = fill  
      filltable[afill] = [run]
    else:
      filltable[afill].append(run)
  #printptable)    
  return ptable    
def parseFirstOrbit(firstorbots):
  ptable = {}  
  for line in firstorbits:
    items = line.split()
    if len(items) != 3:
      continue
    #print(items)
    key = items[0]
    run = items[1]
    if key in ptable:
      ptable[key].add(run)
    else:
      ptable[key] = {run}
  #print(ptable)
  return ptable
def isRunAvailable(runs,fills):
  for fill in fills:
    #print(fills[fill])  
    n = 0  
    for run in fills[fill]:
      #print(run, run in runs)  
      if run in runs: 
        n += 1
    if n == 0:
      print("===> no runs:",fill)
    else:
      print(n,fill)    
if __name__ == "__main__":
  lhc22o = openfile("LHC22o.txt")
  ptex = parseExcel(lhc22o)
  firstorbits = openfile("orbit.log")
  ptfo = parseFirstOrbit(firstorbits)
  print(ptex["LHC22o"])
  print("ptab fo")
  print(ptfo["LHC22o"])
  isRunAvailable(ptfo["LHC22o"],ptex["LHC22o"])

