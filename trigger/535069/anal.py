#!/usr/bin/env python3
import sys
dir = ""
def openfile(filename, max = 0):
  filename = dir+filename
  print("openfile Opening:",filename)
  try:
    with open(filename) as f:
      lines = f.readlines()
  except IOError:
    print("Can not open file, exiting") 
    exit(1)
  trginfo = []  
  print("Opened ok")  
  id = 0
  ic = 0
  #for line in lines:
  MAX = len(lines)
  if max: MAX = max
  for id in range(MAX):
    line = lines[ic]
    ic += 1
    if len(line) < 2: continue  
    #print(line)
    l = line.replace(':',' ')
    items = l.split(" ")
    #print(items[1])
    trg = []
    trg.append(int(items[1]))
    trg.append(int(items[2]))
    trg.append(int(items[3],base=16))
    trg.append(int(items[4],base=16))
    trg.append(id)
    trg.append(0) # not used
    trginfo.append(trg)
    id += 1
  return trginfo  
def anal(us,sk):
  #print(un)
  hist = {}
  print("#us:",len(us)," #sk:",len(sk))
  for i in us:
    print("===>",i)
    n = 0
    for j in sk:
      if i[0] == j[0]:
        print(j, compare(i,j))
        j[5] = 1
        n+= 1
    if n in hist:
      hist[n] += 1
    else:
      hist[n] = 1
  print("freqA:",hist)
  notused = 0
  for i in sk:
    if i[5] ==0 : notused +=1
  print("notused:", notused)
def checkSorted(tr):
  i0 =  tr[0]
  n = 0
  for i in tr:
    if i[0] < i0[0]:
      #print("NOT SOOOORTED")
      n += 1
    i0 = i
  print(len(tr), "sorted violations:",n)    
def compare(tr1,tr2):
  mask = ""
  if tr1[0] == tr2[0] : mask = "bcs "
  if tr1[1] == tr2[1] : mask += "col "
  if tr1[2] == tr2[2] : mask += "trg "
  if tr1[3] == tr2[3] : mask += "sel "
  return mask
def correlateAll(us,sk):
  mask = 1 << 54
  delta = 10
  corr = {}
  for i in us:
    posi = i[0]
    trigi = i[3] & mask
    if trigi:
      for j in sk:
        trigj = j[2] & mask
        if trigi and trigj:
          posj = j[0]
          if posj + delta < posi: continue
          if posi + delta < posj: continue
          dist = posi - posj + delta
          if dist in corr:
            corr[dist] += 1
          else:
            corr[dist] = 1
          print("idtrg ",i[4],j[4],dist)
  print("corr",corr)
def freqBC(tr):
  bcf = {}
  prev = 0
  count = 1 
  for i in tr:
    if i[0] != prev:
      if prev : bcf[prev] =  count
      count = 1
      prev = i[0]
    else: count += 1
  bcf[prev] = count
  #print("bcFreq:",bcf)
  freq = {}
  for key,value in bcf.items():
    if value in freq:
      freq[value] += 1
    else:
      freq[value] = 1
  print("freq:",freq)
if __name__ == "__main__":
  us = openfile("sel_unskimmed.txt")
  sk = openfile("trg_skimmed.txt")
  print("unskimmed")
  checkSorted(us)
  freqBC(us)
  print("skimmed:")
  checkSorted(sk)
  freqBC(sk)
  anal(us,sk)
  #correlateAll(us,sk)
  anal(sk,us)
