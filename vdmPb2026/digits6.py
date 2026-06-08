#!/usr/bin/env python3
import sys,re
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
  orbit0 = 0;
  i = 0;
  ntot = 0
  map = {}
  ntrgmask = 0
  ninpmask = 0
  for line in lines:
    #print(line)
    line = " ".join(line.split())
    if len(line)==0: continue
    if line[0] == "#": continue
    if line.find("Processing") >= 0: continue
    #if line.find("TCR") >= 0: continuea
    flagClass = 0
    if line.find("TCR") >= 0: 
     flagClass = 1
    if line.find("int") >= 0: continue
    if line.find("stop") >= 0: return map
    #items = line.split(" ")
    items = re.split(r"\s*[: ]\s*", line)
    #print(line)
    if len(items) < 4:
     print(line)
     print(len(items))
     exit(1)
    #print("=== ",flagClass, items)
    ipos = 0
    if flagClass: ipos = 1
    orbit = int(items[ipos],16)
    bcid = int(items[ipos + 1],16)
    inpmask = int(items[ipos + 2],16)
    trgmask = int(items[ipos + 3],16)
    gbcid = orbit*3564 + bcid
    if inpmask: ninpmask += 1
    if trgmask: ntrgmask += 1
    map[gbcid] = [bcid,orbit,inpmask,trgmask]
    #print(hex(orbit))
    if orbit < orbit0:
        print("decreasing",orbit,orbit0)
        i += 1
    orbit0 = orbit
    ntot += 1
  #for l in out: print(l)
  print("Tot number:",ntot,"inp:",ninpmask," trg:",ntrgmask)
  #print("map:",len(map))
  #for k in map: print("parser:",hex(k))
  return map
########################################
NCLASS = 64
class Descriptor:
  def __init__(self):  
    self.name = None
    self.mask = None
    self.inps =[]
class CTPConfig:
  def __init__(self):
    self.Inputs = {}
    self.Descs = {}
    self.Clusters = {}
    self.TClasses = {}
    self.ClassCounters = []
    for i in range(NCLASS):
      self.ClassCounters.append(0)
  def parseConfig(self,run):
    filename = dir+run+".rcfg2"
    print("openfile Opening:",filename)
    with open(filename) as f:
      lines = f.readlines()
    print("file opened ok:",len(lines))  
    state = None 
    for rline in lines:
      line=rline.strip(' \n')  
      if len(line) == 0: continue
      #print(state)  
      if line.find("INPUTS") >= 0 :
        state = "inputs"
        continue
      if line.find("BCMASKS") >= 0:
        state = "bcmasks"
        continue
      if line.find("DESCRIPTORS") >= 0:
        state = "descs"
        continue
      if line.find("LTG") >= 0:
        state = "ltg"
      if line.find("cluster") >= 0:
        state = "cluster"  
        continue
      items = line.split(" ")
      if state == "inputs":
        #print(items)
        self.Inputs[items[1].strip()] = int(items[2])
      if state == "descs":
        desc =  Descriptor()
        desc.name = items[0]
        mask = 0
        for d in items[1:]:
          dd = d.strip()
          desc.inps.append(dd)
          if dd in self.Inputs:
            mask += 1 << (self.Inputs[dd]-1)
        desc.mask = mask  
        #print(desc.inps)
        self.Descs[desc.name] = desc 
      if state == "cluster":
        if line[0] != 'C':
          print("Parse ERROR",line)
          exit(1)
        dname = items[2].strip()
        #print(dname)
        self.TClasses[1 << int(items[1])] = self.Descs[dname].mask 
    return 0
  def checkClassVsInput(self,inpmask,classmask,bc,orb):
    nerrs = 0
    nfired = 0
    icls=0
    for cls in self.TClasses:
      inpsok = True  
      isclassfired = classmask & cls
      #print("clsmask,cls:",hex(classmask),hex(cls), isclassfired)
      if isclassfired:
        nfired += 1  
        #print("isfired:",isclassfired,classmask,cls,nfired)
        self.ClassCounters[icls] += 1
        inpsok = inpmask & self.TClasses[cls]
      if not inpsok:
        #if icls > 8 and icls < 21:
        print("class:",hex(cls),hex(self.TClasses[cls]),"in classmask",hex(classmask),"not ok with inps",hex(inpmask)," BC,orb:",hex(bc),hex(orb)) 
        nerrs += 1
      icls += 1
    return nerrs, nfired
  def printConfig(self):
    for inp in self.Inputs:
      print(inp,self.Inputs[inp])
    for d in self.Descs:
      print(d,self.Descs[d].inps,self.Descs[d].mask)  
    for c in self.TClasses:
      print(hex(c),hex(self.TClasses[c]))  
  def printClassCounters(self):
    print("ClassCounters:")
    #print(self.ClassCounters)
    for i in range(len(self.ClassCounters)):
      print(i,self.ClassCounters[i])
###########################      
def glob(bc,orbit):
    return orbit*3564 + bc
def bcorb(glob):
    bc = glob % 3564
    orb = int(glob / 3564)
    return bc,orb
#######################################
lmmask = 0xfff;
l0mask = 0xfff000
l1mask = 0xffffff000000
nlost = 0
def shift(tforbit,key,orb,bc,inp,level,outmap):
  global nlost
  shift = 15
  lxmask = l0mask  
  if level == 1: 
    lxmask = l1mask
    shift = 295+1
  #print("shift:",shift,level)
  #print(outmap)
  if (orb > tforbit) or (bc >= shift):
    bcn,orbn = bcorb(key-shift)
    globn = glob(bcn,orbn)
    if key - globn != shift: print("error 1")
    if globn in outmap:  
      inpn = outmap[globn][2]
      if inpn & lxmask: 
        print("WARNING 1 ",i)
      inpn = inpn | (inp&lxmask)
      outmap[globn][2] = inpn
      #print("case 1:",outmap[globn])
    else:    
      outmap[globn] = [bcn,orbn,inp&lxmask,0]
      #print("case 2:",outmap[globn])
  else:
    #print("LOST",orb,bc)
    nlost += 1
  return 0
###############################################  
nl0 = 0
nl1 = 0
nl01 = 0
def shiftDigits(digits):
  global nl0, nl1, nl01
  tforbit = 0  
  outmap = {}  
  n1e=0
  i = 0
  for key in digits:
    bc = digits[key][0]
    orb = digits[key][1]
    inp = digits[key][2]
    trg = digits[key][3]
    #
    if (orb & 0x1f) == 0:
      tforbit = orb
      #print("new tf:", hex(tforbit))
    #
    lm = (inp & lmmask) > 0
    l0 = (inp & l0mask) > 0
    l1 = (inp & l1mask) > 0
    lut = lm + (l0<<1) + (l1<<2)
    #print("lut:",lut,digits[key])
    if lut==0:
      #print("case No inputs ",lut)
      outmap[key] = [bc,orb,inp,trg]
    elif lut == 1: #LM
      outmap[key] = [bc,orb,inp,trg]
    elif lut == 2: #L0
      shift(tforbit,key,orb,bc,inp,0,outmap)
      if trg: 
        outmap[key] = [bc,orb,0,trg]
        print("=====> L0 error adding trg wo inp:",outmap[key]," ori:",bc,orb)
        nl0 += 1
    elif lut == 4: #L1
      shift(tforbit,key,orb,bc,inp,1,outmap)
      if trg: 
        outmap[key] = [bc,orb,0,trg]
        print("=====> L1 error adding trg wo inp:",outmap[key]," ori:",bc,orb)
        nl1 += 1
    elif lut == 6: # L0 and L1  
      shift(tforbit,key,orb,bc,inp,0,outmap)
      shift(tforbit,key,orb,bc,inp,1,outmap)
      if trg:
        outmap[key] = [bc,orb,0,trg]
        print("=====> L0 and  L1 error adding trg wo inp:",outmap[key]," ori:",bc,orb) 
        nl01 == 1
    elif lut == 3:   # LM and L0
      #print("case LM and L0 ",lut)
      #if trg : print("3 trg",key,orb,bc,trg)
      shift(tforbit,key,orb,bc,inp,0,outmap) #!!
      outmap[key] = [bc,orb, inp&(~l0mask),trg]
    elif lut == 5: # LM and L1
      shift(tforbit,key,orb,bc,inp,1,outmap)
      outmap[key] = [bc,orb, inp&(~l1mask),trg]
    elif lut == 7: # LM, L0 and L1
      shift(tforbit,key,orb,bc,inp,0,outmap)
      shift(tforbit,key,orb,bc,inp,1,outmap)
      outmap[key] = [bc,orb, inp&(lmmask),trg]
    else: print("FATAL internal ERROR")
    #print(i,outmap)  
    i+=1
  return outmap
def shiftGlobal(shift,digits):
  outmap = {}  
  for key in digits:
    bc = digits[key][0]
    orb = digits[key][1]
    inp = digits[key][2]
    trg = digits[key][3]
    if trg:
      bcn,orbn = bcorb(key-shift)
      globn = glob(bcn,orbn)
      if globn in outmap:
        trgn = outmap[globn][3]
        if trgn:
          print("====>error: trgn already there")
        else:  
          outmap[globn][3] = trg
      else:
        outmap[globn] = [bcn,orbn,0,trg] 
      if inp:
        globo = glob(bc,orb)
        outmap[globo] = [bc,orb,inp,0]
    else:
      outmap[key]=[bc,orb,inp,0]   
  print("shiftGlobal done")    
  return outmap    
def printRec(d):
  print(hex(d[1]),hex(d[0]),hex(d[2]),hex(d[3])," (",d[1],d[0],")")
  return
def printMap(m):
  for i in m:
    d = m[i]  
    print(i,hex(d[1]),hex(d[0]),hex(d[2]),hex(d[3])," (",d[1],d[0],")")
    #print(d[1],d[0],d[2],(d[3])," (",d[1],d[0],")")
  return          
def compare(dspy,ds):
  n = 0  
  nmiss = 0
  print(" checking # of k in dspy and not in ds")
  for k in dspy:
    if k in ds:
      if dspy[k] != ds[k] : nmiss += 1
      continue
    else:
      d=dspy[k]  
      print(k,hex(d[1]),hex(d[0]),hex(d[2]),hex(d[3])," (",d[1],d[0],")")
      n+=1
  print("# of k in dspy and not in ds",n, " nmiss:", nmiss)
  print("checking # of k in ds and not in dspy")
  n = 0
  nmiss = 0
  for k in ds:
    if k in dspy:
      if dspy[k] != ds[k] : nmiss += 1
      continue
    else:
      n+=1
      d=ds[k]
      print(k,hex(d[1]),hex(d[0]),hex(d[2]),hex(d[3])," (",d[1],d[0],")") 
  print("# of k in ds and not in dspy",n," nmiss:",nmiss)
def checkClassesVsInputs(ds,cfg):
  nerrs = 0  
  nfired = 0
  ncls = 0
  for key in ds:
     inpmask = ds[key][2]
     classmask = ds[key][3]
     bc = ds[key][0]
     orb = ds[key][1]
     if classmask: ncls+=1
     #prt = 0
     #if classmask == 1024: prt=1
     #if inpmask & 0x302000: prt=1
     #if classmask: prt=1
     #if prt:
       #print("prt",key,hex(inpmask),hex(classmask),ds[key][0],ds[key][1])
       #def checkClassVsInput(self,inpmask,classmask,bc,orb):
     n1,n2 = cfg.checkClassVsInput(inpmask,classmask,bc,orb)
     nerrs += n1
     nfired += n2
  print("Classes versus Inputs errors:",nerrs, "ncls:",ncls, "nfired:",nfired," tot digits:",len(ds))
def CountLxTC(ds):
  nLM = 0
  nL0 = 0
  nL1 = 0
  nTCwI = 0
  nTCwoI = 0
  for d in ds:
    inp = digits[d][2]
    trg = digits[d][3]
    if inp & 0xfff: nLM += 1
    if inp & 0xfff000: nL0 += 1
    if inp & 0xffffff000000: nL1 += 1
    if trg:
      if inp:nTCwI += 1
      else: nTCwoI += 1
  print("LM:",nLM," L0:",nL0," L1:",nL1," T with I:",nTCwI," T wo I:",nTCwoI)
if __name__ == "__main__":
 dWD = openfile("digitsWShift.txt")
 dWOD = openfile("digitsNoShift.txt")
 dWODS = shiftDigits(dWOD)
 #compare(dWD,dWODS)
 #print("NLOST:",nlost," Ntot:",len(dWODS),len(dWOD))
 #print("Warns: L0:",nl0," L1:",nl1," L0L1:",nl01)
