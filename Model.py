# HYCUD
# Copyright (C) 2014 Klama, Frederik and Rezaei-Ghaleh, Nasrollah
#
# This file is part of HYCUD.
#
# HYCUD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HYCUD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HYCUD.  If not, see <http://www.gnu.org/licenses/>.


import os
import io
import re
# import string
import gc
import copy
# import math

from os import path
from multiprocessing    import Process, Queue
from Fragment           import Fragment, FragmentStatistics
from HelperFunctions    import printWarning, removeChainIDs, outSumm
from HelperFunctions    import ANSI_ESC, flush
from StatItem           import StatItem
from Parsers            import modelParser, endmdlParser
from Parsers            import atomParser, pdbResParser, pdbParser
from Parsers            import diffMatrixParser
import numpy as np

# Linear algebra functions from scipy are using BLAS and are faster
# so we try them before using numpy as a fallback
try:
  import scipy.linalg as linalg
except ImportError:
  import numpy.linalg as linalg


class readThread:
  def __init__(self, fileName, bufferSize=(1024**2)*16):
    self.file       = fileName
    self.bufferSize = bufferSize

  def __call__(self, q):
    fd  = open(self.file, mode='r', encoding='utf-8')
    buf = fd.readlines(self.bufferSize)
    while(len(buf) > 0):
      q.put(buf)
      buf = fd.readlines(self.bufferSize)
    q.put([])
    fd.close()


class writeThread:
  def __init__(self, fileName):
    self.file = fileName

  def __call__(self, q):
    fd  = open(self.file, "w", encoding='utf-8')
    inBuf = q.get()
    buf = []
    while inBuf is not None:
      buf.append(inBuf)
      if len(buf) >= 1000:
        fd.writelines(buf)
        buf = []
      inBuf = q.get()
    fd.writelines(buf)
    fd.close()


class Model:
  """Class for the different pdt models"""
  def __init__(self, opt, num):
    self.num = num
    self.opt = opt
    m = num
    n = 1
    modelPath = ""
    while m > 9:
      m -= m % 10
      m /= 10
      n *= 10
      t = "{:09n}".format(m * n)
      modelPath = path.join(t, modelPath)
    self.path        = path.join(opt.tmpPath, modelPath)
    self.basename    = path.join(self.path, "{:09n}".format(num))
    if not path.exists(self.basename):
      os.makedirs(self.basename)
    # self.pdb        = self.basename + '.pdb'
    self.fragments   = []
    self.script      = ""
    self.transScript = ""
    # self.result     = self.basename + '.result'

  def addFrag(self, num, partial=False):
    frag = Fragment(num, self.basename, partial)
    self.fragments.append(frag)
    return frag

  def getPDB(self):
    return self.basename + '.pdb'

  def getResult(self):
    return self.basename + '.result'

  def getTransRes(self):
    return self.basename + '.trans.result'

  def getSprouted(self):
    if self.opt.skipREMO:
      return self.basename + '.pdb'
    else:
      return self.basename + '.sprouted'

  def doneParsing(self):
    for f in self.fragments:
      f.doneParsing()

  def getWeight(self):
    weight = 0.0
    for f in self.fragments:
      weight += f.getWeight()
    return weight


class printFragStatVal:
  """Helper class to print fragment statistics"""
  def __init__(self, subject, avg, unit, stdDev):
    self.s = subject
    self.a = avg
    self.u = unit
    self.d = stdDev

  def __call__(self):
    print("  {}: {:e} {} (+/- {:e})".format(self.s, self.a, self.u, self.d))


class Models:
  def __init__(self):
    self.models = []
    self.statsExist = False
    self.fragStats  = None
    self.transStats = None
    self.maxRes     = 0
    self.residues   = []

  def parsePDT(self, opt, fileName):
    """Function parsing the PDT file, spltting it int PDBs inside the temporary directory
       and populating protein information"""
    # Initializing the readThread class with the file to read
    pTarget = readThread(fileName)

    # Generating a Queue and passing it to the read thread, it will read up to
    # 32 times the buffer size, and output a list of lines into q
    q = Queue(32)
    p = Process(target=pTarget, args=(q,))
    p.start()

    modelNo     = -1
    gotRes      = 0
    output      = None
    # r           = 1

    # Initialize variable line with first line
    if opt.verbose > 0:
      print("Parsing PDT input file and splitting into seperate models")
    lines = q.get()
    # As long as we do not encounter an empty set of lines (indicates EOF)
    while len(lines) > 0:
      for line in lines:
        line  = removeChainIDs(line)

        # Testing for 'MODEL' line
        res   = modelParser.match(line)
        if res:
          # r = 1
          modelNo = int(res.group(1))
          model   = Model(opt, modelNo)
          self.models.append(model)

          # Here we are initiailizing the output thread, with a queue of 1024 lines
          oQ = Queue(1024)
          oT = writeThread(model.getPDB())
          oP = Process(target=oT, args=(oQ,))
          oP.start()
          output = "Thread"

        # Testing for 'ENDMDL' line
        res = endmdlParser.match(line)
        if res:
          if output is not None:
            # Making sure everything has been written into the file
            oQ.put(None)
            oP.join()
            output = None
            modelNo = -1

        # As long we are inside a model
        if modelNo > 0:
          pdbMatch = pdbResParser.match(line)
          if pdbMatch:
            # prefix    = pdbMatch.group(1)
            resNum    = int(pdbMatch.group(2))
            # postfix   = string.rstrip(pdbMatch.group(3))

            if gotRes < resNum:  # i.e. first atom of new Residue
              gotRes = resNum
              opt.resCount += 1

            if resNum > self.maxRes:
              self.maxRes = resNum
            if resNum not in self.residues:
              self.residues.append(resNum)

            # Write the pdbLine into the output file
            oQ.put(line)
          else:
            # Catching invalid PDB atom lines and warning about them
            # this is only a problem if all atoms of a residue can not
            # be parsed by the previous parser
            atomMatch = atomParser.match(line)
            if atomMatch:
              oQ.put(line)
              printWarning("Could not parse the following line:")
              print(line)

      # Read next line to process
      lines = q.get()

    p.join()
    gc.collect()

  def splitFragments(self, opt):
    """This function splits each protein model into fragments as defined by *opt.fragmentation*"""
    # firstFrag = opt.fragmentation.firstFragment
    # fragments = opt.fragmentation.fragments
    fragDef = opt.fragmentation
    self.fragmentation = fragDef
    REMOpdt = ["REMARK REMO processed version of {}\n".format(opt.pdt)]

    for m in self.models:
      if opt.verbose > 0:
        print("{0}2K{0}1GNow splitting fragments for model {1:6n}/{2:n}".format(
          ANSI_ESC(), m.num, len(self.models)), end='')
        flush()

      script  = []
      script.append("#!/bin/bash\n")
      script.append("cd {}\n".format(m.basename))
      script.append("cat > tmp.pl << EOF\n")
      script.append("#!/usr/bin/perl\n")
      script.append("\\$fragR = 0; \\$fragV = 0;\n")
      script.append("while(<>){\n")
      script.append(" if(/^\\s+This\\sfile\\s*:\\s+(\\S+)-res.txt/) {\n")
      script.append("  \\$frag = \\$1; }\n")
      script.append(" if(/^\\s+Radius\\sof\\sgyration:\\s+(\\S+)/) {\n")
      script.append("  print \"{:09n}.\\$frag  Stokes:  \\$1\\n\"; }}\n".format(m.num))
      script.append(" if(/^\\s+Intrinsic\\sviscosity:\\s+(\\S+)\\s+cm/) {\n")
      script.append("  print \"{:09n}.\\$frag  Viscos:  \\$1\\n\"; }}\n".format(m.num))
      script.append(" if(/^\\s+Harm.*time:\\s+(\\S+)/) {\n")
      script.append("  print \"{:09n}.\\$frag  HarmMe:  \\$1\\n\"; }}\n".format(m.num))
      script.append(" if(/^\\s*(-?\\d\\.\\d+[eE][-+]\\d+\\s+.*)\\s*$/) {\n")
      script.append("  print \"{:09n}.\\$frag Mat:  \\$1 \\$2 \\$3 \\$4 \\$5 \\$6\n\"; }}\n".format(m.num))
      script.append("}\nEOF\n")
      script.append("chmod +x tmp.pl\n")

      REMOpdt.append("MODEL {:n}\n".format(m.num))

      with io.open(m.getSprouted(), mode='r', encoding='utf-8') as pdbInFile:
    #     currRes       = 1
    #     currFrag      = 0
    #     fragSum       = firstFrag + fragments[currFrag].num
    #     static        = fragments[currFrag].static
        pdbLines      = pdbInFile.readlines()
        for line in pdbLines:
          REMOpdt.append(line)
        for frag in range(0, fragDef.fragmentCount()):
          fragment      = m.addFrag(frag, fragDef.fragments[frag].partial)
          static        = fragDef.fragments[frag].static
          pdbOut        = []
          for line in pdbLines:
            pdbMatch = pdbParser.match(line)
            if pdbMatch:
              # atomNum   = int(pdbMatch.group(1))
              atomName  = pdbMatch.group(2).strip()
              resName   = pdbMatch.group(3).strip()
              resNum    = int(pdbMatch.group(4))
              atomX     = float(pdbMatch.group(5))
              atomY     = float(pdbMatch.group(6))
              atomZ     = float(pdbMatch.group(7))

              if fragDef.partOfFragment(frag, resNum):
                # if atomName.strip() == "N":
                #   print(line)
                #   print(atomX, atomY, atomZ, np.array([atomX, atomY, atomZ], np.float64))
                fragment.addAtom(atomName, resNum, resName, np.array([atomX, atomY, atomZ], np.float64))
                pdbOut.append(line)

          opt.HydroConf[1]          = "{}\n".format(fragment.basename)
          opt.HydroConf[2]          = "{}.pdb\n".format(fragment.basename)
          opt.HydroConf[10]         = "{:f}\n".format(fragment.getWeight())

          with io.open(fragment.getDat(), mode='w', encoding='utf-8') as confOutFile:
            confOutFile.writelines(opt.HydroConf)

          with io.open(fragment.getPDB(), mode='w', encoding='utf-8') as pdbOutFile:
            pdbOutFile.writelines(pdbOut)

          if static:
            script.append("echo \"{:09n}.{}  Stokes:  0.0\" >> {}\n".format(m.num, fragment.basename, m.getResult()))
            script.append("echo \"{:09n}.{}  Viscos:  {:e}\"  >> {}\n".format(m.num, fragment.basename, fragDef.fragments[frag].viscos, m.getResult()))
            script.append("echo \"{:09n}.{}  HarmMe:  {:e}\"  >> {}\n".format(m.num, fragment.basename, fragDef.fragments[frag].HarmMe, m.getResult()))
          else:
            script.append("ln -fs {} hydropro.dat\n".format(fragment.getDat()))
            script.append("{} > /dev/null 2>&1\n".format(opt.exePath))
            script.append("./tmp.pl < {}/{}-res.txt >> {}\n".format(m.basename, fragment.basename, m.getResult()))

      REMOpdt.append("ENDMDL\n")

      if not opt.keepTemp:
        script.append("rm -Rf tmp.pl\n")
      if opt.verbose > 3:
        print("Finishing off script...", end='')
      script.append("cd ..\n")
      if not opt.keepTemp:
        script.append("rm -Rf {}\n".format(m.basename))
        script.append("rm -Rf {}.pdb\n".format(m.basename))

      scriptOutFile = "{}.script".format(m.basename)
      with io.open(scriptOutFile, mode='w') as scriptOut:
        scriptOut.writelines(script)

      os.system("chmod +x {}".format(scriptOutFile))
      m.script = scriptOutFile
      m.doneParsing()
      gc.collect()
      if opt.verbose > 3:
        print("done")

    if opt.verbose > 0:
      print("")
    if opt.REMOout != "":
      with open(opt.REMOout, "w") as REMOoutFile:
        REMOoutFile.writelines(REMOpdt)


  def translationalDiff(self, opt):
    for m in self.models:
      script = []
      script.append("#!/bin/bash\n")
      script.append("cd {}\n".format(m.basename))
      script.append("cat > tmpTrans.pl << EOF\n")
      script.append("#!/usr/bin/perl\n")
      script.append("while(<>) {\n")
      script.append("  if(/^\\s+Translational\\s+diffusion\\s+coefficient:\\s+(\\S+)/) {")
      script.append("    print \"{:09n}  Trans:  \\$1\\n\"; }}\n".format(m.num))
      script.append("}\nEOF\n")
      script.append("chmod +x tmpTrans.pl\n")

      (sproutedPath, sproutedFile) = path.split(m.getSprouted())

      opt.HydroConf[1]          = "FullModel\n"
      opt.HydroConf[2]          = "FullModel.pdb\n"
      opt.HydroConf[10]         = "{:f}\n".format(m.getWeight())

      datFile    = m.basename + "/hydropro.trans.dat"
      scriptFile = "{}.trans.script".format(m.basename)
      with io.open(datFile, mode='w', encoding='utf-8') as confOutFile:
        confOutFile.writelines(opt.HydroConf)

      script.append("ln -fs hydropro.trans.dat hydropro.dat\n")
      script.append("ln -fs ../{} FullModel.pdb\n".format(sproutedFile))
      script.append("{} > /dev/null 2>&1\n".format(opt.exePath))
      script.append("./tmpTrans.pl < {}/FullModel-res.txt >> {}\n".format(m.basename, m.getTransRes()))
      script.append("rm -Rf tmpTrans.pl\n")
      script.append("cd ..\n")
      with io.open(scriptFile, mode='w', encoding='utf-8') as scriptOut:
        scriptOut.writelines(script)

      os.system("chmod +x {}".format(scriptFile))
      m.transScript = scriptFile

  def parseResults(self, opt):
    """This function parses the values calculated by HydroPro and stores them in the fragments"""
    stokesParser = re.compile("^(\\d+)\\.Frag(\\d+)\\s+Stokes:\\s+(\\S+)\\s*$")
    harmMeParser = re.compile("^(\\d+)\\.Frag(\\d+)\\s+HarmMe:\\s+(\\S+)\\s*$")
    viscosParser = re.compile("^(\\d+)\\.Frag(\\d+)\\s+Viscos:\\s+(\\S+)\\s*$")
    matrixParser = re.compile("^(\\d+)\\.Frag(\\d+)\\s+Mat:\\s+(.*)\\s*$")
    translParser = re.compile("^(\\d+)\\s+Trans:\\s+(\\S+)\\s*$")

    for m in self.models:
      resFile = m.getResult()
      transResFile = m.getTransRes()
      if path.exists(transResFile):
        with io.open(transResFile) as transResFD:
          lines = transResFD.readlines()
          for line in lines:
            if opt.verbose > 3:
              print("Parsing line: \"{}\"".format(line.strip()))
            parsed = False

            res = translParser.match(line)
            if res:
              parsed = True
              m.transDiff = float(res.group(2))

      if not opt.onlyTrans:
        with io.open(resFile) as resFD:
          lines = resFD.readlines()
          for line in lines:
            if opt.verbose > 3:
              print("Parsing line: \"{}\"".format(line.strip()))
            parsed = False

            res = stokesParser.match(line)
            if res:
              # num   = int(res.group(1))
              frag  = int(res.group(2))
              r     = float(res.group(3))
              for fragment in m.fragments:
                if fragment.num == frag:
                  parsed = True
                  fragment.calcValues(radius=r)

            res = harmMeParser.match(line)
            if res:
              # num   = int(res.group(1))
              frag  = int(res.group(2))
              hm    = float(res.group(3))
              for fragment in m.fragments:
                if fragment.num == frag:
                  parsed = True
                  fragment.calcValues(HarmMe=hm)

            res = viscosParser.match(line)
            if res:
              # num   = int(res.group(1))
              frag  = int(res.group(2))
              eta   = float(res.group(3))
              for fragment in m.fragments:
                if fragment.num == frag:
                  parsed = True
                  fragment.calcValues(viscosity=eta)

            res = matrixParser.match(line)
            if res:
              frag    = int(res.group(2))
              matStr  = res.group(3)
              matRes  = diffMatrixParser.match(matStr)
              if matRes:
                for fragment in m.fragments:
                  if fragment.num == frag:
                    parsed = True
                    fragment.diffMat.append([ matRes.group(1)
                                            , matRes.group(2)
                                            , matRes.group(3)
                                            , matRes.group(4)
                                            , matRes.group(5)
                                            , matRes.group(6)
                                            ])

            if not parsed:
              printWarning("Could not parse line: \"{}\"".format(line.strip()))

  def calculateWeightingFactors(self, opt):
    """This function calculates the weighting due to fragment distances"""
    for m in self.models:
      for fragI in m.fragments:
        if not fragI.partial:
          i         = fragI.num
          commonSum = 0

          # NA = 6.02214129e23    # $N_A$
          #
          #       1 molecule         1
          # c   = ----------   = ------------ (mol/Å³)
          #  ij    (   _    )3    N  × 6 r³
          #       (  ³√6 r   )     A      ij
          #        (      ij)
          #
          # 1Å³ = 1e-30 m³ = 1e-27ℓ = 1e-24 cm³
          #
          #                 1
          # c   = ---------------------
          #  ij              -24        (mol / cm³)
          #        6 × N × 10    × r³
          #             A           ij
          #
          #
          #                 1                1
          # constJ  = -------------- = -------------
          #                     -24     3.613284774
          #           6 × N × 10
          #                A
          #

          constJ  = 1 / 3.613284774
          centerI = fragI.center.getCenter()
          consideredFragments = []
          if opt.sdf:
            for frag in m.fragments:
              if not frag.partial and frag.num % opt.fragSize == i % opt.fragSize:
                consideredFragments.append(frag)
            if opt.sdfPartial:
              # Consiering the beginning
              normBegFrag = min(fragI.residues)
              while normBegFrag > 0:
                normBegFrag -= opt.fragSize
              normBegFrag += opt.fragSize
              if opt.verbose > 4:
                print("F =", fragI.num, end='')
                print("   B:", normBegFrag, end='')
              if normBegFrag >= opt.fragSize / 2:
                for f in m.fragments:
                  if f.partial and max(f.residues) == normBegFrag - 1:
                    consideredFragments.append(f)
                    if opt.verbose > 4:
                      print(" add", f.residues, end='')

              # Considering the end
              normEndFrag = max(fragI.residues)
              while normEndFrag - 1 < self.maxRes:
                normEndFrag += opt.fragSize
              normEndFrag -= opt.fragSize
              if opt.verbose > 4:
                print("   E:", normEndFrag, "  ", end='')
              if self.maxRes - normEndFrag >= opt.fragSize / 2:
                for f in m.fragments:
                  if f.partial and min(f.residues) == normEndFrag + 1:
                    consideredFragments.append(f)
                    if opt.verbose > 4:
                      print("add", f.residues, end='')
              if opt.verbose > 4:
                print("")
          else:
            for frag in m.fragments:
              if not frag.partial:
                consideredFragments.append(frag)

          # print("{:d}: fCount:{:d}".format(i, len(consideredFragments)), end="")
          fragI.adjFragCount = len(consideredFragments)
          for fragJ in consideredFragments:
            j = fragJ.num
            if i != j:
              centerJ = fragJ.center.getCenter()
              rij = linalg.norm(centerI - centerJ)
              # print("\n   i: ", centerI, " j: ", centerJ,
              #       " rij:", rij,
              #       " constJ:", constJ,
              #       " wj:", fragJ.getWeight(),
              #       " eta:", fragJ.getEta(),
              #       " cs:", ((constJ * fragJ.getWeight())) / (rij ** 3) * fragJ.getEta()
              #       ,end="")
              commonSum += ((constJ * fragJ.getWeight())) / (rij ** 3) * fragJ.getEta()
          commonSum += 1
          # print(" \n   cs+1:", commonSum)

          fragI.values.corrected   = fragI.values.values.mult(commonSum)
          fragI.weightingFactor    = commonSum

  # def calculateWeightingFactorsMT(self, opt):
  #   from MTFunctions import calcWeightFact
  #   from multiprocessing import Pool
  #   pool = Pool(processes=opt.threads)
  #   pool.map(calcWeightFact, self.models)

  def initFragStats(self, o):
    self.fragStats = FragmentStatistics(o)

  def populateFragStats(self):
    self.fragStats.populate(self)

  def calcStats(self):
    self.statsExist = True

    self.Rsum  = StatItem()  # count
    self.HMsum = StatItem()  # count
    self.eta   = StatItem()  # count
    self.RcW   = StatItem()  # weight
    self.RuW   = StatItem()  # weight
    self.RcP   = StatItem()  # protons
    self.RuP   = StatItem()  # protons
    self.HMcW  = StatItem()  # weight
    self.HMuW  = StatItem()  # weight
    self.HMcP  = StatItem()  # protons
    self.HMuP  = StatItem()  # protons
    for m in self.models:
      for f in m.fragments:
        w   = f.getWeight()
        p   = f.getProtons()
        r   = f.values.getR(corr=False)
        rC  = f.values.getR(corr=True)
        hm  = f.values.getHM(corr=False)
        hmC = f.values.getHM(corr=True)
        eta = f.values.getEta()
        self.Rsum.addValue( rC)
        self.HMsum.addValue(hmC)
        self.eta.addValue(  eta)
        self.RcW.addValue(  rC , w)
        self.RuW.addValue(  r  , w)
        self.RcP.addValue(  rC , p)
        self.RuP.addValue(  r  , p)
        self.HMcW.addValue( hmC, w)
        self.HMuW.addValue( hm , w)
        self.HMcP.addValue( hmC, p)
        self.HMuP.addValue( hm , p)

  def output(self, opt):
    if not self.statsExist:
      self.calcStats()
    print("\n--- S U M M A R Y ---")
    if opt.verbData :
      outSumm("Average corrected  ","radius             ","                               ",self.Rsum.getAvg() ,"cm   ",self.Rsum.getStdDev())
    outSumm(  "Average corrected  ","harmonic mean time ","                               ",self.HMsum.getAvg(),"s    ",self.HMsum.getStdDev())
    if opt.verbData:
      outSumm("Average            ","intrinsic viscosity","                               ",self.eta.getAvg()  ,"cm3/g",self.eta.getStdDev())
      if opt.showWeightedAvg:
        outSumm("Average corrected  ","radius             ","(Norm: Mol Weight)             ",self.RcW.getAvg()  ,"cm   ",self.RcW.getStdDev())
        outSumm("Average uncorrected","radius             ","(Norm: Mol Weight)             ",self.RuW.getAvg()  ,"cm   ",self.RuW.getStdDev())
        outSumm("Average corrected  ","radius             ","(Norm: non-exchangable protons)",self.RcP.getAvg()  ,"cm   ",self.RcP.getStdDev())
        outSumm("Average uncorrected","radius             ","(Norm: non-exchangable protons)",self.RuP.getAvg()  ,"cm   ",self.RuP.getStdDev())
    if opt.showWeightedAvg:
      outSumm(  "Average corrected  ","harmonic mean time ","(Norm: Mol Weight)             ",self.HMcW.getAvg() ,"s    ",self.HMcW.getStdDev())
      outSumm(  "Average uncorrected","harmonic mean time ","(Norm: Mol Weight)             ",self.HMuW.getAvg() ,"s    ",self.HMuW.getStdDev())
      outSumm(  "Average corrected  ","harmonic mean time ","(Norm: non-exchangable protons)",self.HMcP.getAvg() ,"s    ",self.HMcP.getStdDev())
      outSumm(  "Average uncorrected","harmonic mean time ","(Norm: non-exchangable protons)",self.HMuP.getAvg() ,"s    ",self.HMuP.getStdDev())

  def transStatisticss(self):
    """This function does the translational statistics"""
    stat = StatItem()
    for m in self.models:
      stat.addValue(m.transDiff)

    self.transStats = {
      'avg': stat.getAvg(),
      'harmMe': stat.getHarmMe(),
      'stdDev': stat.getStdDev()}

  def outputFragResults(self, opt):
    """Print per fragment statistics"""
    i = 1
    for stat in self.fragStats.stats:
      print("--- Fragment #{:n}  {:n} residues ({}) ---".format(
        i, stat.residues, stat.resPrint))
      printList = []
      avg = stat.avg
      if opt.harmonicMean:
        avg = stat.harmMean
      if opt.verbData:
        printList.append(printFragStatVal("Radius                   ",
                          avg.getR(corr=False),
                          "cm   ", stat.stdDev.getR(corr=False)))
      printList.append(printFragStatVal("Harmonic mean time       ",
                        avg.getHM(corr=False),
                        "s    ", stat.stdDev.getHM(corr=False)))
      printList.append(printFragStatVal("Intrinsic viscosity      ",
                        avg.getEta(),
                        "cm3/g", stat.stdDev.getEta()))
      if(opt.verbData):
        printList.append(printFragStatVal("Corr. Radius             ",
                          avg.getR(corr=True),
                          "cm   ", stat.stdDev.getR(corr=True)))
      printList.append(printFragStatVal("Corr. Harmonic mean time ",
                        avg.getHM(corr=True),
                        "s    ", stat.stdDev.getHM(corr=True)))
      for item in printList:
        item()
      i += 1

  def outputTransFragResults(self, opt):
    avg = 'avg'
    if opt.harmonicMean:
      avg = 'harmMe'
    out = printFragStatVal(
      "Translational diffusion coefficient",
      self.transStats[avg],
      "cm2/s",
      self.transStats['stdDev'])
    out()

  def size(self):
    return len(self.models)

  def subModels(self, size, offset=0):
    endPos = size+offset
    self.models = self.models[offset:endPos]

  def removeModels(self, mNumList):
    out = Models()
    out.statsExist = False
    out.fragmentation = self.fragmentation
    modelNumbers = []
    for m in self.models:
      modelNumbers.append(m.num)
    for exclNum in mNumList:
      if exclNum not in modelNumbers:
        printWarning("Models number {:n} does not exist, ignoring it.".format(exclNum))
    for m in self.models:
      if m.num not in mNumList:
        out.models.append(m)
    return out

  def getStdDevList(self, o):
    out = []
    if self.fragStats is None:
      self.fragStats = FragmentStatistics(o)
      self.fragStats.populate(self)
    for stat in self.fragStats.stats:
      out.append(stat.stdDev.getHM(corr=True))
    return out

  def findOutliers(self):
    out = []
    avgList = []
    for stat in self.fragStats.stats:
      avgList.append(stat.avg.getHM(corr=True))
    for i in range(0, self.fragmentation.fragmentCount()):
      out.append({'num':None, 'val':None})

    for m in self.models:
      for f in range(0, len(m.fragments)):
        dev = m.fragments[f].getR(corr=True) - avgList[f]
        dev *= dev
        if out[f]['val'] is None or out[f]['val'] < dev:
          out[f]['val'] = dev
          out[f]['num'] = m.num

    return out

