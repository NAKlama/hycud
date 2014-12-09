#!/usr/bin/python3

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

import numpy as np
import re
import scipy.optimize as optimize

from multiprocessing import Pool
from sdfFunction import sdfFunction

from StatItem import npStatItem
from Minimizer import Minimize

import copy
# import sys

class boundingFunc:
  def __init__(self, bounds):
    self.bounds = bounds

  def __call__(self, par):
    bAvg    = (self.bounds[:,0] + self.bounds[:,1]) / 2
    bHDiff  = (self.bounds[:,0] - self.bounds[:,1]) / 2
    res     = np.sin(par) * bHDiff + bAvg
    return res


class sdfMinFunc:
  def __init__(self, models, wH, inVals, inWeights, res, threads, bounds):
    self.models     = models
    self.wH         = wH
    self.inVals     = inVals
    self.inWeights  = inWeights
    self.residue    = res
    self.threads    = threads
    self.bounds     = bounds

  def __call__(self, par):
    # Defining bounds
    # boundF  = boundingFunc(self.bounds)
    # usePar  = boundF(par)
    # print("Testing: ", usePar, " (", par, ")", sep='', end='')
    # print("Testing: ", par, end='')

    # Preparing function to Minimize
    task    = sdfFunction(self.wH)
    task.changeParameters(par)
    # task.changeParameters(usePar)
    task.singleResidueSet(self.residue)

    # Minimizing function
    try:
      sdf_pool = Pool(processes=self.threads)
      result   = sdf_pool.map(task, self.models)
      sdf_pool.close()
    except ProcessLookupError:
      pass

    resSDF =  [ { 'res':    self.residue
                , 'wh':     wh
                , 'J':      None
                , 'R':      None
                , 'sigma':  None
                , 'noe':    None
                , 'ratio':  None
                , 'tau':    None
                , 'fragC':  0
                } for wh in self.wH ]


    for i in range(len(result)):
      m = self.models[i]
      for f in m.fragments:
        if not f.partial and self.residue in f.residues:
          f.sdf = []
          for r in result[i]:
            if r['frag'] == f.num:
              f.sdf.append(r)

    # Averaging over fragments within each model
    for m in self.models:
      m.resSDF =  [ { 'res':    self.residue
                    , 'wh':     wh
                    , 'J':      None
                    , 'R':      None
                    , 'sigma':  None
                    , 'noe':    None
                    , 'ratio':  None
                    , 'tau':    None
                    , 'fragC':  0
                    } for wh in self.wH ]
      for f in m.fragments:
        if not f.partial and self.residue in f.residues:
          for sdf in f.sdf:
            r = sdf['res']
            for sdfAvg in m.resSDF:
              if sdfAvg['res'] == r and sdfAvg['wh'] == sdf['wh']:
                if sdfAvg['J'] is None:
                  sdfAvg['J']     = npStatItem()
                  sdfAvg['R']     = npStatItem()
                  sdfAvg['sigma'] = npStatItem()
                  sdfAvg['noe']   = npStatItem()
                  sdfAvg['ratio'] = npStatItem()
                  sdfAvg['tau']   = npStatItem()
                sdfAvg['R'].addValue(sdf['R'])
                sdfAvg['noe'].addValue(np.array([sdf['noe']]))

    # Avereraging over all models
    for sdf in resSDF:
      if sdf['J'] is None:
        sdf['J']      = npStatItem()
        sdf['R']      = npStatItem()
        sdf['sigma']  = npStatItem()
        sdf['noe']    = npStatItem()
        sdf['ratio']  = npStatItem()
        sdf['tau']    = npStatItem()
      for m in self.models:
        for msdf in m.resSDF:
          if msdf['J'] is not None and sdf['res'] == msdf['res'] and sdf['wh'] == msdf['wh']:
            sdf['R'].addStatItem(msdf['R'])
            sdf['noe'].addStatItem(msdf['noe'])

    comp = []
    for wh in self.wH:
      comp.append(  list(
                    [ sdf['R'].getAvg()[0]
                    , sdf['R'].getAvg()[1]
                    , sdf['noe'].getAvg()[0]
                    ] for sdf in resSDF if sdf['wh'] == wh )
                 )
    inVals    = None
    inWeights = None
    for v in self.inVals:
      if v['res'] == self.residue:
        # inVals    = np.array(v['data'])[:,:,0:2]
        # inWeights = np.array(v['weight'])[:,:,0:2]
        inVals    = np.array(v['data'])
        inWeights = np.array(v['weight'])

    diff = (np.array(comp) - inVals) / inWeights

    return np.ravel(diff)


class sdfMinimize:
  def __init__(self, opt, inFile):
    weights       = False
    LarmorF       = []
    inData        = []
    inWeights     = []
    guessTau      = []
    guessS        = []
    fitRes        = []

    optionsParser   = re.compile("^\\s*(\\w+)\\((.*)\\)\\s*")
    frequencyParser = re.compile("^\\s*([+-]?\\d+(?:\\.\\d+)?(?:[eE][+-]?\\d+))\\s*")
    lines = []
    with open(inFile, mode="r", encoding="utf-8") as inF:
      lines = inF.readlines()
    firstLine = lines[0]
    options   = firstLine.split(';')
    for o in options:
      if o.strip() == "weights" or o.strip() == "uncertainties":
        weights       = True
      else:
        resOpt = optionsParser.match(o)
        if resOpt:
          key   = resOpt.group(1)
          param = resOpt.group(2)
          if key == "freq" or key == "larmor":
            freqList = param.split(',')
            for f in freqList:
              resFreq = frequencyParser.match(f)
              if resFreq:
                LarmorF.append(float(resFreq.group(1)))
          if key == "Tau" or key == "tau":
            tList = param.split(',')
            for t in tList:
              guessTau.append(float(t))
          if key == "S" or key == "s":
            sList = param.split(',')
            for s in sList:
              guessS.append(float(s))
    if LarmorF == []:
      LarmorF.append(opt.larmorFreq)

    fieldMult = 3
    dataParserDef = "^\\s*(\\d+)"
    for f in LarmorF:
      dataParserDef += "\\s+(\\S+)"
      dataParserDef += "\\s+(\\S+)"
      dataParserDef += "\\s+(\\S+)"
      if weights:
        fieldMult = 6
        dataParserDef += "\\s+(\\S+)"
        dataParserDef += "\\s+(\\S+)"
        dataParserDef += "\\s+(\\S+)"
    dataParser = re.compile(dataParserDef)

    # for f in LarmorF:
    #   inData.append([])
    #   inWeights.append([])
    inData = []

    for l in lines:
      res = dataParser.match(l)
      if res:
        data    = res.groups()
        resid     = int(data[0])
        index   = 1
        Data    = []
        Weight  = []
        for i in range(len(LarmorF)):
          d = []
          w = []
          if weights:
            d.append( [ float(data[index])
                      , float(data[index+2])
                      , float(data[index+4])
                      ] )
            w.append( [ float(data[index+1])
                      , float(data[index+3])
                      , float(data[index+5])
                      ])
          else:
            d.append( [ float(data[index])
                      , float(data[index+1])
                      , float(data[index+2])
                      ] )
            w.append([1.0,1.0,1.0])
          Data.append(d)
          Weight.append(w)
          index += fieldMult
        inData.append(  { 'res':    resid
                        , 'data':   copy.copy(Data)
                        , 'weight': copy.copy(Weight)
                        } )
    for l in lines:
      res = dataParser.match(l)
      if res:
        data   = res.groups()
        resid  = int(data[0])
        bounds = np.array([ [1e-12,   1e-8]
                          , [1e-5,    0.99999]
                          ])

        minFunc = sdfMinFunc(
          opt.models.models, LarmorF, inData, inWeights, resid, opt.threads,
          bounds)

        result = None
        sB  = []
        tB  = []
        S   = None
        Tau = None
        if len(guessS) >= 2:
          sB  = [ guessS[0] - guessS[1]
                , guessS[0] + guessS[1]
                ]
          if len(guessS) >= 3:
            sB.append(int(guessS[2])-1)
          else:
            sB.append(2)
          S = np.arange(sB[0], sB[1]+ (sB[1] - sB[0])/sB[2], (sB[1] - sB[0])/sB[2])
        else:
          S = np.array([guessS[0]])
        if len(guessTau) >= 2:
          tB  = [ guessTau[0] - guessTau[1]
                , guessTau[0] + guessTau[1]
                ]
          if len(guessS) >= 3:
            tB.append(int(guessTau[2])-1)
          else:
            tB.append(2)
          Tau = np.arange(tB[0], tB[1]+ (tB[1] - tB[0])/tB[2], (tB[1] - tB[0])/tB[2])
        else:
          Tau = np.array(guessTau[0])
        inV = (Tau[0:sB[2]+1], S[0:tB[2]+1])
        inVals = np.array(inV)

          # halfInCount = int(inVals.shape[1]**inVals.shape[0] / 2)
          # print(halfInCount)
        minimize    = Minimize(remFact=10.0, useCount=25)
        result      = minimize(
          minFunc,
          inVals,
          np.array([2e-10, 0.2]),
          bounds
          )
        # result = optimize.leastsq(
        #   minFunc,
        #   np.array([guessTau[0],guessS[0]]),
        #   # epsfcn=1e-5
        #   )
        # result  = optimize.minimize(
        #   minFunc,
        #   np.array([1e-5, 0.5]),
        #   # method="L-BFGS-B",
        #   method="TNC",
        #   # method="SLSQP",
        #   bounds=[(1e-14, 1e-5), (0.1, 0.9)],
        #   options={'disp': True}
        #   )
        # result  = optimize.basinhopping(
        #   minFunc,
        #   np.array([1e-5, 0.5]),
        #   # method="L-BFGS-B",
        #   # method="TNC",
        #   # method="SLSQP",
        #   # bounds=[(1e-14, 1e-5), (0.1, 0.9)],
        #   # options={'disp': True}
        #   )

        bFunc = boundingFunc(bounds)
        print("{}:".format(resid), result)


        # res = bFunc(result[0])
        # print("{}: {:e}  {:e}".format(resid, res[0], res[1]))


