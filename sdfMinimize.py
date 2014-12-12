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
# import scipy.optimize as optimize
import math

from multiprocessing import Pool
from sdfFunction import sdfFunction

from StatItem import npStatItem
from Minimizer import Minimize
from PascalTriangle import *

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
  def __init__(self, models, wH, inVals, inWeights, res, threads, bounds=None):
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


    # Testing for boundary conditions
    if self.bounds is not None:
      outFail = []
      for i in range(len(self.wH)):
        outFail += [1e100, 1e100, 1e100]
      if par[0] < self.bounds[0][0] or par[0] > self.bounds[0][1]:
        return outFail
      if par[1] < self.bounds[1][0] or par[1] > self.bounds[1][1]:
        return outFail


    # Preparing function to Minimize
    task    = sdfFunction(self.wH)
    task.changeParameters(par)
    # task.changeParameters(usePar)
    task.singleResidueSet(self.residue)

    # print(par)
    # print(self.residue)

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

    # print(result)

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
    # print(resSDF)
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
    models        = opt.models.models
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
    residues = []

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
        residues.append(resid)

    residues  = sorted(residues)

    # Split center residues and edges
    centerRes = []
    edgeRes   = []
    results   = []
    for r in residues:
      if r > opt.fragSize and r <= residues[-1] - opt.fragSize:
        centerRes.append(r)
      else:
        edgeRes.append(r)

    # Boundry conditions
    bounds = [[1e-12, 1e-8], [1e-5, 0.99999]]

    centerPoint = 0
    for r in centerRes:
      centerPoint += r
    centerPoint = math.floor(float(centerPoint) / len(centerRes))

    # Initialize Minimization Functions
    minimize = Minimize(gridSpec=(guessTau, guessS))

    # Get starting point for grid search
    initVal = np.array([guessTau[0], guessS[0]])
    optRss  = 1e200
    if minimize.enableGrid:
      print("Grid search for Residue ", centerPoint)
      F = sdfMinFunc(models, LarmorF, inData, inWeights, centerPoint, opt.threads, bounds)
      R = minimize(F, None, mode="grid")

      # Spot Clusters
      clusterCutOff = 0.001 ############
      n   = 0
      while n < len(R):
        val  = R[n]['par']
        Rout = R[:n+1]
        for r in R[n+1:]:
          ratio = np.fabs((r['par'] - val) / val)
          if np.any(np.greater(ratio, clusterCutOff)):
            Rout.append(r)
        R = Rout
        n += 1

      results.append({'res': centerPoint, 'min': R})
      initVal = R[0]['par']
      optRss  = R[0]['rss']

    oddRes  = centerRes[::2]
    evenRes = centerRes[1::2]

    # First minimization on optimal value
    for r in evenRes:
      if r != centerPoint:
        print("Minimizing residue", r)
        F = sdfMinFunc(models, LarmorF, inData, inWeights, r, opt.threads, bounds)
        R = minimize(F, initVal, mode="quick")
        results.append({'res': r, 'min': R})

    # Sort results by rss for comparisons
    sResults = sorted(results, key=lambda v:v['min'][0]['rss'])

    # Trying other minima for those above the cutoff
    rssTrustCutoff = 10 * optRss
    for r in sResults:
      if r['min'][0]['rss'] > rssTrustCutoff:
        print("Residue", r['res'], "is above the cutoff trying different minima")
        for initVal in results[0]['min'][1:]:
          F = sdfMinFunc(models, LarmorF, inData, inWeights, r['res'], opt.threads, bounds)
          R = minimize(F, initVal, mode="quick")
          if R['rss'] <= rssTrustCutoff:
            r['min'] = R
            break

    # Outward moving grid search for those above the cutoff
    for r in sResults:
      if r['min'][0]['rss'] > rssTrustCutoff:
        print("Residue", r['res'], "is above the cutoff trying grid search")
        F = sdfMinFunc(models, LarmorF, inData, inWeights, r['res'], opt.threads, bounds)
        R = minimize(F, r['min'][0]['par'], mode="extend", compare=rssTrustCutoff)
        if len(R) > 0:
          r['min'] = R

    print(results)


    # for r in residues:
    #   # bounds = np.array([ [1e-12,   1e-8]
    #   #                   , [1e-5,    0.99999]
    #   #                   ])

    #   minFunc = sdfMinFunc(
    #     models, LarmorF, inData, inWeights, resid, opt.threads)

    #   result = None
    #   sB  = []
    #   tB  = []
    #   S   = None
    #   Tau = None

    #   minimize    = Minimize()


    #   result      = minimize(
    #     minFunc,
    #     inVals,
    #     np.array([2e-10, 0.2]),
    #     bounds
    #     )


    #   print("{}:".format(resid), result)


