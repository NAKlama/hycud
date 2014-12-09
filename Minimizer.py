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

import math
import numpy as np
from copy import copy, deepcopy

import scipy.linalg as linalg
import scipy.optimize as optimize


class sdfRSS:
  def __init__(self, func, verbose=False):
    self.minF = func
    self.verb = verbose

  def __call__(self, par):
    res = self.minF(par)
    rss = np.sum(res**2)
    if self.verb:
      print(par, " ", end='')
      print(res, " ", end='')
      print(rss)

    return rss


class Simplex:
  def __init__(self, directions=3):
    self.dirN = directions
    if self.dirN < 3:
      self.dirN = 3
    if self.dirN % 2 == 0:
      self.dirN += 1
    self.dir = []
    angles = np.arange(0.0, 2 * np.pi, (2*np.pi)/self.dirN)
    for a in np.nditer(angles):
      rotM = np.array([ [np.cos(a), -np.sin(a)],
                        [np.sin(a), np.cos(a)] ])
      self.dir.append(np.dot(rotM, np.array([0.0, 1.0])))

  def __call__(self, func, par, step, bounds, rss, diffIn):
    points = []
    for d in self.dir:
      p = par + d * step
      if np.all(np.greater_equal(p, bounds[:,0])) and np.all(np.less_equal(p, bounds[:,1])):
        points.append({'par': p, 'rss': None, 'st': step})

    for p in points:
      diff = func(p['par'])
      p['dVect']  = diff
      p['rss']    = np.sum(diff**2)
    points.append({'par': par, 'rss': rss, 'st': step/2, 'dVect': diffIn})
    return points

class GaussNewton:
  def __init__(self, dampMult=2):
    self.dampMult = dampMult

  def __call__(self, func, par, step, bounds, rssIn, diffIn, damp=0.5):
    points  = []
    dampSt  = copy(damp)
    Dir     = [np.array([x,0]) for x in (-1,1)]
    Dir    += [np.array([0,x]) for x in (-1,1)]
    rssV    = []
    # mat   = []
    minDamp = damp
    minRss  = 1e300
    for d in Dir:
      rssV.append(func(par + d * step))
    dX        = rssV[1] - rssV[0]
    dY        = rssV[3] - rssV[2]
    J         = np.matrix([dX, dY])
    JJT       = J * J.T

    invJ      = linalg.inv(JJT + damp * np.diag(JJT))
    matS      = (invJ * J)
    diff      = np.array(np.dot(matS, diffIn))[0]
    newPar    = par - diff * step
    if np.all(np.greater_equal(newPar, bounds[:,0])) and np.all(np.less_equal(newPar, bounds[:,1])):
      res       = func(newPar)
      rss       = np.sum(res**2)
      # print("damp =", damp, "  rss =", rss)
      minDamp   = damp
      minRss    = rss
      points.append({
        'par':    newPar,
        'rss':    rss,
        'st':     step,
        'dVect':  res,
        'damp':   damp})


    # for i in [0.25, 0.5, 2, 4, 8]:
    for i in [0.5, 2, 4]:
      damp    = dampSt * i
      invJ    = linalg.inv(JJT + damp * np.diag(JJT))
      matS    = (invJ * J)
      diff    = np.array(np.dot(matS, diffIn))[0]
      newPar  = par - diff * step
      if np.all(np.greater_equal(newPar, bounds[:,0])) and np.all(np.less_equal(newPar, bounds[:,1])):
        res     = func(newPar)
        rss     = np.sum(res**2)
        if rss < minRss:
          minRss  = rss
          minDamp = damp
        # print("damp =", damp, "  rss =", rss)
        points.append({
          'par':    newPar,
          'rss':    rss,
          'st':     step,
          'dVect':  res,
          'damp':   damp})
          # 'damp':   (damp + dampSt)/2})
        if rss < rssIn:
          break


    # for s in [0.25, 0.5, 2, 4]:
    #   newPar  = par - diff * s * step
    #   res     = func(newPar)
    #   if np.all(np.greater_equal(newPar, bounds[:,0])) and np.all(np.less_equal(newPar, bounds[:,1])):
    #     points.append({
    #       'par':    newPar,
    #       'rss':    np.sum(res**2),
    #       'st':     step*s,
    #       'dVect':  res,
    #       'damp':   minDamp})
    points.append({
      'par':    par,
      'rss':    rssIn,
      'st':     step/4,
      'dVect':  diffIn,
      'damp':   minDamp})
    return points



class Minimize:
  def __init__(self, useCount=0, remFact=None, directions=3, detailed=False):
    self.count   = useCount
    if useCount < 1:
      self.count = 0
    self.directions = directions
    self.remFact = remFact
    self.detail  = detailed

  def __call__(self, func, initVal, step, bounds, tol=1e-6):
    points = []
    diff   = None
    # First loop
    print("Starting grid search")
    v = []
    for x in initVal[0]:
      for y in initVal[1]:
        v = np.array([x, y])
        d = func(v)
        r = np.sum(d**2)
        print(v, "\t", d, "  \t", "{:>14.6f}".format(r), "   \t", step, )
        points.append({'par': v, 'rss': r, 'st': copy(step), 'dVect': d})
    # for i in range(shape[1]**2):
    #   x = i % shape[1]
    #   y = (i - x) / shape[1]
    #   inVal = np.array([initVal[0][x], initVal[1][y]])
    #   d = func(inVal)
    #   r = np.sum(d**2)
    #   print(inVal, "\t", d, "\t", r, "   \t", step, )
    #   points.append({'par': inVal, 'rss': r, 'st': copy(step), 'dVect': d})
    print("Intial Loop done, entering simplex loop.")

    # minRoutine = Simplex(self.directions)
    minRoutine = GaussNewton()

    # Main loop
    # print(points)
    points = sorted(points, key=lambda v:v['rss'])
    count = len(points)
    # stepCurr = step
    acc   = 1e-1
    # firstRun = True
    while diff is None or diff > tol:
      rCount = count
      maxN   = 1
      if len(points) >= 5:
        rCount = 0
        avg = 0.0
        for p in points:
          avg += p['rss']
        avg /= len(points)
        print("avg(rss)    =", avg, end='')
        stdDev = 0.0
        for p in points:
          stdDev += (p['rss']-avg)**2
        stdDev /= len(points)
        stdDev = math.sqrt(stdDev)
        print(" stdDev(rss) =", stdDev)
        cutoff = copy(avg)
        # if len(points) <= 30:
        cutoff += stdDev
        for p in points:
          if p['rss'] <= cutoff:
            rCount += 1
      else:
        if self.remFact is not None:
          rCount = 0
          for p in points:
            if p['rss'] < points[0]['rss'] * self.remFact:
              rCount += 1
      if count == 0:
        count = rCount
      count   = min(count, rCount)

      iPoints = points[0:count]
      # begPoints = copy(iPoints)
      print("\nInput points: ({})".format(count))
      for p in iPoints:
        # p['st'] = copy(stepCurr)
        print("p: {}  d: {}  \trss: {:14.8f}  st: {}".format(
          p['par'], p['dVect'], p['rss'], p['st']))

      points  = []

      if self.detail:
        # Simplex
        rss = 1e100
        # first = True
        i = 0
        oPts = copy(iPoints)
        # while first or 1.0-(rss - points[0]['rss']) / rss >= max(acc, tol):
        for i in range(maxN):
          # first = False
          print("+++ Levenberg-Marquard {} ({:e})+++".format(i, 1.0 - (rss - iPoints[0]['rss']) / rss))
          points = copy(oPts)
          rss = points[0]['rss']
          oPts = []
          for ptIn in points:
            pt  = deepcopy(ptIn)
            lSt = pt['st']
            pnts = [pt]
            n = 0
            # while pnts[0]['rss'] >= pt['rss'] and n < 15 and np.all(np.greater((pnts[0]['st'] / pnts[0]['par']), 1e-10)):
            while pnts[0]['rss'] >= ptIn['rss'] and np.any(np.greater((np.absolute(pnts[0]['st'] / pnts[0]['par'])), 1e-12)):
              pt = pnts[0]
              if 'damp' in pt:
                pnts = sorted(minRoutine(func, pt['par'], lSt, bounds, pt['rss'], pt['dVect'], pt['damp']),
                            key=lambda v:v['rss'])
              else:
                pnts = sorted(minRoutine(func, pt['par'], lSt, bounds, pt['rss'], pt['dVect']),
                            key=lambda v:v['rss'])
              if not pnts[0]['rss'] >= pt['rss']:
                print(chr(27) + "[1m", end='')
                print(pnts[0]['par'], "\t", pnts[0]['dVect'], "   \t{:14.6f}".format(pnts[0]['rss']), "   \t", pnts[0]['st'], end='')
                if 'damp' in pnts[0]:
                  print("\t", pnts[0]['damp'])
                else:
                  print("")
                print(chr(27) + "[m", end='')
              else:
                print(pnts[0]['par'], "\t", pnts[0]['dVect'], "   \t{:14.6f}".format(pnts[0]['rss']), "   \t", pnts[0]['st'], end='')
                if 'damp' in pnts[0]:
                  print("\t", pnts[0]['damp'])
                else:
                  print("")
              lSt = pnts[0]['st']
              n+=1
              if n >= maxN:
                break
            oPts.append(pnts[0])
            rss = pnts[0]['rss']

          points = sorted(oPts, key=lambda v:v['rss'])
          i += 1
          # print("{:e} {} {} {}".format((rss - iPoints[0]['rss']) / rss, rss, iPoints[0]['rss'], max(acc, tol)))
        print("")

        oPoints = sorted(points, key=lambda v:v['rss'])
        W = 0.0
        D = 0.0
        for i in range(len(oPoints)):
          w  = len(oPoints) - i
          d  = ((iPoints[i]['rss'] - oPoints[i]['rss']) / iPoints[i]['rss']) * float(w)
          W += w
          D += d
          # print("{} * {}".format(d, w))

        diff = float(D)/float(W)
        # diff = (begPoints[0]['rss'] - oPoints[0]['rss']) / begPoints[0]['rss']
        print("Diff = {:14.8f}".format(diff))
        points = oPoints
        acc /= 2
        # for p in points:
          # p['st'] *= 2
        if self.remFact is None:
          count = int(float(count) * 0.75)
        # if count < 2:
        #   count = 2
        maxN += 1
        if maxN > 5:
          maxN = 5
      else:
        resultList = []
        for p in iPoints:
          print("Starting optimization on", p['par'])
          minFunc = sdfRSS(func)
          result  = optimize.minimize(
            minFunc,
            np.array(p['par']),
            # method="BFGS",
            method="L-BFGS-B",
            # method="TNC",
            # method="SLSQP",
            bounds=[(1e-14, 1e-5), (0.0001, 0.9999)],
            # tol=1e-4,
            options={'disp': True, 'maxiter': 20}
            )
          res = func(result.x)
          rss = np.sum(res**2)
          points.append({ 'par':    result.x
                        , 'dVect':  res
                        , 'rss':    rss
                        , 'st':     [0,0]
                        })
          print("Minimization Result:", result.x, res, rss)





    print("\nResult:")
    for p in points:
      print(p['par'], "\t", p['rss'], "   \t", p['st'])
    return points[0]['par']









