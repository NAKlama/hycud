#!/usr/bin/python3

# HYCUD
# Copyright (C) 2014 Klama, Nina Alexandra and Rezaei-Ghaleh, Nasrollah
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
import sys
import numpy as np
from copy import copy

# import scipy.linalg as linalg
import scipy.optimize as optimize


class Minimize:
  def __init__(self, maxResult=10, gridSpec=None, verbose=True):
    self.gridSpec   = gridSpec
    self.maxResult  = maxResult
    self.enableGrid = False
    self.verbose    = verbose

    # Calculate exact grid
    self.grid      = []

    gsTau     = self.gridSpec[0]
    gsS       = self.gridSpec[1]
    if len(gsTau) > 1 and len(gsS) > 1:
      self.enableGrid = True
      countTau  = 5
      countS    = 5
      if len(gsTau) > 2:
        countTau = int(gsTau[2])
      if len(gsS) > 2:
        countS = int(gsS[2])
      minTau    =  gsTau[0] - gsTau[1]
      maxTau    = (gsTau[0] + gsTau[1]) * (1+ (1/ (2*countTau)))
      minS      =  gsS[0] - gsS[1]
      maxS      = (gsS[0] + gsS[1])     * (1+ (1/ (2*countS)))
      tau       = np.arange(minTau, maxTau, (gsTau[1] * 2.0) / countTau)
      S         = np.arange(minS,   maxS,   (gsS[1] * 2.0)   / countS)
      for t in np.nditer(tau):
        for s in np.nditer(S):
          self.grid.append( np.array([t, s]) )

      self.dTau = tau[1] - tau[0]
      self.dS   = S[1]   - S[0]
      self.bounds = [ [minTau, maxTau], [minS, maxS] ]



  def levMar(self, func, initVal, limit=None):
    res = None
    if limit is not None:
      res = optimize.leastsq(func, initVal, maxfev=limit)
    else:
      res = optimize.leastsq(func, initVal)
    # print(res)
    diff = func(res[0])
    rss  = np.sum(diff ** 2)
    return (res[1], res[0], diff, rss)

  def weeder(self, points):
    avg    = 0.0
    stdDev = 0.0

    # Calculate Average
    for p in points:
      avg += p['rss']
    avg /= len(points)

    # Calculate StdDev
    for p in points:
      stdDev += (p['rss'] - avg)**2
    stdDev /= len(points)
    stdDev  = math.sqrt(stdDev)

    # Making output list
    out = []
    for p in points:
      if p['rss'] <= avg + stdDev:
        out.append(p)
    return out

  def gridSearch(self, func, grid):
    points = []
    for g in grid:
      diff = func(g)
      rss  = np.sum(diff**2)
      points.append({ 'par':    g
                    , 'diff':   diff
                    , 'rss':    rss
                    })
    return points


  def __call__(self, func, initVal, mode="quick", compare=None):
    # Modes with increasing complexity
    # quick   Only do LevMar
    # extend  Used after quick, initVal is assumed to be optimum found by quick
    #         will do an extending grid search
    # grid    Does the full grid search

    results   = []
    points    = []
    if mode == "grid" and self.enableGrid:

      # Doing initial full grid search
      points = self.gridSearch(func, self.grid)

      levMarSteps = 8

      while len(points) > 0:
        # Sorting with increasing rss
        points  = sorted(points, key=lambda v:v['rss'])
        if self.verbose:
          for p in points:
            print("par:", p['par'], "  diff:", p['diff'], "  rss:", p['rss'])

        # Spotting Clusters
        clusterCutOff = 0.01
        n = 0
        while n < len(points):
          val = points[n]['par']
          out = points[:n+1]
          for r in points[n+1:]:
            ratio = np.fabs((r['par'] - val) / val)
            if np.any(np.greater(ratio, clusterCutOff)):
              out.append(r)
          points = out
          n += 1

        # Weeding out worst results
        work    = self.weeder(points)
        points  = []

        for p in work:
          (ier, parOut, diff, rss) = self.levMar(func, p['par'], levMarSteps)

          if ier >=1 and ier <= 4:
            print("Result: par:", parOut, "  diff:", diff, "  rss:", rss)
            results.append( { 'par':    parOut
                            , 'diff':   diff
                            , 'rss':    rss
                            })
          else:
            points.append({ 'par':    parOut
                          , 'diff':   diff
                          , 'rss':    rss
                          })

        # Increase maximum LevMar function calls
        levMarSteps += 4

    # Extending grid Mode
    elif mode == "extend" and self.enableGrid:
      if compare is None:
        sys.exit(0)
      bestRss = compare * 1e50

      step = np.array([self.dTau, self.dS])

      n = 1
      while bestRss > compare:
        grid = []
        p = initVal

        # Faces first
        for d in -n, n:
          for i in range(-(n-1), n):
            grid.append( np.array([p[0] + i * step[0], p[1] + d * step[1]]) )
        for d in -n, n:
          for i in range(-(n-1), n):
            grid.append( np.array([p[0] + d * step[0], p[1] + i * step[1]]) )

        # Corners
        for d in (n,n), (-n,n), (n,-n), (-n,-n):
          grid.append( p + np.array(d) * step )

        # Cheching bounds
        workGrid = copy(grid)
        grid     = []
        for g in workGrid:
          if g[0] >= self.bounds[0][0] and g[0] <= self.bounds[0][1]:
            if g[1] >= self.bounds[1][0] and g[1] <= self.bounds[1][1]:
              grid.append(g)

        # Abort when all points are out of bounds
        if len(grid) == 0:
          bestRss = 0

        # Grid search and sort
        points = self.gridSearch(func, grid)
        outP   = []
        for p in points:
          (ier, parOut, diff, rss) = self.levMar(func, p['par'])
          if ier >=1 and ier <= 4:
            outP.append( { 'par':    parOut
                         , 'diff':   diff
                         , 'rss':    rss
                         })

        points = sorted(outP, key=lambda v:v['rss'])
        if bestRss > points[0]['rss']:
          bestRss = points[0]['rss']
          results.append(points[0])

      results = [sorted(results, key=lambda v:v['rss'])[0]]

    else: # "quick"
      (ier, parOut, diff, rss) = self.levMar(func, initVal)
      if ier >=1 and ier <= 4:
        results.append( { 'par':    parOut
                        , 'diff':   diff
                        , 'rss':    rss
                        })

    out = sorted(results, key=lambda v:v['rss'])
    if self.verbose:
      print("Results:")
      for p in results:
        print("par:", p['par'], "  diff:", p['diff'], "  rss:", p['rss'])
    return out








    # print("Starting grid search")
    # v = []
    # for x in initVal[0]:
    #   for y in initVal[1]:
    #     v = np.array([x, y])
    #     d = func(v)
    #     r = np.sum(d**2)
    #     print(v, "\t", d, "  \t", "{:>14.6f}".format(r), "   \t", step, )
    #     points.append({'par': v, 'rss': r, 'st': copy(step), 'dVect': d})
    # # for i in range(shape[1]**2):
    # #   x = i % shape[1]
    # #   y = (i - x) / shape[1]
    # #   inVal = np.array([initVal[0][x], initVal[1][y]])
    # #   d = func(inVal)
    # #   r = np.sum(d**2)
    # #   print(inVal, "\t", d, "\t", r, "   \t", step, )
    # #   points.append({'par': inVal, 'rss': r, 'st': copy(step), 'dVect': d})
    # print("Intial Loop done, entering simplex loop.")

    # # minRoutine = Simplex(self.directions)
    # minRoutine = GaussNewton()

    # # Main loop
    # # print(points)
    # points = sorted(points, key=lambda v:v['rss'])
    # count = len(points)
    # # stepCurr = step
    # acc   = 1e-1
    # # firstRun = True
    # while diff is None or diff > tol:
    #   rCount = count
    #   maxN   = 1
    #   if len(points) >= 5:
    #     rCount = 0
    #     avg = 0.0
    #     for p in points:
    #       avg += p['rss']
    #     avg /= len(points)
    #     print("avg(rss)    =", avg, end='')
    #     stdDev = 0.0
    #     for p in points:
    #       stdDev += (p['rss']-avg)**2
    #     stdDev /= len(points)
    #     stdDev = math.sqrt(stdDev)
    #     print(" stdDev(rss) =", stdDev)
    #     cutoff = copy(avg)
    #     # if len(points) <= 30:
    #     cutoff += stdDev
    #     for p in points:
    #       if p['rss'] <= cutoff:
    #         rCount += 1
    #   else:
    #     if self.remFact is not None:
    #       rCount = 0
    #       for p in points:
    #         if p['rss'] < points[0]['rss'] * self.remFact:
    #           rCount += 1
    #   if count == 0:
    #     count = rCount
    #   count   = min(count, rCount)

    #   iPoints = points[0:count]
    #   # begPoints = copy(iPoints)
    #   print("\nInput points: ({})".format(count))
    #   for p in iPoints:
    #     # p['st'] = copy(stepCurr)
    #     print("p: {}  d: {}  \trss: {:14.8f}  st: {}".format(
    #       p['par'], p['dVect'], p['rss'], p['st']))

    #   points  = []

    #   if self.detail:
    #     # Simplex
    #     rss = 1e100
    #     # first = True
    #     i = 0
    #     oPts = copy(iPoints)
    #     # while first or 1.0-(rss - points[0]['rss']) / rss >= max(acc, tol):
    #     for i in range(maxN):
    #       # first = False
    #       print("+++ Levenberg-Marquard {} ({:e})+++".format(i, 1.0 - (rss - iPoints[0]['rss']) / rss))
    #       points = copy(oPts)
    #       rss = points[0]['rss']
    #       oPts = []
    #       for ptIn in points:
    #         pt  = deepcopy(ptIn)
    #         lSt = pt['st']
    #         pnts = [pt]
    #         n = 0
    #         # while pnts[0]['rss'] >= pt['rss'] and n < 15 and np.all(np.greater((pnts[0]['st'] / pnts[0]['par']), 1e-10)):
    #         while pnts[0]['rss'] >= ptIn['rss'] and np.any(np.greater((np.absolute(pnts[0]['st'] / pnts[0]['par'])), 1e-12)):
    #           pt = pnts[0]
    #           if 'damp' in pt:
    #             pnts = sorted(minRoutine(func, pt['par'], lSt, bounds, pt['rss'], pt['dVect'], pt['damp']),
    #                         key=lambda v:v['rss'])
    #           else:
    #             pnts = sorted(minRoutine(func, pt['par'], lSt, bounds, pt['rss'], pt['dVect']),
    #                         key=lambda v:v['rss'])
    #           if not pnts[0]['rss'] >= pt['rss']:
    #             print(chr(27) + "[1m", end='')
    #             print(pnts[0]['par'], "\t", pnts[0]['dVect'], "   \t{:14.6f}".format(pnts[0]['rss']), "   \t", pnts[0]['st'], end='')
    #             if 'damp' in pnts[0]:
    #               print("\t", pnts[0]['damp'])
    #             else:
    #               print("")
    #             print(chr(27) + "[m", end='')
    #           else:
    #             print(pnts[0]['par'], "\t", pnts[0]['dVect'], "   \t{:14.6f}".format(pnts[0]['rss']), "   \t", pnts[0]['st'], end='')
    #             if 'damp' in pnts[0]:
    #               print("\t", pnts[0]['damp'])
    #             else:
    #               print("")
    #           lSt = pnts[0]['st']
    #           n+=1
    #           if n >= maxN:
    #             break
    #         oPts.append(pnts[0])
    #         rss = pnts[0]['rss']

    #       points = sorted(oPts, key=lambda v:v['rss'])
    #       i += 1
    #       # print("{:e} {} {} {}".format((rss - iPoints[0]['rss']) / rss, rss, iPoints[0]['rss'], max(acc, tol)))
    #     print("")

    #     oPoints = sorted(points, key=lambda v:v['rss'])
    #     W = 0.0
    #     D = 0.0
    #     for i in range(len(oPoints)):
    #       w  = len(oPoints) - i
    #       d  = ((iPoints[i]['rss'] - oPoints[i]['rss']) / iPoints[i]['rss']) * float(w)
    #       W += w
    #       D += d
    #       # print("{} * {}".format(d, w))

    #     diff = float(D)/float(W)
    #     # diff = (begPoints[0]['rss'] - oPoints[0]['rss']) / begPoints[0]['rss']
    #     print("Diff = {:14.8f}".format(diff))
    #     points = oPoints
    #     acc /= 2
    #     # for p in points:
    #       # p['st'] *= 2
    #     if self.remFact is None:
    #       count = int(float(count) * 0.75)
    #     # if count < 2:
    #     #   count = 2
    #     maxN += 1
    #     if maxN > 5:
    #       maxN = 5
    #   else:
    #     resultList = []
    #     for p in iPoints:
    #       print("Starting optimization on", p['par'])
    #       minFunc = sdfRSS(func)
    #       result  = optimize.minimize(
    #         minFunc,
    #         np.array(p['par']),
    #         # method="BFGS",
    #         method="L-BFGS-B",
    #         # method="TNC",
    #         # method="SLSQP",
    #         bounds=[(1e-14, 1e-5), (0.0001, 0.9999)],
    #         # tol=1e-4,
    #         options={'disp': True, 'maxiter': 20}
    #         )
    #       res = func(result.x)
    #       rss = np.sum(res**2)
    #       points.append({ 'par':    result.x
    #                     , 'dVect':  res
    #                     , 'rss':    rss
    #                     , 'st':     [0,0]
    #                     })
    #       print("Minimization Result:", result.x, res, rss)





    # print("\nResult:")
    # for p in points:
    #   print(p['par'], "\t", p['rss'], "   \t", p['st'])
    # return points[0]['par']









