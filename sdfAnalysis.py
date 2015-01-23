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
# import math
import copy
import sys
# import time
from multiprocessing import Pool
from sdfFunction import sdfFunction

# Linear algebra functions from scipy are using BLAS and are faster
# so we try them before using numpy as a fallback
# try:
#   import scipy.linalg as linalg
# except ImportError:
#   import numpy.linalg as linalg

from StatItem import npStatItem
from ANSI import *

class sdfAnalysis:
  def __init__(self, opt, models, wH):
    self.models   = models
    self.size     = opt.fragSize
    self.resSDF   = []
    self.verbose  = opt.verbose
    self.threads  = opt.threads
    self.tau      = opt.SDFtau
    self.order    = opt.SDForder
    if type(wH) == type(3.2):
      self.wH     = [wH]
    elif type(wH) == type([3.2]):
      self.wH     = wH
    else:
      print("ERROR!!!")
      sys.exit(1)

  def calc(self):
    sdf_pool = Pool(processes=self.threads)
    task     = sdfFunction(self.wH)
    task.changeParameters((self.tau, self.order))
    result   = sdf_pool.map(task, self.models.models)
    for i in range(len(result)):
      m = self.models.models[i]
      for f in m.fragments:
        f.sdf = []
        for r in result[i]:
          if r['frag'] == f.num:
            f.sdf.append(r)
    # for m in self.models.models:
    #   for frag in m.fragments:
    #     if not frag.partial:
    #       frag.sdf   = []
    #       diffRaw = []
    #       for v in frag.diffMat[3:6]:
    #         diffRaw.append(v[3:6])
    #       diff    = np.array(diffRaw, np.float64) / frag.weightingFactor
    #       if diff.shape != (3,3):
    #         print("ERROR: diffusion matrix does not have the right shape.", diff.shape)
    #         time.sleep(0.01)
    #       Dc, Qc  = linalg.eig(diff)
    #       D       = np.real(Dc)
    #       Q       = np.real(Qc)
    #       Diso    = (D[0] + D[1] + D[2]) / 3.0
    #       L2      = (D[0] * D[1] + D[0] * D[2] + D[1] * D[2]) / 3.0
    #       delta   = (D - Diso) / (math.sqrt(Diso**2 - L2))
    #       tau     = np.array( [ 1.0 / (4 * D[0] + D[1] + D[2])
    #                         ,   1.0 / (D[0] + 4 * D[1] + D[2])
    #                         ,   1.0 / (D[0] + D[1] + 4 * D[2])
    #                         ,   1.0 / (6 * Diso + 6 * math.sqrt(Diso**2 - L2))
    #                         ,   1.0 / (6 * Diso - 6 * math.sqrt(Diso**2 - L2))
    #                         ])
    #       tt      = tau**2
    #       pi      = math.pi
    #       Gn      = -2.712e7              # rad/(T s)
    #       Gh      = 2.6752e8              # rad/(T s)
    #       rhn     = 1.02e-10              # m
    #       h       = 6.62606896e-34        # J s
    #       u0      = 4 * pi * 1e-7         # uo = Bo/H; H = A/m; uo = T m/A
    #       csaN    = -170e-6                  # CSA of 15N
    #       dd      = -math.sqrt(1/10) * (u0/(4*pi))*(h/(2*pi))*Gn*Gh / rhn**3
    #       for wh in self.wH:
    #         # wh      = 600.25e6              # 1H Larmor frequency
    #         B0      = wh * 2e6* pi / Gh
    #         cc      = -math.sqrt(2/15) * csaN * B0 * Gn
    #         wn      = wh * -(Gn / Gh)
    #         f       = np.array([0, wn, wh-wn, wh, wh+wn])
    #         w       = (2e6 * pi * f)**2
    #         q       = np.matrix(w).T * np.matrix(tt)
    #         one     = np.ones((5,5))
    #         j       = np.tile(tau.T, (5,1)) * np.array(one / (one + q))
    #         for r in frag.residues:
    #           SRaw  = frag.protein.getNHvect(r)
    #           if SRaw is not None:
    #             S     = np.matrix(SRaw) * 1e-10
    #             QS    = np.matrix(Q).H * S.T
    #             x     = QS[0,0] / linalg.norm(QS)
    #             y     = QS[1,0] / linalg.norm(QS)
    #             z     = QS[2,0] / linalg.norm(QS)
    #             # print(QS, "  x",x," y",y," z",z)
    #             A     = np.matrix(  [ 3 * y**2 * z**2
    #                               ,   3 * x**2 * z**2
    #                               ,   3 * x**2 * y**2
    #                               ,    (1/4) * (3*(x**4 + y**4 + z**4)-1)\
    #                                 - (1/12) * (delta[0]*(3* x**4 +6* y**2 * z**2 -1)
    #                                            +delta[1]*(3* y**4 +6* x**2 * z**2 -1)
    #                                            +delta[2]*(3* z**4 +6* x**2 * y**2 -1))
    #                               ,    (1/4) * (3*(x**4 + y**4 + z**4)-1)\
    #                                 + (1/12) * (delta[0]*(3* x**4 +6* y**2 * z**2 -1)
    #                                            +delta[1]*(3* y**4 +6* x**2 * z**2 -1)
    #                                            +delta[2]*(3* z**4 +6* x**2 * y**2 -1))
    #                               ]).T
    #             # print("A",A)
    #             J     = np.matrix(j) * A
    #             R1dd  = dd**2       * (3*(J[1,0]+6*J[2,0]+J[4,0]))
    #             R1c   = (cc**2)     * J[1,0]
    #             R1    = R1dd + R1c
    #             R2dd  = (dd**2 / 2) * (4*J[0,0]+3*J[1,0]+6*J[2,0]+6*J[3,0]+J[4,0])
    #             R2c   = (cc**2 / 6) * (4*J[0,0]+3*J[1,0])
    #             R2    = R2dd + R2c
    #             sigma = dd**2 * (6*J[2,0]-J[4,0])
    #             noe   = 1+(sigma*Gh / (R1*Gn))
    #             ratio = R2/R1
    #             R     = np.array([R1, R2])
    #             frag.sdf.append({ 'res':    r
    #                             , 'wh':     wh
    #                             , 'J':      (J[:,0] * 1e9)
    #                             , 'R':      R
    #                             , 'A':      A
    #                             , 'tau':    tau
    #                             , 'sigma':  sigma
    #                             , 'noe':    noe
    #                             , 'ratio':  ratio
    #                             , 'fragC':  frag.adjFragCount
    #                             })

  def average(self):
    # print("Average")
    resSDF = self.resSDF
    for r in self.models.residues:
      for wh in self.wH:
        resSDF.append({ 'res':    r
                      , 'wh':     wh
                      , 'J':      None
                      , 'R':      None
                      , 'sigma':  None
                      , 'noe':    None
                      , 'ratio':  None
                      , 'tau':    None
                      , 'A':      None
                      , 'fragC':  0
                      })
      # print(self.models.residues)
    for m in self.models.models:
      m.resSDF = copy.deepcopy(resSDF)
      for f in m.fragments:
        if not f.partial:
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
                  sdfAvg['A']     = npStatItem()
                sdfAvg['J'].addValue(sdf['J'])
                sdfAvg['R'].addValue(sdf['R'])
                sdfAvg['sigma'].addValue(np.array([sdf['sigma']]))
                sdfAvg['noe'].addValue(np.array([sdf['noe']]))
                sdfAvg['ratio'].addValue(np.array([sdf['ratio']]))
                sdfAvg['tau'].addValue(np.array([sdf['tau']]))
                sdfAvg['A'].addValue(np.array([sdf['A']]))
                sdfAvg['fragC'] += sdf['fragC']

    for sdf in resSDF:
      if sdf['J'] is None:
        sdf['J']      = npStatItem()
        sdf['R']      = npStatItem()
        sdf['sigma']  = npStatItem()
        sdf['noe']    = npStatItem()
        sdf['ratio']  = npStatItem()
        sdf['tau']    = npStatItem()
        sdf['A']    = npStatItem()
      for m in self.models.models:
        for msdf in m.resSDF:
          if msdf['J'] is not None and sdf['res'] == msdf['res']:
            sdf['J'].addStatItem(msdf['J'])
            sdf['R'].addStatItem(msdf['R'])
            sdf['sigma'].addStatItem(msdf['sigma'])
            sdf['noe'].addStatItem(msdf['noe'])
            sdf['ratio'].addStatItem(msdf['ratio'])
            sdf['tau'].addStatItem(msdf['tau'])
            sdf['A'].addStatItem(msdf['A'])
            sdf['fragC'] += msdf['fragC']

    self.fragCmax = 0
    for sdf in resSDF:
      if sdf['fragC'] > self.fragCmax:
        self.fragCmax = sdf['fragC']

  def output(self):
    csi = chr(27) + '['
    twoTables = True
    (width, height) = getTerminalSize()
    if width > 217:
      twoTables = False
    if self.verbose > 0:
      print("Colors show confidence by showing what percentage of max. possible fragments")
      print("where taken into account:")
      print("{}  0% -  24.999%".format(getColor(0.1)))
      print("{} 25% -  49.999%".format(getColor(0.3)))
      print("{} 50% -  79.999%".format(getColor(0.6)))
      print("{} 80% -  98.999%".format(getColor(0.85)))
      print("{} 99% - 100%".format(getColor(1.0)))
      print("\n", resetColor())

    print(" Res  | J(w=0)               J(w=wN)              "
         +"J(w=wH+wN)           J(w=wH)              J(w=wH-wN)"
         , end="")
    if not twoTables:
      print("           R1                   R2                   "
           +"sigma                noe                  R2/R1")
    else:
      print("")
    print("------+----------------------------------------------"
         +"----------------------------------------------------------"
         , end="")
    if not twoTables:
      print("------------------------------------------------"
           +"---------------------------------------------------------")
    else:
      print("")

    for sdf in self.resSDF:
      if sdf['J'].weight > 0.0:
        count     = sdf['J'].weight / len(self.models.models)
        ratio     = count / self.size
        fragRatio = sdf['fragC'] / self.fragCmax
        print("{:>5d} | ".format(sdf['res']), end='')
        avg     = sdf['J'].getAvg()[0]
        stdDev  = sdf['J'].getStdDev()[0]
        # print(avg)
        # print(stdDev)
        for i in range(0, avg.shape[0]):
          print("{}±{}  ".format(
            coloredExp(avg[i], fragRatio), coloredExp(stdDev[i], fragRatio)), end='')
        if not twoTables:
          avg     = sdf['R'].getAvg()
          stdDev  = sdf['R'].getStdDev()
          print("{}±{}  ".format(
            coloredExp(avg[0], fragRatio), coloredExp(stdDev[0], fragRatio)), end='')
          print("{}±{}  ".format(
            coloredExp(avg[1], fragRatio), coloredExp(stdDev[1], fragRatio)), end='')
          print("{}±{}  ".format(
            coloredExp(sdf['sigma'].getAvg()[0], fragRatio),
            coloredExp(sdf['sigma'].getStdDev()[0], fragRatio)), end='')
          print("{}±{}  ".format(
            coloredExp(sdf['noe'].getAvg()[0], fragRatio),
            coloredExp(sdf['noe'].getStdDev()[0], fragRatio)), end='')
          print("{}±{}  ".format(
            coloredExp(sdf['ratio'].getAvg()[0], fragRatio),
            coloredExp(sdf['ratio'].getStdDev()[0], fragRatio)), end='')

          if self.verbose > 2:
            print("ratio: {:.3f}  fragC: {:5d} ({:.3f}) ratio: {:.6f})".format(
              ratio, sdf['fragC'], fragRatio,
              (sdf['fragC'] / self.fragCmax) / (count / self.size)), end="")
          if self.verbose > 3:
            print(" w: ", sdf['J'].weight,
                  " m: ", len(self.models.models),
                  " c: {:.3f}".format(count),
                  end="")
        print()

    if twoTables:
      print("\n{}0m".format(csi))
      print(" Res  | R1                   R2                   "
           +"sigma                noe                  R2/R1")
      print("------+-------------------------------------------"
           +"-------------------------------------------------------------")
      for sdf in self.resSDF:
        if sdf['J'].weight > 0.0:
          count  = sdf['J'].weight / len(self.models.models)
          ratio  = count / self.size
          fragRatio = sdf['fragC'] / self.fragCmax
          print("{:>5d} | ".format(sdf['res']), end='')
          avg     = sdf['R'].getAvg()
          stdDev  = sdf['R'].getStdDev()
          print("{}±{}  ".format(
            coloredExp(avg[0], fragRatio), coloredExp(stdDev[0], fragRatio)), end='')
          print("{}±{}  ".format(
            coloredExp(avg[1], fragRatio), coloredExp(stdDev[1], fragRatio)), end='')
          print("{}±{}  ".format(
            coloredExp(sdf['sigma'].getAvg()[0], fragRatio),
            coloredExp(sdf['sigma'].getStdDev()[0], fragRatio)), end='')
          print("{}±{}  ".format(
            coloredExp(sdf['noe'].getAvg()[0], fragRatio),
            coloredExp(sdf['noe'].getStdDev()[0], fragRatio)), end='')
          print("{}±{}  ".format(
            coloredExp(sdf['ratio'].getAvg()[0], fragRatio),
            coloredExp(sdf['ratio'].getStdDev()[0], fragRatio)), end='')
          if self.verbose > 2:
            print("ratio:", ratio, end="")

          print()


    print("\n{}0m".format(csi))
    # print(" Res  | R1                   R2                   "
    #      +"sigma                noe                  R2/R1")
    # print("------+-------------------------------------------"
    #      +"-------------------------------------------------------------")
    for sdf in self.resSDF:
      if sdf['A'].weight > 0.0:
        count  = sdf['A'].weight / len(self.models.models)
        ratio  = count / self.size
        fragRatio = sdf['fragC'] / self.fragCmax
        print("{:>5d} | ".format(sdf['res']), end='')
        avg     = sdf['A'].getAvg()[0]
        stdDev  = sdf['A'].getStdDev()[0]
        # print(avg)
        # print(stdDev)
        for i in range(0, avg.shape[0]):
          print("{}±{}  ".format(
            coloredExp(avg[i][0], fragRatio), coloredExp(stdDev[i][0], fragRatio)), end='')
        print("|  ", end='')
        avg     = sdf['tau'].getAvg()
        stdDev  = sdf['tau'].getStdDev()
        # print(avg)
        # print(stdDev)
        for i in range(0, avg.shape[0]):
          print("{}±{}  ".format(
            coloredExp(avg[i][0], fragRatio), coloredExp(stdDev[i][0], fragRatio)), end='')

        print()




