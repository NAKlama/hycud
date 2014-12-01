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
import math

try:
  import scipy.linalg as linalg
except ImportError:
  import numpy.linalg as linalg


class sdfFunction:
  def __init__(self, wH):
    if type(wH) == type(3.2):
      self.wH     = np.array([wH])
    elif type(wH) == type([3.2]):
      self.wH     = np.array(wH)
    else:
      self.wH     = wH
    self.tau            = 0.0
    self.S              = 1.0
    self.parChange      = False
    self.singleResidue  = None

  def singleResidueSet(self, res):
    self.singleResidue  = res

  def changeParameters(self, par):
    self.parChange = True
    self.tau = par[0]
    self.S   = par[1]

  def __call__(self, m):
    sdf   = []
    for frag in m.fragments:
      if not frag.partial:
        if self.singleResidue is None or self.singleResidue in frag.residues:
          diffRaw = []
          for v in frag.diffMat[3:6]:
            diffRaw.append(v[3:6])
          diff    = np.array(diffRaw, np.float64) / frag.weightingFactor
          Dc, Qc  = linalg.eig(diff)
          D       = np.real(Dc)
          Q       = np.real(Qc)
          Diso    = (D[0] + D[1] + D[2]) / 3.0
          L2      = (D[0] * D[1] + D[0] * D[2] + D[1] * D[2]) / 3.0
          delta   = (D - Diso) / (math.sqrt(Diso**2 - L2))
          tau     = np.array( [ 1.0 / (4 * D[0] + D[1] + D[2])
                            ,   1.0 / (D[0] + 4 * D[1] + D[2])
                            ,   1.0 / (D[0] + D[1] + 4 * D[2])
                            ,   1.0 / (6 * Diso + 6 * math.sqrt(Diso**2 - L2))
                            ,   1.0 / (6 * Diso - 6 * math.sqrt(Diso**2 - L2))
                            ])
          tt      = tau**2
          pi      = math.pi
          Gn      = -2.712e7              # rad/(T s)
          Gh      = 2.6752e8              # rad/(T s)
          rhn     = 1.02e-10              # m
          h       = 6.62606896e-34        # J s
          u0      = 4 * pi * 1e-7         # uo = Bo/H; H = A/m; uo = T m/A
          csaN    = -170e-6                  # CSA of 15N
          dd      = -math.sqrt(1/10) * (u0/(4*pi))*(h/(2*pi))*Gn*Gh / rhn**3
          for wh in self.wH:
            # wh      = 600.25e6              # 1H Larmor frequency
            B0      = 2e6 * wh * pi / Gh
            cc      = -math.sqrt(2/15) * csaN * B0 * Gn
            wn      = wh * -(Gn / Gh)
            f       = np.array([0, wn, wh-wn, wh, wh+wn])
            w       = (2e6 * pi * f)**2
            q       = np.matrix(w).T * np.matrix(tt)
            one     = np.ones((5,5))
            j       = self.S * np.tile(tau.T, (5,1)) * np.array(one / (one + q))
            if self.parChange:
              tauDash  = 1 / (1 / tau + 1 / self.tau)
              qDash    = np.matrix(w).T * np.matrix(tauDash**2)
              j       += (1.0 - self.S) * np.tile(tauDash.T, (5,1)) * np.array(one / (one + qDash))
            for r in frag.residues:
              if self.singleResidue is None or self.singleResidue in frag.residues:
                SRaw  = frag.protein.getNHvect(r)
                if SRaw is not None:
                  S     = np.matrix(SRaw) * 1e-10
                  QS    = np.matrix(Q).H * S.T
                  x     = QS[0,0] / linalg.norm(QS)
                  y     = QS[1,0] / linalg.norm(QS)
                  z     = QS[2,0] / linalg.norm(QS)
                  # print(QS, "  x",x," y",y," z",z)
                  A     = np.matrix(  [ 3 * y**2 * z**2
                                    ,   3 * x**2 * z**2
                                    ,   3 * x**2 * y**2
                                    ,    (1/4) * (3*(x**4 + y**4 + z**4)-1)\
                                      - (1/12) * (delta[0]*(3* x**4 +6* y**2 * z**2 -1)
                                                 +delta[1]*(3* y**4 +6* x**2 * z**2 -1)
                                                 +delta[2]*(3* z**4 +6* x**2 * y**2 -1))
                                    ,    (1/4) * (3*(x**4 + y**4 + z**4)-1)\
                                      + (1/12) * (delta[0]*(3* x**4 +6* y**2 * z**2 -1)
                                                 +delta[1]*(3* y**4 +6* x**2 * z**2 -1)
                                                 +delta[2]*(3* z**4 +6* x**2 * y**2 -1))
                                    ]).T
                  # print("A",A)
                  J     = np.matrix(j) * A
                  R1dd  = dd**2       * (3*(J[1,0]+6*J[2,0]+J[4,0]))
                  R1c   = (cc**2)     * J[1,0]
                  R1    = R1dd + R1c
                  R2dd  = (dd**2 / 2) * (4*J[0,0]+3*J[1,0]+6*J[2,0]+6*J[3,0]+J[4,0])
                  R2c   = (cc**2 / 6) * (4*J[0,0]+3*J[1,0])
                  R2    = R2dd + R2c
                  sigma = dd**2 * (6*J[2,0]-J[4,0])
                  noe   = 1+(sigma*Gh / (R1*Gn))
                  ratio = R2/R1
                  R     = np.array([R1, R2])
                  sdf.append({ 'res':    r
                             , 'frag':   frag.num
                             , 'wh':     wh
                             , 'J':      (J[:,0] * 1e9)
                             , 'R':      R
                             , 'A':      A
                             , 'tau':    tau
                             , 'sigma':  sigma
                             , 'noe':    noe
                             , 'ratio':  ratio
                             , 'fragC':  frag.adjFragCount
                             })
    return sdf


