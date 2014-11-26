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



# Linear algebra functions from scipy are using BLAS and are faster
# so we try them before using numpy as a fallback
try:
  import scipy.linalg as linalg
except ImportError:
  import numpy.linalg as linalg


def calcWeightFact(m):
  for fragI in m.fragments:
    i         = fragI.num
    commonSum = 0
    cj  = 10
    cj /= 36.12
    centerI = fragI.center.getCenter()
    for fragJ in m.fragments:
      j = fragJ.num
      if i != j:
        centerJ = fragJ.center.getCenter()
        rij = linalg.norm(centerI - centerJ)
        commonSum += ((cj * fragJ.getWeight())) / (rij ** 3) * fragJ.getEta()
    commonSum += 1

    fragI.values.corrected   = fragI.values.values.mult(commonSum)

