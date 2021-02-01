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
import numpy as np
from copy import *

class StatItem:
  """Class for calculation of weighted average and weighted standard deviation for a set of values"""
  def __init__(self):
    self.sum        = 0.0
    self.recipSum   = 0.0
    self.stdDev     = None
    self.values     = []
    self.weight     = 0.0

  def addValue(self, val, weight=1):
    self.values.append(val * weight)
    self.sum      += val * weight
    self.weight   += weight
    if val != 0.0:
      self.recipSum += weight/val

  def getAvg(self):
    return self.sum / self.weight

  def getHarmMe(self):
    if self.weight != 0:
      x = self.recipSum / self.weight
      if x != 0.0:
        return 1/x

    return 0.0

  def getStdDev(self):
    if self.stdDev is None:
      self.stdDev = 0.0
      for v in self.values:
        self.stdDev += (v - (self.sum / self.weight)) ** 2
      self.stdDev /= self.weight
      self.stdDev  = math.sqrt(self.stdDev)
      self.values  = []
    return self.stdDev



class npStatItem:
  def __init__(self):
    self.sum        = None
    self.recipSum   = None
    self.stdDev     = None
    self.shape      = None
    self.values     = []
    self.weight     = 0.0

  def addValue(self, val, weight=1):
    v = np.array(val)
    self.stdDev = None
    self.values.append(v * weight)
    self.weight   += weight
    if self.sum is None:
      self.shape      = deepcopy(v.shape)
      self.recipSum   = np.zeros(self.shape)
      self.sum        = np.zeros(self.shape)
      self.sum       += v * weight
      nonZero         = np.nonzero(v)
      for i in np.nditer(nonZero):
        self.recipSum[i] += weight / v[i]
    elif val.shape == self.shape:
      self.sum       += v * weight
      nonZero         = np.nonzero(v)
      for i in np.nditer(nonZero):
        self.recipSum[i] += weight / v[i]
    else:
      print("Shapes don't match (own:", self.shape, " in:", val.shape, ")")

  def addStatItem(self, inD):
    self.stdDev = None
    if inD is not None:
      if inD.shape == self.shape:
        # print("Input viable and shapes match")
        if self.sum is None:
          self.sum        = np.zeros(inD.shape)
          self.recipSum   = np.zeros(inD.shape)
        self.sum        += inD.sum
        self.recipSum   += inD.recipSum
        for v in inD.values:
          self.values.append(deepcopy(v))
        self.weight     += inD.weight
      elif self.shape is None:
        self.shape    = deepcopy(inD.shape)
        self.sum      = deepcopy(inD.sum)
        self.recipSum = deepcopy(inD.recipSum)
        self.values   = deepcopy(inD.values)
        self.weight   = deepcopy(inD.weight)

  def getAvg(self):
    if self.sum is not None:
      return (self.sum / self.weight).T
    return None

  def getHarmMe(self):
    if self.weight != 0 and self.sum is not None:
      x       = self.recipSum / self.weight
      nonZero = np.nonzero(x)
      for i in np.nditer(nonZero):
        x[i] = 1 / x[i]
      return x
    return np.zeros(self.shape)

  def output(self, text=""):
    print("\n{} npStatItem:".format(text))
    print("  Sum:         ", self.sum)
    print("  recipSum:    ", self.recipSum)
    print("  stdDev:      ", self.stdDev)
    print("  shape:       ", self.shape)
    print("  values.count:", len(self.values))
    print("  weight:      ", self.weight)

  def getStdDev(self):
    if self.weight is None or self.weight == 0.0:
      return 0.0
    if self.stdDev is None:
      # if (self.sum / self.weight).T[0][0] == 0.0:
        # self.output()
        # for v in self.values:
        #   print(v)
      self.stdDev = np.zeros(self.shape)
      for v in self.values:
        self.stdDev += (v - (self.sum / self.weight)) ** 2
      self.stdDev /= self.weight
      self.stdDev  = np.sqrt(self.stdDev)
      # self.values  = []
    return self.stdDev.T
