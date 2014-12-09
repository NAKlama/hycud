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
      self.recipSum += 1/val

  def getAvg(self):
    return self.sum / self.weight

  def getHarmMe(self):
    if self.weight != 0:
      x = (1/self.weight) * self.recipSum
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
