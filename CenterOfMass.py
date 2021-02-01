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



import numpy as np
from copy import deepcopy


class CenterOfMass:
  """Class for calculating the center of Mass"""
  def __init__(self, weight=0.0, x=0.0, y=0.0, z=0.0, point=None):
    self.weight = weight
    if point is None:
      self.center = np.array([x, y, z])
    else:
      self.center = deepcopy(point)

  def getCenter(self):
    return self.center / self.weight

  def addPoint(self, weight, point):
    self.center += point * weight
    self.weight += weight

  def printPoint(self):
    p = self.getCenter()
    print("({:g}, {:g}, {:g}) <= ".format(p[0], p[1], p[2]), end='')
    print("({:g}, {:g}, {:g}) / {:g}".format(self.center[0], self.center[1], self.center[2], self.weight))


class IndexedCoM:
  """
  Class for calculating the center of Mass for many points, takes an index
  and calculates individial centers for each center to reduce rounding error
  """
  def __init__(self):
    self.centerList = []
    self.center     = None

  def addPoint(self, index, weight, point):
    found       = None
    self.center = None
    for i in self.centerList:
      if i['index'] == index:
        found = i
    if found is None:
      self.centerList.append({'index':index, 'center':CenterOfMass(weight, point=point)})
    else:
      found['center'].addPoint(weight, point)

  def getCenter(self):
    if self.center is None:
      self.center = CenterOfMass()
      if self.centerList == []:
        return None
      for i in self.centerList:
        c = i['center']
        self.center.addPoint(c.weight, c.getCenter())
      # self.centerList = []
    return self.center.getCenter()

  def printCenter(self):
    p = self.getCenter()
    print("({:g}, {:g}, {:g})  comes from".format(p[0], p[1], p[2]))
    for c in self.centerList:
      print("{:5n}: ".format(c['index']),
            c['center'].printPoint())
