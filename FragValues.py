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


class CalcValues:
  """Class for handling harmonic mean and stokes radius at the same time (math capable)"""
  def __init__(self, r=0.0, hm=0.0):
    self.r    = r
    self.hm   = hm

  def add(self, cv):
    ret    = CalcValues()
    ret.r  = self.r  + cv.r
    ret.hm = self.hm + cv.hm
    return ret

  def addTo(self, cv):
    self.r  += cv.r
    self.hm += cv.hm

  def sub(self, cv):
    return CalcValues(r=(self.r - cv.r), hm=(self.hm - cv.hm))

  def mult(self, v):
    # return CalcValues(r=(self.r * v), hm=(self.hm * v))
    ret     = CalcValues()
    ret.r   = self.r  * v
    ret.hm  = self.hm * v
    return ret

  def div(self, v):
    ret     = CalcValues()
    ret.r   = self.r  / v
    ret.hm  = self.hm / v
    return ret

  def divTo(self, v):
    self.r  /= v
    self.hm /= v

  def sqr(self):
    ret     = CalcValues()
    ret.r   = self.r  ** 2
    ret.hm  = self.hm ** 2
    return ret

  def sqrt(self):
    ret     = CalcValues()
    ret.r   = math.sqrt(self.r)
    ret.hm  = math.sqrt(self.hm)
    return ret

  def recip(self):
    ret       = CalcValues()
    if self.r != 0.0:
      ret.r   = 1.0 / self.r
    else:
      ret.r   = 0.0
    if self.hm != 0.0:
      ret.hm  = 1.0 / self.hm
    else:
      ret.hm  = 0.0
    return ret

  def out(self, prefix=""):
    print("{} r ={:e}".format(prefix, self.r))
    print("{} hm={:e}".format(prefix, self.hm))

  def getR(self):
    return self.r

  def getHM(self):
    return self.hm


class FragValues:
  """Class holding the calculated values for a fragment (math capable)"""
  def __init__(self):
    self.values             = CalcValues()
    self.corrected          = CalcValues()
    self.eta                = 0.0

  def calcValues(self, viscosity=0.0, HarmMe=0.0, radius=0.0):
    if viscosity != 0.0:
      self.eta        = viscosity
    if HarmMe != 0.0:
      self.values.hm  = HarmMe
    if radius != 0.0:
      self.values.r   = radius

  def addTo(self, fv):
    self.values.addTo(fv.values)
    self.corrected.addTo(fv.corrected)
    self.eta        = self.eta + fv.eta

  def subTo(self, fv):
    self.values     = self.values.sub(fv.values)
    self.corrected  = self.corrected.sub(fv.corrected)
    self.eta        = self.eta - fv.eta

  def sub(self, fv):
    ret             = FragValues()
    ret.values      = self.values.sub(fv.values)
    ret.corrected   = self.corrected.sub(fv.corrected)
    ret.eta         = self.eta - fv.eta
    return ret

  def multTo(self, v):
    self.values     = self.values.mult(v)
    self.corrected  = self.corrected.mult(v)
    self.eta        = self.eta * v

  # def multRet(self, v):
  #   ret = FragValues()
  #   ret.values     = self.values.mult(v)
  #   ret.corrected  = self.corrected.mult(v)
  #   ret.eta        = self.eta * v
  #   return ret

  def divTo(self, v):
    self.values.divTo(v)
    self.corrected.divTo(v)
    self.eta        = self.eta / v

  def sqrtTo(self):
    self.values     = self.values.sqrt()
    self.corrected  = self.corrected.sqrt()
    self.eta        = math.sqrt(self.eta)

  def sqr(self):
    ret             = FragValues()
    ret.values      = self.values.sqr()
    ret.corrected   = self.corrected.sqr()
    ret.eta         = self.eta ** 2
    return ret

  def recip(self):
    ret             = FragValues()
    ret.values      = self.values.recip()
    ret.corrected   = self.corrected.recip()
    if self.eta != 0.0:
      ret.eta       = 1.0 / self.eta
    else:
      ret.eta       = 0.0
    return ret

  def getEta(self):
    return self.eta

  def getR(self, corr=False):
    if corr:
      return self.corrected.getR()
    else:
      return self.values.getR()

  def getHM(self, corr=False):
    if corr:
      return self.corrected.getHM()
    else:
      return self.values.getHM()
