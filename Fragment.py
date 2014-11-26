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


from os import path

from FragValues       import FragValues
from Protein          import Protein
from CenterOfMass     import IndexedCoM
from HelperFunctions  import atomWeight


class FragStatistics:
  """Class for keeping the statistics for a fragment"""
  def __init__(self):
    self.n        = 0
    self.avg      = FragValues()
    self.stdDev   = FragValues()
    self.harmMean = FragValues()
    self.residues = 0


class Fragment:
  """Class describing fragments"""
  def __init__(self, num, modelPath, partial):
    self.num          = num
    self.basename     = "Frag{:04n}".format(num)
    self.basepath     = path.join(modelPath, self.basename)
    self.values       = FragValues()
    self.protein      = Protein()
    self.center       = IndexedCoM()
    self.resCenters   = []
    self.stat         = FragStatistics()
    self.atomCount    = 0
    self.residues     = []
    self.diffMat      = []
    self.partial      = partial

  def calcValues(self, viscosity=0.0, HarmMe=0.0, radius=0.0):
    self.values.calcValues(viscosity=viscosity, HarmMe=HarmMe, radius=radius)

  def addAtom(self, atomName, resNum, resName, point):
    weight = atomWeight(atomName)
    self.center.addPoint(resNum, weight, point=point)
    self.protein.addResidue(resNum, resName)
    self.atomCount   += 1
    if resNum not in self.residues:
      self.residues.append(resNum)
    if resNum > 1:
      if atomName.strip() == "N":
        self.protein.setNpos(resNum, point)
      if atomName.strip() == "H":
        self.protein.setHpos(resNum, point)

  def getWeight(self):
    return self.protein.getWeight()

  def getProtons(self):
    return self.protein.getProtons()

  def hasResidue(self, num):
    if num in self.residues:
      return True
    return False

  def getCenter(self):
    return self.center.getCenter()

  def getEta(self):
    return self.values.getEta()

  def getR(self, corr):
    return self.values.getR(corr)

  def getHM(self, corr):
    return self.values.getHM(corr)

  def getPDB(self):
    return self.basepath + '.pdb'

  def getDat(self):
    return self.basepath + '.dat'

  def doneParsing(self):  # This function simply exists to free memory
    self.protein.done()


class FragmentStatistics:
  def __init__(self, opt):
    self.stats              = []
    self.fragmentation      = opt.fragmentation
    for frag in self.fragmentation.fragments:
      stat = FragStatistics()
      stat.residues = frag.resCount()
      stat.resPrint = frag.resPrint()
      self.stats.append(stat)

  def populate(self, models):
    for m in models.models:
      for frag in m.fragments:
        self.stats[frag.num].avg.addTo(frag.values)
        self.stats[frag.num].n     += 1
        self.stats[frag.num].harmMean.addTo(frag.values.recip())

    for s in self.stats:
      s.avg.divTo(s.n)
      s.harmMean.multTo(1.0 / s.n)
      s.harmMean = s.harmMean.recip()

    for m in models.models:
      for frag in m.fragments:
        self.stats[frag.num].stdDev.addTo(frag.values.sub(self.stats[frag.num].avg).sqr())

    for s in self.stats:
      s.stdDev.divTo(s.n)
      s.stdDev.sqrtTo()
