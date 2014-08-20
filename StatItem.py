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
