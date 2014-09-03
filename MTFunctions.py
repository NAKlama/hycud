#!/usr/bin/python3

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
        rij = centerI.dist(centerJ)
        commonSum += ((cj * fragJ.getWeight())) / (rij ** 3) * fragJ.getEta()
    commonSum += 1

    fragI.values.corrected   = fragI.values.values.mult(commonSum)

