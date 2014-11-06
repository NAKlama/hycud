from HelperFunctions import printError


def residueWeight(resName):
  """Function to translate residue-string to molecular weight"""
  # Source: http://web.expasy.org/findmod/findmod_masses.html#AA
  wgt = {"ALA":  71.08, "GLY":  57.05, "ILE": 113.16, "LEU": 113.16, "MET": 131.19
        ,"PRO":  97.12, "VAL":  99.13, "ASN": 114.10, "CYS": 103.14, "GLN": 128.13
        ,"SER":  87.08, "THR": 101.11, "PHE": 147.18, "TRP": 186.21, "TYR": 163.18
        ,"ASP": 115.09, "GLU": 129.12, "ARG": 156.19, "HIS": 137.14, "LYS": 128.17}

  # wgt = {"ALA":  89.09, "GLY":  75.07, "ILE": 131.17, "LEU": 121.17, "MET": 149.21
  #       ,"PRO": 115.13, "VAL": 117.15, "ASN": 132.12, "CYS": 121.16, "GLN": 146.14
  #       ,"SER": 105.09, "THR": 119.12, "PHE": 165.19, "TRP": 204.23, "TYR": 181.19
  #       ,"ASP": 133.10, "GLU": 147.13, "ARG": 174.20, "HIS": 155.15, "LYS": 146.19}
  if resName in wgt:
    return wgt[resName]
  else:
    printError("Residue not recognized: {}".format(resName))


def residueProtons(resName):
  """Function to translate residue-string to non-exchangable proton-count"""
  prot = {"ALA":  4, "GLY":  2, "ILE": 10, "LEU": 10, "MET":  8
         ,"PRO":  7, "VAL":  8, "ASN":  3, "CYS":  3, "GLN":  5
         ,"SER":  3, "THR":  5, "PHE":  8, "TRP":  8, "TYR":  7
         ,"ASP":  3, "GLU":  5, "ARG":  7, "HIS":  5, "LYS":  9}
  if resName in prot:
    return prot[resName]
  else:
    printError("Residue not recognized: {}".format(resName))


class Residue:
  """Simple class to store Residue infomation, automatically looks up weight and proton count"""
  def __init__(self, num, res):
    self.num      = num
    self.res      = res
    self.weight   = residueWeight(res)
    self.protons  = residueProtons(res)
    self.bbNpos   = None
    self.bbCApos  = None
    self.NCvect   = None

  def getNCvect(self):
    if self.NCvect is not None:
      return self.NCvect
    if self.bbNpos is None or self.bbCApos is None:
      return None
    else:
      self.NCvect = self.bbNpos.sub(self.bbCApos)
      return self.NCvect

  def setNpos(self, point):
    if self.bbNpos is not None:
      self.bbNpos = point

  def setCApos(self, point):
    if self.bbCApos is not None:
      self.bbCApos = point


class Protein:
  """Simple Protein(fragment) class, automatically calculates weight and proton count"""
  def __init__(self):
    self.residues = []
    self.weight   = 0.0
    self.protons  = 0

  def addResidue(self, resNum, resName):
    found = False
    for r in self.residues:
      if r.num == resNum:
        found = True
    if not found:
      self.residues.append(Residue(resNum, resName))
      self.weight  += residueWeight(resName)
      self.protons += residueProtons(resName)

  def getWeight(self):
    return self.weight

  def getProtons(self):
    return self.protons

  def done(self):
    pass
    # self.residues = []

  def setNpos(self, res, point):
    for r in self.residues:
      if r.num == res:
        r.setNpos(point)

  def setCApos(self, res, point):
    for r in self.residues:
      if r.num == res:
        r.setCApos(point)

  def getNCvect(self, res):
    for r in self.residues:
      if r.num == res:
        return r.getNCvect()
    return None

