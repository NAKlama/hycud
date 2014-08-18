import sys

from HelperFunctions import metaConvert


class DataDump:
  """Class defining the data in the DataFile"""
  def __init__(self, fragmentation, model, version):
    self.frag  = fragmentation
    self.model = model
    self.meta  = [{'multiHydroVersion': version}]

  def addMetaData(self, key, value):
    self.meta.append({key: value})

  def getMetaData(self, key):
    meta = metaConvert(self.meta)
    if key in meta:
      return meta[key]
    else:
      return None

  def printInfo(self):
    argumentBlacklist = [ "convertToNewDataformat"
                        , "inData"]
    m = metaConvert(self.meta)
    print("Data version:            {}".format(m['multiHydroVersion']))
    print("Number of models:        {:n}".format(self.model.size()))

    if 'cwd' in m:
      print("Executed from directory: {}".format(m['cwd']))

    if 'args' in m:
      print("\nArguments:")
      for a in m['args'].items():
        if not a[0] in argumentBlacklist:
          print("    {:<30s} : {}".format(a[0], a[1]))

    if 'hydropro.dat' in m:
      print("\nHydroPro configuration used:")
      for line in m['hydropro.dat']:
        print("    ",line,sep='',end='')

    sys.exit()

  def getVerbose(self):
    m = metaConvert(self.meta)
    if 'args' in m:
      for a in m['args'].items():
        if a[0] == 'verbose':
          return a[1]
    return 0
