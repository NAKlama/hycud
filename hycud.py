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

import sys
from os               import path

# Find customized Options.py in home directory
sys.path.insert(0, path.expanduser("~/.hycud"))

from Options import *
import argparse
import io
import os
# import string
import tempfile
import pickle
import gzip
import numpy as np
from time import sleep

from HelperFunctions  import checkPathExists, printError
from HelperFunctions  import vers2Num, metaConvert
from Model            import Models
from Fragmentation    import Fragmentation
from REMO             import sproutModels
from HydroPro         import hydroPro
from DataDump         import DataDump
from OptionsUpdater   import updateUserOptions
from sdfAnalysis      import sdfAnalysis

default_temporaryStorage    = path.abspath(default_temporaryStorage)
version                     = "v3.4.6"

if vers2Num(options_ver) < vers2Num(version):
  updateUserOptions(version)

class Options:
  """This class holds all options defined on the command line,
  after parsing and some processing"""
  def __init__(self):
    self.verbose    = 0
    self.debug      = 0
    self.resCount   = 0

try:
  default_larmor = larmor_freq_proton
except NameError:
  default_larmor = 600.25



#############################
#                           #
#  M A I N   P R O G R A M  #
#                           #
#############################

if __name__ == '__main__':
  # We use the options class to keep track of all options in one place
  o   = Options()

  # Directory the program was called from
  cwd = os.getcwd()
  desc = "HYCUD {}".format(version)
  # Argument parser
  argParser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawDescriptionHelpFormatter,
      epilog='''\
  MANUAL FRAGMENT SPECIFICATION (--fragments | -F)

    [X][#]<n>[([vVhHtT]<value>,[vVhHtT]<value>)]:[#]<n>[([vVhHtT]<value>},[vVhHtT]<value>})]:...

    A sequence of numbers split by ':'. Without a '#' the numbers are a fragment count, with the '#'
    they are residue numbers. Putting a X before the first number will exclude that fragment from the
    calculation. Any residues that are left at the end will be discarded.
    It is possible to suplly hydrodynamic time and viscosity manually by adding them in brackets behind
    the fragment length definition.

    Example: 7:56(v2.746,h2.99e-9):6


  DETAILED FRAGMENT SPECIFICATION (--detailedFrag | -D)

    (<res>[-<res][ <res>[-<res>][...]]][ [vVhHtT]<value>])...

    A sequence of fragments defined by fragment numbers or number ranges, seperated by spaces, in brackets.
    Optionally hydrodynamic time and viscosity can be manually fixed by adding the value prefixed by 'v'
    or 'V' for viscosity and one of 'h', 'H', 't', 'T' for hydrodynamic time.


  ATTENTION: This program assumes your Residues are numbered incrementally starting from 1!
             If this is not the case, please use --fixResidueNumbering or --detailedFrag
    ''')

  argParser.add_argument('--exe',
      default=default_calcHydr, type=str,
      help="Executable file to use")
  argParser.add_argument('--REMOdir', '-R',
      default=default_REMOdir, type=str,
      help="Path to the REMO directory")
  argParser.add_argument('--NoREMO',
      action='store_true',
      help="Disable the REMO step. (Only for full proteins)")
  argParser.add_argument('--REMOout',
      default="", type=str,
      help="Output PDT after REMO")
  argParser.add_argument('--in', '--template', '-I',
      default=default_templateFile, type=str,
      help="Template file for configuration of executable")
  argParser.add_argument('--tmpDir', '-T',
      default=default_temporaryStorage, type=str,
      help="Directory for temporary files")
  argParser.add_argument('--threads', '-t',
      default=default_threads, type=int,
      help="Number of threads to run REMO with (Default: {:n})".format(default_threads))
  if allow_HydroPro_MultiTreaded:
    argParser.add_argument('--hydroProThreads',
      default=1, type=int,
      help="Number of threads to run HydroPro with    USE WITH CAUTION!")
  argParser.add_argument('--nice',
      default=default_niceness, type=int,
      help="Nice level to use (Default: {:n})".format(default_niceness))
  argParser.add_argument('--translation',
      action='store_true',
      help="Calculate translational diffusion coefficient as well")
  argParser.add_argument('--translationOnly',
      action='store_true',
      help="Calculate translational diffusion coefficient only")
  argParser.add_argument('--fragSize', '-s',
      default=(-1), type=int,
      help="Fragment size (all fragments will be sized equally)")
  argParser.add_argument('--fragments', '-F',
      type=str, default="",
      help="Manually specify fragments. (See MANUAL FRAGMENT SPECIFICATION)")
  argParser.add_argument('--detailedFrag', '-D',
      default="", type=str,
      help="Specify fragments in Detail (See DETAILED FRAGMENT SPECIFICATION)")
  argParser.add_argument('--keepTempFiles', action='store_true',
      help="Don't automatically delete temporary files")
  argParser.add_argument('--verbose', '-v', action='count')
  argParser.add_argument('--pdt',
      type=str, default="",
      help="PDT-File with multiple conformations")
  argParser.add_argument('--convertToNewDataformat', action='store_true',
      help="Converts --inData to new Dataformat and writes it to --outData")
  argParser.add_argument('--outputAdditionalData', action='store_true',
      help="Output additional calculated data")
  argParser.add_argument('--outData',
      default="", type=str,
      help="Name of data output file")
  argParser.add_argument('--inData',
      default="", type=str,
      help="Name of data input file to use for analysis")
  argParser.add_argument('--inCount',
      default="", type=str,
      help="How many results to read from data input format <count>[:offset]")
  argParser.add_argument('--inInfo',
      action='store_true',
      help="Print information about the datafile")
  # argParser.add_argument('--analyseDistribution',
  #     action='store_true',
  #     help="Analyse the distrbutions of Averages")
  # argParser.add_argument('--analyseFragmentDistribution',
  #     action='store_true',
  #     help="Analyse the distrbutions of Averages for each fragment")
  # argParser.add_argument('--distributionStep',
  #     default=0.5, type=float,
  #     help="Stepping for --analyseDistribution")
  argParser.add_argument('--filterOutliers',
      default=0.0, type=float,
      help="Remove Outliers above a cutoff value")
  argParser.add_argument('--outputOutliers',
      action="store_true",
      help="Output model numbers of outliers that were removed")
  argParser.add_argument('--removeModels',
      default="", type=str,
      help="Comma seperated list of models to exclude from analysis")
  argParser.add_argument('--dumpDataTable',
      default="", type=str,
      help="Dump fragment Data into text table defined by argument")
  argParser.add_argument('--displayHarmonicMean',
      action='store_true',
      help="Instead of arithmetric averages display harmonic means")
  argParser.add_argument('--spectralDensityFunction', '--sdf',
      action="store_true",
      help="Calculate spectral density function")
  argParser.add_argument('--sdfResidueSkip',
      default=1, type=int,
      help="Reduce Calculation time by skipping more than one residue in fragmentation")
  argParser.add_argument('--larmorFreq',
      default=default_larmor, type=float,
      help="Set the larmor frequency in MHz (needed for spectral density function)")
  argParser.add_argument('--tau',
      default=1e-10, type=float,
      help="Set tau for the SDF calculation")
  argParser.add_argument('--order',
      default=1.0, type=float,
      help="Set the order parameter for the SDF calculation (1.0)")

  if allow_option_weighted_averages and not default_show_weighted_averages:
    argParser.add_argument('--displayWeightedAverages',
      action='store_true',
      help="Display extra averages, weighted by weight or non exchangeable protons in the summary")
  argParser.add_argument('--version',
      action="store_true",
      help="Display HYCUD version number")

  # Convert arguments into a hash datastructure
  args = vars(argParser.parse_args())

  if args['version']:
    print("{}".format(version[1:]))
    sys.exit(0)

  # Checking if multiple fragmentation specifications where given
  if args['fragSize'] != (-1) and (args['fragments'] != "" or args['detailedFrag']):
    msg = "Please do not specify fragments manually when using the --fragSize option."
    printError(msg)

  if args['fragments'] and args['detailedFrag']:
    msg = "Please only use one of --fragments or --detailedFrag, not both"
    printError(msg)

  o.calcHydr  = path.abspath(args['exe'])
  o.REMOPath  = path.abspath(args['REMOdir'])
  o.runDir    = os.getcwd()

  # Generate temprary directory
  if args['inData'] == "":
    checkPathExists(args['tmpDir'])
    o.tmpPath   = tempfile.mkdtemp(dir=path.abspath(args['tmpDir']), prefix="HYCUD_")
  else:
    o.tmpPath   = cwd

  if not path.isdir(o.tmpPath):
    msg = "{} is not a valid directory".format(o.tmpPath)
    printError(msg)
  if args['fragments'] == "" and args['detailedFrag'] == "" and args['fragSize'] == (-1):
    fragSize = default_fragmentSize
  elif args['fragSize'] > 0:
    fragSize = args['fragSize']
  else:
    fragSize = (-1)

  # Writing the options into the Options clasee
  o.exePath         = path.abspath(args['exe'])
  o.keepTemp        = args['keepTempFiles']
  o.threads         = int(args['threads'])
  o.nice            = args['nice']
  o.skipREMO        = args['NoREMO']
  o.outData         = args['outData']
  o.inData          = args['inData']
  o.convert         = args['convertToNewDataformat']
  o.inCount         = args['inCount']
  o.inInfo          = args['inInfo']
  o.detailFrag      = args['detailedFrag']
  o.REMOout         = args['REMOout']
  o.pdt             = args['pdt']
  # o.distAnaly       = args['analyseDistribution']
  # o.fragDisAn       = args['analyseFragmentDistribution']
  # o.distAnSt        = args['distributionStep']
  o.filtOutl        = args['filterOutliers']
  o.outpOutl        = args['outputOutliers']
  o.translation     = args['translation']
  o.onlyTrans       = args['translationOnly']
  o.outDataTable    = args['dumpDataTable']
  o.verbData        = args['outputAdditionalData']
  o.harmonicMean    = args['displayHarmonicMean']
  o.sdfResidueSkip  = args['sdfResidueSkip']
  o.sdf             = args['spectralDensityFunction']
  o.larmorFreq      = args['larmorFreq']
  o.SDFtau          = args['tau']
  o.SDForder        = args['order']

  if allow_option_weighted_averages and not default_show_weighted_averages:
    o.showWeightedAvg = args['displayWeightedAverages']
  else:
    o.showWeightedAvg = default_show_weighted_averages
  if allow_HydroPro_MultiTreaded:
    o.hydroMT     = args['hydroProThreads']
  else:
    o.hydroMT     = 1

  # If we are only claculation translation, we swich on translation calculation
  if o.onlyTrans and not o.translation:
    o.translation = True

  # Split the comma seperated list into a list of integers
  if args['removeModels'] != "":
    o.remove     = list(map(int, args['removeModels'].split(',')))
  else:
    o.remove     = []

  o.fragSize   = fragSize

  o.verbose     = args['verbose']
  o.keepTemp    = args['keepTempFiles']
  o.sdfPartial  = sdf_partial_fragments

  if not o.verbose:
    o.verbose = 0

  # Checking paths
  if o.inData == "":
    checkPathExists(o.exePath)
    checkPathExists(args['in'])
    if not o.skipREMO:
      checkPathExists(path.join(o.REMOPath, 'REMO.pl'))
    checkPathExists(o.tmpPath)
  else:
    checkPathExists(o.inData)


  if o.inData == "":
    # Hydropro configuration file -> memory
    o.HydroConf       = []
    with io.open(args['in']) as template:
      lines = template.readlines()
      for line in lines:
        o.HydroConf.append(line)

    # Parse PDT and initialize models list
    models = Models()
    models.parsePDT(o, args['pdt'])

    o.fragmentation = Fragmentation(o)
    if o.sdf:
      o.skipREMO = True
      fSize = o.fragSize
      skip  = o.sdfResidueSkip
      start = 1
      end   = models.maxRes
      if skip > 1:
        fSizeDiv =  [ d for d in range(2,fSize//2+1) if fSize % d == 0 ]
        if skip not in fSizeDiv:
          print ( "WARNING: The value selected for skipping ("
                + str(skip)
                + ") is not a divisor of fragment size ("
                + str(fSize)
                + ")" )
          print ( "         Values for different residues will be calculated from "
                + "different number of fragments.")
          print ( "         This means that values are not as consistent over "
                + "different residues.")
          sleep(3)
        if skip > 3:
          print ( "WARNING: Using very high skip values reduces accuracy significantly.")
          sleep(2)
        if (fSize / skip) <= 2:
          print ( "WARNING: With a skip this high, residues will only be found is very "
                + "few fragments.")
          print ( "         Please check your options. If you really want this, please "
                + "wait for 15s.")
          sleep(15)
      o.fragmentation.sdfFrags(start, end, fSize, skip)
    else:
      # Determine fragmentation
      if args['fragments'] != "":
        o.fragmentation.fragDef(args['fragments'])
      elif args['detailedFrag'] != "":
        o.fragmentation.detailedFrag(args['detailedFrag'])
      else:
        o.fragmentation.fragSize(fragSize)

    # Use sprouting tool to make sure we have a full protein
    if not o.skipREMO:
      sproutModels(o,  models)
    if o.verbose > 3:
      print("Residue Count  = {:n}".format(o.resCount))
      print("Fragment Count = {:n}".format(o.fragmentation.fragmentCount()))
    if not o.onlyTrans:
      models.splitFragments(o)
    if o.translation:
      models.translationalDiff(o)

    if o.verbose > 4:
      for i in range(0, o.fragmentation.fragmentCount()):
        print(o.fragmentation.fragmentList(i))

    # HydroPro calculation
    hydroPro(o, models)
    models.parseResults(o)
    if not o.onlyTrans:
      models.calculateWeightingFactors(o)

    # Generate output data structure
    data = DataDump(o.fragmentation, models, version)
    data.addMetaData('args', args)
    data.addMetaData('hydropro.dat', o.HydroConf)
    data.addMetaData('cwd', cwd)
    if o.outData != "":
      with gzip.GzipFile(o.outData, "wb", 9) as outFile:
        outFile.write(pickle.dumps(data, -1))

  else: # Results file
    if o.convert:
      from OldFragmentation import FragmentDef, SubFragment, ProtFragment
      from Model            import Model
      from REMO             import runREMO
      from Fragment         import Fragment
      from CenterOfMass     import *
      from Point3D          import *
      from FragValues       import *
      from Protein          import *
      from Fragmentation    import FragmentDefinition
    with gzip.GzipFile(o.inData, "rb", 9) as inFile:
      data            = pickle.loads(inFile.read())
      o.fragmentation = data.frag
      models          = data.model
      o.onlyTrans     = metaConvert(data.meta)['args']['translationOnly']
      o.translation   = metaConvert(data.meta)['args']['translation']
      o.fragSize      = metaConvert(data.meta)['args']['fragSize']
      if metaConvert(data.meta)['args']['spectralDensityFunction'] and o.fragSize == -1:
        o.fragSize = default_fragmentSize
      if o.onlyTrans and not o.translation:
        o.translation = True
    if vers2Num(metaConvert(data.meta)['multiHydroVersion']) < vers2Num("v2.3.0"):
      print("\nWARNING: The datafile has an old data format. Please use an older version to deccode it!\n")
    if o.inInfo:
      data.printInfo()
    if o.convert and o.outData != "":
      if vers2Num(metaConvert(data.meta)['multiHydroVersion']) < vers2Num("v3.0.0"):
        models = Models()
        models.models = data.model

        oldFragmentation = o.fragmentation.fragments
        o.fragmentation.fragments  = []
        for frag in oldFragmentation:
          definition = ""
          for subFrag in frag.frag.subFragments:
            if definition != "":
              definition += ","
            definition   += "{:n}-{:n}".format(subFrag.start, subFrag.end)
          outFragment = FragmentDefinition(rangesString=definition)
          outFragment.static = frag.static
          outFragment.viscos = frag.viscos
          outFragment.HarmMe = frag.HarmMe
          o.fragmentation.fragments.append(outFragment)

      if vers2Num(metaConvert(data.meta)['multiHydroVersion']) < vers2Num("v3.4.0"):
        for m in models:
          for f in m.fragments:
            for c in f.center.centerList:
              d = c.center
              c.center = np.array([d.x, d.y, d.z])

      dataNew = DataDump(o.fragmentation, models, version)
      for m in data.meta:
        print("{!r}".format(m))
        if not 'multiHydroVersion' in m:
          dataNew.meta.append(m)
      with gzip.GzipFile(o.outData, "wb", 9) as outFile:
        outFile.write(pickle.dumps(dataNew, -1))


  if o.sdf:
    models.calculateWeightingFactors(o)
    sdf = sdfAnalysis(opt, models, o.larmorFreq)
    sdf.calc()
    sdf.average()
    sdf.output()
    sys.exit(0)



  # da = DistributionAnalysis(o, models)
  #
  # if o.distAnaly:
  #   da.generateAnalysis()
  #   da.printAnalysis()
  #   sys.exit()

  # --removeModels
  if len(o.remove) > 0:
    models = models.removeModels(o.remove)

  # --inCount
  if len(o.inCount) > 0:
    o.inCount = o.inCount.split(':')
    o.inCount = list(map(int, o.inCount))
    if(len(o.inCount) > 1):
      models.subModels(o.inCount[0], o.inCount[1])
    else:
      models.subModels(o.inCount[0])

  # --filterOutliers
  if o.filtOutl > 0.0:
    removedModels = []
    if o.filtOutl >= 1.0:
      printError("The number provided should be a ratio, so values ≥ 1 do not make sense.")
    elif o.filtOutl <= 0.0:
      printError("Providing a filtering ratio ≤ 0 will result in no filtering. Please remove the option or provide a sensible filtering ratio.")
    origSize = models.size()
    oldDev = models.getStdDevList(o)
    fragNums = []
    for i in range(0, o.fragmentation.fragmentCount()):
      fragNums.append(i)
    outliers = models.findOutliers()
    if(o.verbose > 2):
      print("Deviation before filtering: {!r}".format(oldDev))
      print("Fragment outliers list: {!r}".format(outliers))
    while len(fragNums) > 0:
      removeList = []
      for n in fragNums:
        removeList.append(outliers[n]['num'])
      newModels = models.removeModels(removeList)
      newModels.calcStats()
      newModels.initFragStats(o)
      newModels.populateFragStats()
      newDev = newModels.getStdDevList(o)
      if(o.verbose > 2):
        print("Deviations after probing removal: {!r}".format(newDev))
      if newModels.size() < origSize / 2:
        print("WARNING: Stopped discarding models, since more than half qualified for being discarded.")
        print("         Please choose a more appropriate minimum percentage change.")
        break
      devCalced = []
      for i in range(0, len(newDev)):
        devFactor = 0
        if oldDev[i] > 0.0 or oldDev[i] < 0.0:
          devFactor = (oldDev[i] - newDev[i])/oldDev[i]
        devCalced.append(devFactor)
        if devFactor < o.filtOutl:
          newFragNums = []
          for fN in fragNums:
            if fN != i:
              newFragNums.append(fN)
          fragNums = newFragNums
      if(o.verbose > 2):
        print("Calculated Deviation Ratios: {!r}:".format(devCalced))
      removeList = []
      for n in fragNums:
        removeList.append(outliers[n]['num'])
      if(o.verbose > 2):
        print("Models to be removed: {!r}".format(removeList))
        print("Fragment count not satisfied before round: {:n}".format(len(fragNums)))

      # Optional code that would only remove if overall improvement exists
      # removeFinal = []
      # testList = []
      # for rM in removeList:
      #   tmpList = removeFinal
      #   tmpList.append(rM)
      #   testModel = models.removeModels(tmpList)
      #   testModel.calcStats()
      #   testModel.initFragStats(o)
      #   testModel.populateFragStats()
      #   testDev = testModel.getStdDevList(o)
      #   score = 0.0
      #   for i in range(0, len(testDev)):
      #     val = oldDev[i]/testDev[i] - 1
      #     score += val
      #   if score >= 0:
      #     removeFinal.append(rM)

      for rM in removeList:
        removedModels.append(rM)

      models = models.removeModels(removeList)
      oldDev = newDev

    # --outputOutliers
    if o.outpOutl:
      print("Removed Models: {}".format(removedModels))

  if o.verbose > 0:
    print("Number of models:        {:n}".format(models.size()))

  # Calculate statistics
  if not o.onlyTrans:
    models.initFragStats(o)
    models.populateFragStats()
  if o.translation:
    models.transStatisticss()
    models.outputTransFragResults(o)

  # if o.fragDisAn:
  #   da.generateAnalysis()
  #   da.fragmentAnalysis()
  #   da.printFragmentwiseAnalysis()
  #   sys.exit()

  # --dumpDataTable
  if o.outDataTable != "":
    with io.open(o.outDataTable, "w", encoding='utf-8') as outTable:
      if not o.onlyTrans:
        outLine  = "                    "
        outLine += "                 Center                   "
        outLine += "\n"
        outTable.write(outLine)
        outLine  = " ModelNo    FragNo  "
        outLine += "           X             Y             Z  "
        outLine += "       eta  "
        outLine += "      r-corr       hm-corr  "
        outLine += "           r            hm  "
        outLine += "\n"
        outTable.write(outLine)
        for m in models.models:
          for f in m.fragments:
            outLine  = "{:8n}  {:8n}  ".format(m.num, f.num)
            center   = f.center.getCenter()
            outLine += "{:12f}  {:12f}  {:12f}  ".format(center[0], center[1], center[2])
            outLine += "{:10f}  ".format(f.values.eta)
            outLine += "{:12e}  {:12e}  ".format(f.values.corrected.r, f.values.corrected.hm)
            outLine += "{:12e}  {:12e}  ".format(f.values.values.r, f.values.values.hm)
            outLine += "\n"
            outTable.write(outLine)
      if o.translation and not o.onlyTrans:
        outTable.write("\n")
        outTable.write("\n")
      if o.translation:
        outLine = "{:>8s}  {:s}\n".format("ModelNo", "Trans. Diff. Coeff")
        outTable.write(outLine)
        for m in models.models:
          outLine  = "{:8n}  ".format(m.num)
          outLine += "{:12e}\n".format(m.transDiff)
          outTable.write(outLine)

  if o.onlyTrans:
    sys.exit(0)
  models.outputFragResults(o)
  models.calcStats()
  models.output(o)
  if o.translation:
    models.outputTransFragResults()

  if not o.keepTemp and o.inData == "":
    os.system("rm -Rf {}".format(o.tmpPath))
