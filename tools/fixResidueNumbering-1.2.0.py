#!/usr/bin/python3

import argparse
import re
# import string
from os import path

version = "v1.2.0"
pdbParser          = re.compile("^ATOM  (.....).(....).(...).(.)(....)(....)(........)(........)(........)(.*)$")
terParser          = re.compile("^TER   (.....)......(...).(.)(....)(.*)$")
modelParser        = re.compile("^MODEL\\s+(\\d+)\s*$")


class ArgumentError:
  pass


def checkPathExists(val):
  if not path.exists(val):
    raise ArgumentError("File does not exist: {}".format(val))

if __name__ == '__main__':
  argParser = argparse.ArgumentParser(description=("fixResidueNumbering.py {}".format(version)))
  argParser.add_argument('--removeDNA', action='store_true',
    help="Don't include DNA in output file")
  argParser.add_argument('in', type=str,
    help="Input file")
  argParser.add_argument('out', type=str,
    help="Output file")

  args = vars(argParser.parse_args())

  checkPathExists(args['in'])

  pdbLines = []
  outLines = []

  with open(args['in'], 'r') as pdbFile:
    pdbLines = pdbFile.readlines()

  outRes = 0
  lastRes = -999999
  for l in pdbLines:
    line = l.strip()
    atomMatch   = pdbParser.match(line)
    terMatch    = terParser.match(line)
    modelMatch  = modelParser.match(line)
    if atomMatch:
      atomNum   = int(atomMatch.group(1))
      atomName  = atomMatch.group(2).strip()
      resName   = atomMatch.group(3).strip()
      misc1     = atomMatch.group(4)
      resNum    = int(atomMatch.group(5))
      misc2     = atomMatch.group(6)
      x         = float(atomMatch.group(7))
      y         = float(atomMatch.group(8))
      z         = float(atomMatch.group(9))
      rest      = atomMatch.group(10)

      if resNum != lastRes:
        outRes += 1
        lastRes = resNum

      if atomName in ["CA", "N", "C", "O"]:
        if atomName != "CA":
          atomName += " "
        atomName += " "

      if not args['removeDNA'] or not resName in ["DA", "DG", "DT", "DC"]:
        outLines.append("ATOM  {:5d} {:4s} {:3s} {:1s}{:4d}{:4s}{:8.3f}{:8.3f}{:8.3f}{:s}\n".format(
          atomNum, atomName, resName, misc1, outRes, misc2, x, y, z, rest))

    elif terMatch:
      atomNum = int(terMatch.group(1))
      resName = terMatch.group(2).strip()
      misc    = terMatch.group(3)
      resNum  = int(terMatch.group(4))
      rest    = terMatch.group(5)

      if resNum != lastRes:
        outRes += 1
        lastRes = resNum

      outLines.append("TER   {:5d}      {:3s} {:1s}{:4d}{:s}\n".format(
        atomNum, resName, misc, outRes, rest))

    elif modelMatch:
      outRes = 0
      lastRes = -999999

      outLines.append(line + "\n")

    else:
      outLines.append(line + "\n")

  with open(args['out'], 'w') as outFile:
    outFile.writelines(outLines)
