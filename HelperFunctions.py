import time
import sys

from os import path

from Parsers import pdbParser, versionParser, bracketParser, atomExtractor


def ANSI_ESC():
  return "\x1b["


def flush():
  sys.stdout.flush()


def waitTillFileExists(filename):
  while not path.exists(filename):
    time.sleep(0.1)


def vers2Num(versString):
  vMatch = versionParser.match(versString)
  out = 0
  if vMatch:
    out += int(vMatch.group(1)) * 1000000
    out += int(vMatch.group(2)) * 1000
    out += int(vMatch.group(3))
  return out


def divisors(num):
  divisors = [1]
  for i in range(2, (num - 1)):
    if num % i == 0:
      divisors.append(i)
  return divisors


def splitBrackets(inString, ret=[]):
  match = bracketParser.match(inString)
  ret.append(match.group(1))
  rest = match.group(2)
  if len(rest) > 0 and rest[0] == '(':
    splitBrackets(rest, ret)
  else:
    return ret


def printError(message, errNo=1):
  print("ERROR: {}".format(message))
  sys.exit(errNo)


def printWarning(message):
  print("WARNING: {}".format(message))


def checkPathExists(val):
  if not path.exists(val):
    printError("File does not exist: {}".format(val))


def outSumm(a, b, c, d, e, f):
  """Helper function to shorten lines and better readability in the *Summary* class"""
  print("{} {} {} = {:e} {} (+/- {:e})".format(a, b, c, d, e, f))


# Functions that return atomic weights and proton counts

def atomWeight(atomName):
  """Function returns atomic weight given an atom-name"""
  res = atomExtractor.match(atomName)
  if res:
    atomSymbol = res.group(1)
    if atomSymbol == "H": return  1.00794
    if atomSymbol == "C": return 12.0107
    if atomSymbol == "N": return 14.0067
    if atomSymbol == "O": return 15.9994
    if atomSymbol == "P": return 30.973762
    if atomSymbol == "S": return 32.065
  printWarning("Defaulting to Carbon for unknown atom {}".format(atomName))
  return 12.0107


def removeChainIDs(inLine):
  """This function takes an ATOM line from a PDB file and replaces the chainID by a space"""
  match = pdbParser.match(inLine)
  if match:
    outLine = "ATOM  {} {} {}  {}    {}{}{}{}\r\n".format(
        match.group(1), match.group(2), match.group(3), match.group(4),
        match.group(5), match.group(6), match.group(7), match.group(8))
  else:
    outLine = inLine
  return outLine

def metaConvert(inMeta):
  l = []
  for m in inMeta:
    l.append(list(m.items())[0])
  return dict(l)
