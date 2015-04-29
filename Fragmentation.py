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


import string
import sys

from HelperFunctions  import printWarning, printError, splitBrackets
from Parsers          import rangeParser, numberParser, dataParser
from Parsers          import fragmentSpecParser, valueExtractor


class FragmentDefinition:
  """This class describes a single fragment"""
  def __init__( self, rangesString=None,
                begin=0, count=0,
                valuesSetStr=None,
                partial=False):
    self.static   = False
    self.viscos   = 0.0
    self.HarmMe   = 0.0
    self.partial  = partial
    self.begin    = begin
    if valuesSetStr is not None:
      dMatch = dataParser.match(valuesSetStr)
      if dMatch:
        gotViscos = False
        gotHarmMe = False
        val = [dMatch.group(1), dMatch.group(2)]
        for v in val:
          vMatch = valueExtractor.match(v)
          if vMatch:
            vType = vMatch.group(1)
            value = vMatch.group(2)
            if vType == 'v' or vType == 'V':
              gotViscos = True
              self.viscos = float(value)
            elif vType == 't' or vType == 'T' or vType == 'h' or vType == 'H':
              gotHarmMe = True
              self.HarmMe = float(value)
        if gotViscos and gotHarmMe:
          self.static = True
      else:
        parts = valuesSetStr.split(' ')
        gotViscos = False
        gotHarmMe = False
        for p in parts:
          vMatch = valueExtractor.match(p)
          if vMatch:
            vType = vMatch.group(1)
            value = vMatch.group(2)
            if vType == 'v' or vType == 'V':
              gotViscos = True
              self.viscos = float(value)
            elif vType == 't' or vType == 'T' or vType == 'h' or vType == 'H':
              gotHarmMe = True
              self.HarmMe = float(value)
        if gotViscos and gotHarmMe:
          self.static = True

    self.resNumbers = []
    if rangesString is not None:
      ranges = rangesString.split(' ')
      for r in ranges:
        match = rangeParser.match(r)
        if match:
          begin = int(match.group(1))
          end   = int(match.group(2))
          if end < begin:
            tmp   = end
            end   = begin
            begin = tmp
          for i in range(begin, end+1):
            self.resNumbers.append(i)
        else:
          match = numberParser.match(r)
          if match:
            self.resNumbers.append(int(match.group(1)))
    else:
      for i in range(begin, begin+count):
        self.resNumbers.append(i)

  def isInFragment(self, n):
    if n in self.resNumbers:
      return True
    return False

  def resCount(self):
    return len(self.resNumbers)

  def resPrint(self):
    out = ""
    prev = None
    start = None
    for rn in self.resNumbers:
      if prev is None:
        prev = rn
        start = rn
      elif prev + 1 == rn:
        prev = rn
      elif prev + 1 != rn:
        if start == prev:
          out += "{:n} ".format(start)
          start = prev
          prev  = rn
        else:
          out += "{:n}-{:n} ".format(start, prev)
          start = rn
          prev  = rn
    if start == prev:
      out += "{:n}".format(start)
    else:
      out += "{:n}-{:n}".format(start, prev)
    return out


class Fragmentation:
  """Class that defines the fragmentation of the protein"""
  def __init__(self, opt):
    """Initialize values from options"""
    self.verbose        = opt.verbose
    self.resCount       = opt.resCount
    self.firstFragment  = 0
    self.fragments      = []

  def fragSize(self, fragSize):
    """Fragmentation by evenly fragmenting the protein into fixed size fragments"""
    begin = 1
    for i in range(fragSize, self.resCount, fragSize):
      self.fragments.append(FragmentDefinition(begin=begin, count=fragSize))
      begin += fragSize
    rest = self.resCount % fragSize
    if (float(rest) / float(fragSize)) > 0.5:
      self.fragments.append(FragmentDefinition(begin=begin, count=rest))
    elif self.verbose > 0 and rest > 0:
      printWarning("Discarding {:n} residues at the end of the protein.".format(rest))

  def fragDef(self, fragmentString):
    """Fragment the protein according to a set of fragment sizes given as a parameter"""
    fragmentSpec = fragmentString.split(':')
    currRes      = 0
    for spec in fragmentSpec:
      specMatch = fragmentSpecParser.match(spec)
      if specMatch:
        options   = specMatch.group(1)
        num       = int(specMatch.group(2))
        values    = specMatch.group(3)
        if options.find('x') != -1 or options.find('X') != -1:
          if currRes != 0:
            printError("The option 'X' is only allowed for the first element.")
          else:
            self.firstFrament = num
            if options.find('#') == -1:
              self.firstFragment += 1
        if options.find('#') != -1:
          if num <= currRes:
            printError("Residue numbers have to be strictly in increasing order!")
          num -= currRes
        if currRes < self.resCount:
          if currRes + num <= self.resCount:
            self.fragments.append(
              FragmentDefinition(
                begin=(currRes+1), count=(currRes+num),
                valuesSetStr=values))
            currRes += num
          else:
            self.fragments.append(
              FragmentDefinition(
                begin=(currRes+1), count=(self.resCount - currRes),
                valuesSetStr=values))
      else:
        print("ERROR: Could not match fragment specification: '{}'".format(spec))
        sys.exit(1)

  def detailedFrag(self, detailedFragString):
    """Fragmentation controlled into the last detail, by a more complicated definition string"""
    ranges = []
    splitBrackets(detailedFragString, ranges)
    if self.verbose > 3:
      print("{}".format(ranges))
    for r in ranges:
      parts = r.split(' ')
      if self.verbose > 3:
        print("{}".format(parts))
      rangeOut = ""
      values   = ""
      for p in parts:
        if rangeParser.match(p):
          rangeOut += p + " "
        elif numberParser.match(p):
          rangeOut += p + " "
        elif valueExtractor.match(p):
          values += p + " "
      if self.verbose > 3:
        print("{}\n{}".format(rangeOut, values))
      self.fragments.append(FragmentDefinition(
        rangesString=rangeOut[:-1],
        valuesSetStr=values))

  def sdfFrags(self, start, end, size, skip):
    """Fragmentation for the spectral density function calculation"""
    for i in range(start, end - size + 2, skip):
      self.fragments.append(FragmentDefinition(begin=i, count=size))

    for i in range(0, size):
      if i < size / 2:
        # print("Fstart: {}  {}".format(start, size-i-1))
        self.fragments.append(FragmentDefinition(
          begin=start,
          count=size - i - 1 ,
          partial=True))

    for i in range(0, size):
      if i < size / 2:
        begin = end - ((end+1) % size) + i
        # print("Fend:   {}  {}".format(begin, end-begin))
        self.fragments.append(FragmentDefinition(
          begin=begin,
          count=end - begin,
          partial=True))

  def fragmentCount(self):
    return len(self.fragments)

  def partOfFragment(self, fragNum, resNum):
    return self.fragments[fragNum].isInFragment(resNum)

  def fragmentList(self, fragNum):
    return self.fragments[fragNum].resNumbers
