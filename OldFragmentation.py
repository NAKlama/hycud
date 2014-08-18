from Parsers import dataParser, valueExtractor, rangeParser


class FragmentDef:
  """Class for fragment definition from options"""
  def __init__(self, frag, values):
    self.frag   = frag
    self.static = False
    self.viscos = 0.0
    self.HarmMe = 0.0
    if values != "":
      dMatch = dataParser.match(values)
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
        parts = values.split(' ')
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

  def isInFragment(self, n):
    return self.frag.isInFragment(n)

  def resCount(self):
    return self.frag.resCount()

  def resPrint(self):
    return self.frag.resPrint()


class SubFragment:
  def __init__(self, start, end):
    if start <= end:
      self.start = start
      self.end   = end
    else:
      self.end   = start
      self.start = end

  def isInFragment(self, n):
    if n >= self.start and n <= self.end:
      return True
    return False

  def resCount(self):
    return self.end - self.start + 1

class ProtFragment:
  def __init__(self, rangesString=None, begin=0, count=0):
    self.subFragments = []
    if rangesString is not None:
      ranges = rangesString.split(' ')
      for r in ranges:
        match = rangeParser.match(r)
        if match:
          self.subFragments.append(SubFragment(int(match.group(1)), int(match.group(2))))
        else:
          match = numberParser.match(r)
          if match:
            self.subFragments.append(SubFragment(int(match.group(1)), int(match.group(1))))
    else:
      self.subFragments.append(SubFragment(begin, begin+count))

  def isInFragment(self, n):
    for f in self.subFragments:
      if f.isInFragment(n):
        return True
    return False

  def resCount(self):
    count = 0
    for f in self.subFragments:
      count += f.resCount()
    return count

  def resPrint(self):
    out = ""
    for f in self.subFragments:
      out += "%s-%s " % (f.start, f.end)
    return out[:-1]
