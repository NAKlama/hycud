# HYCUD
# Copyright (C) 2014 Klama, Nina Alexandra and Rezaei-Ghaleh, Nasrollah
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


import re

expParser = re.compile("^\\s*(-?\\d\\.\\d+)([Ee][+-]\\d+)\\s*$")
csi = chr(27) + '['

def getColor(ratio):
  if ratio >= 0.99:
    return csi + "32m"
  elif ratio >= 0.8:
    return csi + "36m"
  elif ratio >= 0.5:
    return csi + "33m"
  elif ratio >= 0.25:
    return csi + "35m"
  else:
    return csi + "31m"


def resetColor():
  return csi + '0m'

def getTerminalSize():
  import os
  env = os.environ
  def ioctl_GWINSZ(fd):
    try:
      import fcntl, termios, struct
      cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,
        '1234'))
    except:
      return
    return cr
  cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
  if not cr:
    try:
      fd = os.open(os.ctermid(), os.O_RDONLY)
      cr = ioctl_GWINSZ(fd)
      os.close(fd)
    except:
      pass
  if not cr:
    cr = (env.get('LINES', 25), env.get('COLUMNS', 80))
  return int(cr[1]), int(cr[0])

def coloredExp(num, ratio, formatStr="{:>11.3e}"):
  inStr = formatStr.format(num)
  res = expParser.match(inStr)
  if res:
    return "{1:}{2:}{3:}{0:}{2:}{4:}{0:}".format(
      resetColor(),
      csi + "1m",
      getColor(ratio),
      res.group(1),
      res.group(2))
  else:
    return inStr
