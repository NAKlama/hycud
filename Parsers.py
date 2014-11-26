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


import re

rangeParser        = re.compile("^(\\d+)-(\\d+)$")
numberParser       = re.compile("^(\\d+)$")
fragmentSpecParser = re.compile("^([xX#]*)(\\d+)(.*)$")
valueExtractor     = re.compile("^([vVhHtT])(\\d+.?\\d*(?:[eE][+-]?\\d+)?)$")
modelParser        = re.compile("^MODEL\\s+(\\d+)\s*$")
endmdlParser       = re.compile("^ENDMDL\\s*$")
atomParser         = re.compile("^ATOM")
pdbResParser       = re.compile("^(ATOM  .{16})(.{4})(.*)[\\n\\r]*$")
pdbParser          = re.compile("^ATOM  (.....).(....).(...)..(....)....(........)(........)(........)(.*)$")
                                #        Atom#  AtomN  ResN    Seq#          X         Y        Z
versionParser      = re.compile("^v(\\d+)\\.(\\d+)\\.(\\d+)$")
bracketParser      = re.compile("^\\(([^\\(\\)]+)\\)(.*)$")
atomExtractor      = re.compile("^\\d*([HCNOPS])\\w*\\d*$")
dataParser         = re.compile("^\\(([vVhHtT].*),([vVhHtT].*)\\)$")

optVerParser			 = re.compile("^\\s*options_ver\\s+=\\s+'(v\\d+\\.\\d+\\.\\d+)'\\s+$")


pattern  = "^\\s*(-?\\d\\.\\d+[eE][-+]\\d+)"
pattern += "\\s+(-?\\d\\.\\d+[eE][-+]\\d+)"
pattern += "\\s+(-?\\d\\.\\d+[eE][-+]\\d+)"
pattern += "\\s+(-?\\d\\.\\d+[eE][-+]\\d+)"
pattern += "\\s+(-?\\d\\.\\d+[eE][-+]\\d+)"
pattern += "\\s+(-?\\d\\.\\d+[eE][-+]\\d+)\\s*$"
diffMatrixParser	 = re.compile(pattern)

