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


from os import path

from HelperFunctions import vers2Num
from Parsers import optVerParser

def updateUserOptions(progVersion):
	userOptFile = path.expanduser("~/.hycud/Options.py")
	if path.exists(userOptFile) and path.isfile(userOptFile):
		optFile = open(userOptFile, "r", encoding="utf-8")
		lines = optFile.readlines()
		optFile.close()
		fileVersion = 0
		for l in lines:
			res = optVerParser.match(l)
			if res:
				versionString = res.group(1)
				fileVersion = vers2Num(versionString)
				break

		for i in range(0, len(lines)):
			if optVerParser.match(lines[i]):
				lines[i] = "options_ver = '{}'\n".format(progVersion)

		if fileVersion < vers2Num('v3.3.0'):
			lines.append('daemon_db_host = "localhost"\n\n')
			lines.append("# DB storage Mode\n# db_mode = 'filesystem'\ndb_mode = 'mysql'\n")

		if fileVersion < vers2Num('v3.4.4'):
			lines.append('\n# In the spectral density function mode, should we try and use partial fragments\n')
			lines.append('# for calculating the correction factor?\n')
			lines.append("# Don't change unless you know what you are doing!\n")
			lines.append('sdf_partial_fragments = True\n')

		if fileVersion < vers2Num('v3.4.6'):
			lines.append('\n# For the spectral density function, the larmor frequency is needed, you can set\n')
			lines.append('# the default here\n')
			lines.append('larmor_freq_proton = 600.25\n')

		with open(userOptFile, "w", encoding='utf-8') as optFile:
			for l in lines:
				optFile.write(l)
