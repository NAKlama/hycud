#!/usr/bin/python3

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

		with open(userOptFile, "w", encoding='utf-8') as optFile:
			for l in lines:
				optFile.write(l)
