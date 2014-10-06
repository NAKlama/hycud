#!/usr/bin/python3

from os import path

from HelperFunctions import vers2Num
from Parsers import optVerParser

def updateUserOptions():
	if path.exists("~/.hycud/Options.py") and path.isfile("~/.hycud/Options.py"):
		optFile = open("~/.hycud/Options.py", "r", encoding="utf-8")
		lines = optFile.readlines()
		optFile.close()
		version = 0
		for l in lines:
			res = optVerParser.match(l)
			if res:
				versionString = res.group(1)
				version = vers2Num(versionString)
				break

		optFile = open("~/.hycud/Options.py", "a", encoding='utf-8')

		if version < vers2Num('v3.3.0'):
			optFile.write('daemon_db_host = "localhost"\n\n')
			optFile.write("# DB storage Mode\n# db_mode = 'filesystem'\ndb_mode = 'mysql'\n")

		optFile.close()
