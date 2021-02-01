#!/usr/bin/python3

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



import argparse
import sys
import os
import re
from os import path

version = "v1.0.1"
dotPDB  = re.compile("^.*\\.pdb$")
modelP  = re.compile("^MODEL")
eModP   = re.compile("^ENDMDL")


if __name__ == '__main__':
	argParser = argparse.ArgumentParser(
		description="makePDT.py {}".format(version))
	argParser.add_argument('--inDir', type=str, default=".",
		help="Directory with input PDT files [.]")
	argParser.add_argument('--outFile', type=str, default="",
		help="Output into file instead of STDOUT")

	args = vars(argParser.parse_args())

	inDir       = args['inDir']
	outFileName = args['outFile']

	if not path.exists(inDir):
		print("The path '{}' does not exist!".format(inDir))
		sys.exit(1)
	else:
		inDir = path.abspath(inDir)

	if outFileName != "":
		(Path, File) = path.split(outFileName)
		if Path != "" and not path.exists(Path):
			print("The directory '{}' for the output file does not exist!".format(Path))
			sys.exit(1)

	outFile = sys.stdout
	if outFileName != "":
		outFile = open(outFileName, "w", encoding="utf-8")

	model = 1

	outFile.write("REMARK Generated with makePDT.py ({}), written by Nina Alexandra Klama\n".format(version))
	outFile.write("REMARK \n")
	outFile.write("REMARK Input directory: {}\n".format(inDir))

	dirList = os.listdir(inDir)
	for f in dirList:
		fileName = path.join(inDir, f)
		if dotPDB.match(fileName):
			outFile.write("MODEL {:n}\n".format(model))
			outFile.write("REMARK Source filename for model: {}\n".format(f))
			model += 1

			with open(fileName, "r", encoding="utf-8") as inPDB:
				pdbLines = inPDB.readlines()
				for l in pdbLines:
					if not modelP.match(l) and not eModP.match(l):
						outFile.write(l)

			outFile.write("ENDMDL\n")



