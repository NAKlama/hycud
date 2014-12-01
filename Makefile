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


VERSION =$(shell python3 hycud.py --version)
PACKDIR =hycud-${VERSION}
TARFILE =hycud-${VERSION}.tar

.DEFAULT: all
all: clean zip tar

.packDirectory:
	echo ${PACKDIR}
	-mkdir ${PACKDIR}
	-cp *.py ${PACKDIR}
	-cp -R tools ${PACKDIR}
	-cp LICENSE ${PACKDIR}
	-touch .packDirectory

.PHONY: zip
zip: .packDirectory
	-zip -r9 hycud-${VERSION}.zip ${PACKDIR}

tar: .packDirectory
	-tar cf ${TARFILE} ${PACKDIR}
	-gzip -9c ${TARFILE} > ${TARFILE}.gz
	-bzip2 -9c ${TARFILE} > ${TARFILE}.bz2
	-xz -e9 ${TARFILE}

clean:
	-rm -Rf .packDirectory
	-rm -Rf __pycache__
	-rm -Rf hycud-*.*.*
