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
