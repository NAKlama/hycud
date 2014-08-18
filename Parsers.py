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
