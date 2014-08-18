import os
import tempfile
import subprocess
import gc

from os import path

from HelperFunctions import waitTillFileExists, ANSI_ESC, flush
from multiprocessing import Pool

esc = ANSI_ESC()

class runREMO:
  """Class returning a function to run REMO for one model"""
  def __init__(self, opt, count):
    self.opt = opt
    self.size = count

  def __call__(self, model):
    tmpDir = ""
    if self.opt.skipREMO:
      # subprocess.check_call(['ln', '-s', model.getPDB(), model.getSprouted()])
      pass
    else:
      if self.opt.keepTemp:
        tmpDir = path.join(self.opt.tmpPath, ("REMO_%09i" % model.num))
        os.makedirs(tmpDir)
      else:
        tmpDir = tempfile.mkdtemp(prefix="REMO_", suffix="", dir=self.opt.tmpPath)
      (PDB_dir, PDB_file) = path.split(model.getPDB())
      tmpPDB  = path.join(tmpDir, PDB_file)
      remoExe = path.join(self.opt.REMOPath, "REMO.pl")
      subprocess.check_call(['ln', '-s', model.getPDB(), tmpPDB])
      waitTillFileExists(tmpPDB)
      os.chdir(tmpDir)
      if self.opt.verbose > 1:
        print("nice -n", str(self.opt.nice), "perl" , remoExe, "0", PDB_file)
      elif self.opt.verbose > 0:
        print("{0}2K{0}1GCalculating REMO {1:6n}/{2:n}".format(
          esc, model.num, self.size), end='')
        flush()
      subprocess.check_output(['nice', '-n', str(self.opt.nice), 'perl', remoExe, "0", PDB_file])
      waitTillFileExists(tmpPDB + ".h")
      subprocess.check_call(['mv', (tmpPDB + '.h'), model.getSprouted()])


def sproutModels(opt, models):
  """Function sprouts full sidechains for a given set of protein models"""
  sprout_pool = Pool(processes=opt.threads)
  task        = runREMO(opt, models.size())
  sprout_pool.map(task, models.models)
  if opt.verbose > 0:
    print("")
  gc.collect()
