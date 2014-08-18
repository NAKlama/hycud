import os

from HelperFunctions import ANSI_ESC, flush


class HydroProThread:
  def __init__(self, opt, translation):
    self.opt   = opt
    self.trans = translation

  def __call__(self, model):
    command = ""
    if self.trans:
      command = "nice -n {:n} {}".format(self.opt.nice, model.transScript)
    else:
      command = "nice -n {:n} {}".format(self.opt.nice, model.script)
    if self.opt.verbose > 1:
      print(command)
    return os.system(command)


def hydroPro(opt, models):
  """This takes a modelset and executes HydroPro for it (in parallel for HydroPro7"""
  taskDesc      = HydroProThread(opt, False)
  transTaskDesc = HydroProThread(opt, True)

  # Although the capability to start Hydropro in parallel threads still exists,
  # it is not used, since HydroPro 10 will use all cores it can get, making multiple
  # instances of HydroPro ineffective
  done  = 0
  count = models.size()
  for m in models.models:
    if opt.verbose > 0:
      done += 1
      print("{0}2K{0}1GCalculating HydroPro {1:6n}/{2:n}".format(
        ANSI_ESC(), done, count), end='')
      flush()
    if opt.translation:
      transTaskDesc(m)
    if not opt.onlyTrans:
      taskDesc(m)
  if opt.verbose > 0:
    print("")
