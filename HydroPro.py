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

  done  = 0
  count = models.size()
  if opt.hydroMT > 1:
    from multprocessing import Pool
    hp_pool = Pool(processes=opt.hydroMT)
    if opt.translation:
      hp_pool.map(taskDesc, models.models)
    if not opt.onlyTrans:
      hp_pool.map(transTaskDesc, models.models)
  else:
    for m in models.models:
      if opt.verbose > 0:
        done += 1
        # Print progress, while staying on one line, using ANSI ESC Sequences
        print("{0}2K{0}1GCalculating HydroPro {1:6n}/{2:n}".format(
          ANSI_ESC(), done, count), end='')
        flush()
      if opt.translation:
        transTaskDesc(m)
      if not opt.onlyTrans:
        taskDesc(m)
    if opt.verbose > 0:
      print("")
