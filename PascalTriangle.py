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


# This function will only output a pascal triagle, with order n
class PascalTriangle:
	def __init__(self, n):
		self.n = n
		self.t = [1,1]
		if n < 0:
			self.t = None
			return
		if n == 0:
			self.t = [1]
			return
		if n == 1:
			return
		for x in range(n-1):
			tOut = [1]
			for i in range(len(self.t)-1):
				tOut.append(self.t[i] + self.t[i+1])
			tOut.append(1)
			self.t = tOut

	def getList(self):
		return self.t

	def getItem(self, n):
		if n < len(self.t):
			return self.t[n]
		return None

	def getCenter(self):
		return (self.n+2.0 / 2.0)

