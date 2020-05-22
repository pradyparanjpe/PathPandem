#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# Copyright 2020 Pradyumna Paranjape
# This file is part of PathPandem.
#
# PathPandem is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PathPandem is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PathPandem.  If not, see <https://www.gnu.org/licenses/>.
'''Class'''

from numpy.random import randint


class person(object):
    '''Person class tracking each Person'''
    def __init__(
            self, parent=None, active: bool=False, recovered: bool=False,
            susceptible: float=None, support = False, health: float=None,
            comorbidity: float=0., progress: float=0., move_per_day: int=0,
            strain: int=None, home: tuple=(0, 0), p_max: int=0,
            rms_v: float=0)-> None:
        '''Initialize a (Null) Person'''
        self.p_max = p_max  # Maximum walking reach (y)
        self.move_per_day = move_per_day  # Random walk edge length
        self.rms_v = rms_v * 1.4142  # Number of random walk vertices per day
        if parent:  # Copy Attributes
            self.p_max = self.p_max or parent.p_max
        self.active: bool = active  # Active Infection
        self.recovered: bool = recovered  # Recovered from Infection
        self.progress: float = progress  # Progress of active infection
        self.strain: int = strain  # Pathogen Strain
        self.comorbidity: float = comorbidity  # Predisposed complications
        self.support: bool = support  # On life support
        if not(self.p_max):
            self.home = (0, 0)
        else:
            # everyone's init position is uniformly randomly guessed
            self.home = (randint(self.p_max), randint(self.p_max))
        if parent:  # To copy attribures, NOT the biological reproduction
            self.active = self.active or parent.active
            self.recovered = self.recovered or parent.recovered
            self.progress = self.progress or parent.progress
            self.move_per_day = self.move_per_day or parent.move_per_day
            self.strain = self.strain or parent.strain
            self.rms_v = self.rms_v or parent.rms_v
            self.p_max = self.p_max or parent.p_max
            self.comorbidity = self.comorbidity or parent.comorbidity
            self.support = self.support or parent.support
            if susceptible == None:
                self.susceptible: float = parent.susceptible
            else:
                self.susceptible:float = susceptible
            if health is None:
                    self.health: float = parent.health
            else:
                self.health: float = health
        else:
            if susceptible == None:
                self.susceptible: float = 1.
            else:
                self.susceptible: float = susceptible
            if health is None:
                    self.health: float = 1.
            else:
                self.health: float = health
        return

    def __copy__(self):
        '''instance copy'''
        return person(parent=self)

