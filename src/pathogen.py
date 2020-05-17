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



class pathogen(object):
    '''Properties of pathogen'''
    def __init__(self, parent=None, raw_cfr=False, cfr: float=0.,
                 day_per_inf: int=0, inf_per_exp: float = 0,
                 persistence = 0)-> None:
        '''Initialize a (Null) Pathogen'''
        self.cfr: float = cfr
        self.inf_per_day: float = 2/day_per_inf  # inverse of days of inf
        self.inf_per_exp: float = inf_per_exp  # Exposure causes infection
        self.persistence: int = persistence  # How long viable in wild
        if parent:
            self.cfr = self.cfr or parent.cfr
            self.inf_per_day = self.inf_per_day or parent.inf_per_day
            self.inf_per_exp = self.inf_per_exp or parent.inf_per_exp
            self.persistence = self.persistence or parent.persistence
        return
