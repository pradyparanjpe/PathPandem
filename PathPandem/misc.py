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
'''Load DB'''


def ih_translate(cfr: float=0, day_per_inf: int=0)-> float:
    '''Convert raw cfr to Irwin Hall calculated Max limit'''
    # "day_per_inf" sums of Uniform random variable should be greater than 1
    # for "cfr" cases, in others, it should be less than 1.
    # Irwin Hall cdf gives the probability that N uniform random numbers
    # sum to less than t threshold. Calculation involves complications
    # Empirically calculated database can be used to set the correct u-bound
    if cfr == 0:
        return 0
    elif cfr >= 1:
        return 0xFFFF  # A Very Large Number
    from . import __file__
    from os import path
    from pickle import load
    dbfile = path.join(path.dirname(__file__),
                       'reverse_cfr_database.pkl')
    try:
        with open(dbfile, "rb") as dbfile_h:
            IH_DB = load(dbfile_h)
            del path, load, dbfile
    except FileNotFoundError:
        try:
            dbfile = 'reverse_cfr_database.pkl'
            with open(dbfile, "rb") as dbfile_h:
                IH_DB = load(dbfile_h)
                del path, load, dbfile
        except:
            print("Could not find Irwin-Hall calculated database file,\
            using THE INCORRECT mathematical formula", flush=True)
            del path, load
            return 0.385 * pow(cfr, 2.2) + 0.115
    lookup = IH_DB[day_per_inf]
    cfr_keys = list(sorted(lookup.keys()))
    for idx, key in enumerate(cfr_keys):
        if key > cfr and idx:
            return lookup[cfr_keys[idx - 1]]
    return 0.385 * pow(cfr, 2.2) + 0.115

