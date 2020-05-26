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
'''Simulate'''

from numpy import array as nparray
from numpy import any as npany
from numpy import append as npappend
from numpy import int16 as npint64
from numpy import random as nprandom
from .definitions import MED_DISCOVERY

def simulate(
        city, logfile, simul_pop, med_eff: float=0.,
        med_recov: float=0, vac_res: float=0, vac_cov: float=0.,
        movement_restrict: int=0, contact_restrict: int=0,
        lockdown_chunk: int=0, lockdown_panic: int=1, seed_inf: int=0,
        zero_lock: bool=False, intervention: bool=False, early_action=False,
        plot_h=None
) -> tuple:
    '''Recursive simulation of each day'''
    vaccine_discovery_date = 0
    drug_discovery_date = 0
    while not vaccine_discovery_date:
        for idx, k in enumerate(nprandom.random(size=(int(5/MED_DISCOVERY)))):
            if k < MED_DISCOVERY:
                vaccine_discovery_date = idx
    while not drug_discovery_date:
        for idx, k in enumerate(nprandom.random(size=(int(5/MED_DISCOVERY)))):
            if k < MED_DISCOVERY:
                drug_discovery_date = idx
    lockdown = 0
    next_lockdown = seed_inf * lockdown_panic
    # Track infection trends
    track: nparray = nparray([[]] * 0, dtype=npint64).reshape((0, 5))
    days = 0
    args = city.survey(simul_pop)
    newsboard = (
        "Total Population %d" % city.pop_size,
        "ICUs Available: %d/%d"
        % (city.infrastructure - args[3], city.infrastructure),
    )
    reaction = (zero_lock, early_action, intervention,
                days > vaccine_discovery_date, days > drug_discovery_date)
    plot_h.update_epidem(days, args, newsboard, lockdown, *reaction)
    track = npappend(track, nparray(args).reshape((1, 5)), axis=0)
    print(*args, file=logfile, flush=True)
    city.pass_day(plot_h)  # IT STARTS!
    while npany(city.space_contam):  # Absent from persons and places
        if days == vaccine_discovery_date:
            city.vaccine_resist = vac_res
            city.vaccine_cov = vac_cov
        if days == drug_discovery_date:
            for idx, pathy in enumerate(city.strain_types):
                if pathy is not None:
                    city.strain_types[idx].inf_per_day /= med_recov
                    city.strain_types[idx].cfr *= med_eff
                    city.inf_per_day /= med_recov
                    city.cfr *= med_eff
        if early_action :
            if not days:
                # Restrict movement
                city.rms_v //= movement_restrict
                city.move_per_day //= contact_restrict
            elif days == zero_lock:
                # End of initial lockdown
                city.rms_v *= movement_restrict
                city.move_per_day *= contact_restrict
        days += 1
        args = city.survey(simul_pop)
        newsboard = (
            "Total Population %d" % city.pop_size,
            "ICUs Available: %d/%d"
            % (city.infrastructure - args[3], city.infrastructure),
        )
        if days > drug_discovery_date:
            newsboard.append("Drug Discovered")
        if days > vaccine_discovery_date:
            newsboard.append("Vaccine Discovered")
        track = npappend(track, nparray(args).reshape((1, 5)), axis=0)
        print(*args, file=logfile, flush=True)
        city.pass_day(plot_h)
        reaction = (zero_lock, early_action, intervention,
                    days > vaccine_discovery_date, days > drug_discovery_date)
        plot_h.update_epidem(days, args, newsboard, lockdown, *reaction)
        if intervention and lockdown == 0 and (args[2] > next_lockdown):
            next_lockdown *= lockdown_panic
            # Panic by infection Spread
            lockdown = 1
            city.rms_v //= movement_restrict
            city.move_per_day //= contact_restrict
        if intervention and lockdown:
            lockdown += 1
        if intervention and lockdown > lockdown_chunk + 1:
            # Business as usual
            city.rms_v *= movement_restrict
            city.move_per_day *= contact_restrict
            lockdown = 0
    args = city.survey(simul_pop)
    newsboard = (
        "Total Population %d" % city.pop_size,
        "ICUs Available: %d/%d"
        % (city.infrastructure - args[3], city.infrastructure),
    )
    if days > drug_discovery_date:
        newsboard.append("Drug Discovered")
    if days > vaccine_discovery_date:
        newsboard.append("Vaccine Discovered")
    track = npappend(track, nparray(args).reshape((1, 5)), axis=0)
    print(*args, file=logfile, flush=True)
    reaction = (zero_lock, early_action, intervention,
                days > vaccine_discovery_date, days > drug_discovery_date)
    plot_h.update_epidem(days, args, newsboard, lockdown, *reaction)
    return

