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
'''Compose a homogenous population'''

from .cli import cli
from .pathogen import pathogen
from .person import person
from .population import population

def compose_homogenous(
        simul_pop: int=0, pop_dense: float=0., infra: float=0.,
        move_per_day: int=0, rms_v: int=0, serious_health: float=0.,
        seed_inf: int=0, feeble_prop: float=0., comorbidity: float=0.,
        resist_prop: float=0., resistance: float=0., cfr: float=0.,
        day_per_inf: int=0, inf_per_exp: float=0., persistence: int=0,
        vac_res: float=0, vac_cov: float=0, resist_def: float=0,
) -> tuple:
    '''A homogenous population'''
    # INITS
    max_space = int(pow(simul_pop/pop_dense, 0.5))  # Sqr_mtr
    infra *= simul_pop
    feeble_prop = round(feeble_prop * simul_pop)
    resist_prop = round(resist_prop * simul_pop)
    ordinary_prop = simul_pop - feeble_prop - resist_prop- seed_inf

    # INIT pathogen, host-type
    pathy = pathogen(cfr=cfr, raw_cfr=True, day_per_inf=day_per_inf,
                     inf_per_exp=inf_per_exp, persistence = persistence)
    ordinary_immun = person(susceptible=1, move_per_day=move_per_day,
                            rms_v=rms_v, p_max=max_space)
    feeble = person(parent=ordinary_immun, comorbidity=comorbidity)
    resest_immun = person(parent=ordinary_immun, susceptible=resist_def)

    # Founder of infection
    founder = person(parent=ordinary_immun, active=True,
                     progress=0.0001, strain=pathy)
    city = population(infrastructure=infra, p_max=max_space,
                      serious_health=serious_health, resist_def=resist_def)
    city.compose_pop(resest_immun, resist_prop)
    city.compose_pop(ordinary_immun, ordinary_prop)
    city.compose_pop(feeble, feeble_prop)
    city.compose_pop(founder, seed_inf)
    simul_pop = city.pop_size
    return city, simul_pop, pathy, max_space
