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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the # GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PathPandem.  If not, see <https://www.gnu.org/licenses/>.
'''Simulate spread of an infectious agent'''


from sys import exit as sysexit
from pickle import load
from .cli import cli
from .population import population
from .pathogen import pathogen
from .person import person
from .misc import ih_translate
from .simul import simulate
from .compose_pop import compose_homogenous
from .plot import plot_wrap


def main():
    SIMUL_POP, POP_DENSE, INFRA, MOVE_PER_DAY, RMS_V, SEED_INF,\
        FEEBLE_PROP, COMORBIDITY, RESIST_PROP, RESISTANCE, CFR,\
        PERSISTENCE, DAY_PER_INF, SERIOUS_HEALTH, INF_PER_EXP,\
        MOVEMENT_RESTRICT, CONTACT_RESTRICT, LOCKDOWN_CHUNK,\
        LOCKDOWN_PANIC, ZERO_LOCK, EARLY_ACTION, INTERVENTION,\
        VAC_RES, VAC_COV, MED_EFF, MED_RECOV, VISUALIZE = cli()

    # INITS
    MAX_SPACE = int(pow(SIMUL_POP/POP_DENSE, 0.5))  # Sqr_mtr

    if SIMUL_POP <= SEED_INF:
        print("Total population must be less than Seed founder population, resetting...")
    if MAX_SPACE < 2:
        print('''Too less space to simulate
        Increase Population size OR
        Decrease Density of population''')
        sysexit(1)

    CFR = ih_translate(cfr=CFR/2, day_per_inf=int(DAY_PER_INF*1.414))
    # Log raw survey numbers
    FNAME_BASE= int(EARLY_ACTION) * "Early_acted_"\
            + int(INTERVENTION) * "Intervened_"\
            + int(not(INTERVENTION or EARLY_ACTION)) * "Uncontrolled_"
    LOGFILE = open("%sdisease_spread.log" %FNAME_BASE, "w")
    print("Active(INST)", "Recovered", "Cases", "Critical(INST)", "Deaths",
          file=LOGFILE, flush=True)

    # INIT pathogen, host-type
    CITY, SIMUL_POP, PATHY, SPACE = compose_homogenous(
        simul_pop=SIMUL_POP, pop_dense=POP_DENSE, infra=INFRA,
        move_per_day=MOVE_PER_DAY, rms_v=RMS_V, serious_health=SERIOUS_HEALTH,
        seed_inf=SEED_INF, feeble_prop=FEEBLE_PROP, comorbidity=COMORBIDITY,
        resist_prop=RESIST_PROP, resistance=RESISTANCE, cfr=CFR,
        day_per_inf=DAY_PER_INF, inf_per_exp=INF_PER_EXP,
        persistence=PERSISTENCE, vac_res=VAC_RES, vac_cov=VAC_COV,
        resist_def=(1 - RESISTANCE),
    )
    PLOT_H = plot_wrap(SPACE, PERSISTENCE, SIMUL_POP, VISUALIZE)
    err = simulate(
        city=CITY, logfile=LOGFILE, simul_pop=SIMUL_POP, med_recov=MED_RECOV, med_eff=MED_EFF,
        vac_res=VAC_RES, vac_cov=VAC_COV, movement_restrict=MOVEMENT_RESTRICT,
        contact_restrict=CONTACT_RESTRICT, lockdown_chunk=LOCKDOWN_CHUNK,
        lockdown_panic=LOCKDOWN_PANIC, seed_inf=SEED_INF, zero_lock=ZERO_LOCK,
        intervention=INTERVENTION, early_action=EARLY_ACTION, plot_h=PLOT_H,
    )

    # Finally, save
    PLOT_H.savefig("%sdisease_plot.jpg"%FNAME_BASE)
    LOGFILE.close()
    sysexit(0)

if __name__ == "__main__":
    main()
