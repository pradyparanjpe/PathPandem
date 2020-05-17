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
'''Simulate spread of an infectious agent'''


from pickle import load
from random import shuffle
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from numpy import round as npround
from numpy import array as nparray
from numpy import float16 as npfloat64
from numpy import intp as npintp
from numpy import int16 as npint64
from numpy import append as npappend
from numpy import delete as npdelete
from numpy import random as nprandom
from numpy import logical_not as npnot
from numpy import logical_and as npand
from numpy import logical_or as npor
from numpy import nonzero as npnonzero
from numpy import any as npany
from matplotlib import pyplot as plt
from gooey import Gooey


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
        if raw_cfr:
            # Correct CFR to adjust for correct probabilities
            # self.cfr =  0.385 * pow(self.cfr, 2.2) + 0.115
            self.cfr = ih_translate(self.cfr/2, day_per_inf)
        return


class person(object):
    '''Person class tracking each Person'''
    def __init__(
            self, parent=None, active: bool=False, recovered: bool=False,
            susceptible: float=1., support = False, health: float=None,
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
            self.home = nparray((0, 0), dtype=npint64).reshape((1, 2))
        else:
            # everyone's init position is uniformly randomly guessed
            self.home = nparray(
                (nprandom.randint(self.p_max), nprandom.randint(self.p_max)),
                dtype=npint64).reshape((1, 2))
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


class population(object):
    '''Population class bearing disease spread'''
    def __init__(self, people: list=[], infrastructure: float=0,
                 pop_size: int=0, p_max=10000,
                 vaccine_resist: float=0, vaccine_cov: float=0):
        self.pop_size = pop_size  # Intermixing Population size
        self.p_max = p_max  # Geographical boundary (x)
        self.strain_types: list = [None]  # To track evolution of pathogen
        self.infrastructure: float = infrastructure  # Available beds
        self.vaccine_resist = vaccine_resist
        self.vaccine_cov = vaccine_cov

        # Use fast numpy ufunc operations on arrays (may be ported to cupy)
        self.active: nparray = nparray([False] * pop_size, dtype=bool)
        self.recovered: nparray = nparray([False] * pop_size, dtype=bool)
        self.susceptible: nparray = nparray([1.] * pop_size, dtype=npfloat64)
        self.health: nparray = nparray([1.] * pop_size, dtype=npfloat64)
        self.support: nparray = nparray([False] * pop_size, dtype=bool)
        self.comorbidity: nparray = nparray([0.] * pop_size, dtype=npfloat64)
        self.progress: nparray = nparray([0.] * pop_size, dtype=npfloat64)
        self.move_per_day: nparray = nparray([], dtype=npint64)
        self.strain: nparray = nparray([0] * pop_size, dtype=npintp)
        self.home: nparray = nparray([[0, 0]] * pop_size,
                                    dtype=npint64).reshape((pop_size, 2))
        self.rms_v: nparray = nparray([0] * pop_size, dtype=npint64)
        self.cfr: nparray = nparray([0] * pop_size, dtype=npfloat64)
        self.inf_per_day: nparray = nparray([0] * pop_size, dtype=npfloat64)

        # Contamination: presence of pathogen in space cell
        # Contamination persists for "int" number of days
        self.space_contam: nparray = nparray([[0] * p_max]
                                             * p_max, dtype=npint64)
        # The pathogen strain (type) present in that strain
        self.space_dep_strain: nparray = nparray([[0] * p_max]
                                                 * p_max, dtype=npint64)
        return

    def __add__(self, indiv: person):
        '''add individuals in population'''
        self.pop_size += 1

        # Append numpy arrays
        self.active = npappend(self.active, indiv.active)
        self.recovered = npappend(self.recovered, indiv.recovered)
        self.susceptible = npappend(self.susceptible, indiv.susceptible)
        self.health = npappend(self.health, indiv.health)
        self.support = npappend(self.support, indiv.support)
        self.comorbidity = npappend(self.comorbidity, indiv.comorbidity)
        self.progress = npappend(self.progress, indiv.progress)
        self.move_per_day = npappend(self.move_per_day, indiv.move_per_day)
        if indiv.strain:
            self.cfr = npappend(self.cfr, indiv.strain.cfr)
            self.inf_per_day = npappend(
                self.inf_per_day, indiv.strain.inf_per_day)
        else:
            self.cfr = npappend(self.cfr, 0)
            self.inf_per_day = npappend(self.inf_per_day, 0)
        if indiv.strain in self.strain_types:
            idx = self.strain_types.index(indiv.strain)
        else:
            self.strain_types.append(indiv.strain)
            idx = len(self.strain_types) - 1
        self.strain = npappend(self.strain, idx)
        self.rms_v = npappend(self.rms_v, indiv.rms_v)
        self.home = npappend(self.home, indiv.home, axis=0)
        return

    def __sub__(self, idx: list):
        '''remove persons by [idx]'''
        idx = list(idx)
        self.pop_size -= len(idx)

        # Delete from numpy array  (We can't track what happened to the dead)
        # Else, remember the dead in a different set of objects
        self.active = npdelete(self.active, idx)
        self.recovered = npdelete(self.recovered, idx)
        self.susceptible = npdelete(self.susceptible, idx)
        self.health = npdelete(self.health, idx)
        self.support = npdelete(self.support, idx)
        self.comorbidity = npdelete(self.comorbidity, idx)
        self.progress = npdelete(self.progress, idx)
        self.move_per_day = npdelete(self.move_per_day, idx)
        self.cfr = npdelete(self.cfr, idx)
        self.inf_per_day = npdelete(self.inf_per_day, idx)
        self.strain = npdelete(self.strain, idx)
        self.home = npdelete(self.home, idx, axis=0)
        self.rms_v = npdelete(self.rms_v, idx)

    def compose_pop(self, person_typ: person=None, person_num: int=0):
        '''Add a "whole-sale" of persons belonging to one type'''
        for _ in range(person_num):
            self + person_typ.__copy__()
        return

    def analyse_person(self, idx: int) -> person:
        '''Extract information from numpy into a person object'''
        p_max = self.p_max

        # From numpy array by index
        active = self.active[idx]
        recovered = self.recovered[idx]
        susceptible = self.susceptible[idx]
        health = self.health[idx]
        support = self.support[idx]
        comorbidity = self.comorbidity[idx]
        progress = self.progress[idx]
        move_per_day = self.move_per_day[idx]
        strain = self.strain_types[self.strain[idx]]
        rms_v = self.rms_v[idx]
        return person(
            active=active, recovered=recovered, susceptible=susceptible,
            health=health, support=support, comorbidity=comorbidity,
            progress=progress, move_per_day=move_per_day, strain=strain,
            p_max=p_max, rms_v=rms_v
        )

    def normalize_pop(self, pop_size=1000000)-> None:
        '''Expand/Shrink population size to pop_size,
        fairly maintaining composition
        '''
        reduction_ratio = int(pop_size / self.pop_size) + 1
        if reduction_ratio > 1:
            for p_num in range(self.pop_size):
                person = self.analyse_person(p_num)
                self.compose_pop(person_typ=person, person_num=reduction_ratio)
        # randomly trim
        trim_num = self.pop_size - pop_size
        if trim_num > 0:
            trim_idx = list(range(self.pop_size))
            shuffle(trim_idx)
            self - trim_idx[:trim_num]
        return

    def calc_exposure(self, indiv, pos)-> None:
        '''Get exposed or expose the space to infection'''
        # Deposit self's pathogen strain at point
        carrier = False  # Not the biological asymptomatic carrier
        i_pos, j_pos = pos[indiv].tolist()
        pathy = self.strain_types[self.strain[indiv]]
        if pathy:  # indiv is carrying infection
            carrier = True
            self.space_contam[i_pos, j_pos] = (
                self.active[indiv] * pathy.persistence)
            self.space_dep_strain[i_pos, j_pos] = (
                self.active[indiv] * self.strain[indiv])
        # Collect pathy
        in_strain = self.strain_types[
            self.space_dep_strain[i_pos, j_pos]]
        if (not(carrier)  # indiv can't renew their infection by depositing
            and self.susceptible[indiv]
            and in_strain):
            if (nprandom.random()
                < self.susceptible[indiv] * in_strain.inf_per_exp):

                # Get infected
                self.active[indiv] = True
                self.progress[indiv] = 0.000001
                self.recovered[indiv] = False
                self.susceptible[indiv] = nprandom.random() * 0.01
                # Some unfortunate indiv still get infected again

                # Possibility of mutation in pathogen
                # (For Future, to simulate evolution of pathogens)
                if nprandom.random() < 0.0001: # Rarely, mutate
                    # Motion and probability of mutation arbitrarily chosen
                    # (Biological cumulative mutation rates are 10^-6to-7)
                    # Cleaner to generate a numpy random array
                    mutations = 1 + nprandom.random(size=4) * 0.02 - 0.01
                    # Then, use each random number in the array
                    mut_cfr = in_strain.cfr * mutations[0]
                    mut_inf_per_day = in_strain.inf_per_day * mutations[1]
                    mut_inf_per_exp = in_strain.inf_per_exp * mutations[2]
                    mut_persistence = in_strain.persistence * mutations[3]

                    # Compose a mutated strain
                    mut_str = pathogen(parent=in_strain, cfr=mut_cfr,
                                       day_per_inf=2/mut_inf_per_day,
                                       inf_per_exp=mut_inf_per_exp,
                                       persistence=mut_persistence)

                    # This strain has now entered the population
                    self.strain_types.append(mut_str)

                    # Indiv is infected by mutated strain
                    self.strain[indiv] = len(self.strain_types) - 1
                    self.cfr[indiv] = mut_cfr
                    self.inf_per_day[indiv] = mut_inf_per_day
                else:

                    # Indiv is infected by old (unmutated strain)
                    self.strain[indiv] = self.strain_types.index(in_strain)
                    self.cfr[indiv] = in_strain.cfr
                    self.inf_per_day[indiv] = in_strain.inf_per_day
        return


    def random_walk(self)-> None:
        '''Let all population walk randomly'''
        walk_left = self.move_per_day.copy()
        # Every day, people start from home
        pos = self.home.copy()
        # Some travel less, some more
        while npany(walk_left):
            walk_left -= 1

            # All randomly move an edge-length
            pos += nparray(
                npround((nprandom.random(size=pos.shape) * 2 - 1)
                        * nparray((self.rms_v, self.rms_v)).T)
                * nparray((walk_left, walk_left), dtype=bool).T,
                dtype=pos.dtype)

            # Can't jump beyond boundary
            pos = pos.clip(min=0, max=self.p_max-1)
            for indiv in range(self.pop_size):
                # TODO: A ufunc would have been faster
                if walk_left[indiv]:
                    self.calc_exposure(indiv, pos)
        return

    def inf_progress(self)-> None:
        '''progress infection every day'''
        # Many logical equations are calculated over numpy ufunc
        # Remember, active, recovered, support are bool

        # Health declines every day
        self.health -= self.active * self.cfr\
            * nprandom.random(size=self.pop_size)
        self.progress += nprandom.random(size=self.pop_size)\
            * self.active * self.inf_per_day
        self.progress = self.progress.clip(min=0, max=1)
        self.recovered = nparray(
            npor(self.progress==1, self.recovered), dtype=bool)
        self.active = nparray(
            self.active * npnot(self.recovered), dtype=bool)
        # If recovered, return to original health
        self.health = nparray(npnot(self.active) * (1 - self.comorbidity),
                              dtype=self.health.dtype)\
                              + nparray(self.active * self.health,
                                        dtype=self.health.dtype)

        # If health below threshold, life support is essential
        self.support = self.health < SERIOUS_HEALTH

        # If support is required but not available, indiv dies
        dead_idx = npnonzero(self.support)[0].tolist()
        shuffle(dead_idx)
        dead_idx = dead_idx[int(self.infrastructure):]

        # If health < 0: death
        dead_idx += npnonzero(self.health <= 0.)[0].tolist()
        dead_idx = list(set(dead_idx))

        # Eliminate dead from population
        if dead_idx:
            self - dead_idx

        # Contamination reduces over time
        self.space_contam -= 1
        self.space_contam = self.space_contam.clip(min=0)
        # Infrastructure may grow, but linearly and very slow
        self.infrastructure = max(min(
            self.infrastructure + 0.2, self.active.sum()/20),
                                  self.infrastructure)

        # Vaccination, when available, happens linearly
        self.susceptible -= nparray(
            self.vaccine_resist * nparray(
                nprandom.random(self.pop_size) < self.vaccine_cov, dtype=bool))
        self.susceptible = self.susceptible.clip(min=0)
        return

    def survey(self, o_size=0) -> tuple:
        '''Testing results: active, recovered, cases, serious, dead
        if original population size(o_size) is provided dead is returned.
        Else, returns negative of current population size.
        '''
        num_active: int = self.active.sum()
        num_recovered: int = self.recovered.sum()
        dead: int = o_size - self.pop_size
        num_cases: int = num_active + num_recovered + dead
        num_serious: int = nparray(self.support).sum()
        return num_active, num_recovered, num_cases, num_serious, dead

    def pass_day(self)-> None:
        '''progress all population and infections'''
        self.random_walk()  # Macro-scale population
        self.inf_progress()  # Micro-scale: infected individual
        return


def update_plot(fig, lines, tp: int, updates: tuple, lockdown: int=0) -> None:
    '''Update'''
    if len(updates) == 5:
        for idx, val in enumerate(updates):
            x, y = lines[idx].get_data()
            x = npappend(x, tp)
            y = npappend(y, val)
            lines[idx].set_data((x, y))
        day, cases = lines[2].get_data()
        if len(day) == 1:
            new_cases = 0
        else:
            new_cases = cases[-1] - cases[-2]
        x, y = lines[5].get_data()
        x = npappend(x, tp)
        y = npappend(y, new_cases)
        lines[5].set_data((x, y))
        bgcolor = 0
        if lockdown or (DAYS < ZERO_LOCK and EARLY_ACTION):
            bgcolor += 0x3F0000
        elif INTERVENTION or EARLY_ACTION:
            bgcolor += 0x1F00
        if VACCINED:
            bgcolor += 0x3F3F00
        if DRUGGED:
            bgcolor += 0x3F
        bgcolor = "#" + "0" * (6 - len(hex(bgcolor)[2:])) + hex(bgcolor)[2:]
        ax.set_facecolor(bgcolor)
        fig.canvas.draw()
        ax.relim();
        ax.autoscale_view(True, True, True)
        plt.xscale("linear")
        plt.yscale("linear")
        plt.pause(0.5)
    return

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
    elif IH_DB is None:
        # Discouraged
        return 0.385 * pow(cfr, 2.2) + 0.115
    else:
        lookup = IH_DB[day_per_inf]
        cfr_keys = list(sorted(lookup.keys()))
        for idx, key in enumerate(cfr_keys):
            if key > cfr and idx:
                return lookup[cfr_keys[idx - 1]]
    return 0.385 * pow(cfr, 2.2) + 0.115


@Gooey
def cli()-> tuple:
    '''cli inputs'''
    parser: ArgumentParser = ArgumentParser(
        description="Simulate spread of a disease",
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-e", "--early-action", action='store_true',
                        help="Community takes early action")
    parser.add_argument("-i", "--intermediate-action", action='store_true',
                        help="Community locks down intermittently")
    parser.add_argument("-P", "--population", type=int, default=5000,
                        help="Population to simulate, sugest: <50000")
    parser.add_argument(
        "-D", "--density", type=float, default=0.0004,
        help='Density of population people/mtr_sqr, to calculate roaming area')
    parser.add_argument("-I", "--infrastructure", default="0.00053",
                        type=float, help="Beds available per person")
    parser.add_argument("-C", "--contacts", default=50, type=int,
                        help="Contacts per person per day")
    parser.add_argument("-T", "--travel", default=40, type=int,
                        help="Distance travelled between two contacts")
    parser.add_argument("-S", "--seed-infect", default=25, type=int,
                        help="Initial infected number that starts spreading")
    parser.add_argument("-F", "--feeble-percent", default=0, type=float,
                        help="Initial % of population with comorbidity")
    parser.add_argument("-m", "--comorbidity", default=65, type=float,
                        help="% Health loss on account of comorbidity")
    parser.add_argument("-R", "--resist-percent", default=0, type=float,
                        help="Initial % of resistant people")
    parser.add_argument("-r","--resistance", default=90, type=float,
                        help="% reduced susceptibility in resistanant people")
    parser.add_argument("-c", "--case-fatality", default=3.0, type=float,
                        help="% of cases turning fatal")
    parser.add_argument("-p", "--persistence", default=2, type=int,
                        help="Days for which pathogen stays viable on surface")
    parser.add_argument("-d", "--days-per-inf", default=10, type=int,
                        help="Average days of infection before clearance")
    parser.add_argument("-H", "--serious-health", default=10, type=float,
                        help="Below this %, life support is essential")
    parser.add_argument("-E", "--efficiency", default=3, type=float,
                        help="(Efficiency) Infection per 100 exposures")
    parser.add_argument("-w", "--worried-movement-ratio", default=5, type=int,
                        help="Fold reduction in movement")
    parser.add_argument("-W", "--worried-contact-ratio", default=5, type=int,
                        help="Fold reduction in contacts")
    parser.add_argument("-L", "--lockdown-chunk", default=14, type=int,
                        help="Number of days per lockdown")
    parser.add_argument("-l", "--lockdown-panic", default=5, type=int,
                        help="Fold change of cases when people panic-lockdown")
    parser.add_argument("-Z", "--zero-lockdown", default=70, type=int,
                        help="Zero lock down period in case of early action")
    parser.add_argument("-V", "--vaccine-resistance", default=90, type=float,
                        help="Resistance due to vaccination")
    parser.add_argument("-v", "--vaccine-coverage", default=3, type=float,
                        help="% Population that's vaccinated per day")
    parser.add_argument("-M", "--medicine-effect", default=80, type=float,
                        help="% reduction in health deterioration due to virus")
    parser.add_argument("-f", "--fast-recover", default=50, type=float,
                        help="% reduction in days on infection")
    args = parser.parse_args()
    return (args.population, args.density, args.infrastructure, args.contacts,
            args.travel, args.seed_infect, args.feeble_percent/100,
            args.comorbidity/100, args.resist_percent/100, args.resistance/100,
            args.case_fatality/100, args.persistence, args.days_per_inf,
            args.serious_health/100, args.efficiency/100,
            args.worried_movement_ratio, args.worried_contact_ratio,
            args.lockdown_chunk, args.lockdown_panic, args.zero_lockdown,
            args.early_action, args.intermediate_action,
            args.vaccine_resistance/100, args.vaccine_coverage/100,
            1 - args.medicine_effect/100, 1 - args.fast_recover/100
    )


if __name__ == "__main__":
    SIMUL_POP, POP_DENSE, INFRA, MOVE_PER_DAY, RMS_V, SEED_INF,\
        FEEBLE_PROP, COMORBIDITY, RESIST_PROP, RESISTANCE, CFR,\
        PERSISTENCE, DAY_PER_INF, SERIOUS_HEALTH, INF_PER_EXP,\
        MOVEMENT_RESTRICT, CONTACT_RESTRICT, LOCKDOWN_CHUNK,\
        LOCKDOWN_PANIC, ZERO_LOCK, EARLY_ACTION, INTERVENTION,\
        VAC_RES, VAC_COV, MED_EFF, MED_RECOV = cli()

    assert all(
        k < 1 for k in
        [FEEBLE_PROP, COMORBIDITY, RESIST_PROP, RESISTANCE, CFR, SERIOUS_HEALTH,
         INF_PER_EXP]
    )
    assert all(0 < k < 1 for k in [VAC_RES, VAC_COV, MED_EFF, MED_RECOV])
    # INITS
    MAX_SPACE = int(pow(SIMUL_POP/POP_DENSE, 0.5))  # Sqr_mtr
    INFRA *= SIMUL_POP
    IH_DB = None
    LOCKDOWN = 0
    NEXT_LOCKDOWN = SEED_INF * LOCKDOWN_PANIC
    ORDINARY_PROP = round(
        (1 - FEEBLE_PROP - RESIST_PROP) * SIMUL_POP) - SEED_INF
    FEEBLE_PROP = round(FEEBLE_PROP * SIMUL_POP)
    RESIST_PROP = round(RESIST_PROP * SIMUL_POP)
    VACCINED = False
    DRUGGED = False

    # Log raw survey numbers
    FNAME_BASE= int(EARLY_ACTION) * "Early_acted_"\
            + int(INTERVENTION) * "Intervened_"\
            + int(not(INTERVENTION or EARLY_ACTION)) * "Uncontrolled_"
    LOGFILE = open("%sdisease_spread.log" %FNAME_BASE, "w")
    try:
        db_file_h = open("reverse_cfr_database.pkl", "rb")
        IH_DB = load(db_file_h)
        db_file_h.close()
        del db_file_h
    except:
        print("Could not find Irwin-Hall calculated database file,\
        using THE INCORRECT mathematical formula", flush=True)
    print("Active(INST)", "Recovered", "Cases", "Critical(INST)", "Deaths",
            file=LOGFILE, flush=True)

    # Init plot
    fig, ax = plt.subplots()
    lines = []
    HOSTTYPE = {"active": "#7F7FFF",
                "recovered": "#7FFF7F",
                "cases": "#FFFF7F",
                "serious/critical": "#FF7F7F",
                "dead": "#FFFFFF",
                "new cases": "#7F7F3F"}
    for names in HOSTTYPE:
        line_n, = ax.plot([],[], label=names, color=HOSTTYPE[names])
        lines.append(line_n)
    ax.legend()
    ax.set_facecolor("#000000")
    ax.grid(color="#7f7f7f", linestyle="dotted", linewidth=1)
    ax.set_xlabel("Days")
    ax.set_ylabel("Persons")
    plt.ion()

    # INIT pathogen, host-type
    PATHY = pathogen(cfr=CFR, raw_cfr=True, day_per_inf=DAY_PER_INF,
                     inf_per_exp=INF_PER_EXP, persistence = PERSISTENCE)
    ORDINARY_IMMUN = person(susceptible=1, move_per_day=MOVE_PER_DAY,
                            rms_v=RMS_V, p_max=MAX_SPACE)
    FEEBLE = person(parent=ORDINARY_IMMUN, comorbidity=COMORBIDITY)
    RESIST_IMMUN = person(parent=ORDINARY_IMMUN, susceptible=(1-RESISTANCE))

    # Founder of infection
    FOUNDER = person(parent=ORDINARY_IMMUN, active=True,
                     progress=0.0001, strain=PATHY)
    CITY = population(infrastructure=INFRA, p_max=MAX_SPACE)
    CITY.compose_pop(ORDINARY_IMMUN, ORDINARY_PROP)
    CITY.compose_pop(RESIST_IMMUN, RESIST_PROP)
    CITY.compose_pop(FEEBLE, FEEBLE_PROP)
    CITY.compose_pop(FOUNDER, SEED_INF)
    SIMUL_POP = CITY.pop_size

    # Track infection trends
    TRACK: nparray = nparray([[]] * 0, dtype=npint64).reshape((0, 5))
    DAYS = 0
    args = CITY.survey(SIMUL_POP)
    update_plot(fig, lines, 0, args, LOCKDOWN)
    TRACK = npappend(TRACK, nparray(args).reshape((1, 5)), axis=0)
    print(*args, file=LOGFILE, flush=True)
    CITY.pass_day()  # IT STARTS!
    while npany(CITY.space_contam):  # Absend from persons and places
        if nprandom.random() < 0.002 and not VACCINED:
            VACCINED= True
            CITY.vaccine_resist = VAC_RES,
            CITY.vaccine_cov = VAC_COV
        if nprandom.random() < 0.002 and not DRUGGED:
            DRUGGED = True
            for idx, pathy in enumerate(CITY.strain_types):
                if pathy is not None:
                    CITY.strain_types[idx].inf_per_day /= MED_RECOV
                    CITY.strain_types[idx].cfr *= MED_EFF
                    CITY.inf_per_day /= MED_RECOV
                    CITY.cfr *= MED_EFF
        if EARLY_ACTION :
            if not DAYS:
                # Restrict movement
                CITY.rms_v //= MOVEMENT_RESTRICT
                CITY.move_per_day //= CONTACT_RESTRICT
            elif DAYS == ZERO_LOCK:
                # End of initial lockdown
                CITY.rms_v *= MOVEMENT_RESTRICT
                CITY.move_per_day *= CONTACT_RESTRICT
        DAYS += 1
        args = CITY.survey(SIMUL_POP)
        TRACK = npappend(TRACK, nparray(args).reshape((1, 5)), axis=0)
        print(*args, file=LOGFILE, flush=True)
        CITY.pass_day()
        update_plot(fig, lines, DAYS, args, LOCKDOWN)
        if INTERVENTION and LOCKDOWN == 0 and (args[2] > NEXT_LOCKDOWN):
            NEXT_LOCKDOWN *= LOCKDOWN_PANIC
            # Panic by infection Spread
            LOCKDOWN = 1
            CITY.rms_v //= MOVEMENT_RESTRICT
            CITY.move_per_day //= CONTACT_RESTRICT
        if INTERVENTION and LOCKDOWN:
            LOCKDOWN += 1
        if INTERVENTION and LOCKDOWN > LOCKDOWN_CHUNK + 1:
            # Business as usual
            CITY.rms_v *= MOVEMENT_RESTRICT
            CITY.move_per_day *= CONTACT_RESTRICT
            LOCKDOWN = 0
    args = CITY.survey(SIMUL_POP)
    TRACK = npappend(TRACK, nparray(args).reshape((1, 5)), axis=0)
    print(*args, file=LOGFILE, flush=True)
    update_plot(fig, lines, DAYS, args, LOCKDOWN)
    plt.savefig("%sdisease_plot.jpg"%FNAME_BASE)
    LOGFILE.close()
    try:
        exit()
    except:
        pass

