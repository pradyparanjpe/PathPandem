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


from random import shuffle
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
from numpy import abs as npabs
from .pathogen import pathogen
from .person import person


class population(object):
    '''Population class bearing disease spread'''
    def __init__(
            self, people: list=[], infrastructure: float=0, pop_size: int=0,
            p_max: int=10000, serious_health: float=0.3, resist_def: float=0,
            vac_resist: float=0, vac_cov: float=0,
    )-> None:
        self.pop_size = pop_size  # Intermixing Population size
        self.p_max = p_max  # Geographical boundary (x)
        self.serious_health = serious_health  # Life support threshold
        self.resist_def = resist_def # Susceptibility below means resistant
        self.strain_types: list = [None]  # To track evolution of pathogen
        self.infrastructure: float = infrastructure  # Available beds
        self.vac_resist = vac_resist
        self.vac_cov = vac_cov

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
        self.home = npappend(self.home, nparray(
            indiv.home, dtype=npint64).reshape((1, 2)), axis=0)
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

    def mutate(self, in_strain):
        '''Mutation'''
        mutations = 1 + nprandom.random(size=3) * 0.02 - 0.01
        # Then, use each random number in the array
        mut_cfr = in_strain.cfr * mutations[0]
        mut_inf_per_day = in_strain.inf_per_day * mutations[1]
        mut_inf_per_exp = in_strain.inf_per_exp * mutations[2]

        # Compose a mutated strain
        mut_str = pathogen(parent=in_strain, cfr=mut_cfr,
                           day_per_inf=2/mut_inf_per_day,
                           inf_per_exp=mut_inf_per_exp,
        )

        # This strain has now entered the population
        self.strain_types.append(mut_str)
        # Indiv is infected by mutated strain
        return len(self.strain_types) - 1,  mut_cfr,  mut_inf_per_day

    def deposit(self, indiv, i_pos, j_pos)-> bool:
        '''Deposit strain if infected'''
        # Deposit self's pathogen strain at point
        pathy = self.strain_types[self.strain[indiv]]
        if not pathy:  # indiv is carrying infection
            return False
        self.space_contam[i_pos, j_pos] = (
            self.active[indiv] * pathy.persistence)
        self.space_dep_strain[i_pos, j_pos] = (
            self.active[indiv] * self.strain[indiv])
        # Can't collect any more
        return True

    def collect(self, indiv, i_pos, j_pos)-> None:
        '''Collect infection'''
        in_strain = self.strain_types[
            self.space_dep_strain[i_pos, j_pos]]
        if not (in_strain and self.susceptible[indiv]):
            return
        if (nprandom.random()
            > self.susceptible[indiv] * in_strain.inf_per_exp):
            return
        # Possibility of mutation in pathogen
        # (For Future, to simulate evolution of pathogens)
        if nprandom.random() < 0.0001: # Rarely, mutate
            # Motion and probability of mutation arbitrarily chosen
            # (Biological cumulative mutation rates are 10^-6to-7)
            # Cleaner to generate a numpy random array
            pathy_attr = self.mutate(in_strain)
        else:
            # Indiv is infected by old (unmutated strain)
            pathy_attr = (self.strain_types.index(in_strain),
                          in_strain.cfr, in_strain.inf_per_day)
        # Get infected
        self.active[indiv] = True
        self.progress[indiv] = 0.000001
        self.recovered[indiv] = False
        # Some unfortunate indiv will still get infected again
        self.susceptible[indiv] = nprandom.random() * 0.01
        self.strain[indiv], self.cfr[indiv], self.inf_per_day[indiv] =\
            pathy_attr
        return

    def calc_exposure(self, indiv, pos)-> None:
        '''Get exposed or expose the space-cell to infection'''
        # Where is indiv standing?
        i_pos, j_pos = pos[indiv].tolist()
        if self.deposit(indiv, i_pos, j_pos):
            return
        # Collect pathy
        self.collect(indiv, i_pos, j_pos)
        return


    def random_walk(self, d=None, plot_h=None)-> None:
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
                        * self.rms_v.T[:, None])
                * (walk_left != 0).T[:, None],
                dtype=pos.dtype)

            # Can't jump beyond boundary
            # So, reflect exploration
            # pos = pos.clip(min=0, max=self.p_max-1)
            pos = nparray(npnot(pos > (self.p_max-1)) * npabs(pos),
                          dtype=pos.dtype)\
                + nparray((pos > (self.p_max-1)) * (2 * (self.p_max - 1) - pos),
                          dtype=pos.dtype)
            for indiv in range(self.pop_size):
                # TODO: A ufunc or async map would have been faster
                if walk_left[indiv]:
                    self.calc_exposure(indiv, pos)
            if plot_h.contam_dots:
                strain_persist = self.strain_types[-1].persistence
                host_types = []
                host_types.append((pos * (npnot(self.active[:, None])
                                          * self.susceptible[:, None]
                                                  > self.resist_def)).tolist())
                host_types.append((pos * self.active[:, None]).tolist())
                host_types.append((pos * (npnot(self.active[:, None])
                                          * (self.susceptible[:, None]
                                             <= self.resist_def))).tolist())
                pathn_pers = []
                for pers in range(int(strain_persist))[::-1]:
                    pathn_pers.append(npnonzero(self.space_contam==(pers+1)))
                plot_h.update_contam(host_types, pathn_pers)
        return

    def inf_progress(self)-> None:
        '''progress infection every day'''
        # Many logical equations are calculated over numpy ufunc
        # Remember, active, recovered, support are bool

        # Health declines every day
        self.health -= nparray(
            self.active
            * nprandom.random(size=self.pop_size) * self.cfr,
            dtype=self.health.dtype)
        self.progress += nprandom.random(size=self.pop_size)\
            * self.active * self.inf_per_day
        self.progress = self.progress.clip(min=0, max=1)
        self.recovered = nparray(
            npor(self.progress==1, self.recovered), dtype=bool)
        self.active = nparray(
            self.active * npnot(self.recovered), dtype=bool)
        # If recovered, return to original health
        self.health = nparray(npnot(self.active) * (1 - self.comorbidity),
                              dtype=self.health.dtype) \
                              + nparray(self.active * self.health,
                                        dtype=self.health.dtype)

        # If health below threshold, life support is essential
        self.support = self.health < self.serious_health

        # Serious patients do not move
        # self.rms_v = nparray(npnot(self.support) * self.rms_v
        #                      + self.support * 0, dtype=self.rms_v.dtype)
        # self.move_per_day = nparray(npnot(self.support) * self.move_per_day
        #                            + self.support * 0,
        #                            dtype=self.move_per_day.dtype)

        dead_idx = []
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
        self.space_dep_strain = nparray(
            self.space_dep_strain * nparray(self.space_contam, dtype=bool),
            dtype=self.space_dep_strain.dtype
        )
        # Infrastructure may grow, but linearly and very slow
        self.infrastructure = max(min(
            self.infrastructure + 0.2, self.active.sum()/5),
                                  self.infrastructure)

        # Vaccination, when available, happens linearly
        self.susceptible -= nparray(
            self.vac_resist * nparray(
                nprandom.random(self.pop_size) < self.vac_cov, dtype=bool))
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

    def pass_day(self, plot_h=None)-> None:
        '''progress all population and infections'''
        self.random_walk(plot_h=plot_h)  # Macro-scale population
        self.inf_progress()  # Micro-scale: infected individual
        return

