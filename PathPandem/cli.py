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
'''UI'''


from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from gooey import Gooey

@Gooey
def cli()-> tuple:
    '''cli inputs'''
    parser: ArgumentParser = ArgumentParser(
        description="Simulate spread of a disease",
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--intermediate-action", action='store_true',
                        help="Community locks down intermittently")
    parser.add_argument("-e", "--early-action", action='store_true',
                        help="Community takes early action")
    parser.add_argument("-g", "--graphical-visualization", action='store_true',
                        help="Visualize population movements, very slow")
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
    parser.add_argument("-F", "--feeble-percent", default=5, type=float,
                        help="Initial % of population with comorbidity")
    parser.add_argument("-m", "--comorbidity", default=65, type=float,
                        help="% Health loss on account of comorbidity")
    parser.add_argument("-R", "--resist-percent", default=0, type=float,
                        help="Initial % of resistant people")
    parser.add_argument("-r","--resistance", default=90, type=float,
                        help="% reduced susceptibility in resistanant people")
    parser.add_argument("-c", "--case-fatality", default=0.3, type=float,
                        help="% of cases turning fatal")
    parser.add_argument("-p", "--persistence", default=3, type=int,
                        help="Days for which pathogen stays viable on surface")
    parser.add_argument("-d", "--days-per-inf", default=10, type=int,
                        help="Average days of infection before clearance")
    parser.add_argument("-H", "--serious-health", default=30, type=float,
                        help="Below this %, life support is essential")
    parser.add_argument("-E", "--efficiency", default=7, type=float,
                        help="(Efficiency) Infection per 100 exposures")
    parser.add_argument("-w", "--worried-movement-ratio", default=3, type=int,
                        help="Fold reduction in movement")
    parser.add_argument("-W", "--worried-contact-ratio", default=3, type=int,
                        help="Fold reduction in contacts")
    parser.add_argument("-L", "--lockdown-chunk", default=5, type=int,
                        help="Number of days per lockdown")
    parser.add_argument("-l", "--lockdown-panic", default=5, type=int,
                        help="Fold change of cases when people panic-lockdown")
    parser.add_argument("-Z", "--zero-lockdown", default=10, type=int,
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
            1 - args.medicine_effect/100, 1 - args.fast_recover/100,
            args.graphical_visualization,
    )

