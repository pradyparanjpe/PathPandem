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
'''Visualization'''


from numpy import append as npappend
from matplotlib import pyplot as plt


def init_plot(space, persistence, visualize=True):
    '''Initiate matplot'''
    fig, ax = plt.subplots(nrows=1, ncols=(visualize + 1))
    fig.set_facecolor("#7F7F7F7F")
    if visualize:
        epidem_ax = ax[0]
        contam_ax = ax[1]
        contam_dots = init_contam(contam_ax, space, persistence)
    else:
        contam_ax, contam_dots, epidem_ax = None, None, ax
    epidem_lines = init_epidem(epidem_ax)
    plt.ion()
    return plt, fig, epidem_ax, contam_ax, epidem_lines, contam_dots


def init_epidem(ax):
    '''Initiate matplot'''
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
    ax.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=2,
        borderaxespad=0, frameon=False)
    # ax.legend()
    ax.set_facecolor("#000000")
    ax.grid(color="#7f7f7f", linestyle="dotted", linewidth=1)
    ax.set_xlabel("Days")
    ax.set_ylabel("Persons")
    return lines


def init_contam(ax, space, persistence):
    '''Initiate space-contamination visualization'''
    hosts = []
    hosttypes = {"unaffected person": "#3F3F3FFF",
                 "carrier person": "#0000FFFF",
                 "resistant person": "#3FFF3FFF"}
    for typ in hosttypes:
        hosts.append(ax.scatter([], [] , s=1, c=hosttypes[typ] , label=typ))
    pathns = []
    for persist in range(persistence):
        colstr = "#FF3F3F" + hex(int(0x7F * (1 - persist/persistence)))[2:]
        if not persist:
            pathns.append(ax.scatter([], [], s=1, c=colstr,
                                     label="Contaminated space"))
        else:
            pathns.append(ax.scatter([], [], s=1, c=colstr))
    ax.set_xlim(0, space)
    ax.set_ylim(0, space)
    ax.set_aspect(1)
    ax.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=2,
        borderaxespad=0, frameon=False)
    # ax.legend()
    ax.set_facecolor("#000000")
    ax.grid(color="#7f7f7f", linestyle="dotted", linewidth=1)
    ax.set_xlabel("East<->West")
    ax.set_ylabel("North<->South")
    return hosts, pathns


def update_epidem(plt, fig, ax, lines, tp: int, updates: tuple, lockdown: int=0,
                days: int=0, zero_lock: bool=False, early_action: bool=False,
                intervention: bool=False, vaccined: bool=False,
                drugged: bool=False) -> None:
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
        if lockdown or (days < zero_lock and early_action):
            bgcolor += 0x3F0000
        elif intervention or early_action:
            bgcolor += 0x1F00
        if vaccined:
            bgcolor += 0x3F3F00
        if drugged:
            bgcolor += 0x3F
        bgcolor = "#" + "0" * (6 - len(hex(bgcolor)[2:])) + hex(bgcolor)[2:]
        ax.set_facecolor(bgcolor)
        fig.canvas.draw()
        ax.relim();
        ax.autoscale_view(True, True, True)
        plt.pause(0.0005)
    return


def update_contam(plt, fig, ax, dots,
                  host_types: list, pathn_pers: list) -> None:
    '''Update'''
    host_scs, pathn_scs = dots
    for idx, typ in enumerate(host_scs):
        typ.set_offsets(host_types[idx])
    for idx, persist in enumerate(pathn_scs):
        pathn = list(zip(*pathn_pers[idx]))
        if pathn:
            persist.set_offsets(pathn)
    fig.canvas.draw()
    plt.pause(0.0005)
    return
