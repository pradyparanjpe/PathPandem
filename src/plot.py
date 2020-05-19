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


def init_plot(space):
    '''Initiate matplot'''
    fig, ax = plt.subplots(2)
    epidem_ax = ax[0]
    epidem_lines = init_epidem(epidem_ax)
    contam_ax = ax[1]
    contam_dots = init_contam(contam_ax, space)
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
    ax.legend()
    ax.set_facecolor("#000000")
    ax.grid(color="#7f7f7f", linestyle="dotted", linewidth=1)
    ax.set_xlabel("Days")
    ax.set_ylabel("Persons")
    return lines


def init_contam(ax, space):
    '''Initiate space-contamination visualization'''
    host = ax.scatter([], [] , s=1, c="#FFFF3FFF", label="Carrier")
    pathn = ax.scatter([], [], s=1, c="#FF3F3F7F", label="Contaminated")
    ax.set_xlim(0, space)
    ax.set_ylim(0, space)
    ax.legend()
    ax.set_facecolor("#000000")
    ax.grid(color="#7f7f7f", linestyle="dotted", linewidth=1)
    ax.set_xlabel("East<->West")
    ax.set_ylabel("North<->South")
    return host, pathn


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


def update_contam(plt, fig, ax, dots, host_pos: list, pathn_pos) -> None:
    '''Update'''
    pathn = list(zip(pathn_pos[0].tolist(), pathn_pos[1].tolist()))
    host_sc, pathn_sc = dots
    pathn_sc.set_offsets(pathn)
    host_sc.set_offsets(host_pos)
    fig.canvas.draw()
    plt.pause(0.0005)
    return
