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


from time import sleep
from numpy import append as npappend
from matplotlib import pyplot as plt
from matplotlib.widgets import CheckButtons as mplCheckButtons


class plot_wrap():
    '''Share variables between axes'''
    def __init__(self, space, persistence, pop_size, visualize=False)-> None:
        '''Initiate matplot'''
        self.linetypes = {"Active": "#7F7FFF",
                         "Recovered": "#7FFF7F",
                         "Cases": "#FFFF7F",
                         "Serious": "#FF7F7F",
                         "Dead": "#FFFFFF",
                         "New cases": "#7F7F3F"}
        self.dottypes = {"Uninfected": "#3F3F3FFF",
                         "Infected": "#0000FFFF",
                         "Resistant": "#3FFF3FFF"}
        self.plt = plt
        self.pop_size = pop_size
        self.fig, self.ax = self.plt.subplots(nrows=1, ncols=(visualize + 1))
        self.plt.subplots_adjust(left=0.2)
        self.fig.set_facecolor("#CFCFCFFF")
        self.epidem_ax = None
        self.lines = None
        self.contam_dots = None
        self.contam_ax = None
        self.track_cb_ax = self.plt.axes([0.05, 0.5, 0.1, 0.4],
                                         facecolor="#00000000")
        self.track_cb_ax.text(
            .5,.95,'Plot Track',
            horizontalalignment='center',
            transform=self.track_cb_ax.transAxes,
            fontweight='bold'
        )
        self.track_cb = mplCheckButtons(
            self.track_cb_ax, [k for k in self.linetypes],
            actives=[True] * len(self.linetypes)
        )
        if visualize:
            self.epidem_ax = self.ax[0]
            self.contam_ax = self.ax[1]
            self.init_contam(space, persistence)
            self.person_cb_ax = self.plt.axes([0.05, 0.25, 0.1, 0.25],
                                             facecolor="#00000000")
            self.person_cb_ax.text(
                .5,.9,'Persons',
                horizontalalignment='center',
                transform=self.person_cb_ax.transAxes,
                fontweight='bold'
            )
            self.person_cb = mplCheckButtons(
                self.person_cb_ax,
                [k for k in self.dottypes],
                actives=[True] * len(self.dottypes)
            )
            self.contam_cb_ax = self.plt.axes([0.05, 0.1, 0.1, 0.15],
                                             facecolor="#00000000")
            self.contam_cb_ax.text(
                .5,.85,'Area',
                horizontalalignment='center',
                transform=self.contam_cb_ax.transAxes,
                fontweight='bold'
            )
            self.contam_cb = mplCheckButtons(
                self.contam_cb_ax,
                ["Contaminated"],
                actives=[True]
            )

        else:
            self.epidem_ax = self.ax
        self.init_epidem()
        self.plt.ion()
        self.plt.show(block=False)
        return


    def init_epidem(self):
        '''Initiate matplot'''
        self.news_text = []
        self.lines = []
        for names in self.linetypes:
            line_n, = self.epidem_ax.plot(
                [],[], label=names, color=self.linetypes[names]
            )
            self.lines.append(line_n)
        self.epidem_ax.legend(
            loc='lower left', bbox_to_anchor= (0.0, 1.01),
            ncol=2, borderaxespad=0, frameon=False
        )
        self.epidem_ax.set_facecolor("#000000")
        self.epidem_ax.grid(color="#7f7f7f", linestyle="dotted", linewidth=1)
        self.epidem_ax.set_xlabel("Days")
        self.epidem_ax.set_ylabel("Persons")
        return


    def init_contam(self, space, persistence):
        '''Initiate space-contamination visualization'''
        hosts = []
        for typ in self.dottypes:
            hosts.append(self.contam_ax.scatter(
                [], [] , s=1, c=self.dottypes[typ] , label=typ)
            )
        pathns = []
        for persist in range(persistence):
            colstr = "#FF3F3F" + hex(
                int(0x7F * (1 - persist/persistence)))[2:]
            if not persist:
                pathns.append(self.contam_ax.scatter(
                    [], [], s=1, c=colstr, label="Contaminated space")
                )
            else:
                pathns.append(self.contam_ax.scatter(
                    [], [], s=1, c=colstr)
                )
        self.contam_ax.set_xlim(0, space)
        self.contam_ax.set_ylim(0, space)
        self.contam_ax.set_aspect(1)
        self.contam_ax.legend(
            loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=2,
            borderaxespad=0, frameon=False)
        self.contam_ax.set_facecolor("#000000")
        self.contam_ax.set_xticks([])
        self.contam_ax.set_yticks([])
        self.contam_ax.set_xlabel("West<->East")
        self.contam_ax.set_ylabel("South<->North")
        self.contam_dots = hosts, pathns
        return


    def update_epidem(self, days: int, updates: tuple, newsboard: list,
                      lockdown: int=0, zero_lock: bool=False,
                      early_action: bool=False, intervention: bool=False,
                      vaccined: bool=False, drugged: bool=False) -> None:
        '''Update'''
        if len(updates) == 5:
            for idx, val in enumerate(updates):
                x, y = self.lines[idx].get_data()
                x = npappend(x, days)
                y = npappend(y, val)
                self.lines[idx].set_data((x, y))
            day, cases = self.lines[2].get_data()
            if len(day) == 1:
                new_cases = 0
            else:
                new_cases = cases[-1] - cases[-2]
            x, y = self.lines[5].get_data()
            x = npappend(x, days)
            y = npappend(y, new_cases)
            self.lines[5].set_data((x, y))
            for idx, news in enumerate(newsboard):
                if idx < len(self.news_text):
                    self.news_text[idx].set_text(news)
                else:
                    self.news_text.append(self.epidem_ax.text(
                        .05,.95 - idx * 0.05,
                        news,
                        horizontalalignment='left',
                        transform=self.epidem_ax.transAxes,
                        fontweight='bold',
                        color="#BF7F7FFF"
                    ))
            bgcolor = 0
            label_on = self.track_cb.get_status()
            for idx in range(len(self.lines)):
                self.lines[idx].set_visible(label_on[idx])
            if lockdown or (days < zero_lock and early_action):
                bgcolor += 0x3F0000
            elif intervention or early_action:
                bgcolor += 0x1F00
            if vaccined:
                bgcolor += 0x3F3F00
            if drugged:
                bgcolor += 0x3F
            bgcolor = "#" + "0" * (6 - len(hex(bgcolor)[2:])) + hex(bgcolor)[2:]
            self.epidem_ax.set_facecolor(bgcolor)
            self.mypause(0.0005)
            self.epidem_ax.relim();
            self.epidem_ax.autoscale_view(True, True, True)
        return


    def update_contam(self, host_types: list, pathn_pers: list) -> None:
        '''Update'''
        label_on = self.person_cb.get_status()
        host_scs, pathn_scs = self.contam_dots
        for idx, persist in enumerate(pathn_scs):
            pathn = list(zip(*pathn_pers[idx]))
            if pathn:
                persist.set_offsets(pathn)
            self.contam_dots[1][idx].set_visible(self.contam_cb.get_status()[0])
        for idx, typ in enumerate(host_scs):
            typ.set_offsets(host_types[idx])
            self.contam_dots[0][idx].set_visible(label_on[idx])
        self.mypause(0.0005)
        return

    def mypause(self, interval):
        # Thanks to pepgma at stackoverflow
        manager = self.plt._pylab_helpers.Gcf.get_active()
        if manager is not None:
            canvas = manager.canvas
            if canvas.figure.stale:
                canvas.draw_idle()
                canvas.start_event_loop(interval)
            else:
                sleep(interval)
        return


    def savefig(self, filehandle):
        '''Save figure'''
        self.plt.savefig(filehandle)
        return
