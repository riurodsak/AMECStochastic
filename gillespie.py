#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# gillespie.py
# Created by Riu on 18/12/2018 for the project AMECStochastic
# This script implements a parallelized Gillespie algorithm to simulate the number of balls in two boxes across time.

import numpy as np

from bokeh.layouts import row
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Panel
from bokeh.plotting import figure

TLIM = 50  # Simulation time limit
dt = 0.002
N = 30  # Number of trayectories
k = 0.06
N1_init = 1000  # Number of balls in box 1
N2_init = 500   # Number of balls in box 2
N_tot = N1_init + N2_init  # Total number of balls constant
lamda = N_tot * k
a_n = lamda

t = np.arange(0, TLIM, dt)
T = len(t)
var_mean_timeframe = int(5.0 / dt)

# Generate dataset
N1 = np.zeros((N, T))  # Number of balls in box 1 across time for each trayectory.
N1[:, 0] = np.repeat(N1_init, N)
prom = np.zeros(T)  # Mean value of balls in box 1 across time.
prom[0] = N1_init
std_dev = np.zeros(T)
tau = np.zeros(N)  # time for next reaction for each trayectory

t_rxn = 0  # relaxation time

for i in range(1, T):
    for n in range(N):
        if t[i] >= tau[n]:  # reaction happens
            if np.random.uniform(0, 1) < N1[n, i - 1] / N_tot:
                N1[n, i] = N1[n, i - 1] - 1  # box 1 -> box 2
            else:
                N1[n, i] = N1[n, i - 1] + 1  # box 2 -> box 1
            tau[n] = t[i] - np.log(np.random.uniform(0, 1)) / a_n  # determine time for next reaction
        else:
            N1[n, i] = N1[n, i - 1]
    prom[i] = np.mean(N1[:, i])
    std_dev[i] = np.sqrt(np.mean(np.power(N1[:, i], 2)) - np.power(prom[i], 2))
    
    if i > var_mean_timeframe:
        coef = np.polyfit(t[i - var_mean_timeframe:i], prom[i - var_mean_timeframe:i], 1)
        if abs(coef[0]) < 0.01 and t_rxn == 0:  # if the variance has stabilized itself, stop the simulation
            t_rxn = i
            t = t[:i]
            TLIM = t[-1]
            N1 = N1[:, :i]
            prom = prom[:i]
            std_dev = std_dev[:i]
            break

print("Gillespie algorithm:")
print("t_rxn: {}".format(t[t_rxn-1]))
print("Mean: {0:5.2f} \t Std: {1:5.2f}".format(prom[-1], std_dev[-1]))
p_eq, edges = np.histogram(N1[:, -1], density=True, bins=range(1, 1500))


# Create data sources
tray_source = ColumnDataSource(data=dict(t=t))
for tray_n in range(N):
    tray_source.data[str(tray_n)] = N1[tray_n]
analysis_source = ColumnDataSource(data=dict(t=t, tdt=t + dt, prom=prom, std_dev=std_dev,
                                             top=prom + std_dev, bottom=prom - std_dev))
hist_source = ColumnDataSource(data=dict(bottom=edges[:-1], top=edges[1:], p_eq=p_eq))

# Plot trayectories
plot_tray = figure(title='Trayectories',
                   x_axis_label='t', y_axis_label='N1(t)',
                   x_range=(0, TLIM), y_range=(600, 1100),
                   output_backend="svg")

plot_tray.line('t', 'top', line_color='black', line_dash='dashed', legend='std', source=analysis_source)
plot_tray.line('t', 'bottom', line_color='black', line_dash='dashed', source=analysis_source)
for tray_n in range(N):  # plot only the first N trayectories
    plot_tray.line('t', str(tray_n), source=tray_source, alpha=0.08)

plot_tray.line('t', 'prom', source=analysis_source, color='red', legend='mean')

# Plot equlibrium histogram
plot_hist = figure(title='Equilibrium distribution',
                   x_axis_label='P_eq(x)', y_axis_label='x', y_range=plot_tray.y_range,
                   output_backend="svg")
plot_hist.quad(right='p_eq', left=0, top='top', bottom='bottom', source=hist_source)

gillespie_tab = Panel(child=row(plot_tray, plot_hist), title='Gillespie')
