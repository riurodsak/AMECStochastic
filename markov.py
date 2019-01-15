#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# markov.py
# Created by Riu on 18/12/2018 for the project AMECStochastic
# Simulation of a Markov process. The simulation is done beforehand. The Bokeh visualization is afterwards.

import numpy as np
from scipy.optimize import curve_fit
from pylab import exp, sqrt, diag

from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Panel
from bokeh.plotting import figure

TLIM = 150  # Simulation time limit
N = 100  # Number of trayectories
dt = 0.01
dt_sqrt = np.sqrt(dt)
xinit = 1.0  # Initial value
a = -1.0  # Sign of the linear term in the potential
d_sqrt = np.sqrt(0.2)
var_mean_timeframe = int(20.0 / dt)  # Number of timeframes to consider when determining the stability of the variance


def a_func(_x, _a):  # Potential
    return -_a * _x - np.power(_x, 3)


def gauss(x, mu, sigma, A):  # Gaussian fit function
    return A * exp(-(x - mu) ** 2 / 2 / sigma ** 2)


def bimodal(x, mu1, sigma1, A1, mu2, sigma2, A2):  # Bimodal fit function
    return gauss(x, mu1, sigma1, A1) + gauss(x, mu2, sigma2, A2)


t = np.arange(0, TLIM, dt)
T = len(t)
xlin = np.linspace(-2, 2, 1000)

# Create dataset
x = np.zeros((N, T))
x[:, 0] = np.repeat(xinit, N)
prom = np.zeros(T)
prom[0] = xinit
var = np.zeros(T)
corr = np.zeros(T)
corr[0] = 1

t_rxn_corr = 0
t_rxn_var = 0

for i in range(1, T):
    x[:, i] = x[:, i - 1] + a_func(x[:, i - 1], a) * dt + d_sqrt * dt_sqrt * np.random.normal(0, 1, N)
    prom[i] = np.mean(x[:, i])
    var[i] = np.mean(np.power(x[:, i], 2)) - np.power(prom[i], 2)
    corr[i] = np.mean(x[:, i] * x[:, 0])
    if corr[i] < 0.0 and t_rxn_corr == 0:
        t_rxn_corr = i
    if i > var_mean_timeframe:
        coef = np.polyfit(t[i - var_mean_timeframe:i], var[i - var_mean_timeframe:i], 1)
        if coef[0] < 0.0 and t_rxn_var == 0:
            t_rxn_var = i
    if t_rxn_corr and t_rxn_var:  # if all relaxation times have been asigned, end simulation
        t = t[:i]
        TLIM = t[-1]
        x = x[:, :i]
        prom = prom[:i]
        var = var[:i]
        corr = corr[:i]
        break

expected = (1, .5, 1.4, -1, .5, 0.9)  # starting coefficients to fir the bimodal function
p_eq, edges = np.histogram(x[:, -1], density=True, bins=20)
params, cov = curve_fit(bimodal, edges[:-1], p_eq, expected)
sigma = sqrt(diag(cov))

print("Markov process:")
print("t_rxn_corr: {} \t t_rxn_var: {}".format(t[t_rxn_corr - 1], t[t_rxn_var - 1]))
print("Mean: {0:1.4f} \t Var: {1:1.4f} \t Corr: {2:1.4f}".format(prom[-1], var[-1], corr[-1]))

# BOKEH Visualization ##################################################################################################
# Create sources
tray_source = ColumnDataSource(data=dict(t=t))

for tray_n in range(N):
    tray_source.data[str(tray_n)] = x[tray_n]
analysis_source = ColumnDataSource(data=dict(t=t, prom=prom, var=var, corr=corr))
hist_source = ColumnDataSource(data=dict(bottom=edges[:-1], top=edges[1:], p_eq=p_eq))
gauss1_source = ColumnDataSource(data=dict(x=xlin, y=gauss(xlin, params[0], params[1], params[2])))
gauss2_source = ColumnDataSource(data=dict(x=xlin, y=gauss(xlin, params[3], params[4], params[5])))

# Plots
plot_tray = figure(title='Trayectories',
                   x_axis_label='t', y_axis_label='X(t)',
                   x_range=(0, TLIM), output_backend="svg")
for tray_n in range(20):  # Plot only the first 100 trayectories
    plot_tray.line('t', str(tray_n), source=tray_source, alpha=min(1.0, 10 / 100))
plot_tray.line('t', 'prom', source=analysis_source, color='red', legend='promedio')

plot_hist = figure(title='Equilibrium distribution',
                   x_axis_label='P_eq(x)', y_axis_label='x', y_range=plot_tray.y_range,
                   output_backend="svg")
plot_hist.quad(right='p_eq', left=0, top='top', bottom='bottom', source=hist_source)
plot_hist.line(x='y', y='x', source=gauss1_source, color='red')
plot_hist.line(x='y', y='x', source=gauss2_source, color='green')

plot_var = figure(title='Variance',
                  x_axis_label='t', y_axis_label='var(X(t))', x_range=plot_tray.x_range,
                  output_backend="svg")
plot_var.line('t', 'var', source=analysis_source)

plot_corr = figure(title='Correlation',
                   x_axis_label='t', y_axis_label='corr(X(t))', x_range=plot_tray.x_range,
                   output_backend="svg")
plot_corr.line('t', 'corr', source=analysis_source)


markov_tab = Panel(child=gridplot(children=[[plot_tray, plot_hist, plot_corr],
                                            [plot_var, None, None]],
                                  plot_width=400,
                                  plot_height=400),
                   title='Markov')
