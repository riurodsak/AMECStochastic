#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# markov_dynamic.py
# Created by Riu on 18/12/2018 for the project AMECStochastic
# Simulation of a Markov process. The simulation is done dynamically with user interaction within Bokeh.

import numpy as np
from scipy.stats import norm

from bokeh.layouts import column, gridplot, widgetbox
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Panel, Button, TextInput
from bokeh.plotting import figure

from utils import LatexLabel

TLIM = 10  # Simulation time limit
N = 20  # Number of trayectories
dt = 0.01
dt_sqrt = np.sqrt(dt)
xinit = 1.0
a = -1.0
d_sqrt = np.sqrt(0.2)
var_mean_timeframe = int(2.0 / dt)  # Number of timeframes to consider when determining the stability of the variance


def a_func(_x, _a):
    return -_a * _x - np.power(_x, 3)


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

mean, std = norm.fit(x[:, -1])
p_eq, edges = np.histogram(x[:, -1], density=True, bins=20)

print("Markov process:")
print("t_rxn_corr: {0:3.2f} \t t_rxn_var: {1:3.2f}".format(t[t_rxn_corr - 1], t[t_rxn_var - 1]))
print("Mean: {0:1.4f} \t Var: {1:1.4f} \t Corr: {2:1.4f}".format(prom[-1], var[-1], corr[-1]))

# BOKEH ################################################################################################################
# Create sources
tray_source = ColumnDataSource(data=dict(t=t))

for tray_n in range(N):
    tray_source.data[str(tray_n)] = x[tray_n]
analysis_source = ColumnDataSource(data=dict(t=t, prom=prom, var=var, corr=corr))
hist_source = ColumnDataSource(data=dict(bottom=edges[:-1], top=edges[1:], p_eq=p_eq))
gauss_source = ColumnDataSource(data=dict(x=xlin, y=norm.pdf(xlin, mean, std)))


# Update plot routine
def update():
    a = float(a_input.value)
    d = float(d_input.value)
    dt = float(dt_input.value)
    xinit = float(xinit_input.value)
    dt_sqrt = np.sqrt(dt)
    d_sqrt = np.sqrt(d)

    x[:, 0] = np.repeat(xinit, N)
    for i in range(1, T):
        x[:, i] = x[:, i - 1] + a_func(x[:, i - 1], a) * dt + d_sqrt * dt_sqrt * np.random.normal(0, 1, N)
        prom[i] = np.mean(x[:, i])
        var[i] = np.mean(np.power(x[:, i], 2)) - np.power(prom[i], 2)
        corr[i] = np.mean(x[:, i] * x[:, 0])
    p_eq, edges = np.histogram(x[:, -1], density=True, bins=20)

    for tray_n in range(N):
        tray_source.data[str(tray_n)] = x[tray_n]
    analysis_source.data['prom'] = prom
    analysis_source.data['var'] = var
    analysis_source.data['corr'] = corr
    hist_source.data['bottom'] = edges[:-1]
    hist_source.data['top'] = edges[1:]
    hist_source.data['p_eq'] = p_eq


# Plots
plot_tray = figure(title='Trayectories',
                   x_axis_label='t', y_axis_label='X(t)',
                   x_range=(0, TLIM), output_backend="svg")
for tray_n in range(N):
    plot_tray.line('t', str(tray_n), source=tray_source, alpha=min(1.0, 10 / N))
plot_tray.line('t', 'prom', source=analysis_source, color='red', legend='mean')

plot_hist = figure(title='Equilibrium distribution',
                   x_axis_label='P_eq(x)', y_axis_label='x', y_range=plot_tray.y_range,
                   output_backend="svg")
plot_hist.quad(right='p_eq', left=0, top='top', bottom='bottom', source=hist_source)
plot_hist.line(x='y', y='x', source=gauss_source, color='red',
               legend="mean: {0:.3f}  std: {1:.3f}".format(float(mean), float(std)))

plot_var = figure(title='Variance',
                  x_axis_label='t', y_axis_label='var(X(t))', x_range=plot_tray.x_range,
                  output_backend="svg")
plot_var.line('t', 'var', source=analysis_source)

plot_corr = figure(title='Correlation',
                   x_axis_label='t', y_axis_label='corr(X(t))', x_range=plot_tray.x_range,
                   output_backend="svg")
plot_corr.line('t', 'corr', source=analysis_source)

# Control panel
a_input = TextInput(title='a', value='1.0', placeholder='1.0')
d_input = TextInput(title='D', value='0.2', placeholder='0.2')
dt_input = TextInput(title='dt', value='0.01', placeholder='0.01')
xinit_input = TextInput(title='X(0)', value='1.0', placeholder='1.0')
update_button = Button(label='Update')

latex = LatexLabel(text="A(x)=-ax-x^3",
                   x=180, y=-80, x_units='screen', y_units='screen',
                   render_mode='css', text_font_size='16pt',
                   background_fill_color='#ffffff')
latex2 = LatexLabel(text="D(x)=D",
                    x=180, y=-150, x_units='screen', y_units='screen',
                    render_mode='css', text_font_size='16pt',
                    background_fill_color='#ffffff')
latex3 = LatexLabel(text="\\Xi(t)=A(X(t))dt+\\sqrt{D}\\sqrt{dt}\\mathcal{N}(0,1)",
                    x=180, y=-220, x_units='screen', y_units='screen',
                    render_mode='css', text_font_size='16pt',
                    background_fill_color='#ffffff')

plot_hist.add_layout(latex)
plot_hist.add_layout(latex2)
plot_hist.add_layout(latex3)

update_button.on_click(update)

controlpanel = column(widgetbox(a_input, width=20),
                      widgetbox(d_input, width=20),
                      widgetbox(dt_input, width=20),
                      widgetbox(xinit_input, width=20),
                      widgetbox(update_button, width=20))

markov_dynamic_tab = Panel(child=gridplot(children=[[plot_tray, plot_hist, plot_corr],
                                                    [plot_var, controlpanel, None]],
                                          plot_width=400,
                                          plot_height=400),
                           title='Markov dynamic')
