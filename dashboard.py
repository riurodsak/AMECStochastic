#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# dashboard.py
# Created by Riu on 18/12/2018 for the project AMECStochastic
# Bokeh dashboard groups all simulation visualizations into one webpage.


from bokeh.io import curdoc
from bokeh.models.widgets import Tabs

# from markov import markov_tab
from markov_dynamic import markov_dynamic_tab
# from gillespie import gillespie_tab


tabs = Tabs(tabs=[markov_dynamic_tab])  # add tabs as required

curdoc().add_root(tabs)
curdoc().title = "AMEC Stochastic"
