# AMECStochastic
Simulation of two stochastic algorithms implemented in Python 3: Gillespie algorithm and Markov process. Interactive web-based visualization in Bokeh.

## Usage
To run, call `Scripts/bokeh-script.py` in your Python environment, with parameters `serve --show \<path-to-repository>\AMECStochastic\dashboard.py --args myargs`.

`dashboard` can be edited to include the desired simulation tabs. Each simulation is contained within it's tab, all in a single webpage.

## Dependencies
* Python 3
* [Bokeh](https://bokeh.pydata.org/en/latest/)

## Examples
Example of the Gillespie algorithm:
![Gillespie algorithm](/Images/gillespie.png)

Example of the Markov process (non-interactive).
![Markov process](/Images/markov_static.png)

Example of the interactive Markov process. Parameters can be changed to update the simulation conditions.
![Interactive Markov process](/Images/markov_dynamic.png)

