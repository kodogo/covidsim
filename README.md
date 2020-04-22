This repo contains a bit of Python code attempting to reproduce some of the results in the Otago modelling of the 2020 COVID-19 pandemic in NZ.

The model itself seems to be the same as v1.0 of the one here: https://covidsim.eu

It is a pretty standard SIR type model with the following compartments:

S - Susceptible
E - Exposed
P - Pre-symptomatic
I - Infectious
R - Recovered or removed

Details are in the two Otago papers in the doc folder. The version of the model currently (23/4/2020) on the site is v1.1 which differs slightly. Details are in a paper in the doc folder - I have not bothered to change the program to reflect the changes.

One slight complication over simple SIR type models is that the E, P, and I phases are modelled as an Erlang distribution using 16 stages -  a couple of papers in the doc section explain this sort of thing. In practice it doesn’t seem to make a big difference, especially given the uncertainties in the underlying quantities.

The code itself is in /src. params.py has the assumptions, and covidsim.py just plugs these into the psych ode solver. Not fully tested, but the results were close enough to the published ones to convince me it was about right.