from params import *
import math
import numpy as np
from scipy import integrate
import matplotlib as mpl
import matplotlib.pyplot as plt

def c_reduce(t):
	if t < t_cont1 or t > t_cont2:
		return 0
	else:
	    return c_cont
	    
def i_iso(i, t):

	if t < t_iso1 or t > t_iso2:
		return 0
	else:
		num = p_sick * i
		if num > q_max:
			num = q_max
		return num
		
def i_home(i, t):
	if t < t_iso1 or t > t_iso2:
		return 0
	else:
		num = p_sick * i
		if num > q_max:
			return num - q_max
		else:
			return 0

def beta_i(t):
	# Contact rate at time t for infected patients
	return beta * (1.0 + a*math.cos(2.0 * math.pi * (t - a_max)/365.25))
	
def beta_p(t):
	# Contact rate at time t for prodromal patients
	return i_p * beta_i(t)
	
def sprime(s, p, i, t):
	# Force of infection at time t.
	return ((-s/N) * (1.0 - c_reduce(t)) * 
			((p * beta_p(t)) + 
			(beta_i(t) *(i - i_iso(i, t) - (i_home(i, t) * c_home))) +
			psi))

def derivs(y, t):
	s = y[0]
	e = sum(y[1:17])
	p = sum(y[17:33])
	i = sum(y[33:49])
	yp = np.zeros(51)
	yp[0] = sprime(s, p, i, t)
	yp[1] = -yp[0] - (epsilon * y[1])
	for i in range(2,17):
		yp[i] = epsilon * (y[i-1] - y[i])
	yp[17] = (epsilon * y[16]) - (phi * y[17])
	for i in range(18,33):
		yp[i] = phi * (y[i-1] - y[i])
	yp[33] = (phi * y[32]) - (gamma * y[33])
	for i in range(34,49):
		yp[i] = gamma * (y[i-1] - y[i])
	yp[49] = gamma * (1.0 - p_sick * p_death) * y[48]
	yp[50] = gamma * p_sick * p_death * y[48]

	return yp


y = np.zeros(51)
y[0] = 5.0e6 - 100
y[1] = 100.0
t = np.arange(0, 501)
xy_t = integrate.odeint(derivs, y, t)
print(xy_t.shape)

s = xy_t[:, 0]
e = xy_t[:, 1:17].sum(axis=1)
p = xy_t[:, 17:33].sum(axis=1)
i = xy_t[:, 33:49].sum(axis=1)
r = xy_t[:, 49]
d = xy_t[:, 50]

print("Max deaths:", d.max())
# Plot the data on six separate curves
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, s/1000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, e/1000, 'y', alpha=0.5, lw=2, label='Exposed')
ax.plot(t, p/1000, 'y', alpha=0.5, lw=2, label='Prodromal')
ax.plot(t, i/1000, 'r', alpha=0.5, lw=2, label='Infected')
ax.plot(t, r/1000, 'g', alpha=0.5, lw=2, label='Recovered')
ax.plot(t, d/1000, 'b', alpha=0.5, lw=2, label='Dead')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number (1000s)')
ax.set_ylim(0,6000)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()

