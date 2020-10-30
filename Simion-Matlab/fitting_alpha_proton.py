import numpy as np
from scipy.optimize import curve_fit
from scipy.special import gamma, factorial
import csv

with open('fitting_data_maxwell_SIMION_alpha.txt') as csvfile:
	readCSV = csv.reader(csvfile, delimiter=',')
	data = list(readCSV)
	velocities_maxwell = []
	f_maxwell = []
	for i in range (0,len(data[0])):
		velocities_maxwell.append(float(data[0][i]))
		f_maxwell.append(float(data[1][i]))

def gaussian(x, n, cen, vt):
	return (n/(np.pi**(3/2)*vt**3)) * np.exp(-(x-cen)**2/vt**2)

x = velocities_maxwell
y = f_maxwell
position = f_maxwell.index(max(f_maxwell))
init_vals = [1e6, x[position], 1e4]  # for [amp, cen, vt]
best_vals, covar = curve_fit(gaussian, x, y, p0=init_vals)
print('best_vals: {}'.format(best_vals))
perr = np.sqrt(np.diag(covar))
print(perr)

