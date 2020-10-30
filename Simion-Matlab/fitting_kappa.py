import numpy as np
from scipy.optimize import curve_fit
from scipy.special import gamma, factorial
import csv

with open('fitting_data_kappa_SIMION.txt') as csvfile:
	readCSV = csv.reader(csvfile, delimiter=',')
	data = list(readCSV)
	velocities_kappa = []
	f_kappa = []
	for i in range (0,len(data[0])):
		velocities_kappa.append(float(data[0][i]))
		f_kappa.append(float(data[1][i]))

def kappa(x, n, kappa, cen, vt):
    	return (n/vt**3) * (2/(np.pi*(2*kappa-3)))**(3/2) * (gamma(kappa+1)/gamma(kappa-1/2))  * (1+(2/(2*kappa-3))*((x-cen) **2)/(vt**2))**(-kappa-1)

x = velocities_kappa
y = f_kappa
position = f_kappa.index(max(f_kappa))
init_vals = [1e6, 1.7, x[position], 1e4]	# for [density, kappa, cen, vt]
best_vals, covar = curve_fit(kappa, x, y, p0=init_vals)
print('best_vals: {}'.format(best_vals))
perr = np.sqrt(np.diag(covar))
print(perr)

