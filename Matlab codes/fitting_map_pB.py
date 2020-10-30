import numpy as np
from scipy.optimize import curve_fit
from scipy.special import gamma, factorial
import csv


list_a = ["0.01", "0.04", "0.07", "0.1", "0.13", "0.16", "0.19", "0.22", "0.25", "0.28", "0.31", "0.34", "0.37", "0.4", "0.43", "0.46", "0.49", "0.52"]
list_l = ["0.02", "0.015", "0.01", "0.0075", "0.005", "0.002", "0.0015", "0.001", "0.00075", "0.0005", "0.0002", "0.0001", "7.5e-05", "5e-05", "2.5e-05"]

density_maxwell = np.zeros((len(list_a), len(list_l)))
thermal_speed_maxwell = np.zeros((len(list_a), len(list_l)))
bulk_speed_maxwell = np.zeros((len(list_a), len(list_l)))

density_kappa = np.zeros((len(list_a), len(list_l)))
thermal_speed_kappa = np.zeros((len(list_a), len(list_l)))
bulk_speed_kappa = np.zeros((len(list_a), len(list_l)))
kappa_param = np.zeros((len(list_a), len(list_l)))

def open_maxwell(file):

	with open(file) as csvfile:
		readCSV = csv.reader(csvfile, delimiter=',')
		data = list(readCSV)
		velocities_maxwell = []
		f_maxwell = []
		for i in range (0,len(data[0])):
			velocities_maxwell.append(float(data[0][i]))
			f_maxwell.append(float(data[1][i]))
	return velocities_maxwell, f_maxwell

def gaussian(x, n, cen, vt):
	return (n/(np.pi**(3/2)*vt**3)) * np.exp(-(x-cen)**2/vt**2)

def open_kappa(file):

	with open(file) as csvfile:
		readCSV = csv.reader(csvfile, delimiter=',')
		data = list(readCSV)
		velocities_kappa = []
		f_kappa = []
		for i in range (0,len(data[0])):
			velocities_kappa.append(float(data[0][i]))
			f_kappa.append(float(data[1][i]))
	return velocities_kappa, f_kappa

def kappa(x, n, kappa, cen, vt):
    	return (n/vt**3) * (2/(np.pi*(2*kappa-3)))**(3/2) * (gamma(kappa+1)/gamma(kappa-1/2))  * (1+(2/(2*kappa-3))*((x-cen) **2)/(vt**2))**(-kappa-1)

for i in range(0,len(list_a)):
	for j in range(0,len(list_l)):
		file_maxwell = 	"fitting_data_maxwell_error_v_B_" + list_a[i] + "_l_" + list_l[j] + ".txt"
		file_kappa = 	"fitting_data_kappa_error_v_B_" + list_a[i] + "_l_" + list_l[j] + ".txt" 
		x_maxwell, y_maxwell = open_maxwell(file_maxwell)
		x_kappa, y_kappa = open_kappa(file_kappa)
		
		if (len(y_kappa) > 4 and len(y_maxwell) > 0):
	
			position = y_maxwell.index(max(y_maxwell))
			init_vals = [1e6, x_maxwell[position], 1e4]  # for [n, cen, vt]
			try:
				best_vals, covar = curve_fit(gaussian, x_maxwell, y_maxwell, p0=init_vals)
				density_maxwell[i][j] = best_vals.item(0)
				bulk_speed_maxwell[i][j] = best_vals.item(1)
				thermal_speed_maxwell[i][j] = best_vals.item(2)
			except (RuntimeError):
				print("Runtime Error")
				pass

			position = y_kappa.index(max(y_kappa))
			init_vals = [1e6, 3, x_kappa[position], 1e4]  # for [amp, kappa, cen, vt]
			try:
				best_vals, covar = curve_fit(kappa, x_kappa, y_kappa, p0=init_vals)
				density_kappa[i][j] = best_vals.item(0)
				kappa_param[i][j] = best_vals.item(1)
				bulk_speed_kappa[i][j] = best_vals.item(2)
				thermal_speed_kappa[i][j] = best_vals.item(3)
			except(RuntimeError):
				print("Runtime Error")
				pass
			
			init_vals = []

			print(best_vals.item(1))
			print(list_a[i], list_l[j])

print(density_maxwell)
print(" ")
print(bulk_speed_maxwell)
print(" ")
print(thermal_speed_maxwell)
print(" ")

print(density_kappa)
print(" ")
print(kappa_param)
print(" ")
print(bulk_speed_kappa)
print(" ")
print(thermal_speed_kappa)


np.savetxt("density_maxwell.csv", density_maxwell, delimiter=",")
np.savetxt("bulk_speed_maxwell.csv", bulk_speed_maxwell, delimiter=",")
np.savetxt("thermal_speed_maxwell.csv", thermal_speed_maxwell, delimiter=",")
np.savetxt("density_kappa.csv", density_kappa, delimiter=",")
np.savetxt("kappa_param.csv", kappa_param, delimiter=",")
np.savetxt("bulk_speed_kappa.csv", bulk_speed_kappa, delimiter=",")
np.savetxt("thermal_speed_kappa.csv", thermal_speed_kappa, delimiter=",")










