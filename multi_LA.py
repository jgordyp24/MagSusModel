import numpy as np
import matplotlib.pyplot as plt
import lightFunctions as lf
from scipy.optimize import curve_fit
import os


#Determine fit params for magnetic field######
sensor_pos = 0.015115 #position of sensor
sensor_w = 0.001 #Estimate sampling width at center distance
res = 1000 #data resolution

z = np.linspace(sensor_pos-sensor_w, sensor_pos+sensor_w, res)
B_r = 1.48 #Tesla
#magnet dimensions in meters
L = 0.0254 #1 inch in m
W = 0.0254 #1 inch in m
T = 0.009525 #3/8th inch in m (thickness)

B_s = lf.B_field(z, B_r, L, W, T)

(popt_B, pcov_B) = curve_fit(lf.B_LinFit, z, B_s, p0=None, sigma=None, bounds=(-np.inf,np.inf))

A, b = popt_B

print("A = {}\nb = {}".format(A, b))
########

### initial guess for fit function
eps = -0.05 #-0.75
S1 = -176#-0.001#-0.056 #-0.1
S2 = -0.00729#-0.001#-1189485 #-0.2
omega = -1 #-1.12
guesses = np.array([eps,S1,S2,omega])
###
density = 5240

# fileNames = np.array(["co1_13.txt","co1_23.txt","co1_33.txt","co2_13.txt",
# 							"co2_33.txt","co3_13.txt","co3_33.txt","co4_33.txt",])

path = "lightData/slitData/3_16inS"
fileNames = os.listdir(path)


file = open("MagSusVals.txt", "w")
file.write("Data for 3/8 inch magnets\n")
file.write("Name\tX\t\teps\t\tS1\t\tS2\t\tomega\n")

X_array = np.zeros(len(fileNames))
S1_array = np.zeros(len(X_array))
S2_array = np.zeros(len(X_array))

# plt.figure()
for i, name in enumerate(fileNames):
	#Determine fit params for light curves#########
	t_raw, T_raw = np.loadtxt(path+"/{}".format(name), skiprows=3, delimiter=None,unpack=True)

	T_RAWlog = lf.matchEXP(T_raw)
	t_T, T_logT = lf.trunc(15, 300, t_raw, T_RAWlog) 
	t_F, T_logF = lf.dataAdj(t_T, T_logT)
	###

	###Curve fitting trasm function with EXP data##
	(popt_T, pcov_T) = curve_fit(lf.transm, t_F, T_logF, p0=guesses, maxfev=10000, bounds=(-np.inf,np.inf))
	#variables are redefined here. Ensure this doesnt cause issues later
	eps, S1, S2, omega = popt_T

	alpha = lf.alpha(S1, S2)
	beta = lf.beta(S1, S2)
	X = lf.mag_sus(density, A, S1, S2)

	X_array[i] = X
	S1_array[i] = S1
	S2_array[i] = S2

	file.write(name+"\t{:.2}\t\t{:.2}\t\t{:.2}\t\t{:.2}\t\t{:.2}\n".format(X, *popt_T))

	time = np.linspace(0, 285, res)

	plt.figure()
	plt.plot(t_F, T_logF)
	plt.plot(time, lf.transm(time, *popt_T), label=name+"\nS1: {:.4}\nS2: {:.4}\nX: {:.3}".format(S1, S2, X))
	plt.legend()
	plt.xlabel("time (s)")
	plt.ylabel("log(1/T)")
	plt.savefig("plots/{}.png".format(name), dpi=80)
	# plt.show()
	plt.close()

print("Avg. S1*S2 Vals: ", np.mean(S1_array*S2_array))

print("X Vals: ", X_array)
avgX = np.mean(X_array)

file.write("AVERAGE Mag Sus:\n")
file.write("{}".format(avgX))

file.close()

print("Average Sus.: X = {:.3}\n".format(avgX))





