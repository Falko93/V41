import math              # Matheumgebung für Funktionen und Konstanten
import numpy as np       # Matheumgebung für Funktionen und Konstanten
import matplotlib as mlp        #plotten
import matplotlib.pyplot as plt #plotten
import scipy.constants as const #fitten
from scipy.optimize import curve_fit
from uncertainties import ufloat, umath #ab hier Pakete für Fehlerrechnung.
from uncertainties import unumpy as unp
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
from scipy.optimize import curve_fit
from scipy.stats import sem
from scipy import integrate

###############################################################################
#Metall
###############################################################################

###############################################################################
#Strukturfaktoren berechnen
###############################################################################

# lambda_1 = 1.5417e-10
# Radius_k = 5.73 #in cm
# radius_rohr = 0.4e-3


# radius_m = np.array([4.5, 6.25, 8.2, 9.8, 11.5]) # in cm

# #Bogenlängenformel b=alpha*r*pi /180

# theta_4 = radius_m * 180 / (np.pi * Radius_k)
# theta = theta_4 * np.pi / (4 * 180) # in radianten

# print(' Theta=', theta)
# theta_sin = np.sin(theta)
# print('sin(\theta) = ', theta_sin)


# ######### Bragg Bedingung ##################


# def bragg(n, lambdaa, theta):
#     return n * lambdaa / (2 * np.sin(theta))

# netzebenenabstand_1 = bragg(1, lambda_1, theta[0])



#print('Netzebenenabstand',netzebenenabstand_1)

######### Streuamplitude ######################


############ Metall ################

### Nicht verschwindende Reflexe werden mit der Strukturamplitude berechnet. Es werden alle hkl Kombinationen von 100 bis 999 betrachtet wobei daruf geachtet wird, dass m = h^2 + k^2+ l^2 sich nicht doppelt.

def Strukturamplitude(f = 1, gitter = 'SC'): # gitter = 'SC', 'BCC', 'FCC', 'Diamant'
	print('Strukturamplitude berechnet für ' + gitter + '.')
	if gitter == 'SC':
		r = np.array([[0,0,0]])
	if gitter == 'BCC':
		r = np.array([[0,0,0],[.5,.5,.5]])
	if gitter == 'FCC':
		r = np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5]])
	if gitter == 'Diamant':
		r = np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5],[0.25,.25,.25],[.75,.75,.25],[.75,.25,.75],[.25,.75,.75]])

	m = []
	temp = 0
	S = 0
	reflexe = []
	for h in range(0,10):
		for k in range(0,h+1):
			for l in range(0,k+1):
				temp = h**2 + k**2 + l**2
				if temp not in m:
					S = 0
					for i in range(0,len(r)):
						S = S + f * np.exp(-2j * np.pi * (r[i][0] * h + r[i][1] * k + r[i][2] * l))
					m.append(h**2 + k**2 + l**2)
					if 0.001 < abs(S):
						reflexe.append([h, k, l])
				else:
					pass

	return reflexe

def main():
	print(Strukturamplitude(gitter = 'SC'))

main()



















