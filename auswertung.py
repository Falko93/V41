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

lambda_1 = 1.5417e-10
Radius_k = 5.73 #in cm
radius_rohr = 0.4e-3


radius_m = np.array([4.5, 6.25, 8.2, 9.8, 11.5]) # in cm

#Bogenlängenformel b=alpha*r*pi /180

theta_4 = radius_m*180/(np.pi*Radius_k)
theta = theta_4*np.pi/(4*180) # in radianten

print(' Theta=', theta)
theta_sin = np.sin(theta)
print('sin(\theta) = ', theta_sin)


######### Bragg Bedingung ##################


def bragg(n, lambdaa, theta):
    return n*lambdaa/(2*np.sin(theta))

netzebenenabstand_1 = bragg(1, lambda_1, theta[0])



#print('Netzebenenabstand',netzebenenabstand_1)

######### Streuamplitude ######################
