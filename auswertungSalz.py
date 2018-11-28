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

####################################################################################################
# alle Konstanten werden hier definiert 
####################################################################################################

lambda_1 = 1.5417e-10 # in m
Radius_kamera = 0.0573 #in m
radius_rohr = 0.4e-3
Abstand_FP = 0.130 # in m -> Abstand Fokus Probe
umfang_kamera = 2*np.pi*Radius_kamera # Umfang
Fehler_radius = 0.001 # in m -> der Fehler der auf jeden Radius Wert durch das Lineal drauf gerechnet werden muss

####################################################################################################
# alle globalen Variablen werden hier definiert
####################################################################################################

gitter_moegl = ['Zinkblende', 'Steinsalz', 'Caesiumchlorid', 'Fluorit']


####################################################################################################
# Funktionen zur Strukturbestimmung
####################################################################################################

def bragg(lambdaa, theta):
    return lambdaa / (2 * np.sin(theta))

def findStructure(m, d):
	verhaeltnis_m = []
	verhaeltnis_d = []

	for i in range(0,len(d)):
		verhaeltnis_m.append(np.sqrt(m[i] / m[0]))
		verhaeltnis_d.append(d[0] / d[i])

	verhaeltnisse = np.array([verhaeltnis_m, verhaeltnis_d])

	return verhaeltnisse

def abweichung(verhaeltnis_m, verhaeltnis_d):
	Xi2 = sum((verhaeltnis_d - verhaeltnis_m)**2 / verhaeltnis_m)

	return Xi2

def gitterkonstanteBragg(m, d):
	a = []

	for i in range(0,len(d)):
		a.append(np.sqrt(m[i]) * d[i])

	return a

def Strukturamplitude(f = 1, gitter = 'SC'): # gitter = 'SC', 'BCC', 'FCC', 'Diamant'
	print('Strukturamplitude wird berechnet für ' + gitter + '.')
	if gitter == 'SC':
		r = np.array([[0,0,0]])
	elif gitter == 'BCC':
		r = np.array([[0,0,0],[.5,.5,.5]])
	elif gitter == 'FCC':
		r = np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5]])
	elif gitter == 'Diamant':
		r = np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5],[0.25,.25,.25],[.75,.75,.25],[.75,.25,.75],[.25,.75,.75]])

	m = []
	m_reflexe = []
	temp = 0
	S = 0
	reflexe = []
	for h in range(0,10):
		for k in range(0,h+1):
			for l in range(0,k+1):
				temp = h**2 + k**2 + l**2
				if temp not in m and temp != 0:
					S = 0
					for i in range(0,len(r)):
						S = S + f * np.exp(-2j * np.pi * (r[i][0] * h + r[i][1] * k + r[i][2] * l))
					m.append(h**2 + k**2 + l**2)
					if 0.001 < abs(S): # es wird nie ganz Null weil Computaional Physics
						reflexe.append([h, k, l])
						m_reflexe.append(h**2 + k**2 + l**2)
				else:
					pass
	infos = [reflexe, m_reflexe]

	return infos

####################################################################################################
# Allgemeine Funktionen
####################################################################################################

def lin(x, a, b):
	return a * x + b

def theta_radiant(radius): # diese Funktion gibt einen uarray raus in dem Theta steht; umgerechnet aus dem gemessenen Radius
	radius_unp = unp.uarray(radius, Fehler_radius)
	theta = radius_unp / (4*Radius_kamera)
	
	return np.sort(theta) # in radianten + Fehler

def radius():
	radius_vorne, radius_rueck = np.genfromtxt('SalzRadius.txt', unpack = True) # das ist ein 2D array, der in der ersten Spalte die Radien weg von der Röntgenquelle hat und in der zweiten zur Röntgenquelle
	radius_vorne = radius_vorne[:-3]
	radius_rueck = umfang_kamera - radius_rueck
	radius = np.concatenate((radius_vorne, radius_rueck), axis = 0)

	return radius

def main():

	print('\n#################### Analyse für Salz ####################\n')

	Xi2_best = ['SC', 20] # SC ist in dem Fall ein Platzhalter. Die 20 garantiert, dass ein jedes Xi zunächst kleiner ist.

	theta = theta_radiant(radius())
	print('Theta mit Fehler: ', theta)
	netzebenenabstand = bragg(lambda_1, noms(theta))

	sThetaLambda = unp.sin(theta/lambda_1)

	print(sThetaLambda)

main()




















