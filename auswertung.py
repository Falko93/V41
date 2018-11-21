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


lambda_1 = 1.5417e-10
Radius_k = 5.73 #in cm
radius_rohr = 0.4e-3


radius_m = np.array([4.5, 6.25, 8.2, 9.8, 11.5]) # Radius der Kreibögen in cm für Metall

### Nicht verschwindende Reflexe werden mit der Strukturamplitude berechnet. Es werden alle hkl Kombinationen von 100 bis 999 betrachtet wobei daruf geachtet wird, dass m = h^2 + k^2+ l^2 sich nicht doppelt.

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

######### Bragg Bedingung ##################


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

def lin(x, a, b):
	return a * x + b


def main():
	# print(Strukturamplitude(gitter = 'SC'))
	# print(Strukturamplitude(gitter = 'FCC'))
	# print(Strukturamplitude(gitter = 'BCC'))
	gitter_moegl = ['SC','FCC','BCC','Diamant']
	Xi2_best = ['SC', 10]

	######################## Metall ##########################
	#Bogenlängenformel b=alpha*r*pi /180
	theta_4 = radius_m * 180 / (np.pi * Radius_k)
	theta = theta_4 * np.pi / (4 * 180) # in radianten
	theta_sin = np.sin(theta)
	netzebenenabstand_metall = bragg(lambda_1, theta)

	for gitter in gitter_moegl:
		infos = Strukturamplitude(gitter = gitter)
		reflexe = infos[0]
		m = infos[1]
	
		verhaeltnisse = findStructure(m, netzebenenabstand_metall)
		verhaeltnis_m = verhaeltnisse[0]
		verhaeltnis_d = verhaeltnisse[1]
	
		print('Verhältnisse für die m Werte: ', verhaeltnis_m)
		print('Verhältnisse für die d Werte: ', verhaeltnis_d)
		print('Abweichung Xi^2 für die ' + gitter + ' Sturktur: ', abweichung(verhaeltnis_m, verhaeltnis_d))

		if abweichung(verhaeltnis_m, verhaeltnis_d) < Xi2_best[1]:
			Xi2_best = [gitter, abweichung(verhaeltnis_m, verhaeltnis_d)]

	print('Struktur mit der kleinsten Abweichung: ', Xi2_best[0])
	print('Abweichung Xi^2: ', Xi2_best[1])
		
	m = Strukturamplitude(gitter = Xi2_best[0])[1]
	a = np.array(gitterkonstanteBragg(m, netzebenenabstand_metall))

	print('Gitterkonstanten für ' + Xi2_best[0] + ' Struktur: ', a)

	params, cov = curve_fit(lin,np.cos(theta)**2,a)
	err = np.sqrt(np.diag(cov))

	a_extrp = ufloat(params[1], err[1])

	print('Extrapolierte Gitterkonstante: ', a_extrp)

	cos2Theta_fit = np.linspace(0, 1)

	plt.plot(np.cos(theta)**2, a, 'x', label = 'Daten')
	plt.plot(cos2Theta_fit, lin(cos2Theta_fit, *params), label = 'Fit')
	plt.legend(loc = 'best')
	plt.xlabel('cos$^2(\Theta)$')
	plt.ylabel('$a$ in Angtröm')
	plt.xlim(0.6,1)
	plt.tight_layout()
	plt.grid()
	plt.savefig('Plots/Metall_Fit.pdf')


main()



















