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

gitter_moegl = ['SC','FCC','BCC','Diamant']


####################################################################################################
# Funktionen zur Strukturbestimmung
####################################################################################################

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

def bragg(lambdaa, theta):
    return lambdaa / (2 * unp.sin(theta))

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
	radius_vorne, radius_rueck = np.genfromtxt('MetallRadius.txt', unpack = True) # das ist ein 2D array, der in der ersten Spalte die Radien weg von der Röntgenquelle hat und in der zweiten zur Röntgenquelle
	radius_rueck = umfang_kamera - radius_rueck
	radius_rueck = radius_rueck[::-1]
	radius = np.concatenate((radius_vorne, radius_rueck), axis = 0)

	return radius

####################################################################################################
# Und los
####################################################################################################

def main():

	print('\n#################### Analyse für Metall ####################\n')
	print('radius', radius())

	Xi2_best = ['SC', 20] # SC ist in dem Fall ein Platzhalter. Die 20 garantiert, dass ein jedes Xi zunächst kleiner ist.
	Xi_Fehler = ['SC', 20]

	daten = np.array([44.5, 64.5, 81.5, 98.5, 114.5, 134.5])
	daten = unp.uarray(daten, 1)
	daten = (daten * np.pi) / (360)
	
	# theta = theta_radiant(radius())
	theta =daten
	print('Theta mit Fehler: ', theta)
	netzebenenabstand = bragg(lambda_1, noms(theta))

	reflexe_SC = []
	reflexe_FCC = []
	reflexe_BCC = []
	reflexe_Diamant = []

	verhaeltnisse_temp = []

	for gitter in gitter_moegl:
		infos = Strukturamplitude(gitter = gitter)
		reflexe = np.array(infos[0])
		# print(gitter +': ',reflexe[np.argsort(infos[1])])
		m = infos[1]
		m = np.sort(m)

		verhaeltnisse = findStructure(m, netzebenenabstand)
		verhaeltnisse_temp = verhaeltnisse[0]
		# verhaeltnis_m = np.sort(verhaeltnisse[0])
		verhaeltnis_m = verhaeltnisse[0]
		if gitter == 'SC':
			reflexe_SC = verhaeltnis_m
		elif gitter == 'FCC':
			reflexe_FCC = verhaeltnis_m
		elif gitter == 'BCC':
			reflexe_BCC = verhaeltnis_m
		elif gitter == 'Diamant':
			reflexe_Diamant = verhaeltnis_m

		verhaeltnis_d = np.sort(verhaeltnisse[1])

		print('sqrt(m_i/m_1): ', verhaeltnis_m)
		print('d_1/d_i: ', verhaeltnis_d)

		print('Verhältnisse für die m Werte: ', verhaeltnis_m)
		print('Verhältnisse für die d Werte: ', verhaeltnis_d)
		print('Abweichung Xi^2 für die ' + gitter + ' Sturktur: ', abweichung(verhaeltnis_m, verhaeltnis_d))

		if abweichung(verhaeltnis_m, verhaeltnis_d) < Xi2_best[1]:
			Xi2_best = [gitter, abweichung(verhaeltnis_m, verhaeltnis_d)]

		
		infos = Strukturamplitude(gitter = gitter)
		reflexe = np.array(infos[0])
		m = infos[1]
		m = np.sort(m)
		netzebenenabstand_Fehler = bragg(lambda_1, theta)
		verhaeltnisse_Fehler = findStructure(m, netzebenenabstand_Fehler)
		verhaeltnis_m_Fehler = verhaeltnisse_Fehler[0]
		verhaeltnis_d_Fehler = np.sort(verhaeltnisse_Fehler[1])

		print('Abweichung Xi^2 Fehler für die ' + gitter + ' Sturktur: ', abweichung(verhaeltnis_m_Fehler, verhaeltnis_d_Fehler))

		if abweichung(verhaeltnis_m_Fehler, verhaeltnis_d_Fehler) < Xi_Fehler[1]:
			Xi_Fehler = [gitter, abweichung(verhaeltnis_m_Fehler, verhaeltnis_d_Fehler)]

	print('Struktur mit der kleinsten Abweichung: ', Xi2_best[0])
	print('Abweichung Xi^2: ', Xi2_best[1])
	print('Abweichung Xi^2 Fehler: ', Xi_Fehler[1])

	# m = Strukturamplitude(gitter = Xi2_best[0])[1]
	m = Strukturamplitude(gitter = 'BCC')[1]
	a = np.array(gitterkonstanteBragg(m, netzebenenabstand))

	print('Gitterkonstanten für ' + Xi2_best[0] + ' Struktur: ', a)

	####################################################################################################
	# Systematischer Fehler
	####################################################################################################

	# systematischer Fehler Absorption der Röntgenstrahlung

	DeltaA = (radius_rohr / (2 * Radius_kamera)) * (1 - Radius_kamera / Abstand_FP) * (np.cos(noms(theta))**2 / noms(theta)) * a

	a_mitFehler = unp.uarray(a, DeltaA)
	print('Gitterkonstanten mit Fehler durch die Absorption: ', a_mitFehler)

	####################################################################################################
	# curve_fit zeugs
	####################################################################################################


	# linearer fit für a gegen cos^2
	# params, cov = curve_fit(lin,np.cos(noms(theta))**2,noms(a_mitFehler), sigma = stds(a_mitFehler))
	params, cov = curve_fit(lin,np.cos(noms(theta))**2,noms(a_mitFehler))
	err = np.sqrt(np.diag(cov))
	a_extrp = ufloat(params[1], err[1])
	print('Extrapolierte Gitterkonstante: ', a_extrp)

	####################################################################################################
	# Plots
	####################################################################################################

	cos2Theta = unp.cos(theta)**2
	cos2Theta_fit = np.linspace(0, 1)

	####### ab hier der a gegen cos^2 Plot ########
	plt.errorbar(np.cos(noms(theta))**2, a, xerr=stds(cos2Theta), yerr=DeltaA, fmt='x', label = 'Daten')
	plt.plot(cos2Theta_fit, lin(cos2Theta_fit, *params), label = 'Fit')
	plt.legend(loc = 'best')
	plt.xlabel('cos$^2(\Theta)$')
	plt.ylabel('$a$ in Angtröm')
	# plt.xlim(0.6,1)
	# plt.ylim(4.8e-10,5.8e-10)
	plt.tight_layout()
	plt.grid()
	plt.savefig('Plots/Metall_Fit.pdf')
	plt.close()

	i = np.linspace(1,6,6)

	verhaeltnisse_daten = []
	for j in range(0,len(daten)):	
		verhaeltnisse_daten.append(unp.sin(daten[j]) / unp.sin(daten[0]))

	plt.plot(i, reflexe_SC, 'v', label = 'SC')
	plt.plot(i, reflexe_FCC, 's', label = 'FCC')
	plt.plot(i, reflexe_BCC, '^', label = 'BCC')
	plt.plot(i, reflexe_Diamant, 'o', label = 'Diamant')
	plt.plot(i, noms(verhaeltnisse_daten), '<', label = 'data')
	plt.xlabel('i')
	plt.xlim(0,7)
	plt.ylim(0,4)
	plt.ylabel(r'$\frac{m_i}{m_1}$')
	plt.legend(loc = 'best')
	plt.tight_layout()
	plt.grid()
	plt.savefig('Plots/verhaeltnisse.pdf')
	plt.close()

main()
