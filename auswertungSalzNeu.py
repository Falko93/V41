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

gitter_moegl = ['Zinkblende','Steinsalz','Caesiumchlorid','Fluorit']


####################################################################################################
# Funktionen zur Strukturbestimmung
####################################################################################################

def Strukturamplitude_Salz(f1 = 1, f2 = 1, gitter = 'Zinkblende'): # gitter = 'SC', 'BCC', 'FCC', 'Diamant' ; bei Fluorit f1=nicht Fluorit und f2=Fluorit 
	print('Strukturamplitude wird berechnet für ' + gitter + '.')
	if gitter == 'Zinkblende':
		r = np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5],[0.25,.25,.25],[.75,.75,.25],[.75,.25,.75],[.25,.75,.75]])
	elif gitter == 'Steinsalz':
		r = np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5],[0.5,.5,.5],[1.,1.,.5],[1.,.5,1.],[.5,1.,1.]])
	elif gitter == 'Caesiumchlorid':
		r = np.array([[0,0,0],[0.5,.5,.5]])
	elif gitter == 'Fluorit':
		r = np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5],[0.25,.25,.25],[.75,.75,.25],[.75,.25,.75],[.25,.75,.75],[0.75,.75,.75],[.25,.25,.75],[.25,.75,.25],[.75,.25,.25]])	

	m = []
	m_reflexe = []
	temp = 0
	S = 0
	reflexe = []
	if gitter == 'Fluorit': #bei Fluorit gibt es 2 mal ein Fluorti Atome und damit 3 Atome in der Struktur insgesamt; nicht 2
		for h in range(0,10):
			for k in range(0,h+1):
				for l in range(0,k+1):
					temp = h**2 + k**2 + l**2
					if temp not in m and temp != 0:
						S = 0
						for i in range(0,4):
							S = S + f1 * np.exp(-2j * np.pi * (r[i][0] * h + r[i][1] * k + r[i][2] * l))
						for i in range(4,12):
							S = S + f2 * np.exp(-2j * np.pi * (r[i][0] * h + r[i][1] * k + r[i][2] * l))
						m.append(h**2 + k**2 + l**2)
						if 0.001 < abs(S): # es wird nie ganz Null weil Computaional Physics
							reflexe.append([h, k, l])
							m_reflexe.append(h**2 + k**2 + l**2)
					else:
						pass
		infos = [reflexe, m_reflexe]
	
		return infos
	elif gitter == 'Caesiumchlorid':
		for h in range(0,10):
			for k in range(0,h+1):
				for l in range(0,k+1):
					temp = h**2 + k**2 + l**2
					if temp not in m and temp != 0:
						S = 0
						for i in range(0,1):
							S = S + f1 * np.exp(-2j * np.pi * (r[i][0] * h + r[i][1] * k + r[i][2] * l))
						for i in range(1,2):
							S = S + f2 * np.exp(-2j * np.pi * (r[i][0] * h + r[i][1] * k + r[i][2] * l))
						m.append(h**2 + k**2 + l**2)
						if 0.001 < abs(S): # es wird nie ganz Null weil Computaional Physics
							reflexe.append([h, k, l])
							m_reflexe.append(h**2 + k**2 + l**2)
					else:
						pass
		infos = [reflexe, m_reflexe]
	
		return infos
	else:
		for h in range(0,10):
			for k in range(0,h+1):
				for l in range(0,k+1):
					temp = h**2 + k**2 + l**2
					if temp not in m and temp != 0:
						S = 0
						for i in range(0,4):
							S = S + f1 * np.exp(-2j * np.pi * (r[i][0] * h + r[i][1] * k + r[i][2] * l))
						for i in range(4,8):
							S = S + f2 * np.exp(-2j * np.pi * (r[i][0] * h + r[i][1] * k + r[i][2] * l))
						m.append(h**2 + k**2 + l**2)
						if 0.001 < abs(S): # es wird nie ganz Null weil Computaional Physics
							reflexe.append([h, k, l])
							m_reflexe.append(h**2 + k**2 + l**2)
					else:
						pass
		infos = [reflexe, m_reflexe]
	
		return infos

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

	print('\n#################### Analyse für Salz ####################\n')

	Xi2_best = ['SC', 20, 'CuF'] # SC ist in dem Fall ein Platzhalter. Die 20 garantiert, dass ein jedes Xi zunächst kleiner ist.

	# daten1, daten2 = np.genfromtxt('SalzRadius.txt', unpack=True)
	# daten = np.array([23, 34.50, 41, 43, 51, 57, 58.5, 65, 70, 78, 83, 84.5, 90,95.5, 103, 108,110,116, 123, 133, 116 +43])
	daten = np.array([23, 34.50, 43, 51, 57, 65, 78, 83, 90,95.5, 103, 108, 123, 133, 116 +43])
	# daten = np.array([23, 34.50, 43, 51, 57, 65, 78, 83,95.5, 103, 108, 123, 133, 116 +43]) #9. raus

	# maske = [True,True,False,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True]
	# daten = daten[maske]
	daten = daten + 4.5
	# daten = daten[daten != 0]
	daten = (daten * np.pi) / (360)

	print('Daten: ', daten)

	# theta = theta_radiant(radius())
	theta = daten

	print('Theta mit Fehler: ', theta)
	netzebenenabstand = bragg(lambda_1, noms(theta))

	reflexe_Zinkblende = []
	reflexe_Steinsalz = []
	reflexe_Caesiumchlorid = []
	reflexe_Fluorit = []

	verhaeltnisse_temp = []

	for gitter in gitter_moegl:
		if gitter == 'Zinkblende':
			############## CuF ##################
			argg = 'CuF'
			print('Berechnung '+gitter+ ' für ' + argg + ':')
			infos = Strukturamplitude_Salz(28, 10, gitter = gitter)
			reflexe = np.array(infos[0])
			# print(gitter +': ',reflexe[np.argsort(infos[1])])
			m = infos[1]
			m = np.sort(m)
	
			verhaeltnisse = findStructure(m, netzebenenabstand)
			verhaeltnisse_temp = verhaeltnisse[0]
			# verhaeltnis_m = np.sort(verhaeltnisse[0])
			verhaeltnis_m = verhaeltnisse[0]
			if gitter == 'Zinkblende':
				reflexe_Zinkblende = verhaeltnis_m
			elif gitter == 'Steinsalz':
				reflexe_Steinsalz = verhaeltnis_m
			elif gitter == 'Caesiumchlorid':
				reflexe_Caesiumchlorid = verhaeltnis_m
			elif gitter == 'Fluorit':
				reflexe_Fluorit = verhaeltnis_m
	
			verhaeltnis_d = np.sort(verhaeltnisse[1])
	
			print('sqrt(m_i/m_1): ', verhaeltnis_m)
			print('d_1/d_i: ', verhaeltnis_d)
	
			print('Verhältnisse für die m Werte: ', verhaeltnis_m)
			print('Verhältnisse für die d Werte: ', verhaeltnis_d)
			print('Abweichung Xi^2 für die ' + gitter + ' Struktur für ' + argg + ': ', abweichung(verhaeltnis_m, verhaeltnis_d))
	
			if abweichung(verhaeltnis_m, verhaeltnis_d) < Xi2_best[1]:
				Xi2_best = [gitter, abweichung(verhaeltnis_m, verhaeltnis_d), argg]

			############## CuCl ##################

			argg = 'CuCl'
			print('Berechnung '+gitter+ ' für ' + argg + ':')
			infos = Strukturamplitude_Salz(28, 18, gitter = gitter)
			reflexe = np.array(infos[0])
			# print(gitter +': ',reflexe[np.argsort(infos[1])])
			m = infos[1]
			m = np.sort(m)
	
			verhaeltnisse = findStructure(m, netzebenenabstand)
			verhaeltnisse_temp = verhaeltnisse[0]
			# verhaeltnis_m = np.sort(verhaeltnisse[0])
			verhaeltnis_m = verhaeltnisse[0]
			if gitter == 'Zinkblende':
				reflexe_Zinkblende = verhaeltnis_m
			elif gitter == 'Steinsalz':
				reflexe_Steinsalz = verhaeltnis_m
			elif gitter == 'Caesiumchlorid':
				reflexe_Caesiumchlorid = verhaeltnis_m
			elif gitter == 'Fluorit':
				reflexe_Fluorit = verhaeltnis_m
	
			verhaeltnis_d = np.sort(verhaeltnisse[1])
	
			print('sqrt(m_i/m_1): ', verhaeltnis_m)
			print('d_1/d_i: ', verhaeltnis_d)
	
			print('Verhältnisse für die m Werte: ', verhaeltnis_m)
			print('Verhältnisse für die d Werte: ', verhaeltnis_d)
			print('Abweichung Xi^2 für die ' + gitter + ' Struktur für ' + argg + ': ', abweichung(verhaeltnis_m, verhaeltnis_d))
	
			if abweichung(verhaeltnis_m, verhaeltnis_d) < Xi2_best[1]:
				Xi2_best = [gitter, abweichung(verhaeltnis_m, verhaeltnis_d), argg]
		elif gitter == 'Steinsalz':
			############## KF ##################
			argg = 'KF'
			print('Berechnung '+gitter+ ' für ' + argg + ':')
			infos = Strukturamplitude_Salz(18, 10, gitter = gitter)
			reflexe = np.array(infos[0])
			# print(gitter +': ',reflexe[np.argsort(infos[1])])
			m = infos[1]
			m = np.sort(m)
	
			verhaeltnisse = findStructure(m, netzebenenabstand)
			verhaeltnisse_temp = verhaeltnisse[0]
			# verhaeltnis_m = np.sort(verhaeltnisse[0])
			verhaeltnis_m = verhaeltnisse[0]
			if gitter == 'Zinkblende':
				reflexe_Zinkblende = verhaeltnis_m
			elif gitter == 'Steinsalz':
				reflexe_Steinsalz = verhaeltnis_m
			elif gitter == 'Caesiumchlorid':
				reflexe_Caesiumchlorid = verhaeltnis_m
			elif gitter == 'Fluorit':
				reflexe_Fluorit = verhaeltnis_m
	
			verhaeltnis_d = np.sort(verhaeltnisse[1])
	
			print('sqrt(m_i/m_1): ', verhaeltnis_m)
			print('d_1/d_i: ', verhaeltnis_d)
	
			print('Verhältnisse für die m Werte: ', verhaeltnis_m)
			print('Verhältnisse für die d Werte: ', verhaeltnis_d)
			print('Abweichung Xi^2 für die ' + gitter + ' Struktur für ' + argg + ': ', abweichung(verhaeltnis_m, verhaeltnis_d))
	
			if abweichung(verhaeltnis_m, verhaeltnis_d) < Xi2_best[1]:
				Xi2_best = [gitter, abweichung(verhaeltnis_m, verhaeltnis_d), argg]

			############## KCl ##################

			argg = 'KCl'
			print('Berechnung '+gitter+ ' für ' + argg + ':')
			infos = Strukturamplitude_Salz(18, 18, gitter = gitter)
			reflexe = np.array(infos[0])
			# print(gitter +': ',reflexe[np.argsort(infos[1])])
			m = infos[1]
			m = np.sort(m)
	
			verhaeltnisse = findStructure(m, netzebenenabstand)
			verhaeltnisse_temp = verhaeltnisse[0]
			# verhaeltnis_m = np.sort(verhaeltnisse[0])
			verhaeltnis_m = verhaeltnisse[0]
			if gitter == 'Zinkblende':
				reflexe_Zinkblende = verhaeltnis_m
			elif gitter == 'Steinsalz':
				reflexe_Steinsalz = verhaeltnis_m
			elif gitter == 'Caesiumchlorid':
				reflexe_Caesiumchlorid = verhaeltnis_m
			elif gitter == 'Fluorit':
				reflexe_Fluorit = verhaeltnis_m
	
			verhaeltnis_d = np.sort(verhaeltnisse[1])
	
			print('sqrt(m_i/m_1): ', verhaeltnis_m)
			print('d_1/d_i: ', verhaeltnis_d)
	
			print('Verhältnisse für die m Werte: ', verhaeltnis_m)
			print('Verhältnisse für die d Werte: ', verhaeltnis_d)
			print('Abweichung Xi^2 für die ' + gitter + ' Struktur für ' + argg + ': ', abweichung(verhaeltnis_m, verhaeltnis_d))
	
			if abweichung(verhaeltnis_m, verhaeltnis_d) < Xi2_best[1]:
				Xi2_best = [gitter, abweichung(verhaeltnis_m, verhaeltnis_d), argg]
		elif gitter == 'Caesiumchlorid':
			############## CsJ ##################
			argg = 'CsJ'
			print('Berechnung '+gitter+ ' für ' + argg + ':')
			infos = Strukturamplitude_Salz(54, 54, gitter = gitter)
			reflexe = np.array(infos[0])
			# print(gitter +': ',reflexe[np.argsort(infos[1])])
			m = infos[1]
			m = np.sort(m)
	
			verhaeltnisse = findStructure(m, netzebenenabstand)
			verhaeltnisse_temp = verhaeltnisse[0]
			# verhaeltnis_m = np.sort(verhaeltnisse[0])
			verhaeltnis_m = verhaeltnisse[0]
			if gitter == 'Zinkblende':
				reflexe_Zinkblende = verhaeltnis_m
			elif gitter == 'Steinsalz':
				reflexe_Steinsalz = verhaeltnis_m
			elif gitter == 'Caesiumchlorid':
				reflexe_Caesiumchlorid = verhaeltnis_m
			elif gitter == 'Fluorit':
				reflexe_Fluorit = verhaeltnis_m
	
			verhaeltnis_d = np.sort(verhaeltnisse[1])
	
			print('sqrt(m_i/m_1): ', verhaeltnis_m)
			print('d_1/d_i: ', verhaeltnis_d)
	
			print('Verhältnisse für die m Werte: ', verhaeltnis_m)
			print('Verhältnisse für die d Werte: ', verhaeltnis_d)
			print('Abweichung Xi^2 für die ' + gitter + ' Struktur für ' + argg + ': ', abweichung(verhaeltnis_m, verhaeltnis_d))
	
			if abweichung(verhaeltnis_m, verhaeltnis_d) < Xi2_best[1]:
				Xi2_best = [gitter, abweichung(verhaeltnis_m, verhaeltnis_d), argg]

			############## CsCl ################## -> könnte komisch sein

			argg = 'CsCl'
			print('Berechnung '+gitter+ ' für ' + argg + '(könnte komisch sein):')
			infos = Strukturamplitude_Salz(54, 18, gitter = gitter)
			reflexe = np.array(infos[0])
			# print(gitter +': ',reflexe[np.argsort(infos[1])])
			m = infos[1]
			m = np.sort(m)
	
			verhaeltnisse = findStructure(m, netzebenenabstand)
			verhaeltnisse_temp = verhaeltnisse[0]
			# verhaeltnis_m = np.sort(verhaeltnisse[0])
			verhaeltnis_m = verhaeltnisse[0]
			if gitter == 'Zinkblende':
				reflexe_Zinkblende = verhaeltnis_m
			elif gitter == 'Steinsalz':
				reflexe_Steinsalz = verhaeltnis_m
			elif gitter == 'Caesiumchlorid':
				reflexe_Caesiumchlorid = verhaeltnis_m
			elif gitter == 'Fluorit':
				reflexe_Fluorit = verhaeltnis_m
	
			verhaeltnis_d = np.sort(verhaeltnisse[1])
	
			print('sqrt(m_i/m_1): ', verhaeltnis_m)
			print('d_1/d_i: ', verhaeltnis_d)
	
			print('Verhältnisse für die m Werte: ', verhaeltnis_m)
			print('Verhältnisse für die d Werte: ', verhaeltnis_d)
			print('Abweichung Xi^2 für die ' + gitter + ' Struktur für ' + argg + ': ', abweichung(verhaeltnis_m, verhaeltnis_d))
	
			if abweichung(verhaeltnis_m, verhaeltnis_d) < Xi2_best[1]:
				Xi2_best = [gitter, abweichung(verhaeltnis_m, verhaeltnis_d), argg]
		elif gitter == 'Fluorit':
			############## CaF ##################
			argg = 'CaF'
			print('Berechnung '+gitter+ ' für ' + argg + ':')
			infos = Strukturamplitude_Salz(18, 10, gitter = gitter)
			reflexe = np.array(infos[0])
			# print(gitter +': ',reflexe[np.argsort(infos[1])])
			m = infos[1]
			m = np.sort(m)
	
			verhaeltnisse = findStructure(m, netzebenenabstand)
			verhaeltnisse_temp = verhaeltnisse[0]
			# verhaeltnis_m = np.sort(verhaeltnisse[0])
			verhaeltnis_m = verhaeltnisse[0]
			if gitter == 'Zinkblende':
				reflexe_Zinkblende = verhaeltnis_m
			elif gitter == 'Steinsalz':
				reflexe_Steinsalz = verhaeltnis_m
			elif gitter == 'Caesiumchlorid':
				reflexe_Caesiumchlorid = verhaeltnis_m
			elif gitter == 'Fluorit':
				reflexe_Fluorit = verhaeltnis_m
	
			verhaeltnis_d = np.sort(verhaeltnisse[1])
	
			print('sqrt(m_i/m_1): ', verhaeltnis_m)
			print('d_1/d_i: ', verhaeltnis_d)
	
			print('Verhältnisse für die m Werte: ', verhaeltnis_m)
			print('Verhältnisse für die d Werte: ', verhaeltnis_d)
			print('Abweichung Xi^2 für die ' + gitter + ' Struktur für ' + argg + ': ', abweichung(verhaeltnis_m, verhaeltnis_d))
	
			if abweichung(verhaeltnis_m, verhaeltnis_d) < Xi2_best[1]:
				Xi2_best = [gitter, abweichung(verhaeltnis_m, verhaeltnis_d), argg]

			############## CaCl ################## -> könnte komisch sein

			argg = 'CaCl'
			print('Berechnung '+gitter+ ' für ' + argg + '(könnte komisch sein):')
			infos = Strukturamplitude_Salz(18, 18, gitter = gitter)
			reflexe = np.array(infos[0])
			# print(gitter +': ',reflexe[np.argsort(infos[1])])
			m = infos[1]
			m = np.sort(m)
	
			verhaeltnisse = findStructure(m, netzebenenabstand)
			verhaeltnisse_temp = verhaeltnisse[0]
			# verhaeltnis_m = np.sort(verhaeltnisse[0])
			verhaeltnis_m = verhaeltnisse[0]
			if gitter == 'Zinkblende':
				reflexe_Zinkblende = verhaeltnis_m
			elif gitter == 'Steinsalz':
				reflexe_Steinsalz = verhaeltnis_m
			elif gitter == 'Caesiumchlorid':
				reflexe_Caesiumchlorid = verhaeltnis_m
			elif gitter == 'Fluorit':
				reflexe_Fluorit = verhaeltnis_m
	
			verhaeltnis_d = np.sort(verhaeltnisse[1])
	
			print('sqrt(m_i/m_1): ', verhaeltnis_m)
			print('d_1/d_i: ', verhaeltnis_d)
	
			print('Verhältnisse für die m Werte: ', verhaeltnis_m)
			print('Verhältnisse für die d Werte: ', verhaeltnis_d)
			print('Abweichung Xi^2 für die ' + gitter + ' Struktur für ' + argg + ': ', abweichung(verhaeltnis_m, verhaeltnis_d))
	
			if abweichung(verhaeltnis_m, verhaeltnis_d) < Xi2_best[1]:
				Xi2_best = [gitter, abweichung(verhaeltnis_m, verhaeltnis_d), argg]


	print('Struktur mit der kleinsten Abweichung: ', Xi2_best[0], ' für ', Xi2_best[2])
	print('Abweichung Xi^2: ', Xi2_best[1])

	m = Strukturamplitude_Salz(18, 18,gitter = Xi2_best[0])[1]
	a = np.array(gitterkonstanteBragg(m, netzebenenabstand))

	print('Gitterkonstanten für ' + Xi2_best[0] + ' Struktur: ', ' für ', Xi2_best[2], a)

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
	plt.savefig('Plots/Salz_Fit_Egor.pdf')
	plt.close()

	i = np.linspace(1,len(daten),len(daten))

	verhaeltnisse_daten = []
	for j in range(0,len(daten)):	
		verhaeltnisse_daten.append(unp.sin(daten[j]) / unp.sin(daten[0]))

	plt.plot(i, reflexe_Zinkblende, 'x', label = 'Zinkblende')
	plt.plot(i, reflexe_Steinsalz, 'x', label = 'Steinsalz')
	plt.plot(i, reflexe_Caesiumchlorid, 'x', label = 'Caesiumchlorid')
	plt.plot(i, reflexe_Fluorit, 'o', label = 'Fluorit')
	plt.plot(i, noms(verhaeltnisse_daten), 'x', label = 'data')
	plt.xlabel('i')
	# plt.xlim(0,7)
	# plt.ylim(0,4.5)
	plt.ylabel('verhältnisse')
	plt.legend(loc = 'best')
	plt.tight_layout()
	plt.grid()
	plt.savefig('Plots/verhaeltnisse_Salz.pdf')
	plt.close()

main()


