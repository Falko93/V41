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


lambda_1 = 1.5417e-10 # in m
Radius_k = 0.0573 #in m
radius_rohr = 0.4e-3
Abstand_FP = 0.130 # in m -> Abstand Fokus Probe
umfang_kamera = 2*np.pi*Radius_k
 

# radius_m = np.array([4.5, 6.25, 8.2, 9.8, 11.5]) # Radius der Kreibögen in cm für Metall
# radius_m = np.array([0.045, 0.0625, 0.082, 0.098, 0.115]) # Radius der Kreibögen in m für Metall # 8.2 und neuer 4.5 rückwertig
radius_m = np.array([0.065, 0.098, 0.115]) # Radius der Kreibögen in m für Metall # 8.2 und neuer 4.5 rückwertig
radius_m_rueck_temp = np.array([0.045, 0.065, 0.082])
radius_m_rueck = umfang_kamera - radius_m_rueck_temp
# durchmesser_s = np.array([0.027, 0.032, 0.056, 0.06, 0.062, 0.079, 0.082, 0.095, 0.103, 0.112, 0.114, 0.125, 0.129, 0.14, 0.15, 0.154, 0.164]) # Radius der Kreibögen in m für Salz -> 9.5, 12.5, 12.9 und 15.4 cm scheint mit Ring zu sein.. ; 2,7 und 3,2 und 6,0 und 6,2 und 8,2 rückwertig
durchmesser_s = np.array([0.056, 0.079, 0.095, 0.112, 0.125, 0.15]) # Radius der Kreibögen in m für Salz -> 9.5, 12.5, 12.9 und 15.4 cm scheint mit Ring zu sein..
durchmesser_s_rueck = np.array([0.0295, 0.061, 0.082, 0.103, 0.114, 0.129, 0.14, 0.154, 0.164])
radius_s = durchmesser_s / 2
radius_s_rueck_temp = durchmesser_s_rueck / 2
radius_s_rueck = umfang_kamera - radius_s_rueck_temp
print(radius_s)
print(radius_s_rueck_temp)
Fehler_radius = 0.001 # in m

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


	# for h in range(0,10):
	# 	for k in range(0,h+1):
	# 		for l in range(0,k+1):
	# 			temp = h**2 + k**2 + l**2
	# 			if temp not in m and temp != 0:
	# 				S = 0
	# 				for i in range(0,len(r)):
	# 					S = S + f * np.exp(-2j * np.pi * (r[i][0] * h + r[i][1] * k + r[i][2] * l))
	# 				m.append(h**2 + k**2 + l**2)
	# 				if 0.001 < abs(S): # es wird nie ganz Null weil Computaional Physics
	# 					reflexe.append([h, k, l])
	# 					m_reflexe.append(h**2 + k**2 + l**2)
	# 			else:
	# 				pass
	# infos = [reflexe, m_reflexe]

	# return infos

def formfaktor(thet, lam, ion): # ion = 'F', 'Cl', 'K', 'Ca', 'Cu', 'J' oder 'Cs' 
	# was passiert bei Lücken in der Tabelle?
	sin_theta = np.sin(thet / lam)
	txtDatei = 'atomfaktoren.txt'
	formfaktor = np.genfromtxt(txtDatei , unpack = True)
	index_formfaktor = 0
	if ion == 'F':
		index_formfaktor = 0
	elif ion == 'Cl':
		index_formfaktor = 1
	elif ion == 'K':
		index_formfaktor = 2
	elif ion == 'Ca':
		index_formfaktor = 3
	elif ion == 'Cu':
		index_formfaktor = 4
	elif ion == 'J':
		index_formfaktor = 5
	elif ion == 'Cs':
		index_formfaktor = 6

	thetaSin1 = np.linspace(0.0, 0.4, 9)
	thetaSin2 = np.linspace(0.5, 0.7, 3)

	if sin_theta <= 0.4:
		for i in range(0,9):
			if sin_theta < thetaSin1[i] + (0.05/2) and sin_theta >= thetaSin1[i] - (0.05/2):
				# output = [ion, formfaktor[index_formfaktor][i]]
				output = ormfaktor[index_formfaktor][i]
				return output
			else:
				pass
	elif sin_theta > 0.4 and sin_theta <= 0.7:
		for i in range(0,3):
			if sin_theta < thetaSin2[i] + (0.1/2) and sin_theta >= thetaSin2[i] - (0.1/2):
				# output = [ion, formfaktor[index_formfaktor][i+9]]
				output = formfaktor[index_formfaktor][i+9]
				return output
			else:
				pass
	else:
		print('ERROR: sin(theta/lam) nicht bekannt!')
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

	####################################################################################################
	# Metall 
	####################################################################################################
	print('\n#################### Analyse für Metall ####################\n')
	# print('\n')

	#Bogenlängenformel b=alpha*r*pi /180
	# theta_4 = radius_m * 180 / (np.pi * Radius_k)
	theta_4_temp = radius_m * 180 / (np.pi * Radius_k) # in Grad!
	theta_4_rueck = (radius_m_rueck * 180 / (np.pi * Radius_k)) #+ 180 # in Grad!
	theta_4 = np.concatenate((theta_4_temp, theta_4_rueck), axis = 0)
	print('Hier!!!!!!!!!!!:' ,theta_4)
	theta = theta_4 * np.pi / (4 * 180) # in radianten
	theta = np.sort(theta)
	theta_sin = np.sin(theta)

	# hier für die Fehlerfortpfalnzung aufgrund der ungenauigkeit des Lineals
	radius_metall_temp = unp.uarray(radius_m, Fehler_radius)
	radius_metall_rueck = unp.uarray(radius_m_rueck, Fehler_radius)
	radius_metall = np.concatenate((radius_metall_temp, radius_metall_rueck), axis = 0)
	Radius_kamera = ufloat(Radius_k, 0.0)
	theta_4_unp = radius_metall * 180 / (np.pi * Radius_kamera)
	theta_unp = theta_4_unp * np.pi / (4 * 180) # in radianten
	theta_unp = np.sort(theta_unp)

	print('Theta mit Fehler: ', theta_unp)
	
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

	# linearer fit für a gegen cos^2
	params, cov = curve_fit(lin,np.cos(theta)**2,a)
	# params, cov = curve_fit(lin,np.cos(theta),a)
	err = np.sqrt(np.diag(cov))

	a_extrp = ufloat(params[1], err[1])

	# systematischer Fehler Absorption der Röntgenstrahlung

	DeltaA = (radius_rohr / (2 * Radius_k)) * (1 - Radius_k / Abstand_FP) * (np.cos(theta)**2 / theta) * a
	# cos2Theta = []
	# for i in range(0,len(theta)):
	# 	cos2Theta.append(ufloat(theta[i], Fehler_Theta[i]))
	cos2Theta = unp.cos(theta_unp)**2

	a_mitFehler = unp.uarray(a, DeltaA)
	print('Gitterkonstanten mit Fehler durch die Absorption: ', a_mitFehler)

	print('Extrapolierte Gitterkonstante: ', a_extrp)

	cos2Theta_fit = np.linspace(0, 1)

	####### ab hier der a gegen cos^2 Plot ########
	# plt.plot(np.cos(theta)**2, a, 'x', label = 'Daten')
	plt.errorbar(np.cos(theta)**2, a, xerr=stds(cos2Theta), yerr=DeltaA, fmt='x', label = 'Daten')
	# plt.plot(np.cos(theta), a, 'x', label = 'Daten')
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

	####################################################################################################
	# Salz 
	####################################################################################################
	print('\n#################### Analyse für Salz ####################\n')

	gitter_moegl_Salz = ['Zinkblende', 'Steinsalz', 'Caesiumchlorid', 'Fluorit']
	Xi2_best_Salz = ['Zinkblende', 10, 'CuCl']

	ion_zink = [['Cu', 'Cl'], ['Cs', 'J'], ['Cu', 'J']]

	#Bogenlängenformel b=alpha*r*pi /180
	# theta_4 = radius_s * 180 / (np.pi * Radius_k)
	theta_4_temp = radius_s * 180 / (np.pi * Radius_k) # in Grad!
	theta_4_rueck = (radius_s_rueck * 180 / (np.pi * Radius_k))# + 180 # in Grad!
	theta_4 = np.concatenate((theta_4_temp, theta_4_rueck), axis = 0)
	print('Hier!!!!!!!!!!!:' ,theta_4)
	theta = theta_4 * np.pi / (4 * 180) # in radianten
	theta = np.sort(theta)
	theta_sin = np.sin(theta)

	# hier für die Fehlerfortpfalnzung aufgrund der ungenauigkeit des Lineals
	# radius_salz = unp.uarray(radius_s, Fehler_radius)
	radius_salz_temp = unp.uarray(radius_s, Fehler_radius)
	radius_salz_rueck = unp.uarray(radius_s_rueck, Fehler_radius)
	radius_salz = np.concatenate((radius_salz_temp, radius_salz_rueck), axis = 0)
	print('Hier!!!???!?!?!?!: ', radius_salz)
	Radius_kamera = ufloat(Radius_k, 0.0)
	theta_4_unp = radius_salz * 180 / (np.pi * Radius_kamera)
	theta_unp = theta_4_unp * np.pi / (4 * 180) # in radianten
	theta_unp = np.sort(theta_unp)

	print('Theta mit Fehler: ', theta_unp)
	
	netzebenenabstand_salz = bragg(lambda_1, theta)

	for gitter in gitter_moegl_Salz:
		# reflexe = []
		# m = []
		# if gitter == 'Zinkblende':
		# 	for ion in ion_zink:
		# 		temp_Salz = Strukturamplitude_Salz(formfaktor(, lambda_1,ion[0]), formfaktor(ion[1]), 'Zinkblende') # wie mit theta????
		# 		reflexe = temp_Salz[0]
		# 		m = temp_Salz[1]
		# 		verhaeltnisse = findStructure(m, netzebenenabstand_salz)
		# 		verhaeltnis_m = verhaeltnisse[0]
		# 		verhaeltnis_d = verhaeltnisse[1]
			
		# 		print('Verhältnisse für die m Werte: ', verhaeltnis_m)
		# 		print('Verhältnisse für die d Werte: ', verhaeltnis_d)
		# 		print('Abweichung Xi^2 für die ' + gitter + ' Sturktur mit' + ion[0] + ion[1] + ' : ', abweichung(verhaeltnis_m, verhaeltnis_d))
		
		# 		if abweichung(verhaeltnis_m, verhaeltnis_d) < Xi2_best_Salz[1]:
		# 			Xi2_best_Salz = [gitter, abweichung(verhaeltnis_m, verhaeltnis_d), ion[0] + ion[1]]
		# elif gitter == 'Steinsalz':
		# 	temp_Salz = Strukturamplitude_Salz(ion[0], ion[1], 'Zinkblende')
		# 	reflexe = temp_Salz[0]
		# 	m = temp_Salz[1]
		# 	verhaeltnisse = findStructure(m, netzebenenabstand_salz)
		# 	verhaeltnis_m = verhaeltnisse[0]
		# 	verhaeltnis_d = verhaeltnisse[1]
		
		# 	print('Verhältnisse für die m Werte: ', verhaeltnis_m)
		# 	print('Verhältnisse für die d Werte: ', verhaeltnis_d)
		# 	print('Abweichung Xi^2 für die ' + gitter + ' Sturktur mit' + ion[0] + ion[1] + ' : ', abweichung(verhaeltnis_m, verhaeltnis_d))
	
		# 	if abweichung(verhaeltnis_m, verhaeltnis_d) < Xi2_best_Salz[1]:
		# 		Xi2_best_Salz = [gitter, abweichung(verhaeltnis_m, verhaeltnis_d), ion[0] + ion[1]]

		infos = Strukturamplitude_Salz(f1 = 1, f2 = 2, gitter = gitter)
		reflexe = infos[0]
		m = infos[1]
	
		verhaeltnisse = findStructure(m, netzebenenabstand_salz)
		verhaeltnis_m = verhaeltnisse[0]
		verhaeltnis_d = verhaeltnisse[1]
	
		print('Verhältnisse für die m Werte: ', verhaeltnis_m)
		print('Verhältnisse für die d Werte: ', verhaeltnis_d)
		print('Abweichung Xi^2 für die ' + gitter + ' Sturktur: ', abweichung(verhaeltnis_m, verhaeltnis_d))

		if abweichung(verhaeltnis_m, verhaeltnis_d) < Xi2_best_Salz[1]:
			Xi2_best_Salz = [gitter, abweichung(verhaeltnis_m, verhaeltnis_d)]

	print('Struktur mit der kleinsten Abweichung: ', Xi2_best_Salz[0])
	print('Abweichung Xi^2: ', Xi2_best_Salz[1])
		
	m = Strukturamplitude_Salz(gitter = Xi2_best_Salz[0])[1]
	a = np.array(gitterkonstanteBragg(m, netzebenenabstand_salz))

	print('Gitterkonstanten für ' + Xi2_best_Salz[0] + ' Struktur: ', a)

	#linearer fit für a gegen cos^2
	params, cov = curve_fit(lin,np.cos(theta)**2,a)
	# params, cov = curve_fit(lin,np.cos(theta),a)
	err = np.sqrt(np.diag(cov))

	a_extrp = ufloat(params[1], err[1])

	# systematischer Fehler Absorption der Röntgenstrahlung

	DeltaA = (radius_rohr / (2 * Radius_k)) * (1 - Radius_k / Abstand_FP) * (np.cos(theta)**2 / theta) * a
	# cos2Theta = []
	# for i in range(0,len(theta)):
	# 	cos2Theta.append(ufloat(theta[i], Fehler_Theta[i]))
	cos2Theta = unp.cos(theta_unp)**2

	a_mitFehler = unp.uarray(a, DeltaA)
	print('Gitterkonstanten mit Fehler durch die Absorption: ', a_mitFehler)

	print('Extrapolierte Gitterkonstante: ', a_extrp)

	cos2Theta_fit = np.linspace(0, 1)

	####### ab hier der a gegen cos^2 Plot ########
	# plt.plot(np.cos(theta)**2, a, 'x', label = 'Daten')
	plt.errorbar(np.cos(theta)**2, a, xerr=stds(cos2Theta), yerr=DeltaA, fmt='x', label = 'Daten')
	# plt.plot(np.cos(theta), a, 'x', label = 'Daten')
	plt.plot(cos2Theta_fit, lin(cos2Theta_fit, *params), label = 'Fit')
	plt.legend(loc = 'best')
	plt.xlabel('cos$^2(\Theta)$')
	plt.ylabel('$a$ in Angtröm')
	# plt.xlim(0.5,1)
	# plt.ylim(4.8e-10,5.8e-10)
	plt.tight_layout()
	plt.grid()
	plt.savefig('Plots/Salz_Fit.pdf')
	plt.close()


main()



















