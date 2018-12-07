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
    radius_vorne = 2*radius_vorne[:-3]
    radius_rueck = umfang_kamera - 2*radius_rueck
    radius_rueck = radius_rueck[::-1]
    radius = np.concatenate((radius_vorne, radius_rueck), axis = 0)

    return radius



def Millerindizes(): # gitter = 'SC', 'BCC', 'FCC', 'Diamant' ; bei Fluorit f1=nicht Fluorit und f2=Fluorit
    	m = []
    	m_reflexe = []
    	temp = 0
    	reflexe = []
    	for h in range(0,10):
    		for k in range(0,h+1):
    			for l in range(0,k+1):
    				temp = h**2 + k**2 + l**2
    				if temp not in m and temp != 0:
    					m.append(h**2 + k**2 + l**2)
    					reflexe.append([h, k, l])
    					m_reflexe.append(h**2 + k**2 + l**2)
    				else:
    					pass
    	infos = [reflexe, m_reflexe]

    	return infos




####################################################################################################
# und los
####################################################################################################

def main():


    print('\n#################### Analyse für Salz ####################\n')

    daten1, daten2 = np.genfromtxt('SalzRadius.txt', unpack=True)
    daten = np.append(daten1, daten2)
    daten = daten[daten != 0]
    daten = (daten * np.pi) / (360)

    print('Daten: ', daten)

    # theta = theta_radiant(radius())
    theta = daten
    print('Theta mit Fehler: ', theta)
    netzebenenabstand = bragg(lambda_1, noms(theta))
    print('Netzebenenabstand', netzebenenabstand)
    print('Radius des Kegels:', radius())


    miller = Millerindizes()
    miller = miller[0]
    print('millerindizes', miller)

    print('\n#################### Verhältnisse d1/di ####################\n')

    verhaeltnis_d = netzebenenabstand[0]/netzebenenabstand
    verhaeltnis_d=np.array(verhaeltnis_d)
    print('d1/di:', verhaeltnis_d)



    print('\n#################### Cäsiumchlorid ####################\n')
    print('Starke Reflexe normiert auf Wurzel2 aus h,k,l = 1,1,0, denn erster Starker reflex')

    miller = np.array(miller[1:]) # denn erster starker reflex bei 1,1,0
#    print('miller', miller)


    info =[]
    miller_stark=[]
    miller_schwach=[]

    m_stark = []
    m_schwach = []

    info_stark =[]
    info_schwach = []

    for hkl in miller:
        if np.sum(hkl)%2 == 0:
            miller_stark.append(hkl) #starke reflexe'
            m_stark.append(np.sqrt(hkl[0]**2+hkl[1]**2+hkl[2]**2))
        else:
            miller_schwach.append(hkl)# schwache reflexe
            m_schwach.append(np.sqrt(hkl[0]**2+hkl[1]**2+hkl[2]**2))


    miller_stark = np.array(miller_stark)
    miller_schwach = np.array(miller_schwach)

    m_stark = np.array(m_stark)
    m_schwach = np.array(m_schwach)

    m_stark_ver = m_stark/m_stark[0] #für normierung
    m_schwach_ver = m_schwach/m_stark[0] #für normierung


    for i in range(0, len(miller_stark)):
        info_stark.append([miller_stark[i], m_stark_ver[i], 'stark'])



    for j in range(0, len(miller_schwach)):
        info_schwach.append([miller_schwach[j], m_schwach_ver[j], 'schwach'])


    info = info_stark+info_schwach
    info.sort(key = lambda x: x[1]) # sortieren der verhältnisse nach Größe
    info = np.array(info)


    for p in range(0, len(info)):
        if info[p][1]< verhaeltnis_d[len(verhaeltnis_d)-1]+0.1:
            p=p+1
        else:
            info = info[:p][:]
            break

    print('Infos Ceasium ', info)

    print('\n#################### Steinsalz ####################\n')
    print('Starke Reflexe normiert auf Wurzel2 aus h,k,l = 2,0,0, denn erster Starker reflex')

    miller = np.array(miller[2:]) # denn erster starker reflex bei 1,1,0
#    print('miller', miller)


    info =[]
    miller_stark=[]
    miller_schwach=[]

    m_stark = []
    m_schwach = []

    info_stark =[]
    info_schwach = []

    verbot_steinsalz = np.arange(1,29,2)

    for hkl in miller:
        if (hkl % 2 != 0).all() or np.sum(hkl) in verbot_steinsalz:
            miller_schwach.append(hkl)# schwache reflexe
            m_schwach.append(np.sqrt(hkl[0]**2+hkl[1]**2+hkl[2]**2))
        else:
            miller_stark.append(hkl) #starke reflexe'
            m_stark.append(np.sqrt(hkl[0]**2+hkl[1]**2+hkl[2]**2))

    miller_stark = np.array(miller_stark)
    miller_schwach = np.array(miller_schwach)

    m_stark = np.array(m_stark)
    print('m_stark=', m_stark[0])
    m_schwach = np.array(m_schwach)

    m_stark_ver = m_stark/m_stark[0] #für normierung
    m_schwach_ver = m_schwach/m_stark[0] #für normierung


    for i in range(0, len(miller_stark)):
        info_stark.append([miller_stark[i], m_stark_ver[i], 'stark'])



    for j in range(0, len(miller_schwach)):
        info_schwach.append([miller_schwach[j], m_schwach_ver[j], 'schwach'])


    info = info_stark+info_schwach
    info.sort(key = lambda x: x[1]) # sortieren der verhältnisse nach Größe
    info = np.array(info)


    for p in range(0, len(info)):
        if info[p][1]< verhaeltnis_d[len(verhaeltnis_d)-1]+0.1:
            p=p+1
        else:
            info = info[:p][:]
            break

#    print('Infos Steinsalz ', info)


    print('\n#################### Flourit ####################\n')
    print('Starke Reflexe normiert auf Wurzel1=1 aus h,k,l = 1,0,0, denn erster Starker reflex')

    miller = Millerindizes()
    miller = np.array(miller[0]) # denn erster starker reflex bei 1,0,0
#    print('miller', miller)


    info =[]
    miller_stark=[]
    miller_schwach=[]

    m_stark = []
    m_schwach = []

    info_stark =[]
    info_schwach = []

    verbot_flourit = np.arange(4,29,6)

    for hkl in miller:
        if (hkl % 2 == 0).all() or np.sum(hkl) in verbot_flourit:
            miller_schwach.append(hkl)# schwache reflexe
            m_schwach.append(np.sqrt(hkl[0]**2+hkl[1]**2+hkl[2]**2))
        else:
            miller_stark.append(hkl) #starke reflexe'
            m_stark.append(np.sqrt(hkl[0]**2+hkl[1]**2+hkl[2]**2))

    miller_stark = np.array(miller_stark)
    miller_schwach = np.array(miller_schwach)

    m_stark = np.array(m_stark)
    print('m_stark=', m_stark[0])
    m_schwach = np.array(m_schwach)

    m_stark_ver = m_stark/m_stark[0] #für normierung
    m_schwach_ver = m_schwach/m_stark[0] #für normierung


    for i in range(0, len(miller_stark)):
        info_stark.append([miller_stark[i], m_stark_ver[i], 'stark'])



    for j in range(0, len(miller_schwach)):
        info_schwach.append([miller_schwach[j], m_schwach_ver[j], 'schwach'])


    info = info_stark+info_schwach
    info.sort(key = lambda x: x[1]) # sortieren der verhältnisse nach Größe
    info = np.array(info)


    for p in range(0, len(info)):
        if info[p][1]< verhaeltnis_d[len(verhaeltnis_d)-1]+0.1:
            p=p+1
        else:
            info = info[:p][:]
            break


    print('Infos Flourit ', info)


############################ Fit und Plot #####################################
    print('\n#################### Flourit ####################\n')
############################ da flourit erster reflex [100] folgt: muss nicht mehr multiplizieren mit der Norm da Norm =1


    m = []
    for q in range(0, len(info)):
        m.append(info[q][1])

############################ m*d rechnen für Gitterkonstante a
    d = lambda_1/(2*unp.sin(theta))
    a = np.array(m)*d
    print('Gitterkonstante', a)


    # linearer fit für a gegen cos^2
    params, cov = curve_fit(lin,np.cos(noms(theta))**2,noms(a))
    # params, cov = curve_fit(lin,np.cos(theta),a)
    err = np.sqrt(np.diag(cov))

    a_extrp = ufloat(params[1], err[1])

    # systematischer Fehler Absorption der Röntgenstrahlung

    DeltaA = (radius_rohr / (2 * Radius_kamera)) * (1 - Radius_kamera / Abstand_FP) * (unp.cos(theta)**2 / theta) * a

    # cos2Theta = []
    # for i in range(0,len(theta)):
    # 	cos2Theta.append(ufloat(theta[i], Fehler_Theta[i]))
    cos2Theta = unp.cos(theta)**2

    a_mitFehler = unp.uarray(noms(a), noms(DeltaA))
    print('Gitterkonstanten mit Fehler durch die Absorption: ', a_mitFehler)

    print('Extrapolierte Gitterkonstante: ', a_extrp)






    cos2Theta_fit = np.linspace(0, 1)

    # plt.plot(np.cos(theta)**2, a, 'x', label = 'Daten')
    plt.errorbar(np.cos(noms(theta))**2, noms(a), xerr=stds(cos2Theta), yerr=noms(DeltaA), fmt='x', label = 'Daten')
    # plt.plot(np.cos(theta), a, 'x', label = 'Daten')
    plt.plot(cos2Theta_fit, lin(cos2Theta_fit, *params), label = 'Fit')
    plt.legend(loc = 'best')
    plt.xlabel('cos$^2(\Theta)$')
    plt.ylabel('$a$ in Angtröm')
    # plt.xlim(0.6,1)
    # plt.ylim(4.8e-10,5.8e-10)
    plt.tight_layout()
    plt.grid()
    plt.savefig('Plots/Salz2_Fit.pdf')
    plt.close()



main()
