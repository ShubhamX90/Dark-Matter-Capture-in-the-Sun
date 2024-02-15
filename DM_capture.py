import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.style

from DarkCapPy.Configure.Constants  import *
from DarkCapPy.Configure.AtomicData import *
from DarkCapPy.Configure.PlanetData  import *
from DarkCapPy.Configure.Conversions import amu2GeV

def helm_formFactor(element, E_R):
    E_N = 0.114/((atomicNumbers[element])**(5./3))
    FN2 = np.exp(-E_R/E_N)
    
    return FN2

def diff_crossSec(element, E_R, m_X):
    m_N = amu2GeV(atomicNumbers[element])
    FN2 = helm_formFactor(element, E_R)
    mu_N = (m_X*m_N) / (m_X + m_N)

    sigma_O = 10**(-44)

    c_Sec = (m_N/(2*(mu_N**2))) * sigma_O * FN2

    return c_Sec


def E_min(u, m_X):
    e_min = 0.5 * m_X * (u**2)
    
    return e_min


def E_max(element, m_X, u, index_esVel):
    m_N = amu2GeV(atomicNumbers[element])
    mu_N = (m_N*m_X) / (m_N + m_X)

    V_escape2 = (escVel2_List[index_esVel])

    e_max = (2 * (mu_N**2) / m_N) * ((u**2) + (V_escape2))
    
    return e_max


def intersection_Vel(element, m_X, index_esVel):
    m_N = amu2GeV(atomicNumbers[element])
    mu_N = (m_N*m_X) / (m_N+m_X)

    sqrt_V_escape2 = np.sqrt(escVel2_List[index_esVel])

    A = 0.5 * m_X 
    B = 2 * (mu_N**2) / m_N
    uInt = np.sqrt( ( B ) / (A-B) ) * sqrt_V_escape2

    return uInt


def Velocity_Energy_Integral(element, m_X, index_esVel):
    
	def integrand(E_R, u):
		fu = fCrossInterp(u)
		integrand = diff_crossSec(element, E_R, m_X) * u * fu

		return integrand

	# Calculate the intersection uInt of eMin and eMax given a specific rIndex
	uInt = intersection_Vel(element, m_X, index_esVel)

	uLow = 0
	uHigh = min(uInt, V_gal) # We take the minimal value between the intersection velocity and galactic escape velocity
	eLow = lambda u: E_min(u, m_X)
	eHigh = lambda u: E_max(element, m_X, index_esVel, u)
	integral = integrate.dblquad(integrand, uLow, uHigh, eLow, eHigh)[0]
	return integral

def sumOverR(element, m_X):

	tempSum = 0
    
	for i in range(0, len(radius_List)):
		r = radius_List[i]
		deltaR = deltaR_List[i]

		n_N = numDensity_Func(element)[i]

		summand = n_N * r**2 * Velocity_Energy_Integral(element, m_X, i) * deltaR
		tempSum += summand

	return tempSum


def singleElementCap(element, m_X):

	Z_N = nProtons[element]
	m_N = amu2GeV(atomicNumbers[element])
	n_X = 0.3/m_X # GeV/cm^3

	conversion = (5.06e13)**-3 * (1.52e24) # Conversion to seconds (cm^-3)(GeV^-2) -> (s^-1)
	prefactors = (4*np.pi)**2
	single_cap = n_X * conversion * prefactors * sumOverR(element, m_X)
	return single_cap


def cCap(m_X):

	totalCap = 0
	for element in element_List:
		elementCap = singleElementCap(element, m_X)
		totalCap += elementCap 
	return totalCap

def plot_capture_rate():
    m_X_values = np.linspace(1, 100, 100)  # Set the range of m_X values for the plot
    total_cap_values = []

    for m_X in m_X_values:
        total_cap = cCap(m_X)
        total_cap_values.append(total_cap)

    # Plotting
    plt.plot(m_X_values, total_cap_values, label='Total Capture Rate')
    plt.xlabel('m_X (GeV)')
    plt.ylabel('Total Capture Rate')
    plt.title('Dark Matter Capture Rate in the Sun')
    plt.legend()
    plt.grid(True)
    plt.show()
