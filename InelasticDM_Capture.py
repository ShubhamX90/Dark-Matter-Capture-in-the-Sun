
################################################################
# Import Python Libraries
################################################################
import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interpolate

from DarkCapPy.Configure.Constants  import *
from DarkCapPy.Configure.AtomicData import *
from DarkCapPy.Configure.PlanetData  import *
from DarkCapPy.Configure.Conversions import amu2GeV


################################################################
# Capture Rate Functions
################################################################


########################
# Nuclear Form Factor
########################

def formFactor2(element, E_R):

    E_N = 0.114/((atomicNumbers[element])**(5./3))
    FN2 = np.exp(-E_R/E_N)
    return FN2




########################
# Photon Scattering Cross Sections
########################


def crossSection(element, E_R): # returns 1/GeV^3

	FN2 = formFactor2(element, E_R)
	
	function = ( ( FN2 ) * (atomicNumbers[element]**2) )    ##### excluded w^-2  and  m_N and mu_XP and sigma_O_SI
	return function






########################
# Kinematics
########################
def eCap(u, m_X, delta):

	function = ( (0.5) * m_X * u**2 ) - delta
	return function


def eMin(element, m_X, rIndex, u, delta):
	
	m_N = amu2GeV(atomicNumbers[element])
	mu = m_N*m_X / (m_N + m_X)
	vCross2 = (escVel2_List[rIndex])
	value1 = np.sqrt( ( 1  -  ( (2*delta) / (mu * (u**2 + vCross2)))))
	value2 = (mu * delta) / m_N
	
	function = ( (mu**2 * (u**2 + vCross2) / m_N)  *  (1 - value1 ) ) - value2
	return function


def eMax(element, m_X, rIndex, u, delta):
	
	m_N = amu2GeV(atomicNumbers[element])
	mu = m_N*m_X / (m_N + m_X)
	vCross2 = (escVel2_List[rIndex])
	value1 = np.sqrt( ( 1  -  ( (2*delta) / (mu * (u**2 + vCross2)))))
	value2 = (mu * delta) / m_N
	
	function = ( (mu**2 * (u**2 + vCross2) / m_N)  *  (1 + value1 ) ) - value2
	return function







##########################
# Intersection Velocities 
##########################


def EminEcapIntersection(element, m_X, rIndex, delta):
	
	m_N = amu2GeV(atomicNumbers[element])
	sqrtvCross2 = np.sqrt(escVel2_List[rIndex])
	
	A = np.sqrt(2* m_X * m_N) /  np.abs(m_X -m_N) * sqrtvCross2 
	B = delta * ( (m_X-m_N) / ( (m_X ** 2) * (sqrtvCross2 ** 2) ) )
	C = 2 * delta * ( (m_X-m_N) / ( m_X * m_N * (sqrtvCross2 ** 2) ) )

	
	uIntMin = A * np.sqrt( 1  -  B  - np.sqrt( 1  -  C)  )

	return uIntMin




def EmaxEcapIntersection(element, m_X, rIndex, delta):
	
	m_N = amu2GeV(atomicNumbers[element])
	sqrtvCross2 = np.sqrt(escVel2_List[rIndex])
	
	A = np.sqrt(2* m_X * m_N) /  np.abs(m_X -m_N) * sqrtvCross2 
	B = delta * ( (m_X-m_N) / ( (m_X ** 2) * (sqrtvCross2 ** 2) ) )
	C = 2 * delta * ( (m_X-m_N) / ( m_X * m_N * (sqrtvCross2 ** 2) ) )

	
	uIntMax = A * np.sqrt( 1  -  B  + np.sqrt( 1  -  C)  )

	return uIntMax






########################################
# Photon Velocity and Energy Integration
########################################

def intDuDEr(element, m_X, rIndex, delta):
    
	def integrand(E_R, u):                                    ####### E_R and u coming from where ????
		fu = fCrossInterp(u)
		integrand = crossSection(element, m_X, E_R) * u * fu

		return integrand

	# Calculate the intersection uInt of eMin and eMax given a specific rIndex
	uIntMin = EminEcapIntersection(element, m_X, rIndex, delta)
	uIntMax = EmaxEcapIntersection(element, m_X, rIndex, delta)

	uLow = uIntMin
	uHigh = min(uIntMax, V_gal) # We take the minimal value between the intersection velocity and galactic escape velocity
	
	eLow = lambda u: eCap(u, m_X, delta)
	eHigh = lambda u: eMax(element, m_X, rIndex, u, delta)
	
	integral = integrate.dblquad(integrand, uLow, uHigh, eLow, eHigh)[0]
	return integral







################
# Sum Over Radii
################

def sumOverR(element, m_X):

	tempSum = 0
    
	for i in range(0, len(radius_List)):
		r = radius_List[i]
		deltaR = deltaR_List[i]
		n_N = numDensity_Func(element)[i]
		summand = n_N * r**2 * intDuDEr(element, m_X, i) * deltaR
		tempSum += summand


	return tempSum






#############################
# Single Element Capture Rate
#############################

def singleElementCap(element, m_X):

	n_X = 0.3/m_X # GeV/cm^3
	m_P = amu2GeV(1)

	conversion = (5.06e13)**-3 * (1.52e24) # Conversion to seconds (cm^-3)(GeV^-2) -> (s^-1)  ###### WHY NEEDED ????
	prefactors = (4*np.pi)**2     ######### DOUBT 
	mu_XP = m_P*m_X / (m_P + m_X)
	
	function = n_X * conversion * (10^(-44)) * prefactors * sumOverR(element, m_X) /   (2*((mu_XP)**2))  
	return function







###################
# Full Capture Rate
###################

def cCap(m_X):

	totalCap = 0
	for element in element_List:
		elementCap = singleElementCap(element, m_X)
		# print ('Element:', element,',' 'Cap: ', elementCap)
		totalCap += elementCap 
	return totalCap





print ("Dark Photon Module Imported")