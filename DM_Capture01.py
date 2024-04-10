import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interpolate

from DarkCapPy.Configure.Constants  import *
from DarkCapPy.Configure.AtomicData import *
from DarkCapPy.Configure.PlanetData  import *
from DarkCapPy.Configure.Conversions import amu2GeV

def eMax(element, m_X, rIndex, u): 
    m_A = amu2GeV(atomicNumbers[element])
    mu = (m_A*m_X) / (m_A + m_X)
    vCross2 = (escVel2_List[rIndex])
    E_max = 2 * mu**2 * (u**2 + vCross2) / m_A
    return E_max
    
def eMin(u, m_X):
    E_min = (0.5) * m_X * u**2
    return E_min



def formFactorIntegral(element, m_X, rIndex, u):
    E_max = eMax(element, m_X, rIndex, u)
    E_cap = eMin(u, m_X)
    if(atomicNumbers[element]== 1):
        return (E_max - E_cap)
    else :
        E_A = 0.114/((atomicNumbers[element])**(5./3))  
        FN_SI = E_A * (np.exp(((-1) * E_cap)/E_A) - np.exp(((-1) * E_max)/E_A))
        return FN_SI


def EminEmaxIntersection(element, m_X, rIndex):
	m_N = amu2GeV(atomicNumbers[element])
	mu = (m_N*m_X)/(m_N+m_X)

	sqrtvCross2 = np.sqrt(escVel2_List[rIndex])
    
	A = m_X/2. 
	B = 2. * mu**2 / m_N
	uInt = np.sqrt( ( B ) / (A-B) ) * sqrtvCross2

	return uInt


def velocityIntegral(element, m_X, rIndex, u):
     
     F_A = formFactorIntegral(element, m_X, rIndex, u)
     def integrand(E,u):
        fu = fCrossInterp(u)
        integrand = u * fu * F_A
        return integrand
     
     uInt = EminEmaxIntersection(element, m_X, rIndex)
     uLow = 0
     uHigh = min(uInt, V_gal)
     integral = integrate.quad(integrand, uLow, uHigh)
     return integral


def sumOverR(element, m_X, rIndex):
     tempSum = 0
     for i in range(0, len(radius_List)):
        r = radius_List[i]
        deltaR = deltaR_List[i]
        n_N = numDensity_Func(element)[i]
        summand = n_N * r**2 * velocityIntegral(element, m_X, rIndex, i) * deltaR
        tempSum += summand
     return tempSum

def sumOverElements(m_X, rIndex):
    totalCap =0
    for element in element_List:
        elementCap = (atomicNumbers[element] * 10**(-44) )**2 * sumOverR(element, m_X, rIndex)
        totalCap += elementCap 

    return totalCap

def totalCaptureRate()

     

