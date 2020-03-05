import numpy as np
import pandas as pd
from scipy.integrate import quad

def A_func(D):
	return 0.25 * np.pi * D**2

def Cp_func(T):
	return 2.224*10**(-10)*T**4-8.91*10**(-7)\
		       *T**3+1.218*10**(-3)*T**2-4.632*10**\
		       (-1)*T+1056
		
def mu_func(T):
	return -3.245*10**-18*T**4+1.65*10**-14*T**3\
		   -3.57*10**-11*T**2+6.32*10**-8*T+2.67*10**-6
		       
def R_func(M):
	return 8.314*10**3/M
		
def perfectGas_func(P,R,T):
	return P*10**6/(R*T)
		
def v_func(mu,Re,rho,D):
	return mu * Re / rho / D

def m_func(rho,v,A):
	return rho * v * A
	
def AL_func(D,L):
	return np.pi * D * L
	
def q_func(Q,AL):
	return Q/AL

if __name__ == '__main__':
	
	D = 0.00424 #m
	L = 0.76 #m
	P_in = 0.477 #Mpa
	T_in = 578.7 #k
	T_out = 1222 #k
	M = 28.7
	Re = 15000
	
	Ad = A_func(D)
	mu_in = mu_func(T_in)
	R = R_func(M)
	rho_in = perfectGas_func(P_in,R,T_in)
	v_in = v_func(mu_in,Re,rho_in,D)
	m = m_func(rho_in,v_in,Ad)
	f = lambda T : Cp_func(T)*m
	Q,er = quad(f,T_in,T_out)
	AL = AL_func(D,L)
	q = q_func(Q,AL)
	
	
	print(q)
