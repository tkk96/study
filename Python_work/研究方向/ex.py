import pandas as pd
import re
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
sys.setrecursionlimit(1000000)

class Paramters():
	
	def __init__(self):
		self.physicalPropertiesOfTestSection()
		self.physicalPropertiesOfPreHeater()
		
	def physicalPropertiesOfTestSection(self):
		self.Pin = 1 #MPa
		self.deltaP = 3000 #Pa
		self.Tout = 773 #K
		self.De = 0.01 #m 管壁内径
		self.Dout = 0.012 #inconcel管壁外径
		self.LDratio = 60 #长经比
		self.L = self.L_func(self.LDratio,self.De) #试验段长度
		self.Cp_out = self.cp_func(self.Tout) #出口定压比热
		self.mu_out = self.mu_func(self.Tout) #出口动力粘度
		self.M = 28.959 #空气相对分子量
		self.R = self.R_func(self.M) #空气理想气体常数
		self.Re_aim = 10000 #目标雷诺数
		self.deltaT = 300 #试验段进出口空气温差
		self.Tin = self.T_func(self.Tout,self.deltaT) #试验段进口空气温度
		self.Kair_out = self.K_func(self.Tout) #出口空气导热率
		self.kt = 0.155 #inconcel导热率
		self.Pr = self.Pr_func(self.mu_out,self.Cp_out,self.Kair_out)
	
	def physicalPropertiesOfPreHeater(self):
		self.P_pre = self.Pin * 1.1
		self.Tin_pre = 298
		self.Tout_pre = self.Tin * 1.1
	
	#def physicalPropertiesOfSurgeTank(self):
	
	def T_func(self,Tout,deltaT):
		return Tout - deltaT
	
	def Pr_func(self,mu,Cp,K):
		return mu * Cp / K
	
	def K_func(self,T):
		return -4.676*(10**-15)*(T**4)+2.257*(10**-11)*(T**3)\
		       -4.617*(10**-8)*(T**2)+9.25*(10**-5)*T+0.00217 #W/(m.k)
	
	def perfectGas_func(self,P,R,T):
		return P/(R*T)
	
	def Pout_func(self,Pin,deltaP):
		return Pin*10**6 - deltaP
	
	def L_func(self,LDratio,D):
		return LDratio * D
	
	def cp_func(self,T):
		return 2.224*10**(-10)*T**4-8.91*10**(-7)\
		       *T**3+1.218*10**(-3)*T**2-4.632*10**\
		       (-1)*T+1056
		             
	def mu_func(self,T):
		return -3.245*10**-18*T**4+1.65*10**-14*T**3\
		       -3.57*10**-11*T**2+6.32*10**-8*T+2.67*10**-6
	
	def R_func(self,M):
		return 8.314*10**3/M
	
	def f_func(self,Re):
		return 0.079*4/(Re**0.25)
	
	def Re_func(self,rho,v,D,mu):
		return rho*v*D/mu
	
	def v_func_fromP(self,deltaP,LDratio,f,rho):
		return np.sqrt(2*deltaP/LDratio/f/rho)
	
	def v_func_fromRe(self,Re,rho,D,mu):
		return Re*mu/rho/D
	
	def Pf_func(self,f,LDratio,rho,v):
		return f*LDratio*rho*(v**2)/2
	
	def A_func(self,D):
		return 1/4*np.pi*(D**2)
	
	def Tf_func(self,Tout,deltaT):
		return Tout - deltaT / 2
	
	def Tw_func(self,Tf,Q,h,A):
		return Tf + Q / h / A
	
	def Tt_func(self,Ttin,Q,R):
		'''Tout:inconel管外表面温度;R:inconel导热热阻'''
		return Ttin - Q * R
	
	def Rt_func(self,dout,din,kt,L):
		return np.log(dout/din)/(2*np.pi*kt*L)
	
	def DB_func(self,Re,Pr):
		return 0.023 * (Re ** 0.8) * (Pr ** 0.4)
	
	def h_func(self,Nu,K,D):
		return Nu * K / D
	
	def Re_iter(self,Re_init):
		f_init = self.f_func(Re_init)
		self.Pout = self.Pout_func(self.Pin,self.deltaP)
		self.rho_out = self.perfectGas_func(self.Pout,self.R,self.Tout)
		v_init = self.v_func_fromP(self.deltaP,self.LDratio,f_init,self.rho_out)
		Re_new = self.Re_func(self.rho_out,v_init,self.De,self.mu_out)
		if abs(Re_new - Re_init) < 10**-4:
			return Re_new
		else:
			return self.Re_iter(Re_new)
	
	def Re_deltaP_iter(self,Re=1000):
		'''get min deltaP'''
		Re0 = self.Re_iter(Re)
		if Re0 - self.Re_aim < 0:
			self.deltaP *= 2
			return self.Re_deltaP_iter()
		elif 0 <= Re0 - self.Re_aim < 10:
			return Re0
		elif Re0 - self.Re_aim > 10:
			self.deltaP *= 3/4
			return self.Re_deltaP_iter() 
		
		
class Main(Paramters):
	
	def __init__(self):
		Paramters.__init__(self)
	
	def __repr__(self,Re_init=10000):
		#试验段
		Re_max = self.Re_deltaP_iter(Re_init)
		v_out = self.v_func_fromRe(Re_max,self.rho_out,self.De,self.mu_out)
		A = self.A_func(self.De)
		Cp_out = self.Cp_out
		Q_out = Cp_out * self.rho_out * v_out * A * self.deltaT
		T_f = self.Tf_func(self.Tout,self.deltaT)
		Nu_out = self.DB_func(Re_max,self.Pr)
		h_out = self.h_func(Nu_out,self.Kair_out,self.De)
		T_w = self.Tw_func(T_f,Q_out,h_out,A)
		Rt = self.Rt_func(self.Dout,self.De,self.kt,self.L) #管壁导热热阻
		Tt_out = self.Tt_func(T_w,Q_out,Rt) #管壁外表面温度
		#预热段
		Tm_pre = (self.Tin_pre + self.Tout_pre) / 2
		Cp_pre = self.cp_func(Tm_pre)
		rho_pre = self.perfectGas_func(self.P_pre*10**6,self.R,Tm_pre)
		W_air = v_out * A
		deltaT_pre = self.Tout_pre - self.Tin_pre
		Q_pre = Cp_pre * rho_pre * W_air * deltaT_pre
		W = rho_pre * W_air
		#Surge Tank
		self.t_oper = 1800 #s
		V_change = W_air * self.t_oper
		V_tank = 10 *  V_change
		#fmt
		width = 50
		
		params_width = 10
		value_width = 14
		space_width = 2
		unit_width = width - params_width - value_width - 2 * space_width
		normal_digit = 4
		lower_digit = 3
		llower_digit = 2
		lllower_digit = 1
		large_digit = 8
		llarge_digit = 6
		
		type_fmt = '+{{:=^{}}}+'.format(width)
		cutting_fmt = '+{{:-^{}}}+'.format(width)
		title_fmt = '+{{:<{}}}{{:<{}}}+'.format(space_width,width-space_width)
		normalValue_fmt = '+{{0:{0}}}{{1:<{1}}}{{2:>{2}.{3}f}}{{3:>{4}}}{{0:{0}}}+'.format(space_width,params_width,value_width,normal_digit,unit_width)
		lower_fmt = '+{{0:{0}}}{{1:<{1}}}{{2:>{2}.{3}f}}{{3:>{4}}}{{0:{0}}}+'.format(space_width,params_width,value_width,lower_digit,unit_width)
		llower_fmt = '+{{0:{0}}}{{1:<{1}}}{{2:>{2}.{3}f}}{{3:>{4}}}{{0:{0}}}+'.format(space_width,params_width,value_width,llower_digit,unit_width)
		lllower_fmt = '+{{0:{0}}}{{1:<{1}}}{{2:>{2}.{3}f}}{{3:>{4}}}{{0:{0}}}+'.format(space_width,params_width,value_width,lllower_digit,unit_width)
		large_fmt = '+{{0:{0}}}{{1:<{1}}}{{2:>{2}.{3}f}}{{3:>{4}}}{{0:{0}}}+'.format(space_width,params_width,value_width,large_digit,unit_width)
		llarge_fmt = '+{{0:{0}}}{{1:<{1}}}{{2:>{2}.{3}f}}{{3:>{4}}}{{0:{0}}}+'.format(space_width,params_width,value_width,llarge_digit,unit_width)
		
		
		
		return type_fmt.format('Test Section')+'\n'\
		       +title_fmt.format('','Material Paramters')+'\n'\
		       +cutting_fmt.format('') +'\n'\
		       +lower_fmt.format('','De',self.De,'m')+'\n'\
		       +lllower_fmt.format('','L',self.L,'m')+'\n'\
		       +large_fmt.format('','A',A,'m^2')+'\n'\
		       +normalValue_fmt.format('','Rt',Rt,'K/W')+'\n'\
		       +llarge_fmt.format('','W',W,'kg/s')+'\n'\
		       +normalValue_fmt.format('','Q',Q_out,'W')+'\n'\
		       +cutting_fmt.format('') +'\n'\
		       +title_fmt.format('','Inlet Paramters')+'\n'\
		       +cutting_fmt.format('') +'\n'\
		       +normalValue_fmt.format('','Pin',self.Pin,'MPa')+'\n'\
		       +llower_fmt.format('','Tin',self.Tin,'K')+'\n'\
		       +llower_fmt.format('','Tt_in',T_w,'K')+'\n'\
		       +cutting_fmt.format('') +'\n'\
		       +title_fmt.format('','Outlet Paramters')+'\n'\
		       +cutting_fmt.format('') +'\n'\
		       +normalValue_fmt.format('','Pout',self.Pout/10**6,'MPa')+'\n'\
		       +normalValue_fmt.format('','Re_out',Re_max,'')+'\n'\
		       +normalValue_fmt.format('','K_air',self.Kair_out,'W/(m.K)')+'\n'\
		       +normalValue_fmt.format('','h_air',h_out,'W/(m.K)')+'\n'\
		       +normalValue_fmt.format('','rho_out',self.rho_out,'kg/(m^3)')+'\n'\
		       +llarge_fmt.format('','mu_out',self.mu_out,'Pa.s')+'\n'\
		       +normalValue_fmt.format('','Cp_out',self.Cp_out,'J/(kg.K)')+'\n'\
		       +llower_fmt.format('','Tout',self.Tout,'K')+'\n'\
		       +lower_fmt.format('','Pr_out',self.Pr,'')+'\n'\
		       +normalValue_fmt.format('','Nu_out',Nu_out,'')+'\n'\
		       +normalValue_fmt.format('','v_out',v_out,'m/s')+'\n'\
		       +llower_fmt.format('','Tt_out',Tt_out,'K')+'\n'\
		       +llarge_fmt.format('','V_out',W_air,'m^3/s')+'\n'\
		       +cutting_fmt.format('') +'\n'\
		       +title_fmt.format('','Difference/Mean')+'\n'\
		       +cutting_fmt.format('') +'\n'\
		       +llower_fmt.format('','deltaT',self.deltaT,'K')+'\n'\
		       +llower_fmt.format('','Tf',T_f,'K')+'\n'\
		       +normalValue_fmt.format('','deltaP',self.deltaP,'Pa')+'\n'\
		       +'\n'\
		       +type_fmt.format('PreHeater')+'\n'\
		       +'\n'\
		       +normalValue_fmt.format('','rho_pre',rho_pre,'kg/(m^3)')+'\n'\
		       +normalValue_fmt.format('','P_pre',self.P_pre,'MPa')+'\n'\
		       +llower_fmt.format('','Tin_pre',self.Tin_pre,'K')+'\n'\
		       +llower_fmt.format('','Tout_pre',self.Tout_pre,'K')+'\n'\
		       +llower_fmt.format('','Tm_pre',Tm_pre,'K')+'\n'\
		       +normalValue_fmt.format('','Q_pre',Q_pre,'W')+'\n'\
		       +'\n'\
		       +type_fmt.format('','V after half one hour')+'\n'\
		       +'\n'\
		       +llarge_fmt.format('','V_change',V_change,'m^3')+'\n'\
		       +'\n'\
		       +type_fmt.format('Surge Tank')+'\n'\
		       +'\n'\
		       +llarge_fmt.format('','V_tank',V_tank,'m^3')+'\n'\
		       +cutting_fmt.format('') +'\n'\
	


if __name__ == '__main__':
	t = Main()
	print(t)
