#coding = utf-8
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import os
import sys
sys.setrecursionlimit(10000)
import glob

class DF():
	
	def __init__(self):
		self.D = 0.00424 #m
		self.Tin = 578.7 #k
		self.q = 100000 #w/m^2
		self.path_m = ''
		self.path_heatedwalls = ''
		#self.m_mean = self.get_m_mean()
		
		
		
	def x_func(self,z):
		return z-6*self.D
	
	def xD_func(self,x):
		return x/self.D
	
	def Cp_func(self,T):
		return 1.2286*10**(-10)*T**4-5.2906*10**(-7)\
		       *T**3+0.7365*10**(-3)*T**2-1.8626*10**\
		       (-1)*T+997.5834
	
		
	def Tb_func(self,x,Tb0=600):
		f = lambda T : self.Cp_func(T)*self.m_mean
		Q = self.q*np.pi*self.D*x
		Q_inter,er = quad(f,self.Tin,Tb0)
		if Q - Q_inter > 10:
			Tb0 *= 2
			return self.Tb_func(x,Tb0)
		elif 0 < Q - Q_inter < 10:
			return Tb0
		elif Q - Q_inter < 0:
			Tb0 *= 3/4
			return self.Tb_func(x,Tb0)
	
	def k_func(self,T):
		return -2.6005*(10**-15)*(T**4)+1.5104*(10**-11)*(T**3)\
		       -3.6349*(10**-8)*(T**2)+8.6896*(10**-5)*T+0.003337
	
	def h_func(self,df):
		return self.q/(df.Tw - df.Tb)
	
	def Nu_func(self,df):
		return df.h*self.D/df.k
	
	def mu_func(self,T):
		return -2.2810*10**-18*T**4+1.303*10**-14*T**3\
		       -3.1078*10**-11*T**2+6.0519*10**-8*T+3.2333*10**-6
	
	def Re_func(self,mu):
		return self.m_mean*self.D/((0.25*np.pi*self.D**2)*mu)
	
	def Pr_func(self,mu,Cp,K):
		return mu * Cp / K
	
	def Nuc_func(self,Re,Pr,Tw,Tb,xD):
		return 0.023*Re**0.8*Pr**0.4*(Tw/Tb)**-(0.57-1.59/xD)
	
	def get_df_final(self):
		
		df = pd.read_csv(self.path_heatedwalls,skiprows=2,sep='\s+',header=None)
		
		df_sort = df.sort_values(by=2)
		
		df_sort.reset_index(drop=True, inplace=True)
	
		df_group = df_sort.groupby(by=2,as_index=False)
		n = df_group[3].count().iloc[0,1] 
		df_final = df_group.aggregate({3:lambda x: sum(x)/n})
		df_final.columns = ['z','Tw']
		df_m = pd.read_csv(self.path_m,skiprows=5,sep='\s+',header=None)
		self.m_mean = df_m.iloc[:,1].sum()/len(df_m.iloc[:,1])


		df_final['x'] = df_final['z'].apply(self.x_func)
		df_final['xD'] = df_final['x'].apply(self.xD_func)
		df_final['Tb'] = df_final.apply(lambda df:self.Tb_func(df.x),axis=1) 
		df_final['k'] = df_final['Tb'].apply(self.k_func)
		df_final['h'] = df_final.apply(self.h_func,axis=1)
		df_final['Nu'] = df_final.apply(self.Nu_func,axis=1)
		df_final['mu'] = df_final['Tb'].apply(self.mu_func)
		df_final['Cp'] = df_final['Tb'].apply(self.Cp_func)
		df_final['Re'] = df_final['mu'].apply(self.Re_func)
		df_final['Pr'] = df_final.apply(lambda df:self.Pr_func(df.mu,df.Cp,df.k),axis=1)
		df_final['Nuc'] = df_final.apply(lambda df:self.Nuc_func
								(df.Re,df.Pr,df.Tw,df.Tb,df.xD),axis=1)
		
		return df_final


if __name__ == '__main__':
	
	path_current = os.getcwd()
	os.chdir(path_current)
	#grid1
	path_m1 = glob.glob('./*1/f*/*/*')[0]
	path_h_dir1 = glob.glob('./*1/h*')[0]
	path_h_list1 = os.listdir(path_h_dir1)
	
	path_h1 = list(map(float,path_h_list1))
	path_h1.sort()
	#path_h11 = glob.glob(path_h_dir1 + '/'+str(path_h1[-2])+'/*')[0]
	path_h11 = glob.glob(path_h_dir1 + '/'+str(path_h1[-1])+'/*')[0]
	
	#grid2
	path_m2 = glob.glob('./*2/f*/*/*')[0]
	path_h_dir2 = glob.glob('./*2/h*')[0]
	path_h_list2 = os.listdir(path_h_dir2)
	
	path_h2 = list(map(float,path_h_list2))
	path_h2.sort()
	#path_h12 = glob.glob(path_h_dir2 + '/'+str(path_h2[-2])+'/*')[0]
	path_h12 = glob.glob(path_h_dir2 + '/'+str(path_h2[-1])+'/*')[0]
	
	t1 = DF()
	t1.path_m = path_m1
	t1.path_heatedwalls = path_h11
	df1 = t1.get_df_final()
	x1 = df1.loc[:,'xD']
	Nu1 = df1.loc[:,'Nu']
	
	t2 = DF()
	t2.path_m = path_m2
	t2.path_heatedwalls = path_h12
	df2 = t2.get_df_final()
	x2 = df2.loc[:,'xD']
	Nu2 = df2.loc[:,'Nu']
	
	plt.plot(x1,Nu1)
	plt.plot(x2,Nu2)
	ax = plt.gca()
	ax.set_yscale('log')
	plt.xlabel('x/D')
	plt.ylabel('Nu')
	plt.legend(['500w axial7146 radial700','555w axial7940 radial700'])
	plt.show()	

