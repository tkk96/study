import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from tkinter import *
from tkinter import messagebox
from tkinter.scrolledtext import ScrolledText
from tkinter.ttk import *
from tkinter.filedialog import asksaveasfilename,askopenfilename,askdirectory


class CalculateEquations():
	
	#Calculate functions
	def k_func(self,T):
		return -4.676*(10**-15)*(T**4)+2.257*(10**-11)*(T**3)\
		       -4.617*(10**-8)*(T**2)+9.25*(10**-5)*T+0.00217 #W/(m.k)
	
	def R_func(self,M):
		return 8.314*10**3/M
	
	def mu_func(self,T):
		return -3.245*10**-18*T**4+1.65*10**-14*T**3\
		       -3.57*10**-11*T**2+6.32*10**-8*T+2.67*10**-6
	
	def nu_func(self,mu,rho):
		return mu / rho
	
	def perfectGas_func(self,P,R,T):
		return P*10**6/(R*T)
	
	def Cp_func(self,T):
		return 2.224*10**(-10)*T**4-8.91*10**(-7)\
		       *T**3+1.218*10**(-3)*T**2-4.632*10**\
		       (-1)*T+1056
	
	def Pr_func(self,mu,Cp,K):
		return mu * Cp / K
	
	def A_func(self,D):
		return 0.25 * np.pi * D**2
	
	def DB_func(self,Re,Pr):
		return 0.023 * (Re ** 0.8) * (Pr ** 0.4)
	
	def h_func(self,Nu,k,D):
		return Nu * k / D
	
	def v_func(self,mu,Re,rho,D):
		return mu * Re / rho / D
	
	def m_func(self,rho,v,A):
		return rho * v * A
	
	def QL_func(self,Cp,m,deltaT):
		return Cp * m * deltaT
	
	def ql0_func(self,Q_L,L):
		return Q_L / 2 * np.pi / L
		
	def Tf_func(self,ql_0,L,z,Cp,m,Tf_in):
		return Tf_in + ql_0 * L / np.pi * (np.sin(np.pi*z/L) + 1) / (Cp * m) 
	
	def Tw_func(self,Tf,ql_0,z,L,D,h):
		return Tf + ql_0 * np.cos(np.pi*z/L) / (np.pi * D * h)
	
	def plot(self):
		#paramters
		ql_0 = eval(self.ql_0value.get())
		L = eval(self.Lvalue_entry.get())
		Cp = eval(self.Cpvalue.get())
		m = eval(self.mvalue.get())
		Tf_in = eval(self.Tinvalue_entry.get())
		D = eval(self.Devalue_entry.get())
		h = eval(self.hvalue.get())
		z = np.linspace(-L/2,L/2,100)
			
		#T
		Tf = self.Tf_func(ql_0,L,z,Cp,m,Tf_in)
		Tw = self.Tw_func(Tf,ql_0,z,L,D,h)
			
		#Plot
		deltaTfw = list(map(lambda x,y:x-y,Tw,Tf))
		max_DTfw = max(deltaTfw)
		max_index = deltaTfw.index(max_DTfw)
		z_max = z[max_index]
		Tmin = Tf[max_index]
		Tmax = Tw[max_index]
		text_x = 0.9 * z_max + 0.05 * L
		text_y = (Tmin + Tmax) // 2
		max_DTfw_fmt = '{:.2f}'.format(max_DTfw)
		text = r'$\leftarrow$'+r'$\Delta T_{max}$'+' = '+max_DTfw_fmt+'K'
		
		self.ax.clear()
		
		self.ax.text(text_x,text_y,text)
		self.ax.vlines(z_max,Tmin,Tmax,colors='black',ls='--')
		self.ax.plot(z,Tf,color='black',ls='-.')
		self.ax.plot(z,Tw,color='black',ls='-')
		self.ax.legend(['Tf','Tw']) 
		self.ax.set_xlabel('L (m)')
		self.ax.set_ylabel('T (K)') 
		self.ax.set_title('Cosine Heating')
		
		self.canvas.draw()
	
class ThermoPhysicalProperties(Frame,CalculateEquations):
	
	def __init__(self):
		self.root = Tk()
		self.setupGUI()
		Frame.__init__(self)
	
	def setupGUI(self):
		self.root.title('')
		self.width = 10
		self.getBestGeometry()
		self.setupNotebook()
		self.setupCanvas()
	
	def getBestGeometry(self):
		self.sw = self.root.winfo_screenwidth()
		self.sh = self.root.winfo_screenheight()
		w = 1000
		h = 600
		#窗口摆放位置
		x = (self.sw - w) / 2
		y = (self.sh - h) / 2
		self.root.geometry('%dx%d+%d+%d'%(w,h,x,y))
	
	def setupNotebook(self):
		notebook = Notebook(self.root)
		
		self.TPP_frame = Frame() #ThermoPhysicalProperties Frame
		notebook.add(self.TPP_frame,text='空气热物性')
		
		self.createTPPInput()
		self.createTPPButton()
		self.createTPPOutput()
		
		self.THC_frame = Frame() #Thermal Hydraulic Calculation Frame
		notebook.add(self.THC_frame,text='热工水力计算')
		
		self.createTHCInput()
		self.createTHCButton()
		self.createTHCOutput()
		
		notebook.pack(side='left',padx=5,pady=5,fill='both',expand=True)
		
	def createTPPInput(self):
		#Temperature
		T_label = Label(self.TPP_frame,text='温度T')
		T_label.grid(row=0,column=0,padx=5,pady=5,sticky=W)
		
		self.Tvalue_entry = Entry(self.TPP_frame,width=self.width)
		self.Tvalue_entry.insert(0,'950')
		self.Tvalue_entry.grid(row=0,column=1,padx=5,pady=5,sticky=E)
		
		Tunit_label = Label(self.TPP_frame,text='K')
		Tunit_label.grid(row=0,column=2,padx=5,pady=5,sticky=W)
		#Pressure
		P_label = Label(self.TPP_frame,text='压力P')
		P_label.grid(row=1,column=0,padx=5,pady=5,sticky=W)
		
		self.Pvalue_entry = Entry(self.TPP_frame,width=self.width)
		self.Pvalue_entry.insert(0,'0.5')
		self.Pvalue_entry.grid(row=1,column=1,padx=5,pady=5,sticky=W)
		
		Punit_label = Label(self.TPP_frame,text='Mpa')
		Punit_label.grid(row=1,column=2,padx=5,pady=5,sticky=W)
	
	def createTPPButton(self):
		#Button functions
		def clear():
			self.Mvalue.set('')
			self.kvalue.set('')
			self.Rvalue.set('')
			self.μvalue.set('')
			self.vvalue.set('')
			self.rhovalue.set('')
			self.Cpvalue.set('')
			self.Prvalue.set('')
			
			self.ax.clear()
			
			self.canvas.draw()
		
		def calculate():
			self.Mvalue.set('28.7')
			
			k = self.k_func(eval(self.Tvalue_entry.get()))
			k_fmt = '{:.4f}'
			self.kvalue.set(k_fmt.format(k))
			
			R = self.R_func(eval(self.Mvalue.get()))
			R_fmt = '{:.2f}'
			self.Rvalue.set(R_fmt.format(R))
			
			rho = self.perfectGas_func(eval(self.Pvalue_entry.get()),
		                                   eval(self.Rvalue.get()),
		                                   eval(self.Tvalue_entry.get()))
			rho_fmt = '{:.4f}'
			self.rhovalue.set(rho_fmt.format(rho))
			
			μ = self.mu_func(eval(self.Tvalue_entry.get()))
			μ_fmt = '{:.8f}'
			self.μvalue.set(μ_fmt.format(μ))
			
			v = self.nu_func(eval(self.μvalue.get()),
							eval(self.rhovalue.get()))
			v_fmt = '{:.8f}'
			self.vvalue.set(v_fmt.format(v))
			
			Cp = self.Cp_func(eval(self.Tvalue_entry.get()))
			Cp_fmt = '{:.4f}'
			self.Cpvalue.set(Cp_fmt.format(Cp))
			
			Pr = self.Pr_func(eval(self.μvalue.get()),
		                          eval(self.Cpvalue.get()),
		                          eval(self.kvalue.get()))
			Pr_fmt = '{:.3f}'
			self.Prvalue.set(Pr_fmt.format(Pr))
					
		#Clear
		Clear_button = Button(self.TPP_frame,text='清空',command=clear)
		Clear_button.grid(row=2,column=0,padx=5,pady=5,sticky=W)
		
		#Calculate
		Calcu_button = Button(self.TPP_frame,text='计算',command=calculate)
		Calcu_button.grid(row=2,column=2,padx=5,pady=5,sticky=E)
		
	def createTPPOutput(self):
		#M
		M_label = Label(self.TPP_frame,text='摩尔质量M')
		M_label.grid(row=3,column=0,padx=5,pady=5,sticky=W)
		
		self.Mvalue = StringVar()
		self.Mvalue.set('')
		Mvalue_label = Label(self.TPP_frame,textvariable=self.Mvalue,anchor='center')
		Mvalue_label.grid(row=3,column=1,padx=5,pady=5,sticky=W)
		
		Munit_label = Label(self.TPP_frame,text=' ')
		Munit_label.grid(row=3,column=2,padx=5,pady=5,sticky=W)
		
		#k
		k_label = Label(self.TPP_frame,text='导热率k')
		k_label.grid(row=4,column=0,padx=5,pady=5,sticky=W)
		
		self.kvalue = StringVar()
		self.kvalue.set('')
		kvalue_label = Label(self.TPP_frame,textvariable=self.kvalue,anchor='center')
		kvalue_label.grid(row=4,column=1,padx=5,pady=5,sticky=W)
		
		kunit_label = Label(self.TPP_frame,text='W/(m.K)')
		kunit_label.grid(row=4,column=2,padx=5,pady=5,sticky=W)
		
		#R
		R_label = Label(self.TPP_frame,text='气体常数Rg')
		R_label.grid(row=5,column=0,padx=5,pady=5,sticky=W)
		
		self.Rvalue = StringVar()
		self.Rvalue.set('')
		Rvalue_label = Label(self.TPP_frame,textvariable=self.Rvalue,anchor='center')
		Rvalue_label.grid(row=5,column=1,padx=5,pady=5,sticky=W)
		
		Runit_label = Label(self.TPP_frame,text=' ')
		Runit_label.grid(row=5,column=2,padx=5,pady=5,sticky=W)
		
		#ρ
		rho_label = Label(self.TPP_frame,text='密度ρ')
		rho_label.grid(row=6,column=0,padx=5,pady=5,sticky=W)
		
		self.rhovalue = StringVar()
		self.rhovalue.set('')
		rhovalue_label = Label(self.TPP_frame,textvariable=self.rhovalue,anchor='center')
		rhovalue_label.grid(row=6,column=1,padx=5,pady=5,sticky=W)
		
		rhounit_label = Label(self.TPP_frame,text='kg/m^3')
		rhounit_label.grid(row=6,column=2,padx=5,pady=5,sticky=W)
		
		#μ
		μ_label = Label(self.TPP_frame,text='动力粘度μ')
		μ_label.grid(row=7,column=0,padx=5,pady=5,sticky=W)
		
		self.μvalue = StringVar()
		self.μvalue.set('')
		μvalue_label = Label(self.TPP_frame,textvariable=self.μvalue,anchor='center')
		μvalue_label.grid(row=7,column=1,padx=5,pady=5,sticky=W)
		
		μunit_label = Label(self.TPP_frame,text='Pa.s')
		μunit_label.grid(row=7,column=2,padx=5,pady=5,sticky=W)
		
		#v
		v_label = Label(self.TPP_frame,text='运动粘度v')
		v_label.grid(row=8,column=0,padx=5,pady=5,sticky=W)
		
		self.vvalue = StringVar()
		self.vvalue.set('')
		vvalue_label = Label(self.TPP_frame,textvariable=self.vvalue,anchor='center')
		vvalue_label.grid(row=8,column=1,padx=5,pady=5,sticky=W)
		
		vunit_label = Label(self.TPP_frame,text='m^2/s')
		vunit_label.grid(row=8,column=2,padx=5,pady=5,sticky=W)
		
		#Cp
		Cp_label = Label(self.TPP_frame,text='定压比热Cp')
		Cp_label.grid(row=9,column=0,padx=5,pady=5,sticky=W)
		
		self.Cpvalue = StringVar()
		self.Cpvalue.set('')
		Cpvalue_label = Label(self.TPP_frame,textvariable=self.Cpvalue,anchor='center')
		Cpvalue_label.grid(row=9,column=1,padx=5,pady=5,sticky=W)
		
		Cpunit_label = Label(self.TPP_frame,text='J/(kg.K)')
		Cpunit_label.grid(row=9,column=2,padx=5,pady=5,sticky=W)
		
		#Pr
		Pr_label = Label(self.TPP_frame,text='普朗特数Pr')
		Pr_label.grid(row=10,column=0,padx=5,pady=5,sticky=W)
		
		self.Prvalue = StringVar()
		self.Prvalue.set('')
		Prvalue_label = Label(self.TPP_frame,textvariable=self.Prvalue,anchor='center')
		Prvalue_label.grid(row=10,column=1,padx=5,pady=5,sticky=W)
		
		Prunit_label = Label(self.TPP_frame,text=' ')
		Prunit_label.grid(row=10,column=2,padx=5,pady=5,sticky=W)
	
	def createTHCInput(self):
		#D
		D_label = Label(self.THC_frame,text='直径D')
		D_label.grid(row=0,column=0,padx=5,pady=5,sticky=W)
		
		self.Devalue_entry = Entry(self.THC_frame,width=self.width)
		self.Devalue_entry.insert(0,'0.01')
		self.Devalue_entry.grid(row=0,column=1,padx=5,pady=5,sticky=W)
		
		Deunit_label = Label(self.THC_frame,text='m')
		Deunit_label.grid(row=0,column=2,padx=5,pady=5,sticky=W)
		
		#L
		L_label = Label(self.THC_frame,text='管长L')
		L_label.grid(row=1,column=0,padx=5,pady=5,sticky=W)
		
		self.Lvalue_entry = Entry(self.THC_frame,width=self.width)
		self.Lvalue_entry.insert(0,'1')
		self.Lvalue_entry.grid(row=1,column=1,padx=5,pady=5,sticky=W)
		
		Lunit_label = Label(self.THC_frame,text='m')
		Lunit_label.grid(row=1,column=2,padx=5,pady=5,sticky=W)
		
		#Tin
		Tin_label = Label(self.THC_frame,text='进口温度Tin')
		Tin_label.grid(row=2,column=0,padx=5,pady=5,sticky=W)
		
		self.Tinvalue_entry = Entry(self.THC_frame,width=self.width)
		self.Tinvalue_entry.insert(0,'800')
		self.Tinvalue_entry.grid(row=2,column=1,padx=5,pady=5,sticky=W)
		
		Tinunit_label = Label(self.THC_frame,text='K')
		Tinunit_label.grid(row=2,column=2,padx=5,pady=5,sticky=W)
		
		#Tout
		Tout_label = Label(self.THC_frame,text='出口温度Tout')
		Tout_label.grid(row=3,column=0,padx=5,pady=5,sticky=W)
		
		self.Toutvalue_entry = Entry(self.THC_frame,width=self.width)
		self.Toutvalue_entry.insert(0,'1100')
		self.Toutvalue_entry.grid(row=3,column=1,padx=5,pady=5,sticky=W)
		
		Toutunit_label = Label(self.THC_frame,text='K')
		Toutunit_label.grid(row=3,column=2,padx=5,pady=5,sticky=W)
		
		#Re
		Re_label = Label(self.THC_frame,text='Re')
		Re_label.grid(row=4,column=0,padx=5,pady=5,sticky=W)
		
		self.Revalue_entry = Entry(self.THC_frame,width=self.width)
		self.Revalue_entry.insert(0,'10000')
		self.Revalue_entry.grid(row=4,column=1,padx=5,pady=5,sticky=W)
		
		Reunit_label = Label(self.THC_frame,text='')
		Reunit_label.grid(row=4,column=2,padx=5,pady=5,sticky=W)
		
		#Cos Heating
		self.cos_heating_var = BooleanVar()
		self.cos_heating_var.set(True)
		Cos_ckb = Checkbutton(self.THC_frame,text='余弦加热',variable=self.cos_heating_var)
		Cos_ckb.grid(row=5,column=0,padx=5,pady=5,sticky=W)
	
	def createTHCButton(self):
		#Button functions
		def clear():
			self.Avalue.set('')
			self.Nuvalue.set('')
			self.hvalue.set('')
			self.Vmeanvalue.set('')
			self.mvalue.set('')
			self.Q_allvalue.set('')
			self.ql_0value.set('')
			#self.Prvalue.set('')
			
			self.ax.clear()
			
			self.canvas.draw()
		
		def calculate():
			self.De = eval(self.Devalue_entry.get())
			A = self.A_func(self.De)
			A_fmt = '{:.8f}'
			self.Avalue.set(A_fmt.format(A))
			
			Nu = self.DB_func(eval(self.Revalue_entry.get()),
								eval(self.Prvalue.get()))
			Nu_fmt = '{:.2f}'
			self.Nuvalue.set(Nu_fmt.format(Nu))
			
			self.h = self.h_func(eval(self.Nuvalue.get()),
								eval(self.kvalue.get()),
								eval(self.Devalue_entry.get()))
			h_fmt = '{:.2f}'
			self.hvalue.set(h_fmt.format(self.h))
			
			Vmean = self.v_func(eval(self.μvalue.get()),
		                                   eval(self.Revalue_entry.get()),
		                                   eval(self.rhovalue.get()),
		                                   eval(self.Devalue_entry.get()))
			Vmean_fmt = '{:.4f}'
			self.Vmeanvalue.set(Vmean_fmt.format(Vmean))
			
			m = self.m_func(eval(self.rhovalue.get()),
								eval(self.Vmeanvalue.get()),
								eval(self.Avalue.get()),)
			m_fmt = '{:.6f}'
			self.mvalue.set(m_fmt.format(m))
			
			deltaT = eval(self.Toutvalue_entry.get()) - eval(self.Tinvalue_entry.get())
			Q_all = self.QL_func(eval(self.Cpvalue.get()),
									eval(self.mvalue.get()),
									deltaT)
			Q_all_fmt = '{:.4f}'
			self.Q_allvalue.set(Q_all_fmt.format(Q_all))
			
			ql_0 = self.ql0_func(eval(self.Q_allvalue.get()),
									eval(self.Lvalue_entry.get()))
			ql_0_fmt = '{:.4f}'
			self.ql_0value.set(ql_0_fmt.format(ql_0))
			
		#Clear
		Clear_button = Button(self.THC_frame,text='清空',command=clear)
		Clear_button.grid(row=6,column=0,padx=5,pady=5,sticky=W)
		
		#Calculate
		Calcu_button = Button(self.THC_frame,text='计算',command=calculate)
		Calcu_button.grid(row=6,column=1,padx=5,pady=5,sticky=E)
		
		#Plot
		Plot_button = Button(self.THC_frame,text='温度分布',command=self.plot)
		Plot_button.grid(row=6,column=2,padx=5,pady=5,sticky=E)
	
	def createTHCOutput(self):
		#A
		A_label = Label(self.THC_frame,text='管子截面积A')
		A_label.grid(row=7,column=0,padx=5,pady=5,sticky=W)
		
		self.Avalue = StringVar()
		self.Avalue.set('')
		Avalue_label = Label(self.THC_frame,textvariable=self.Avalue,anchor='center')
		Avalue_label.grid(row=7,column=1,padx=5,pady=5,sticky=W)
		
		Aunit_label = Label(self.THC_frame,text='m^2')
		Aunit_label.grid(row=7,column=2,padx=5,pady=5,sticky=W)
		
		#Nu
		Nu_label = Label(self.THC_frame,text='努塞尔特数Nu')
		Nu_label.grid(row=8,column=0,padx=5,pady=5,sticky=W)
		
		self.Nuvalue = StringVar()
		self.Nuvalue.set('')
		Nuvalue_label = Label(self.THC_frame,textvariable=self.Nuvalue,anchor='center')
		Nuvalue_label.grid(row=8,column=1,padx=5,pady=5,sticky=W)
		
		Nuunit_label = Label(self.THC_frame,text='')
		Nuunit_label.grid(row=8,column=2,padx=5,pady=5,sticky=W)
		
		#h
		h_label = Label(self.THC_frame,text='换热系数h')
		h_label.grid(row=9,column=0,padx=5,pady=5,sticky=W)
		
		self.hvalue = StringVar()
		self.hvalue.set('')
		hvalue_label = Label(self.THC_frame,textvariable=self.hvalue,anchor='center')
		hvalue_label.grid(row=9,column=1,padx=5,pady=5,sticky=W)
		
		hunit_label = Label(self.THC_frame,text='W/(m.K)')
		hunit_label.grid(row=9,column=2,padx=5,pady=5,sticky=W)
		
		#Vmean
		Vmean_label = Label(self.THC_frame,text='流速v')
		Vmean_label.grid(row=10,column=0,padx=5,pady=5,sticky=W)
		
		self.Vmeanvalue = StringVar()
		self.Vmeanvalue.set('')
		Vmeanvalue_label = Label(self.THC_frame,textvariable=self.Vmeanvalue,anchor='center')
		Vmeanvalue_label.grid(row=10,column=1,padx=5,pady=5,sticky=W)
		
		Vmeanunit_label = Label(self.THC_frame,text='m/s')
		Vmeanunit_label.grid(row=10,column=2,padx=5,pady=5,sticky=W)
		
		#m
		m_label = Label(self.THC_frame,text='质量流量m')
		m_label.grid(row=11,column=0,padx=5,pady=5,sticky=W)
		
		self.mvalue = StringVar()
		self.mvalue.set('')
		mvalue_label = Label(self.THC_frame,textvariable=self.mvalue,anchor='center')
		mvalue_label.grid(row=11,column=1,padx=5,pady=5,sticky=W)
		
		munit_label = Label(self.THC_frame,text='kg/s')
		munit_label.grid(row=11,column=2,padx=5,pady=5,sticky=W)
		
		#Q_all
		Q_all_label = Label(self.THC_frame,text='总功率Q')
		Q_all_label.grid(row=12,column=0,padx=5,pady=5,sticky=W)
		
		self.Q_allvalue = StringVar()
		self.Q_allvalue.set('')
		Q_allvalue_label = Label(self.THC_frame,textvariable=self.Q_allvalue,anchor='center')
		Q_allvalue_label.grid(row=12,column=1,padx=5,pady=5,sticky=W)
		
		Q_allunit_label = Label(self.THC_frame,text='W')
		Q_allunit_label.grid(row=12,column=2,padx=5,pady=5,sticky=W)
		
		#ql_0
		ql_0_label = Label(self.THC_frame,text='线功率密度ql0')
		ql_0_label.grid(row=13,column=0,padx=5,pady=5,sticky=W)
		
		self.ql_0value = StringVar()
		self.ql_0value.set('')
		ql_0value_label = Label(self.THC_frame,textvariable=self.ql_0value,anchor='center')
		ql_0value_label.grid(row=13,column=1,padx=5,pady=5,sticky=W)
		
		ql_0unit_label = Label(self.THC_frame,text='W')
		ql_0unit_label.grid(row=13,column=2,padx=5,pady=5,sticky=W)
		
	def setupCanvas(self):
		#Canvas
		self.f = Figure()
		self.ax = self.f.add_subplot(111)
		self.canvas = FigureCanvasTkAgg(self.f,master=self.root)
		self.canvas.get_tk_widget().pack(side='top',padx=5,pady=5,fill='both',expand=True)
		#Toolbal
		toolbal = NavigationToolbar2Tk(self.canvas,self.root)
		toolbal.pack(side='bottom',padx=5,pady=5,fill='both',expand=True)
		toolbal.update()
		
if __name__ == '__main__':
	
	app = ThermoPhysicalProperties()
	app.mainloop()
