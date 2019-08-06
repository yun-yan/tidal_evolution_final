#encoding=utf-8  
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import os
import shutil


#自定义函数 e指数形式
def func(x, a, b, c):
    return a*x**2+b*x**1.5+c

def func2(x):
    return -PlutoCharon.mu()*(PlutoCharon.Rp*x)**2/PlutoCharon.Cp()+PlutoCharon.mu()*PlutoCharon.Rp**2*17**0.5/PlutoCharon.Cp()*x**1.5-PlutoCharon.Cs()/PlutoCharon.Cp()*(10-1)


def fitting(wcvalue=2,evalue=0.5):

    #导入数据及x、y散点坐标
    data = pd.read_csv("output_a_big_data__wc{0}_e{1}_exactly.csv".format(wcvalue,evalue))
    x = data['orbital semimajor axis_ini']
    y = data['pluto_spin_ini']

    #非线性最小二乘法拟合
    guess = (-2.5,1,0)
    popt, pcov = curve_fit(func, x, y)
    
    #获取popt里面是拟合系数
    a = popt[0] 
    b = popt[1]
    c = popt[2]
    yvals1 = func(x,a,b,c) #拟合y值
    
     
    #绘图
    plot1 = plt.plot(x, y, 's',label=r'$\alpha_c$={0},e={1}'.format(wcvalue,evalue))
    plot2 = plt.plot(x, yvals1,label='{0}x^2+{1}x^1.5+({2})'.format(round(a,3),round(b,3),round(c,3)))
    plt.ylabel('spin angular velocities of Pluto (mean motion)',fontsize=80)
    plt.xlabel('orbital semimajor axis (Pluto radius)',fontsize=80)
    plt.xticks(fontsize=80)
    plt.yticks(fontsize=80)   
    plt.legend(loc=1,fontsize=30,fancybox=True, framealpha=0.1) #指定legend的位置右上角
    plt.savefig('new_no_fit')
    

  
    
    return

fig=plt.figure(figsize=(54,36))
for i in [10,20,30,40]:
    fitting(i,0)