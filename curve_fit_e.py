#encoding=utf-8  
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import os

os.chdir("fitting_e_to_a")


wpvalue=8.5


#导入数据及x、y散点坐标
data = pd.read_csv("wp_{0}wc_2.csv".format(wpvalue))
x = data['orbital semimajor axis_ini']
y = data['eccentricity_ini']


 
#自定义函数 e指数形式
def func(x, a, b, c):
    return a*x**0.5+b*x**1+c

def func2(x):
    return -PlutoCharon.mu()*(PlutoCharon.Rp*x)**2/PlutoCharon.Cp()+PlutoCharon.mu()*PlutoCharon.Rp**2*17**0.5/PlutoCharon.Cp()*x**1.5-PlutoCharon.Cs()/PlutoCharon.Cp()*(10-1)


#非线性最小二乘法拟合
guess = (-2.5,1,0)
popt, pcov = curve_fit(func, x, y)

#获取popt里面是拟合系数
a = popt[0] 
b = popt[1]
c = popt[2]
yvals1 = func(x,a,b,c) #拟合y值
print ('系数a:', a)
print ('系数b:', b)
print ('系数c:', c)

yvals2 = func2(x) 
 
#绘图
fig=plt.figure(figsize=(12,9))
plot1 = plt.plot(x, y, 's',label='original values')
#plot2 = plt.plot(x, yvals1, 'r',label='{0}x^2+{1}x^1.5+({2})'.format(round(a,3),round(b,3),round(c,3)))
#plot2 = plt.plot(x, yvals2, 'r',label='{0}*x**2+{1}*x**1.5+{2}'.format(func2(1)-round(np.real(func2(1)-func2(-1)),3)-func2(0),round(np.real(func2(1)-func2(-1)),3),func2(0)))
plt.ylabel('orbital eccentricity')
plt.xlabel('orbital semimajor axis (Pluto radius)')
plt.legend(loc=1) #指定legend的位置右下角
plt.title(r'$\alpha_p={0},\alpha_c=2$'.format(wpvalue))
#plt.show()
plt.savefig('wp_{0}wc_2.csv.png'.format(wpvalue))
