#encoding=utf-8  
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import os
import shutil

wcvalue=2
evalue=0.5


#导入数据及x、y散点坐标
data = pd.read_csv("output_a_big_data__wc{0}_e{1}_exactly.csv".format(wcvalue,evalue))
x = data['orbital semimajor axis_ini']
y = data['pluto_spin_ini']


 
#自定义函数 e指数形式
def func(x, a, b, c):
    return a*x**2+b*x**1.5+c

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
plot2 = plt.plot(x, yvals1, 'r',label='{0}x^2+{1}x^1.5+({2})'.format(round(a,3),round(b,3),round(c,3)))
#plot2 = plt.plot(x, yvals2, 'r',label='{0}*x**2+{1}*x**1.5+{2}'.format(func2(1)-round(np.real(func2(1)-func2(-1)),3)-func2(0),round(np.real(func2(1)-func2(-1)),3),func2(0)))
plt.ylabel('spin angular velocities of Pluto (mean motion)')
plt.xlabel('orbital semimajor axis (Pluto radius)')
plt.legend(loc=1) #指定legend的位置右下角
plt.title(r'$\alpha_c$={0},e={1}'.format(wcvalue,evalue))
#plt.show()
plt.savefig('wc={0}_e={1}.png'.format(wcvalue,evalue))

try:
    file_dir="fitting_ini_data_re"
    os.mkdir(file_dir)
except:
    pass

shutil.copy("output_a_big_data__wc{0}_e{1}_exactly.csv".format(wcvalue,evalue),file_dir)
os.remove("output_a_big_data__wc{0}_e{1}_exactly.csv".format(wcvalue,evalue))
shutil.copy("wc={0}_e={1}.png".format(wcvalue,evalue) ,file_dir)
os.remove("wc={0}_e={1}.png".format(wcvalue,evalue))