import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import os





def func(x, a, b, c):
    return a*x**2+b*x**1.5+c


def fit_curve_wp_a(x,y,x0,fin_a,object):
    popt, pcov = curve_fit(func, x, y)
    
    a = popt[0] 
    b = popt[1]
    c = popt[2]
    yvals = func(x,a,b,c)
    
    a_minic=-object.mu()*object.Rp**2/object.Cp()
    b_minic=(object.mu()*object.Rp**2*fin_a**0.5+object.Cp()*fin_a**(-1.5)+object.Cs()*fin_a**(-1.5))/object.Cp()
    c_minic=-object.Cs()/object.Cp()
    
    yvals_minic = a_minic*x**2+b_minic*x**1.5+c_minic
    fig=plt.figure(figsize=(32,24))
    plt.grid(True)
    plot1 = plt.plot(x, y, 's',label='original values')
    plot2 = plt.plot(x, yvals, 'b',label='{0}x^2+{1}x^1.5+({2})'.format(round(a,3),round(b,3),round(c,3)))
    plot3 = plt.plot(x, yvals_minic, 'r',linestyle="-.",label='{0}x^2+{1}x^1.5+({2})'.format(round(a_minic,3),round(b_minic,3),round(c_minic,3)))
    plt.ylabel('spin angular velocities of Pluto (mean motion)',fontsize=30)
    plt.xlabel('orbital semimajor axis (Pluto radius)',fontsize=30)
    plt.legend(loc=1,fontsize=30) #指定legend的位置右下角
    plt.title("Initial conditon({},{},{},{})".format(x0[0],x0[1],x0[2],x0[3]),fontsize=30)
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40) 
    
    plt.savefig("{}_{}_{}_{}_fit.png".format(x0[0],x0[1],x0[2],x0[3]))


    return 


def fit_curve_wp_a_ini(para,value=10):
    
    
    
    os.chdir("fitting_data_change_{0}".format(para))



    data = pd.read_csv("{0}{1}_afin17.csv".format(para,value))
    x = data['orbital semimajor axis_ini']
    y = data['pluto_spin_ini']

    popt, pcov = curve_fit(func, x, y)

    a = popt[0] 
    b = popt[1]
    c = popt[2]
    yvals = func(x,a,b,c)


    fig=plt.figure(figsize=(12,9))
    plot1 = plt.plot(x, y, 's',label='original values')
    plot2 = plt.plot(x, yvals, 'r',label='{0}x^2+{1}x^1.5+({2})'.format(round(a,3),round(b,3),round(c,3)))
    plt.ylabel('spin angular velocities of Pluto (mean motion)')
    plt.xlabel('orbital semimajor axis (Pluto radius)')
    plt.legend(loc=1)
    plt.title(r'$\alpha_c$={0},e=0'.format(value))
    #plt.show()
    plt.savefig('{0}{1}_afin17.png'.format(para,value))  
    os.chdir('..')
    
    