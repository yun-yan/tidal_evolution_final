#This file is to define all variable 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import ode
from TenTai import TenTai 
import run_ode
from plot_evolution import plot
import os
import curve_fit


#=============variable===================================================

G=6.674e-11

#   Earth and Moon
EARTH_MASS=59760e20;        #kg
MOON_MASS=737e20;           #kg
EARTH_RADIUS=6378e3;        #meter
MOON_RADIUS=1738e3;         #meter
EARTH_DELTAT=600;           #s
MOON_DELTAT=600;            #s
EARTH_LOVE_NUMBER=0.299;
MOON_LOVE_NUMBER=0.03;
EARTH_DISSIPATION_FUNCTION=12;
MOON_DISSIPATION_FUNCTION=27;


#  Pluto and caron
PLUTO_MASS=870.3/G*1e9;     #kg
CHARON_MASS=101.4/G*1e9;    #kg
PLUTO_RADIUS=1153e3;        #meter
CHARON_RADIUS=606e3;        #meter
PLUTO_DELTAT=600;           #s
CHARON_DELTAT=600;          #s
PLUTO_LOVE_NUMBER=0.058;
CHARON_LOVE_NUMBER=0.006;
PLUTO_DISSIPATION_FUNCTION=100;
CHARON_DISSIPATION_FUNCTION=100;

#------------------------------alteration-------------------------------------

#   Earth and Moon

MOON_ADt=1;
MOON_AQ=0.5;            
EarthMoon_Initial= (26,0,60,0)   # Earth spin/ Moon spin/ distance/ eccentricity 
EarthMoon_start= 1e-2            # year
EarthMoon_end= 1e12              # year

#  Pluto and Charon

CHARON_ADt=10*2;
CHARON_AQ=1.15; 
PlutoCharon_Initial= (5.5,2,4,0.1)   # Pluto spin/ Charon spin/ distance/ eccentricity 
PlutoCharon_start= 1e-2            # year
PlutoCharon_end= 1e8              # year


#==================================================================================

EarthMoon=TenTai(Mp=EARTH_MASS,Ms=MOON_MASS,
                   Rp=EARTH_RADIUS,Rs=MOON_RADIUS,
                   Dtp=EARTH_DELTAT,
                   k2p=EARTH_LOVE_NUMBER,k2s=MOON_LOVE_NUMBER,
                   Qp=EARTH_DISSIPATION_FUNCTION,
                   ADt=MOON_ADt,
                   AQ=MOON_AQ                                      
                   )


PlutoCharon=TenTai(Mp=PLUTO_MASS,Ms=CHARON_MASS,
                   Rp=PLUTO_RADIUS,Rs=CHARON_RADIUS,
                   Dtp=PLUTO_DELTAT,
                   k2p=PLUTO_LOVE_NUMBER,k2s=CHARON_LOVE_NUMBER,
                   Qp=PLUTO_DISSIPATION_FUNCTION,
                   ADt=CHARON_ADt,
                   AQ=CHARON_AQ
                   )



#====================================get_path==================================

PATH=os.getcwd()

#===============================================================================


#def figure_c(x0=(2,2,2,0)):
#    
#    PlutoCharon_start= 1e-2            # year
#    PlutoCharon_end= 5e7              # year
#       
#    #fig=plt.figure(figsize=(32,24))
#    t_xs=run_ode.run_Dt(x0,PlutoCharon_start,PlutoCharon_end,PlutoCharon.Dt_model_exactly)
#    plot(t_xs,x0,'Pluto','Charon',4)
#    
#    os.chdir(PATH)
#    
#    return


#fig=plt.figure(figsize=(32,24))
#for CHARON_ADt in [5,10,20]:
#    PlutoCharon=TenTai(Mp=PLUTO_MASS,Ms=CHARON_MASS,
#                   Rp=PLUTO_RADIUS,Rs=CHARON_RADIUS,
#                   Dtp=PLUTO_DELTAT,
#                   k2p=PLUTO_LOVE_NUMBER,k2s=CHARON_LOVE_NUMBER,
#                   Qp=PLUTO_DISSIPATION_FUNCTION,
#                   ADt=CHARON_ADt,
#                   AQ=CHARON_AQ
#                   )
#    figure_c((5.3,2,3.9,0.2))

def data(x0=(2,2,2,0)):
    
    PlutoCharon_start= 1e-2            # year
    PlutoCharon_end= 0.5e8              # year
       
    
    t_xs=run_ode.run_Dt_find_t(x0,PlutoCharon_start,PlutoCharon_end,PlutoCharon.Dt_model_exactly)
    np.save('test_1',t_xs)
    
    
    return t_xs


def figure(x0=(2,2,2,0)):
       
#    fig=plt.figure(figsize=(32,24))
#    t_xs=run_ode.run_Dt(x0,PlutoCharon_start,PlutoCharon_end,PlutoCharon.Dt_model_exactly)
#    plot(t_xs,x0,'Pluto','Charon',4)
    

#    t_xs=run_ode.run_Dt(x0,PlutoCharon_start,PlutoCharon_end,PlutoCharon.Q_model_exactly_continuous)
#    plot(t_xs,x0,'Pluto','Charon',4)
    
    t_xs=run_ode.run_Q_exactly(x0,PlutoCharon_start,PlutoCharon_end,PlutoCharon.Q_model_exactly)
    for i in range(len(t_xs[:,1])):
        if t_xs[i,0]> 2000:
            box = np.ones(300)/300
            break
    a = np.convolve(t_xs[i:,2], box)
    a=a[299:]
    t_xs[i:,2]=a
    
    
    for i in range(len(t_xs[:,1])):
        if t_xs[i,0]> 5e7:
            break
    t_xs=t_xs[:i]
    plot(t_xs,x0,'Pluto','Charon',4)
    
    
    #t_xs=run_ode.run_Dt(x0,PlutoCharon_start,PlutoCharon_end,PlutoCharon.Dt_model)
    #plot(t_xs,x0,'Pluto','Charon',4)

#    t_xs=run_ode.run_Q_exactly(x0,PlutoCharon_start,PlutoCharon_end,PlutoCharon.Q_model_exactly)
#    plot(t_xs,x0,'Pluto','Charon',4)
   
    #t_xs=run_ode.run_Q(x0,PlutoCharon_start,PlutoCharon_end,PlutoCharon.Q_model)
    #plot(t_xs,x0,'Pluto','Charon',4)
#    
    
    #t_xs=run_ode.run_Dt(x0,EarthMoon_start,EarthMoon_end,EarthMoon.Dt_model_exactly)
    #plot(t_xs,x0,'Earth','Moon',4)
    
    
    os.chdir(PATH)
    
    return

def figure_c(x0=(2,2,2,0)):
    
    PlutoCharon_start= 1e-2            # year
    PlutoCharon_end= 5e7              # year
       
    #fig=plt.figure(figsize=(32,24))
    t_xs=run_ode.run_Dt(x0,PlutoCharon_start,PlutoCharon_end,PlutoCharon.Dt_model_exactly)
    plot(t_xs,x0,'Pluto','Charon',4)
    
    os.chdir(PATH)
    
    return

def figure_i(x0=(2,2,2,0,0)):
    
    PlutoCharon_start= 1e-2            # year
    PlutoCharon_end= 0.5e8              # year
       
#    fig=plt.figure(figsize=(32,36))
    t_xs=run_ode.run_Dt(x0,PlutoCharon_start,PlutoCharon_end,PlutoCharon.i_Dt_model)
    t_xs[:,5]=t_xs[:,5]*180/np.pi
    plot(t_xs,x0,'Pluto','Charon',5)
    
    os.chdir(PATH)
    
    return


#fig=plt.figure(figsize=(32,24))
#fig=plt.figure(figsize=(32,36))
#for i in [0,0.1,0.2,0.3]:
#for i in [0,np.pi/6,np.pi/3,np.pi/2]:
#    figure_i((5.3,2,3.9,i,np.pi/3))


fig=plt.figure(figsize=(32,24))
figure_c((5.3,2,3.9,0.1))
figure_c((5.3,2,2.9,0.1))
figure_c((5.3,2,3.9,0.4))
figure_c((3.3,2,3.9,0.1))
figure_c((5.3,4,3.9,0.1))


def figure_wp_to_a(x0=(2,0,2,0),fit=True):
    
    PlutoCharon_start= 1e-2            # year
    PlutoCharon_end= 0.5e8              # year
    
    t_xs=run_ode.run_Dt(x0,PlutoCharon_start,PlutoCharon_end,PlutoCharon.Dt_model_exactly)
    
    
    fig=plt.figure(figsize=(32,24))
    plt.clf
    
    plt.scatter(t_xs[:,3],t_xs[:,1]) 
    plt.grid(True)
    # label axes
    plt.ylabel('spin angular velocities of Pluto (mean motion)',fontsize=30)
    plt.xlabel('orbital semimajor axis (Pluto radius)',fontsize=30)
    plt.title('Initial condition({},{},{},{})'.format(x0[0],x0[1],x0[2],x0[3]),fontsize=30)
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)    
    
    try:
        file_dir="figure_wp_to_a"
        os.mkdir(file_dir)
    except:
        pass
    os.chdir(file_dir)
    
    plt.savefig("{}_{}_{}_{}.png".format(x0[0],x0[1],x0[2],x0[3]) )
    
    if fit==True:
        curve_fit.fit_curve_wp_a(t_xs[:,3],t_xs[:,1],x0,t_xs[:,3][-1],PlutoCharon)
    else:
        pass
    
    os.chdir(PATH)
    
    return 

def figure_max_to_e():
    
    PlutoCharon_start= 1e-2            # year
    PlutoCharon_end= 0.5e8              # year
    
    e_initial_list=[i/100 for i in range(70)]
    
    
    max_result=[0,0]
    
    for i in e_initial_list:
        t_xs=run_ode.run_Dt((5.3,2,3.9,i),PlutoCharon_start,PlutoCharon_end,PlutoCharon.Dt_model_exactly)
        max_result=np.row_stack((max_result,[i,np.amax(t_xs[:,3])]))
        
        
    r'$\alpha_c$={0},e={1}'
#    fig=plt.figure(figsize=(32,24))
    plt.ylabel('max X',fontsize=30)
    plt.xlabel('initial eccentricity',fontsize=30)
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)    
    
    plt.scatter(max_result[1:,0],max_result[1:,1],label=r'$A\Delta t$={0}'.format(CHARON_ADt)) 
    plt.legend(loc=1,fontsize=40)
    try:
        file_dir="figure_max_to_e"
        os.mkdir(file_dir)
    except:
        pass
    os.chdir(file_dir)
    
    plt.savefig("{0}.png".format(CHARON_ADt))
    
     
    os.chdir(PATH)

    return 

#fig=plt.figure(figsize=(32,24))

#
#for CHARON_ADt in [6,8,10]:
#
#    PlutoCharon.ADt=CHARON_ADt
#
#    figure_max_to_e()

def figure_e_to_a(x0=(2,2,2,0.24),fit=False):
    
    PlutoCharon_start= 1e-2            # year
    PlutoCharon_end= 1e8              # year
    
    t_xs=run_ode.run_Dt_find_t(x0,PlutoCharon_start,PlutoCharon_end,PlutoCharon.Dt_model)
    
    
    fig=plt.figure(figsize=(32,24))
    plt.clf
    
    plt.scatter(t_xs[:,3],t_xs[:,4]) 
    plt.grid(True)
    # label axes
    plt.ylabel('orbital eccentricity',fontsize=30)
    plt.xlabel('orbital semimajor axis (Pluto radius)',fontsize=30)
    plt.title('Initial condition({},{},{},{})'.format(x0[0],x0[1],x0[2],x0[3]),fontsize=30)
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)    
    
    try:
        file_dir="figure_e_to_a"
        os.mkdir(file_dir)
    except:
        pass
    os.chdir(file_dir)
    
    plt.savefig("{}_{}_{}_{}.png".format(x0[0],x0[1],x0[2],x0[3]) )
    
    if fit==True:
        curve_fit.fit_curve_wp_a(t_xs[:,3],t_xs[:,1],x0)
    else:
        pass
    
    os.chdir(PATH)
    
    
    return 


def figure_a_to_e2(x0=(2,2,2,0.24),fit=False):
    
    PlutoCharon_start= 1e-2            # year
    PlutoCharon_end= 1e8              # year
    
    t_xs=run_ode.run_Dt_find_t(x0,PlutoCharon_start,PlutoCharon_end,PlutoCharon.Dt_model)
    #t_xs=run_ode.run_Dt(x0,PlutoCharon_start,PlutoCharon_end,PlutoCharon.Dt_model)
    x0_o=x0[:-1]+(0,)
    t_xs_o=run_ode.run_Dt_find_t(x0_o,PlutoCharon_start,PlutoCharon_end,PlutoCharon.Dt_model)
    
    
    fig=plt.figure(figsize=(32,24))
    plt.clf
    
    
    #plt.scatter(t_xs[:,4]**2,t_xs[:,3]) 
    plt.scatter(t_xs[:,4]-t_xs_o[:,4],t_xs[:,3]-t_xs_o[:,3]) 
    #plt.scatter(t_xs[:,0],(t_xs[:,3]-t_xs_o[:,3])**2-200*(t_xs[:,4]-t_xs_o[:,4])**4) 
    
    plt.grid(True)
    # label axes
    plt.xlabel('orbital eccentricity',fontsize=30)
    plt.ylabel('orbital semimajor axis (Pluto radius)',fontsize=30)
    plt.title('Initial condition({},{},{},{})'.format(x0[0],x0[1],x0[2],x0[3]),fontsize=30)
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)    
    
    try:
        file_dir="figure_a_to_e2"
        os.mkdir(file_dir)
    except:
        pass
    os.chdir(file_dir)
    
    plt.savefig("{}_{}_{}_{}.png".format(x0[0],x0[1],x0[2],x0[3]) )
    
    if fit==True:
        curve_fit.fit_curve_wp_a(t_xs[:,3],t_xs[:,1],x0)
    else:
        pass
    
    os.chdir(PATH)
    
    
    return 

#figure_a_to_e2((5,2,5.2,0.34))
#figure_a_to_e2((8,2,11.7,0.5))
#figure_a_to_e2((8.5,2,9.5,0.01))

def write_ini_wp_to_a(wcc=2,ee=0):
    
    
    
    log_data=np.array(["pluto_spin_ini","charon_spin_ini","orbital semimajor axis_ini","eccentricity_ini"
                   ,"pluto_spin_fin","charon_spin_fin","orbital semimajor axis_fin","eccentricity_fin"])

    
    
    wp_list=[i/10 for i in range(150,130,-1)]
    a_list=[i/10 for i in range(163,158,-1)]
    
       
    
    for j in a_list:
        for k in wp_list:
    #        try:
            t_xs=run_ode.run_Dt((k,wcc,j,ee),PlutoCharon_start,PlutoCharon_end,PlutoCharon.Dt_model_exactly)
    #        except:
    #            pass
            if 16.9<t_xs[:,3][-1]<17.1:
                try:
                    log_data=np.row_stack((log_data,[k,wcc,j,ee,t_xs[:,1][-1],t_xs[:,2][-1],t_xs[:,3][-1],t_xs[:,4][-1]]))
        
                    out_put_data=pd.DataFrame(log_data[1:],columns=log_data[0])
                    out_put_data.to_csv('output_a_big_data__wc{0}_e{1}_exactly.csv'.format(wcc,ee),index=False)
                
                    print("sucess:{0},{1},{2},{3}".format(k,wcc,j,ee))
                except:
                    pass
            else:
                print("fail:{0},{1},{2},{3}".format(k,wcc,j,ee))
                   
    
    

    return

def write_change_e(wp=5):
    
    
    log_data=np.array(["pluto_spin_ini","charon_spin_ini","orbital semimajor axis_ini","eccentricity_ini"
                   ,"pluto_spin_fin","charon_spin_fin","orbital semimajor axis_fin","eccentricity_fin"])
       
    a_list=[i/10 for i in range(201,11,-1)]
    e_list=[i/40 for i in range(0,40)]
    
    for j in a_list:
        for k in e_list:
    #        try:
            t_xs=run_ode.run_Dt((wp,2,j,k),PlutoCharon_start,PlutoCharon_end,PlutoCharon.Dt_model)
    #        except:
    #            pass
            if 16.9<t_xs[:,3][-1]<17.1:
                try:
                    log_data=np.row_stack((log_data,[wp,2,j,k,t_xs[:,1][-1],t_xs[:,2][-1],t_xs[:,3][-1],t_xs[:,4][-1]]))
                    print("sucess:{0},{1},{2},{3}".format(wp,2,j,k))
                except:
                    pass
            else:
                print("fail:{0},{1},{2},{3}".format(wp,2,j,k))
                   
    
    out_put_data=pd.DataFrame(log_data[1:],columns=log_data[0])
    out_put_data.to_csv('output_a_big_data_wp{}_wc2_exactly.csv'.format(wp),index=False)

    
    return

