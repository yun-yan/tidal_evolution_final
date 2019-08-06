import scipy.integrate
import numpy as np

import run_ode



#z=((1+0.1*np.cos(x))/(1-0.1^2))

f = lambda x: ((1+0.1*np.cos(x))/(1-0.1**2))**6
i=scipy.integrate.quad(f,0,2*np.pi)[0]
print(i)

'''

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
MOON_AQ=1;            
EarthMoon_Initial= (26,0,60,0)   # Earth spin/ Moon spin/ distance/ eccentricity 
EarthMoon_start= 1e-2            # year
EarthMoon_end= 1e12              # year

#  Pluto and Charon

CHARON_ADt=10;
CHARON_AQ=1.15; 
PlutoCharon_Initial= (5.5,2,4,0.1)   # Pluto spin/ Charon spin/ distance/ eccentricity 
PlutoCharon_start= 1e-2            # year
PlutoCharon_end= 1e8              # year


Mp=PLUTO_MASS
Ms=CHARON_MASS
Rp=PLUTO_RADIUS
Rs=CHARON_RADIUS
Dtp=PLUTO_DELTAT
k2p=PLUTO_LOVE_NUMBER
k2s=CHARON_LOVE_NUMBER
Qp=PLUTO_DISSIPATION_FUNCTION
ADt=CHARON_ADt
AQ=CHARON_AQ



Cs=0.4*Ms*Rs**2
Cp=0.328*Mp*Rp**2
Qs=Qp/AQ*(k2s/k2p)*(Mp/Ms)**2*(Rs/Rp)**5


def Q_model_exactly(t,x):
    Eqs= np.zeros((4)) 
    
#        feight = lambda x: ((1+x[3]*np.cos(x))/(1-x[3]**2))**8
#        fsix = lambda x: ((1+x[3]*np.cos(x))/(1-x[3]**2))**6
#        fthree = lambda x: ((1+x[3]*np.cos(x))/(1-x[3]**2))**3
    
#        
#        feight = lambda x: ((1+0.1*np.cos(x))/(1-0.1**2))**8
#        fsix = lambda x: ((1+0.1*np.cos(x))/(1-0.1**2))**6
#        fthree = lambda x: ((1+0.1*np.cos(x))/(1-0.1**2))**3
#        
#        avsix=scipy.integrate.quad(fsix,0,2*np.pi)
#        avthree=scipy.integrate.quad(fthree,0,2*np.pi)
#        aveight=scipy.integrate.quad(feight,0,2*np.pi)
    
    
    
    avsix=1
    avthree=2
    aveight=3
    
    
    
    
    Eqs[0]= (-3*(G)/((G*(Mp+Ms)/(x[2]*Rp)**3)**(0.5))/(2*Cp*x[2]**6)*k2p/Qp*Ms**2*(Rp**-1)*(((x[0]-1)*avsix)/(x[0]-avthree))
            +3/2*x[0]*(3*(G*(Mp+Ms)/(x[2]*Rp)**3)**(0.5)*k2p/Qp*(Ms/Mp)/x[2]**5*((((x[0]-1)*aveight)/(x[0]-avthree))+AQ*(((x[0]-1)*aveight)/(x[0]-avthree)))))*31536000
    Eqs[1]=(-3*(G)/((G*(Mp+Ms)/(x[2]*Rp)**3)**(0.5))/(2*Cs()*x[2]**6)*k2s/Qs*Mp**2*(Rs**5/Rp**6)*(((x[0]-1)*avsix)/(x[0]-avthree))
            +3/2*x[1]*(3*(G*(Mp+Ms)/(x[2]*Rp)**3)**(0.5)*k2p/Qp*(Ms/Mp)/x[2]**5*((((x[0]-1)*aveight)/(x[0]-avthree))+AQ*(((x[0]-1)*aveight)/(x[0]-avthree)))))*31536000
    Eqs[2]= x[2]*(3*(G*(Mp+Ms)/(x[2]*Rp)**3)**(0.5)*k2p/Qp*(Ms/Mp)/x[2]**5*((((x[0]-1)*aveight)/(x[0]-avthree))+AQ*(((x[0]-1)*aveight)/(x[0]-avthree))))*31536000
    Eqs[3]= x[3]*((G*(Mp+Ms)/(x[2]*Rp)**3)**(0.5)*k2p/Qp*(Ms/Mp)/x[2]**5*(57/8+AQ*57/8))*31536000
      
    return Eqs 

#t_xs=run_ode.run_Q((5,2,4,0),1e-2,1e8,Q_model_exactly)
    

def model(t,x):
    Eqs= np.zeros((2)) 
    
    Eqs[0] = 1 * x[0]
    Eqs[1] = 1 * x[0]
    return Eqs


t = np.linspace(0,20)

y0 = 5

#y = scipy.integrate.odeint(Q_model_exactly,y0,t)
y = scipy.integrate.odeint(model,y0,t)

'''