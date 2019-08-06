import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from collections import OrderedDict


y=np.linspace(1,18,100)
x=np.linspace(1,18,100)
X,Y=np.meshgrid(x,y)
def consevation_angular(xa,ap):
    return (-ap-0.3181499860691723 *xa**2+1.3265927734795924*xa**(1.5))/0.03925018904571284


Z=[]
for i in x:
    z=[]
    for j in y:
        ze=consevation_angular(j,i)
        if ze<1:
            z.append(-1)
        else:
            z.append(ze)
    Z.append(z)
            

x_minic=np.linspace(1,17,100)
y_minic = -0.3181499860691723 *x_minic**2+1.3265927734795924*x_minic**1.5-0.03925018904571284*10
    
plt.axes(xlim=(1, 18), ylim=(1, 18))
plt.pcolormesh(X,Y,Z)
plt.plot(x_minic, y_minic,'w',label=r'$\alpha_c$=10,e=0')
legend=plt.legend(loc=1,fancybox=True, framealpha=0.1)
plt.setp(legend.get_texts(), color='w')

plt.colorbar()
plt.savefig('fit_show.png')

