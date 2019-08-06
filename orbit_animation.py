import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import patches
import main

#tx=main.data((5.3,2,3.9,0.2))
tx=np.load('test_1.npy')

G=6.674e-11
Mp=870.3/G*1e9
Ms=101.4/G*1e9


tlist=tx[:,0]
xdata=[]
ydata=[]
pdata=[0]
cdata=[0]
ndata=[]


for i in range(len(tlist)-1):
    n=(G*(Mp+Ms)/(tx[:,3][i]*1153e3)**3)**0.5*31536000/10000
    xdata.append(tx[:,3][i]*np.cos(n*tlist[i]))
    ydata.append(tx[:,3][i]*np.sin(n*tlist[i]))
    pdata.append(pdata[i]+tx[:,1][i]*n*(tlist[i+1]-tlist[i])*180/np.pi)
    cdata.append(cdata[i]+tx[:,2][i]*n*(tlist[i+1]-tlist[i])*180/np.pi)
    ndata.append(n*tlist[i]*180/np.pi)






fig=plt.figure()
plt.rcParams['figure.figsize'] = (8.0, 8.0)
ax = plt.axes(xlim=(-30, 30), ylim=(-30, 30))

pluto = plt.Circle((0, 0), 1, color='r')
pluto_ellipse =patches.Ellipse((0, 0), 3, 2,
                     angle=0, linewidth=1, fill=False, zorder=2)
pluto_face =patches.Ellipse((0, 0), 3, 2,
                     angle=0, linewidth=1, color='r',fill=False, zorder=2)
charon = plt.Circle((4, 0), 0.5, color='b')
charon_ellipse =patches.Ellipse((4, 0), 1.5, 1,
                     angle=0, linewidth=1, fill=False, zorder=2)
charon_face =patches.Ellipse((4, 0), 1.5, 1,
                     angle=0, linewidth=1, color='b',fill=False, zorder=2)



def init():
#    pluto.center=(data[0,0],1)
#    pluto_ellipse.angle=(data[0,0])
#    charon.center=(data[0,0],4)
#    charon_ellipse.angle=(data[0,0])

    ax.add_patch(pluto)
    ax.add_patch(pluto_ellipse)
    ax.add_patch(pluto_face)
    ax.add_patch(charon)
    ax.add_patch(charon_ellipse)
    ax.add_patch(charon_face)

    return pluto, pluto_ellipse,pluto_face, charon, charon_ellipse, charon_face

def animate(i):
    # for state in data:
    j=i+49999000
    #j=i
    pluto_ellipse.angle=(ndata[j])
    pluto_face.angle=(pdata[j])
    charon.center=(xdata[j],ydata[j])
    charon_ellipse.angle=(ndata[j])
    charon_ellipse.center=(xdata[j],ydata[j])
    charon_face.angle=(cdata[j])
    charon_face.center=(xdata[j],ydata[j])
    return pluto, pluto_ellipse,pluto_face, charon, charon_ellipse, charon_face


anim=animation.FuncAnimation(fig,animate,init_func=init,frames=1000,blit=True)

#pyplot.show()
anim.save('final_evolution.mp4', fps=60)

