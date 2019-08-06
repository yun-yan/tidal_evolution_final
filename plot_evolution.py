#plot

import numpy as np
import matplotlib.pyplot as plt
import os 

def plot(t_xs,x0,P,S,canvas_number):
     
    plt.clf
    
    plt.subplot(round((canvas_number+1)/2),2,3) 
    plt.plot(t_xs[:,0],t_xs[:,1]) 
    plt.xscale('log')
    plt.grid(True)
    # label axes
    plt.ylabel('spin angular velocities of {0} (mean motion)'.format(P),fontsize=30)
    plt.xlabel('time (year)',fontsize=30)
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)    
    
    
    plt.subplot(round((canvas_number+1)/2),2,4)
    plt.plot(t_xs[:,0],t_xs[:,2])
    plt.xscale('log')
    plt.grid(True)
    plt.ylabel('spin angular velocities of {0} (mean motion)'.format(S),fontsize=30)
    plt.xlabel('time  (year)',fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    
    
    plt.subplot(round((canvas_number+1)/2),2,1)
    plt.plot(t_xs[:,0],t_xs[:,3])
    plt.xscale('log')
    # Show the major grid lines with dark grey lines
    plt.grid(b=True, which='major', color='#666666', linestyle='-')

    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
    plt.ylabel('orbital semimajor axis ({0} radius)'.format(P),fontsize=30)
    plt.xlabel('time  (year)',fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)    
    
    plt.subplot(round((canvas_number+1)/2),2,2)
    plt.plot(t_xs[:,0],t_xs[:,4],label='{0}'.format(x0[3]))
    #plt.plot(t_xs[:,0],t_xs[:,4])
    plt.xscale('log')
    plt.grid(True)
    plt.ylabel('orbital eccentricity',fontsize=30)
    plt.xlabel('time  (year)',fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    
    if canvas_number==5:
        plt.subplot(round((canvas_number+1)/2),2,5)
        plt.plot(t_xs[:,0],t_xs[:,5])
        plt.xscale('log')
        plt.grid(True)
        plt.ylabel('inclination(degree)',fontsize=30)
        plt.xlabel('time  (year)',fontsize=30)
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)
    else:
        pass
    try:
        file_dir="fig_outputt_at_paper"
        os.mkdir(file_dir)
    except:
        pass
    plt.savefig(os.getcwd()+"\\"+file_dir+"\\{}_{}_{}_{}.png".format(x0[0],x0[1],x0[2],x0[3]) )
    
    return 
        
    
    
