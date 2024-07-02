import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
try:
    data=np.loadtxt("fort.60")
    data2=np.loadtxt("fort.50")
    sw_spa=True
except OSError:
    sw_spa=False
try:
    data=np.loadtxt('fort.210')
    sw_map=True
except OSError:
    sw_map=False

cdict = {'blue': ((0.0, 0.0, 0.0),
                  (0.25, 0.8, 0.8),
                  (0.377, 1, 1),
                  (0.67, 0, 0),
                  (1, 0,0)),
         'green': ((0.0, 0, 0),
                   (0.174, 0, 0),
                   (0.410, 1, 1),
                   (0.66, 1, 1),
                   (0.915, 0, 0),
                   (1, 0, 0)),
         'red': ((0.0, 0, 0),
                 (0.387, 0, 0),
                 (0.680, 1, 1),
                 (0.896, 1, 1),
                 (1.0, 0.65, 0.65))}
cmapcolor = colors.LinearSegmentedColormap('color3',cdict, 10024)
#cmapcolor = plt.jet()

if sw_spa:
    d1=int(data[:,0].max())+1
    d2=int(len(data[:,0])/d1)
    x=data[:,0].reshape(d1,d2)
    y=data[:,1].reshape(d1,d2)*1.e3
    z=data[:,2].reshape(d1,d2)

    x2=data2[:,0].reshape(d1,d2)
    y2=data2[:,1].reshape(d1,d2)*1.e3
    z2=data2[:,2].reshape(d1,d2)

    ymax=y.max()
    #ymax=150
    #pnum=121
    pnum=111

    fig=plt.figure()
    ax=fig.add_subplot(pnum,xlim=[0,x.max()],ylim=[0,ymax],title='$\chi_s$ map')
    ax.set_xticks([0,x.max()/2,x.max()])
    #ax.set_xticklabels(['0','$2\pi$/c'])
    ax.set_xticklabels(['(0.5,0.5,0)','(0.5,0.5,1)','(0.5,0.5,2)'])
    ax.set_ylabel('$\omega$ meV')
    ax.set_xlabel('$(0.5,0.5,H/(2\pi/c))$')
    cont=ax.contourf(x,y,z,100,cmap=cmapcolor)
    fig.colorbar(cont)

    if False:
        ax1=fig.add_subplot(222,xlim=[0,y.max()],ylim=[0,z.max()*1.1],title='$\chi_s(Q=(\pi,0,k_z))$')
        ax1.set_ylabel('Intensity')
        ax1.set_xlabel('$\omega$ meV')
        for x1,y1,z1 in zip(x,y,z):
            ax1.plot(y1,z1,c=cm.jet(x1[0]/x.max()))
        ax2=fig.add_subplot(224,xlim=[0,y2.max()],ylim=[0,z2.max()*1.1],title='$\chi_0(Q=(\pi,0,k_z))$')
        ax2.set_ylabel('Intensity')
        ax2.set_xlabel('$\omega$ meV')
        for x3,y3,z3 in zip(x2,y2,z2):
            ax2.plot(y3,z3,c=cm.jet(x3[0]/x2.max()))
    plt.show()

if sw_map:
    xmax=int(data[:,0].max())+1
    ymax=int(data[:,1].max())+1
    x0=data[:,0].reshape(xmax,ymax)
    y0=data[:,1].reshape(xmax,ymax)
    chi0=data[:,3].reshape(xmax,ymax)

    x1=np.hstack((x0-xmax,x0))
    x=np.vstack((x1,x1)) #+xmax
    y1=np.hstack((y0,y0))
    y=np.vstack((y1-ymax,y1)) #+ymax
    chi1=np.hstack((chi0,chi0))
    chi=np.vstack((chi1,chi1))

    xlim_list=[0,xmax]
    ylim_list=[0,ymax]
    #xlim_list=[-xmax,xmax]
    #ylim_list=[-ymax,ymax]
    #axislabel=['-1','-0.5','0','0.5','1']
    axislabel=['0','0.5','1']
    fig=plt.figure()
    ax=fig.add_subplot(111,xlim=[0,xmax],ylim=[0,ymax],title='$\chi_s$ map')
    #ax.set_xticks([-xmax,-int(xmax*.5),0,int(xmax*.5),xmax])
    ax.set_xticks([0,int(xmax*.5),xmax])
    ax.set_xticklabels(axislabel)
    #ax.set_yticks([-ymax,-int(ymax*.5),0,int(ymax*.5),ymax])
    ax.set_yticks([0,int(ymax*.5),ymax])
    ax.set_yticklabels(axislabel)
    cont=ax.contourf(x,y,chi,100,cmap=cmapcolor)
    fig.colorbar(cont)
    plt.show()
