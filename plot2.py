import numpy as np
import matplotlib.pyplot as plt

try:
    data=np.loadtxt('fort.210')
    xmax=int(data[:,0].max())+1
    ymax=int(data[:,1].max())+1
    x0=data[:,0].reshape(xmax,ymax)
    y0=data[:,1].reshape(xmax,ymax)
    chi0=data[:,3].reshape(xmax,ymax)

    x1=np.hstack((x0,x0+xmax))
    x=np.vstack((x1,x1))[:xmax+1,:ymax+1]
    y1=np.hstack((y0,y0))
    y=np.vstack((y1,y1+ymax))[:xmax+1,:ymax+1]
    chi1=np.hstack((chi0,chi0))
    chi=np.vstack((chi1,chi1))[:xmax+1,:ymax+1]

    xlim_list=[0,xmax]
    ylim_list=[0,ymax]
    axislabel=['0','0.5','1']
    fig=plt.figure()
    ax=fig.add_subplot(111,xlim=xlim_list,ylim=ylim_list,title='$\chi_s$ map')
    ax.set_xticks([0,int(xmax*.5),xmax])
    ax.set_xticklabels(axislabel)
    ax.set_yticks([0,int(ymax*.5),ymax])
    ax.set_yticklabels(axislabel)
    cont=ax.contourf(x,y,chi,100,cmap=plt.jet())
    fig.colorbar(cont)
    plt.show()
except OSError:
    pass
