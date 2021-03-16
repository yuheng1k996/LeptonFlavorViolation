#!/usr/bin/env python
print "Loading Libraries"
import re
import os, sys
import pickle, cPickle, sys
import pprint
from math import fabs, log, exp
from numpy import log10
from numpy import *

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
#from matplotlib.mlab import griddata
import numpy as np
from scipy.interpolate import griddata
from math import ceil, sqrt
from copy import deepcopy

from pylab import genfromtxt 
from matplotlib.pyplot import text

def plot_logx_logy_logz(xvalues, xlabel, yvalues, ylabel, zvalues, zlabel, title, filename, linecolor, contourstyle, contourlevel, contourlinewidth,clabel_positions,clabelformat):
    print "Plotting '"+title+"'"

    # remove invalid entries which have either x or z <= 0
    newxvalues, newyvalues, newzvalues = xvalues[:], yvalues[:], zvalues[:]
    npopped = 0
    newX, newY, newZ =  list(), list(), list()
    for i in range(len(xvalues)):
           newX.append(xvalues[i])
           newY.append(yvalues[i])
           if zvalues[i] == 0:
              newZ.append(1E-99)
           else:
              newZ.append(zvalues[i])
	   
    x = np.array(newX)
    y = np.array(map(log10, newY))
    z = np.array(newZ)
    if len(x) == 0:
       return
    plt.xlabel(xlabel,fontsize = 16)
    plt.ylabel(ylabel,fontsize = 16)
    plt.title(title,fontsize = 16)
    plt.axis([0.54, 1.645, 1E-8, 5E-4])
    zmin, zmax =  1E-3, 1E3 #min([z[i] for i in range(len(z)) if x[i] >= -9 and x[i] <= -5 and y[i] >= -9 and y[i] <= -5 and z[i] != 0.])

    ax.set_xscale('linear')
    ax.set_yscale('log')

    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)

    
    n=200
    xg = np.linspace(x.min(),x.max(),n)
    yg = np.linspace(y.min(),y.max(),n)
    X,Y = np.meshgrid(xg,yg)

    # interpolate Z values on defined grid
    Z = griddata(np.vstack((x.flatten(),y.flatten())).T, \
    np.vstack(z.flatten()),(X,Y),method='linear').reshape(X.shape)
    # mask nan values, so they will not appear on plot
    Zm = np.ma.masked_where(np.isnan(Z),Z)

    X = np.array(map(lambda k : k, X))    
    Y = np.array(map(lambda k : 10.**k, Y))
    
    colmap = mpl.cm.Greys_r
    colmap.set_bad('k',1)
    #ax.pcolormesh(X,Y,Zm,vmin=zmin,vmax=zmax,shading='gouraud', norm=mpl.colors.LogNorm(),cmap = colmap)



    cons2 = plt.contour(X, Y, Zm, [contourlevel] ,colors=linecolor,locator=mpl.ticker.LogLocator(), linewidths=contourlinewidth, linestyles = contourstyle)
#
# Label levels with specially formatted floats
    if (clabelformat == "log"):
    	#fmt = mpl.ticker.LogFormatterMathtext()
    	#fmt.create_dummy_axis()
    	fmt = {}
	for l in cons2.levels:
   		fmt[l] = r'$10^{ ' + str("%.f" % log10(l)) + r'}$ m'
    	ax.clabel(cons2, cons2.levels, inline=True, fmt=fmt, fontsize=8,manual=clabel_positions)
    elif (clabelformat == '3signal'):
    	fmt = {}
	for l in cons2.levels:
		#fmt[l] = str("%.f" % l)+r' signal events'
		fmt[l] = r'$N_S=3$'
    	ax.clabel(cons2, cons2.levels, inline=True, fmt=fmt, fontsize=8,manual=clabel_positions)
    elif (clabelformat == '21signal'):
    	fmt = {}
	for l in cons2.levels:
		#fmt[l] = str("%.f" % l)+r' signal events'
		fmt[l] = r'$N_S=21$'
    	ax.clabel(cons2, cons2.levels, inline=True, fmt=fmt, fontsize=8,manual=clabel_positions)
    elif (clabelformat == 'tauBR'):
	#fmt = r'$\delta{\cal B}(\tau^-\rightarrow \pi^- \nu_\tau)=0.0005$'
	fmt = r'${\mathcal B(\tau\to \pi \tilde{\chi}_1^0 )}=2\sigma_{\mathcal B(\tau\to \pi \nu_\tau)}$'
    	ax.clabel(cons2, cons2.levels, inline=True, fmt=fmt, fontsize=8,manual=clabel_positions)


    plt.scatter(map(lambda k: 10**k, x), map(lambda k: 10**k, y), s=-20, c=z, vmin=zmin, vmax=zmax, cmap = mpl.cm.Greys_r,norm=mpl.colors.LogNorm())


##############
#tau decay width limit 
sigmaBrtaupinu=0.057e-2
##############

########################################################################################################
########################################################################################################
########################################################################################################
print "Reading in Data"

fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1) 
massstring = 'test'
for i in range(len(sys.argv))[1:]:
  resdircollection = sys.argv[i]
  print ""
  for resdir in resdircollection.split(":"):
    #os.chdir(resdir)
        data=np.zeros((1,10))
        print " - "+resdir
        #counter = 0
        for filename in os.listdir(resdir):
            newdata=np.zeros((1,10))
            #print(filename)
            #counter = counter + 1
            #print(counter)
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "mN1:" in l), None)
                 newdata[0,0]=line[line.find("mN1:")+len("mN1:"):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "lambOverMSq:" in l), None)
                 newdata[0,1]=line[line.find("lambOverMSq:")+len("lambOverMSq:"):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "reallyobservedN1:" in l), None)
                 newdata[0,2]=line[line.find("reallyobservedN1:")+len("reallyobservedN1:"):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "reallyvisibleN1:" in l), None)
                 newdata[0,3]=line[line.find("reallyvisibleN1:")+len("reallyvisibleN1:"):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "ctau [m]:" in l), None)
                 newdata[0,4]=line[line.find("ctau [m]:")+len("ctau [m]:"):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "BRtauTon1pim:" in l), None)
                 newdata[0,5]=line[line.find("BRtauTon1pim:")+len("BRtauTon1pim:"):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "BRtauTon1rhom:" in l), None)
                 newdata[0,6]=line[line.find("BRtauTon1rhom:")+len("BRtauTon1rhom:"):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "BRtauTon1a1m:" in l), None)
                 newdata[0,7]=line[line.find("BRtauTon1a1m:")+len("BRtauTon1a1m:"):]
                 newdata[0,8]=newdata[0,5]+newdata[0,6]+newdata[0,7]
                 newdata[0,9]=3./newdata[0,2]*newdata[0,8]
            data=np.vstack((data,newdata))
    #os.chdir(resdir)
        data=np.delete(data,0,0)

#converting nan to zero
where_are_NaNs = isnan(data)
data[where_are_NaNs] = 0


#export data
np.savetxt("plots/data/data_311.dat", data, fmt="%s")



# os.chdir("..")	   
print "Plotting"

font = {'size'   : 18}
mpl.rc('font', family='serif')
plt.rc("text", usetex=True)
fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1) 
plt.rc("text", usetex=True)
         
#title = r'$\epsilon_{\textrm{det.}}=8.4\%$' 
title = r'$\lambda^\prime_{311} \neq 0$'
xlabel = r'$m_{\tilde{\chi}^0_1}~[\mathrm{GeV}]$'
ylabel = r'$\lambda^\prime_{311}/m^2_{\tilde{f}}~[\mathrm{GeV}^{-2}]$'
zlabel = r'$three_event$'
filename = "lambOverMSq_vs_mass_311"

#plot signal
#plot_logx_logy_logz(data[:,0], xlabel, data[:,1], ylabel, data[:,3], zlabel, title, filename, 'k', 'solid', 3 , 1 ,[(0.8,5e-8)],'3signal')


#use data files expoted by Mathematica where MovingAverage was used:
signal311NoBG = genfromtxt("plots/data/signal_311_NoBG.dat");
signal311NoBG1InvAb = genfromtxt("plots/data/signal_311_NoBG_1InvAb.dat");

plt.plot(signal311NoBG[:,0], signal311NoBG[:,1],  color='blue',linewidth=1,linestyle='solid' ,label=r"$N_S = ~ \,$3, 50 ab$^{-1}$ (Belle II)" );
plt.plot(signal311NoBG1InvAb[:,0], signal311NoBG1InvAb[:,1],  color='darkturquoise',linewidth=1,linestyle='solid' ,label=r"$N_S = ~ \,$3, 1 ab$^{-1}$ (Belle)" );

plt.legend(loc='upper right',prop={'size': 9})

#plot ctau
for exponent in range(-7,5,2):
	plot_logx_logy_logz(data[:,0], xlabel, data[:,1], ylabel, data[:,4], zlabel, title, filename, 'darkorange', 'dashed', 10**exponent,0.5,[(1.4,2e-4)],'log')
for exponent in range(5,6):
	plot_logx_logy_logz(data[:,0], xlabel, data[:,1], ylabel, data[:,4], zlabel, title, filename, 'darkorange', 'dashed', 10**exponent,0.5,[(0.5,2e-9)],'log')


#plot tau decay BR constraint   removed because it is always for ctau < 1 m
#plot_logx_logy_logz(data[:,0], xlabel, data[:,1], ylabel, data[:,8], zlabel, title, filename, 'red', 'dashdot', 2*sigmaBrtaupinu , 0.5, [(1.35,8e-6)],'tauBR')   


#plot LQD coupling constraint
msfermion=[1e3,5e3]
for mass in msfermion:
	limit = 0.20/(1e3*mass)+0.046/(mass**2)
	plt.plot([0,2],[limit,limit],linestyle='dashed',color='red',linewidth=0.5)
	plt.text(0.59, limit*0.92, r"$m_{\tilde{d}_{R}}=" + str(int(mass/1000))  + "$ TeV", color='red',ha='left',va='top', fontsize = '8')


ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.05))
ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

plt.grid(which='major',axis='x',alpha=0.5,linestyle='--',linewidth=0.5)
plt.grid(which='major',axis='y',alpha=0.5,linestyle='--',linewidth=0.5)
plt.savefig("plots/"+filename+".png",bbox_inches='tight')


#br vs ctau plot
plt.clf()#clear the current figure
#plt.cla()#clear the current axes

fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1)

# os.chdir("..")	   
print "Plotting"

font = {'size'   : 18}
mpl.rc('font', family='serif')
plt.rc("text", usetex=True)
fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1) 
plt.rc("text", usetex=True)
title = r'$\lambda^\prime_{311} \neq 0$, 50 ab$^{-1}$ (Belle II)' 
xlabel = r'$c\tau_{\tilde{\chi}_1^0}~[\mathrm{m}]$'
#ylabel = r'${\cal B}(\tau\rightarrow \tilde{\chi}_1^0+X)\cdot{\cal B}(\tilde{\chi}_1^0\rightarrow$vis.$)\cdot \epsilon_{\textrm{det.}}$'
ylabel = r'${\cal B}(\tau\rightarrow \tilde{\chi}_1^0+X)$'
zlabel = r'$three_event$'
filename = "br_vs_ctau_311"


data0p6GeV = np.zeros((1,10))
for x in data:
	if (x[0]==0.6  and x[9]<1):
		data0p6GeV=np.vstack((data0p6GeV,x))
data0p6GeV=np.delete(data0p6GeV,0,0)
data0p6GeV = sorted(data0p6GeV, key = lambda x: x[4])
#print (np.transpose(data0p6GeV)[4],np.transpose(data0p6GeV)[9] )

data1GeV = np.zeros((1,10))
for x in data:
	if (x[0]==1  and x[9]<1):
		data1GeV=np.vstack((data1GeV,x))
data1GeV=np.delete(data1GeV,0,0)
data1GeV = sorted(data1GeV, key = lambda x: x[4])
#print (np.transpose(data1GeV)[4],np.transpose(data1GeV)[9] )

data1p5GeV = np.zeros((1,10))
for x in data:
	if (x[0]==1.5 and x[9]<1):
		data1p5GeV=np.vstack((data1p5GeV,x))
data1p5GeV=np.delete(data1p5GeV,0,0)
data1p5GeV = sorted(data1p5GeV, key = lambda x: x[4])
#print (np.transpose(data1p5GeV)[4],np.transpose(data1p5GeV)[9] )

np.savetxt("plots/data/data_311_1p5GeV_tobechanged.dat", data1p5GeV, fmt="%s")

plt.xlabel(xlabel,fontsize = 16)
plt.ylabel(ylabel,fontsize = 16)
plt.title(title,fontsize = 16)
plt.axis([5e-4, 2e6, 2E-12, 5E-3])

plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)

plt.loglog(np.transpose(data0p6GeV)[4],np.transpose(data0p6GeV)[9]  , label = r'0.6 GeV, $N_S=3$', linewidth=1  )
plt.loglog(np.transpose(data1GeV)[4],np.transpose(data1GeV)[9]  , label = r'1.0 GeV, $N_S=3$', linewidth=1  )
#plt.loglog(np.transpose(data1p5GeV)[4],np.transpose(data1p5GeV)[9]  , label = r'1.5 GeV, $N_S=3$' , linewidth=1 )
signal311MI1p5GeV = genfromtxt("plots/data/data_311_1p5GeV.dat");
plt.loglog(signal311MI1p5GeV[:,0],signal311MI1p5GeV[:,1] , label = r'1.5 GeV, $N_S=3$' , linewidth=1)


#br(tau->hadron) limit
#plt.loglog([1e0,1e10],[0.0005,0.0005],color='red', linestyle='dashdot', label = r'$\delta{\cal B}(\tau^-\rightarrow \pi^- \nu_\tau)=0.05\%$',linewidth=0.5)
plt.loglog([1e0,1e0,1e10],[1,2*sigmaBrtaupinu,2*sigmaBrtaupinu],color='red', linestyle='dashdot', #label = r'$\delta{\cal B}(\tau^-\rightarrow \pi^- \nu_\tau)=0.05\%$',linewidth=0.5)
label = r'${\mathcal B(\tau\to \pi \tilde{\chi}_1^0 )}=2\sigma_{\mathcal B(\tau\to \pi \nu_\tau)}$',linewidth=0.5)


plt.text(0.1, 5e-5, r"${\cal B}(\tilde{\chi}_1^0\rightarrow$vis.$)\cdot \epsilon_{\textrm{det.}}=1$", color='black',ha='left',va='top', fontsize = '9')


plt.legend(loc='lower right',prop={'size': 9})
#if the range of x-axis is large enough, the intervals become 100 instead of 10 and the minor ticks are gone. Fix below, which was found online:
#label the major tickers properly:
locmaj = mpl.ticker.LogLocator(base=10.0, subs=(0.1,1.0, ), numticks=100)
ax.xaxis.set_major_locator(locmaj)
#label the minor tickers properly:
locmin = mpl.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1 ),numticks=100) 
ax.xaxis.set_minor_locator(locmin)
ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
#which must be added here instead of inside the plotting funciont plot_logx_logy_logz defined above.

#if the range of y-axis is large enough, the intervals become 100 instead of 10 and the minor ticks are gone. Fix below, which was found online:
#label the major tickers properly:
locmaj = mpl.ticker.LogLocator(base=10.0, subs=(0.1,1.0, ), numticks=100)
ax.yaxis.set_major_locator(locmaj)
#label the minor tickers properly:
locmin = mpl.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1 ),numticks=100) 
ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
#which must be added here instead of inside the plotting funciont plot_logx_logy_logz defined above.


plt.grid(which='major',axis='x',alpha=0.5,linestyle='--',linewidth=0.5)
plt.grid(which='major',axis='y',alpha=0.5,linestyle='--',linewidth=0.5)

plt.savefig("plots/"+filename+".png",bbox_inches='tight')
#3* ([8]/[2])   vs [4]
