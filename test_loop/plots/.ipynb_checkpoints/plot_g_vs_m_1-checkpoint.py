#!/usr/bin/env python
print ("Loading Libraries")
import matplotlib
matplotlib.use('Agg')
import re
import os, sys
#clean the storage folder before reading .out files
os.system("rm -rf storage/result_g_m_1/.ipynb_checkpoints")
import pickle, cPickle
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

#This script needs to be executed in python2

def plot_logx_logy_logz(xvalues, xlabel, xmin, xmax, yvalues, ylabel, ymin, ymax, zvalues, zlabel, title, filename, linecolor, contourstyle, linelabel, contourlevel, contourlinewidth, clabel_positions, clabelformat):
    print ("Plotting '"+title+"'")

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
    
    #scenario_024
    plt.axis([xmin, xmax, ymin, ymax])
    #scenario_13
    #plt.axis([1E-9, 1e-1, 1E-10, 1e-6])
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


    cons2 = plt.contour(X, Y, Zm, [contourlevel] ,colors=linecolor,locator=mpl.ticker.LogLocator(), linewidths=contourlinewidth, linestyles=contourstyle)
#
    labels = [linelabel]
    for i in range(len(labels)):
	cons2.collections[i].set_label(labels[i])
    
    plt.legend(loc='lower left',prop={'size': 9}, ncol = 1)

#    plt.scatter(map(lambda k: k, x), map(lambda k: 10**k, y), s=-20, c=z, vmin=zmin, vmax=zmax, cmap = mpl.cm.Greys_r,norm=mpl.colors.LogNorm())


##############
#tau decay width limit 
sigmaBrtaupinu=0.057e-2
##############

def Read(fname, key = "(ind. and cum.):"):
    f = open(fname)
    data = f.read()
    sd = data.split(key)[1:]
    Ind, Cum = [], []
    for s in sd:
        s = s[:s.find("\n")]
        s = s.split(' ')
        Ind.append(float(s[1]))
        Cum.append(float(s[2]))
    return Ind, Cum

def ReadLine(resdir, filename, key):
#     key = key + ":"
    with open(resdir+"/"+filename,'r') as f: 
        line = next((l for l in f if key in l), None)
        return line[line.find(key)+len(key):]

########################################################################################################
########################################################################################################
########################################################################################################
print ("Reading in Data")

fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1) 
massstring = 'test'
keys = ["mX:", "g_EE_L:", "g_TE_L:","BRtau2xmu:", "BRtau2xe:", "BRtau2xmu + BRtau2xe:", "produced tau excluding those daughters are also tau lepton:", "BRx2tautau:", "BRx2mumu:", "BRx2ee:", "BRx2gmgm:", "Total Gamma of X [GeV]: ", "ctau [m]:", "produced X:", "n_prompt_smallr: ", "n_prompt_d0z0: ", "n_prompt_baselineeff: ", "prompt small_r efficiencies (ind. and cum.): ", "prompt d0z0 efficiencies (ind. and cum.): ", "prompt baseline efficiencies (ind. and cum.): ", "n_disp_fidvol: ", "n_disp_baselineeff: ", "n_disp_DispTrackEff: ", "displaced fiducial volume efficiencies (ind. and cum.): ", "displaced baseline efficiencies (ind. and cum.): ", "displaced displaced-tracking efficiencies (ind. and cum.): ", "n_mod_DispTrackEff: ", "modified efficiencies: ", "reallyProducedX: ", "Prompt_reallyobservedX: ", "Prompt_reallyvisibleX: ", "Displaced_reallyobservedX: ", "Displaced_reallyvisibleX: ", "Modified_reallyobservedX: ", "Modified_reallyvisibleX: "]
for i in range(len(sys.argv))[1:]:
    resdircollection = sys.argv[i]
    print ("")
    for resdir in resdircollection.split(":"):
    #os.chdir(resdir)
        data=np.zeros((1,41))
        print (" - "+resdir)
        #counter = 0
        for filename in os.listdir(resdir):
            newdata=np.zeros((1,41))
            
            ndata = 0
            KEYS = []
            for k in keys:
                if "(ind. and cum.)" in k:
                    test = ReadLine(resdir, filename, k)
                    print(filename)
                    newdata[0,ndata] = test.split()[0]
                    KEYS.append(k[:k.find("(ind. and cum.)")] + "(ind.)")
                    newdata[0,ndata+1] = test.split()[1]
                    KEYS.append(k[:k.find("(ind. and cum.)")] + "(cum.)")
                    ndata +=2
                else:
                    newdata[0,ndata] = ReadLine(resdir, filename, k)
                    #print(filename)
                    KEYS.append(k)
                    ndata +=1
            #print(filename)
            #counter = counter + 1
            #print(counter)
            #print (resdir+"/"+filename)
           
            #data = data.astype(np.double)
            data=np.vstack((data,newdata))
    #os.chdir(resdir)
        data=np.delete(data,0,0)
    if i==1:
        data1=data
   
        
#converting nan to zero
where_are_NaNs = isnan(data1)
data1[where_are_NaNs] = 0

#export data
#np.savetxt("plots/data/data_311.dat", data, fmt="%s")



# os.chdir("..")	   
print ("Plotting")

font = {'size'   : 18}
mpl.rc('font', family='serif')
plt.rc("text", usetex=True)
fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1) 
plt.rc("text", usetex=True)
         
#title = r'$\epsilon_{\textrm{det.}}=8.4\%$' 

#scenario_0
title = r'$g \,$ vs. $m_{X}$ $(g = g_{\tau e} = g_{ee})$'
xlabel = r'$m_{X}$ [GeV]'
ylabel = r'$g$ [GeV$^{-1}$]'
zlabel = r'$three_event$'
filename = "g_vs_m_1"

#plot signal
#[1E-12, 1e-4, 1E-12, 1e-6]

plot_logx_logy_logz(data1[:,KEYS.index("mX:")], xlabel, -0.05, 1.85, data1[:,KEYS.index("g_TE_L:")], ylabel, 1e-9, 1e1, data1[:,KEYS.index("Prompt_reallyvisibleX: ")]*0.988/50, zlabel, title, filename, 'c','-',r'Prompt Belle', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data1[:,KEYS.index("mX:")], xlabel, -0.05, 1.85, data1[:,KEYS.index("g_TE_L:")], ylabel, 1e-9, 1e1, data1[:,KEYS.index("Displaced_reallyvisibleX: ")]*0.988/50, zlabel, title, filename, 'c','--',r'Displaced Belle', 3 , 1 ,[(0.8,5e-8)],'3signal')

plot_logx_logy_logz(data1[:,KEYS.index("mX:")], xlabel, -0.05, 1.85, data1[:,KEYS.index("g_TE_L:")], ylabel, 1e-9, 1e1, data1[:,KEYS.index("Prompt_reallyvisibleX: ")], zlabel, title, filename, 'k','-',r'Prompt Belle II', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data1[:,KEYS.index("mX:")], xlabel, -0.05, 1.85, data1[:,KEYS.index("g_TE_L:")], ylabel, 1e-9, 1e1, data1[:,KEYS.index("Displaced_reallyvisibleX: ")], zlabel, title, filename, 'k','--',r'Displaced Belle II', 3 , 1 ,[(0.8,5e-8)],'3signal')


#use data files expoted by Mathematica where MovingAverage was used:
#signal311NoBG = genfromtxt("plots/data/signal_311_NoBG.dat");
#signal311NoBG1InvAb = genfromtxt("plots/data/signal_311_NoBG_1InvAb.dat");

#plt.plot(signal311NoBG[:,0], signal311NoBG[:,1],  color='blue',linewidth=1,linestyle='solid' ,label=r"$N_S = ~ \,$3, 50 ab$^{-1}$ (Belle II)" );
#plt.plot(signal311NoBG1InvAb[:,0], signal311NoBG1InvAb[:,1],  color='darkturquoise',linewidth=1,linestyle='solid' ,label=r"$N_S = ~ \,$3, 1 ab$^{-1}$ (Belle)" );

#plt.legend(loc='upper right',prop={'size': 9})
plt.legend(loc='upper center', shadow=True, fontsize='small')

#plot ctau
#for exponent in range(-7,5,2):
#	plot_logx_logy_logz(data[:,0], xlabel, data[:,1], ylabel, data[:,4], zlabel, title, filename, 'darkorange', 'dashed', 10**exponent,0.5,[(1.4,2e-4)],'log')
#for exponent in range(5,6):
#	plot_logx_logy_logz(data[:,0], xlabel, data[:,1], ylabel, data[:,4], zlabel, title, filename, 'darkorange', 'dashed', 10**exponent,0.5,[(0.5,2e-9)],'log')


#plot tau decay BR constraint   removed because it is always for ctau < 1 m
#plot_logx_logy_logz(data[:,0], xlabel, data[:,1], ylabel, data[:,8], zlabel, title, filename, 'red', 'dashdot', 2*sigmaBrtaupinu , 0.5, [(1.35,8e-6)],'tauBR')   


#plot LQD coupling constraint
#msfermion=[1e3,5e3]
#for mass in msfermion:
#	limit = 0.20/(1e3*mass)+0.046/(mass**2)
#	plt.plot([0,2],[limit,limit],linestyle='dashed',color='red',linewidth=0.5)
#	plt.text(0.59, limit*0.92, r"$m_{\tilde{d}_{R}}=" + str(int(mass/1000))  + "$ TeV", color='red',ha='left',va='top', fontsize = '8')


#ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
#ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.05))
#ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

plt.grid(which='major',axis='x',alpha=0.5,linestyle='--',linewidth=0.5)
plt.grid(which='major',axis='y',alpha=0.5,linestyle='--',linewidth=0.5)

plt.savefig("plots/"+filename+".png",bbox_inches='tight')
plt.show()

plt.clf()#clear the current figure

##################################################################################################################################################################

print ("Plotting")

font = {'size'   : 18}
mpl.rc('font', family='serif')
plt.rc("text", usetex=True)
fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1) 
plt.rc("text", usetex=True)
         
#title = r'$\epsilon_{\textrm{det.}}=8.4\%$' 

#scenario_0
title = r'$g \,$ vs. $m_{X}$ $(g = g_{\tau e} = g_{ee})$'
xlabel = r'$m_{X}$ [GeV]'
ylabel = r'$g$ [GeV$^{-1}$]'
zlabel = r'$three_event$'
filename = "g_vs_m_mod_1"

#plot signal
#[1E-12, 1e-4, 1E-12, 1e-6]

plot_logx_logy_logz(data1[:,KEYS.index("mX:")], xlabel, -0.05, 1.85, data1[:,KEYS.index("g_TE_L:")], ylabel, 1e-9, 1e1, data1[:,KEYS.index("Prompt_reallyvisibleX: ")]*0.988/50, zlabel, title, filename, 'c','-',r'Prompt Belle', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data1[:,KEYS.index("mX:")], xlabel, -0.05, 1.85, data1[:,KEYS.index("g_TE_L:")], ylabel, 1e-9, 1e1, data1[:,KEYS.index("Modified_reallyvisibleX: ")]*0.988/50, zlabel, title, filename, 'c','--',r'Displaced Belle', 3 , 1 ,[(0.8,5e-8)],'3signal')

plot_logx_logy_logz(data1[:,KEYS.index("mX:")], xlabel, -0.05, 1.85, data1[:,KEYS.index("g_TE_L:")], ylabel, 1e-9, 1e1, data1[:,KEYS.index("Prompt_reallyvisibleX: ")], zlabel, title, filename, 'k','-',r'Prompt Belle II', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data1[:,KEYS.index("mX:")], xlabel, -0.05, 1.85, data1[:,KEYS.index("g_TE_L:")], ylabel, 1e-9, 1e1, data1[:,KEYS.index("Modified_reallyvisibleX: ")], zlabel, title, filename, 'k','--',r'Displaced Belle II', 3 , 1 ,[(0.8,5e-8)],'3signal')
                    
plt.legend(loc='upper center', shadow=True, fontsize='small')

plt.grid(which='major',axis='x',alpha=0.5,linestyle='--',linewidth=0.5)
plt.grid(which='major',axis='y',alpha=0.5,linestyle='--',linewidth=0.5)

plt.savefig("plots/"+filename+".png",bbox_inches='tight')
plt.show()

plt.clf()#clear the current figure
