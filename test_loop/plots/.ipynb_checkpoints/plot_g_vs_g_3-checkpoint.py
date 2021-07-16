#!/usr/bin/env python
print ("Loading Libraries")
import matplotlib
matplotlib.use('Agg')
import re
import os, sys
#clean the storage folder before reading .out files
os.system("rm -rf storage/results_g_g_3/3_max/.ipynb_checkpoints")
os.system("rm -rf storage/results_g_g_3/3_med/.ipynb_checkpoints")
os.system("rm -rf storage/results_g_g_3/3_min/.ipynb_checkpoints")
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
	   
    x = np.array(map(log10, newX))
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

    ax.set_xscale('log')
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

    X = np.array(map(lambda k : 10.**k, X))    
    Y = np.array(map(lambda k : 10.**k, Y))
    
    colmap = mpl.cm.Greys_r
    colmap.set_bad('k',1)
    #ax.pcolormesh(X,Y,Zm,vmin=zmin,vmax=zmax,shading='gouraud', norm=mpl.colors.LogNorm(),cmap = colmap)


    ax.plot(xg,-1*yg,linecolor,label=linelabel,linestyle=contourstyle)
    cons2 = plt.contour(X, Y, Zm, [contourlevel] ,colors=linecolor,locator=mpl.ticker.LogLocator(), linewidths=contourlinewidth, linestyles=contourstyle)
#


    plt.scatter(map(lambda k: 10**k, x), map(lambda k: 10**k, y), s=-20, c=z, vmin=zmin, vmax=zmax, cmap = mpl.cm.Greys_r,norm=mpl.colors.LogNorm())


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
keys = ["mX:", "g_EE_L:", "g_TM_L:","BRtau2xmu:", "BRtau2xe:", "BRtau2xmu + BRtau2xe:", "produced tau excluding those daughters are also tau lepton:", "BRx2tautau:", "BRx2mumu:", "BRx2ee:", "BRx2gmgm:", "Total Gamma of X [GeV]: ", "ctau [m]:", "produced X:", "n_prompt_smallr: ", "n_prompt_d0z0: ", "n_prompt_baselineeff: ", "prompt small_r efficiencies (ind. and cum.): ", "prompt d0z0 efficiencies (ind. and cum.): ", "prompt baseline efficiencies (ind. and cum.): ", "n_disp_fidvol: ", "n_disp_baselineeff: ", "n_disp_DispTrackEff: ", "displaced fiducial volume efficiencies (ind. and cum.): ", "displaced baseline efficiencies (ind. and cum.): ", "displaced displaced-tracking efficiencies (ind. and cum.): ", "n_mod_DispTrackEff: ", "modified efficiencies: ", "reallyProducedX: ", "Prompt_reallyobservedX: ", "Prompt_reallyvisibleX: ", "Displaced_reallyobservedX: ", "Displaced_reallyvisibleX: ", "Modified_reallyobservedX: ", "Modified_reallyvisibleX: "]
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
                    newdata[0,ndata] = test.split()[0]
                    KEYS.append(k[:k.find("(ind. and cum.)")] + "(ind.)")
                    newdata[0,ndata+1] = test.split()[1]
                    KEYS.append(k[:k.find("(ind. and cum.)")] + "(cum.)")
                    ndata +=2
                else:
                    newdata[0,ndata] = ReadLine(resdir, filename, k)
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
    elif i==2:
        data2=data
    elif i==3:
        data3=data
    
    
        
#converting nan to zero
where_are_NaNs = isnan(data1)
data1[where_are_NaNs] = 0
where_are_NaNs = isnan(data2)
data2[where_are_NaNs] = 0
where_are_NaNs = isnan(data3)
data3[where_are_NaNs] = 0


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

#scenario_3
title = r'$g_{\tau \mu}$ vs. $g_{ee} (Belle II)$'
xlabel = r'$g_{ee}$ [GeV$^{-1}$]'
ylabel = r'$g_{\tau \mu}$ [GeV$^{-1}$]'
zlabel = r'$three_event$'
filename = "g_vs_g_3"

#plot signal
#[1E-9, 1e-1, 1E-10, 1e-6]


plot_logx_logy_logz(data3[:,KEYS.index("g_EE_L:")], xlabel, 1e-9, 1e-1, data3[:,KEYS.index("g_TM_L:")], ylabel, 1e-10, 1e-6, data3[:,KEYS.index("Prompt_reallyvisibleX: ")], zlabel, title, filename, 'c','-',r'Prompt $\:m_{X} = 0.5 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data2[:,KEYS.index("g_EE_L:")], xlabel, 1e-9, 1e-1, data2[:,KEYS.index("g_TM_L:")], ylabel, 1e-10, 1e-6, data2[:,KEYS.index("Prompt_reallyvisibleX: ")], zlabel, title, filename, 'r','-',r'Prompt $\: m_{X} = 1.0 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data1[:,KEYS.index("g_EE_L:")], xlabel, 1e-9, 1e-1, data1[:,KEYS.index("g_TM_L:")], ylabel, 1e-10, 1e-6, data1[:,KEYS.index("Prompt_reallyvisibleX: ")], zlabel, title, filename, 'k','-',r'Prompt $\: m_{X} = 1.5 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')

plot_logx_logy_logz(data3[:,KEYS.index("g_EE_L:")], xlabel, 1e-9, 1e-1, data3[:,KEYS.index("g_TM_L:")], ylabel, 1e-10, 1e-6, data3[:,KEYS.index("Displaced_reallyvisibleX: ")], zlabel, title, filename, 'c','--',r'Displaced $\: m_{X} = 0.5 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data2[:,KEYS.index("g_EE_L:")], xlabel, 1e-9, 1e-1, data2[:,KEYS.index("g_TM_L:")], ylabel, 1e-10, 1e-6, data2[:,KEYS.index("Displaced_reallyvisibleX: ")], zlabel, title, filename, 'r','--',r'Displaced $\: m_{X} = 1.0 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data1[:,KEYS.index("g_EE_L:")], xlabel, 1e-9, 1e-1, data1[:,KEYS.index("g_TM_L:")], ylabel, 1e-10, 1e-6, data1[:,KEYS.index("Displaced_reallyvisibleX: ")], zlabel, title, filename, 'k','--',r'Displaced $\: m_{X} = 1.5 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')


#plt.legend(loc='upper right',prop={'size': 9})
plt.legend(loc='lower left', shadow=True, fontsize='small')


#plot ctau
plt.grid(which='major',axis='x',alpha=0.5,linestyle='--',linewidth=0.5)
plt.grid(which='major',axis='y',alpha=0.5,linestyle='--',linewidth=0.5)

plt.savefig("plots/"+filename+".png",bbox_inches='tight')
plt.show()

#br vs ctau plot
plt.clf()#clear the current figure

print ("Plotting")

font = {'size'   : 18}
mpl.rc('font', family='serif')
plt.rc("text", usetex=True)
fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1) 
plt.rc("text", usetex=True)

#scenario_3
title = r'${\rm Br}(\tau \to X \, \mu)$ vs. $c\tau(X \longrightarrow e \, \bar{e}) (Belle II)$'
xlabel = r'$c\tau(X \longrightarrow e \, \bar{e})$ [m]'
ylabel = r'${\rm Br}(\tau \longrightarrow X \, \mu)$'
zlabel = r'$three_event$'
filename = "br_vs_ct_3"

#plt.cla()#clear the current axes
#[1E-6, 1e+2, 1E-10, 1e-3]
plot_logx_logy_logz(data3[:,KEYS.index("ctau [m]:")], xlabel, 1e-6, 1e+2, data3[:,KEYS.index("BRtau2xmu + BRtau2xe:")], ylabel, 1e-10, 1e-3, data3[:,KEYS.index("Prompt_reallyvisibleX: ")], zlabel, title, filename, 'c','-',r'Prompt$ \: m_{X} = 0.5 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data2[:,KEYS.index("ctau [m]:")], xlabel, 1e-6, 1e+2, data2[:,KEYS.index("BRtau2xmu + BRtau2xe:")], ylabel, 1e-10, 1e-3, data2[:,KEYS.index("Prompt_reallyvisibleX: ")], zlabel, title, filename, 'r','-',r'Prompt$ \: m_{X} = 1.0 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data1[:,KEYS.index("ctau [m]:")], xlabel, 1e-6, 1e+2, data1[:,KEYS.index("BRtau2xmu + BRtau2xe:")], ylabel, 1e-10, 1e-3, data1[:,KEYS.index("Prompt_reallyvisibleX: ")], zlabel, title, filename, 'k','-',r'Prompt$ \: m_{X} = 1.5 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')

plot_logx_logy_logz(data3[:,KEYS.index("ctau [m]:")], xlabel, 1e-6, 1e+2, data3[:,KEYS.index("BRtau2xmu + BRtau2xe:")], ylabel, 1e-10, 1e-3, data3[:,KEYS.index("Displaced_reallyvisibleX: ")], zlabel, title, filename, 'c','--',r'Displaced$ \:m_{X} = 0.5 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data2[:,KEYS.index("ctau [m]:")], xlabel, 1e-6, 1e+2, data2[:,KEYS.index("BRtau2xmu + BRtau2xe:")], ylabel, 1e-10, 1e-3, data2[:,KEYS.index("Displaced_reallyvisibleX: ")], zlabel, title, filename, 'r','--',r'Displaced$ \:m_{X} = 1.0 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data1[:,KEYS.index("ctau [m]:")], xlabel, 1e-6, 1e+2, data1[:,KEYS.index("BRtau2xmu + BRtau2xe:")], ylabel, 1e-10, 1e-3, data1[:,KEYS.index("Displaced_reallyvisibleX: ")], zlabel, title, filename, 'k','--',r'Displaced$ \: m_{X} = 1.5 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')

plt.legend(loc='upper center', shadow=True, fontsize='small')

plt.grid(which='major',axis='x',alpha=0.5,linestyle='--',linewidth=0.5)
plt.grid(which='major',axis='y',alpha=0.5,linestyle='--',linewidth=0.5)

plt.savefig("plots/"+filename+".png",bbox_inches='tight')

plt.clf()#clear the current figure

##################################################################################################################################################################################
# the modified plot

print ("Plotting")

font = {'size'   : 18}
mpl.rc('font', family='serif')
plt.rc("text", usetex=True)
fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1) 
plt.rc("text", usetex=True)
         
#title = r'$\epsilon_{\textrm{det.}}=8.4\%$' 

#scenario_3
title = r'$g_{\tau \mu}$ vs. $g_{ee} (Belle II)$'
xlabel = r'$g_{ee}$ [GeV$^{-1}$]'
ylabel = r'$g_{\tau \mu}$ [GeV$^{-1}$]'
zlabel = r'$three_event$'
filename = "g_vs_g_3_mod"

#plot signal
#[1E-9, 1e-1, 1E-10, 1e-6]


plot_logx_logy_logz(data3[:,KEYS.index("g_EE_L:")], xlabel, 1e-9, 1e-1, data3[:,KEYS.index("g_TM_L:")], ylabel, 1e-10, 1e-6, data3[:,KEYS.index("Prompt_reallyvisibleX: ")], zlabel, title, filename, 'c','-',r'Prompt $\:m_{X} = 0.5 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data2[:,KEYS.index("g_EE_L:")], xlabel, 1e-9, 1e-1, data2[:,KEYS.index("g_TM_L:")], ylabel, 1e-10, 1e-6, data2[:,KEYS.index("Prompt_reallyvisibleX: ")], zlabel, title, filename, 'r','-',r'Prompt $\: m_{X} = 1.0 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data1[:,KEYS.index("g_EE_L:")], xlabel, 1e-9, 1e-1, data1[:,KEYS.index("g_TM_L:")], ylabel, 1e-10, 1e-6, data1[:,KEYS.index("Prompt_reallyvisibleX: ")], zlabel, title, filename, 'k','-',r'Prompt $\: m_{X} = 1.5 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')

plot_logx_logy_logz(data3[:,KEYS.index("g_EE_L:")], xlabel, 1e-9, 1e-1, data3[:,KEYS.index("g_TM_L:")], ylabel, 1e-10, 1e-6, data3[:,KEYS.index("Modified_reallyvisibleX: ")], zlabel, title, filename, 'c','--',r'Displaced $\: m_{X} = 0.5 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data2[:,KEYS.index("g_EE_L:")], xlabel, 1e-9, 1e-1, data2[:,KEYS.index("g_TM_L:")], ylabel, 1e-10, 1e-6, data2[:,KEYS.index("Modified_reallyvisibleX: ")], zlabel, title, filename, 'r','--',r'Displaced $\: m_{X} = 1.0 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data1[:,KEYS.index("g_EE_L:")], xlabel, 1e-9, 1e-1, data1[:,KEYS.index("g_TM_L:")], ylabel, 1e-10, 1e-6, data1[:,KEYS.index("Modified_reallyvisibleX: ")], zlabel, title, filename, 'k','--',r'Displaced $\: m_{X} = 1.5 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')


#plt.legend(loc='upper right',prop={'size': 9})
plt.legend(loc='lower left', shadow=True, fontsize='small')



#plot ctau

plt.grid(which='major',axis='x',alpha=0.5,linestyle='--',linewidth=0.5)
plt.grid(which='major',axis='y',alpha=0.5,linestyle='--',linewidth=0.5)

plt.savefig("plots/"+filename+".png",bbox_inches='tight')
plt.show()

#br vs ctau plot
plt.clf()#clear the current figure

print ("Plotting")

font = {'size'   : 18}
mpl.rc('font', family='serif')
plt.rc("text", usetex=True)
fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1) 
plt.rc("text", usetex=True)

#scenario_3
title = r'${\rm Br}(\tau \to X \, \mu)$ vs. $c\tau(X \to e \, \bar{e}) (Belle II)$'
xlabel = r'$c\tau(X \longrightarrow e \, \bar{e})$ [m]'
ylabel = r'${\rm Br}(\tau \longrightarrow X \, \mu)$'
zlabel = r'$three_event$'
filename = "br_vs_ct_3_mod"

#plt.cla()#clear the current axes
#[1E-6, 1e+2, 1E-10, 1e-3]
plot_logx_logy_logz(data3[:,KEYS.index("ctau [m]:")], xlabel, 1e-6, 1e+2, data3[:,KEYS.index("BRtau2xmu + BRtau2xe:")], ylabel, 1e-10, 1e-3, data3[:,KEYS.index("Prompt_reallyvisibleX: ")], zlabel, title, filename, 'c','-',r'Prompt$ \: m_{X} = 0.5 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data2[:,KEYS.index("ctau [m]:")], xlabel, 1e-6, 1e+2, data2[:,KEYS.index("BRtau2xmu + BRtau2xe:")], ylabel, 1e-10, 1e-3, data2[:,KEYS.index("Prompt_reallyvisibleX: ")], zlabel, title, filename, 'r','-',r'Prompt$ \: m_{X} = 1.0 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data1[:,KEYS.index("ctau [m]:")], xlabel, 1e-6, 1e+2, data1[:,KEYS.index("BRtau2xmu + BRtau2xe:")], ylabel, 1e-10, 1e-3, data1[:,KEYS.index("Prompt_reallyvisibleX: ")], zlabel, title, filename, 'k','-',r'Prompt$ \: m_{X} = 1.5 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')

plot_logx_logy_logz(data3[:,KEYS.index("ctau [m]:")], xlabel, 1e-6, 1e+2, data3[:,KEYS.index("BRtau2xmu + BRtau2xe:")], ylabel, 1e-10, 1e-3, data3[:,KEYS.index("Modified_reallyvisibleX: ")], zlabel, title, filename, 'c','--',r'Displaced$ \:m_{X} = 0.5 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data2[:,KEYS.index("ctau [m]:")], xlabel, 1e-6, 1e+2, data2[:,KEYS.index("BRtau2xmu + BRtau2xe:")], ylabel, 1e-10, 1e-3, data2[:,KEYS.index("Modified_reallyvisibleX: ")], zlabel, title, filename, 'r','--',r'Displaced$ \:m_{X} = 1.0 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data1[:,KEYS.index("ctau [m]:")], xlabel, 1e-6, 1e+2, data1[:,KEYS.index("BRtau2xmu + BRtau2xe:")], ylabel, 1e-10, 1e-3, data1[:,KEYS.index("Modified_reallyvisibleX: ")], zlabel, title, filename, 'k','--',r'Displaced$ \: m_{X} = 1.5 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')

plt.legend(loc='upper center', shadow=True, fontsize='small')

plt.grid(which='major',axis='x',alpha=0.5,linestyle='--',linewidth=0.5)
plt.grid(which='major',axis='y',alpha=0.5,linestyle='--',linewidth=0.5)

plt.savefig("plots/"+filename+".png",bbox_inches='tight')

#################################################################################################################################################################################

fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1)
