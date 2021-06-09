#!/usr/bin/env python
print ("Loading Libraries")
import matplotlib
matplotlib.use('Agg')
import re
import os, sys
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

def plot_logx_logy_logz(xvalues, xlabel, yvalues, ylabel, zvalues, zlabel, title, filename, linecolor, contourstyle, linelabel, contourlevel, contourlinewidth, clabel_positions, clabelformat):
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
    plt.axis([1E-4, 1e+2, 1E-10, 1e-3])
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


    ax.plot(xg,-1*yg,linecolor,label=linelabel)
    cons2 = plt.contour(X, Y, Zm, [contourlevel] ,colors=linecolor,locator=mpl.ticker.LogLocator(), linewidths=contourlinewidth, linestyles = contourstyle)
#


    plt.scatter(map(lambda k: 10**k, x), map(lambda k: 10**k, y), s=-20, c=z, vmin=zmin, vmax=zmax, cmap = mpl.cm.Greys_r,norm=mpl.colors.LogNorm())


##############
#tau decay width limit 
sigmaBrtaupinu=0.057e-2
##############

########################################################################################################
########################################################################################################
########################################################################################################
print ("Reading in Data")

fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1) 
massstring = 'test'
for i in range(len(sys.argv))[1:]:
  resdircollection = sys.argv[i]
  print ("")
  for resdir in resdircollection.split(":"):
    #os.chdir(resdir)
        data=np.zeros((1,21))
        print (" - "+resdir)
        #counter = 0
        for filename in os.listdir(resdir):
            newdata=np.zeros((1,30))
            #print(filename)
            #counter = counter + 1
            #print(counter)
            #print (resdir+"/"+filename)
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "mX:" in l), None)
                 newdata[0,0]=line[line.find("mX:")+len("mX:"):]
            
            #scenario_0
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "g_CC_L:" in l), None)
                 newdata[0,1]=line[line.find("g_CC_L:")+len("g_CC_L:"):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "g_AB_L:" in l), None)
                 newdata[0,2]=line[line.find("g_AB_L:")+len("g_AB_L:"):]
            
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "BRtau2xmu:" in l), None)
                 newdata[0,3]=line[line.find("BRtau2xmu:")+len("BRtau2xmu:"):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "BRtau2xe:" in l), None)
                 newdata[0,4]=line[line.find("BRtau2xe:")+len("BRtau2xe:"):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "BRtau2xmu + BRtau2xe:" in l), None)
                 newdata[0,5]=line[line.find("BRtau2xmu + BRtau2xe:")+len("BRtau2xmu + BRtau2xe:"):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "produced tau excluding those daughters are also tau lepton:" in l), None)
                 newdata[0,6]=line[line.find("produced tau excluding those daughters are also tau lepton:")+len("produced tau excluding those daughters are also tau lepton:"):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "BRx2tautau:" in l), None)
                 newdata[0,7]=line[line.find("BRx2tautau:")+len("BRx2tautau:"):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "BRx2mumu:" in l), None)
                 newdata[0,8]=line[line.find("BRx2mumu:")+len("BRx2mumu:"):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "BRx2ee:" in l), None)
                 newdata[0,9]=line[line.find("BRx2ee:")+len("BRx2ee:"):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "BRx2gmgm:" in l), None)
                 newdata[0,10]=line[line.find("BRx2gmgm:")+len("BRx2gmgm:"):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "Total Gamma of X [GeV]: " in l), None)
                 newdata[0,11]=line[line.find("Total Gamma of X [GeV]: ")+len("Total Gamma of X [GeV]: "):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "ctau [m]:" in l), None)
                 newdata[0,12]=line[line.find("ctau [m]:")+len("ctau [m]:"):]
                
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "produced X:" in l), None)
                 newdata[0,13]=line[line.find("produced X:")+len("produced X:"):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "n_prompt_smallr: " in l), None)
                 newdata[0,14]=line[line.find("n_prompt_smallr: ")+len("n_prompt_smallr: "):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "n_prompt_d0z0: " in l), None)
                 newdata[0,15]=line[line.find("n_prompt_d0z0: ")+len("n_prompt_d0z0: "):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "n_prompt_baselineeff: " in l), None)
                 newdata[0,16]=line[line.find("n_prompt_baselineeff: ")+len("n_prompt_baselineeff: "):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "prompt small_r efficiencies (ind. and cum.): " in l), None)
                 newdata[0,17]=line[line.find("prompt small_r efficiencies (ind. and cum.): ")+len("prompt small_r efficiencies (ind. and cum.): "):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "prompt d0z0 efficiencies (ind. and cum.): " in l), None)
                 newdata[0,18]=line[line.find("prompt d0z0 efficiencies (ind. and cum.): ")+len("prompt d0z0 efficiencies (ind. and cum.): "):]
                
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "n_disp_fidvol: " in l), None)
                 newdata[0,19]=line[line.find("n_disp_fidvol: ")+len("n_disp_fidvol: "):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "n_disp_baselineeff: " in l), None)
                 newdata[0,20]=line[line.find("n_disp_baselineeff: ")+len("n_disp_baselineeff: "):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "n_disp_DispTrackEff: " in l), None)
                 newdata[0,21]=line[line.find("n_disp_DispTrackEff: ")+len("n_disp_DispTrackEff: "):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "displaced fiducial volume efficiencies (ind. and cum.): " in l), None)
                 newdata[0,22]=line[line.find("displaced fiducial volume efficiencies (ind. and cum.): ")+len("displaced fiducial volume efficiencies (ind. and cum.): "):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "displaced baseline efficiencies (ind. and cum.): " in l), None)
                 newdata[0,23]=line[line.find("displaced baseline efficiencies (ind. and cum.): ")+len("displaced baseline efficiencies (ind. and cum.): "):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "displaced displaced-tracking efficiencies (ind. and cum.): " in l), None)
                 newdata[0,24]=line[line.find("displaced displaced-tracking efficiencies (ind. and cum.): ")+len("displaced displaced-tracking efficiencies (ind. and cum.): "):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "reallyProducedX: " in l), None)
                 newdata[0,25]=line[line.find("reallyProducedX: ")+len("reallyProducedX: "):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "Prompt_reallyobservedX: " in l), None)
                 newdata[0,26]=line[line.find("Prompt_reallyobservedX: ")+len("Prompt_reallyobservedX: "):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "Prompt_reallyvisibleX: " in l), None)
                 newdata[0,27]=line[line.find("Prompt_reallyvisibleX: ")+len("Prompt_reallyvisibleX: "):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "Displaced_reallyobservedX: " in l), None)
                 newdata[0,28]=line[line.find("Displaced_reallyobservedX: ")+len("Displaced_reallyobservedX: "):]
            with open(resdir+"/"+filename,'r') as f: 
                 line = next((l for l in f if "Displaced_reallyvisibleX: " in l), None)
                 newdata[0,29]=line[line.find("Displaced_reallyvisibleX: ")+len("Displaced_reallyvisibleX: "):]
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
data2[where_are_NaNs] = 0


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
title = r'${\rm Br}(l_{\alpha} \longrightarrow X + l_{\beta})$ vs. $c\tau(X \longrightarrow l_{\alpha} + \overline{l}_{\alpha})$'
xlabel = r'$c\tau(X \longrightarrow l_{\alpha} + \overline{l}_{\alpha})$ [m]'
ylabel = r'${\rm Br}(l_{\alpha} \longrightarrow X + l_{\beta})$'
zlabel = r'$three_event$'
filename = "br_vs_ct_0_prompt"

#plot signal
plot_logx_logy_logz(data1[:,9], xlabel, data1[:,12], ylabel, data1[:,27], zlabel, title, filename, 'k','-', r'$m_{X} = 1.5 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data2[:,9], xlabel, data2[:,12], ylabel, data2[:,27], zlabel, title, filename, 'r','-', r'$m_{X} = 1.0 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')
plot_logx_logy_logz(data3[:,9], xlabel, data3[:,12], ylabel, data3[:,27], zlabel, title, filename, 'c','-', r'$m_{X} = 0.5 \: \textrm{GeV}$', 3 , 1 ,[(0.8,5e-8)],'3signal')


#use data files expoted by Mathematica where MovingAverage was used:
#signal311NoBG = genfromtxt("plots/data/signal_311_NoBG.dat");
#signal311NoBG1InvAb = genfromtxt("plots/data/signal_311_NoBG_1InvAb.dat");

#plt.plot(signal311NoBG[:,0], signal311NoBG[:,1],  color='blue',linewidth=1,linestyle='solid' ,label=r"$N_S = ~ \,$3, 50 ab$^{-1}$ (Belle II)" );
#plt.plot(signal311NoBG1InvAb[:,0], signal311NoBG1InvAb[:,1],  color='darkturquoise',linewidth=1,linestyle='solid' ,label=r"$N_S = ~ \,$3, 1 ab$^{-1}$ (Belle)" );

plt.legend(loc='lower left', shadow=True, fontsize='large')

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


#br vs ctau plot
plt.clf()#clear the current figure
#plt.cla()#clear the current axes

fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1)
