#!/usr/bin/env python
print ("Loading Libraries")
import matplotlib
matplotlib.use('Agg')
import re
import os, sys
os.system("rm -rf storage/results_g_g_0/0_max/.ipynb_checkpoints")
os.system("rm -rf storage/results_g_g_0/0_med/.ipynb_checkpoints")
os.system("rm -rf storage/results_g_g_0/0_min/.ipynb_checkpoints")
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

from operator import itemgetter


#This script needs to be executed in python2

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
        #print len(line[line.find(key)+len(key):].split(" "))
        #print line[line.find(key)+len(key):].split(" ")
        return line[line.find(key)+len(key):].split(" ")

def ZipAndSort(data, index1, index2):
    l1=[]
    zipped_lists = zip(data[:,index1], data[:,index2])
    sorted_zipped_lists = sorted(zipped_lists)
    sorted_list2 = [element for _, element in sorted_zipped_lists]
    l1=sorted(data[:,index1])
    return [l1, sorted_list2]

########################################################################################################
########################################################################################################
########################################################################################################
print ("Reading in Data")

fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1) 
#keys = ["mX:", "g_CC_L:", "g_AB_L:","BRtau2xmu:", "BRtau2xe:", "BRtau2xmu + BRtau2xe:", "produced tau excluding those daughters are also tau lepton:", "BRx2tautau:", "BRx2mumu:", "BRx2ee:", "BRx2gmgm:", "Total Gamma of X [GeV]: ", "ctau [m]:", "produced X:", "n_prompt_smallr: ", "n_prompt_d0z0: ", "n_prompt_baselineeff: ", "prompt small_r efficiencies (ind.): ", "prompt small_r efficiencies (cum.): ", "prompt d0z0 efficiencies (ind.): ", "prompt d0z0 efficiencies (cum.): ", "prompt baseline efficiencies (ind.): ", "prompt baseline efficiencies (cum.): ", "n_disp_fidvol: ", "n_disp_baselineeff: ", "n_disp_DispTrackEff: ", "displaced fiducial volume efficiencies (ind.): ", "displaced fiducial volume efficiencies (cum.): ", "displaced baseline efficiencies (ind.): ", "displaced baseline efficiencies (cum.): ", "displaced displaced-tracking efficiencies (ind.): ", "displaced displaced-tracking efficiencies (cum.): ", "reallyProducedX: ", "Prompt_reallyobservedX: ", "Prompt_reallyvisibleX: ", "Displaced_reallyobservedX: ", "Displaced_reallyvisibleX: "]
keys = ["mX:  ", "g_CC_L:  ", "g_AB_L:  ","BRtau2xmu: ", "BRtau2xe: ", "BRtau2xmu + BRtau2xe: ", "produced tau excluding those daughters are also tau lepton: ", "BRx2tautau: ", "BRx2mumu: ", "BRx2ee: ", "BRx2gmgm: ", "Total Gamma of X [GeV]: ", "ctau [m]: ", "produced X: ", "n_prompt_smallr: ", "n_prompt_d0z0: ", "n_prompt_baselineeff: ", "prompt small_r efficiencies (ind. and cum.): ", "prompt d0z0 efficiencies (ind. and cum.): ", "prompt baseline efficiencies (ind. and cum.): ", "n_disp_fidvol: ", "n_disp_baselineeff: ", "n_disp_DispTrackEff: ", "displaced fiducial volume efficiencies (ind. and cum.): ", "displaced baseline efficiencies (ind. and cum.): ", "displaced displaced-tracking efficiencies (ind. and cum.): ", "n_mod_DispTrackEff: ", "modified efficiencies: ", "reallyProducedX: ", "Prompt_reallyobservedX: ", "Prompt_reallyvisibleX: ", "Displaced_reallyobservedX: ", "Displaced_reallyvisibleX: ", "Modified_reallyobservedX: ", "Modified_reallyvisibleX: "]

for i in range(len(sys.argv))[1:]:
    resdircollection = sys.argv[i]
    print ("")
    for resdir in resdircollection.split(":"):
    #os.chdir(resdir)
        data=np.zeros((1,53))
        print (" - "+resdir)
        #counter = 0
        for filename in os.listdir(resdir):
            newdata=np.zeros((1,53))
            #print filename
            ndata = 0
            KEYS = []
            for k in keys:
                #if "(ind. and cum.)" in k:
                    #test = ReadLine(resdir, filename, k)
                    #newdata[0,ndata] = test.split()[0]
                    #KEYS.append(k[:k.find("(ind. and cum.)")] + "(ind.): ")
                    #newdata[0,ndata+1] = test.split()[1]
                    #KEYS.append(k[:k.find("(ind. and cum.)")] + "(cum.): ")
                    #ndata +=2
                #else:
                lengthReadLine = len(ReadLine(resdir, filename, k) ) 
                for n in range(ndata, ndata+lengthReadLine  ):
                    #print "k: "+str(k)+", "+"n: "+str(n)
                    newdata[0,n] = float(ReadLine(resdir, filename, k)[n-ndata])
                    if ("(ind. and cum.)" in k)   and   (n==ndata)   :
                        KEYS.append( k.replace("ind. and cum.","ind.")   )
                    elif ("(ind. and cum.)" in k) and (n==ndata+1)   :
                        KEYS.append( k.replace("ind. and cum.","cum.")   )
                    else: KEYS.append(k)
                    #print KEYS
                ndata += lengthReadLine
                    
            #print(filename)
            #counter = counter + 1
            #print(counter)
            #print (resdir+"/"+filename)
           
            # data = data.astype(np.double)
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
#np.savetxt("data1.dat", data1, fmt="%s")
#np.savetxt("data2.dat", data2, fmt="%s")
#np.savetxt("data3.dat", data3, fmt="%s")

# os.chdir("..")	   

#print(KEYS)

########################################################################################################################################################################################################################################################################################################################################################

#prom ind
print ("Plotting")

font = {'size'   : 18}
mpl.rc('font', family='serif')
plt.rc("text", usetex=True)
fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1) 
plt.rc("text", usetex=True)
         
#title = r'$\epsilon_{\textrm{det.}}=8.4\%$' 

#scenario_0
title = r'Prompt individual $\epsilon \,$ vs. $c\tau$'
xlabel = r'$c\tau(X \longrightarrow l_{\alpha} + \overline{l}_{\alpha})$ [m]'
ylabel = r'Efficiency $\: \epsilon$'
zlabel = r'$three_event$'
filename = "ef_vs_ct_prom_ind_0"

#plot signal
#[1E-12, 1e-4, 1E-12, 1e-6]

plt.xscale("log")
plt.yscale("log")
plt.xlabel(xlabel,fontsize = 16)
plt.ylabel(ylabel,fontsize = 16)
plt.title(title,fontsize = 16)
#plt.axis([1e-3, 1e2, 4e-4, 1])
plt.xlim(1e-4,6e2)
plt.ylim(1e-4,1.5e0)

plt.plot(ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("prompt small_r efficiencies (ind.): "))[0], ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("prompt small_r efficiencies (ind.): "))[1],'c-', label = '$r_{min} \: m_{X} = 0.5$GeV')
plt.plot(ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("prompt small_r efficiencies (ind.): "))[0], ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("prompt small_r efficiencies (ind.): "))[1],'r-', label = '$r_{min} \: m_{X} = 1.0$GeV')
plt.plot(ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("prompt small_r efficiencies (ind.): "))[0], ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("prompt small_r efficiencies (ind.): "))[1],'k-', label = '$r_{min} \: m_{X} = 1.5$GeV')

plt.plot(ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("prompt d0z0 efficiencies (ind.): "))[0], ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("prompt d0z0 efficiencies (ind.): "))[1],'c--', label = '$d_{0}$-$z_{0} \: m_{X} = 0.5$GeV')
plt.plot(ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("prompt d0z0 efficiencies (ind.): "))[0], ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("prompt d0z0 efficiencies (ind.): "))[1],'r--', label = '$d_{0}$-$z_{0} \: m_{X} = 1.0$GeV')
plt.plot(ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("prompt d0z0 efficiencies (ind.): "))[0], ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("prompt d0z0 efficiencies (ind.): "))[1],'k--', label = '$d_{0}$-$z_{0} \: m_{X} = 1.5$GeV')

plt.plot(ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("prompt baseline efficiencies (ind.): "))[0], ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("prompt baseline efficiencies (ind.): "))[1],'c:', label = 'base. $m_{X} = 0.5$GeV')
plt.plot(ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("prompt baseline efficiencies (ind.): "))[0], ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("prompt baseline efficiencies (ind.): "))[1],'r:', label = 'base. $m_{X} = 1.0$GeV')
plt.plot(ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("prompt baseline efficiencies (ind.): "))[0], ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("prompt baseline efficiencies (ind.): "))[1],'k:', label = 'base. $m_{X} = 1.5$GeV')

plt.legend(loc='lower left', shadow=True, fontsize='xx-small')

plt.grid(which='major',axis='x',alpha=0.5,linestyle='--',linewidth=0.5)
plt.grid(which='major',axis='y',alpha=0.5,linestyle='--',linewidth=0.5)

plt.savefig("plots/"+filename+".png",bbox_inches='tight')
plt.show()

plt.clf()#clear the current figure
############################################################################################################################################################################
#pron cum
print ("Plotting")

font = {'size'   : 18}
mpl.rc('font', family='serif')
plt.rc("text", usetex=True)
fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1) 
plt.rc("text", usetex=True)

#scenario_0
title = r'Prompt cumulative $\epsilon \,$ vs. $c\tau$'
xlabel = r'$c\tau(X \longrightarrow l_{\alpha} + \overline{l}_{\alpha})$ [m]'
ylabel = r'Efficiency $\: \epsilon$'
zlabel = r'$three_event$'
filename = "ef_vs_ct_prom_cum_0"

#plot signal
#[1E-12, 1e-4, 1E-12, 1e-6]

plt.xscale("log")
plt.yscale("log")
plt.xlabel(xlabel,fontsize = 16)
plt.ylabel(ylabel,fontsize = 16)
plt.title(title,fontsize = 16)
#plt.axis([1e-3, 1e2, 4e-4, 1])
plt.xlim(1e-4,6e2)
plt.ylim(1e-4,1.5e0)

plt.plot(ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("prompt small_r efficiencies (cum.): "))[0], ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("prompt small_r efficiencies (cum.): "))[1],'c-', label = '$r_{min} \: m_{X} = 0.5$GeV')
plt.plot(ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("prompt small_r efficiencies (cum.): "))[0], ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("prompt small_r efficiencies (cum.): "))[1],'r-', label = '$r_{min} \: m_{X} = 1.0$GeV')
plt.plot(ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("prompt small_r efficiencies (cum.): "))[0], ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("prompt small_r efficiencies (cum.): "))[1],'k-', label = '$r_{min} \: m_{X} = 1.5$GeV')

plt.plot(ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("prompt d0z0 efficiencies (cum.): "))[0], ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("prompt d0z0 efficiencies (cum.): "))[1],'c--', label = '$d_{0}$-$z_{0} \: m_{X} = 0.5$GeV')
plt.plot(ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("prompt d0z0 efficiencies (cum.): "))[0], ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("prompt d0z0 efficiencies (cum.): "))[1],'r--', label = '$d_{0}$-$z_{0} \: m_{X} = 1.0$GeV')
plt.plot(ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("prompt d0z0 efficiencies (cum.): "))[0], ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("prompt d0z0 efficiencies (cum.): "))[1],'k--', label = '$d_{0}$-$z_{0} \: m_{X} = 1.5$GeV')

plt.plot(ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("prompt baseline efficiencies (cum.): "))[0], ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("prompt baseline efficiencies (cum.): "))[1],'c:', label = 'base. $m_{X} = 0.5$GeV')
plt.plot(ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("prompt baseline efficiencies (cum.): "))[0], ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("prompt baseline efficiencies (cum.): "))[1],'r:', label = 'base. $m_{X} = 1.0$GeV')
plt.plot(ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("prompt baseline efficiencies (cum.): "))[0], ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("prompt baseline efficiencies (cum.): "))[1],'k:', label = 'base. $m_{X} = 1.5$GeV')

plt.legend(loc='lower left', shadow=True, fontsize='xx-small')

plt.grid(which='major',axis='x',alpha=0.5,linestyle='--',linewidth=0.5)
plt.grid(which='major',axis='y',alpha=0.5,linestyle='--',linewidth=0.5)

plt.savefig("plots/"+filename+".png",bbox_inches='tight')
plt.show()

plt.clf()#clear the current figure

############################################################################################################################################################################
#disp ind
print ("Plotting")

font = {'size'   : 18}
mpl.rc('font', family='serif')
plt.rc("text", usetex=True)
fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1) 
plt.rc("text", usetex=True)

#scenario_0
title = r'Displaced individual $\epsilon \,$ vs. $c\tau$'
xlabel = r'$c\tau(X \longrightarrow l_{\alpha} + \overline{l}_{\alpha})$ [m]'
ylabel = r'Efficiency $\: \epsilon$'
zlabel = r'$three_event$'
filename = "ef_vs_ct_disp_ind_0"

#plot signal
#[1E-12, 1e-4, 1E-12, 1e-6]

plt.xscale("log")
plt.yscale("log")
plt.xlabel(xlabel,fontsize = 16)
plt.ylabel(ylabel,fontsize = 16)
plt.title(title,fontsize = 16)
#plt.axis([1e-3, 1e2, 4e-4, 1])
plt.xlim(1e-4,1e4)
plt.ylim(1e-4,1.5e0)

plt.plot(ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("displaced fiducial volume efficiencies (ind.): "))[0], ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("displaced fiducial volume efficiencies (ind.): "))[1],'c-', label = 'f-vol. $m_{X} = 0.5$GeV')
plt.plot(ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("displaced fiducial volume efficiencies (ind.): "))[0], ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("displaced fiducial volume efficiencies (ind.): "))[1],'r-', label = 'f-vol. $m_{X} = 1.0$GeV')
plt.plot(ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("displaced fiducial volume efficiencies (ind.): "))[0], ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("displaced fiducial volume efficiencies (ind.): "))[1],'k-', label = 'f-vol. $m_{X} = 1.5$GeV')

plt.plot(ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("displaced baseline efficiencies (ind.): "))[0], ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("displaced baseline efficiencies (ind.): "))[1],'c--', label = 'base. $m_{X} = 0.5$GeV')
plt.plot(ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("displaced baseline efficiencies (ind.): "))[0], ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("displaced baseline efficiencies (ind.): "))[1],'r--', label = 'base. $m_{X} = 1.0$GeV')
plt.plot(ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("displaced baseline efficiencies (ind.): "))[0], ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("displaced baseline efficiencies (ind.): "))[1],'k--', label = 'base. $m_{X} = 1.5$GeV')

plt.plot(ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("displaced displaced-tracking efficiencies (ind.): "))[0], ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("displaced displaced-tracking efficiencies (ind.): "))[1],'c:', label = 'track. $m_{X} = 0.5$GeV')
plt.plot(ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("displaced displaced-tracking efficiencies (ind.): "))[0], ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("displaced displaced-tracking efficiencies (ind.): "))[1],'r:', label = 'track. $m_{X} = 1.0$GeV')
plt.plot(ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("displaced displaced-tracking efficiencies (ind.): "))[0], ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("displaced displaced-tracking efficiencies (ind.): "))[1],'k:', label = 'track. $m_{X} = 1.5$GeV')

plt.legend(loc='lower center', shadow=True, fontsize='xx-small')

plt.grid(which='major',axis='x',alpha=0.5,linestyle='--',linewidth=0.5)
plt.grid(which='major',axis='y',alpha=0.5,linestyle='--',linewidth=0.5)

plt.savefig("plots/"+filename+".png",bbox_inches='tight')
plt.show()

plt.clf()#clear the current figure
############################################################################################################################################################################
#disp cum
print ("Plotting")

font = {'size'   : 18}
mpl.rc('font', family='serif')
plt.rc("text", usetex=True)
fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1) 
plt.rc("text", usetex=True)

#scenario_0
title = r'Displaced cumulative $\epsilon \,$ vs. $c\tau$'
xlabel = r'$c\tau(X \longrightarrow l_{\alpha} + \overline{l}_{\alpha})$ [m]'
ylabel = r'Efficiency $\: \epsilon$'
zlabel = r'$three_event$'
filename = "ef_vs_ct_disp_cum_0"

#plot signal
#[1E-12, 1e-4, 1E-12, 1e-6]

plt.xscale("log")
plt.yscale("log")
plt.xlabel(xlabel,fontsize = 16)
plt.ylabel(ylabel,fontsize = 16)
plt.title(title,fontsize = 16)
#plt.axis([1e-3, 1e2, 4e-4, 1])
plt.xlim(1e-4,1e4)
plt.ylim(1e-4,1.5e0)

plt.plot(ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("displaced fiducial volume efficiencies (cum.): "))[0], ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("displaced fiducial volume efficiencies (cum.): "))[1],'c-', label = 'f-vol. $m_{X} = 0.5$GeV')
plt.plot(ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("displaced fiducial volume efficiencies (cum.): "))[0], ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("displaced fiducial volume efficiencies (cum.): "))[1],'r-', label = 'f-vol. $m_{X} = 1.0$GeV')
plt.plot(ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("displaced fiducial volume efficiencies (cum.): "))[0], ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("displaced fiducial volume efficiencies (cum.): "))[1],'k-', label = 'f-vol. $m_{X} = 1.5$GeV')

plt.plot(ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("displaced baseline efficiencies (cum.): "))[0], ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("displaced baseline efficiencies (cum.): "))[1],'c--', label = 'base. $m_{X} = 0.5$GeV')
plt.plot(ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("displaced baseline efficiencies (cum.): "))[0], ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("displaced baseline efficiencies (cum.): "))[1],'r--', label = 'base. $m_{X} = 1.0$GeV')
plt.plot(ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("displaced baseline efficiencies (cum.): "))[0], ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("displaced baseline efficiencies (cum.): "))[1],'k--', label = 'base. $m_{X} = 1.5$GeV')

plt.plot(ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("displaced displaced-tracking efficiencies (cum.): "))[0], ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("displaced displaced-tracking efficiencies (cum.): "))[1],'c:', label = 'track. $m_{X} = 0.5$GeV')
plt.plot(ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("displaced displaced-tracking efficiencies (cum.): "))[0], ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("displaced displaced-tracking efficiencies (cum.): "))[1],'r:', label = 'track. $m_{X} = 1.0$GeV')
plt.plot(ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("displaced displaced-tracking efficiencies (cum.): "))[0], ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("displaced displaced-tracking efficiencies (cum.): "))[1],'k:', label = 'track. $m_{X} = 1.5$GeV')

plt.legend(loc='upper right', shadow=True, fontsize='xx-small')

plt.grid(which='major',axis='x',alpha=0.5,linestyle='--',linewidth=0.5)
plt.grid(which='major',axis='y',alpha=0.5,linestyle='--',linewidth=0.5)

plt.savefig("plots/"+filename+".png",bbox_inches='tight')
plt.show()

plt.clf()#clear the current figure
################################################################################################################################################################################
#disp modified
print ("Plotting")

font = {'size'   : 18}
mpl.rc('font', family='serif')
plt.rc("text", usetex=True)
fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1) 
plt.rc("text", usetex=True)

#scenario_0
title = r'Displaced modified $\epsilon \,$ vs. $c\tau$'
xlabel = r'$c\tau(X \longrightarrow l_{\alpha} + \overline{l}_{\alpha})$ [m]'
ylabel = r'Efficiency $\: \epsilon$'
zlabel = r'$three_event$'
filename = "ef_vs_ct_mod_0"

#plot signal
#[1E-12, 1e-4, 1E-12, 1e-6]

plt.xscale("log")
plt.yscale("log")
plt.xlabel(xlabel,fontsize = 16)
plt.ylabel(ylabel,fontsize = 16)
plt.title(title,fontsize = 16)
#plt.axis([1e-3, 1e2, 4e-4, 1])
plt.xlim(1e-4,1e4)
plt.ylim(1e-4,1.5e0)

plt.plot(ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("modified efficiencies: "))[0], ZipAndSort(data3,KEYS.index("ctau [m]: "),KEYS.index("modified efficiencies: "))[1],'c-', label = 'mod. $m_{X} = 0.5$GeV')
plt.plot(ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("modified efficiencies: "))[0], ZipAndSort(data2,KEYS.index("ctau [m]: "),KEYS.index("modified efficiencies: "))[1],'r-', label = 'mod. $m_{X} = 1.0$GeV')
plt.plot(ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("modified efficiencies: "))[0], ZipAndSort(data1,KEYS.index("ctau [m]: "),KEYS.index("modified efficiencies: "))[1],'k-', label = 'mod. $m_{X} = 1.5$GeV')

plt.legend(loc='upper right', shadow=True, fontsize='xx-small')

plt.grid(which='major',axis='x',alpha=0.5,linestyle='--',linewidth=0.5)
plt.grid(which='major',axis='y',alpha=0.5,linestyle='--',linewidth=0.5)

plt.savefig("plots/"+filename+".png",bbox_inches='tight')
plt.show()

plt.clf()#clear the current figure


fig = plt.figure(num=None, figsize=(5,4), dpi=300, facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1)
