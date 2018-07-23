#
# Readout of membrane dynamics
#
import numpy as np
from math import *
import matplotlib.pyplot as plt
import random as pyrnd

clrs = ['blue','green']

T = 1010 #2000

dN = 200

Min = 200
Kin = 5
Nin = Min*Kin

RA = 100
gmax = 0.0025
gI = 0.00075
uk_min = 0.001
dendmodel = 'active'
learningrule = 'opt'
wiring = 'rewired_dropped_'
sd_bias = 1.0
release_prob = 1.0
ik = 0

tos = [0,100]
for toidx in range(len(tos)):
    strmbp = ''
    if tos[toidx] == 0:
        strmbp = 'data/nrn_simul_mbps_K' + str(Kin) + "_RA" + str(RA) \
                 + '_gm' + str(gmax) + '_gi' + str(gI) \
                 + '_dm_' + str(dendmodel) + '_lr_' + str(learningrule) \
                 + '_cw_' + str(wiring) + '_ukm_' + str(uk_min) \
                 + '_sdb_' + str(sd_bias) + '_pr_' + str(release_prob)\
                 + '_ik' + str(ik) + '.txt'
    else:
        strmbp = 'data/nrn_simul_mbpf_K' + str(Kin) + "_RA" + str(RA) \
                 + '_gm' + str(gmax) + '_gi' + str(gI) \
                 + '_dm_' + str(dendmodel) + '_lr_' + str(learningrule) \
                 + '_cw_' + str(wiring) + '_ukm_' + str(uk_min) \
                 + '_sdb_' + str(sd_bias) + '_pr_' + str(release_prob)\
                 + '_ik' + str(ik) + '.txt'

    tss = []; vss = []
    for nidx in range(dN):
        vss.append([])
    for line in open(strmbp,'r'):
        vtmps = line.split(" ")
        vidx = int(vtmps[0])
        if vidx == 0:
            tss.append( float(vtmps[1]) )
        vss[vidx].append( float( vtmps[2].strip("\n") ) )
        
    v1ss = []; v2ss = []
    for tidx in range(len(tss)):
        v1ss.append(0.0); v2ss.append(0.0)
    v1cnt = 0.0; v2cnt = 0.0
    for nidx in range(dN):
        #if prespikes[nidx] == prespike_const:
        if len(vss[nidx]) > 0 and max(vss[nidx]) < -30:
            if nidx < dN/2:
                v1cnt += 1
            else:
                v2cnt += 1

    for nidx in range(dN):
        #if prespikes[nidx] == prespike_const:
        if len(vss[nidx]) > 0 and max(vss[nidx]) < -30:
            if nidx < dN/2:
                for tidx in range(len(vss[nidx])):
                    v1ss[tidx] += vss[nidx][tidx]/float(v1cnt)
            else:
                for tidx in range(len(vss[nidx])):
                    v2ss[tidx] += vss[nidx][tidx]/float(v2cnt)

    #Plotting
    fig = plt.subplot(1,2,toidx+1)
    fig.spines['right'].set_visible(False)
    fig.spines['top'].set_visible(False)
    fig.yaxis.set_ticks_position('left')
    fig.xaxis.set_ticks_position('bottom')
    for nidx in range(dN):
        if pyrnd.random() < 0.1 and len(vss[nidx]) > 0:# and max(vss[nidx]) < -30:
            plt.plot(tss,vss[nidx],color=clrs[nidx/(dN/2)],alpha=0.4)
    plt.plot(tss,v2ss,color=clrs[1],linewidth=5.0)
    plt.plot(tss,v1ss,color=clrs[0],linewidth=5.0)
    plt.xlim(0,100)
    plt.xticks([10,35,60,85],[0,25,50,75],fontsize=20)
    plt.ylim(-80,-55)
    if toidx == 0:
        plt.yticks([-80,-75,-70,-65,-60],fontsize=20)
    else:
        plt.yticks([])
plt.show()
