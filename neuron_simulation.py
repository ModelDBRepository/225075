# ----------------------------------------------------------
# Main simulation control
#
# Naoki Hiratani, N.Hiratani@gmail.com
#
# Based on the model by Dr. Tiago Branco (Smith et al., 2013)
# 
# ----------------------------------------------------------

import numpy as np
import pickle
import time

import libcell as lb
#import saveClass as sc

from math import *
import matplotlib.pyplot as plt
import random as pyrnd
import sys

#----------------------------------------------------------------------------
# Functions and Classes

#calculation of maximum membrane potential difference
def calc_vks(v,Nin):
    vkmin = min(v)
    vkmaxs = []
    vks = []
    for tidx in range(2,len(v)):
        if (v[tidx-2] - v[tidx-1])*(v[tidx-1] - v[tidx]) < 0.0:
            if v[tidx-1] > vkmin + 0.01:
                vks.append( v[tidx-1] - v[-1] )
    if len(vks) != Nin:
        print "the program failed to estimate vks"
        
    return vks

#calculation of unit EPSPs
def estimate_vkappas(model):
    locs = lb.calc_locs_kappa(model)
    lb.add_Estims(model,locs)
    lb.add_AMPAsyns(model,locs,model.gmax)
    model.tstop = 100
    strvkp = "data/spf4d4_vkappas_RA" + str(model.RA) + "_gm" + str(model.gmax) + "_dm_" + str(model.dendmodel) + ".txt"
    fvkp  = open(strvkp,'w')
    vkappas = []
    for lidx in range(len(locs)):
        for l2idx in range(len(locs)):
            model.Estimlist[l2idx].number = 0
        model.Estimlist[lidx].number = 1
        ttmp, vtmp = lb.simulate(model)
        dvtmp = max(vtmp) - vtmp[-1]
        print lidx,dvtmp
        vkappas.append( dvtmp )
        fvkp.write( str(locs[lidx][0]) + " " + str(locs[lidx][1]) + " " + str(dvtmp) + "\n" )
        fvkp.flush()
    plt.hist(vkappas)
    plt.show()

def estimate_vks(model, vks):
    model.sinterval = 1000.0
    for i in range(model.Nin):
        model.Estimlist[i].interval = model.sinterval
        model.Estimlist[i].number = 1
        model.Estimlist[i].start = 500 + i*model.sinterval
        model.Estimlist[i].noise = 0
    model.tstop = model.Nin*model.sinterval + 1000
    ttmp, vtmp = lb.simulate(model)

    return calc_vks(vtmp,model.Nin)

def pickup_vks(model,locs): #pick up the value of vk from the file 'strvkp'
    strvkp = "data/nrn_simul_vkappas_RA" + str(model.RA) + "_gm" + str(model.gmax) + "_dm_" + str(model.dendmodel) + ".txt"
    model.dendidxs = []
    model.dendpos = []
    model.vkps = []
    for line in open(strvkp,'r'):
        vktmp = line.split(" ")
        model.dendidxs.append( int(vktmp[0]) )
        model.dendpos.append( float(vktmp[1]) )
        model.vkps.append( float(vktmp[2].strip("\n")) )

    vks = []
    for nidx in range(len(locs)):
        lidx = locs[nidx][0]
        ptmp = locs[nidx][1]
        dptmp = 1.0
        vks.append(0.0)
        for n2idx in range(len(model.dendidxs)):
            if model.dendidxs[n2idx] == lidx:
                if abs( ptmp - model.dendpos[n2idx] ) <= dptmp:
                    vks[nidx] = model.vkps[n2idx]
                    dptmp = abs( ptmp - model.dendpos[n2idx] )
        if vks[nidx] == 0.0:
            print "an error occurred at nidx = " + str(nidx) + ", while reading the vkappa file."
    #print vks
    return vks

def pick_a_vk(model,loc): 
    vktmp = 0.0
    lidx = loc[0]
    ptmp = loc[1]
    dptmp = 1.0
    for nidx in range(len(model.dendidxs)):
        if model.dendidxs[nidx] == lidx:
            if abs( ptmp - model.dendpos[nidx] ) <= dptmp:
                vktmp = model.vkps[nidx]
                dptmp = abs( ptmp - model.dendpos[nidx] )
    return vktmp

def estimate_dks(model,locs): #distance from the soma to the synapses
    dks = []
    for lidx in range(len(locs)):
        loc = locs[lidx]
        dks.append( model.dists[loc[0]] + loc[1]*model.dends[loc[0]].L )
    return dks    

def calc_pvks(vks): #distribution of vk
    pvks = []
    vkmin = min(vks)
    vkmax = max(vks)
    dvk = (vkmax - vkmin)/10.0
    for i in range(len(vks)):
        pvks.append(0.0)
        for i2 in range(len(vks)):
            if vks[i]-0.5*dvk < vks[i2] and vks[i2] < vks[i]+0.5*dvk:
                pvks[i] += 1.0/float(len(vks)*dvk)
    return pvks

def initialize_ukgks(model,pvks):
    uks = []
    gks = []
    Zuks = []
    for j in range(model.Min):
        Zuks.append(0.0)
    for i in range(model.Nin):
        uks.append(1.0/float(pvks[i]))
        Zuks[ model.preidx[i] ] += uks[i]
        gks.append(1.0)
    for i in range(model.Nin):
        uks[i] = uks[i]/Zuks[ model.preidx[i] ]
    return uks,gks

#update of spine sizes with near-optimal learning
def update_uks(model, vks, uks, pre_spikes, rewiring, dropping):
    gamma = model.gamma
    uk_min = model.uk_min; uk_max = model.uk_max
    vkmin = model.vkmin
    min_rate = model.min_rate
    prev_uks = uks

    connec_exists = []
    for j in range(model.Min):
        connec_exists.append(0)
    for i in range(model.Nin):
        jidx = model.preidx[i];
        if uks[i] > 0.0:
            connec_exists[jidx] += 1

    Zuks = np.zeros((model.Min))
    for i in range(model.Nin):
        jidx = model.preidx[i]
        if connec_exists[jidx] > 0:
            tmp_rate = min_rate*exp( gamma*vks[i] )
            uks[i] *= (tmp_rate**pre_spikes[i])*exp(-tmp_rate)/factorial( pre_spikes[i] )
            Zuks[jidx] += uks[i]
    for i in range(model.Nin):
        jidx = model.preidx[i]
        if connec_exists[jidx] > 0:
            if Zuks[jidx] > 0.0:
                uks[i] = uks[i]/Zuks[jidx]
            else:
                uks[i] = prev_uks[i]
                if uks[i] > 0.0:
                    print model.preidx[i], i, vks[i], uks[i]
            if prev_uks[i] > 0.0 and uks[i] == 0.0:
                uks[i] = 0.1*uk_min
            if (not rewiring) and (not dropping) and uks[i] < uk_min:
                uks[i] = uk_min
            if uks[i] > uk_max:
                uks[i] = uk_max
    return uks

#update of spine sizes without pre-neuron specific normalization
def update_uks_apr(model, vks, uks, prespikes, rewiring, dropping):
    Nin = model.Nin; Min = model.Min; Kin = model.Kin
    gamma = model.gamma
    uk_min = model.uk_min; uk_max = model.uk_max
    vkmin = model.vkmin; vo = model.vo
    min_rate = model.min_rate

    duks = []
    Zdtot = 0.0
    for i in range(Nin):
        jidx = model.preidx[i]
        rhojk = min_rate*exp(gamma*vks[i])

        vojk = (1.0 - uks[i])*vo + uks[i]*vks[i]
        rhoo = min_rate*exp(gamma*vojk)
        duks.append( ( (rhojk/rhoo)**prespikes[i] )*exp( -(rhojk-rhoo) ) )
        Zdtot += log(duks[i])/float(Nin)

    epZdtot = exp(Zdtot)
    for i in range(Nin):
        uks[i] = uks[i]*duks[i]/epZdtot
        if uks[i] > 0.5*uk_max:
            uks[i] = 0.5*uk_max
        if (not rewiring) and (not dropping) and uks[i] < uk_min:
            uks[i] = uk_min
    return uks

#calculation of EPSP area
def calc_integ_v(vtmp,drt):
    integ_v = 0.0
    vmin = vtmp[400] #min(vtmp) # membrane potential at t=10ms
    for tidx in range(len(vtmp)):
        integ_v += (vtmp[tidx] - vmin)*drt
    return integ_v

def evaluate_perf(model,prespikes,ttmp,vtmp):
    epsilon = model.epsilon
    opt_pv = 0.0; opt_nv = 0.0;
    for j in range(model.Min):
        opt_pv += prespikes[j]*(model.log_p_rates[j] + log(epsilon))
        opt_nv += prespikes[j]*(model.log_n_rates[j] + log(epsilon))

    drt = ttmp[1] - ttmp[0]
    integ_v = calc_integ_v(vtmp,drt)
    depol_v = max(vtmp) - vtmp[400]

    total_prespikes = 0
    for prespike in prespikes:
        total_prespikes += prespike

    perftmp = [opt_pv, opt_nv, integ_v, depol_v, total_prespikes]
    return perftmp

def is_match(t,trs):
    tof = False
    for tr in trs:
        if t == tr:
            tof = True
    return tof

def simul(Kin, gmax, gI, uk_min, sd_bias, release_prob, ik): 
    model = lb.L23() # Layer 2/3 model is selected
    dendActive = True # with passive/active dendrrite
    opt_update = True  #whether update rule is optimal or approximate
    rewiring = True  #if the model includes rewiring or not
    dropping = True #if the model includes uncompensated elimination 
    membrane_recording = True #if membrane dynamics is recorded or not (generate a heavy file if true)
    prespike_recording = False #if the prespikes are recorded or not

    # active channels are distributed over both axon and dendrite
    lb.init_active(model, axon=True, soma=True, dend=dendActive, dendNa=True, dendCa=True)
    lb.init_params(model, Kin, gmax, gI, uk_min, sd_bias, release_prob)

    model.dendmodel = 'active' if dendActive else 'passive'
    model.learningrule = 'opt' if opt_update else 'apr'
    if rewiring or dropping:
        model.connections = ''
        if rewiring:
            model.connections = model.connections + 'rewired_'
        if dropping:
            model.connections = model.connections + 'dropped_'
    else:
        model.connections = 'fixed'

    #CAUTION: 'estimate_vkappas' resets vkappa
    #estimate_vkappas(model)
    
    #locs = lb.calc_locs(model) #uniformly locate synapses
    locs = lb.calc_locs_random(model) # randomly locate synapses
    lb.allocate_syns(model)
    #lb.plot_morphology(model,locs)

    lb.add_Estims(model,locs)
    lb.add_AMPAsyns(model, locs, model.gmax)

    inh_locs = lb.calc_inh_locs_uniform(model)
    lb.add_Istims(model, inh_locs)
    lb.add_GABAsyns(model, inh_locs, model.gI)

    dks = estimate_dks(model,locs) #estimate the distances from the soma
    vks = pickup_vks(model,locs) #readout values of vks
    mean_pre_rates = np.full((model.Nin), model.spn_rate)
    
    pvks = calc_pvks(vks) #prior distribution of vks

    uks, gks = initialize_ukgks(model,pvks)
    lb.generate_pre_tuning(model)
    lb.generate_input_rates(model)

    model.vkmax = max(vks); model.vkmin = min(vks)
    model.gamma = max(model.log_p_rates)/model.vkmax
    model.vo = 1.5*model.vkmin#
    
    fstrtmp = 'K' + str(model.Kin) + "_RA" + str(model.RA) + '_gm' + str(model.gmax) + '_gi' + str(model.gI) + \
              '_dm_' + str(model.dendmodel) + '_lr_' + str(model.learningrule) + '_cw_' + str(model.connections) \
              + '_ukm_' + str(uk_min) + '_sdb_' + str(sd_bias) + '_pr_' + str(release_prob) + '_ik' + str(ik) + '.txt'
    strperf = 'data/nrn_simul_perf_' + fstrtmp
    fperf = open(strperf,'w')
    strduvk = 'data/nrn_simul_duvk_' + fstrtmp
    fduvk = open(strduvk,'w')
    strtune = 'data/nrn_simul_tune_' + fstrtmp
    ftune = open(strtune,'w')
    for j in range(model.Min):
        ftune.write(str(j) + ' ' + str(model.phis[j]) + ' ' + str(model.thetas[j]) + ' ' + str(model.rfdists[j]) + ' ' + str(model.log_p_rates[j]/model.gamma) + ' ' + str(model.log_n_rates[j]/model.gamma) + '\n')
    ftune.flush()
    if prespike_recording:
        strprespike = 'data/nrn_simul_prespike_' + fstrtmp
        fprespike = open(strprespike,'w')

    model.tstop = 100
    trs = [0,100]
    perfs = []
    for q in range(3):
        perfs.append([])
    for t in range(model.trials):
        if is_match(t,trs):#evaluation
            if membrane_recording:
                if t == trs[0]:
                    strmbps = 'data/nrn_simul_mbps_' + fstrtmp
                    fmbps = open(strmbps,'w')
                elif t == trs[-1]:
                    strmbpf = 'data/nrn_simul_mbpf_' + fstrtmp
                    fmbpf = open(strmbpf,'w')
            
            for i in range(model.Nin):
                fduvk.write(str(locs[i][0]) + " " + str(locs[i][1]) + " " + str(model.preidx[i]) + " " + str(dks[i]) + " " + str(uks[i]) + " " + str(vks[i]) + "\n")
                            
            dist_bt_syns = lb.calc_dist_bt_syns(model, locs)
            for didx in range(len(dist_bt_syns)):
                fduvk.write( str(dist_bt_syns[didx]) + ' ' )
            fduvk.write( '\n' ); fduvk.flush()
            
            #test phase
            for t2 in range(200): 
                true_pre_spikes = []; pre_spikes = []
                if t2 < 100:
                    true_pre_spikes, pre_spikes = lb.generate_prespikes(model, model.p_rates)
                else:
                    true_pre_spikes, pre_spikes = lb.generate_prespikes(model, model.n_rates)
                if prespike_recording:
                    for j in range(model.Min):
                        fprespike.write( str(t) + " " + str(t2) + " " + str(j) + " " + str(true_pre_spikes[j]) + "\n" )
                    fprespike.flush()

                lb.generate_inputs(model, uks, pre_spikes)
                lb.generate_inh_inputs(model, true_pre_spikes)
                
                ttmp, vtmp = lb.simulate(model)
                
                perftmp = evaluate_perf(model,true_pre_spikes,ttmp,vtmp)
                strtmp = str(t) + ' ' + str(perftmp[0]) + ' ' + str(perftmp[1]) + ' ' + str(perftmp[2]) + ' ' + str(perftmp[3]) + ' ' + str(perftmp[4]) + '\n'
                fperf.write(strtmp)
                if membrane_recording:
                    if t == trs[0]:
                        for t3 in range(len(ttmp)):
                            strtmp = str(t2) + ' ' + str(ttmp[t3]) + ' ' + str(vtmp[t3]) + '\n'
                            fmbps.write( strtmp )
                        fmbps.flush()
                    if t == trs[-1]:
                        for t3 in range(len(ttmp)):
                            strtmp = str(t2) + ' ' + str(ttmp[t3]) + ' ' + str(vtmp[t3]) + '\n'
                            fmbpf.write( strtmp )
                        fmbpf.flush()
            
        #training phase
        true_pre_spikes, pre_spikes = lb.generate_prespikes(model, model.p_rates)
        mean_pre_rates = np.multiply(1.0 - 1.0/model.taumr, mean_pre_rates) + np.multiply(1.0/model.taumr, pre_spikes)
        if opt_update:
            uks = update_uks(model, vks, uks, pre_spikes, rewiring, dropping)
        else:
            uks = update_uks_apr(model, vks, uks, pre_spikes, rewiring, dropping)
        if prespike_recording:
            for j in range(model.Min):
                fprespike.write( str(t) + " " + str(j) + " " + str(true_pre_spikes[j]) + "\n" )
            fprespike.flush()
   
        if rewiring: #rewiring of connections with small spine sizes
            rwcnt = 0
            loctmp = []
            for i in range(model.Nin):
                if uks[i] < model.ukthreshold and uks[i] != 0.0:
                    if np.random.random() < model.rewiring_freq:
                        loctmp = lb.restricted_resampling(model, locs, i)
                                
                        locs[i][0] = loctmp[0]; locs[i][1] = loctmp[1] 
                        vks[i] = pick_a_vk(model,locs[i])
                        uks[i] = model.ukinit

                        dks[i] = model.dists[locs[i][0]] + locs[i][1]*model.dends[locs[i][0]].L
                        lb.rewire_synapse(model,locs[i],i)
                        rwcnt += 1
        if dropping: #elimination of inactive connections
            for i in range(model.Nin):
                if mean_pre_rates[i] < model.pre_rates_th and np.random.random() < model.rewiring_freq:
                    uks[i] = 0.0
        for i in range(model.Nin):
            if uks[i] != 0.0 and uks[i] < model.ukthreshold:
                uks[i] = model.uk_min
            
if __name__ == "__main__":
    #def main(args=None):
    #"""Main"""
    param = sys.argv
    Kin = int(param[1]) #total redundancy of synapses (Kin=5)
    gmax = float(param[2]) #standard conductance [muS] (gmax=0.0025)
    gI = float(param[3]) #inhibitory conductance [muS] (gI = 0.00075)
    uk_min = float(param[4]) #threshold of rewiring (ukth=0.001)
    sd_bias = float(param[5]) #bias in synaptic distribution (sd_bias=1.0)
    release_prob = float(param[6]) #Release probability of synaptic vesicle (release_prob = 1.0)
    ik = 0
    
    simul(Kin, gmax, gI, uk_min, sd_bias, release_prob, ik)
















