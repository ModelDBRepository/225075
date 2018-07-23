# ----------------------------------------------------------
# Library of cell classes and functions
#
# Naoki Hiratani, N.Hiratani@gmail.com
#
# Based on the model by Dr. Tiago Branco
# ----------------------------------------------------------

import numpy as np
import neuron

from neuron import h
from neuron import load_mechanisms
from neuron import gui

from math import *
import random as pyrnd
import matplotlib.pyplot as plt
import scipy.special as scisp

load_mechanisms('/directory_where_mod_files_have_been_compiled')
h('objref nil')

# ----------------------------------------------------------
# MODELS
class L23(object):
    def __init__(self):
        h('xopen("./L23.hoc")')
        props(self)
        self._geom()
        self._topol()
        self._biophys()
        self._dist_from_soma()

    def _geom(self):
        self.axon = h.Section()
        self.axon.L = 300
        self.axon.diam = 1

    def _topol(self):            
        self.soma = h.soma
        self.axon.connect(self.soma,1,0)
        #self.axon = h.axon
        self.dends = [] #list of dendritic section
        for sec in h.allsec():
            self.dends.append(sec)
        self.dends.remove(self.soma)   # Remove soma from the list
        self.dends.remove(self.axon)   # and the Axon

        dnames = ['dend1_1','dend1_11','dend1_111','dend1_1111','dend1_1112','dend1_112','dend1_1121','dend1_1122','dend1_12','dend1_121','dend1_1211','dend1_1212','dend1_122','dend1_1221','dend1_1222','dend1_12221','dend1_12222','dend2_1','dend2_11','dend2_111','dend2_1111','dend2_1112','dend2_112','dend2_1121','dend2_1122','dend2_12','dend2_121','dend2_1211','dend2_12111','dend2_121111','dend2_121112','dend2_12112','dend2_121121','dend2_1211211','dend2_1211212','dend2_12112121','dend2_12112122','dend2_121122','dend2_1212','dend2_12121','dend2_121211','dend2_121212','dend2_12122','dend2_122','dend2_1221','dend2_12211','dend2_12212','dend2_1222','dend2_12221','dend2_12222','dend3_1','dend3_11','dend3_111','dend3_1111','dend3_1112','dend3_11121','dend3_11122','dend3_112','dend3_1121','dend3_1122','dend3_12','dend3_121','dend3_1211','dend3_1212','dend3_12121','dend3_12122','dend3_121221','dend3_121222','dend3_1212221','dend3_1212222','dend3_12122221','dend3_12122222','dend3_121222221','dend3_1212222211','dend3_12122222111','dend3_12122222112','dend3_1212222212','dend3_121222222','dend3_1212222221','dend3_1212222222','dend3_12122222221','dend3_12122222222','dend3_122','dend3_1221','dend3_12211','dend3_12212','dend3_1222','dend4_1','dend4_11','dend4_111','dend4_1111','dend4_1112','dend4_11121','dend4_11122','dend4_112','dend4_1121','dend4_1122','dend4_12','dend4_121','dend4_122','dend4_1221','dend4_1222','dend4_12221','dend4_12222']
    
    def _biophys(self):
        for sec in h.allsec():
            sec.cm = self.CM
            sec.insert('pas')
            sec.e_pas = self.E_PAS
            sec.g_pas = 1.0/self.RM
            sec.Ra = self.RA
            sec.nseg = int(ceil(sec.L/5.0)) #5um segmentation

    def _dist_from_soma(self):
        self.dists = [] #distance from the soma to dendrites
        h('soma distance(0.0,0.5)')
        lidx = 0
        for sec in h.allsec():
            if lidx < len(self.dends):
                h('''proc calc_distance(){
                    disttmp = distance(1.0,0.0)}''')
                h.calc_distance()
                self.dists.append(h.disttmp)
            lidx += 1
            
        self.branchdists = [] #dendritic distance between two branches
        l1idx = 0
        for sec1 in h.allsec():
            if l1idx < len(self.dends):
                self.branchdists.append([])
                h('distance()') 
                l2idx = 0
                for sec2 in h.allsec():
                    if l2idx < len(self.dends):
                        h('''proc calc_distance(){
                            disttmp = distance(1.0,0.0)}''')
                        h.calc_distance()
                        self.branchdists[l1idx].append( h.disttmp )
                    l2idx += 1
            l1idx += 1

# ----------------------------------------------------------
# INSTRUMENTATION FUNCTIONS
def props(model):

    # Passive properties
    model.CM = 1.0
    model.RM = 7000.0
    model.RA = 100
    model.E_PAS = -75
    model.CELSIUS = 35

    # Active properties
    model.Ek = -90
    model.Ena = 60
    model.Eca = 140
    
    model.gna_axon = 1000
    model.gkv_axon = 100
    
    model.gna_soma = 1000
    model.gkv_soma = 100 
    model.gkm_soma = 2.2 
    model.gkca_soma = 3 
    model.gca_soma = 0.5 
    model.git_soma = 0.0003 
    
    model.gna_dend = 27#80
    model.gna_dend_hotSpot = 600
    model.gkv_dend = 1#3
    model.gkm_dend = 0.3#1
    model.gkca_dend = 3
    model.gca_dend = 0.5
    model.git_dend = 0.00015
    model.gh_dend = 0

def init_active(model, axon=False, soma=False, dend=True, dendNa=False,
                dendCa=False):
    if axon:
        model.axon.insert('na'); model.axon.gbar_na = model.gna_axon
        model.axon.insert('kv'); model.axon.gbar_kv = model.gkv_axon
        model.axon.ena = model.Ena
        model.axon.ek = model.Ek

    if soma:
        model.soma.insert('na'); model.soma.gbar_na = model.gna_soma
        model.soma.insert('kv'); model.soma.gbar_kv = model.gkv_soma
        model.soma.insert('km'); model.soma.gbar_km = model.gkm_soma
        model.soma.insert('kca'); model.soma.gbar_kca = model.gkca_soma
        model.soma.insert('ca'); model.soma.gbar_ca = model.gca_soma
        model.soma.insert('it'); model.soma.gbar_it = model.git_soma
        #model.soma.insert('cad');
        model.soma.ena = model.Ena
        model.soma.ek = model.Ek
        model.soma.eca = model.Eca

    if dend:
        for d in model.dends:
            d.insert('na'); d.gbar_na = model.gna_dend*dendNa
            d.insert('kv'); d.gbar_kv = model.gkv_dend
            d.insert('km'); d.gbar_km = model.gkm_dend
            d.insert('kca'); d.gbar_kca = model.gkca_dend
            d.insert('ca'); d.gbar_ca = model.gca_dend*dendCa
            d.insert('it'); d.gbar_it = model.git_dend*dendCa
            #d.insert('cad')
            d.ena = model.Ena
            d.ek = model.Ek
            d.eca = model.Eca

def init_params(model, Kin, gmax, gI, uk_min, sd_bias, release_prob):
    #simulation control
    model.trials = 1001# #trial number

    #synapses
    model.Min = 200 #number of input neurons
    model.Kin = Kin #redanduncy in synaptic connections
    model.Nin = model.Kin*model.Min #number of synaptic inputs
    
    model.gmax = gmax#0.0015 #standard conductance [muS]
    model.sd_bias = sd_bias #bias in synaptic distribution    

    model.uk_min = uk_min #lower boundary condition of uks
    model.uk_max = 1.0 #upper boundary condition of uks

    #rewiring parameters
    model.ukthreshold = model.uk_min #spine cutoff threshold
    model.ukinit = 1.0/float(Kin) #model.uk_min #size of a newly created spine

    model.taumr = 10.0 #the time constant for mean firing rate
    model.rewiring_freq = 0.2 #relative frequency of rewiring
    model.pre_rates_th = 0.05 #threshold for the presynaptic firing rate in supervised trials
    
    #input structure
    model.spn_rate = 1.5*pi #spontaneous firing rate
    model.epsilon = 0.01 #the minimum relative rate
    model.min_rate = model.epsilon*model.spn_rate
    model.release_prob = release_prob# presynaptic release probability

    model.p_theta = 0.0 #preferred orientation
    model.n_theta = model.p_theta + pi/2.0 #non-preferred orientation

    model.rfdist_zero = 1.0 # standard distance between somatic and synaptic receptive fields
    model.rfdist_min = 0.01 # minimum distance between somatic and synaptic receptive fields
    model.rfdist_max = 3.0 # maximum distance between somatic and synaptic receptive fields
    
    model.kappa_zero = 2.0 # orientation selectivity
    model.kappa_phi = 4.0 # association field selectivity
    model.krfdist_const = exp(model.kappa_phi)/100.0 #constant for rfdist in calculation of k_i (to avoid numerical instability)

    #inhibitory synapses
    model.inhN = 200 #Number of inhibitory synapses
    model.gI = gI #inhibitory conductance
    model.Inh_thr = -90.0 #mV

    #stimulation protocol
    model.sstart = 20.0 #starting timing of the stimulation (from the initiation of the simulation)
    model.sduration = 20.0 #stimulation duration[ms]
    model.snoise = 0.0 #noise in spike train

def plot_morphology(model,locs):
    h('access soma')
    h('objref sh')
    h('sh = new PlotShape(0)')
    h('sh.size(-300,300,-299.522,299.522)')

    h('''proc plot_synapse(){
            xdtmp = x3d($1)
            ydtmp = y3d($1)
            ctmp = 7
            sh.mark(xdtmp,ydtmp,"o",6,ctmp,4)
        }''')
    lidx = 0
    for sec in h.allsec():
        if lidx < len(model.dends):
            for i in range(len(locs)):
                loc = locs[i]
                if loc[0] == lidx:
                    h('n3dtmp = n3d()')
                    ndtmp = int(floor(h.n3dtmp*loc[1]))
                    h.plot_synapse(ndtmp,model.preidx[i])

        lidx += 1
    h('sh.view(-300, -400.522, 600, 900.043, 265, 450, 200.64, 400.32)')

def add_Estims(model,locs,sstart=20.0,sinterval=20,snumber=1,snoise=0):
    model.Estimlist = []
    lidx = 0
    for loc in locs:
        Estim = h.NetStim()
        Estim.interval = sinterval
        Estim.number = snumber
        Estim.start = sstart
        Estim.noise = snoise
        model.Estimlist.append(Estim)
        lidx += 1

def add_Istims(model,inh_locs,sstart=20.0,sinterval=20,snumber=1,snoise=0):
    model.Istimlist = []
        
    lidx = 0
    for loc in inh_locs:
        Istim = h.NetStim()
        Istim.interval = sinterval
        Istim.number = snumber
        Istim.start = sstart
        Istim.noise = snoise
        model.Istimlist.append(Istim)
        lidx += 1

#select dendritic locations for calculation of the unit EPSPs 
def calc_locs_kappa(model):
    dL = 5.0
    locs = []
    lidx = 0
    nidx = 0
    for sec in model.dends:
        nsec = int(ceil(sec.L/dL))
        for sidx in range( nsec ):
            locs.append([])
            locs[nidx].append( lidx )
            locs[nidx].append( (sidx+0.5)/float(nsec) )
            if locs[nidx][1] < 0.0 or locs[nidx][1] > 1.0:
                print lidx,nidx,locs[nidx][0],locs[nidx][1]
            nidx += 1
        lidx += 1
    #print '#kappa_locs',nidx
    return locs

#Uniform selection of dendritic locations of the N synaptic contacts
def calc_locs_uniform(model):
    locs = []
    Ltot = 0.0
    for sec in model.dends:
        Ltot += sec.L
    dL = Ltot/(model.Nin+1.0)
    Ltmp = 0.0
    nidx = 0
    lidx = 0
    for sec in model.dends:
        Ltmp += sec.L
        while nidx*dL < Ltmp and Ltmp <= (nidx+1)*Ltmp:
            if nidx > 0 and nidx <= model.Nin:
                locs.append([])
                locs[nidx-1].append(lidx)
                locs[nidx-1].append( 1.0 - (Ltmp-nidx*dL)/(sec.L) )
            nidx += 1
        lidx += 1
    #print Ltot, dL, len(locs)
    return locs

#Random selection of dendritic locations of the N synaptic contacts from a Beta distribution characterized by (bsd, 2-bsd)
def calc_locs_random(model):
    locs = []
    Lcums = []
    Lcums.append(0.0)
    for sec in model.dends:
        Lcums.append(Lcums[-1] + sec.L)
    Ltot = Lcums[-1]
    dmax = 500.0#max(model.dists)

    bsda = model.sd_bias; bsdb = 2.0 - model.sd_bias #parameters for Beta distribution
    Zbeta = 10.0*scisp.beta(bsda,bsdb)
    #print dmax, bsda, bsdb, Zbeta

    nidx = 0
    while(nidx < model.Nin):
        Ltmp = Ltot*pyrnd.random()        
        for j in range(len(model.dends)):
            if Lcums[j] <= Ltmp and Ltmp < Lcums[j+1]:
                rdtmp = (model.dists[j] + (Ltmp-Lcums[j]))/dmax
                if pyrnd.random() < pow(rdtmp,bsda-1.0)*pow(1.0-rdtmp,bsdb-1.0)/Zbeta:
                    locs.append([])
                    locs[nidx].append(j)
                    locs[nidx].append( (Ltmp-Lcums[j])/(Lcums[j+1]-Lcums[j]) )
                    nidx += 1

    #print Ltot, len(locs)
    return locs 

#Uniform selection of dendritic locations of the inhN synaptic contacts
def calc_inh_locs_uniform(model):
    inh_locs = []
    Ltot = 0.0
    for sec in model.dends:
        Ltot += sec.L
    dL = Ltot/(model.inhN+1.0)
    Ltmp = 0.0
    nidx = 0
    lidx = 0
    for sec in model.dends:
        Ltmp += sec.L
        while nidx*dL < Ltmp and Ltmp <= (nidx+1)*Ltmp:
            if nidx > 0 and nidx <= model.inhN:
                inh_locs.append([])
                inh_locs[nidx-1].append(lidx)
                inh_locs[nidx-1].append( 1.0 - (Ltmp-nidx*dL)/(sec.L) )
            nidx += 1
        lidx += 1
    print Ltot, dL, len(inh_locs)
    return inh_locs

#restricted resampling: the posisiton of the new synapse is restricted to the set of dendritic branches where original connections were made  
def restricted_resampling(model, locs, i): 
    loctmp = []
    
    Lcums = []
    Lcums.append(0.0)
    for j in range(model.Nin):
        if model.preidx[j] == model.preidx[i]: 
            Lcums.append(Lcums[-1] + model.dends[locs[j][0]].L)
    Ltot = Lcums[-1]
    
    Ltmp = Ltot*pyrnd.random()
    lidx = 0
    for j in range(model.Nin):
        if model.preidx[j] == model.preidx[i]: 
            if Lcums[lidx] <= Ltmp and Ltmp < Lcums[lidx+1]:
                loctmp.append( locs[j][0])
                loctmp.append( (Ltmp-Lcums[lidx])/(Lcums[lidx+1]-Lcums[lidx]) )
            lidx += 1
    return loctmp

def calc_dist_bt_syns(model,locs):#calculate the distance between synapses
    disttmps = []
    for i1 in range(model.Nin):
        jidx = model.preidx[i1]
        d1 = locs[i1][0]
        for i2 in range(i1+1,model.Nin):
            if model.preidx[i2] == jidx:
                d2 = locs[i2][0]
                disttmp = model.branchdists[d1][d2] + model.dends[d1].L*locs[i1][1] + model.dends[d2].L*locs[i2][1] 
                disttmps.append( disttmp )
    return disttmps
    
#presynaptic index allocation
def allocate_syns(model):
    rndidx = range(model.Nin)
    pyrnd.shuffle(rndidx)
    model.preidx = []
    for i in range(model.Nin):
        model.preidx.append(0)
    for iidx in range(model.Nin):
        model.preidx[rndidx[iidx]] = iidx/model.Kin

#add AMPA synapses
def add_AMPAsyns(model, locs=[[0, 0.5]], gmax=0.5, tau1=0.5, tau2=2.5):
    model.AMPAlist = []
    model.ncAMPAlist = []

    for lidx in range(len(locs)):
        loc = locs[lidx]
        AMPA = h.Exp2Syn(float(loc[1]), sec=model.dends[int(loc[0])]) 
        AMPA.tau1 = tau1
        AMPA.tau2 = tau2
        #NC = h.NetCon(h.nil, AMPA, 0, 0, gmax)
        NC = h.NetCon(model.Estimlist[lidx], AMPA, 0, 0, gmax)
        model.AMPAlist.append(AMPA)
        model.ncAMPAlist.append(NC)

def add_GABAsyns(model, inh_locs=[[0, 0.5]], gI=0.5, tau1=0.5, tau2=2.5):
    model.GABAlist = []
    model.ncGABAlist = []

    for lidx in range(len(inh_locs)):
        inh_loc = inh_locs[lidx]
        GABA = h.Exp2Syn(float(inh_loc[1]), sec=model.dends[int(inh_loc[0])]) 
        GABA.tau1 = tau1
        GABA.tau2 = tau2
        GABA.e = model.Inh_thr
        NC = h.NetCon(model.Istimlist[lidx], GABA, 0, 0, model.gI)
        model.GABAlist.append(GABA)
        model.ncGABAlist.append(NC)

def rewire_synapse(model,loc,lidx,tau1=0.5,tau2=2.5):
    AMPAtmp = h.Exp2Syn(float(loc[1]), sec=model.dends[int(loc[0])]) 
    AMPAtmp.tau1 = tau1
    AMPAtmp.tau2 = tau2
    NCtmp = h.NetCon(model.Estimlist[lidx], AMPAtmp, 0, 0, model.gmax)
    model.AMPAlist[lidx] = AMPAtmp
    model.ncAMPAlist[lidx] = NCtmp

def generate_pre_tuning(model):
    model.phis = [] #polar direction of RF
    model.thetas = [] #preferred orientation
    model.rfdists = [] #Receptive field distance
    rfdist_zero = model.rfdist_zero
    rfdist_min = model.rfdist_min
    rfdist_max = model.rfdist_max
    for j in range(model.Min):
        model.phis.append( 2.0*pi*np.random.random() )
        model.thetas.append( pi*np.random.random() )
        model.rfdists.append( rfdist_min + (rfdist_max-rfdist_min)*np.random.random() )

def generate_log_rates(model, thetao):
    rfdist_zero = model.rfdist_zero 
    krfdist_const = model.krfdist_const
    
    kappa_zero = model.kappa_zero
    kappa_phi = model.kappa_phi
    
    min_rate = model.min_rate
    spn_rate = model.spn_rate
    
    rates = []; log_rates = []
    for j in range(model.Min):
        phi = model.phis[j]; rfdist = model.rfdists[j]; theta = model.thetas[j]
        kappa_r = rfdist_zero*exp(kappa_phi*cos(2.0*(phi-thetao)))/(rfdist + krfdist_const)
        ktmp = sqrt( kappa_zero*kappa_zero + kappa_r*kappa_r + 2*kappa_zero*kappa_r*cos( 2.0*(theta-thetao) ) )
        rates.append( spn_rate*exp(-rfdist/rfdist_zero)*scisp.i0(ktmp)/(2*pi*scisp.i0(kappa_zero)*scisp.i0(kappa_r)) )
        if rates[-1] > min_rate:
            log_rates.append( log(rates[-1]/min_rate) )
        else:
            log_rates.append( 0.0 )
    return rates, log_rates

def generate_input_rates(model):
    # spontaneous firing rates
    model.s_rates = [] 
    for j in range(model.Min):
        model.s_rates.append( model.spn_rate )

    # firing rate for the preferred orientation
    p_rates, log_p_rates = generate_log_rates(model, model.p_theta)
    model.p_rates = p_rates
    model.log_p_rates = log_p_rates

    # firing rate for the non-preferred orientation
    n_rates, log_n_rates = generate_log_rates(model, model.n_theta)
    model.n_rates = n_rates
    model.log_n_rates = log_n_rates
    #print "p_rate, n_rate: ", np.average(p_rates), np.average(n_rates)
    
def generate_prespikes(model, pre_rates):
    release_prob = model.release_prob
    Min = model.Min; Nin = model.Nin
    true_pre_spikes = []; pre_spikes = []
    for j in range(Min):
        #Poisson pre-spikes
        true_pre_spikes.append( np.random.poisson( pre_rates[j] ) )
    for i in range(Nin):
        jidx = model.preidx[i]
        pre_spikes.append( np.random.binomial( true_pre_spikes[jidx], release_prob ) )

    return true_pre_spikes, pre_spikes

def generate_inputs(model, uks, pre_spikes):
    release_prob = model.release_prob
    for i in range(model.Nin):
        if pre_spikes[i] > 0:
            model.Estimlist[i].interval = model.sduration/float(pre_spikes[i])
        else:
            model.Estimlist[i].interval = model.sduration            
        model.Estimlist[i].number = pre_spikes[i]
        model.Estimlist[i].start = model.sstart + (model.Estimlist[i].interval)*np.random.random()
        model.Estimlist[i].noise = model.snoise

        weight_tmp = model.gmax*uks[i]/release_prob
        #weight_tmp = model.gmax*uks[i]
        NC = h.NetCon(model.Estimlist[i], model.AMPAlist[i], 0, 0, weight_tmp)
        model.ncAMPAlist[i] = NC

def generate_inh_inputs(model, true_pre_spikes):
    inhN = model.inhN
    Nin = model.Nin; Min = model.Min

    Erate = 0.0
    for i in range(Nin):
        Erate += true_pre_spikes[ model.preidx[i] ]/float(Nin)
    Irate = ( float(Nin)/float(inhN) )*Erate
    Ispikes = np.random.poisson( Irate, (inhN) )
    
    #print Ispikes
    for i in range(inhN):
        if Ispikes[i] > 0:
            model.Istimlist[i].interval = model.sduration/float(Ispikes[i])
        else:
            model.Istimlist[i].interval = model.sduration            
        model.Istimlist[i].number = Ispikes[i]
        model.Istimlist[i].start = model.sstart + (model.Istimlist[i].interval)*np.random.random()
        model.Istimlist[i].noise = model.snoise
        
        NC = h.NetCon(model.Istimlist[i], model.GABAlist[i], 0, 0, model.gI)
        model.ncGABAlist[i] = NC
    

# ----------------------------------------------------------
# SIMULATION RUN
def simulate(model):
    trec, vrec = h.Vector(), h.Vector()
    trec.record(h._ref_t)
    vrec.record(model.soma(0.5)._ref_v)

    h.celsius = model.CELSIUS
    #h.FInitializeHandler(1, initSpikes2)
    h.finitialize(model.E_PAS)
    neuron.run(model.tstop)

    return np.array(trec), np.array(vrec)

