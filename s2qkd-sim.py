'''
    This is a temporary branch created only for generateind a few plots.
    
    s2qkd-sim.py:   This interactive program estimates and visualizes Qubit error
                    rate, classical communication rate and private key rate as a
                    function of the Transmission rate and other parameters.

    input:          Sliders are used to to variate the following parameters
                    Dark count rate                     =   Dark Count on Alice or
                                                             Bob's side per second.
                    Coincidence window                  =   Coincidence window size
                                                            in nano seconds
                    Visibility                          =   The visibility of
                                                            entanglement. (this
                                                            parameter quantifies
                                                            the quality of 
                                                            entanglement)
                    Entangled pair generation rate      =   The number of entangled
                                                            pairs generated from
                                                            the source per 
                                                            second.   
                    Detector temperature (C)            =   Given in Celsius. The
                                                            detector temperature
                                                            drives the dark
                                                            count rate.  
                    reset                               =   On click the sliders
                                                            are reset to the
                                                            default values.
                   
    output:         outputs a QBER.png file containing the default plots.
                   
    Usages:         To use this program run:    python s2qkd-sim.py
   
    Prerequisites:  This program is tested for python2.7
                    You need the following libraries
                  
                    numpy
                    Scipy
                    matplotlib 
                   
                    To install them on Debian or Ubuntu run:
                        
                    sudo apt-get install python-numpy python-scipy python-matplotlib 
                    
                    Alternatively you can also use pip to install them. Please look
                    up the pip online documentation for pip instructions.
                        
    Copyright (C) 2016 Tanvirul Islam, National University
                         of Singapore <tanvirulbd@gmail.com>

    This source code is free software; you can redistribute it and/or
    modify it under the terms of the GNU Public License as published 
    by the Free Software Foundation; either version 3 of the License,
    or (at your option) any later version.

    This source code is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    Please refer to the GNU Public License for more details.

    You should have received a copy of the GNU Public License along with
    this source code; if not, see: <https://www.gnu.org/licenses/gpl.html>

'''

import math
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy.optimize import curve_fit



#load temperature vs dark count data

TM = []
DC = []
with open("new_apd.csv") as f:
    for line in f:
        tmp, dc = line.split()
        TM.append(float(tmp)) 
        DC.append(float(dc)*4) # there are 4 APDs

#curve fit

def func(x, a, b, c):
    return a * np.exp(b * x) + c

#x = np.linspace(0,4,50)
#y = func(x, 2.5, 1.3, 0.5)
x = np.array(TM)
yn = np.array(DC)

print len(TM)
print len(DC)
#exit()
#yn = y + 20*np.random.normal(size=len(x))

popt, pcov = curve_fit(func, x, yn)




# default values
tmp0 = 10
V0 = 0.97  # the default visibility
dc0 = func(tmp0, *popt) #150e3  # the default dark count reate
tau_c0 = 1  # coincidance window size in ns
r_pair0 = 1e6  # default rate of pair generation

axcolor = 'lightgoldenrodyellow'

# global parameters

# Transmission rate
T = 1

# pairs generated in the sweetspot of the source
r_pair = r_pair0  # per second

# Dark counts are the same for both Alice and Bob
dc = dc0

# coincidence window
tau_c = tau_c0  # nano seconds

# visibility
V = V0  # 97% of the paris are good!

'''
    function:   c_si() computes the detector counts of the sender    
    inputs: r_pair = the Entangled Pair generation rate
            dc  = Dark count reate at the senders side. Cumulative of all the 4 detectors.
    output: si = the singles count on the senders side
'''


def c_si(r_pair, dc):  # singles on idler Alice
    return r_pair + 4 * r_pair + dc  # 20% of the singles are coincidence
    # will change depending on the window
    # size and other parameters.


'''
    function:   c_ss() computes the detector counts of the receiver    
    inputs: r_pair = the Entangled Pair generation rate
            dc  = Dark count reate at the receiver's side. Cumulative of all the 4 detectors.
    output: si = the singles count on the receiver's side
'''


def c_ss(r_pair, T, dc):
    return r_pair * T + dc  # Singles on Bob is dark count + the portion of the pairs successfully transmitted.


# signal coincidence rate after loss
def c_rc(r_pair, T):
    return r_pair * T


# compute the accidental coincedence rate per second.
def c_ra(si, ss, tau_c):
    return si * ss * tau_c  # the transmission ratio is considered inside the ss


# intrinsic Qubit Error Rate
def c_qi(V):
    return (1 - V) / 2.0


# signal rate is 1/2 of the coincidence rate.
def c_rsig(rc):
    return rc * 0.5


# compute overall QBER

''' function: c_QBER computes the overall qubit error rate
    input:  r_pair = pair rate /s
            dc = Dark count /s
            tau_c = coincidence window s
            V =  visibility 
            T = Transmission factor
    output: QBER =  the overall quantum bit error rate.
'''


def c_QBER(r_pair, dc, tau_c, V, T):
    ss = c_ss(r_pair, T, dc)
    si = c_si(r_pair, dc)
    ra = c_ra(si, ss, tau_c)
    rc = c_rc(r_pair, T)
    qi = c_qi(V)
    rsig = c_rsig(rc)

    QBER = (1 / (rsig + ra)) * (qi * rsig + 0.5 * ra)
    return QBER


qt = c_QBER(r_pair, dc, tau_c, V, T)  # for test

'''
    function:   c_private() computes the number of private keys generated per second.
    input:  r_pair = pair rate /s
            dc = Dark count /s
            tau_c = coincidence window s
            V =  visibility 
            T = Transmission factor
            
    output:
            r_private = private key rate in kbPS (kilo bits per seconds)
'''


def c_private(r_pair, dc, tau_c, V, T):
    rc = c_rc(r_pair, T)  # actual coincidance
    rsig = c_rsig(rc)
    # ra = c_ra(si,ss,tau_c)

    qt = c_QBER(r_pair, dc, tau_c, V, T)

    # number of raw key bits lost to error correction 
    rlost_err = rsig * qt * math.log(1 / qt, 2)

    # print rlost_err

    # corrected raw keys
    rcorr = rsig - rlost_err

    # print rcorr

    # after privacy amplification the corrected key is shrinked by half

    r_private = rcorr / 2

    return r_private / (1000)  # bits to kilo Bytes


print 'key'
# print c_private(qt,r_pair,T)

print 'QBER', qt

'''
    compute the classical link rate in mega bits per secont Mbps (10**6 bits = 1 mega bits)
    function:   c_ccr() computes the classical communication required to transmit the tags from Bob to Alice (receiver to sender)
    input:  r_pair = pair rate /s
            dc = Dark count /s
            tau_c = coincidence window s
            V =  visibility 
            T = Transmission factor

    output: ckBps = Number of classical bits exchanged to transfer the compressed time tags. In MBps (mega bytes per second)
'''


def c_ccr(r_pair, dc, tau_c, V, T):
    raw_tag_size = 32.0
    # T = 0.01 #test only
    ss = c_ss(r_pair, T, dc)  # Bob's count rate

    # time per event
    tpe = 1 / ss

    # if the clock runs at 1 tick each nanosecond then
    # ticks per event
    ticks_pe = tpe / 1e-9

    # MSB position per event tag
    pmsb = math.ceil(math.log(ticks_pe + 1, 2))

    # number of 1 in the binary tag on average
    num_of_1 = pmsb / 2  # assume on avarage half of the bits after msb are 1
    # print 'pmsb', pmsb
    p1 = num_of_1 / raw_tag_size  # assuming 32 bit tags
    p0 = (raw_tag_size - num_of_1) / raw_tag_size  # if not 1 it is 0

    # shanon ent per bit
    h = -p1 * math.log(p1, 2) - p0 * math.log(p0, 2)
    # print 'shanon', h
    # lassical bits needed for tagging ss number of events
    total_tag_bits = ss * raw_tag_size

    # total compressed tag bits
    c_t_bits = total_tag_bits * h
    # print total_tag_bits/(8*1024)

    # print c_t_bits*1.2/(8*1024)
    ckBps = c_t_bits * 1.2 / (10 ** 6)

    return ckBps  # 1.2 is the implementation ineffieiency of the compression algorithm


print 'ccr', c_ccr(r_pair, dc, tau_c, V, T)  # test

# T_list = [x*0.0001 for x in range(0,10001)]
T_db = [x * (-0.01) for x in reversed(range(1500, 6001))]  # transmission factor in dB
# print T_list

print len(T_db), 'tdb'
T_list = [10 ** (x / 10) for x in T_db]  # transmission factor
print len(T_list), 'tlist'

# QBER_list = [c_QBER(r_pair, dc, tau_c * 1e-9, V, x) for x in T_list]
# # print QBER_list

# # compute the private key rate
# keyr_list = [c_private(r_pair, dc, tau_c * 1e-9, V, x) for x in T_list]

# # print keyr_list
# # the classicla link rate Bob to Alice
# ccr_list = [c_ccr(r_pair, dc, tau_c * 1e-9, V, x) for x in T_list]

# # print ccr_list
# # def c_private(qt,rc,T):


'''
plotting code below
'''
# 10 
gs = gridspec.GridSpec(11, 1,
                       height_ratios=[12, 1, 12, 1, 12, 4, 1, 1, 1, 1,1],
                       width_ratios=[35, 1, 2, 1, 2, 1, 1, 1, 1, 1,1],
                       )

gs.update(left=0.30, right=0.95, wspace=0.05)

'''The figure drawing'''

fig = plt.figure(figsize=(12, 6))
fig.canvas.set_window_title('S^2QKD TEst')

#generate all the plot data for the following temperature values
#tmptmp = [-20,-15,-10,-5,0,5,10,15,20]
tmptmp = [-10,-5,0,5,10,15,20]
#tmptmp = [-20,-10,0]

qball = {} #dictionary

for t in tmptmp:
        dc = func(t, *popt)
        qball[t] = [c_QBER(r_pair, dc, tau_c * 1e-9, V, x) for x in T_list]
        #qball[t] = [c_QBER(r_pair, dc, tau_c * 1e-9, V, x) for x in T_list]

#plt.plot(x,y,label="old noise data");
#plt.plot(x,z,label="old noise data");
#plt.plot(T_db, QBER_list, lw=2, color='red')

plt.axhline(y=0.11, xmin=0, xmax=1, hold=None)
pos_T_db = [-1*vx for vx in T_db]

#for t in reversed(tmptmp):
        #plt.plot(pos_T_db, qball[t], lw=2,label="APD temp. = "+str(t)+"$ ^\circ$C")
#        plt.plot(pos_T_db, qball[t], lw=2,label=str(t)+"$ ^\circ$C")
        #plt.plot(pos_T_db, qball[t], lw=2,label=str(t+273.15)+"K")

plt.plot(pos_T_db, qball[20], lw=2,label=str(20)+"$ ^\circ$C", linestyle='-.')
plt.plot(pos_T_db, qball[15], lw=2,label=str(15)+"$ ^\circ$C", linestyle='--')
plt.plot(pos_T_db, qball[10], lw=2,label=str(10)+"$ ^\circ$C", dashes=[3,2,3,4])
plt.plot(pos_T_db, qball[5], lw=2,label=str(5)+"$ ^\circ$C", linestyle=':')
plt.plot(pos_T_db, qball[0], lw=2,label=str(0)+"$ ^\circ$C", dashes=[1,2,3,4])
plt.plot(pos_T_db, qball[-5], lw=2,label=str(-5)+"$ ^\circ$C", dashes=[3,4])
plt.plot(pos_T_db, qball[-10], lw=2,label=str(-10)+"$ ^\circ$C", linestyle='-')
#inserting degree celsius symbol
#http://stackoverflow.com/questions/19926246/inserting-a-degree-symbol-into-python-plot


plt.xlabel("Optical losses (dB)")
plt.ylabel("QBER")
plt.title("QBER vs. Optical losses at different APD temperatures")

plt.legend(loc=0)
plt.savefig("QBER-TempPlot-bw-friendly.pdf")

plt.show()
