'''
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
with open("apd_dc.csv") as f:
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
V0 = 0.97  # the default visibility
dc0 = func(12, *popt) #150e3  # the default dark count reate
tau_c0 = 2  # coincidance window size in ns
r_pair0 = 1e6  # default rate of pair generation
tmp0 = 12
axcolor = 'lightgoldenrodyellow'

# global parameters

# Transmission rate
T = 1

# pairs generated in the sweetspot of the source
r_pair = 1e6  # per second

# Dark counts are the same for both Alice and Bob
dc = 150e3

# coincidence window
tau_c = 1  # seconds

# visibility
V = 0.97  # 97% of the paris are good!

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

    qt = c_QBER(r_pair, dc, tau_c * 1e-9, V, T)

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
T_db = [x * (-0.01) for x in reversed(range(0, 6001))]  # transmission factor in dB
# print T_list

print len(T_db), 'tdb'
T_list = [10 ** (x / 10) for x in T_db]  # transmission factor
print len(T_list), 'tlist'

QBER_list = [c_QBER(r_pair, dc, tau_c * 1e-9, V, x) for x in T_list]
# print QBER_list

# compute the private key rate
keyr_list = [c_private(r_pair, dc, tau_c * 1e-9, V, x) for x in T_list]

# print keyr_list
# the classicla link rate Bob to Alice
ccr_list = [c_ccr(r_pair, dc, tau_c * 1e-9, V, x) for x in T_list]

# print ccr_list
# def c_private(qt,rc,T):


'''
plotting code below
'''
# 10 
gs = gridspec.GridSpec(11, 1,
                       height_ratios=[12, 1, 12, 1, 12, 4, 1, 1, 1, 1,1],
                       width_ratios=[35, 1, 2, 1, 2, 1, 1, 1, 1, 1,1],
                       )

gs.update(left=0.30, right=0.95, wspace=0.05)

# the main figure object
fig = plt.figure(figsize=(12, 10))
fig.canvas.set_window_title('S^2QKD simulations')
axQBER = plt.subplot(gs[0])
axQBER.set_ylabel('QBER')
# axQBER.set_xlabel('Transmission factor T (dB)')
plt.axhline(y=0.11, xmin=0, xmax=1, hold=None)

lqb, = axQBER.plot(T_db, QBER_list, lw=2, color='red')
axQBER.axis([-60, 0, 0, 0.5])

# make these tick labels invisible
# plt.setp(axQBER.get_xticklabels(), visible=False)


axtau = plt.subplot(gs[5 + 2], axisbg=axcolor)
stau = Slider(axtau, 'Coincidence window (ns)', 0.5, 5, valinit=tau_c0)

axdc = plt.subplot(gs[4 + 2], axisbg=axcolor)
sdc = Slider(axdc, 'Dark count rate at each side', 0, 400e3, valinit=dc0)

axrpair = plt.subplot(gs[7 + 2], axisbg=axcolor)
srpair = Slider(axrpair, 'Entengled pair generation rate', 0, 2e6, valinit=r_pair0)

axtmp = plt.subplot(gs[10], axisbg=axcolor)
stmp = Slider(axtmp, 'Detector temperature (C) ' , -25, 25, valinit=tmp0)


axvis = plt.subplot(gs[6 + 2], axisbg=axcolor)
svis = Slider(axvis, 'Visibility', 0.78, 1, valinit=V0)

axccr = plt.subplot(gs[2], sharex=axQBER)
axccr.set_ylabel('Classical bit rate (Mbps)')
axccr.axis([-60, 0, 0, 6 * 8])

lcr, = axccr.plot(T_db, ccr_list, lw=2, color='green')

# axccr.set_xlabel('Transmission factor T (dB)')
# make these tick labels invisible
# plt.setp(axccr.get_xticklabels(), visible=False)


axkeyr = plt.subplot(gs[4], sharex=axQBER)
axkeyr.set_ylabel('Private key rate (kbps)')
axkeyr.set_xlabel('Transmission factor T (dB)')
axkeyr.axis([-60, 0, 0, 50 * 8])

lkr, = axkeyr.plot(T_db, keyr_list, lw=2, color='blue')


# the onclick update module for the slider
def update(val):
    global V
    global dc
    global tau_c
    global QBER_list
    global keyr_list
    global ccr_list
    global r_pair

    tau_c = stau.val
    V = svis.val
    dc = sdc.val
    r_pair = srpair.val

    QBER_list = [c_QBER(r_pair, dc, tau_c * 1e-9, V, x) for x in T_list]
    lqb.set_ydata(QBER_list)

    keyr_list = [c_private(r_pair, dc, tau_c * 1e-9, V, x) for x in T_list]
    lkr.set_ydata(keyr_list)

    ccr_list = [c_ccr(r_pair, dc, tau_c * 1e-9, V, x) for x in T_list]
    lcr.set_ydata(ccr_list)

    fig.canvas.draw_idle()




def update_temp(val):
    global dc
    dc = func(stmp.val, *popt)
    sdc.set_val(dc)
    
    
    
svis.on_changed(update)
sdc.on_changed(update)
stau.on_changed(update)
srpair.on_changed(update)
stmp.on_changed(update_temp)

resetax = plt.axes([0.85, 0.93, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

    
# the reset function for the reset button
def reset(event):
    svis.reset()
    sdc.reset()
    stau.reset()
    srpair.reset()
    stmp.reset()


button.on_clicked(reset)

plt.savefig("QBER.png")

plt.show()
