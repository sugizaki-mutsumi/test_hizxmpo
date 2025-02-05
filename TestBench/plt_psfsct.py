#!/usr/bin/env python

#import os
import matplotlib.pyplot as plt
import numpy as np

import uproot

if __name__=="__main__" :

    ### for HiZ
    imgz = 300
    
    #rootfname  = "test_gamma1kev.root"
    #rootfname  = "gamma1kev_1e4.root"
    rootfname  = "gamma_1000ev_x0y0_n100000.root"
    
    rootfile = uproot.open(rootfname)
    tevthits = rootfile["tEventHits"]
    
    vposx = tevthits["posX"].array()
    vposy = tevthits["posY"].array()
    vposz = tevthits["posZ"].array()
    vnhits = tevthits["NumHitData"].array()

    vsel_focal = np.fabs(vposz-imgz)<0.1
    vcut_infov = np.logical_and(np.fabs(vposx)<62.0, np.fabs(vposy)<62.0)
    vcut = np.logical_and(vsel_focal, vcut_infov)
    ncut = np.sum(vcut*1)
    #ncut = int((vcut*1).sum())

    sigma = 0.0
    #sigma = 0.371 ## 10 arcmin at 300 mm focal length
    #sigma = 0.185 ### 5 arcmin at 300 mm focal length
    vpsfx = vposx[vcut] + np.random.normal(0.0, sigma, ncut)
    vpsfy = vposy[vcut] + np.random.normal(0.0, sigma, ncut)
    vpsfnhits = vnhits[vcut]

    
    fig, ax = plt.subplots(figsize=(6,6))
    #ax.set_xlim(-5,5)
    #ax.set_ylim(-5,5)
    
    nhits_list = [3,5,7,9,11,13]
    color_list = ['k', 'r', 'g', 'b', 'm']
    #grphlist = []
    num = len(nhits_list)
    #for nhits in nhits_list :
    for idx in range(num) :
        nhits = nhits_list[idx]

        if idx<len(color_list) :
            color = color_list[idx]
        else :
            color = color_list[-1]
        #color = "C%d"%(idx)
        
        #vcuthit = np.logical_and(vcut, vnhits==nhits)
        #ax.plot(vposx[vcuthit],vposy[vcuthit], '.', ms=0.5, color=color)
        vcuthit = vpsfnhits==nhits
        ax.plot(vpsfx[vcuthit],vpsfy[vcuthit], '.', ms=0.5, color=color, alpha=0.5)
        #ax.plot(vpsfx[vcuthit],vpsfy[vcuthit], '.', ms=0.5, color='k', alpha=0.1)

    plt.ion()
    plt.show()
    
        
