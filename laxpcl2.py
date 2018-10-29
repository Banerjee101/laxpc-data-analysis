#############################################################################################################################################################
#This code was originally developed by Vikas Chand (TIFR) during his PhD Thesis work. 
#My work on this code was to improve on it as I see Fit
#If you want to edit/improve it, please leave this section intact
#also, for any suggestions, email ankushabanerjee@gmail.com
#############################################################################################################################################################
from __future__ import division
import os
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.io import fits
from pylab import *
from astropy.io.fits import update



def LAXE(E_start, E_stop, E_rsp1):
	for i in range(len(E_rsp1)):
		if ( E_rsp1[i] <= E_start ) and (E_rsp1[i+1] > E_start ):
			E_start_n = i
			print i
	for i in range(len(E_rsp1)):
		if (E_rsp1[i] <= E_stop ) and ( E_rsp1[i+1] > E_stop ):
			E_stop_n = i
			print i
	print "Start Energy is: ", E_rsp1[E_start_n]
	print "Stop Energy is: ", E_rsp1[E_stop_n]
	print "Start Channel is: ", E_start_n
	print "Stop Channel is: ", E_stop_n
	return E_start_n,E_stop_n, E_rsp1[E_start_n], E_rsp1[E_stop_n]




def pcu_chan_sel(pcu_no):
    if pcu_no==1:
        N_Chns = 512
    elif pcu_no==2:
        N_Chns = 256
    return N_Chns


pcu_no = 2					#PCU you need
time_bin = 1					#Time Bin in seconds
anod = 0					#Anode selection, refer the laxpc readme file

minimum_energy = 3.0				#energy range selecton. Need to be float
maximum_energy = 50.0

uld_bin = -1 					#change as required, ref the laxpc readme file for reference
Erth_occ = 1 
evt_flg = 2

if pcu_no == 1:
	rsp_file = ascii.read("lx10cshm07v1.0_t1.txt", header_start = None , data_start = 0)			#change this each time u use a new observation/response
	print rsp_file
elif pcu_no ==2:
	rsp_file = ascii.read("lx20cshp05v1.0_t1.txt", header_start = None , data_start = 0)
	print rsp_file
else:
	print "Check pcu_unit, code wants 1 or 2"

Ersp1 = rsp_file['col2']
Ersp2 = rsp_file['col3']



ch_min, ch_max, E_min, E_max = LAXE(minimum_energy, maximum_energy,  Ersp1)
print ch_min, ch_max, E_min, E_max

ascii.write([[ch_min], [ch_max],[E_min], [E_max]], "Selected_channels_for_this_run.txt", overwrite=True)

pcu_nos=pcu_chan_sel(pcu_no)


 
    
lx_cmd = """printf \""""+str(pcu_no)+""" """+str(time_bin)+""" """+str(anod)+"""\n""" \
         +str(ch_min)+""" """+str(ch_max)+"""\n"""+str(pcu_nos)+"""\n""" \
         +str(uld_bin)+""" """+str(Erth_occ)+""" """+ str(evt_flg)+"""\n\" | ./laxpcl1"""
print lx_cmd
os.system(lx_cmd)


rename_lc = "mv lxp"+str(pcu_no)+"level2.lcurv"+" lxp"+str(pcu_no)+str(anod)+"level2_"+str(time_bin)+"s_"+str(int(round(E_min)))+"_"+str(int(round(E_max)))+".lcurv"
os.system(rename_lc)

mv_results = "mv ./"+"lxp"+str(pcu_no)+str(anod)+"level2_"+str(time_bin)+"s_"+str(int(round(E_min)))+"_"+str(int(round(E_max)))+".lcurv "+"Products"
os.system(mv_results)
	








