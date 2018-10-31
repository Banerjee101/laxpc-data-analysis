###################################################################################################################################################################################################################################################################################################################################
#This code was originally developed by Vikas Chand (TIFR) during his PhD Thesis work. 
#My work on this code was to improve on it as I see Fit
#If you want to edit/improve it, please leave this section intact
#also, for any suggestions, email ankushabanerjee@gmail.com
###################################################################################################################################################################################################################################################################################################################################
from __future__ import division
import os
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.io import fits
from pylab import *
from astropy.io.fits import update

###################################################################################################################################################################################################################################################################################################################################

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

def pcu_chan_bin(pcu_no, ch_min, ch_max):
	if pcu_no==1:
        	ch_min_binned,ch_max_binned = ch_min*2, ch_max*2
		N_Chns = 512
    	elif pcu_no==2:
        	ch_min_binned,ch_max_binned = ch_min*4, ch_max*4
		N_Chns = 256
    	return ch_min_binned,ch_max_binned,N_Chns

###################################################################################################################################################################################################################################################################################################################################

pcu_no = 1					#PCU you need
time_bin = 16					#Time Bin in seconds
anod = 1					#Anode selection, refer the laxpc readme file (ian)

date_obs_yy = 17				#date of observation, change this as required
date_obs_mm = 4
date_obs_dd = 8

minimum_energy = 3.0				#energy range selecton. Need to be float
maximum_energy = 21.0				#For HID and/or CCD remember to use 1:2 2:4 ratio for the energies of soft and hard colors (3-5 5-9 9-13 13-21 3-9 9-21 3-21)
						#Also, when plotting the HID, you need to take the intensity of the LC with the energy ranges spanning the extreme ends of the HID energies that you are considering
						#meaning, when plotting HID of 3-5 and 5-9 you take the intensity axis from 3-9

uld_bin = -1 					#change as required, ref the laxpc readme file for reference
Erth_occ = 1 
evt_flg = 2

directory_source_loc = "/home/ankush/Desktop/Workspace_project/1A_2466-588/LAXPC/outputdata/laxpcsoftv2.5_6Aug2018/"		#Change this as needed. This is where your data would be generated and later stored
directory_dest_loc = "./Productions/"

rmf_list_lxp1 = "lx10cshm07v1.0_t1.txt"					#change this each time u use a new observation/response. You would have to open the respons file and save it as an ascii table using fv (heasoft)
rmf_list_lxp2 = "lx20cshp05v1.0_t1.txt"

###################################################################################################################################################################################################################################################################################################################################

if pcu_no == 1:
	rsp_file = ascii.read(rmf_list_lxp1, header_start = None , data_start = 0)
	print rsp_file
elif pcu_no ==2:
	rsp_file = ascii.read(rmf_list_lxp2, header_start = None , data_start = 0)
	print rsp_file
else:
	print "Check pcu_unit, code wants 1 or 2"

Ersp1 = rsp_file['col2']
Ersp2 = rsp_file['col3']

ch_min, ch_max, E_min, E_max = LAXE(minimum_energy, maximum_energy,  Ersp1)
print ch_min, ch_max, E_min, E_max

ch_min_binned,ch_max_binned,ch_nos = pcu_chan_bin(pcu_no, ch_min, ch_max)

ascii.write([[ch_min_binned], [ch_max_binned],[E_min], [E_max]], "Selected_channels_for_this_run.txt", overwrite=True)

###################################################################################################################################################################################################################################################################################################################################
 
#This is where the magic happens, the code calls in laxpc software and passes the necessary arguments to it. The same is done with backshiftv2 and lcmath to generate the shifted bg and subrtact from the lightcurve

#One problem with this section is lcmath gives asynchronous warning for 1 second tbin, this does not happen with other time bins

lx_cmd = "printf \""+str(pcu_no)+" "+str(time_bin)+" "+str(anod)+"\n"+str(ch_min_binned)+" "+str(ch_max_binned)+"\n"+str(ch_nos)+"\n"+str(uld_bin)+" "+str(Erth_occ)+" "+ str(evt_flg)+"\n\" | ./laxpcl1"
#print lx_cmd
os.system(lx_cmd)

bg_cmd = "printf \""+str(pcu_no)+" "+str(anod)+"\n"+str(time_bin)+" "+str(ch_min_binned)+" "+str(ch_max_binned)+" "+str(uld_bin)+"\n"+"lxp"+str(pcu_no)+"level2back.pha"+"\n"+str(date_obs_yy)+" "+str(date_obs_mm)+" "+str(date_obs_dd)+"\n"+"/"+"\n\" | ./backshiftv2"
#print bg_cmd
os.system(bg_cmd)

lcmath_cmd = "printf \""+"lxp"+str(pcu_no)+"level2.lcurv"+"\n"+"lxp"+str(pcu_no)+"level2_shifted.lcbk\n"+"level2_"+str(pcu_no)+"_"+str(time_bin)+"s_"+str(int(E_min))+"_"+str(int(E_max))+"_shifted_backsub.lcurv\n 1.\n 1.\n no\n\" | lcmath"
#print lcmath_cmd
os.system(lcmath_cmd)


###################################################################################################################################################################################################################################################################################################################################
#This section moves all the essential files to another directory

lc_move = "mv "+directory_source_loc+"lxp"+str(pcu_no)+"level2.lcurv"+" "+directory_dest_loc+"lxp"+str(pcu_no)+"_"+str(anod)+"level2_"+str(time_bin)+"s_"+str(int(round(E_min)))+"_"+str(int(round(E_max)))+".lcurv"
os.system(lc_move)
backlc_move = "mv "+directory_source_loc+"lxp"+str(pcu_no)+"level2.lcbk"+" "+directory_dest_loc+"lxp"+str(pcu_no)+"_"+str(anod)+"level2_"+str(time_bin)+"s_"+str(int(round(E_min)))+"_"+str(int(round(E_max)))+".lcbk"
os.system(backlc_move)
gti_move = "mv "+directory_source_loc+"lxp"+str(pcu_no)+"level2.gti"+" "+directory_dest_loc+"lxp"+str(pcu_no)+"_"+str(anod)+"level2_"+str(time_bin)+"s_"+str(int(round(E_min)))+"_"+str(int(round(E_max)))+".gti"
os.system(gti_move)
pha_move = "mv "+directory_source_loc+"lxp"+str(pcu_no)+"level2.pha"+" "+directory_dest_loc+"lxp"+str(pcu_no)+"_"+str(anod)+"level2_"+str(time_bin)+"s_"+str(int(round(E_min)))+"_"+str(int(round(E_max)))+".pha"
os.system(pha_move)
backpha_move = "mv "+directory_source_loc+"lxp"+str(pcu_no)+"level2back.pha"+" "+directory_dest_loc+"lxp"+str(pcu_no)+"_"+str(anod)+"level2_"+str(time_bin)+"s_"+str(int(round(E_min)))+"_"+str(int(round(E_max)))+"back.pha"
os.system(backpha_move)
FITS_event_move = "mv "+directory_source_loc+"lxp"+str(pcu_no)+"level2.evn"+" "+directory_dest_loc+"lxp"+str(pcu_no)+"_"+str(anod)+"level2_"+str(time_bin)+"s_"+str(int(round(E_min)))+"_"+str(int(round(E_max)))+".evn"
os.system(FITS_event_move)
spectrum_move = "mv "+directory_source_loc+"lxp"+str(pcu_no)+"level2.spec"+" "+directory_dest_loc+"lxp"+str(pcu_no)+"_"+str(anod)+"level2_"+str(time_bin)+"s_"+str(int(round(E_min)))+"_"+str(int(round(E_max)))+".spec"
os.system(spectrum_move)

shift_lcbk_move = "mv "+directory_source_loc+"lxp"+str(pcu_no)+"level2_shifted.lcbk"+" "+directory_dest_loc+"lxp"+str(pcu_no)+"_"+str(anod)+"level2_"+str(time_bin)+"s_"+str(int(round(E_min)))+"_"+str(int(round(E_max)))+"back_shifted.lcbk"
os.system(shift_lcbk_move)
shift_pha_move = "mv "+directory_source_loc+"lxp"+str(pcu_no)+"level2back_shifted.pha"+" "+directory_dest_loc+"lxp"+str(pcu_no)+"_"+str(anod)+"level2_"+str(time_bin)+"s_"+str(int(round(E_min)))+"_"+str(int(round(E_max)))+"back_shifted.pha"
os.system(shift_pha_move)
shift_lc_move = "mv "+directory_source_loc+"lxp"+str(pcu_no)+"level2_shifted.lc"+" "+directory_dest_loc+"lxp"+str(pcu_no)+"_"+str(anod)+"level2_"+str(time_bin)+"s_"+str(int(round(E_min)))+"_"+str(int(round(E_max)))+"back_shifted.lc"
os.system(shift_lc_move)
#shift_spec_move = "mv "+directory_source_loc+"lxp"+str(pcu_no)+"level2back_shifted.spec"+" "+directory_dest_loc+"lxp"+str(pcu_no)+"_"+str(anod)+"level2_"+str(time_bin)+"s_"+str(int(round(E_min)))+"_"+str(int(round(E_max)))+"back_shifted.spec"
#os.system(shift_spec_move)
corr_lcsr_move = "mv "+directory_source_loc+"lxp"+str(pcu_no)+"level2_corr.lcsr"+" "+directory_dest_loc+"lxp"+str(pcu_no)+"_"+str(anod)+"level2_"+str(time_bin)+"s_"+str(int(round(E_min)))+"_"+str(int(round(E_max)))+"_corr.lcsr"
os.system(corr_lcsr_move)

shift_backsub_lc_move = "mv "+directory_source_loc+"level2_"+str(pcu_no)+"_"+str(time_bin)+"s_"+str(int(E_min))+"_"+str(int(E_max))+"_shifted_backsub.lcurv"+" "+directory_dest_loc+"lxp"+str(pcu_no)+"_"+str(anod)+"level2_"+str(time_bin)+"s_"+str(int(round(E_min)))+"_"+str(int(round(E_max)))+"_shifted_backsub.lcurv"
os.system(shift_backsub_lc_move)

###################################################################################################################################################################################################################################################################################################################################






