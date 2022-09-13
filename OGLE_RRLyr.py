#!/usr/bin/env python3

'''
Add a description of the code later
'''

# Standard library imports
from datetime import datetime
from pathlib import Path, os
import warnings, re, sys
from os import listdir
from os.path import isfile, join

# Third-party imports
from astropy import units as u
from astropy.coordinates import SkyCoord, Distance
from astropy.io import ascii
from astropy.table import Table, Column, vstack, hstack
from astroquery.simbad import Simbad
import inquirer as iq
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.patches as patch
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from sklearn.mixture import GaussianMixture
from sklearn import mixture, preprocessing
import scipy.stats
from scipy.stats import multivariate_normal
from termcolor import colored, cprint

# Local application imports
#from local_module import local_class
#from test_emcee_mod import log_prior	# testing! (get() functions...)

#################################  ==========================
### OGLE_RRLyr.py: TO-DO-LIST ###	Deadline: [...]
#################################  ==========================
#--------------------------------------------------
# 1) Finish membership() function			OK! ###
# 2) Compute <V>/<I> weighted by the memb.	OK!	###
# 3) Improve Fig. RA-DEC (circle & Gaia)	OK!	###
# 4) Improve PMs figure (2D-Gaussian)		**	### (FIX details?...)
# 5) Choose a new cluster from Harris		OK! ### Paper-done
# 6) Reload fits catalogues from Python		OK!	###	Paper-done
# 7) Astropy.Table() to save/read cats		OK!	### Paper-done
# 8) NGC6266: problem with RA/DEC clem		OK!	### Paper-done
# 9) Apply the Gaia quality-flags!!!		OK!	### Paper-done
# 10) Problem with prev-results (Vmag)		OK! ### Paper-done
# 11) Insert Tables in each function called from main()... --> remove
# 12) Memb. prob. from 2D-Gaussians			NO!	### ---> URGENT!!!
# 13) err_field is missing in GMM-2D			###	
# 14) Write the 5+ get functions in external!!	### (not needed for Paper)
# 15) Missing: Gaia catalogue of RRLs*			### (not needed for Paper)
# 16) Work in the case of 0 stars in OGLE*		### (not needed for Paper)
# 17) Test Astropy's sigma_clip 				### (not needed for Paper)
# 18) Cross-match with best photometry?			### (not needed for Paper)
#-- 	[...]	[...]	[...]	[...]			###
# 19) Save the PM-cleaned Gaia photometry		###
# 20) RDP with Gaia: AngSep and King profile	###
# 21) Put the alpha in prev file (review)		###
# 22) Recheck HP1 and other with few stars		### L: field has too big weight
# 23) Fix the membership scale in the VPD		###
#--------------------------------------------------

### DECIDE ABOUT RUWE AND NO_EXCESS FLAGS ###
#			RUWE < 2.0 is good for now!

def main():
	
	#.Required to import local modules:
	path = os.path.dirname(os.path.realpath(__file__))
	sys.path.insert(0, os.path.join(path))
	sys.path.insert(0, os.path.join(path, '..'))

	# To-do: Make local imports with clem, ogle, gaia and ir!!
	global NGC, c1, memb, tel, reHeader, ncat, p_member
	NGC, c1, RA_DEC, new_gc = choose_cluster()
	if not os.path.isdir(f'{path}/{NGC}'):
		os.mkdir(f'{path}/{NGC}')
		os.mkdir(f'{path}/{NGC}/Aux_Cats')

	ncat = iq.prompt({iq.Confirm('t', message='Download new catalogues?', \
					default=False)})['t']
	rrls, gaia_full, gaia_filt, rrl_gaia, rrl_gaiac, ir_filt, rrl_irc,\
		tel, rr_mag, rad = read_cats_Xmatch(NGC, c1, new_gc)
	
	if tel == 'vvv': Jname, Hname, Ksname = 'Jmag3','Hmag3','Ksmag3'
	else: Jname, Hname, Ksname = 'Jmag','Hmag','Kmag'
	Imag = rrl_gaiac['Imag']
	
	# FROM HERE: do not use read_catalogues with several global variables,
	# and retrieve only the useful variables from tables!
	
	# DELETE funct. read_catalogues(), after removing all other calls!!!
	#read_catalogues(NGC)
	memb = iq.prompt({iq.Confirm('memb', message='Use previous membership'+\
					' results?', default=False)})['memb']
	options = ['King profile (half-light radius)', 'Sigma-clipping (<V>,'+\
				' <I> & <Ks>)', 'Proper motions (from Gaia EDR3)']
	reHeader, head, alpha = not memb, not memb, 2.0
	pars = {}

	if not memb:
		recalc = iq.prompt({iq.Checkbox('rec', message='What do you want to'+\
					' recalculate? ', choices=options)})['rec']
		if options[0] in recalc:
			RA_DEC_Plane(memb, head, c1, RA_DEC, rad, gaia_filt, rrl_gaiac)

		if options[1] in recalc:
			### ARRUMAR PARA OS QUE TEM FILTRO B!::: URGENTE!
			if rr_mag == 'I': rr_mag = 'Ic'
			list_V = mag_vs_period(memb, rrl_gaiac, rr_mag, 2.5)
			if len(Imag) != 0 and Imag[0] != -99.99:
				list_I = mag_vs_period(memb, rrl_gaiac, 'I', alpha)		
			list_Ks = mag_vs_period(memb, rrl_irc, Ksname, alpha)
			list_mags = [list_V, list_I, list_Ks]
			pars['list_mags'] = [list_V, list_I, list_Ks]

		if options[2] in recalc:
			plot_CMDs(memb, tel, gaia_filt, rrl_gaiac, ir_filt, rrl_irc)
			ans_VPD = VPD_select(memb, gaia_filt, rrl_gaiac, c1, tel)
			c_rad, pms_RA_1d, pms_DE_1d, pms_2d, f_pms_2d = ans_VPD
			pars['c_rad'] = c_rad
			pars['pms_1d_2d'] = [pms_RA_1d, pms_DE_1d, pms_2d, f_pms_2d]
			### PAREI AQUI: TERMINAR LOGO PARA COMPLETAR O ARTIGO!!!
			if options[1] not in recalc:
				pars = save_read(memb, False, recalc, pars, alpha)
				list_mags = pars['list_mags']
			pars = membership(gaia_filt, rrl_gaiac, c1, pars, alpha)
		else:
			pars = save_read(memb, False, recalc, pars, alpha)
			pars = membership(gaia_filt, rrl_gaiac, c1, pars, alpha)
		
		pars = save_read(memb, True, recalc, pars, alpha)
		memb = True
	
	pars = save_read(memb, False, options, pars, alpha)
	p_member, list_mags = pars['memb_bellini'], pars['list_mags']
	#list_V, list_I, list_Ks = (list_mags)
	ans_VPD = VPD_select(memb, gaia_filt, rrl_gaiac, c1, tel, pars)
	if not memb: c_rad, pms_RA_1d, pms_DE_1d, pms_2d, f_pms_2d = ans_VPD
	else: c_rad = ans_VPD
	RA_DEC_Plane(memb, head, c1, RA_DEC, rad, gaia_filt, rrl_gaiac)

	mag_vs_period(memb, rrl_gaiac, rr_mag, 2.5, list_mags[0])
	if len(Imag) != 0 and Imag[0] != -99.99:
		mag_vs_period(memb, rrl_gaiac, 'I', alpha, list_mags[1])
	#mag_vs_period(memb, Ksmag, Ksname, alpha, list_mags[2])	### dimension problem!
	plot_CMDs(memb, tel, gaia_filt, rrl_gaiac, ir_filt, rrl_irc, p_member)
	print()

def sig_clip(mag, filt, alpha, numb, std=False):
	'''
	Computes the mean/error of <mag> (vs Period), using a sig-clipping method
	to avoid the effect of outliers (field RRLs). It's made through a series of
	iterations, computing med and removing stars outside med+/-alpha*sig.
	'''
	m1 = mag != -99.99
	stds,med_old = [],round(np.median(mag[m1]),3)
	std_old = round(np.std(mag[m1]),3) if len(mag[m1]) > 1 else 0.25
	stds.append(std_old)
	num = 2 if filt in ['V','B','Ic'] else 3
	print ('\033[1m  %d) <%s> vs. Period:\033[0m %.3f +/- %.3f (with outliers)'\
			%(num, filt, med_old, std_old))
	
	if len(mag[m1]) == 1:
		mean, med_new = np.mean(mag[m1]), np.median(mag[m1])
		err, std_new  = 0.25, 0.25
	
	else:
		for j in range(numb+1):
			if j > 0: med_old, std_old = med_new, std_new
			m2 = (mag > med_old - alpha*std_old) & (mag < med_old + alpha*std_old)
			med_new, std_new = np.median(mag[m2]), np.std(mag[m2])
			stds.append(round(std_new,3))
			
			if j == numb:
				if std: print (np.array(stds))
				mean, err  = np.mean(mag[m2]), np.std(mag[m2])/(len(mag[m2])**0.5)
				if not memb: print (u'    \u25B8 \u03C3-clipping: %.3f +/- %.3f '%(med_new,std_new)+\
						'(%.3f +/- %.3f)%s'%(mean,err,colored(u' \033[1m\u2713\033[0m','green')))

	return mean,err,med_new,std_new

def auto_ticks(ax):
	'''
	This aux function places ticks inside the plot in both axis, and sets major/
	minor ticks depending on the entries (xmaj,xmin...).
	'''
	ax.tick_params(direction='in',which='both')
	ax.xaxis.set_ticks_position('both')
	ax.yaxis.set_ticks_position('both')
	ax.xaxis.set_major_locator(ticker.AutoLocator())
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	ax.yaxis.set_major_locator(ticker.AutoLocator())
	ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

def choose_cluster():
	'''
	Write the description later!!!
	This function reads the cats properties and the number of stars in order to
	make a summary of these cats, to be printed in the terminal.
	'''
	warnings.simplefilter(action='ignore', category=FutureWarning)
	u_coord = (u.si.hourangle, u.si.deg)
	#os.system('printf "\n";python -V')
	#print (sys.version_info)
	print (f'{sys.version} --> {os.path.basename(__file__)} \n')
	
	new_gc = iq.prompt({iq.Confirm('new', message='New cluster or update'+\
						' catalogues?', default=False)})['new']
	path = os.path.dirname(os.path.realpath(__file__))
	harris = Table.read(f'{path}/harris_mwgc.fits', format='fits')
	IDs, names = np.array(harris['ID']), np.array(harris['Name'])
	RAs, DECs = np.array(harris['RA']), np.array(harris['DEC'])
	IDs = np.array([i.decode('utf-8').strip() for i in IDs])
	names = np.array([n.decode('utf-8').strip() for n in names])
	RAs = np.array([ra.decode('utf-8').strip() for ra in RAs])
	DECs = np.array([dec.decode('utf-8').strip() for dec in DECs])
	
	#.Insert cluster, validating from Harris catalogue
	if new_gc:
		IDs_case = [i.casefold() for i in IDs]
		names_case = [n.casefold() for n in names]
		check = list(IDs) + [n.casefold() for n in names if n != ''] + \
				IDs_case + [n for n in names if n != '']
		NGC = iq.prompt([iq.Text('c', message='Cluster name', validate= \
				lambda _,x: x.replace(' ','').casefold() in check)])['c']
		NGC = NGC.replace(' ','').casefold()
		if NGC in IDs_case: NGC = IDs[IDs_case.index(NGC)]
		else: NGC = IDs[names_case.index(NGC)]
	
	#.Choose previous cluster
	elif not new_gc:
		fold = ['__pycache__','Clement','GaiaEDR3','L_versus_B','OGLE-CVS',\
				'Paper','References','testCodes']
		GCs = [n for n in os.listdir(f'{path}') if '.' not in n and n not in fold]
		NGC = iq.prompt([iq.List(name='ngc', message='Choose the cluster', \
					choices=sorted(GCs), carousel=True)])['ngc']
	
	#.Get coordinates from Simbad or Harris
	t = Simbad.query_object(NGC)
	if len(t) >= 1:		# Simbad
		c = SkyCoord(t['RA'][0], t['DEC'][0], unit=u_coord)
	else:	# Harris
		RA, DEC = RAs[list(IDs).index(NGC)], DECs[list(IDs).index(NGC)]
		c = SkyCoord(RA, DEC, unit=u_coord)	
	RA_DEC = c.to_string('hmsdms')
	
	#.Print header [name(s), coords] in terminal
	if names[list(IDs).index(NGC)] != '':
		head = ' {} ({}) '.format(NGC, names[list(IDs).index(NGC)])
	else: head = ' {} '.format(NGC)
	l, s = u'\u2500'*(len(head)), ' '*4
	r, d = RA_DEC.split()[0], RA_DEC.split()[1]
	print (u'\n\033[1m{}\u250C{}\u2510\n{}\u2502{}\u2502'.format(s,l,s,head),
			u'{}RA {} ; DEC {}\n{}\u2514{}\u2518\n\33[0m'.format(s,r,d,s,l))

	return NGC, c, RA_DEC, new_gc

def get_clement(ngc, center):	# 23.mar.21: OK!
	'''
	Add description later...
	ESO 452-SC11 ==> 1636-283 in H10
	ESO 455-SC11 ==> HP1 in H10 (HauteProvince)
	ESO 333-SC016 ==> Ton2 in H10 (Tonantzintla)
	There are 6 lines of diff: 8 only H, 2 only C!!!
	'''
	path = os.path.dirname(os.path.realpath(__file__))
	c_clusters = Table.read(f'{path}/Clement/c1-clusters.fit', format='fits')
	c_radec = SkyCoord(ra=c_clusters['RAJ2000'], dec=c_clusters['DEJ2000'])
	c_rads = center.separation(c_radec)
	c_name = c_clusters['Name'][c_rads.argmin()].split(',')[0].replace(' ','')
	c_name = c_name.replace('Ruprecht','Rup').replace('2MASS','2MS')
	c_NED = c_clusters['NEDname'][c_rads.argmin()].replace(' ','')
	c_NED = c_NED.replace('Palomar', 'Pal')

	# Correct HP1 for Haute Provence 1!!!
	H10_notC = '(Ko1-2, Whiting1, BH176, FSR1735, Bh261, GLIMPSE1-2)'
	if ngc == c_NED or ngc == c_name.replace('-','') or c_name in \
		['ESO 452-SC11','ESO 455-SC11','ESO 333-SC016']:

		c_file = c_clusters['File'][c_rads.argmin()].strip()
		c_upd = '({})'.format(c_clusters['Update'][c_rads.argmin()].strip())
		c_Nv = c_clusters['Nv'][c_rads.argmin()]
		
		c_vars = Table.read(f'{path}/Clement/c2-variables.fit', format='fits')
		#3111-> 39, 1892,  76  ,   4   ,  26  , 986 ,  50  ,  4  ,  6   , 28
		ts = ['RR','RR0','RR01','RR01?','RR0?','RR1','RR1?','RR2','RR2?','RR?']
		t_flag = np.array([typ in ts for typ in c_vars['Type']])
		f_flag = np.array([f.strip() == c_file for f in c_vars['File']])
		flags = (t_flag & f_flag)
		rrl_epochs = list(c_vars['Epoch'][flags])
		if not(rrl_epochs and all(rrl == 'R0' for rrl in rrl_epochs)):
			print ('Not all RRLs have J2000 coordinates! Stop.')
		if c_Nv != len(c_vars[f_flag]):
			print ('Number of variables does not match! Stop.')
		
		rrl_ids = np.array(c_vars['Var'][flags], dtype='str')
		rrl_ras = c_vars['RAJ2000'][flags]
		rrl_decs = c_vars['DEJ2000'][flags]
		rrl_pers = np.array(c_vars['Per'][flags])
		rrl_mags = np.array(c_vars['__mag_'][flags])
		rrl_amps = np.array(c_vars['Amp'][flags])
		rrl_filts = np.array(c_vars['Filt'][flags], dtype='str')
		rrl_filts = np.where(rrl_filts==' ', '', rrl_filts)
		rrl_types = np.array(c_vars['Type'][flags], dtype='str')

		#..NGC6266: RA-DEC offset
		if ngc == 'NGC6266':
			rrl_ras = rrl_ras - 2.1e-4
			rrl_decs = rrl_decs + 1.12e-3

		#..NGC6441 (V36 - V104): ' ' -> 'V' (Layden+1999 & Pritzl+2001)
		if ngc == 'NGC6441': rrl_filts[:47] = 'V'

		rrls = (rrl_ids, rrl_ras, rrl_decs, rrl_pers, rrl_mags, rrl_amps, \
				rrl_filts, rrl_types)
		lab = ['ID','RA', 'DEC', 'cPeriod', 'mag', 'Ampl', 'Filt', 'cType']
		dt = [(np.unicode_,8), np.float64, np.float64, np.float32, np.float32,\
				np.float32, (np.unicode_,8), (np.unicode_,8)]
		fname = f'{path}/{ngc}/{ngc}_clem.fits'
		out = Table(rrls, names=lab, dtype=dt)
		out.replace_column('RA', np.around(rrl_ras,decimals=6))
		out.replace_column('DEC', np.around(rrl_decs,decimals=6))
		out.write(fname, format='fits', overwrite=True)
		return out, c_upd

	print ("Cluster is NOT in Clement's catalogue:\n{}.".format(H10_notC))
	return Table(), ''

def get_ogle(ngc, center):		# 23.mar.21: OK!
	'''
	Reads the OGLE catalogues, returning the RRLs inside a given radius.
	'''
	ns = ['Name', 'oType', 'RA', 'DEC']
	path = os.path.dirname(os.path.realpath(__file__))
	ident = ascii.read(f'{path}/L_versus_B/ident_edit.dat', names=ns, data_start=0)
	ident2 = ascii.read(f'{path}/L_versus_B/ident2_edit.dat', names=ns, data_start=0)
	ogleIV, c_unit = vstack([ident, ident2]), (u.si.hourangle, u.si.deg)
	o_coord = SkyCoord(ra=ogleIV['RA'], dec=ogleIV['DEC'], unit=c_unit)
	o_rad = center.separation(o_coord)

	use_Rt = iq.prompt({iq.Confirm('Rt', message='Use King tidal radius?',\
					default=False)})['Rt']
	if use_Rt:
		#..Call King function (Francisco)	#!!!
		pass
	else:
		valid, m = r'^[0-9]*(?:\.(\d|\d\d|\d\d\d))?$', 'Insert radius '
		lim = iq.prompt([iq.Text('lim', message=m+'(arcmin)', validate= \
				lambda _, x: re.match(valid, x) and float(x) <= 20)])['lim']
		lim = float(lim)
		inside = o_rad <= lim*u.si.arcmin
		print()

	#.Selecting stars inside the given radius
	ogle_in = ogleIV[inside]
	fname = '{}_ocvs_{}arcmin.fits'.format(ngc, str(lim).replace('.','p'))

	#.No OGLE stars
	if len(ogle_in) == 0:
		return Table(), lim

	#.Transforming RA-DEC (hmsdms->deg, 6 dec-digits)
	o_ras_deg = np.array(np.around(o_coord[inside].ra, decimals=6))
	o_decs_deg = np.array(np.around(o_coord[inside].dec, decimals=6))
	ns, dt = ('RA','DEC'), list(np.repeat([np.float64],2))
	coords = Table(np.c_[o_ras_deg, o_decs_deg], names=ns, dtype=dt)
	ogle_in.remove_columns(['RA','DEC'])
	ogle_in.add_columns(list(coords.columns.values()))	# corrected 20.fev.2021 --> continue from here

	#.Extra columns with photometry and period
	abc = ['ID', 'Imag', 'Vmag', 'Period', 'uPeriod', 'Time-max', 'Iamp',\
			'R_21', 'phi_21', 'R_31', 'phi_31']
	d = abc + ['Period2', 'uPeriod2', 'Time-max2', 'Iamp2', 'R_21f', \
			'phi_21f','R_31f','phi_31f']
	blg_ab = ascii.read(f'{path}/OGLE-CVS/blg/RRab.dat', names=abc, data_start=0)
	blg_c = ascii.read(f'{path}/OGLE-CVS/blg/RRc.dat', names=abc, data_start=0)
	blg_d = ascii.read(f'{path}/OGLE-CVS/blg/RRd.dat', names=d, data_start=0)
	blg_ad = ascii.read(f'{path}/OGLE-CVS/blg/aRRd.dat', names=d, data_start=0)
	gd_ab = ascii.read(f'{path}/OGLE-CVS/gd/RRab.dat', names=abc, data_start=0)
	gd_c = ascii.read(f'{path}/OGLE-CVS/gd/RRc.dat', names=abc, data_start=0)
	gd_d = ascii.read(f'{path}/OGLE-CVS/gd/RRd.dat', names=d, data_start=0)
	gd_ad = ascii.read(f'{path}/OGLE-CVS/gd/aRRd.dat', names=d, data_start=0)
	
	#.Declaring the phot arrays
	o_imags = np.zeros(len(ogle_in))
	o_vmags, o_pers = np.copy(o_imags), np.copy(o_imags)
	o_upers, o_amps = np.copy(o_imags), np.copy(o_imags)
	name = {'BLG-RRab': blg_ab, 'BLG-RRc': blg_c, 'BLG-RRd': blg_d,\
			'BLG-aRRd': blg_ad, 'GD-RRab': gd_ab, 'GD-RRc': gd_c, \
			'GD-RRd': gd_d, 'GD-aRRd': gd_ad}

	#.Looping to add phot columns
	for (idx, rrl) in enumerate(ogle_in):
		key = '{}-{}'.format(rrl['Name'].split('-')[1], rrl['oType'])
		good = name[key]['ID'] == rrl['Name']
		o_imags[idx] = float(name[key]['Imag'][good])
		try:	# lines with '-' --> NaN
			o_vmags[idx] = name[key]['Vmag'][good]
		except ValueError:
			o_vmags[idx] = np.nan
		o_pers[idx] = name[key]['Period'][good]
		o_upers[idx] = name[key]['uPeriod'][good]
		o_amps[idx] = name[key]['Iamp'][good]
	
	#.Appending new columns and save
	cols = np.c_[o_imags, o_vmags, o_pers, o_upers, o_amps]
	ns = ('Imag','o_Vmag','o_Period','o_uPeriod','Iamp')
	dt = list(np.repeat([np.float32],5))
	phot = Table(cols, names=ns, dtype=dt)
	ogle_in.add_columns(list(phot.columns.values()))
	ogle_in.write(f'{path}/{ngc}/{fname}', format='fits', overwrite=True)

	return ogle_in, lim

def clean_gaia2(ngc, rad, full_t, bulge=True):	# 18.mar.21: OK!
	'''
	WRITE DOCS LATER!!!
	It is separated from get_gaia, so we can use it in other codes!!!
	Read Arenou+2018. fullt needs to be an AstroPy table!
	phot_bp_rp_excess_factor is the ratio of the sum of G_BP and G_RP fluxes
	over the G flux and should be around one for normal stars.
	'''
	#.Defining the main variables
	chi2, nu = full_t['chi2AL'], full_t['NgAL'] - 5
	u_filt = np.sqrt(chi2/nu)
	E_fac = full_t['E_BR_RP_']

	#.Defining the flags from Equations 1, 2 and 3
	flag1 = np.ones(len(full_t),dtype=bool)
	flag2a, flag2b, flag3 = np.copy(flag1), np.copy(flag1), np.copy(flag1)
	for (i, star) in enumerate(full_t):
		flag1[i] = u_filt[i] < 1.2*max(1, np.exp(-0.2*(star['Gmag']-19.5)))
		flag2a[i] = E_fac[i] > 1.0 + 0.015*(star['BPmag']-star['RPmag'])**2
		flag2b[i] = E_fac[i] < 1.3 + 0.06*(star['BPmag']-star['RPmag'])**2
		flag3[i] = star['Nper'] > 8

	#.Normal: Equations 1 and 2
	if not bulge: flags = flag1 & flag2a & flag2b
	#.Preferable for bulge proper motions
	elif bulge: flags = flag1 & flag3

	#.Returning the output table
	#print ('  > Gaia ({}\'): Before {} stars / After {} stars'.format(rad,\
	#		len(full_t), len(full_t[flags])) + ' (Arenou+2018)')
	return full_t[flags]

def get_gaia2(ngc, center, radius):		# 18.mar.21: OK!
	'''
	This	(I/345/gaia2)		WRITE DOCS LATER!!!
	For those with a well-measured parallax, it is possible to convert the
	epochs from 2015.5 to J2000 (January 1st, 2020).
	There is one line more in astroquery.Vizier?
	Method 1 transformed to J2000 is exactly equal the Vizier RAJ2000 !!! OK!
	(But since there are some stars without parralax, I will use Vizier)
	'''

	#print('Entrando em get_gaia2')
	path = os.path.dirname(os.path.realpath(__file__))
	rad, pm_u = radius*u.si.arcmin, u.si.mas/u.si.yr

	## Method 1: astroquery.Gaia (slow, epoch 2015.5 to J2000)
	# Links: https://astroquery.readthedocs.io/en/latest/gaia/gaia.html
	# https://docs.astropy.org/en/stable/coordinates/apply_space_motion.html
	# [Bailer-Jones (2015, 1507.02105); Bailer-Jones+ (2018, 1804.10121)
	'''from astroquery.gaia import Gaia
	from astropy.time import Time
	j = Gaia.cone_search_async(center, radius)
	tab = j.get_results()
	c0 = SkyCoord(ra=tab['ra'][0]*u.si.deg, dec=tab['dec'][0]*u.si.deg,\
				distance=Distance(parallax=tab['parallax'][0] * u.si.mas),\
				pm_ra_cosdec=tab['pmra'][0]*pm_u,pm_dec=tab['pmdec'][0]*pm_u,\
				obstime=Time(tab['ref_epoch'][0], format='decimalyear'))
	epoch_J2000 = Time('2000-01-01')
	c_2mass_epoch = c0.apply_space_motion(epoch_J2000) # Checked!'''

	## Method 2: Get also RAJ2000 and DEJ2000 (fast, Vizier)
	# Link: https://astroquery.readthedocs.io/en/latest/vizier/vizier.html
	# I tried to change e-/s, but it's not good (TOPCAT)
	from astropy.io.fits.verify import VerifyWarning
	warnings.simplefilter(action='ignore', category=u.UnitsWarning)
	warnings.simplefilter(action='ignore', category=VerifyWarning)

	from astroquery.vizier import Vizier
	#rad = 1 * u.si.arcmin		### change later: only testing!!!
	std_cols = Vizier(row_limit=1000000).query_region(center, radius=rad,\
				catalog='I/345/gaia2')
	v1 = Vizier(columns=['RAJ2000','e_RAJ2000','DEJ2000','e_DEJ2000'],\
				row_limit=1000000)
	v2 = Vizier(columns=['NgAL','chi2AL','Nper','phot_bp_rp_excess_factor'],
				row_limit=1000000)
	J2000 = v1.query_region(center, radius=rad, catalog='I/345/gaia2')
	flags = v2.query_region(center, radius=rad, catalog='I/345/gaia2')
	gaia2 = hstack([std_cols[0], J2000[0], flags[0]])

	#.Reordering the columns (J2000, flags, std_cols)
	order = ['RAJ2000','e_RAJ2000','DEJ2000','e_DEJ2000','Source','NgAL',\
			'chi2AL','Nper','E_BR_RP_','Plx','e_Plx','pmRA','e_pmRA','pmDE',\
			'e_pmDE','Dup','FG','e_FG','Gmag','e_Gmag','FBP','e_FBP','BPmag',\
			'e_BPmag','FRP','e_FRP','RPmag','e_RPmag','BP-RP','RV','e_RV',\
			'Teff','AG','E_BP-RP_','Rad','Lum','RA_ICRS','e_RA_ICRS',\
			'DE_ICRS','e_DE_ICRS']
	gaia2_full = gaia2[order]
	
	#.Changing masked values to np.NaN -> Method 1 (quick!)
	def fill_cols(tbl, fill=np.nan, kind='f'):
		"""
		In-place fill of ``tbl`` columns which have dtype ``kind``
		with ``fill`` value.
		"""
		for col in tbl.itercols():
			if col.dtype.kind == kind:
				col[...] = col.filled(fill)
	fill_cols(gaia2_full, fill=np.nan, kind='f')	### IMPORTANT!!!
	fill_cols(gaia2_full, fill=-99, kind='i')

	#.Changing masked values to np.NaN -> Method 2 (slower)
	# for col in gaia2_full.colnames:
	# 	for (idx,elem) in enumerate(gaia2_full[col]):
	# 		if elem is np.ma.masked or elem == 1.0E20:
	# 			gaia2_full[col][idx] = np.NaN

	#.Write full table (before and after Arenou+18 filtering)
	fname = '{}_gaia2_{}arcmin.fits'.format(ngc,str(radius).replace('.','p'))
	gaia2_full.write(f'{path}/{ngc}/{fname}', format='fits', overwrite=True)
	gaia2_filt = clean_gaia2(ngc, radius, gaia2_full)
	fname = fname.replace('gaia2','gaia2c')
	gaia2_filt.write(f'{path}/{ngc}/{fname}', format='fits', overwrite=True)

	return gaia2_full, gaia2_filt

def get_gaia3(ngc, center, radius):		# 18.mar.21: OK!

	#print('Entrando em get_gaia3')
	from astroquery.gaia import Gaia
	from GaiaEDR3 import GaiaEDR3_data as gaiaedr3

	cwd = os.getcwd()
	path = os.path.dirname(os.path.realpath(__file__))
	path = os.path.join(path, 'GaiaEDR3')
	Gaia.login(credentials_file=f'{path}/my_credentials.txt')
	ngc, c, rad, out, outname, lind20 = gaiaedr3.initial_inputs(ngc, radius)
	gaia3, gaia3_valid = gaiaedr3.astroquery_adql(ngc, c, rad, outname)
	#success = gaiaedr3.upload_to_gaia(gaia3, 'tablename', outname)
	gaia3_corr = gaiaedr3.gaia_codes(gaia3, gaia3_valid, lind20, outname)
	#print (center == c1)

	return gaia3, gaia3_corr

def get_vvv_2mass(ngc, center, radius, full=False):		# 23.mar.21: OK!
	'''
	WRITE DESCRIPTION (DOCS) LATER
	Get the IR data from VVV (preferable) or 2MASS (ABDC quality flags)
	Don't save/return the unfiltered catalogue, only the filtered one!!!
	The 'full' kwarg determines if the ir_full Table will be saved.
	'''
	from astroquery.vizier import Vizier
	path = os.path.dirname(os.path.realpath(__file__))
	rad = radius * u.si.arcmin

	#.Check if cluster is inside VVV-DR2 covered area
	col1 = ['RAJ2000','DEJ2000','iauname','mClass','Z-Y','e_Z-Y','J-H',\
			'e_J-H','H-Ks','e_H-Ks','Zmag3','e_Zmag3','Zperrbits','Ymag3',\
			'e_Ymag3','Yperrbits']
	col2 = ['Jmag3','e_Jmag3','Jperrbits','Hmag3','e_Hmag3','Hperrbits',\
			'Ksmag3','e_Ksmag3','Ksperrbits']
	v0 = Vizier(row_limit=1000000, columns=col1, catalog='II/348/vvv2')
	v1 = Vizier(row_limit=1000000, columns=col2, catalog='II/348/vvv2')
	try:
		tab0 = v0.query_region(center, radius=rad)['II/348/vvv2']
		tab1 = v1.query_region(center, radius=rad)['II/348/vvv2']
		ir_full = hstack([tab0, tab1])
		tel, f = 'vvv', 'perrbits'
	except TypeError:
		tel = '2mass'
	
	#.VVV-DR2 photometry tables (with flags)
	if tel == 'vvv':
		fmag = ['Zmag3', 'Ymag3', 'Jmag3', 'Hmag3', 'Ksmag3']
		fcol = ['Z-Y', 'J-H', 'H-Ks']
		ir_filt = Table(ir_full)
		for star in ir_filt:
			for (idx, mag) in enumerate(fmag):
				if star[mag[:-4]+f] > 255 or star[mag] is np.ma.masked:
					star[mag], star['e_'+mag] = np.nan, np.nan
			for (idx, col) in enumerate(fcol):
				if star[col] is np.ma.masked:
					star[col], star['e_'+col] = np.nan, np.nan

	#.2MASS photometry tables (with flags)
	elif tel == '2mass':
		c = ['RAJ2000','DEJ2000','2MASS','Jmag','e_Jmag','Hmag','e_Hmag',\
			'Kmag','e_Kmag','Qflg','Rflg','Bflg','Cflg','Xflg','Aflg']
		v2 = Vizier(row_limit=1000000, columns=c, catalog='II/246/out')
		ir_full = v2.query_region(center, radius=rad)['II/246/out']
		ir_filt, val_Qflg = Table(ir_full), ['A', 'B', 'C', 'D']
		fmag = ['Jmag', 'Hmag', 'Kmag']
		for star in ir_filt:
			for (idx, mag) in enumerate(fmag):
				if star['Qflg'][idx] not in val_Qflg:
					star[mag], star['e_'+mag] = np.nan, np.nan
		ir_filt.meta['description'] = ir_filt.meta['description'][:55]
	
	#.Write/Return ONLY the filtered IR table (VVV or 2MASS)
	rad_str = str(radius).replace('.','p')
	fname = f'{ngc}_{tel}c_{rad_str}arcmin.fits'
	ir_filt.write(f'{path}/{ngc}/{fname}', format='fits', overwrite=True)
	if full:
		fname = fname.replace(tel+'c', tel)
		ir_full.write(f'{path}/{ngc}/{fname}', format='fits', overwrite=True)

	return ir_filt, tel

def read_cats_Xmatch(NGC, c, new_gc):	# 23.mar.21: OK!
	'''
	Write a description (and docs) later.
	Reads all catalogues and make the cross-matches!!!
	Clement, OGLE, Xmatch, Gaia(Xmatch), IR (Xmatch)
	Return: Xmatch, gaia_ [...]
	'''
	from astropy.io.fits.verify import VerifyWarning
	warnings.simplefilter(action='ignore', category=u.UnitsWarning)
	warnings.simplefilter(action='ignore', category=VerifyWarning)

	global n_gaiac	# remove!
	u_coord = (u.si.deg, u.si.deg)
	path = os.path.dirname(os.path.realpath(__file__))

	#.Calling the get_clement() function
	if new_gc:
		c_rrls, c_upd = get_clement(NGC, c)
	elif Path(f'{path}/{NGC}/{NGC}_clem.fits').is_file():
		c_rrls = Table.read(f'{path}/{NGC}/{NGC}_clem.fits', format='fits')
	else: c_rrls = Table()

	#.Calling the get_ogle() function
	fs = [f for f in listdir(f'{path}/{NGC}') if 'ocvs' in f and '.fits' in f]
	if new_gc and ncat:
		o_rrls, o_lim = get_ogle(NGC, c)
	elif len(fs) > 0:
		o_rrls = Table.read(f'{path}/{NGC}/{fs[0]}', format='fits')
		o_lim = float(fs[0].split('_')[-1].split('arc')[0].replace('p','.'))
	else:
		o_rrls, o_lim = Table(), 0

	#.Cross-match OGLE (source) vs. Clement
	breakpoint()
	if len(o_rrls) == 0:
		print('\n\033[1m * Cluster is not in the OGLE footprint *\033[0m\n')
		sys.exit()
	c_ogle = SkyCoord(o_rrls['RA'], o_rrls['DEC'], unit=u_coord)
	c_clem = SkyCoord(c_rrls['RA'], c_rrls['DEC'], unit=u_coord)
	idx, d2d, d3d = c_ogle.match_to_catalog_sky(c_clem)
	max_sep = 4.5 if NGC=='NGC6266' else 1.0	# important!
	sep_constraint = d2d < max_sep * u.si.arcsec
	clem_matches = c_rrls[idx[sep_constraint]]

	#.Adding OGLE stars without match
	c_empty_row = ['', np.nan, np.nan, np.nan, np.nan, np.nan, '', '']
	n_new, sort_idxs = 0, np.zeros(len(c_ogle), dtype=int)
	for i in range(len(c_ogle)):
		if not sep_constraint[i]:
			clem_matches.add_row(c_empty_row)
			sort_idxs[len(clem_matches)-1] = i
			n_new += 1
		else: sort_idxs[i-n_new] = i
	c_mixed = clem_matches[np.argsort(sort_idxs)]
	
	#.Adding Clement stars without match
	for (i, item) in enumerate(c_rrls['ID']):
		if item not in c_mixed['ID']:
			c_mixed.add_row(c_rrls[i])
	
	#.OGLE+Clement tables and write output
	o_empty_row = ['', ''] + list(np.repeat(np.nan,7))
	n_clem_out = len(c_rrls)-len(c_rrls[idx[sep_constraint]])
	for i in range(n_clem_out):
		out = [c_mixed['RA'][-n_clem_out+i], c_mixed['DEC'][-n_clem_out+i]]
		o_empty_row = ['', ''] + out + list(np.repeat(np.nan,5))
		o_rrls.add_row(o_empty_row)
	if o_rrls: del c_mixed['RA', 'DEC']
	rrls = hstack([o_rrls, c_mixed], join_type='exact')
	rrls.write(f'{path}/{NGC}/{NGC}_Xmatch.fits', format='fits', overwrite=True)

	#.Organize prints in the terminal
	n_ogle, n_tot = len(c_ogle), len(c_rrls)+n_new
	if not c_rrls: print('  > Clement: No RR Lyrae')
	if not new_gc: c_upd = ''
	else: print(f'  > Clement: {len(c_rrls)} RR Lyrae {c_upd}')
	if o_rrls:
		col = colored(f'{n_tot} in total!\033[0m', 'yellow')
		print(f'  > OGLE ({str(o_lim)}\'): {n_ogle} RR Lyrae\033[1m',
				f'({n_new} new) ==> {col}\n')
	else:
		print(f'There are no OGLE RRLs inside rad={str(o_lim)}.\n')

	#.Calling the get_gaia2() or get_gaia3() function
	dr = 'gaia3'
	rad_str = str(o_lim).replace('.','p') + 'arcmin'
	cat = [f for f in listdir(f'{path}/{NGC}') if f'{dr}_c{rad_str}' in f]
	source = {'gaia2': 'Source', 'gaia3': 'source_id'}
	
	if new_gc and ncat:
		if dr == 'gaia3':
			gaia3_orig, gaia_full = get_gaia3(NGC, c, o_lim)
		else: gaia_full, gaia_filt = get_gaia2(NGC, c, o_lim)	# DR2
	elif cat:
		filt = cat[0].replace(f'{dr}', f'{dr}c')
		full = cat[0] if dr == 'gaia2' else filt
		gaia_full = Table.read(f'{path}/{NGC}/{full}', format='fits')
		gaia_filt = Table.read(f'{path}/{NGC}/{filt}', format='fits')
	else:
		gaia_full, gaia_filt = Table(), Table()
	
	if dr == 'gaia3':
		# gaia_filt = gaia_full[gaia3_full['all_flags'] == 1]
		# gaia_filt = gaia_full[gaia_full['flag_noExcess'] == 1]
		flag_noExcess_ruwe2p0 = (gaia_full['flag_ruwe_2p0'] == 1) & \
								(gaia_full['flag_gmag'] == 1) & \
								(gaia_full['flag_rpmag'] == 1)			### TESTAR SE DEU CERTO!!!
		gaia_filt = gaia_full[flag_noExcess_ruwe2p0]
	print (f'\n  > {dr} ({str(o_lim)}\'): Before {len(gaia_full)} stars',\
			f'/ After {len(gaia_filt)} stars (cat. validation)')

	# breakpoint()

	def remove_duplicates(tab1, tab2, source_str):
		unique = np.unique(tab2[source_str], return_counts=True)[0]
		counts = np.unique(tab2[source_str], return_counts=True)[1]
		duplic = []

		for i in range(len(counts)):	# remove duplicated entries: ok
			if counts[i] > 1: duplic.append(unique[i])
		for dup in duplic:
			idxs = np.nonzero(tab2[source_str] == dup)[0]
			r = 0 if d2d[sep_const][idxs[0]] > d2d[sep_const][idxs[1]] else 1
			tab1.remove_row(idxs[r])
			tab2.remove_row(idxs[r])

		return hstack([tab1, tab2], join_type='exact')

	#.Cross-match OGLE+Clem vs. Gaia cleaned
	c_match = SkyCoord(rrls['RA'], rrls['DEC'], unit=u_coord)
	c_gaiac = SkyCoord(gaia_filt['RAJ2000'], gaia_filt['DEJ2000'], unit=u_coord)
	idx, d2d, d3d = c_match.match_to_catalog_sky(c_gaiac)
	sep_const = d2d < 1.0 * u.si.arcsec
	rrlc_matches = rrls[sep_const]
	gaiac_matches = gaia_filt[idx[sep_const]]
	rrl_gaiac = remove_duplicates(rrlc_matches, gaiac_matches, source[dr])
	fnamec = f'{path}/{NGC}/{NGC}_Xmatch_wGaiac.fits'
	rrl_gaiac.write(fnamec, format='fits', overwrite=True)

	#breakpoint()
	
	#.Cross-match OGLE+Clem vs. Gaia full
	c_gaia = SkyCoord(gaia_full['RAJ2000'] ,gaia_full['DEJ2000'], unit=u_coord)
	idx, d2d, d3d = c_match.match_to_catalog_sky(c_gaia)
	sep_const = d2d < 1.0 * u.si.arcsec
	rrl_matches = rrls[sep_const]
	gaia_matches = gaia_full[idx[sep_const]]
	rrl_gaia = remove_duplicates(rrl_matches, gaia_matches, source[dr])
	fname = f'{path}/{NGC}/{NGC}_Xmatch_wGaia.fits'
	rrl_gaia.write(fname, format='fits', overwrite=True)

	#.Gaia vs RR Lyrae cross-match: Final prints
	n_gaia, n_gaiac = len(rrl_gaia), len(rrl_gaiac)
	n_miss = n_tot-len(rrl_gaiac)
	final = f'({n_miss} missing)' if n_miss else '(all)'
	print (f'  > Gaia vs RRLs: {n_gaia} / {n_gaiac} RR Lyrae {final}')
	# breakpoint()

	#.Calling the get_vvv_2mass() function
	filt_f = [f for f in listdir(f'{path}/{NGC}') if ('vvvc' in f or \
				'2massc' in f) and 'arcmin.fits' in f]
	if new_gc and ncat:
		ir_filt, tel = get_vvv_2mass(NGC, c, o_lim)
	elif len(filt_f) > 0:
		ir_filt = Table.read(f'{path}/{NGC}/{filt_f[0]}', format='fits')
		tel = filt_f[0].split('_')[1][:-1]
	else:
		ir_filt, tel = Table()	, '2mass'

	#.Cross-match OGLE+Clem vs. IR data (clean)
	c_match = SkyCoord(rrls['RA'], rrls['DEC'], unit=u_coord)
	c_ir = SkyCoord(ir_filt['RAJ2000'], ir_filt['DEJ2000'], unit=u_coord)
	idx, d2d, d3d = c_match.match_to_catalog_sky(c_ir)
	sep_const = d2d < 1.0 * u.si.arcsec
	rrl_matches = rrls[sep_const]
	ir_matches = ir_filt[idx[sep_const]]

	#.OGLE+Clem + IRc data, write output and print
	label = 'iauname' if tel == 'vvv' else '_2MASS'
	rrl_irc = remove_duplicates(rrl_matches, ir_matches, label)
	fname = f'{path}/{NGC}/{NGC}_Xmatch_w{tel}c.fits'
	rrl_irc.write(fname, format='fits', overwrite=True)
	Kname = 'Ksmag3' if tel == 'vvv' else 'Kmag'
	Kgood = len(ir_filt[~np.isnan(ir_filt[Kname])])
	print (f'\n  > {tel.upper()} ({str(o_lim)}\'): Before {len(ir_filt)}',
			f'stars / After {Kgood} stars (phot-flags)')
	n_irc, n_miss = len(rrl_irc), n_tot-len(rrl_irc)
	final = f'({n_miss} missing)' if n_miss else '(all)'
	print (f'  > {tel.upper()} vs RRLs: {n_irc} RR Lyrae {final}')

	#.Print missing and returning all the important tables
	input(colored('\n%s <<< Press ENTER to start >>> '%(' '*19),'green'))
	from collections import Counter
	ct = Counter(rrls['Filt'][rrls['Filt'] != ''])
	rr_mag, n_mag = ct.most_common()[0][0], ct.most_common()[0][1]
	if rr_mag == 'P':
		rr_mag, n_mag = ct.most_common()[1][0], ct.most_common()[0][1]
	print (f'\n  > {n_tot-n_mag} RR Lyr without <{rr_mag}> value (Clement)')
	print (f'  > {n_tot-n_ogle} RR Lyr without <I> value (OGLE)')
	Kmiss = n_tot - len(rrl_irc[Kname][~np.isnan(rrl_irc[Kname])])
	print (f'  > {Kmiss} RR Lyr without Ks value ({tel.upper()})\n')

	return rrls, gaia_full, gaia_filt, rrl_gaia, rrl_gaiac, ir_filt,\
			rrl_irc, tel, rr_mag, o_lim

def RA_DEC_Plane(memb, head, cen, RA_DEC, rad, gaia_filt, rrl_gaiac):	# 24.mar.21: OK!
	'''
	Plots a RA-DEC plane with the spatial distribution of the catalogued RRLs.
	There are two types: with or without membership colors. The SOAR version has
	the extensive computations of the Circle, with SkyCoord.
	'''
	if head:
		print(f"\033[1m  1) RA-DEC plane:\033[0m RR Lyrae (r={rad}')", end='')
		print(colored(u' \033[1m\u2713\033[0m','green'))
	if not memb: fig = plt.figure(figsize=(6.19,6))
	else:
		# fig = plt.figure(figsize=(6.19,6.63))
		fig = plt.figure(figsize=(6.1,6))
	ax = plt.gca()
	
	cen_ra, cen_dec = cen.ra.deg, cen.dec.deg
	ra, dec = rrl_gaiac['ra'], rrl_gaiac['dec']
	cividisBig = cm.get_cmap('plasma_r', 512)
	newcmp = ListedColormap(cividisBig(np.linspace(0.1, 0.95, 256)))
	if not memb:
		ax.scatter(ra, dec, s=10, color='red', marker='o', zorder=3,\
					label='RR Lyrae')
	else:
		pcm = plt.scatter(ra, dec, s=25, c=p_member, cmap=newcmp, label=\
						'RR Lyrae', zorder=3)		# norm=colors.PowerNorm(gamma=1/8)
		#divider = make_axes_locatable(ax)
		#cax = divider.append_axes('top', size='2%', pad=0.3)
		#cbar = plt.colorbar(pcm, cax=cax, orientation='horizontal', ticks=\
		#		[0.0,0.2,0.4,0.6,0.8,1.0])
		#plt.clim(0, 1)
		#cax.xaxis.set_ticks_position('top')
		#cbar.set_label('Membership', fontsize=11)
	
	ax.scatter(cen_ra, cen_dec, s=150, color='black', marker='+', zorder=4, \
				label='Center', alpha=0.8)
	ax.scatter(gaia_filt['ra'], gaia_filt['dec'], s=0.15, marker='x',\
				color='gainsboro', label=r'$Gaia_{\rm filt}$')
	
	# !!! CURIOUS !!! Circle vs. Ellipse
	#circ = plt.Circle((cen_ra, cen_DE),rad/60,ls='--',color='k',zorder=5,fill=False)
	#ax.add_artist(circ)

	u, v = cen_ra, cen_dec				# x and y-positions of the center
	a = (rad/60)/np.cos(cen_dec*np.pi/180)		#radius on the x-axis
	b = rad/60    							#radius on the y-axis
	t = np.linspace(0, 2*np.pi, 200)
	ax.plot(u + a*np.cos(t), v + b*np.sin(t), ls='--', color='green', label=\
			r'r = $%d\,^{\prime}$'%rad)
	#for it in t: plt.plot([u, u + a*np.cos(it)],[v, v + b*np.sin(it)])
	#ellipse = patch.Ellipse((u,v), 2*a, 2*b, angle=0, ls='--', fill=False)
	#ax.add_artist(ellipse)

	#angles = np.linspace(0,360,25)
	#for angle in angles:
	#	l_angle = angle * np.pi/180
	#	DE_l = DEc + (rad/60)*np.sin(l_angle)
	#	RA_l = RAc + (rad/60)*np.cos(l_angle)/np.cos(DE_l*np.pi/180)
	#	plt.plot([RAc,RA_l],[DEc,DE_l],zorder=1000)

	# ax.update({'xlabel':'Right ascension (deg)', 'ylabel':'Declination (deg)'})
	ax.set_xlabel('Right ascension (deg)', fontsize=15.5)
	ax.set_ylabel('Declination (deg)', fontsize=15.5)
	ax.tick_params(axis='both', which='major', labelsize=14)
	auto_ticks(ax)
	ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
	ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.02))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
	ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.02))
	ax.invert_xaxis()
	ax.legend(loc=1, fontsize=14.5)
	plt.title(r'%s$\,$%s: %s$^{\rm{h}}$%s$^{\rm{m}}$%s$^{\rm{s}}$ %s'%\
				(NGC[:3],NGC[3:],RA_DEC.split()[0][:2],RA_DEC.split()[0][3:5],\
				RA_DEC.split()[0][6:-1], RA_DEC.split()[1][:3]) + r'$^{\circ}$'\
				r'%s$^{\prime}$%s$^{\prime\prime}$'%(RA_DEC.split()[1][4:6],\
				RA_DEC.split()[1][7:-1]),fontsize=17)
	plt.tight_layout()
	plt.subplots_adjust(bottom=0.098, right=0.967, left=0.158, top=0.927) # bottom=0.073
	plt.show()

def mag_vs_period(memb, rrls, filt, alpha, list_med=[]):		# 29.mar.21: OK!
	'''
	This function generates the <mag> vs Period plot, using the sig-clipping
	function above, and can be applied to any magnitude.
	'''
	#tel, Jname, Hname, Ksname = 'VVV','Jmag3','Hmag3','Ksmag3'
	#else: tel, Jname, Hname, Ksname = '2mass','Jmag','Hmag','Kmag'
	#mag_dict = {'V':'mag', 'B':'mag', 'Ic': 'mag', 'I':'Imag', 'Ksmag3':filt, 'Kmag':filt}
	num = 3 if filt == 'I' else 2	# falta aqui tbm
	mag_dict = {('V','B','Ic'): 'mag', ('I'):'Imag', ('Ksmag3','Kmag'):filt}
	for (i, item) in enumerate([*mag_dict]):
		if filt in item: break
	mag = rrls[mag_dict[[*mag_dict][i]]]
	if filt in ['Ksmag3', 'Kmag']: filt = 'Ks'
	period = rrls['o_Period']
	ra = rrls['ra'] if filt != 'Ks' else rrls['RA']
	dec = rrls['dec'] if filt != 'Ks' else rrls['DEC']

	#.Correction if there is more than one Filt
	if filt in ('V','B','Ic'):
		f = rrls['Filt'] == filt[0]
		mag, ra, dec, period = mag[f], ra[f], dec[f], period[f]
	
	if not memb:
		list_med = []
		mean, err, med, std = sig_clip(mag[~np.isnan(mag)], filt, alpha, 20)
		list_med.extend([mean, err, med, std])
		#print (list_med)
	else:
		mean, err = list_med[0], list_med[1]
		med, std = list_med[2], list_med[3]

		if reHeader:
			print(f'\033[1m  {num}) <{filt}> vs. Period:\033[0m ', end='')
			print(f'{mean:.3f} +/- {err:.3f} ({med:.3f} +/- {std:.3f})', end='')
			print(colored(u' \033[1m\u2713\033[0m','green'))
	
	if memb:
		f = rrls['Filt'] == filt[0]
		pmemb = np.array(p_member)[f] if filt in ('V','B','Ic') else np.array(p_member)
	if NGC == 'NGC6717':
		err = 0.100 if filt in ('V','B','Ic') else 0.080

	if filt == 'Ic': filt = 'I'
	flag = abs(mag-med) <= alpha*std
	ra.unit, dec.unit = u.si.deg, u.si.deg

	c2 = SkyCoord(ra, dec, frame='icrs')
	sep = c1.separation(c2).arcmin
	size = np.array([15 if s <= 2.8 else 5 for s in sep])
	mag_in, period_in, size_in = mag[flag], period[flag], size[flag]
	mag_out, period_out, size_out = mag[~flag], period[~flag], size[~flag]
	color = ['red' if p >= 0.4 else 'blue' for p in period_in]

	plt.figure(figsize=(7.1,2.5)) # x_old = 6, y = 2.65
	gs = gridspec.GridSpec(1, 2, width_ratios=[14,1])
	ax1, ax2 = plt.subplot(gs[0]), plt.subplot(gs[1])
	ax1.tick_params(axis='both', which='major', labelsize=13.5)
	ax1.set_xlabel('Period (days)', fontsize=13.5)
	if filt == 'Ks': filt = r'K_{\rm{S}}'
	ax1.set_ylabel(r'$\langle %s\rangle$'%filt, fontsize=13.5)
	auto_ticks(ax1)

	# if NGC[0:3] == 'NGC':
	# 	ax1.set_title(r'%s$\,$%s: $\langle %s\rangle$ vs. Period'%(NGC[:3],NGC[3:],filt))
	# else: ax1.set_title(r'%s$\,$%s: $\langle %s\rangle$ vs. Period'%(NGC[:-1],NGC[-1],filt))
	ax1.text(0.77, mean-2+0.76, r'%s$\,$%s'%(NGC[:3],NGC[3:]), \
			color='black',fontsize=17, family='sans-serif')	# max(period)+0.05-0.135
	
	#ax1.set_xlim(min(period)-0.05,max(period)+0.05)	# normal x-range
	ax1.set_xlim(0.21, 0.94)							# panels
	ax1.set_ylim(mean-2.0, mean+2.0)
	if not memb:
		for x, y, s, c in zip(period_in, mag_in, size_in, color):
			ax1.scatter(x, y, s=s, zorder=5, color=c)
		ax1.scatter(period_out, mag_out, s=size_out, zorder=4, color='silver')
	else:
		cividisBig = cm.get_cmap('plasma_r', 512)
		newcmp = ListedColormap(cividisBig(np.linspace(0.1, 0.95, 256)))
		print(len(mag), len(mag[flag]))
		pcm = ax1.scatter(period[flag], mag[flag], s=28, c=pmemb[flag], cmap=\
				newcmp, zorder=4)#, norm=colors.PowerNorm(gamma=1))	# 1/8
		if len(mag[~flag]) != 0: 
			pcm2 = ax1.scatter(period[~flag], mag[~flag], s=28, c=pmemb[~flag], \
					cmap=newcmp, zorder=4)#, norm=colors.PowerNorm(gamma=1/4))
			pcm2.set_facecolor('none')
		#ax1.scatter(period, mag, s=15, c='silver', zorder=5)
	
	print(alpha)
	label = r'%.3f$\,\pm\,$%.3f'%(mean,err)
	ax1.axhline(mean, color='k', ls='--', lw=1.2, zorder=1, label=label)
	ax1.axhline(mean+err, color='gray', ls=':', lw=1, zorder=2)
	ax1.axhline(mean-err, color='gray', ls=':', lw=1, zorder=2)
	ax1.axhline(med+alpha*std, color='gray', ls=':', lw=1, zorder=2)
	ax1.axhline(med-alpha*std, color='gray', ls=':', lw=1, zorder=2)
	ax1.legend(loc=2, fontsize=13.5)
	ax1.invert_yaxis()
	ax2.axis('off')
	ax2.set_ylim(mean-2.0,mean+2.0)
	ax2.invert_yaxis()
	ax2.hist(mag,bins=30,range=[mean-3.0,mean+3.0],histtype='step',orientation='horizontal',density=1)
	plt.tight_layout()
	plt.subplots_adjust(wspace=0.05)
	plt.show()
	
	if not memb: return list_med

def plot_CMDs(memb, tel, gaia_filt, rrl_gaiac, ir_filt, rrl_irc, p_member=[]):	# 29.mar.21: OK!
	'''
	Function to plot 2 CMDs (VVV and Gaia), also showing the RR Lyrae.
	'''

	if memb:
		Gmag, period = rrl_gaiac['corr_phot_g_mean_mag'], rrl_gaiac['o_Period']
		Gmag_lim = Gmag[(Gmag > 15) & (p_member > 1e-02)]
		period_lim = period[(Gmag > 15) & (p_member > 1e-02)]
		p_member_lim = p_member[(Gmag > 15) & (p_member > 1e-02)]
		#print (np.sort(p_member))
		#print (np.sort(p_member_lim))
		mean = np.mean(Gmag_lim)
		err = np.std(Gmag_lim)/(len(Gmag_lim)**(0.5))
		print(f'  \033[1m4) <Gmag> vs. Period:\033[0m {mean:.3f} +/- {err:.3f}')
		print (len(Gmag_lim), len(Gmag))
		plt.figure(figsize=(7,3)); ax = plt.gca()
		pcm = ax.scatter(period_lim, Gmag_lim, s=15, c=p_member_lim, cmap=\
						'plasma_r', norm=colors.PowerNorm(gamma=1/8), label=\
						'RR Lyrae', zorder=3)
		ax.axhline(mean, ls='--', color='k', label=r'%.2f$\pm$%.2f'%(mean,err), zorder=1)
		ax.axhline(mean+err, color='gray', ls=':', lw=0.8, zorder=2)
		ax.axhline(mean-err, color='gray', ls=':', lw=0.8, zorder=2)
		plt.xlim(0.2,0.75)
		plt.ylim(mean-0.40, mean+0.40)
		plt.gca().invert_yaxis()
		plt.ylabel(r'$\langle$Gmag$\rangle$')
		plt.xlabel('Period (days)')
		plt.legend(loc=2)
		plt.tight_layout()
		plt.show()
	
	rad = 15 if NGC == 'NGC6266' else 10
	n = ['Jmag3','Ksmag3'] if tel == 'vvv' else ['Jmag','Kmag']
	if not memb or reHeader: print(f'\033[1m  4) CMD Ks vs. J-Ks:\033[0m %d stars %s'%\
				(len(ir_filt[n[1]]),colored(u' \033[1m\u2713\033[0m','green')))
	if not memb: fig,(ax1,ax2) = plt.subplots(1,2,figsize=(7.5,6.2))
	else: fig,(ax1,ax2) = plt.subplots(1,2,figsize=(8,7))
	ax1.set_xlabel(r'$J-K_{\rm{S}}$')
	ax1.set_ylabel(r'$K_{\rm{S}}$')
	ax2.set_xlabel(r'$\rm{BP}-\rm{RP}$')
	ax2.set_ylabel(r'Gmag')
	ax1.scatter(ir_filt[n[0]]-ir_filt[n[1]], ir_filt[n[1]], s=0.1, color='gray',\
				marker='.', label=r'%d$\,$arcmin'%rad)
	ax2.scatter(gaia_filt['bp_rp'], gaia_filt['corr_phot_g_mean_mag'], s=0.1,\
				color='gray', marker='.', label=r'%d$\,$arcmin'%rad)
	
	if not memb: 
		ax1.scatter(rrl_irc[n[0]]-rrl_irc[n[1]], rrl_irc[n[1]], s=7, color='red',\
					marker='o', label='RR Lyrae')
		ax2.scatter(rrl_gaiac['bp_rp'], rrl_gaiac['corr_phot_g_mean_mag'], s=7,\
					color='red', marker='o', label='RR Lyrae')
	else:
		ax1.scatter(rrl_irc[n[0]]-rrl_irc[n[1]], rrl_irc[n[1]], s=7, color='red',\
					marker='o',label='RR Lyrae')
		#pcm = ax1.scatter(Jmag-Ksmag,Ksmag,s=7,c=p_member_K,cmap='plasma_r',norm=\
		#				colors.PowerNorm(gamma=1/8),label='RR Lyrae',zorder=3)
		pcm2 = ax2.scatter(rrl_gaiac['bp_rp'], rrl_gaiac['corr_phot_g_mean_mag'],\
						s=7, c=p_member, cmap='plasma_r', label='RR Lyrae', norm=\
						colors.PowerNorm(gamma=1/8), zorder=3)
		divider = make_axes_locatable(ax1)
		cax = divider.append_axes('top', size='2%', pad=0.3)
		cbar = plt.colorbar(pcm,cax=cax, orientation='horizontal',\
							ticks=[0.0,0.001,0.01,0.1,1.0])
		cax.xaxis.set_ticks_position('top')
		cbar.set_label('Membership', fontsize=10)
		divider2 = make_axes_locatable(ax2)
		cax2 = divider2.append_axes('top', size='2%', pad=0.3)
		cbar2 = plt.colorbar(pcm2,cax=cax2, orientation='horizontal',\
							ticks=[0.0,0.001,0.01,0.1,1.0])
		cax2.xaxis.set_ticks_position('top')
		cbar2.set_label('Membership', fontsize=10)
	
	m1, c1 = np.array(ir_filt[n[1]]), np.array(ir_filt[n[0]] - ir_filt[n[1]])
	m2, c2 = np.array(gaia_filt['corr_phot_g_mean_mag']), np.array(gaia_filt['bp_rp'])
	m1, c1 = m1[~np.isnan(m1)], c1[~np.isnan(c1)]
	m2, c2 = m2[~np.isnan(m2)], c2[~np.isnan(c2)]
	mad1 = scipy.stats.median_abs_deviation(c1)
	mad2 = scipy.stats.median_abs_deviation(c2)
	ax1.set_xlim(np.median(c1)-6.5*mad1, np.median(c1)+6.5*mad1)
	ax1.set_ylim(18.5, m1.min()-0.5)
	ax2.set_xlim(np.median(c2)-6.5*mad2, np.median(c2)+6.5*mad2)
	ax2.set_ylim(21, m1.min()+1.0)
	if not memb:
		ax1.set_title(r'%s$\,$%s: $K_{\rm{S}}$ vs. $J-K_{\rm{S}}$ (%s)'%\
						(NGC[:3], NGC[3:], tel.upper()))
		ax2.set_title(r'%s$\,$%s: Gmag vs. $\rm{BP}-\rm{RP}$ (Gaia)'%\
						(NGC[:3], NGC[3:]))
	auto_ticks(ax1); auto_ticks(ax2)
	ax1.legend(loc=2); ax2.legend(loc=2)
	plt.tight_layout()
	plt.tight_layout()
	plt.show()

	# Scatter plot with density plot behind (https://stackoverflow.com/questions/
	# 20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib)
	'''from scipy.stats import gaussian_kde
	x = Jmag_all-Ksmag_all
	y = Ksmag_all
	xy = np.vstack([x,y])
	z = gaussian_kde(xy)(xy)
	idx = z.argsort()
	x, y, z = x[idx], y[idx], z[idx]
	ax1.scatter(x,y,c=z,s=0.5,edgecolor='')'''

def VPD_select(memb, gaia_filt, rrl_gaiac, c1, tel, pars={}):	# 31.mar.21: OK!
	'''
	This function plots a vector-point diagram, marking the stars within certain
	radius, to apply the Gaussian Mixed Models. Thus the code gets proper-motion
	in right ascension and declination, and the sigmas and errors.
	####
	Need to insert pm/e_pm < 20??? If using Arenou, it's not necessary!!!
	'''
	def GMM_1d(data, k=2, xmin=-15, xmax=15):
		gmm = GaussianMixture(n_components=k, covariance_type='full',tol=1e-10,\
				 n_init=5, max_iter=1000)
		gmm = gmm.fit(np.expand_dims(data, 1))
		gmm_x = np.linspace(xmin, xmax, 1000)
		gmm_y = np.exp(gmm.score_samples(gmm_x.reshape(-1, 1)))
		#score = gmm.score(gmm_x)	# average Likelihood
		#print (score)
		
		# Getting the mean values, weights, covariance
		if k == 2:
			m1, m2 = gmm.means_
			c1, c2 = gmm.covariances_
			w1, w2 = gmm.weights_
		elif k == 3:
			m1, m2, m3 = gmm.means_
			c1, c2, c3 = gmm.covariances_
			w1, w2, w3 = gmm.weights_
		c1 = np.array((c1))
		print ('1D: means, covariances and weights')
		print (m1,m2)
		print (c1,c2)
		print (w1,w2)
		
		# Contribution of each Gaussian multiplied by the weight
		g1 = np.array(scipy.stats.norm.pdf(gmm_x, m1, c1**0.5)*w1)
		g2 = np.array(scipy.stats.norm.pdf(gmm_x, m2, c2**0.5)*w2)
		g1 = g1.reshape(len(gmm_x))
		g2 = g2.reshape(len(gmm_x))
		g  = [g1, g2]
		if k == 3:
			g3 = np.array(scipy.stats.norm.pdf(gmm_x, m3, c3**0.5)*w3)
			g3 = g3.reshape(len(gmm_x))
			g.append(g3)
		
		global weight_c, weight_f
		ind_min_cov = np.argmin(gmm.covariances_)
		pm = gmm.means_[ind_min_cov]
		sig = np.sqrt(gmm.covariances_[ind_min_cov])
		weight_c = gmm.weights_[ind_min_cov]
		gamma = scipy.stats.norm.pdf(gmm_x,pm,sig) * weight_c
		print (gamma.shape, gmm_y.shape)
		S = sum(np.array(gamma)[0]/gmm_y)
		cov = (gmm.covariances_[ind_min_cov])
		error = (cov/S)**0.5
		print (S, error)
		
		ind_max_cov = np.argmax(gmm.covariances_)
		pm_field = gmm.means_[ind_max_cov]
		sig_field = np.sqrt(gmm.covariances_[ind_max_cov])
		weight_f = gmm.weights_[ind_max_cov]
		gammaf = scipy.stats.norm.pdf(gmm_x, pm_field, sig_field) * weight_f
		Sf = sum(np.array(gammaf)[0]/gmm_y)
		cov_field = (gmm.covariances_[ind_max_cov])
		error_field = (cov_field/Sf)**0.5
		print (Sf, error_field)
		#print (weight_c, weight_f)		# check the adequate weight (average?)

		pms_1d = (pm, sig, error, pm_field, sig_field, error_field)
		return g, gmm_x, gmm_y, pms_1d
	
	def GMM_2d(dataX, dataY, data_out, k=2, xmin=-15, xmax=15):	### MISSING: fix err_field
		data = np.c_[dataX, dataY]
		#Y = np.random.randint(xmin, xmax, size=(1, 2))
		GMM = GaussianMixture(n_components=k, covariance_type='full',tol=1e-15,\
				 n_init=10, max_iter=1000).fit(data)
		print('\nConverged:',GMM.converged_)

		#gmm_x_1d = np.linspace(xmin, xmax, 1000)
		gmm_x = np.linspace((xmin, xmin),(xmax, xmax), 1000)
		x,y = np.meshgrid(np.sort(gmm_x[:,0]),np.sort(gmm_x[:,1]))
		XY = np.array([x.ravel(),y.ravel()]).T
		gmm_y = np.exp(GMM.score_samples(gmm_x))	# must be 2d !!!
		print (data.shape, gmm_x.shape, XY.shape, gmm_y.shape)

		if k == 2:
			m1, m2 = GMM.means_
			c1, c2 = GMM.covariances_
			w1, w2 = GMM.weights_
		elif k == 3:
			m1, m2, m3 = GMM.means_
			c1, c2, c3 = GMM.covariances_
			w1, w2, w3 = GMM.weights_
		
		#.calculating cluster PM and sigma (RA & DEC)
		ind_min_cov = np.argmin(GMM.covariances_[:,0][:,0])
		pm_2d = GMM.means_[ind_min_cov]
		RApm_2d, DEpm_2d = pm_2d
		covm_2d = GMM.covariances_[ind_min_cov]
		pmRAsig_2d, pmDEsig_2d = np.sqrt(covm_2d.diagonal())
		wc_2d = GMM.weights_[ind_min_cov]
		multi_normal = multivariate_normal(mean=pm_2d,cov=covm_2d)
		gamma = multi_normal.pdf(XY).reshape(len(gmm_x),len(gmm_x)) * wc_2d
		# !!! gamma.shape = (1000,1000) ; gmm_y.shape = (1000,) !!!		or 1000*1000???
		# STUDY: https://scikit-learn.org/stable/modules/generated/sklearn.mixture.GaussianMixture.html#sklearn.mixture.GaussianMixture.score_samples
		#https://scikit-learn.org/stable/auto_examples/mixture/plot_gmm_pdf.html
		sum_list = [sum(np.array(gamma[i])/gmm_y) for i in range(len(gmm_y))]
		S_2d = np.mean(sum_list)
		#print (gamma.shape, gmm_y.shape)
		#S_2d = sum(np.array(gamma)/gmm_y)
		errRA_2d, errDE_2d = np.array([pmRAsig_2d, pmDEsig_2d]) / (S_2d**0.5)
		print ('\nCluster :::')
		print ('pmRA, pmDE =', RApm_2d, DEpm_2d)
		print ('pmRAsig, pmDEsig =', pmRAsig_2d, pmDEsig_2d)
		print ('weight =', wc_2d)
		print (S_2d)
		print ('errRA, errDE =', errRA_2d, errDE_2d)
		
		#.calculating field PM and sigma (RA & DEC)
		ind_max_cov = np.argmax(GMM.covariances_[:,0][:,0])
		pm_field_2d = GMM.means_[ind_max_cov]
		RApm_field_2d, DEpm_field_2d = pm_field_2d
		covm_field_2d = GMM.covariances_[ind_max_cov]
		pmRAsig_field_2d, pmDEsig_field_2d = np.sqrt(covm_field_2d.diagonal())
		wf_2d = GMM.weights_[ind_max_cov]
		multi_normalf = multivariate_normal(mean=pm_field_2d,cov=covm_field_2d)
		gammaf = multi_normalf.pdf(XY).reshape(len(gmm_x),len(gmm_x)) * wf_2d
		sum_listf = [sum(np.array(gammaf[i])/gmm_y) for i in range(len(gmm_y))]
		Sf_2d = np.mean(sum_listf)
		errRA_field_2d, errDE_field_2d = np.array([pmRAsig_field_2d, pmDEsig_field_2d]) / (Sf_2d**0.5)
		print ('\nField :::')
		print ('pmRA, pmDE =', RApm_field_2d, DEpm_field_2d)
		print ('pmRAsig, pmDEsig =', pmRAsig_field_2d, pmDEsig_field_2d)
		print ('weight =', wf_2d)
		print (covm_field_2d)
		print (Sf_2d)
		print ('errRA, errDE =', errRA_field_2d, errDE_field_2d) ### WRONG: FIX IT!!! 

		if not memb: plt.figure(figsize=(5,5))
		else: fig = plt.figure(figsize=(6.2,5.2))
		ax = plt.gca()
		for m, c in zip(GMM.means_, GMM.covariances_):
			multi_normal = multivariate_normal(mean=m,cov=c)
			ax.contour(np.sort(gmm_x[:,0]), np.sort(gmm_x[:,1]),\
						multi_normal.pdf(XY).reshape(len(gmm_x),len(gmm_x)), colors=\
						'black',alpha=0.7,linewidths=1,zorder=4,levels=8,
						linestyles='dashed')	# levels=8
			ax.scatter(m[0],m[1],c='black',marker='+',zorder=10,s=150)
		plt.scatter(data_out[0], data_out[1], s=0.8, marker='x', c='gainsboro', label= \
				r'Gaia EDR3 (%d$^{\prime}$)'%rad,zorder=1)
		plt.scatter(data[:, 0], data[:, 1], s=1, marker='x', c='c',zorder=2,\
					label=r'r $\leq$ %.1f$^{\prime}$'%(c_rad*60))
		ax.set_xlabel('pmRA (mas/yr)', fontsize=12.5)
		ax.set_ylabel('pmDE (mas/yr)', fontsize=12.5)
		ax.tick_params(axis='both', which='major', labelsize=13.5)
		#ax.update({'xlabel': , 'ylabel': 'pmDE (mas/yr)'})
		#ax.update({'xlim': (-15,15), 'ylim': (-15,15)})
		if not memb:
			#plt.title(r'VPD: %s$\,$%s (Gaia vs. %s)'%(NGC[:3],NGC[3:],tel))
			plt.scatter(pmRA,pmDE,s=8,marker='o',c='red',label='RR Lyrae',zorder=3)

		# https://matplotlib.org/3.1.0/tutorials/colors/colormap-manipulation.html
		# from matplotlib import cm
		# from matplotlib.colors import ListedColormap, LinearSegmentedColormap
		# # cividis = cm.get_cmap('gist_earth_r', 256)
		# newcolors = cividis(np.linspace(0, 1, 256))
		# pink = np.array([248/256, 24/256, 148/256, 1])
		# newcolors[:25, :] = pink
		# newcmp = ListedColormap(newcolors)

		cividisBig = cm.get_cmap('plasma_r', 512)
		newcmp = ListedColormap(cividisBig(np.linspace(0.1, 0.95, 256)))
		# plot_examples([cividis, newcmp])

		# good options: plasma_r (original), viridis, cividis, gist_earth
		# viridis: 0.15 to 0.95 is very good! (o azul confunde)
		if memb:
			# 'plasma_r'
			# pcm = plt.scatter(pmRA, pmDE, s=25, c=p_member, cmap='cividis_r',\
			# 					label='RR Lyrae', zorder=3)
			pcm = plt.scatter(pmRA, pmDE, s=32, c=p_member, cmap=newcmp,\
								label='RR Lyrae', zorder=3)
			divider = make_axes_locatable(ax)
			cax = divider.append_axes('right', size='3%', pad=0.3)
			cbar = fig.colorbar(pcm, cax=cax, orientation='vertical', ticks=\
								[0.0,0.2,0.4,0.6,0.8,1.0])
			plt.clim(0,1)
			cax.xaxis.set_ticks_position('top')
			cbar.set_label('Membership', fontsize=11)
			plt.tick_params(axis='both', which='major', labelsize=11)
		
		print(NGC)
		lims = {('NGC6266','NGC6441','NGC6626'): [1.75,0.49,1.43],
				('NGC6638'): [1.25,0.34,1.02], ('NGC6642'): [1.05,0.29,0.85],
				('NGC6717'): [0.85,0.05,0.52]}
		for (i, item) in enumerate([*lims]):
			if NGC in item: break
		lim = lims[[*lims][i]]
		
		if NGC != 'NGC6717':
			ax.set_xlim(RApm_field_2d-lim[0]*pmRAsig_field_2d,RApm_field_2d+lim[0]*pmRAsig_field_2d)
			ax.set_ylim(DEpm_field_2d-lim[0]*pmDEsig_field_2d,DEpm_field_2d+lim[0]*pmDEsig_field_2d)
		else:
			ax.set_xlim(-4.23, 0.96)
			ax.set_ylim(-6.84, -1.39)
		text = r'%s$\,%s$'%(NGC[:3],NGC[3:])
		ax.text(RApm_field_2d+lim[1]*pmRAsig_field_2d, DEpm_field_2d+lim[2]*pmDEsig_field_2d,\
				text, fontsize=20, color='black', family='sans-serif')
		ax.legend(loc=2,fontsize=12)
		plt.tight_layout()
		plt.tight_layout()
		plt.show()

		#return GMM
		# still missing the corrected err_pm...
		pms = (RApm_2d, DEpm_2d, pmRAsig_2d, pmDEsig_2d)
		f_pms = (RApm_field_2d, DEpm_field_2d, pmRAsig_field_2d, pmDEsig_field_2d)

		return pms, f_pms
			
	pmRA, pmDE = rrl_gaiac['corr_pmra'], rrl_gaiac['corr_pmdec']
	pmRA_all, e_pmRA_all = gaia_filt['corr_pmra'], gaia_filt['pmra_error']
	pmDE_all, e_pmDE_all = gaia_filt['corr_pmdec'], gaia_filt['pmdec_error']
	
	#global RApm,pmRAsig,RApm_f,pmRAsig_f,DEpm, pmDEsig, DEpm_f, pmDEsig_f
	if not memb: 
		c_rad = float(input('\033[1m  5) VPD:\033[0m PM-sel. with Rad(deg) = '))
		plt.figure(figsize=(5,5))
		plt.title(r'VPD: %s$\,$%s (Gaia vs. %s)'%(NGC[:3], NGC[3:], tel.upper()))
	else:
		c_rad = pars['c_rad']
		p_member = pars['memb_bellini']
		plt.figure(figsize=(4.75,5))
	ax = plt.gca()
	auto_ticks(ax)
	ax.update({'xlabel': 'pmRA (mas/yr)', 'ylabel': 'pmDE (mas/yr)'})
	#lim = (-15,15) if NGC != 'NGC6266' else (-15,15)
	ax.update({'xlim': (-15,15), 'ylim': (-15,15)})

	#.Ang-sep of all Gaia stars (original is on SOAR version...)
	c_gaia = SkyCoord(ra=gaia_filt['ra'], dec=gaia_filt['dec'])
	sep_all = (c1.separation(c_gaia)).value	# deg
	# b_all = np.array([1 if sep <= c_rad else 0 for sep in sep_all])
	# out = np.c_[gaia_filt['ra'], gaia_filt['dec'], sep_all, b_all, \
	# 			gaia_filt['corr_phot_g_mean_mag'], gaia_filt['bp_rp'],\
	# 			gaia_filt['corr_pmra'], gaia_filt['corr_pmdec']]
	# lab = ['RA_icrs','DEC_icrs','rad','flag','Gmag','BP-RP','pmRA','pmDE']
	# dtypes = np.repeat('f4', out.shape[1])
	# f = Table(out, names=lab, dtype=dtypes)
	# f.write(f'{path}/{NGC}/Aux_Cats/{NGC}_AngSep.fits',format='fits', overwrite=True)

	f_pm = (~np.isnan(pmRA_all)) & (~np.isnan(pmDE_all))
	f_cen = (sep_all <= c_rad) & (~np.isnan(gaia_filt['corr_phot_g_mean_mag']))
	# f_cen = (sep_all <= c_rad) & (Gmag_all < 20.6) & (Gmag_all > -99.99)
	pmRA_cen, e_pmRA_cen = pmRA_all[f_cen], e_pmRA_all[f_cen]
	pmDE_cen, e_pmDE_cen = pmDE_all[f_cen], e_pmDE_all[f_cen]
	pmRA_out, pmDE_out = pmRA_all[f_pm & ~f_cen], pmDE_all[f_pm & ~f_cen]
	data_out = [pmRA_out, pmDE_out]
	n_pmGaia = len(sep_all[f_pm])

	rad = 10 if NGC != 'NGC6266' else 15
	plt.scatter(pmRA_out, pmDE_out, s=0.8, marker='x', c='silver', label= \
				r'Gaia DR2 (%d$^{\prime}$)'%rad, zorder=1)
	plt.scatter(pmRA_cen, pmDE_cen, s=0.5, marker='o', c='c', label= \
				r'Center$\sim$%.1f$^{\prime}$'%(c_rad*60), zorder=2)
	if not memb:
		plt.scatter(pmRA,pmDE,s=10,marker='o',c='red',label='RR Lyrae',zorder=3)
	else:
		pcm = plt.scatter(pmRA, pmDE, s=10, c=p_member, cmap='plasma_r', norm=\
						colors.PowerNorm(gamma=1/8), label='RR Lyrae', zorder=3)
		divider = make_axes_locatable(ax)
		cax = divider.append_axes('top', size='2%', pad=0.3)
		cbar = plt.colorbar(pcm, cax=cax, orientation='horizontal', ticks=\
							[0.0,0.001,0.01,0.1,1.0])
		cax.xaxis.set_ticks_position('top')
		cbar.set_label('Membership', fontsize=10)
	ax.legend(loc=2,fontsize=9.8)
	plt.tight_layout()
	plt.show()

	m3 = (abs(pmDE_cen) <= 15.0) & (e_pmRA_cen < 0.20) & \
		 (abs(pmRA_cen) <= 15.0) & (e_pmDE_cen < 0.20) # quality-flag
	pmRA_clean, pmDE_clean = pmRA_cen[m3], pmDE_cen[m3]

	if memb:	# -->> WRONG!?! pm/e_pm < 0.20
		pms_2d, f_pms_2d = GMM_2d(pmRA_clean, pmDE_clean, data_out)
	
	n_pmGaia_RR = len(pmRA[(~np.isnan(pmRA)) & (~np.isnan(pmDE))])
	if not memb:
		print ('%s%d stars with good measurement of pmRA and pmDE%s'%(' '*6, \
				n_pmGaia, colored(u' \033[1m\u2713\033[0m','green')))
		print ('%s%d/%d RRLs with good measurement'%(' '*6,n_pmGaia_RR,len(rrl_gaiac)))

		# pmRAsig = Gaussian dispersion ; pmRAerr = uncert. in PM determination
		GMM_k = int(input('   > GMM components (2 or 3) = '))
		RAg, RAgmm_x, RAgmm_y, pms_RA_1d = GMM_1d(pmRA_clean, GMM_k)
		DEg, DEgmm_x, DEgmm_y, pms_DE_1d = GMM_1d(pmDE_clean, GMM_k)
		
		RApm, pmRAsig, pmRAerr, RApm_f, pmRAsig_f, pmRAerr_f = pms_RA_1d
		DEpm, pmDEsig, pmDEerr, DEpm_f, pmDEsig_f, pmDEerr_f = pms_DE_1d

		pms_2d, f_pms_2d = GMM_2d(pmRA_clean, pmDE_clean, data_out, GMM_k)
		RApm_2d, DEpm_2d, pmRAsig_2d, pmDEsig_2d = pms_2d
		RApm_field_2d, DEpm_field_2d, pmRAsig_field_2d, pmDEsig_field_2d = f_pms_2d
		
		'''print (RAgmm_x[RAg[0] == max(RAg[0])], max(RAg[0]))
		print (RAgmm_x[RAg[1] == max(RAg[1])], max(RAg[1]))
		print (DEgmm_x[DEg[0] == max(DEg[0])], max(DEg[0]))
		print (DEgmm_x[DEg[1] == max(DEg[1])], max(DEg[1]))
		input()'''
		print ('   > RA-cluster:', RApm, pmRAsig, pmRAerr)
		print ('   > DE-cluster:', DEpm, pmDEsig, pmDEerr)
		print ('   > RA-field:', RApm_f, pmRAsig_f, pmRAerr_f)
		print ('   > DE-field:', DEpm_f, pmDEsig_f, pmDEerr_f)
		
		fig, (ax1, ax2) = plt.subplots(1,2,figsize=(8,4),sharex=True)
		ax1.plot(RAgmm_x, RAgmm_y, color='crimson', lw=2, label='GMM')
		ax1.plot(RAgmm_x, RAg[0], ':')
		ax1.plot(RAgmm_x, RAg[1], ':')
		if GMM_k == 3: ax1.plot(RAgmm_x, RAg[2],':')
		ax1.hist(pmRA_clean, bins=20, density=True, histtype='step')
		ax1.update({'xlabel': 'pmRA (mas/yr)', 'ylabel': 'Frequency'})
		ax1.set_title(r'$\langle pmRA\rangle$ = %2.3f; $\sigma_{pmRA}$ = %2.3f'\
						%(RApm[0],pmRAsig[0][0]))
		ax1.set_xlim(-15,15)
		ax2.plot(DEgmm_x, DEgmm_y, color='crimson', lw=2, label='GMM')
		ax2.plot(DEgmm_x,DEg[0],':')
		ax2.plot(DEgmm_x,DEg[1],':')
		if GMM_k == 3: ax2.plot(DEgmm_x,DEg[2],':')
		ax2.hist(pmDE_clean, bins=20, density=True, histtype='step')
		ax2.update({'xlabel': 'pmDE (mas/yr)', 'ylabel': 'Frequency'})
		ax2.set_title(r'$\langle pmDE\rangle$ = %2.3f; $\sigma_{pmDE}$ = %2.3f'\
						%(DEpm[0],pmDEsig[0][0]))
		ax2.set_xlim(-15,15)
		plt.tight_layout()
		plt.show()
	
		return (c_rad, pms_RA_1d, pms_DE_1d, pms_2d, f_pms_2d)
	return c_rad
	
def membership(gaia_filt, rrl_gaiac, c1, pars, alpha):	# 01.abr.21: OK!
	'''
	Calculates the final membership.
	'''
	#global p_member, p_member_field, p_member_norm, p_member_norm2
	#global memb_bellini
	#p_member_norm, p_member_norm2 = [],[]
	print ()
	#print (len(Vmag),len(Imag),len(Ksmag))
	#print (len(pmRA),len(pmRA[(pmRA != -99.99) & (pmDE != -99.99)]))
	#print (len(pmDE),len(pmDE[(pmRA != -99.99) & (pmDE != -99.99)]))

	ra, dec = np.array(rrl_gaiac['ra']), np.array(rrl_gaiac['dec'])
	pmRA, pmDE = np.array(rrl_gaiac['corr_pmra']), np.array(rrl_gaiac['corr_pmdec'])
	e_pmRA, e_pmDE = np.array(rrl_gaiac['pmra_error']), np.array(rrl_gaiac['pmdec_error'])
	pmRA_all, e_pmRA_all = gaia_filt['corr_pmra'], gaia_filt['pmra_error']
	pmDE_all, e_pmDE_all = gaia_filt['corr_pmdec'], gaia_filt['pmdec_error']
	
	# ask if 1d or 2d gaussian?
	pms_2d, f_pms_2d = pars['pms_1d_2d'][-2], pars['pms_1d_2d'][-1]
	RApm, DEpm, pmRAsig, pmDEsig = pms_2d
	RApm_f, DEpm_f, pmRAsig_f, pmDEsig_f = f_pms_2d
	
	import warnings
	warnings.filterwarnings("ignore", category=RuntimeWarning)
	
	r_h = {'NGC6266':0.92,'NGC6352':2.05,'NGC6626':1.97,'NGC6638':0.51,
			'NGC6642':0.73,'NGC6717':0.68}
	c2 = SkyCoord(ra*u.si.deg, dec*u.si.deg, frame='icrs')
	sep_RRgaia = c1.separation(c2).arcmin

	#.Calculating the old p_member:
	chi_pmRA = np.exp(-(pmRA - RApm)**2/(2*pmRAsig**2))#/max(np.exp(-(pmRA - RApm)**2/(2*pmRAsig**2)))
	chi_pmDE = np.exp(-(pmDE - DEpm)**2/(2*pmDEsig**2))#/max(np.exp(-(pmDE - DEpm)**2/(2*pmDEsig**2)))
	#chi_Radi = np.exp(-(sep_RRgaia**2)/(2*(r_h[NGC]**2)))/max(np.exp(-(sep_RRgaia**2)/(2*(r_h[NGC]**2))))
	p_member = chi_pmRA * chi_pmDE

	#.Calculating the old p_member_field
	chi_pmRA_field = np.exp(-(pmRA - RApm_f)**2/(2*pmRAsig_f**2))#/max(np.exp(-(pmRA - RApm)**2/(2*pmRAsig**2)))
	chi_pmDE_field = np.exp(-(pmDE - DEpm_f)**2/(2*pmDEsig_f**2))#/max(np.exp(-(pmDE - DEpm)**2/(2*pmDEsig**2)))
	p_member_field = chi_pmRA_field * chi_pmDE_field

	#.Calculating the new memb_bellini (from Bellini et al. 2009)
	#if pmRA[i] != -99.99 and pmDE[i] != -99.99 and e_pmRA[i] != -99.99 and e_pmDE[i] != -99.99:
	gamma_bellini = ((pmRA-RApm_f)*(pmDE-DEpm_f))/(pmRAsig_f*pmDEsig_f)
	phic_bellini = (np.exp(-0.5*(((pmRA-RApm)**2/(pmRAsig**2+e_pmRA**2))+ \
					(((pmDE-DEpm)**2/(pmDEsig**2+e_pmDE**2))))))/(2*np.pi* \
					((pmRAsig**2 + e_pmRA**2)**(1/2))*((pmDEsig**2 + e_pmDE**2)**(1/2)))
	phif_bellini = (np.exp((-1/(2*(1-gamma_bellini**2)))*(((pmRA-RApm_f)**2/(pmRAsig_f**2+e_pmRA**2))-((2*gamma_bellini*(pmRA-RApm_f)*(pmDE-DEpm_f))/(((pmRAsig_f**2+e_pmRA**2)**(1/2))*((pmDEsig_f**2+e_pmDE**2)**(1/2))))+(((pmDE-DEpm_f)**2)/(DEpm_f**2+e_pmDE**2)))))/(2*np.pi*((1-gamma_bellini**2)**(1/2))*((pmRAsig_f**2 + e_pmRA**2)**(1/2))*((pmDEsig_f**2 + e_pmDE**2)**(1/2)))
	phi_bellini = (weight_c*phic_bellini) + (weight_f*phif_bellini)
	memb_bellini = weight_c * (phic_bellini/phi_bellini)
	membf_bellini = weight_f * (phif_bellini/phi_bellini)

	'''
	for i in range(len(ra)):
		a = 1
		#chi_pmRA = np.exp(-(pmRA[i] - RApm[0])**2/(2*pmRAsig[0][0]**2))/max(np.exp(-(pmRA - RApm[0])**2/(2*pmRAsig[0][0]**2)))
		#chi_pmDE = np.exp(-(pmDE[i] - DEpm[0])**2/(2*pmDEsig[0][0]**2))/max(np.exp(-(pmDE - DEpm[0])**2/(2*pmDEsig[0][0]**2)))
		chi_pmRA = np.exp(-(pmRA[i] - RApm[0])**2/(2*pmRAsig[0][0]**2))#/max(np.exp(-(pmRA - RApm[0])**2/(2*pmRAsig[0][0]**2)))
		chi_pmDE = np.exp(-(pmDE[i] - DEpm[0])**2/(2*pmDEsig[0][0]**2))#/max(np.exp(-(pmDE - DEpm[0])**2/(2*pmDEsig[0][0]**2)))
		#chi_Radi = np.exp(-(sep_RRgaia[i]**2)/(2*(r_h[NGC]**2)))/max(np.exp(-(sep_RRgaia**2)/(2*(r_h[NGC]**2))))

		### USE chi_Vmag??? NO!
		#if Vmag[i] != -99.99:
		#	chi_Vmag = np.exp(-(Vmag[i] - Vmag_mean)**2/(2*Vmag_std**2))/max(np.exp(-(Vmag - Vmag_mean)**2/(2*Vmag_std**2)))
		#else: chi_Vmag = 1e-5
		#if Ksmag[i] != -99.99:		# useless (for the cases with 2MASS, it's very few stars...)
		#	chi_Kmag = np.exp(-(Ksmag[i] - Ksmag_mean)**2/(2*Ksmag_std**2))/max(np.exp(-(Ksmag - Ksmag_mean)**2/(2*Ksmag_std**2)))
		#else: chi_Kmag = 1e-5
		#p_member.append(chi_pmRA * chi_pmDE * chi_Radi * chi_Vmag * chi_Kmag)	# Very old
		#p_member.append(chi_pmRA * chi_pmDE * chi_Vmag * chi_Kmag)				# Old
		p_member.append(chi_pmRA * chi_pmDE)
		
		chi_pmRA_field = np.exp(-(pmRA[i] - RApm_f[0])**2/(2*pmRAsig_f[0][0]**2))#/max(np.exp(-(pmRA - RApm[0])**2/(2*pmRAsig[0][0]**2)))
		chi_pmDE_field = np.exp(-(pmDE[i] - DEpm_f[0])**2/(2*pmDEsig_f[0][0]**2))#/max(np.exp(-(pmDE - DEpm[0])**2/(2*pmDEsig[0][0]**2)))
		p_member_field.append(chi_pmRA_field * chi_pmDE_field)
		
		# Old membership
		#p_member_norm.append(p_member[i]/(p_member[i]+p_member_field[i]))
		#p_member_norm2.append((chi_pmRA/(chi_pmRA+chi_pmRA_field))*(chi_pmDE/(chi_pmDE+chi_pmDE_field)))
		#p_member = np.array(p_member)#/max(p_member)
		#p_member_field = np.array(p_member_field)#/max(p_member_field)
		#p_member_norm2 = np.array(p_member_norm2)
		
		# pmRAall = np.array(pmRA_all)
		# pmDEall = np.array(pmDE_all)
		# e_pmRAall = np.array(e_pmRA_all)
		# e_pmDEall = np.array(e_pmDE_all)
		# mask1 = (pmRAall > -99.99) & (pmDEall > -99.99)
		# pmRAall, e_pmRAall = pmRAall[mask1], e_pmRAall[mask1]
		# pmDEall, e_pmDEall = pmDEall[mask1], e_pmDEall[mask1]
		
		# Bellini et al. (2009) ===>> Very good!
		# if pmRA[i] != -99.99 and pmDE[i] != -99.99 and e_pmRA[i] != -99.99 and e_pmDE[i] != -99.99:
		# 	gamma_bellini = ((pmRA[i]-RApm_f[0])*(pmDE[i]-DEpm_f[0]))/(pmRAsig_f[0][0]*pmDEsig_f[0][0])
		# 	phic_bellini = (np.exp(-0.5*(((pmRA[i]-RApm[0])**2/(pmRAsig[0][0]**2+e_pmRA[i]**2))+(((pmDE[i]-DEpm[0])**2/(pmDEsig[0][0]**2+e_pmDE[i]**2))))))/(2*np.pi*((pmRAsig[0][0]**2 + e_pmRA[i]**2)**(1/2))*((pmDEsig[0][0]**2 + e_pmDE[i]**2)**(1/2)))
		# 	phif_bellini = (np.exp((-1/(2*(1-gamma_bellini**2)))*(((pmRA[i]-RApm_f[0])**2/(pmRAsig_f[0][0]**2+e_pmRA[i]**2))-((2*gamma_bellini*(pmRA[i]-RApm_f[0])*(pmDE[i]-DEpm_f[0]))/(((pmRAsig_f[0][0]**2+e_pmRA[i]**2)**(1/2))*((pmDEsig_f[0][0]**2+e_pmDE[i]**2)**(1/2))))+(((pmDE[i]-DEpm_f[0])**2)/(DEpm_f[0]**2+e_pmDE[i]**2)))))/(2*np.pi*((1-gamma_bellini**2)**(1/2))*((pmRAsig_f[0][0]**2 + e_pmRA[i]**2)**(1/2))*((pmDEsig_f[0][0]**2 + e_pmDE[i]**2)**(1/2)))
		# 	phi_bellini = (weight_c*phic_bellini) + (weight_f*phif_bellini)
		# 	membership_bellini = weight_c*phic_bellini/phi_bellini
		# 	membf_bellini = weight_f*phif_bellini/phi_bellini
		# else: membership_bellini = 0.0
		# memb_bellini.append(membership_bellini)
		
		# Bellini et al. (2009) ==>>> Deu certo!!!
		# print ('Gamma-Bellini =',gamma_bellini)
		# print ('Phic-Bellini =',phic_bellini)
		# print ('Phif-Bellini =',phif_bellini)
		# print ('\033[1mMEMB-Bellini =',membership_bellini,'\033[0m')
		# print ('\033[1mMEMB =',membf_bellini+membership_bellini,'\033[0m')
		# input()
				
		#print (chi_pmRA , chi_pmDE, chi_pmRA_field, chi_pmDE_field)
		# print ('p_Member-RA:',chi_pmRA)
		# print ('p_Field-RA:',chi_pmRA_field)
		# print ('Sum-RA:',chi_pmRA+chi_pmRA_field)
		# print ('Norm-RA:', chi_pmRA/(chi_pmRA+chi_pmRA_field))
		# print ('---')
		# print ('p_Member-DE:',chi_pmDE)
		# print ('p_Field-DE:',chi_pmDE_field)
		# print ('Sum-DE:',chi_pmDE+chi_pmDE_field)
		# print ('Norm-DE:', chi_pmDE/(chi_pmDE+chi_pmDE_field))
		# print ('---')
		# print ('Norm-FINAL:', (p_member_norm2[i]))
		# print ()
		# print ()
		# input ()
	'''
	
	memb_bellini = np.nan_to_num(memb_bellini)
	membf_bellini = np.nan_to_num(membf_bellini)
	
	if not memb:
		print ('\033[1mmemb_bellini (sorted) =\033[0m')
		print (np.sort(memb_bellini)[::-1])

	plt.figure(figsize=(6,3.5))
	abc = []
	print ('\n   > %d RR Lyrae with membership = 0.00'%len(memb_bellini[memb_bellini == 0.0]))
	print ('   > %d RR Lyrae with membership > 0.10'%len(memb_bellini[memb_bellini > 0.1]))
	memb_limit = input('\n   > Insert the membership limit for %s [~0.5]: '%NGC)
	if memb_limit != '': memb_limit = float(memb_limit)
	else: memb_limit = 0.50
	print ('   > Number of stars with membership > %.2f: %d'%(memb_limit, \
			len(memb_bellini[memb_bellini > memb_limit])))
	input()
	
	for j in range(len(memb_bellini)):
		if memb_bellini[j] > memb_limit: abc.append(memb_bellini[j])
	y, x, _ = plt.hist(memb_bellini,bins=25,color='blue',label='All')
	plt.plot([memb_limit,memb_limit],[0,y.max()+4],ls='--',color='k',label='Limit = %.2f'%memb_limit)
	plt.xlabel('Membership',fontsize=11)
	plt.ylabel('Frequency',fontsize=11)
	plt.xlim(-0.03,1.03)
	plt.ylim(0, y.max()+2)
	plt.legend(loc='best')
	plt.tight_layout()
	plt.show()
	
	Vmag = np.nan_to_num(rrl_gaiac['mag'], nan=-99.99)
	Imag = np.nan_to_num(rrl_gaiac['Imag'], nan=-99.99)
	for i in range(len(memb_bellini)):
		print (memb_bellini[i], Vmag[i], Imag[i])
	#Vmag, Imag = np.nan_to_num(Vmag, nan=-99.99), Imag[~np.isnan(Imag)]
	V_mean, V_err, V_med, V_std = sig_clip(Vmag, 'V', alpha, 20, True)
	I_mean, I_err, I_med, I_std = sig_clip(Imag, 'I', alpha, 20, True)
	print([V_mean, V_err, V_med, V_std])
	print([I_mean, I_err, I_med, I_std])
	
	flag_sigV, flag_sigI = abs(Vmag-V_med) < 2*V_std, abs(Imag-I_med) < 2*I_std
	Vmag_clean, membV_clean = Vmag[flag_sigV], memb_bellini[flag_sigV]
	Imag_clean, membI_clean = Imag[flag_sigI], memb_bellini[flag_sigI]
	print ('Vmag:',len(Vmag[~np.isnan(Vmag)]),len(Vmag_clean),len(membV_clean))
	print ('Imag:',len(Vmag[~np.isnan(Imag)]),len(Imag_clean),len(membI_clean))
	
	#global Wmean_Vmag, err_Vmag, Wmean_Imag, err_Imag

	mean_Vmag  = np.mean(Vmag_clean)
	err_Vmag = np.std(Vmag_clean)/((len(membV_clean))**(1/2))
	Wmean_Vmag = (sum(Vmag_clean*membV_clean))/(sum(membV_clean))
	Werr_Vmag  = (sum(membV_clean))**(-0.5)
	print ()
	print ('\033[1mNormal <V> mean:\033[0m',mean_Vmag,'+/-',err_Vmag)
	print ('\033[1mWeighted <V> mean:\033[0m',Wmean_Vmag,'+/-',Werr_Vmag)
	
	mean_Imag  = np.mean(Imag_clean)
	err_Imag = np.std(Imag_clean)/((len(membI_clean))**(1/2))
	Wmean_Imag = (sum(Imag_clean*membI_clean))/(sum(membI_clean))
	Werr_Imag  = (sum(membI_clean))**(-0.5)
	print ()
	print ('\033[1mNormal <I> mean:\033[0m',mean_Imag,'+/-',err_Imag)
	print ('\033[1mWeighted <I> mean:\033[0m',Wmean_Imag,'+/-',Werr_Imag)

	weighted_V = [Wmean_Vmag, Werr_Vmag, mean_Vmag, err_Vmag]
	weighted_I = [Wmean_Imag, Werr_Imag, mean_Imag, err_Imag]
	
	input ('\n============ STOP HERE! ==========\n')

	pars['memb_bellini'] = memb_bellini
	pars['weighted_V'], pars['weighted_I'] = weighted_V, weighted_I
	return pars
	
def save_read(memb, writ, rec, pars, alpha=2.0):	# 01.abr.21: OK!
	'''
	Saves and/or reads the quantities needed in the membership function.
	'''
	print()
	path = os.path.dirname(os.path.realpath(__file__))

	if 'list_mags' in pars.keys():
		list_V, list_I, list_Ks = pars['list_mags']
		Vmag_mean, Vmag_err, Vmag_med, Vmag_std = list_V
		Imag_mean, Imag_err, Imag_med, Imag_std = list_I
		Ksmag_mean, Ksmag_err, Ksmag_med, Ksmag_std = list_Ks
	
	if 'c_rad' in pars.keys(): c_rad = pars['c_rad']

	if 'pms_1d_2d' in pars.keys():
		pms_RA_1d, pms_DE_1d, pms_2d, f_pms_2d = pars['pms_1d_2d']
		RApm, pmRAsig, pmRAerr, RApm_f, pmRAsig_f, pmRAerr_f = pms_RA_1d
		DEpm, pmDEsig, pmDEerr, DEpm_f, pmDEsig_f, pmDEerr_f = pms_DE_1d
		RApm_2d, DEpm_2d, pmRAsig_2d, pmDEsig_2d = pms_2d
		RApm_field_2d, DEpm_field_2d, pmRAsig_field_2d, pmDEsig_field_2d = f_pms_2d
	
	if 'memb_bellini' in pars.keys():
		memb_bellini = pars['memb_bellini']
		weighted_V, weighted_I = pars['weighted_V'], pars['weighted_I']

	filepath = f'{path}/{NGC}/Aux_Cats/{NGC}_prevMemb.txt'
	if not os.path.isfile(filepath):
		fil = open(filepath, "w").close()
	with open(filepath, 'r+') as fil:
		lines_save = fil.readlines()
	for i in range(len(lines_save)):
		if lines_save[i] == 'p_member = \n': starting = i+1
		elif lines_save[i][-2:] == ']\n': 
			finishing = i
			break
	
	### PAREI AQUIII::: TERMINAR URGENTE!!! E RODAR PARA OS 6 AGLOMERADOS!
	#.Reading previous values (if memb is True)
	if memb or rec == ['King profile (half-light radius)']:
		c_rad = float(lines_save[0].split('=')[1])
		alpha = float(lines_save[1].split('=')[1])
		p_member = []
		for i in range(starting, finishing+1):
			lista = lines_save[i].split()
			if lista[0][0] == '[': lista[0] = lista[0][1:]
			if lista[-1][-1] == ']': lista[-1] = lista[-1][:-1]
			for i in range(len(lista)): 
				lista[i] = float(lista[i])
				p_member.append(lista[i])
		p_member = np.array(p_member)

		for i in range(len(lines_save)):
			if lines_save[i][:8] == '1D: RApm': break
			
		#.1D: RApm & DEpm (cluster and field)
		RApm, pmRAsig, pmRAerr = lines_save[i+1].split(',')
		DEpm, pmDEsig, pmDEerr = lines_save[i+3].split(',')
		RApm, pmRAsig = [float(RApm[1:-1])], [[float(pmRAsig[3:-2])]]
		DEpm, pmDEsig = [float(DEpm[1:-1])], [[float(pmDEsig[3:-2])]]
		pmRAerr, pmDEerr = [[float(pmRAerr[3:-3])]], [[float(pmDEerr[3:-3])]]

		RApm_f, pmRAsig_f, pmRAerr_f = lines_save[i+6].split(',')
		DEpm_f, pmDEsig_f, pmDEerr_f = lines_save[i+8].split(',')
		RApm_f, pmRAsig_f = [float(RApm_f[1:-1])], [[float(pmRAsig_f[3:-2])]]
		DEpm_f, pmDEsig = [float(DEpm_f[1:-1])], [[float(pmDEsig_f[3:-2])]]
		pmRAerr_f, pmDEerr_f = [[float(pmRAerr_f[3:-3])]], [[float(pmDEerr_f[3:-3])]]

		pms_RA_1d = RApm, pmRAsig, pmRAerr, RApm_f, pmRAsig_f, pmRAerr_f
		pms_DE_1d = DEpm, pmDEsig, pmDEerr, DEpm_f, pmDEsig_f, pmDEerr_f
		
		#.2D: RApm & DEpm (cluster and field)
		RApm_2d, pmRAsig_2d = lines_save[i+10].split('=')[1].split(',')
		RApm_2d, pmRAsig_2d = float(RApm_2d), float(pmRAsig_2d)
		DEpm_2d, pmDEsig_2d = lines_save[i+11].split('=')[1].split(',')
		DEpm_2d, pmDEsig_2d = float(DEpm_2d), float(pmDEsig_2d)
		RApm_field_2d, pmRAsig_field_2d = lines_save[i+13].split('=')[1].split(',')
		RApm_field_2d, pmRAsig_field_2d = float(RApm_field_2d), float(pmRAsig_field_2d)
		DEpm_field_2d, pmDEsig_field_2d = lines_save[i+14].split('=')[1].split(',')
		DEpm_field_2d, pmDEsig_field_2d = float(DEpm_field_2d), float(pmDEsig_field_2d)

		pms_2d = RApm_2d, DEpm_2d, pmRAsig_2d, pmDEsig_2d
		f_pms_2d = RApm_field_2d, DEpm_field_2d, pmRAsig_field_2d, pmDEsig_field_2d

		Vmag_mean = float(lines_save[i+16].split()[2])
		Vmag_err = float(lines_save[i+16].split()[4])
		Vmag_med = float(lines_save[i+17].split()[2])
		Vmag_std = float(lines_save[i+17].split()[6])
		Imag_mean = float(lines_save[i+19].split()[2])
		Imag_err = float(lines_save[i+19].split()[4])
		Imag_med = float(lines_save[i+20].split()[2])
		Imag_std = float(lines_save[i+20].split()[6])
		Ksmag_mean = float(lines_save[i+22].split()[2])
		Ksmag_err = float(lines_save[i+22].split()[4])
		Ksmag_med = float(lines_save[i+23].split()[2])
		Ksmag_std = float(lines_save[i+23].split()[6])
		
		list_V = [Vmag_mean, Vmag_err, Vmag_med, Vmag_std]
		list_I = [Imag_mean, Imag_err, Imag_med, Imag_std]
		list_Ks = [Ksmag_mean, Ksmag_err, Ksmag_med, Ksmag_std]

		pms_RA_1d = (RApm, pmRAsig, pmRAerr, RApm_f, pmRAsig_f, pmRAerr_f)
		pms_DE_1d = (DEpm, pmDEsig, pmDEerr, DEpm_f, pmDEsig_f, pmDEerr_f)
		pms_2d = (RApm_2d, DEpm_2d, pmRAsig_2d, pmDEsig_2d)
		f_pms_2d = (RApm_field_2d, DEpm_field_2d, pmRAsig_field_2d, pmDEsig_field_2d)

		pars['c_rad'] = c_rad
		pars['list_mags'] = [list_V, list_I, list_Ks]
		pars['pms_1d_2d'] = [pms_RA_1d, pms_DE_1d, pms_2d, f_pms_2d]
		pars['memb_bellini'] = p_member
		pars['weighted_V'], pars['weighted_I'] = list_V, list_I	# is it right?
		return pars
	
	#.Load only the not marked rec_values
	for i in range(len(lines_save)):
		if lines_save[i][:8] == '1D: RApm': break
	
	if 'Sigma-clipping (<V>, <I> & <Ks>)' not in rec:
		Vmag_mean = float(lines_save[i+16].split()[2])
		Vmag_err = float(lines_save[i+16].split()[4])
		Vmag_med = float(lines_save[i+17].split()[2])
		Vmag_std = float(lines_save[i+17].split()[6])
		Imag_mean = float(lines_save[i+19].split()[2])
		Imag_err = float(lines_save[i+19].split()[4])
		Imag_med = float(lines_save[i+20].split()[2])
		Imag_std = float(lines_save[i+20].split()[6])
		Ksmag_mean = float(lines_save[i+22].split()[2])
		Ksmag_err = float(lines_save[i+22].split()[4])
		Ksmag_med = float(lines_save[i+23].split()[2])
		Ksmag_std = float(lines_save[i+23].split()[6])
		pars['list_mags'] = [list_V, list_I, list_Ks]
		
	elif 'Proper motions (from Gaia EDR3)' not in rec:
		c_rad = float(lines_save[0].split('=')[1].split('\n')[0])
		p_member = []
		for j in range(starting,finishing+1):
			lista = lines_save[j].split()
			if lista[0][0] == '[': lista[0] = lista[0][1:]
			if lista[-1][-1] == ']': lista[-1] = lista[-1][:-1]
			for k in range(len(lista)): 
				lista[k] = float(lista[k])
				p_member.append(lista[k])
		p_member = np.array(p_member)
		
		#.1D: RApm & DEpm (cluster and field)
		RApm, pmRAsig, pmRAerr = lines_save[i+1].split(',')
		DEpm, pmDEsig, pmDEerr = lines_save[i+3].split(',')
		RApm, pmRAsig = [float(RApm[1:-1])], [[float(pmRAsig[3:-2])]]
		DEpm, pmDEsig = [float(DEpm[1:-1])], [[float(pmDEsig[3:-2])]]
		pmRAerr, pmDEerr = [[float(pmRAerr[3:-3])]], [[float(pmDEerr[3:-3])]]

		RApm_f, pmRAsig_f, pmRAerr_f = lines_save[i+6].split(',')
		DEpm_f, pmDEsig_f, pmDEerr_f = lines_save[i+8].split(',')
		RApm_f, pmRAsig_f = [float(RApm_f[1:-1])], [[float(pmRAsig_f[3:-2])]]
		DEpm_f, pmDEsig = [float(DEpm_f[1:-1])], [[float(pmDEsig_f[3:-2])]]
		pmRAerr_f, pmDEerr_f = [[float(pmRAerr_f[3:-3])]], [[float(pmDEerr_f[3:-3])]]

		pms_RA_1d = RApm, pmRAsig, pmRAerr, RApm_f, pmRAsig_f, pmRAerr_f
		pms_DE_1d = DEpm, pmDEsig, pmDEerr, DEpm_f, pmDEsig_f, pmDEerr_f
		
		#.2D: RApm & DEpm (cluster and field)
		RApm_2d, pmRAsig_2d = lines_save[i+10].split('=')[1].split(',')
		RApm_2d, pmRAsig_2d = float(RApm_2d), float(pmRAsig_2d)
		DEpm_2d, pmDEsig_2d = lines_save[i+11].split('=')[1].split(',')
		DEpm_2d, pmDEsig_2d = float(DEpm_2d), float(pmDEsig_2d)
		RApm_field_2d, pmRAsig_field_2d = lines_save[i+13].split('=')[1].split(',')
		RApm_field_2d, pmRAsig_field_2d = float(RApm_field_2d), float(pmRAsig_field_2d)
		DEpm_field_2d, pmDEsig_field_2d = lines_save[i+14].split('=')[1].split(',')
		DEpm_field_2d, pmDEsig_field_2d = float(DEpm_field_2d), float(pmDEsig_field_2d)

		pms_2d = RApm_2d, DEpm_2d, pmRAsig_2d, pmDEsig_2d
		f_pms_2d = RApm_field_2d, DEpm_field_2d, pmRAsig_field_2d, pmDEsig_field_2d

		pars['c_rad'] = c_rad
		pars['memb_bellini'] = p_member
		pars['pms_1d_2d'] = [pms_RA_1d, pms_DE_1d, pms_2d, f_pms_2d]
	
	if not writ:
		return pars
	
	#.Write output in prevMemb.txt file
	with open(filepath, 'r+') as fil:
		fil.truncate()
		fil.write(f'center_radius = {c_rad:.3f}\nalpha = {alpha:.2f}\n\n')
		fil.write('p_member = \n[')
		for i in range(len(memb_bellini)):
			if (i+1)%4 == 0 and i != len(memb_bellini)-1 :
				fil.write('%.8E\n'%memb_bellini[i])
			elif i == len(memb_bellini)-1:
				fil.write('%.8E]\n\n'%memb_bellini[i])
			else: fil.write('%.8E\t'%memb_bellini[i])
		s = f'[{pms_RA_1d[0][0]}], [{pms_RA_1d[1][0]}], [{pms_RA_1d[2][0]}]'
		fil.write(f'1D: RApm, pmRAsig, pmRAerr = \n{s}\n')
		s = f'[{pms_DE_1d[0][0]}], [{pms_DE_1d[1][0]}], [{pms_DE_1d[2][0]}]'
		fil.write(f'1D: DEpm, pmDEsig, pmDEerr = \n{s}\n\n')
		s = f'[{pms_RA_1d[3][0]}], [{pms_RA_1d[4][0]}], [{pms_RA_1d[5][0]}]'
		fil.write(f'1D-field: RApm_f, pmRAsig_f, pmRAerr_f = \n{s}\n')
		aux = f'[{pms_DE_1d[3][0]}], [{pms_DE_1d[4][0]}], [{pms_DE_1d[5][0]}]'
		fil.write(f'1D-field: DEpm_f, pmDEsig_f, pmDEerr_f = \n{aux}\n\n')

		fil.write(f'2D: RApm, pmRAsig = {pms_2d[0]}, {pms_2d[2]}\n')
		fil.write(f'2D: DEpm, pmDEsig = {pms_2d[1]}, {pms_2d[3]}\n\n')
		fil.write(f'2D-field: RApm, pmRAsig = {f_pms_2d[0]}, {f_pms_2d[2]}\n')
		fil.write(f'2D-field: DEpm, pmDEsig = {f_pms_2d[1]}, {f_pms_2d[3]}\n\n')

		err = Vmag_err if NGC != 'NGC6717' else 0.10000
		fil.write(f'<V> = {weighted_V[0]:.5f} +/- {err:.5f}\n')
		fil.write(f'V_med = {list_V[2]:.5f} ; V_std = {list_V[3]:.5f}\n\n')
		fil.write(f'<I> = {weighted_I[0]:.5f} +/- {list_I[1]:.5f}\n')
		fil.write(f'I_med = {list_I[2]:.5f} ; I_std = {list_I[3]:.5f}\n\n')
		fil.write(f'Ks = {list_Ks[0]:.5f} +/- {list_Ks[1]:.5f}\n')
		fil.write(f'Ks_med = {list_Ks[2]:.5f} ; Ks_std = {list_Ks[3]:.5f}')
	
	return pars

if __name__ == '__main__':
	main()
	
########################################################################
############################## EOF #####################################
########################################################################