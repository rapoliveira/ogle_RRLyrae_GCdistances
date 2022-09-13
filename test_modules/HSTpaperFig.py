import numpy as np, matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from termcolor import colored, cprint
from astropy.coordinates import SkyCoord
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker

global NGC
NGC = ['NGC6362','NGC6717','NGC6723']
RA = {'NGC6362':[],'NGC6717':[],'NGC6723':[]}
DEC = {'NGC6362':[],'NGC6717':[],'NGC6723':[]}
Period = {'NGC6362':[],'NGC6717':[],'NGC6723':[]}
magV = {'NGC6362':[],'NGC6717':[],'NGC6723':[]}
memb = {'NGC6362':[],'NGC6717':[],'NGC6723':[]}
Wmean_Vmag = {'NGC6362':[],'NGC6717':[],'NGC6723':[]}
std_Vmag = {'NGC6362':[],'NGC6717':[],'NGC6723':[]}
mag_lines, mag_titles, memb_lines, memb_titles = [],[],[],[]

for i in range(len(NGC)):
	
	print (i)
	with open('%s/Exported_Planes/txt_format/%s_Clement_vs_Gaia.txt'% \
			(NGC[i],NGC[i]),'r') as mags:
		mag_lines.append(mags.readlines())
		mag_titles.append(mag_lines[i][0].split('\t'))
	for j in range(2,len(mag_lines[i])):
		RA[NGC[i]].append(mag_lines[i][j].split('\t')[mag_titles[i].index('RA')])		
		DEC[NGC[i]].append(mag_lines[i][j].split('\t')[mag_titles[i].index('DEC')])
		Period[NGC[i]].append(float(mag_lines[i][j].split('\t')[mag_titles[i].index('Period')]))
		magV[NGC[i]].append(float(mag_lines[i][j].split('\t')[mag_titles[i].index('<V>')]))
	RA[NGC[i]], DEC[NGC[i]], Period[NGC[i]], magV[NGC[i]] = np.array(RA[NGC[i]]),\
			np.array(DEC[NGC[i]]),np.array(Period[NGC[i]]),np.array(magV[NGC[i]])
				
	with open('%s/Aux_Cats/%s_prevMemb.txt'%(NGC[i],NGC[i]),'r') as membs:
		memb_lines.append(membs.readlines())
	for k in range(len(memb_lines[i])):
		if memb_lines[i][k][0] == '[': start = k
		if memb_lines[i][k][-2:] == ']\n':
			end = k
			break
	for k in range(start,end+1):
		if k == start:
			aux = memb_lines[i][k].split('\t')
			for m in range(len(aux)):
				if aux[m][0] == '[' and NGC[i] != 'NGC6717': memb[NGC[i]].append(float(aux[m][1:]))
				elif aux[m][0] == '[' and NGC[i] == 'NGC6717': memb[NGC[i]].append(float(aux[m][1:-2]))
				else: memb[NGC[i]].append(float(aux[m]))
		elif k != start and k != end:
			for m in range(4):
				memb[NGC[i]].append(float(memb_lines[i][k].split('\t')[m]))
		elif k == end:
			aux = memb_lines[i][k].split('\t')
			for m in range(len(aux)):
				if aux[m][-2:] != ']\n': memb[NGC[i]].append(float(aux[m]))
				else: memb[NGC[i]].append(float(aux[m][:-2]))
	for k in range(len(memb_lines[i])):
		if memb_lines[i][k][0:3] == '<V>':
			Wmean_Vmag[NGC[i]].append(float(memb_lines[i][k].split()[2]))
			std_Vmag[NGC[i]].append(float(memb_lines[i][k].split()[4]))
	
	memb[NGC[i]], Wmean_Vmag[NGC[i]] = np.array(memb[NGC[i]]), np.array(Wmean_Vmag[NGC[i]])
	std_Vmag[NGC[i]] = np.array(std_Vmag[NGC[i]])

fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(5.5,5.5),sharex=True)
ax3.set_xlabel('Period (days)',fontsize=13)
ax1.set_ylabel(r'$\langle V\rangle$',fontsize=13)
ax2.set_ylabel(r'$\langle V\rangle$',fontsize=13)
ax3.set_ylabel(r'$\langle V\rangle$',fontsize=13)

ax1.tick_params(direction='in',which='both')
ax2.tick_params(direction='in',which='both')
ax3.tick_params(direction='in',which='both')
ax1.xaxis.set_ticks_position('both')
ax1.yaxis.set_ticks_position('both')
ax2.xaxis.set_ticks_position('both')
ax2.yaxis.set_ticks_position('both')
ax3.xaxis.set_ticks_position('both')
ax3.yaxis.set_ticks_position('both')
ax1.xaxis.set_major_locator(ticker.AutoLocator())
ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax1.yaxis.set_major_locator(ticker.AutoLocator())
ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax2.xaxis.set_major_locator(ticker.AutoLocator())
ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax2.yaxis.set_major_locator(ticker.AutoLocator())
ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax3.xaxis.set_major_locator(ticker.AutoLocator())
ax3.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax3.yaxis.set_major_locator(ticker.AutoLocator())
ax3.yaxis.set_minor_locator(ticker.AutoMinorLocator())

#ax1.set_title(r'%s$\,$%s: $\langle %s\rangle$ vs. Period'%(NGC[:3],NGC[3:],filt))
#ax1.scatter(Period['NGC6362'],magV['NGC6362'],s=10)
ax1.set_xlim(min(Period['NGC6362'])-0.05,max(Period['NGC6362'])+0.05)
ax1.set_ylim(round(Wmean_Vmag['NGC6362'][0],2)+0.5,round(Wmean_Vmag['NGC6362'][0],2)-0.5)
#ax2.set_xlim(min(Period['NGC6717'])-0.05,max(Period['NGC6717'])+0.05)
ax2.set_ylim(Wmean_Vmag['NGC6717'][0]+0.5,Wmean_Vmag['NGC6717'][0]-0.5)
#ax3.set_xlim(min(Period['NGC6723'])-0.05,max(Period['NGC6723'])+0.05)
ax3.set_ylim(Wmean_Vmag['NGC6723'][0]+0.5,Wmean_Vmag['NGC6723'][0]-0.5)
pcm = ax1.scatter(Period['NGC6362'],magV['NGC6362'],s=17,c=memb['NGC6362'],cmap='Greys',norm=colors.Normalize(vmin=0, vmax=1.0),zorder=3)
print (memb['NGC6717'])
pcm = ax2.scatter(Period['NGC6717'],magV['NGC6717'],s=17,c=memb['NGC6717'],cmap='Greys',norm=colors.Normalize(vmin=0, vmax=1.0),zorder=3)
pcm = ax3.scatter(Period['NGC6723'],magV['NGC6723'],s=17,c=memb['NGC6723'],cmap='Greys',norm=colors.Normalize(vmin=0, vmax=1.0),zorder=3)

#cbar = plt.colorbar(pcm,cax=cax, orientation='vertical',ticks=[0.0,0.25,0.50,0.75,1.0])
#cax.xaxis.set_ticks_position('top')
cbar_ax = fig.add_axes([0.84, 0.075, 0.02, 0.87])
#divider = make_axes_locatable(cbar_ax)
#cax = divider.append_axes('right', size='1%', pad=0.25)
#cbar_ax.set_label('Membership', fontsize=10)
cbar = fig.colorbar(pcm, cax=cbar_ax,orientation='vertical',ticks=[0.0,0.25,0.50,0.75,1.00])
cbar.set_label('Membership', fontsize=13)
cbar.ax.tick_params(labelsize=11)

#cb_ax = fig.add_axes([0.83, 0.1, 0.02, 0.8])
#cbar = fig.colorbar(im, cax=cb_ax)

ax1.tick_params(axis='both', which='major', labelsize=12)
ax2.tick_params(axis='both', which='major', labelsize=12)
ax3.tick_params(axis='both', which='major', labelsize=12)
deviation1 = std_Vmag['NGC6362'][0]#*np.sqrt(len(Period['NGC6362']))
deviation2 = std_Vmag['NGC6717'][0]#*np.sqrt(len(Period['NGC6717']))
deviation3 = std_Vmag['NGC6723'][0]#*np.sqrt(len(Period['NGC6723']))
print (deviation1,deviation2,deviation3)
ax1.plot([min(Period['NGC6362'])-0.05,max(Period['NGC6362'])+0.05],[Wmean_Vmag['NGC6362'],Wmean_Vmag['NGC6362']],color='k',ls='--',lw=1.5,label=r'$%.3f\pm%.3f$'%(Wmean_Vmag['NGC6362'],deviation1))
ax1.plot([min(Period['NGC6362'])-0.05,max(Period['NGC6362'])+0.05],[Wmean_Vmag['NGC6362']+deviation1,Wmean_Vmag['NGC6362']+deviation1],color='gray',ls=':',lw=1)
ax1.plot([min(Period['NGC6362'])-0.05,max(Period['NGC6362'])+0.05],[Wmean_Vmag['NGC6362']-deviation1,Wmean_Vmag['NGC6362']-deviation1],color='gray',ls=':',lw=1)
ax2.plot([min(Period['NGC6362'])-0.05,max(Period['NGC6362'])+0.05],[Wmean_Vmag['NGC6717'],Wmean_Vmag['NGC6717']],color='k',ls='--',lw=1.5,label=r'$%.3f\pm%.3f$'%(Wmean_Vmag['NGC6717'],deviation2))
ax2.plot([min(Period['NGC6362'])-0.05,max(Period['NGC6362'])+0.05],[Wmean_Vmag['NGC6717']+deviation2,Wmean_Vmag['NGC6717']+deviation2],color='gray',ls=':',lw=1)
ax2.plot([min(Period['NGC6362'])-0.05,max(Period['NGC6362'])+0.05],[Wmean_Vmag['NGC6717']-deviation2,Wmean_Vmag['NGC6717']-deviation2],color='gray',ls=':',lw=1)
ax3.plot([min(Period['NGC6362'])-0.05,max(Period['NGC6362'])+0.05],[Wmean_Vmag['NGC6723'],Wmean_Vmag['NGC6723']],color='k',ls='--',lw=1.5,label=r'$%.3f\pm%.3f$'%(Wmean_Vmag['NGC6723'],deviation3))
ax3.plot([min(Period['NGC6362'])-0.05,max(Period['NGC6362'])+0.05],[Wmean_Vmag['NGC6723']+deviation3,Wmean_Vmag['NGC6723']+deviation3],color='gray',ls=':',lw=1)
ax3.plot([min(Period['NGC6362'])-0.05,max(Period['NGC6362'])+0.05],[Wmean_Vmag['NGC6723']-deviation3,Wmean_Vmag['NGC6723']-deviation3],color='gray',ls=':',lw=1)

plt.subplots_adjust(top=0.965,bottom=0.095,left=0.155,right=0.810,hspace=0.110)
ax1.legend(loc=2,fontsize=12)
ax2.legend(loc=2,fontsize=12)
ax3.legend(loc=2,fontsize=12)
ax1.text(0.595,15.71,r'NGC$\,$6362',fontsize=13)#,weight='heavy')
ax2.text(0.595,16.12,r'NGC$\,$6717',fontsize=13)
ax3.text(0.595,15.88,r'NGC$\,$6723',fontsize=13)
#plt.tight_layout()
plt.show()


'''
def mag_vsPeriod(memb,mag,per,filt,alpha,list_med):
	
	#mag_mean,mag_err,mag_med,mag_std = sig_clip(mag,filt,alpha,20,False)
	#list_med.extend([mag_mean,mag_err,mag_med,mag_std])
			
	#mag_inside_lims, per_inside_lims, size_inside, size = [],[],[],[]
	#mag_outside_lims, per_outside_lims, size_outside = [],[],[]
	#m2 = (Vmag > Vmag_med - alpha*Vmag_std) & (Vmag < Vmag_med + alpha*Vmag_std)
	for i in range(len(mag)):
		if mag_med - alpha*mag_std < mag[i] < mag_med + alpha*mag_std:
			mag_inside_lims[NGC].append(round(mag[i],3))
			per_inside_lims[NGC].append(per[i])
			#size_inside.append(size[i])
		else:
			mag_outside_lims[NGC].append(round(mag[i],3))
			per_outside_lims[NGC].append(per[i])
			#size_outside.append(size[i])

	fig = plt.figure(figsize=(7.8,2.5)) # x_old = 6
	gs = matplotlib.gridspec.GridSpec(1,2,width_ratios=[3.5,0.9])
	ax1, ax2 = plt.subplot(gs[0]), plt.subplot(gs[1])
	ax1.set_xlabel('Period (days)',fontsize=11)
	if filt == 'Ks': filt = r'K_{\rm{S}}'
	ax1.set_ylabel(r'$\langle %s\rangle$'%filt,fontsize=11)
	ax1.tick_params(direction='in',which='both')
	ax1.xaxis.set_ticks_position('both')
	ax1.yaxis.set_ticks_position('both')
	ax1.xaxis.set_major_locator(ticker.AutoLocator())
	ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	ax1.yaxis.set_major_locator(ticker.AutoLocator())
	ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())

	ax1.set_title(r'%s$\,$%s: $\langle %s\rangle$ vs. Period'%(NGC[:3],NGC[3:],filt))
	ax1.set_xlim(min(per)-0.05,max(per)+0.05)
	ax1.set_ylim(mag_mean-2.1,mag_mean+2.1)
	
	for i in range(len(mag_inside_lims[NGC])):
		if per_inside_lims[NGC][i] >= 0.4000: ax1.scatter(per_inside_lims[NGC][i],mag_inside_lims[NGC][i],zorder=5,s=40,color='red',edgecolors = 'black')		# edgecolors = 'black' #################
		else: ax1.scatter(per_inside_lims[NGC][i],mag_inside_lims[NGC][i],zorder=5,s=40,color='blue',edgecolors = 'black')
	for i in range(len(mag_outside_lims[NGC])):
		ax1.scatter(per_outside_lims[NGC][i],mag_outside_lims[NGC][i],zorder=4,color='silver')
	ax1.plot([min(per)-0.05,max(per)+0.05],[mag_mean,mag_mean],color='k',ls='--',lw=1,zorder=1,label=r'%.3f$\,\pm\,$%.3f'%(mag_mean,mag_err))
	ax1.plot([min(per)-0.05,max(per)+0.05],[mag_mean+mag_err,mag_mean+mag_err],color='gray',ls=':',lw=0.8,zorder=2)
	ax1.plot([min(per)-0.05,max(per)+0.05],[mag_mean-mag_err,mag_mean-mag_err],color='gray',ls=':',lw=0.8,zorder=3)
	ax1.legend(loc=1,fontsize=11)
	ax1.invert_yaxis()
	ax2.axis('off')
	ax2.set_ylim(mag_mean-2.1,mag_mean+2.1)
	ax2.invert_yaxis()
	ax2.hist(mag,bins=25,range=[mag_mean-3.0,mag_mean+3.0],histtype='step',orientation='horizontal',density=1)
	plt.tight_layout()
	plt.subplots_adjust(wspace=0.05)
	plt.show()
	
	return list_med
	
#print (magV['NGC6304'],len(magV['NGC6304']))
#print (Period['NGC6304'],len(Period['NGC6304']))
#input()

#memb, alpha, list_med_V = False, 2.5, {'NGC6304':[],'NGC6362':[],'NGC6652':[],'NGC6723':[]}
#global mag_inside_lims, per_inside_lims, mag_outside_lims, per_outside_lims
#mag_inside_lims = {'NGC6304':[],'NGC6362':[],'NGC6652':[],'NGC6723':[]}
#per_inside_lims = {'NGC6304':[],'NGC6362':[],'NGC6652':[],'NGC6723':[]}
#mag_outside_lims = {'NGC6304':[],'NGC6362':[],'NGC6652':[],'NGC6723':[]}
#per_outside_lims = {'NGC6304':[],'NGC6362':[],'NGC6652':[],'NGC6723':[]}

#for k in range(len(NGCs)):
#	NGC = NGCs[k]
#	list_med_V[NGC].append(mag_vsPeriod(memb,magV[NGC],Period[NGC],'V',alpha,list_med_V[NGC]))
#	Vmag_mean,Vmag_err,Vmag_med,Vmag_std = list_med_V[NGC][0],list_med_V[NGC][1],list_med_V[NGC][2],list_med_V[NGC][3]
#	print (Vmag_mean, Vmag_err, Vmag_med, Vmag_std)
#print(mag_inside_lims['NGC6304'],len(mag_inside_lims['NGC6304']))
#print(per_inside_lims['NGC6304'],len(per_inside_lims['NGC6304']))
#input()

#print (list_med_V['NGC6304'])
#print (list_med_V['NGC6362'])
#print (list_med_V['NGC6652'])
#print (list_med_V['NGC6723'])
# list_med_V =>>>> Vmag_mean, Vmag_err, Vmag_med, Vmag_std

filt = 'V'
fig = plt.figure(figsize=(7.8,6))
ax1 = plt.subplot2grid((4, 4), (0, 0), colspan=3)
ax2 = plt.subplot2grid((4, 4), (0, 3))
ax3 = plt.subplot2grid((4, 4), (1, 0), colspan=3)
ax4 = plt.subplot2grid((4, 4), (1, 3))
ax5 = plt.subplot2grid((4, 4), (2, 0), colspan=3)
ax6 = plt.subplot2grid((4, 4), (2, 3))
ax7 = plt.subplot2grid((4, 4), (3, 0), colspan=3)
ax8 = plt.subplot2grid((4, 4), (3, 3))
axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]

for l in np.arange(0,2*len(NGCs),2):	
	
	nnn = int(l*0.5)
	NGC = NGCs[nnn]
	#abc = "41"+str(+1)
	#ax = fig.add_subplot(int(abc))
	
	per = Period[NGC]
	mag_mean, mag_err, mag_med, mag_std = list_med_V[NGC][0],list_med_V[NGC][1],list_med_V[NGC][2],list_med_V[NGC][3]
	
	gs = matplotlib.gridspec.GridSpec(4, 2,
                       width_ratios=[3, 1],
                       height_ratios=[4, 1]
                       )
	
	ax1, ax2 = plt.subplot(gs[0]), plt.subplot(gs[1])
	ax3, ax4 = plt.subplot(gs[2]), plt.subplot(gs[3])
	ax5, ax6 = plt.subplot(gs[4]), plt.subplot(gs[5])
	ax7, ax8 = plt.subplot(gs[6]), plt.subplot(gs[7])

	#ax1.set_xlabel('Period (days)',fontsize=11)
	#if filt == 'Ks': filt = r'K_{\rm{S}}'
	if NGC != 'NGC6652':
		axes[l].set_ylabel(r'$\langle %s\rangle$'%filt,fontsize=11)
	else: axes[l].set_ylabel(r'$\langle B\rangle$',fontsize=11)
	axes[l].tick_params(direction='in',which='both',labelsize=11)
	axes[l].xaxis.set_ticks_position('both')
	axes[l].yaxis.set_ticks_position('both')
	axes[l].xaxis.set_major_locator(ticker.AutoLocator())
	axes[l].xaxis.set_minor_locator(ticker.AutoMinorLocator())
	axes[l].yaxis.set_major_locator(ticker.AutoLocator())
	axes[l].yaxis.set_minor_locator(ticker.AutoMinorLocator())

	#ax1.set_title(r'%s$\,$%s: $\langle %s\rangle$ vs. Period'%(NGC[:3],NGC[3:],filt))
	#axes[l].set_xlim(min(per)-0.05,max(per)+0.05)
	axes[l].set_xlim(0.2,0.9)
	if NGC in ['NGC6304','NGC6652']:
		axes[l].set_ylim(mag_mean-2.1,mag_mean+2.1)
	else: axes[l].set_ylim(mag_mean-0.75,mag_mean+0.75)

	print (len(magV[NGC]),len(per_inside_lims[NGC]), len(mag_inside_lims[NGC]))
	for i in range(len(per_inside_lims[NGC])):
		if per_inside_lims[NGC][i] >= 0.4000: axes[l].scatter(per_inside_lims[NGC][i],mag_inside_lims[NGC][i],zorder=5,s=40,color='red',edgecolors = 'black')		# edgecolors = 'black' #################
		else: axes[l].scatter(per_inside_lims[NGC][i],mag_inside_lims[NGC][i],zorder=5,s=40,color='blue',edgecolors = 'black')
	for i in range(len(mag_outside_lims[NGC])):
		axes[l].scatter(per_outside_lims[NGC][i],mag_outside_lims[NGC][i],zorder=4,color='silver')
	axes[l].plot([0.2,0.9],[mag_mean,mag_mean],color='k',ls='--',lw=1,zorder=1,label=r'%.3f$\,\pm\,$%.3f'%(mag_mean,mag_err))
	axes[l].plot([0.2,0.9],[mag_mean+mag_err,mag_mean+mag_err],color='gray',ls=':',lw=0.8,zorder=2)
	axes[l].plot([0.2,0.9],[mag_mean-mag_err,mag_mean-mag_err],color='gray',ls=':',lw=0.8,zorder=3)
	if NGC == 'NGC6304': axes[l].legend(loc=2,fontsize=11)
	else: axes[l].legend(loc=1,fontsize=11)
	axes[l].invert_yaxis()
	axes[l+1].axis('off')
	if NGC in ['NGC6304','NGC6652']:
		axes[l+1].set_ylim(mag_mean-2.1,mag_mean+2.1)
		axes[l+1].hist(magV[NGC],bins=10,range=[mag_mean-2.5,mag_mean+2.5],histtype='step',orientation='horizontal',density=1)
	else:
		axes[l+1].set_ylim(mag_mean-0.75,mag_mean+0.75)
		axes[l+1].hist(magV[NGC],bins=15,range=[mag_mean-0.8,mag_mean+0.8],histtype='step',orientation='horizontal',density=1)
	axes[l+1].invert_yaxis()
	if l == 2*len(NGCs)-2:
		axes[l].set_xlabel('Period (days)',fontsize=11)
	else:
		axes[l].set_xticklabels([])

plt.tight_layout()
plt.subplots_adjust(wspace=0.07,hspace=0.1)
plt.show()


#list_med_Ks = mag_vsPeriod(memb,Ksmag,'Ks',alpha,list_med_Ks)
#Ksmag_mean,Ksmag_err,Ksmag_med,Ksmag_std = list_med_Ks[0],list_med_Ks[1],list_med_Ks[2],list_med_Ks[3]

'''
