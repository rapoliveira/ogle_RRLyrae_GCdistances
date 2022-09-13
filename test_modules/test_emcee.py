import numpy as np, inquirer as iq, matplotlib, os
import matplotlib.pyplot as plt
import scipy.optimize as optimization

with open('Cassisi_zahb_T385aes.txt','r') as Cas:
	lines = Cas.readlines()[4:]
	titles = lines[0].split()[1:]

#[M/H] = [Fe/H] + np.log10(0.638 * 10[α/Fe] + 0.362)	=> My thesis
#alpha = np.log10((10**(MH - FEH) - 0.362)/0.638)		  (inverted)

MH = np.zeros(len(lines)-2)
FEH, Mv, Mi = np.copy(MH), np.copy(MH), np.copy(MH)
for i in range(len(MH)):
	MH[i]  = float(lines[i+2].split()[0])
	FEH[i] = float(lines[i+2].split()[0]) - np.log10(0.694 * 10**0.445 + 0.306)
	Mv[i]  = float(lines[i+2].split()[4]) - 0.08
	Mi[i]  = float(lines[i+2].split()[4]) - float(lines[i+2].split()[7]) - 0.08

###.Fitting a model to data (emcee docs)
#https://emcee.readthedocs.io/en/stable/tutorials/line/

#.Least squares estimation
f_true = 0.1
x, y, yerr = np.copy(FEH), np.copy(Mv), np.repeat(0.05,len(MH))
y += np.abs(f_true * y) * np.random.randn(len(MH))
A = np.vander(x, 2)
C = np.diag(yerr * yerr)
ATA = np.dot(A.T, A / (yerr ** 2)[:, None])
cov = np.linalg.inv(ATA)
w = np.linalg.solve(ATA, np.dot(A.T, y / yerr ** 2))
print("\nLeast-squares estimates:")
print("m = {0:.3f} ± {1:.3f}".format(w[0], np.sqrt(cov[0, 0])))
print("b = {0:.3f} ± {1:.3f}".format(w[1], np.sqrt(cov[1, 1])))


#.Maximum likelihood estimation
def log_likelihood(theta, x, y, yerr):
    m, b, log_f = theta
    model = m * x + b
    sigma2 = yerr ** 2 + model ** 2 * np.exp(2 * log_f)
    return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))

from scipy.optimize import minimize
m_true, b_true = 0.3, 1
np.random.seed(42)
nll = lambda *args: -log_likelihood(*args)
initial = np.array([m_true, b_true, np.log(f_true)]) + 0.1 * np.random.randn(3)
soln = minimize(nll, initial, args=(x, y, yerr))
m_ml, b_ml, log_f_ml = soln.x
print("\nMaximum likelihood estimates:")
print("m = {0:.3f}".format(m_ml))
print("b = {0:.3f}".format(b_ml))
print("f = {0:.3f}".format(np.exp(log_f_ml)))


#.Marginalization & uncertainty estimator
def log_prior(theta):
    m, b, log_f = theta
    if 0 < m < 1 and 0.0 < b < 2 and -3 < log_f < 0:
        return 0.0
    return -np.inf

def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)

import emcee
pos = soln.x + 1e-4 * np.random.randn(50, 3)
nwalkers, ndim = pos.shape
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(x, y, yerr))
sampler.run_mcmc(pos, 5000, progress=True)

fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
labels = ["m", "b", "log(f)"]
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)
axes[-1].set_xlabel("step number")

tau = sampler.get_autocorr_time()
print(tau)
flat_samples = sampler.get_chain(discard=2500, thin=15, flat=True)
print(flat_samples.shape)

import corner
fig = corner.corner(
    flat_samples, labels=labels, truths=[m_true, b_true, np.log(f_true)])
plt.show()

for i in range(ndim):
    mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
    q = np.diff(mcmc)
    txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
    txt = txt.format(mcmc[1], q[0], q[1], labels[i])
    print (txt)


plt.figure(figsize=(5.8,4.3))
FEH_range = np.arange(-2.7,0.2,0.01)
Gaia = 0.214*FEH_range + 0.88
Murav_lin = 0.34*FEH_range + 1.17
test_RAP = w[0]*FEH_range + w[1]
m_em, b_em = 0.219, 0.863
#Murav_quad = 0.02*(FEH_range**2) + 0.39*FEH_range + 1.19
Cassisi = 0.23*FEH_range + 0.81
plt.plot(FEH_range, Cassisi,label='Cassisi-Sergio')
plt.plot(FEH_range, np.dot(np.vander(FEH_range, 2), [m_ml, b_ml]), ":k", label="ML")
plt.plot(FEH_range, np.dot(np.vander(FEH_range, 2), [m_em, b_em]), ":g", label="emcee")
#plt.plot(FEH_range,test_RAP,label='Raphael')
plt.plot(FEH_range,Gaia,ls='--',label='GaiaCollab+17')
plt.plot(FEH_range,Murav_lin,ls='--',label='Muraveva+18')
#plt.plot(FEH_range,Murav_quad,ls='--',label='Muraveva_quadratic')
plt.scatter(MH,Mv,s=20,color='#CECECE',facecolor='none',label='Data with [M/H]')
plt.scatter(FEH,Mv,s=20,color='black',label='Data with [Fe/H]')
plt.gca().invert_yaxis()
plt.xlabel(r'[Fe/H]',fontsize=11)
plt.ylabel(r'M$_V$',fontsize=11)
plt.legend(loc=3)
plt.tight_layout()
plt.show()

	
