import numpy as np
from matplotlib import pyplot as plt
from math import pi, cos, sin

cen_RA = 255.30249999999995
cen_DE = -30.112361111111113
u = cen_RA      #x-position of the center
v = cen_DE    #y-position of the center
a = (15/60)/np.cos(cen_DE*np.pi/180)     #radius on the x-axis
b = 15/60   #radius on the y-axis
t_rot=pi/4 #rotation angle

t = np.linspace(0, 2*pi, 100)	#std=50
Ell = np.array([a*np.cos(t) , b*np.sin(t)])  
	#u,v removed to keep the same center location
R_rot = np.array([[cos(t_rot) , -sin(t_rot)],[sin(t_rot) , cos(t_rot)]])  
    #2-D rotation matrix

Ell_rot = np.zeros((2,Ell[:].shape[1]))
for i in range(Ell[:].shape[1]):
    Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])

import matplotlib.patches as patch
ellipse = patch.Ellipse((u,v), 2*a, 2*b, angle=0, ls='--', fill=False)
plt.gca().add_artist(ellipse)
plt.plot( u+Ell[0,:] , v+Ell[1,:] )     					#initial ellipse
plt.plot( u+Ell_rot[0,:] , v+Ell_rot[1,:],'darkorange' )	#rotated ellipse
plt.grid(color='lightgray',linestyle='--')
plt.show()