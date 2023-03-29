# this python interface uses the DPStokes repo to compute pair mobilities
import numpy as np
from mobility_matrix import *
import scipy.io

device = 'cpu'
domType = 'DPSC'

eta = 1/6/np.sqrt(np.pi)
radP = 2.0
xmin = 0.0; xmax = 76.8
ymin = 0.0; ymax = 76.8
zmin = 0.0; zmax = 19.2

solver = FCMJoint(device); nP = 2
M = np.zeros((3*nP, 3*nP), dtype = np.double)

xsep = 4
xPs = np.zeros(3*nP, dtype = np.double)
xPs = np.array([xmax/2+xsep/2, ymax/2, 4, xmax/2-xsep/2, ymax/2, 7], dtype = np.double)

solver.Initialize(numberParticles=nP, hydrodynamicRadius=radP, kernType=0,
                  domType=domType, has_torque=False,
                  xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
                  viscosity=eta, optInd=0, ref=False)

M, NormMat = mobility_matrix(np.reshape(xPs, (-1, 3*nP)), solver)
M = np.reshape(M, (3*nP,3*nP))
solver.Clean()

if domType == 'DPBW':
  outname = 'pair_mobility_bw.mat'
elif domType == 'DPSC':
  outname = 'pair_mobility_sc.mat'
  
scipy.io.savemat(outname, dict(M=M, xPs=xPs))

