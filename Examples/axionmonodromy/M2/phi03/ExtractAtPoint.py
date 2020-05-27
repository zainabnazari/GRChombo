import yt
from yt import derived_field
import numpy as np
#from numpy.linalg import inv
#from scipy.interpolate import interp1d
#from scipy.optimize import fsolve
import time
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#yt.enable_parallelism()

#Loading dataset
loc = '/home/fmuia/GRChombo/Examples/Tanh/M014/phi03/T014_phi03_P_00**'
ds = yt.load(loc)

#from mpi4py import MPI
#comm = MPI.COMM_WORLD

#initial_location
location =  np.array([0,0,0])
#location =  ds[0].domain_right_edge/2
#variable = 'VofPhi'
variable = 'chi'
time_data = []
var_data = []
CycleData = []
#filename = 'VofPhiT02.txt'
filename = 'chi03T014.txt'

last_current_time = 0
proper_time = 0
for i in ds:
  print("Now at ", i.current_time)

  # find the center of the BH
  #value, location = i.find_min("chi")
  #center = [float(location[0]), float(location[1]), float(location[2])]
  #print ("New center ", center)

  point = i.point(location)
  lapse = point["lapse"]
  dt = i.current_time - last_current_time
  last_current_time = i.current_time
  proper_time += lapse * dt
  var_out = point[variable]
#  if(comm.rank==0) :
  datafile=open(filename,'a')
  datafile.write("%f    %.8f \n" % (proper_time, var_out))
  datafile.close()
