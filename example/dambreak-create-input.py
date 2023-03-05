# --- GraSPH Input Demo Python Script
# Script shows how the example hdf5 input dambreak.h5 can be generated using 
# Python and h5py. Each dataset is gzip compressed.h5. Note that Python stores
# arrays in row-major, whereas Fortran stores arrays in column-major.

import h5py
import numpy as np

# constants
dim = 3
dx = 0.5
rho = 1000

# Opening new file to write
with h5py.File('example/dambreak.h5','w') as f:

    # Create real group to write real particle data to
    g_r = f.create_group("real")

    # Preparing real particles' coordinate data
    xp = 50 # number of particles in x direction
    yp = 25 #                        y
    zp = 50 #                        z
    x = []
    for i in range(xp):
        for j in range(yp):
            for k in range(zp):
                x.append([(i + 0.5)*dx, (j + 0.5)*dx, (k + 0.5)*dx])

    nreal = len(x)

    # write number of real particles as attribute
    g_r.attrs.create("n", data=nreal, dtype="i")

    # write coordinate data
    g_r.create_dataset("x", 
                       data=x, 
                       compression="gzip", 
                       dtype="f8")
    
    # write initial velocity for each particle
    g_r.create_dataset("v", 
                       data=np.zeros(shape=(nreal, dim)), # zero initial velocity
                       compression="gzip", 
                       dtype="f8")
    
    # write initial pressure for each particle
    g_r.create_dataset("p", 
                       data=np.zeros(shape=(nreal)), # zero initial pressure
                       compression="gzip", 
                       dtype="f8")

    # write initial density for each particle
    g_r.create_dataset("rho", 
                       data=np.full(shape=(nreal), fill_value=rho), 
                       compression="gzip", 
                       dtype="f8")

    # write type index for each particle
    g_r.create_dataset("type", 
                       data=np.ones(shape=(nreal)), 
                       compression="gzip", 
                       dtype="i")

    # write global particle ID
    g_r.create_dataset("ind", 
                       data=np.arange(1, nreal+1), # Fortran indexing from 1
                       compression="gzip", 
                       dtype="i")

    # Create virtual group to write virtual particle data to
    g_v = f.create_group("virt")

    # Preparing coordinate data and type data
    # Particles are given type IDs based on whether they are a bottom, top, north, south, east, or west face, or corner
    xp = 150 # number of particles in x direction
    yp = 25 #                         y
    zp = 100 #                        z
    nlayer = 4
    x = []
    itype = []
    for i in range(-nlayer, xp+nlayer):
        for j in range(-nlayer, yp+nlayer):
            for k in range(-nlayer, zp+nlayer):
                if i < 0 or j < 0 or k < 0 or i >= xp or j >= yp or k >= zp:
                    x.append([(i + 0.5)*dx, (j + 0.5)*dx, (k + 0.5)*dx])
                    if (i < 0 and j < 0) or (i < 0 and j >= yp) or (i >= xp and j < 0) or (i >= xp and j >= yp):
                        itype.append(-5)
                    elif j < 0 or j >= yp:
                        itype.append(-4)
                    elif i < 0 or i >= xp:
                        itype.append(-3)
                    elif k >= zp:
                        itype.append(-2)
                    elif k < 0:
                        itype.append(-1)

    nvirt = len(x)

    # write number of virtual particles as attribute
    g_v.attrs.create("n", data=nvirt, dtype="i")

    # write writing virtual position data
    g_v.create_dataset("x", 
                       data=x, 
                       compression="gzip", 
                       dtype="f8")
    
    # write initial velocity for each particle
    g_v.create_dataset("v", 
                       data=np.zeros(shape=(nvirt, dim)), # zero initial velocity
                       compression="gzip", 
                       dtype="f8")
    
    # write initial pressure for each particle
    g_v.create_dataset("p", 
                       data=np.zeros(shape=(nvirt)), # zero initial pressure
                       compression="gzip", 
                       dtype="f8")

    # write initial density for each particle
    g_v.create_dataset("rho", 
                       data=np.full(shape=(nvirt), fill_value=rho), 
                       compression="gzip", 
                       dtype="f8")

    # write type index for each particle
    g_v.create_dataset("type", 
                       data=itype, 
                       compression="gzip", 
                       dtype="i")

    # write global particle ID
    g_v.create_dataset("ind", 
                       data=np.arange(1, nvirt+1), # Fortran indexing from 1
                       compression="gzip", 
                       dtype="i")