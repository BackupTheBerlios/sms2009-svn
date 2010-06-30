import pyglim
import pyglimThck
import numpy as np
import pylab
import Scientific.IO.NetCDF as nc

nx = 31
ny = 31

ism = pyglim.IceSheetModel2D()
ism.setThicknessSolver(pyglimThck.SIAsolver)
a=np.zeros([nx,ny])
a[:]=0.3

f=nc.NetCDFFile("testpyglim.nc","w")
f.createDimension("x0",nx)
f.createDimension("y0",ny)
f.createDimension("x1",nx-1)
f.createDimension("y1",ny-1)
f.createDimension("t",None)

thckvar=f.createVariable("thck","d",("t","x0","y0"))
diffvar=f.createVariable("diff","d",("t","x1","y1"))

for i in range(300):
    print 'Time: %i, Volume: %f'%(i*100,np.sum(ism.thck)*(25e8))
    ism.timestep(a,10)
    thckvar[i,:,:] = ism.thck
    diffvar[i,:,:] = ism.diffu

ism.finish()
f.close()
