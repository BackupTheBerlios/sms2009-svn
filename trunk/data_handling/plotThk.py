import pylab
import numpy
import Scientific.IO.NetCDF

time=-1

# load file
aFile = Scientific.IO.NetCDF.NetCDFFile('fenscan.2000a.nc','r')

# create a mask based on ice thicknesses
mask = numpy.array(aFile.variables['thk'][time,:,:])<1.e-5

# create array holding ice thickness
thk = numpy.ma.array(aFile.variables['thk'][time,:,:],mask=mask)

pylab.imshow(thk,origin='lower')
pylab.show()
