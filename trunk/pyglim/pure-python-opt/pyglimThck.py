# SIA thickness solver for pyglim

import numpy as np
import pyglimGrid as pg
import pysparse.spmatrix as pysp
import pysparse.itsolvers as its
import Scientific.IO.NetCDF as nc
import pyglimConstants as pygc
import math

diffuconst = -2.0*math.pow((pygc.rhoi*pygc.grav),pygc.gn) / (pygc.gn + 2)

def thckcrit(ca,cb):

    if np.logical_or(ca > 0.0,cb > 0.0).any():
       return True
    else:
       return False

def glide_maskthck(thck,acab,mask,grid):
    
    # Calculates the contents of the mask array.
    mask[:] = 0
    covtot  = 0

    empty = True

    for ns in range(grid.ny):
        for ew in range(grid.nx):
            if thckcrit(thck[max(0,ew-1):min(grid.nx,ew+1),max(0,ns-1):min(grid.ny,ns+1)],acab[ew,ns]):
                covtot += 1
                mask[ew,ns] = covtot 
                if (empty): empty = False

    return covtot, empty

class SIAsolver(object):

    def __init__(self,mainGrid,veloGrid,output=False,dt=10.0):
        self.mainGrid = mainGrid
        self.veloGrid = veloGrid

        # Declare some internal arrays
        self.stagthck = pg.newArray(self.veloGrid) # Thickness on velo grid
        self.dusrfdew = pg.newArray(self.veloGrid) # Horizontal derivs of usrf
        self.dusrfdns = pg.newArray(self.veloGrid) 
        self.dthckdew = pg.newArray(self.veloGrid) # Horizontal derivs of thck
        self.dthckdns = pg.newArray(self.veloGrid)
        self.mask     = pg.newArray(self.mainGrid,dtype=np.int32) # Mask

        # Parameters
        self.pmax = 50     # Maximum number of Picard iterations
        self.linear = True # Linear model
        self.alpha = 0.5   # Still not sure what this does
        self.fc2_3 = (1.0 - self.alpha) / self.alpha
        self.fc2_4 =  1.0 / self.alpha

        # Set timestep
        self.setTimestep(dt)

        # Diagnostic output
        self.output = output
        if self.output:
            self.f=nc.NetCDFFile("SIAsolver.nc","w")
            self.f.createDimension("x0",self.mainGrid.nx)
            self.f.createDimension("y0",self.mainGrid.ny)
            self.f.createDimension("x1",self.veloGrid.nx)
            self.f.createDimension("y1",self.veloGrid.ny)
            self.f.createDimension("t",None)

            self.maskvar = self.f.createVariable("mask","i",("t","x0","y0"))
            self.stagvar = self.f.createVariable("stagthck","d",("t","x1","y1"))
            self.du1var  = self.f.createVariable("dusrfdew","d",("t","x1","y1"))
            self.du2var  = self.f.createVariable("dusrfdns","d",("t","x1","y1"))
            self.dt1var  = self.f.createVariable("dthckdew","d",("t","x1","y1"))
            self.dt2var  = self.f.createVariable("dthckdns","d",("t","x1","y1"))

            self.nccount = 0

    def setTimestep(self,dt):
        self.dt = dt
        self.fc2_1 = self.alpha * dt / (2.0 * self.mainGrid.dx * self.mainGrid.dx)
        self.fc2_5 = self.alpha * dt / (2.0 * self.mainGrid.dy * self.mainGrid.dy)

    def _calcDiffu(self,flwa,diffu):

        diffu[:] = np.where(self.stagthck!=0.0,
                                 diffuconst * flwa * np.power(self.stagthck,pygc.p4) 
                                 * np.power(np.sqrt(np.power(self.dusrfdew,2)
                                                    +np.power(self.dusrfdns,2)),pygc.p2),
                                 0.0)

    def _thckEvolve(self,rhs,ans,calc_rhs,old_thck,new_thck,diffu,lsrf,acab):

        matrix = pysp.ll_mat(self.totpts,self.totpts)

        # Boundary Conditions ---------------------------------------------------------------
        # lower and upper BC
        for bc in [0,self.mainGrid.ny-1]:

            rawIdx = np.nonzero( self.mask[:,bc])
            idx = self.mask[ rawIdx, bc].ravel()-1
            matrix[idx,idx] = 1.0
            if calc_rhs:
                rhs[idx] = old_thck[ rawIdx, bc].ravel()
            ans[idx] = new_thck[ rawIdx, bc].ravel()

        # left and right BC (N.B. Periodic BCs are not implemented)
        for bc in [0,self.mainGrid.nx-1]:
            rawIdx = np.nonzero( self.mask[bc,:])
            idx = self.mask[ bc,rawIdx].ravel()-1
            matrix[idx,idx] = 1.0
            if calc_rhs:
                rhs[idx] = old_thck[ rawIdx, bc].ravel()
            ans[idx] = new_thck[ bc,rawIdx].ravel()
            
        # Ice body
        sumd = np.zeros(5,dtype=np.float)

        rawIdx = np.nonzero(self.mask[1:-1,1:-1])
        idx  = (self.mask[1:-1,1:-1])[rawIdx] -1

        idx2 = (self.mask[0:-2,1:-1])[rawIdx] -1
        sum1 = (self.fc2_1 * (diffu[0:-1,0:-1] + diffu[0:-1,1:]))
        matrix.put(sum1.ravel(), idx, idx2)

        idx2 = (self.mask[2:,1:-1])[rawIdx] -1
        sum2 = (self.fc2_1 * (diffu[1:,0:-1] + diffu[1:,1:]))
        matrix.put(sum2.ravel(), idx, idx2)

        idx2 = (self.mask[1:-1,0:-2])[rawIdx] -1
        sum3 = (self.fc2_5 * (diffu[0:-1,0:-1] + diffu[1:,0:-1]))
        matrix.put(sum3.ravel(), idx, idx2)

        idx2 = (self.mask[1:-1,2:])[rawIdx] -1
        sum4 = (self.fc2_5 * (diffu[0:-1,1:] + diffu[1:,1:]))
        matrix.put(sum4.ravel(), idx, idx2)

        sum5 = -(sum1 + sum2 + sum3 + sum4)
        matrix.put(1.+sum5.ravel(),idx,idx)

        # RHS vector
        if calc_rhs:
            rhs[idx] = (old_thck[1:-1,1:-1] * (1.0 - self.fc2_3 * sum5)
                        - self.fc2_3 * (old_thck[0:-2,1:-1]   * sum1
                                        + old_thck[2:,1:-1]   * sum2
                                        + old_thck[1:-1,0:-2] * sum3
                                        + old_thck[1:-1,2:]   * sum4)
                        - self.fc2_4 * (lsrf[1:-1,1:-1]       * sum5
                                        + lsrf[0:-2,1:-1]     * sum1
                                        + lsrf[2:,1:-1]       * sum2
                                        + lsrf[1:-1,0:-2]     * sum3
                                        + lsrf[1:-1,2:]       * sum4)
                        + acab[1:-1,1:-1] * self.dt)[rawIdx]
             

        ans[idx] = (new_thck[1:-1,1:-1])[rawIdx]

        # Solve system
        info, iter, relres = its.pcg(matrix, rhs, ans, 1e-10, 100)
        print info, iter, relres

        # Rejig the solution onto a 2D array
        rawIdx = np.nonzero(self.mask)
        idx  = self.mask[np.nonzero(self.mask)]-1
        new_thck.flat[idx] = ans[idx]

        # Remove negatives
        new_thck[:] = np.maximum(0.0,new_thck)

    def timestep(self,ISM,acab):

        # Check inputs are correct format
        assert isinstance(acab,np.ndarray),"Mass-balance (acab) must be an array"
        assert len(acab.shape)==2,"Mass-balance (acab) must be 2D"
        assert (acab.shape[0]==ISM.mainGrid.nx and 
                acab.shape[1]==ISM.mainGrid.ny), "Mass-balance (acab) must be correct size"

        # Staggered fields and derivatives
        pg.stagvarb(ISM.thck,self.stagthck,self.mainGrid)
        pg.df_field_2d_staggered(ISM.usrf,self.dusrfdew,self.dusrfdns,self.mainGrid)
        pg.df_field_2d_staggered(ISM.thck,self.dthckdew,self.dthckdns,self.mainGrid)

        # Do masking
        self.totpts, empty = glide_maskthck(ISM.thck,acab,self.mask,self.mainGrid)

        # RHS and unknown of matrix solution
        rhs = np.zeros(self.totpts,dtype=np.float)
        ans = np.zeros(self.totpts,dtype=np.float)

        # Timestep-dependent constants

        if empty:
            ISM.thck[:] = np.maximum(0.0,ISM.thck + acab * self.dt)
        else:
            first_p = True
            oldthck = ISM.thck.copy()

            # Picard iterations
            for p in range(self.pmax):

                # For non-linear computations only
                if not self.linear: oldthck2 = thck.copy()
                if not self.linear and p>0:
                    pg.stagvarb(ISM.thck,self.stagthck,self.mainGrid)
                    pg.df_field_2d_staggered(ISM.usrf,self.dusrfdew,self.dusrfdns,self.mainGrid)
                    pg.df_field_2d_staggered(ISM.thck,self.dthckdew,self.dthckdns,self.mainGrid)

                # Calculate diffusivity
                self._calcDiffu(ISM.flwa,ISM.diffu)

                # Evolve thickness
                self._thckEvolve(rhs,ans,first_p,oldthck,ISM.thck,ISM.diffu,ISM.lsrf,acab)

                first_p = False
                if not self.linear: 
                    residual = np.maximum(np.abs(thck-oldthck2))
                    if residual < tol: break
                else:
                    break

        ISM.usrf = ISM.topg + ISM.thck

        # Write some output
        if self.output:
            self.maskvar[self.nccount,:,:] = self.mask
            self.du1var [self.nccount,:,:] = self.dusrfdew
            self.du2var [self.nccount,:,:] = self.dusrfdns
            self.dt1var [self.nccount,:,:] = self.dthckdew
            self.dt2var [self.nccount,:,:] = self.dthckdns
            self.stagvar[self.nccount,:,:] = self.stagthck
            self.nccount += 1
        
    def finish(self):
        if self.output: self.f.close()
                
