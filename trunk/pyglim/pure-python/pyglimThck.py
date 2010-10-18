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
        self.tol = 1e-13

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

    def _findSums(self,diffu,sumd,ewm,ew,nsm,ns):
        sumd[0] = self.fc2_1 * ((diffu[ewm,nsm] + diffu[ewm,ns]))
        sumd[1] = self.fc2_1 * ((diffu[ew,nsm]  + diffu[ew,ns]))
        sumd[2] = self.fc2_5 * ((diffu[ewm,nsm] + diffu[ew,nsm]))
        sumd[3] = self.fc2_5 * ((diffu[ewm,ns]  + diffu[ew,ns]))
        sumd[4] = - (sumd[0] + sumd[1] + sumd[2] + sumd[3])

    def _applyULBCs(self,matrix,old_thck,new_thck,mask,calc_rhs,
                    rhs,ans):
        # lower and upper BC
        for ew in range(self.mainGrid.nx):
            ns = 0
            if mask[ew,ns] != 0: 
                matrix[int(mask[ew,ns]-1),int(mask[ew,ns]-1)] = 1.0
                if calc_rhs: rhs[mask[ew,ns]-1] = old_thck[ew,ns] 
                ans[mask[ew,ns]-1] = new_thck[ew,ns]

            ns=self.mainGrid.ny-1
            if mask[ew,ns] != 0:
                matrix[int(mask[ew,ns]-1),int(mask[ew,ns]-1)] = 1.0
                if calc_rhs: rhs[mask[ew,ns]-1] = old_thck[ew,ns] 
                ans[mask[ew,ns]-1] = new_thck[ew,ns]

    def _applyLRBCs(self,matrix,old_thck,new_thck,mask,calc_rhs,
                    rhs,ans):

        # left and right BC (N.B. Periodic BCs are not implemented)
        for ns in range(1,self.mainGrid.ny-1):
            ew=0
            if mask[ew,ns] != 0:
                matrix[int(mask[ew,ns]-1),int(mask[ew,ns]-1)] = 1.0
                if calc_rhs: rhs[mask[ew,ns]-1] = old_thck[ew,ns] 
                ans[mask[ew,ns]-1] = new_thck[ew,ns]

            ew=self.mainGrid.nx-1
            if mask[ew,ns] != 0:
                matrix[int(mask[ew,ns]-1),int(mask[ew,ns]-1)] = 1.0
                if calc_rhs: rhs[mask[ew,ns]-1] = old_thck[ew,ns] 
                ans[mask[ew,ns]-1] = new_thck[ew,ns]

    def _calcRHS(self,old_thck,lsrf,acab,sumd,ew,ns):
        return (old_thck[ew,ns] * (1.0 - self.fc2_3 * sumd[4])
                - self.fc2_3 * (old_thck[ew-1,ns]   * sumd[0]      
                                + old_thck[ew+1,ns] * sumd[1]
                                + old_thck[ew,ns-1] * sumd[2]
                                + old_thck[ew,ns+1] * sumd[3]) 
                - self.fc2_4 * (lsrf[ew,ns]     * sumd[4] 
                                + lsrf[ew-1,ns] * sumd[0]
                                + lsrf[ew+1,ns] * sumd[1]
                                + lsrf[ew,ns-1] * sumd[2]
                                + lsrf[ew,ns+1] * sumd[3]) 
                + acab[ew,ns] * self.dt)

    def _fillMatrix(self,matrix,mask,sumd,ew,ns):
        matrix[int(mask[ew,ns]-1),int(mask[ew-1,ns]-1)] = sumd[0]
        matrix[int(mask[ew,ns]-1),int(mask[ew+1,ns]-1)] = sumd[1]
        matrix[int(mask[ew,ns]-1),int(mask[ew,ns-1]-1)] = sumd[2]
        matrix[int(mask[ew,ns]-1),int(mask[ew,ns+1]-1)] = sumd[3]
        matrix[int(mask[ew,ns]-1),int(mask[ew,ns]-1)]   = 1.0 + sumd[4]

    def _solveSystem(self,matrix,rhs,ans):
        info, iter, relres = its.pcg(matrix, rhs, ans, 1e-10, 100)
        #print info, iter, relres

    def _regrid(self,new_thck,mask,ans):

        for ns in range(self.mainGrid.ny):
            for ew in range(self.mainGrid.nx):
                if mask[ew,ns] != 0: new_thck[ew,ns] = ans[mask[ew,ns]-1]

    def _thckEvolve(self,rhs,ans,mask,calc_rhs,old_thck,new_thck,diffu,lsrf,acab):

        matrix = pysp.ll_mat(self.totpts,self.totpts)

        # Boundary Conditions ---------------------------------------------------------------
        self._applyULBCs(matrix,old_thck,new_thck,self.mask,calc_rhs,rhs,ans)
        self._applyLRBCs(matrix,old_thck,new_thck,self.mask,calc_rhs,rhs,ans)

        # Ice body
        sumd = np.zeros(5,dtype=np.float)

        for ns in range(1,self.mainGrid.ny-1):
            for ew in range(1,self.mainGrid.nx-1):
                if mask[ew,ns] != 0:
                    self._findSums(diffu,sumd,ew-1,ew,ns-1,ns)
                    # Matrix elements
                    self._fillMatrix(matrix,self.mask,sumd,ew,ns)
                    # RHS vector
                    if calc_rhs:
                        rhs[mask[ew,ns]-1] = self._calcRHS(old_thck,lsrf,acab,sumd,ew,ns)
                    ans[mask[ew,ns]-1] = new_thck[ew,ns]

        # Solve system
        self._solveSystem(matrix,rhs,ans)

        # Rejig the solution onto a 2D array
        self._regrid(new_thck,self.mask,ans)

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
                if not self.linear: oldthck2 = ISM.thck.copy()
                if not self.linear and p>0:
                    pg.stagvarb(ISM.thck,self.stagthck,self.mainGrid)
                    pg.df_field_2d_staggered(ISM.usrf,self.dusrfdew,self.dusrfdns,self.mainGrid)
                    pg.df_field_2d_staggered(ISM.thck,self.dthckdew,self.dthckdns,self.mainGrid)

                # Calculate diffusivity
                self._calcDiffu(ISM.flwa,ISM.diffu)

                # Evolve thickness
                self._thckEvolve(rhs,ans,self.mask,first_p,oldthck,ISM.thck,ISM.diffu,ISM.lsrf,acab)

                first_p = False
                if not self.linear: 
                    residual = np.maximum(np.abs(ISM.thck-oldthck2))
                    if residual < self.tol: break
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
                
