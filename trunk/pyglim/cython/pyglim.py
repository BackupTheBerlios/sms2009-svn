# Top-level ISM module
import numpy as np
import pyglimCGrid as pgc

# Class for 2D ice sheet model
class IceSheetModel2D(object):

    def __init__(self,nx=31,ny=31,dx=50e3,dy=50e3):
        
        # Set up two grids
        self.mainGrid = pgc.Grid(nx,  ny,  dx,dy)
        self.veloGrid = pgc.Grid(nx-1,ny-1,dx,dy,dx/2.0,dy/2.0)

        # Create main arrays
        self.thck  = pgc.newArray(self.mainGrid) # Thickness
        self.topg  = pgc.newArray(self.mainGrid) # Basal topg
        self.usrf  = pgc.newArray(self.mainGrid) # Surface elevation
        self.lsrf  = pgc.newArray(self.mainGrid) # Lower surface elevation
        self.flwa  = pgc.newArray(self.veloGrid) # Vertically-integrated flow-law
        self.diffu = pgc.newArray(self.veloGrid) # Diffusivity
        self.uvel  = pgc.newArray(self.veloGrid) # x-velocity (depth-mean)
        self.vvel  = pgc.newArray(self.veloGrid) # y-velocity (depth-mean)

        # For now, set flowlaw to a constant
        self.flwa[:] = 1e-16

        # Solvers
        self.thckSolver = None
        self.tempSolver = None

    def setThicknessSolver(self,thckSolver):
        """Set which thickness solver to use"""

        # Checks need to be added to make sure this
        # is a suitable class

        self.thckSolver = thckSolver(self.mainGrid,self.veloGrid)

    def setTemperatureSolver(self,tempSolver):
        """Set which temperature solver to use"""

        # Checks need to be added to make sure this
        # is a suitable class

        self.tempSolver = tempSolver(self.mainGrid,self.veloGrid)

    def timestep(self,acab,nt=1):
        """Perform timesteps"""

        for t in range(nt): 
            if self.thckSolver: self.thckSolver.timestep(self,acab)
            if self.tempSolver: self.tempSolver.timestep(self,acab)

    def finish(self):
        if self.thckSolver: self.thckSolver.finish()
        if self.tempSolver: self.tempSolver.finish()
