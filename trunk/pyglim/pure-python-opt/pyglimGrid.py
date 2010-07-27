
import numpy as np

class Grid(object):
    def __init__(self,nx,ny,dx,dy,x0=0.0,y0=0.0):
        self.nx = nx
        self.ny = ny
        self.dx = dx
        self.dy = dy
        self.x0 = x0
        self.y0 = y0

def newArray(grid,nz=0,dtype=np.float):
    if nz==0:
        return np.zeros([grid.nx,grid.ny],dtype=dtype)
    else:
        return np.zeros([grid.nx,grid.ny,nz],dtype=dtype)

   
def stagvarb(input,output,grid):
    """Interpolates a two-dimensional field from the main grid to the velocity grid"""

    output[0:grid.nx-1,0:grid.ny-1] = (input[1:grid.nx,    0:grid.ny-1] 
                                       + input[0:grid.nx-1,1:grid.ny] 
                                       + input[1:grid.nx,  1:grid.ny]   
                                       + input[0:grid.nx-1,0:grid.ny-1]) / 4.0



def df_field_2d_staggered(f, out_dfdx, out_dfdy, grid):
    xfactor = 1/(2.0*grid.dx)
    yfactor = 1/(2.0*grid.dy)
    
    out_dfdx[0:grid.nx-1,0:grid.ny-1] = (f[1:grid.nx,0:grid.ny-1] + f[1:grid.nx,1:grid.ny] -
                                         f[0:grid.nx-1,0:grid.ny-1] - f[0:grid.nx-1,1:grid.ny])*xfactor
    out_dfdy[0:grid.nx-1,0:grid.ny-1] = (f[0:grid.nx-1,1:grid.ny] + f[1:grid.nx,1:grid.ny] -
                                         f[0:grid.nx-1,0:grid.ny-1] - f[1:grid.nx,0:grid.ny-1])*yfactor


