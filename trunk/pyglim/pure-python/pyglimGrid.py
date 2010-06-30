
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


def dfdx_2d_stag(f, i, j, delta):
    return (f[i+1, j] + f[i+1, j+1] - f[i, j] - f[i, j+1])/(2.0*delta) 

def dfdy_2d_stag(f, i, j, delta):
    return (f[i, j+1] + f[i+1, j+1] - f[i,j] - f[i+1, j])/(2.0*delta)

def df_field_2d_staggered(f, out_dfdx, out_dfdy, grid):
    # For now, we'll use the function calls defined above.
    # Later on we might want to refactor?    
    for x in range(grid.nx - 1): # We go to nx - 1 because we're using a staggered grid
       for y in range(grid.ny - 1):
          out_dfdx[x,y] = dfdx_2d_stag(f, x, y, grid.dx)
          out_dfdy[x,y] = dfdy_2d_stag(f, x, y, grid.dy)

