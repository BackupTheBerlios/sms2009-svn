import numpy as np
cimport numpy as np

def calcRHS(np.ndarray[np.float64_t, ndim=2] old_thck,
            np.ndarray[np.float64_t, ndim=2] lsrf,
            np.ndarray[np.float64_t, ndim=2] acab,
            np.ndarray[np.float64_t, ndim=1] sumd,
            Py_ssize_t ew,
            Py_ssize_t ns,
            fc2_3,fc2_4,dt):
        return (old_thck[ew,ns] * (1.0 - fc2_3 * sumd[4])
                - fc2_3 * (old_thck[ew-1,ns]   * sumd[0]      
                                + old_thck[ew+1,ns] * sumd[1]
                                + old_thck[ew,ns-1] * sumd[2]
                                + old_thck[ew,ns+1] * sumd[3]) 
                - fc2_4 * (lsrf[ew,ns]     * sumd[4] 
                                + lsrf[ew-1,ns] * sumd[0]
                                + lsrf[ew+1,ns] * sumd[1]
                                + lsrf[ew,ns-1] * sumd[2]
                                + lsrf[ew,ns+1] * sumd[3]) 
                + acab[ew,ns] * dt)

def fillMatrix(matrix,
	       np.ndarray[np.int64_t,   ndim=2] mask,
	       np.ndarray[np.float64_t, ndim=1] sumd,
	       Py_ssize_t ew,
               Py_ssize_t ns):
        matrix[int(mask[ew,ns]-1),int(mask[ew-1,ns]-1)] = sumd[0]
        matrix[int(mask[ew,ns]-1),int(mask[ew+1,ns]-1)] = sumd[1]
        matrix[int(mask[ew,ns]-1),int(mask[ew,ns-1]-1)] = sumd[2]
        matrix[int(mask[ew,ns]-1),int(mask[ew,ns+1]-1)] = sumd[3]
        matrix[int(mask[ew,ns]-1),int(mask[ew,ns]-1)]   = 1.0 + sumd[4]

def findSums(np.ndarray[np.float64_t, ndim=2] diffu,
	     np.ndarray[np.float64_t, ndim=1] sumd,
	     Py_ssize_t ewm,
	     Py_ssize_t ew,
	     Py_ssize_t nsm,
	     Py_ssize_t ns,
	     np.float64_t fc2_1,
	     np.float64_t fc2_5):
        sumd[0] = fc2_1 * ((diffu[ewm,nsm] + diffu[ewm,ns]))
        sumd[1] = fc2_1 * ((diffu[ew,nsm]  + diffu[ew,ns]))
        sumd[2] = fc2_5 * ((diffu[ewm,nsm] + diffu[ew,nsm]))
        sumd[3] = fc2_5 * ((diffu[ewm,ns]  + diffu[ew,ns]))
        sumd[4] = - (sumd[0] + sumd[1] + sumd[2] + sumd[3])

def thckcrit(np.ndarray[np.float64_t, ndim=2] ca, np.float64_t cb):
	
    cdef np.int64_t count
    cdef Py_ssize_t i,j

    count = 0

    cdef int ii = ca.shape[0]
    cdef int jj = ca.shape[1]

    for i in range(ii):
       for j in range(jj):
          if ca[i,j]>0.0: count +=1

    if (count > 0) or (cb > 0.0):
       return True
    else:
       return False

def glide_maskthck(np.ndarray[np.float64_t, ndim=2] thck,
                   np.ndarray[np.float64_t, ndim=2] acab,
                   np.ndarray[np.int64_t,   ndim=2] mask,
                   grid):
    
    # Calculates the contents of the mask array.
    mask[:] = 0
    covtot  = 0

    empty = True

    cdef Py_ssize_t ns
    cdef Py_ssize_t ew

    for ns in range(grid.ny):
        for ew in range(grid.nx):
            if thckcrit(thck[max(0,ew-1):min(grid.nx,ew+1),max(0,ns-1):min(grid.ny,ns+1)],acab[ew,ns]):
                covtot += 1
                mask[ew,ns] = covtot 
                if (empty): empty = False

    return covtot, empty

def regrid(np.ndarray[np.float64_t, ndim=2] new_thck,
           np.ndarray[np.int64_t,   ndim=2] mask,
           np.ndarray[np.float64_t, ndim=1] ans,
           mainGrid):

    cdef Py_ssize_t ns
    cdef Py_ssize_t ew

    for ns in range(mainGrid.ny):
        for ew in range(mainGrid.nx):
            if mask[ew,ns] != 0: new_thck[ew,ns] = ans[mask[ew,ns]-1]
