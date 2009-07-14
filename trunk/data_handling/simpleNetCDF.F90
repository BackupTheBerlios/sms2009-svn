!> simple program loading a netCDF file
!! \author Magnus Hagdorn
!! \date July 2009
program simple
  use netcdf
  implicit none
  
  integer status, ncid

  status = nf90_open("example.nc",NF90_NOWRITE,ncid)
  call checkError(status,__FILE__,__LINE__)
  
  status = nf90_close(ncid)
  call checkError(status,__FILE__,__LINE__)

contains
  subroutine checkError(status,fname,line)
    implicit none
    integer, intent(in) :: status                   !< the return value of the netCDF call
    character(len=*), optional, intent(in) :: fname !< name of the file where error occured
    integer, optional, intent(in) :: line           !< line number where error occured
    
    if (status/=NF90_NOERR) then
       write(*,'(a)',advance='no') 'Error '
       if (present(fname) .and. present(line)) then
          write(*,'(a,i0,a)',advance='no') '('//trim(fname)//',',line,') '
       end if
       write(*,'(a)') nf90_strerror(status)
       stop
    end if
  end subroutine checkError
end program simple
