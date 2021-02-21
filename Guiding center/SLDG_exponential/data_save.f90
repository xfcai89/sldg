    subroutine data_save
    implicit none
    integer :: kk0,kk1

    integer :: lx,ly  

    open(unit=120,file='restart',form='unformatted', status='unknown')
    rewind 120
    write(120)  (((( ref(kk0,kk1,lx,ly),kk0=1,nx),kk1=1,ny),lx=1,6),ly=1,6)
    close(120)

    return
    end subroutine data_save
    !**************************************************************************
    subroutine data_read
    implicit none
    integer :: irest
    integer :: kk0,kk1

    integer :: lx,ly     

    irest = 1
    if(irest.eq.1) then
        open(unit=120,file='restart',form='unformatted',status='unknown')
        rewind 120
        read(120) (((( ref(kk0,kk1,lx,ly),kk0=1,nx),kk1=1,ny),lx=1,6),ly=1,6)
        close(120)
    end if

    end subroutine data_read