    subroutine  id_get17(xy,dxy,id17,ixy)
    implicit none
    real,intent(in) :: xy,dxy
    integer,intent(in) :: ixy
    integer,intent(out) :: id17
    integer :: i00

    integer :: id1703
    if(ixy == 1)then
        id1703 = id_get( (xy-xleft)/dxy )
        do i00 = id1703 +1-2, id1703 + 1+2
            if( abs( xy - xgrid(i00) ) .le. dxy + 1d-14 .and. xy-xgrid(i00) .le. 0. )then !must be 1d-14
                id17 = i00-1
            endif
        enddo
    else
        id1703 = id_get( (xy-yleft)/dxy )
        !print *,id1703,xy
        do i00 = id1703 +1-2, id1703 +1 +2
            !print*,ygrid(i00)
            if( abs( xy - ygrid(i00) ) .le. dxy + 1d-14 .and. xy-ygrid(i00) .le. 0. )then
                id17 = i00-1
                !print *,123
            endif
        enddo
        !print *,id17
    endif


    end subroutine id_get17

    integer function id_get(x16)
    !*******************************************************************************
    !
    !   Purpose : getting the location index of point in a Eulerian cell
    !
    !   Notice  : Background cell and Eulerian cell is different in the implementation
    !             Background grid : x(i) + 1d-14
    !             Eulerian cell : x(i)
    !
    !*******************************************************************************
    implicit none
    real,intent(in) :: x16

    if(x16<0.)then
        id_get = int(x16)
    else
        id_get = int(x16) + 1
    endif

    end function id_get