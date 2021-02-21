    subroutine get_norm
    implicit none
    real :: norm_1
    real :: xrg,yrg
    integer :: lx,ly

    !**********************************************************Mass
    norm_1 = 0.D0

    do i = 1 , nx
        do j = 1 , ny
            !
            do lx = 1,6
                do ly = 1,6
                    xrg =  x(i)+xg(lx)*dx
                    yrg =  y(j)+xg(ly)*dy
                    norm_1 = norm_1 + polynomial2d( element(i,j,io)%umodal(1:n_moment),xrg,x(i),dx,yrg,y(j),dy,n_moment )*wg(lx)*wg(ly)
                    
                enddo
            enddo

        enddo
    enddo
    norm_1 = norm_1 * dx*dy
    if(irelative==1)then
        r_norm1 = norm_1
    endif
    irelative = 0
    write(111,*) time,(norm_1-r_norm1)/r_norm1

    end subroutine get_norm