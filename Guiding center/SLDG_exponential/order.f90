    subroutine order
    !*******************************************************************************
    !
    !   Purpose  : accuracy and order of average
    !
    !*******************************************************************************
    implicit none
    integer :: kk0,kk1
    real,allocatable :: uct(:,:)
    real :: den
    integer :: lx,ly
    real :: xrg,yrg
    real :: error1,error2,error3
    real :: rr1,rr2,rr3

    allocate( uct( 1-nghost:nx+nghost , 1-nghost:nx+nghost) )

    !open(111,file='_'//char(48)//char(48+kkkk)//'.plt')
    write(111,*)'zone    ','i=',nx,',    j=',ny


    DO  J=1,NY
        DO  I=1 ,NX
            WRITE(111,*) X(I),Y(J), element(i,j,0)%umodal(1)
        enddo
    enddo

    CLOSE(111)

    !open(222,file='y=0_'//char(48)//char(48+kkkk)//'.plt')
    DO  I=1 ,NX
        WRITE(222,*) X(I), element(i,1,0)%umodal(1)
    enddo


    CLOSE(222)

    do  i=1,nx
        do  j=1,ny
            uct(i,j)=0.0
            do lx=1,6
                do ly=1,6
                    uct(i,j)=uct(i,j)+exact(x(i)+xg(lx)*dx ,y(j)+xg(ly)*dy,0.) *wg(lx)*wg(ly)
                enddo
            enddo
        enddo
    enddo

    error1=0.0
    error2=0.0
    error3=0.0

    do kk0=1,nx
        do kk1=1,ny
            den=abs( element(kk0,kk1,0)%umodal(1) -uct(kk0,kk1) )
            error1=error1+den*dx*dy
            error2=error2+den*den*dx*dy
            error3=max(error3,den)

        enddo
    enddo
    error1=error1/( (xright-xleft)*(ytop-ybottom) )
    error2=sqrt(error2/( (xright-xleft)*(ytop-ybottom) ) )
    if(kkkk.eq.1) write(101,103) nx,ny,error1,error2,error3
    write(*,*) error1,error2,error3
    if(kkkk.gt.1) then
        !rr1=log(er11/error1)/log(real(nx)/real(norder(kkkk-1) ))
        !rr2=log(er22/error2)/log(real(nx)/real(norder(kkkk-1) ))
        !rr3=log(er33/error3)/log(real(nx)/real(norder(kkkk-1) ))

        rr1=log(er11/error1)/log( real(kkkk)/real(kkkk-1) )
        rr2=log(er22/error2)/log( real(kkkk)/real(kkkk-1) )
        rr3=log(er33/error3)/log( real(kkkk)/real(kkkk-1) )

        write(101,102) nx,ny,error1,rr1,error2, rr2,error3, rr3
        write(*,*) nx,ny,rr1,rr2,rr3
    endif
    er11=error1
    er22=error2
    er33=error3

    deallocate(uct)

    close(111)


111 format(4(1x,f12.4))
102 format(i6,'*',i6,1x,3('&',1x, es12.2e2,1x,'&' 1x,f8.2 ,1x),'&',i8,'&',i8,'\\',1x,'\hline')
103 format(i6,'*',i6,1x,3('&',1x,es12.2E2,1x,'&',1x) ,'&',i8,'&',i8,'\\',1x,'\hline')


    return
    end subroutine order