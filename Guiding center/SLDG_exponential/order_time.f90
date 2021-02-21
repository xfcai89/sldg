    subroutine order_time
    !*******************************************************************************
    !
    !   Purpose  : temporal accuracy and order of DG solution
    !
    !
    !*******************************************************************************
    implicit none
    integer :: kk0,kk1
    real :: den
    integer :: lx,ly
    integer :: inf_x,inf_y
    real :: xrg,yrg
    real :: rr1,rr2,rr3,error1,error2,error3


    error1=0.0
    error2=0.0
    error3=0.0

    if(kkkk==0)then
        do kk0=1,nx
            do kk1 = 1,ny
                do lx=1,6
                    do ly =1,6
                        ref(kk0,kk1,lx,ly) = &
                            polynomial2d( element(kk0,kk1,0 )%umodal(1:n_moment), x(kk0)+xg(lx)*dx,x(kk0),dx,y(kk1)+xg(ly)*dy,y(kk1),dy,n_moment )
                    enddo
                enddo
            enddo
        enddo
    endif

    !call data_save
    !call data_read

    do kk0=1,nx
        do kk1=1,ny
            do lx=1,6
                do ly=1,6
                    den=abs( ref(kk0,kk1,lx,ly)  &
                        - polynomial2d( element(kk0,kk1,0)%umodal(1:n_moment), x(kk0)+xg(lx)*dx,x(kk0),dx,y(kk1)+xg(ly)*dy,y(kk1),dy,n_moment ) )
                    error1=error1+den*wg(lx)*wg(ly)*dx*dy
                    error2=error2+den*den*wg(lx)*wg(ly)*dx*dy
                enddo
            enddo
            do lx=1,6
                do ly=1,6
                    xrg =  x(kk0)+xg(lx)*dx
                    yrg =  y(kk1)+xg(ly)*dy
                    error3=max(error3,&
                        abs( ref(kk0,kk1,lx,ly)  &
                        - polynomial2d( element(kk0,kk1,0)%umodal(1:n_moment),xrg,x(kk0),dx,yrg,y(kk1),dy,n_moment )  ))

                enddo
            enddo

        enddo
    enddo
    error1=error1/( (xright-xleft)*(ytop-ybottom) )
    error2=sqrt(error2/( (xright-xleft)*(ytop-ybottom) ) )
    if(kkkk.eq.1) write(133,103) nx,ny,error1,error2,error3
    write(*,*) error1,error2,error3
    if(kkkk.gt.1) then
        rr1=log(er111/error1)/log( real(kkkk-1)/real(kkkk) )
        rr2=log(er222/error2)/log( real(kkkk-1)/real(kkkk) )
        rr3=log(er333/error3)/log( real(kkkk-1)/real(kkkk) )
        write(133,102) nx,ny,error1,rr1,error2, rr2,error3, rr3
        write(*,*) nx,ny,rr1,rr2,rr3
    endif
    er111=error1
    er222=error2
    er333=error3

    close(111)


111 format(4(1x,f12.4))
102 format(i6,'*',i6,1x,3('&',1x, es12.2e2,1x,'&' 1x,f8.2 ,1x),'&',i8,'&',i8,'\\',1x,'\hline')
103 format(i6,'*',i6,1x,3('&',1x,es12.2E2,1x,'&',1x) ,'&',i8,'&',i8,'\\',1x,'\hline')

    !***********************************************************************************************
    !open(121,file='temp_1705.plt')
    !write(121,*)'zone    ','i=',nx*6,',    j=',ny*6
    !
    !DO  J=1,NY
    !    do ly = 1,6
    !        DO  I=1 ,NX
    !            do lx = 1,6
    !                xrg =  x(i)+xg(lx)*dx
    !                yrg =  y(j)+xg(ly)*dy
    !                WRITE(121,*) xrg,yrg,  ref( i,j,lx,ly)  &
    !                    - polynomial( element(i,j)%umodal(1:n_moment),xrg,x(i),dx,yrg,y(j),dy,n_moment )
    !            enddo
    !        enddo
    !    enddo
    !enddo
    !
    !
    !CLOSE(121)

    end subroutine order_time