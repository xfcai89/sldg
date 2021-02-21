    subroutine order_DG
    implicit none
    !  real,allocatable :: uct(:,:)
    real :: den

    real :: xrg,yrg
    real :: error1,error2,error3,rr1,rr2,rr3
    integer :: i,j,n1,n2

    error1=0.0
    error2=0.0
    error3=0.0

    do j=1,nel_y
        do i=1,nel_x
            do n2=1,n_g
                do n1=1, n_g
                    !den=abs( sin( gx(n1,n2,i,j) + gy(n1,n2,i,j) - 2.*  tprint )  &
                    !    -  elem(i,j)%psi( n1,n2 )    )

                    !den=abs( sin( gx(n1,n2,i,j)  -   tprint )  &
                    ! -  elem(i,j)%psi( n1,n2 )    )
                    
                    den=abs( exact(  gx(n1,n2,i,j) , gy(n1,n2,i,j), tprint)  &
                        -  elem(i,j)%psi( n1,n2 )    )
                    error1=error1+den*w_g(n1)* w_g(n2)*dx*dy
                    error2=error2+den*den*w_g(n1)*w_g(n2)*dx*dy
                enddo
            enddo


            do n2=1,n_g
                do n1=1,n_g
                    !error3=max(error3,abs( sin( gx(n1,n2,i,j) + gy(n1,n2,i,j) -   2.*tprint )  &
                    !    -  elem(i,j)%psi( n1,n2 )    )    )
                    !    error3=max(error3,abs( sin( gx(n1,n2,i,j)  -   tprint )  &
                    ! -  elem(i,j)%psi( n1,n2 )    )    )
                    !error3=max(error3,abs( sin( gy(n1,n2,i,j)  -   tprint )  &
                    !  -  elem(i,j)%psi( n1,n2 )    )    )
                    
                     error3=max(error3,abs( exact(  gx(n1,n2,i,j) , gy(n1,n2,i,j), tprint)   &
                        -  elem(i,j)%psi( n1,n2 )    )    )                   
                enddo
            enddo
            !   error3=max(error3,den)

        enddo
    enddo
    error1=error1/( (xright-xleft)*(yright-yleft) )
    error2=sqrt(error2/( (xright-xleft)*(yright-yleft) ))
    if(kkkk.eq.1) write(123,103) nel_x,nel_y,error1,error2,error3
    write(*,*) error1,error2,error3
    if(kkkk.gt.1) then
        rr1=log(er11/error1)/log(real(nel_x)/real(norder(kkkk-1)  ))
        rr2=log(er22/error2)/log(real(nel_x)/real(norder(kkkk-1)  ))
        rr3=log(er33/error3)/log(real(nel_x)/real(norder(kkkk-1)  ))
        write(123,102) nel_x,nel_y,error1,rr1,error2, rr2,error3, rr3
        write(*,*) nel_x,nel_y,rr1,rr2,rr3
    endif
    er11=error1
    er22=error2
    er33=error3


    close(111)


111 format(4(1x,f12.4))
102 format(i6,'*',i6,1x,3('&',1x, es12.2e2,1x,'&' 1x,f8.2 ,1x),'&',i8,'&',i8,'\\',1x,'\hline')
103 format(i6,'*',i6,1x,3('&',1x,es12.2E2,1x,'&',1x) ,'&',i8,'&',i8,'\\',1x,'\hline')

    end subroutine order_DG
