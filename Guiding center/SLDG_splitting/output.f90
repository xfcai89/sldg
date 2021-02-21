    subroutine output
    implicit none
    integer :: ic,jc
    real :: xrg,yrg


    open(121,file='temp_16.plt')
    write(121,*)'zone    ','i=',nel_x*n_g,',    j=',nel_y*n_g

    DO  J=1,Nel_Y
        do jc = 1,n_g
            DO  I=1 ,Nel_X
                do ic = 1,n_g
                    xrg =  x(i)+ x_g(ic)*dx
                    yrg =  y(j)+ x_g(jc)*dy
                    WRITE(121,*) xrg,yrg,  elem(i,j)%psi( ic,jc )
                enddo
            enddo
        enddo
    enddo

    CLOSE(121)

    end subroutine output

    subroutine output_temp
    implicit none
    integer :: lx,ly,ic,jc
    real :: xrg,yrg

    open(121,file='temp_17.plt')
    write(121,*)'zone    ','i=',nel_x*3,',    j=',nel_y*3

    DO  J=1,Nel_Y
        do jc = 1,3
            DO  I=1 ,Nel_X
                do ic = 1,3
                    xrg =  x(i)+ (ic-2.) *dx*0.25
                    yrg =  y(j)+ (jc-2.) *dy*0.25
                    WRITE(121,*) xrg,yrg,  polynomial2d( temp11(i,j,1:n_moment_2d),xrg,x(i),dx,yrg,y(j),dy,n_moment_2d )
                enddo
            enddo
        enddo
    enddo

    CLOSE(121)

    end subroutine output_temp