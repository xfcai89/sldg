    subroutine output
    implicit none
    integer :: lx,ly,ic,jc
    real :: xrg,yrg

    open(125,file='temp_16.plt')
    write(125,*)'zone    ','i=',nx*3,',    j=',ny*3

    DO  J=1,NY
        do jc = 1,3
            DO  I=1 ,NX
                do ic = 1,3
                    xrg =  x(i)+ (ic-2.) *dx*0.25
                    yrg =  y(j)+ (jc-2.) *dy*0.25                
                    WRITE(125,*) xrg,yrg,  polynomial2d( element(i,j,io)%umodal(1:n_moment),xrg,x(i),dx,yrg,y(j),dy,n_moment )
                enddo
            enddo
        enddo
    enddo


    CLOSE(125) 




    end subroutine output