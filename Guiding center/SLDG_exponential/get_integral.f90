    subroutine get_integral
    implicit none
    com_mass(1:nx,1:ny) = 0.
    do i = 1 , nx
        do j = 1 , ny
            element_star(i,j,io)%vertex5 => nodex_star(i,j)
            element_star(i,j,io)%vertex6 => nodey_star(i+1,j)
            element_star(i,j,io)%vertex7 => nodex_star(i,j+1)
            element_star(i,j,io)%vertex8 => nodey_star(i,j)
            element_star(i,j,io)%vertex9 => nodec_star(i,j )
            call green(umod_t(i,j,1:n_moment) )
        enddo
    enddo

    !open(121,file='ave16.plt')
    !write(121,*)'zone    ','i=',nx,',    j=',ny
    !do i = 1 , nx
    !    do j = 1 , ny
    !        if( iadditive(io) ==0 )then
    !            WRITE(121,*) X(I),Y(J),element(i,j,0)%umodal(1) -com_mass(i,j)
    !            if( abs(element(i,j,0)%umodal(1) - com_mass(i,j))>5e-3)then
    !                print *,element(i,j,0)%umodal(1) - com_mass(i,j),i,j ,io
    !                pause
    !            endif
    !        else
    !            WRITE(121,*) X(I),Y(J),element(i,j,io)%umodal(1) -com_mass(i,j)
    !            if( abs(element(i,j,io)%umodal(1) - com_mass(i,j))>5e-3)then
    !                print *,element(i,j,io)%umodal(1) - com_mass(i,j),i,j ,io
    !                pause
    !            endif
    !        endif
    !
    !    enddo
    !enddo
    !
    !close(121)
    !pause

    do i = 1 , nx
        do j = 1 , ny
            element(i,j,io+1)%umodal(1:n_moment) = umod_t(i,j,1:n_moment)
        enddo
    enddo

    end subroutine get_integral

    !