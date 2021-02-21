    subroutine get_subfaces(p)
    implicit none
    type(type_face),pointer :: p
    real :: v1(1:2),v2(1:2),v3(1:2)
    real :: a_trans,b_trans,c_trans,d_trans

    real :: xt1,yt1,xf1,yf1
    real :: xt2,yt2,xf2,yf2
    type :: type_extreme
        sequence
        real kexi
        integer itype
    end type

    type(type_extreme) extreme(1:2),extreme_t

    real :: eh,del

    v1(1:2) = p%point_origin%coor(1:2)
    v2(1:2) = p%point_midt%coor(1:2)
    v3(1:2) = p%point_end%coor(1:2)

    call trans(v1,v3,a_trans,b_trans,c_trans,d_trans)

    p%et2(1) = trans_to_ex( a_trans,b_trans,c_trans,v2 )
    p%et2(2) = trans_to_ey( a_trans,b_trans,d_trans,v2 )

    ! for x(\xi)
    p%xa = ( v3(2)-v1(2) )/2.*p%et2(2)/(p%et2(1)**2-1)
    p%xb = ( v3(1)-v1(1) )/2.
    p%xc = ( v3(1)+v1(1) )/2.-p%xa
    ! for y(\xi)
    p%ya = -( v3(1)-v1(1) )/2.*p%et2(2)/(p%et2(1)**2-1)
    p%yb = ( v3(2)-v1(2) )/2.
    p%yc = ( v3(2)+v1(2) )/2. - p%ya


    p%nsubface = 1


    if(  abs(p%xa) <=eps )then
        ! do nothing
    elseif(  abs(p%xa) >eps )then
        eh = -p%xb/(2*p%xa)
        !del =   p%xb**2-4*p%xa*p%xc  
        !if( eh>-1. .and. eh<1. .and. del>eps .and. -del/4./p%xa<-eps  )then
        if( eh>-1. .and. eh<1.    )then
            extreme( p%nsubface )%kexi = eh
            extreme( p%nsubface )%itype = 1
            p%nsubface = p%nsubface + 1
        endif
    endif


    if(  abs(p%ya) <=eps )then
        ! do nothing
    elseif(  abs(p%ya) >eps )then
        eh = -p%yb/(2*p%ya)
        !del =   p%yb**2-4*p%ya*p%yc 
        !if( eh>-1. .and. eh<1.  .and. del>eps .and. -del/4./p%ya<-eps   )then
        if( eh>-1. .and. eh<1.       )then
            extreme( p%nsubface )%kexi = eh
            extreme( p%nsubface )%itype = 2
            p%nsubface = p%nsubface + 1
        endif
    endif
    
    !if( abs(  p%et2(2)  )<=1000*eps .and. p%nsubface>1 )then
    !    print *,io,i,j,p%et2(1:2),p%nsubface
    !    pause
    !endif
    

    
    
    p%ist = 0
    !if( abs(  p%et2(2)  )<=100*eps   )then
    !    p%nsubface = 1
    !    p%ist = 1
    !endif 
    
    !if( abs(  p%et2(2)  )<=100*eps .or. abs(p%et2(1)**2-1) <1000*eps  )then
    !    pause
    !endif
     
    !p%nsubface = 1
    !*******************************************************
    ! Assemble the extreme points of a quadratic-curved face
    !*******************************************************
    if( p%nsubface==1 )then
        ! do nothing!
    elseif( p%nsubface==2 )then
        xt1 = quadratic_function( p%xa,p%xb,p%xc,extreme(1)%kexi )
        yt1 = quadratic_function( p%ya,p%yb,p%yc,extreme(1)%kexi )

        p%point_mid1%xxii = extreme(1)%kexi
 
        if( extreme(1)%itype==1 )then   
            xf1 = xt1
            ! confine y
            call confine( v1(2),v3(2),yt1,yf1 )
        elseif( extreme(1)%itype==2 )then
            ! confine x
            call confine( v1(1),v3(1),xt1,xf1)
            yf1 = yt1
        endif
        
        p%point_mid1%coor(1) = xf1
        p%point_mid1%coor(2) = yf1

        p%point_mid1%id(1) = ceiling( (xf1-xleft)/dx )
        p%point_mid1%id(2) = ceiling( (yf1-ybottom)/dy )
  
        !print *,extreme(1)
         !pause
    elseif( p%nsubface==3 )then
        if( extreme(1)%kexi>extreme(2)%kexi )then
            extreme_t = extreme(1)
            extreme(1) = extreme(2)
            extreme(2) = extreme_t
        endif
        xt1 = quadratic_function( p%xa,p%xb,p%xc,extreme(1)%kexi )
        yt1 = quadratic_function( p%ya,p%yb,p%yc,extreme(1)%kexi )
        xt2 = quadratic_function( p%xa,p%xb,p%xc,extreme(2)%kexi )
        yt2 = quadratic_function( p%ya,p%yb,p%yc,extreme(2)%kexi )
 
        
        if( extreme(1)%itype==1 .and. extreme(2)%itype==2 )then
            xf1 = xt1
            yf2 = yt2
            
            call confine3( v1(2),v3(2),yf2, yt1,yf1 )
            call confine3( v1(1),v3(1),xf1, xt2,xf2 )
        elseif( extreme(1)%itype==2 .and. extreme(2)%itype==1 )then
            yf1 = yt1
            xf2 = xt2 
            
            call confine3( v1(2),v3(2),yf1, yt2,yf2 )
            call confine3( v1(1),v3(1),xf2, xt1,xf1 )            
        endif
        
        p%point_mid1%xxii = extreme(1)%kexi
        p%point_mid1%coor(1) = xf1
        p%point_mid1%coor(2) = yf1
        p%point_mid1%id(1) = ceiling( (xf1-xleft)/dx )
        p%point_mid1%id(2) = ceiling( (yf1-ybottom)/dy )



        p%point_mid2%xxii = extreme(2)%kexi
        p%point_mid2%coor(1) = xf2
        p%point_mid2%coor(2) = yf2
        p%point_mid2%id(1) = ceiling( (xf2-xleft)/dx )
        p%point_mid2%id(2) = ceiling( (yf2-ybottom)/dy )

    endif


    end subroutine get_subfaces
    !****************************************************
    subroutine trans( v1,v3,a_trans,b_trans,c_trans,d_trans )
    implicit none
    real,intent(in) :: v1(1:2),v3(1:2)
    real,intent(out) :: a_trans,b_trans,c_trans,d_trans
    real :: x1,y1,x2,y2

    x1 = v1(1); y1 = v1(2); x2 = v3(1); y2 = v3(2)

    if( abs(x2-x1)< abs(y2-y1) )then
        a_trans = 2.*( (x2-x1)/(y2-y1)/(y2-y1) ) / (   ( (x2-x1)/(y2-y1) )**2  + 1.    )
        b_trans = 2./( y2-y1 )  /  (   ( (x2-x1)/(y2-y1) )**2  + 1.    )
        c_trans = (  (x1-x2)/(y1-y2)*(x1+x2)/(y1-y2)  + (y1+y2)/(y1-y2)  ) /  (   ( (x2-x1)/(y2-y1) )**2  + 1.    )
        d_trans = (  (x1+x2)/(y1-y2)  +  (y1+y2)/(y2-y1)*(x2-x1)/(y2-y1)   )  /  (   ( (x2-x1)/(y2-y1) )**2  + 1.    )
    else
        a_trans = 2./(x2-x1) / (   1.+(   (y2-y1)/(x2-x1)  )**2    )
        b_trans = 2.*(y2-y1)/(x2-x1)/(x2-x1)   / (   1.+(   (y2-y1)/(x2-x1)  )**2    )
        c_trans = (  (x1+x2)/(x1-x2) + (y1-y2)/(x1-x2)*(y1+y2)/(x1-x2)  )  / (   1.+(   (y2-y1)/(x2-x1)  )**2    )
        d_trans = (  (x1+x2)/(x1-x2)*(y1-y2)/(x1-x2)  + (y1+y2)/(x2-x1)  )  / (   1.+(   (y2-y1)/(x2-x1)  )**2    )
    endif
    end subroutine trans
    !**********************************************************
    real function trans_to_ex( a,b,c,v )
    implicit none
    real,intent(in) :: a,b,c,v(1:2)

    trans_to_ex = a*v(1) + b*v(2) + c

    end function trans_to_ex
    !**********************************************************
    real function trans_to_ey( a,b,d,v )
    implicit none
    real,intent(in) :: a,b,d,v(1:2)

    trans_to_ey = b*v(1) - a*v(2) + d

    end function trans_to_ey
    !**********************************************************
    real function quadratic_function(a,b,c,x)
    implicit none
    real,intent(in) :: a,b,c,x

    quadratic_function = a*x**2+b*x+c

    end function quadratic_function