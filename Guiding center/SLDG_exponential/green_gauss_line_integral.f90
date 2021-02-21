    subroutine segment_int(nm,v1,v2,x_base,y_base,idx,idy,centrl,dl,de,sol,s_ab)
    implicit none
    integer,intent(in) :: nm
    real,intent(in) :: v1(1:2),v2(1:2),x_base,y_base,centrl(1:2),dl(1:2),de(1:2)
    integer,intent(in) :: idx,idy
    real,intent(in) :: sol(1:nm)
    real,intent(out) :: s_ab(1:nm)
    real :: x1,y1,x2,y2
    real :: xratio,yratio
    real :: xshift,yshift
    real :: c1_00
    real :: c1_10,c1_01
    real :: c1_20,c1_11,c1_02
    real :: xcc,ycc ,x2m1,y2m1
    real :: xgs(3),ygs(3)
    real :: xg1_2,xg2_2,yg1_2,yg2_2
    real :: xg3_2,yg3_2
    real :: xg1_4,xg2_4,xg3_4
    real :: yg1_4,yg2_4,yg3_4
    real :: xy1,xy2,xy3
    real :: slope,slope_x,slope_y
    real :: green00
    real :: green10,green01,green20,green11,green02
    real :: green30,green21,green12,green03
    real :: green40,green31,green22,green13,green04
    real :: temp2,temp3

    integer :: i_degree

    x1 = (v1(1) -centrl(1) )/dl(1) ;
    y1 = (v1(2) -centrl(2) )/dl(2) ;
    x2 = (v2(1) -centrl(1) )/dl(1) ;
    y2 = (v2(2) -centrl(2) )/dl(2) ;

    xratio = dl(1)/de(1)
    yratio = dl(2)/de(2)

    xshift = ( centrl(1) - x_base ) /de(1)
    yshift = ( centrl(2) - y_base ) /de(2)
    !

    x2m1 = x2-x1
    y2m1 = y2-y1

    if( abs(x2m1)< 1d-13 .and. abs(y2m1)<1d-13 )then
        s_ab(1:nm) = 0.
    else
        if( abs(x2m1) > abs(y2m1) )then
            xcc = ( x1+x2 )*0.5
            ycc = ( y1+y2 )*0.5

            slope = (y2m1)/(x2m1)

            slope_x = slope*x2m1

            green00 = xcc* slope_x

            if(nm>1)then
                xgs(1) = xcc + x2m1 * gau2(1,1)
                xgs(2) = xcc + x2m1 * gau2(2,1)
                ygs(1) = ycc + y2m1 * gau2(1,1)
                ygs(2) = ycc + y2m1 * gau2(2,1)
                xg1_2 = xgs(1)*xgs(1)
                xg2_2 = xgs(2)*xgs(2)
                yg1_2 = ygs(1)*ygs(1)
                yg2_2 = ygs(2)*ygs(2)

                green10 = ( xg1_2 + xg2_2 )*0.25 *slope_x
                green01 = -( yg1_2 + yg2_2 )*0.25*x2m1
                green20 = ( xg1_2*xgs(1) + xgs(2)*xg2_2 )/6. *slope_x
                green11 = ( xg1_2*ygs(1) + xg2_2*ygs(2) )*0.25 *slope_x
                green02 = -( ygs(1)*yg1_2 + ygs(2)*yg2_2 )/6. *x2m1
            endif

            if(nm>3)then

                xgs(1) = xcc + x2m1 * gau3(1,1)
                xgs(2) = xcc + x2m1 * gau3(2,1)
                xgs(3) = xcc + x2m1 * gau3(3,1)
                ygs(1) = ycc + y2m1 * gau3(1,1)
                ygs(2) = ycc + y2m1 * gau3(2,1)
                ygs(3) = ycc + y2m1 * gau3(3,1)

                xg1_2 = xgs(1)*xgs(1)
                xg2_2 = xgs(2)*xgs(2)
                xg3_2 = xgs(3)*xgs(3)
                yg1_2 = ygs(1)*ygs(1)
                yg2_2 = ygs(2)*ygs(2)
                yg3_2 = ygs(3)*ygs(3)

                xy1 = xgs(1)*ygs(1)
                xy2 = xgs(2)*ygs(2)
                xy3 = xgs(3)*ygs(3)
                xg1_4 = xg1_2*xg1_2
                xg2_4 = xg2_2 *xg2_2
                xg3_4 = xg3_2 *xg3_2

                yg1_4 = yg1_2 * yg1_2
                yg2_4 = yg2_2 * yg2_2
                yg3_4 = yg3_2 * yg3_2

                green30=(xg1_4*gau3(1,2)+xg2_4*gau3(2,2)+xg3_4*gau3(3,2))/4.*slope_x
                green21=(xg1_2*xy1*gau3(1,2)+xg2_2*xy2*gau3(2,2)+xg3_2*xy3*gau3(3,2))/3.*slope_x
                green12=-(xy1*yg1_2*gau3(1,2)+xy2*yg2_2*gau3(2,2)+xy3*yg3_2*gau3(3,2))/3.*x2m1
                green03=-(yg1_4*gau3(1,2)+yg2_4*gau3(2,2)+yg3_4*gau3(3,2))/4.*x2m1
                !
                green40=(xg1_4*xgs(1)*gau3(1,2)+xg2_4*xgs(2)*gau3(2,2)+xg3_4*xgs(3)*gau3(3,2))/5.*slope_x
                green31=(xg1_4*ygs(1)*gau3(1,2)+xg2_4*ygs(2)*gau3(2,2)+xg3_4*ygs(3)*gau3(3,2))/4.* slope_x
                green22=(xgs(1)*xg1_2*yg1_2*gau3(1,2)+xgs(2)*xg2_2*yg2_2*gau3(2,2)+xgs(3)*xg3_2*yg3_2*gau3(3,2))/3.*slope_x
                green13=-(yg1_4*xgs(1)*gau3(1,2)+yg2_4*xgs(2)*gau3(2,2)+yg3_4*xgs(3)*gau3(3,2))/4.*x2m1
                green04=-(ygs(1)*yg1_4*gau3(1,2)+ygs(2)*yg2_4*gau3(2,2)+ygs(3)*yg3_4*gau3(3,2))/5.*x2m1
            endif


        else
            xcc = ( x1+x2 )*0.5
            ycc = ( y1+y2 )*0.5

            slope = (x2m1)/(y2m1)
            slope_y = slope* y2m1

            green00 = xcc* y2m1

            if(nm>1)then
                xgs(1) = xcc + x2m1 * gau2(1,1)
                xgs(2) = xcc + x2m1 * gau2(2,1)
                ygs(1) = ycc + y2m1 * gau2(1,1)
                ygs(2) = ycc + y2m1 * gau2(2,1)
                xg1_2 = xgs(1)*xgs(1)
                xg2_2 = xgs(2)*xgs(2)
                yg1_2 = ygs(1)*ygs(1)
                yg2_2 = ygs(2)*ygs(2)

                green10 = ( xg1_2 + xg2_2 )*0.25 *y2m1
                green01 = -( yg1_2 + yg2_2 )*0.25*slope_y
                green20 = ( xg1_2*xgs(1) + xgs(2)*xg2_2  )/6. *y2m1
                green11 = ( xg1_2*ygs(1) + xg2_2*ygs(2) )*0.25 *y2m1
                green02 = -( ygs(1)*yg1_2 + ygs(2)*yg2_2 )/6.*slope_y
            endif

            if(nm>3)then
                xgs(1) = xcc + x2m1 * gau3(1,1)
                xgs(2) = xcc + x2m1 * gau3(2,1)
                xgs(3) = xcc + x2m1 * gau3(3,1)
                ygs(1) = ycc + y2m1 * gau3(1,1)
                ygs(2) = ycc + y2m1 * gau3(2,1)
                ygs(3) = ycc + y2m1 * gau3(3,1)

                xg1_2 = xgs(1)*xgs(1)
                xg2_2 = xgs(2)*xgs(2)
                xg3_2 = xgs(3)*xgs(3)
                yg1_2 = ygs(1)*ygs(1)
                yg2_2 = ygs(2)*ygs(2)
                yg3_2 = ygs(3)*ygs(3)

                xy1 = xgs(1)*ygs(1)
                xy2 = xgs(2)*ygs(2)
                xy3 = xgs(3)*ygs(3)
                xg1_4 = xg1_2*xg1_2
                xg2_4 = xg2_2 *xg2_2
                xg3_4 = xg3_2 *xg3_2

                yg1_4 = yg1_2 * yg1_2
                yg2_4 = yg2_2 * yg2_2
                yg3_4 = yg3_2 * yg3_2

                green30=(xg1_4*gau3(1,2)+xg2_4*gau3(2,2)+xg3_4*gau3(3,2))/4.*y2m1
                green21=(xg1_2*xy1*gau3(1,2)+xg2_2*xy2*gau3(2,2)+xg3_2*xy3*gau3(3,2))/3.*y2m1
                green12=-(xy1*yg1_2*gau3(1,2)+xy2*yg2_2*gau3(2,2)+xy3*yg3_2*gau3(3,2))/3.*slope_y
                green03=-(yg1_4*gau3(1,2)+yg2_4*gau3(2,2)+yg3_4*gau3(3,2))/4.*slope_y
                !
                green40=(xg1_4*xgs(1)*gau3(1,2)+xg2_4*xgs(2)*gau3(2,2)+xg3_4*xgs(3)*gau3(3,2))/5.*y2m1
                green31=(xg1_4*ygs(1)*gau3(1,2)+xg2_4*ygs(2)*gau3(2,2)+xg3_4*ygs(3)*gau3(3,2))/4.*y2m1
                green22=(xgs(1)*xg1_2*yg1_2*gau3(1,2)+xgs(2)*xg2_2*yg2_2*gau3(2,2)+xgs(3)*xg3_2*yg3_2*gau3(3,2))/3.*y2m1
                green13=-(yg1_4*xgs(1)*gau3(1,2)+yg2_4*xgs(2)*gau3(2,2)+yg3_4*xgs(3)*gau3(3,2))/4.*slope_y
                green04=-( ygs(1)*yg1_4*gau3(1,2)+ygs(2)*yg2_4*gau3(2,2)+ygs(3)*yg3_4*gau3(3,2))/5.*slope_y
            endif

        endif

        if(nm ==1 )then
            c1_00 = sol(1)
            s_ab(1) = c1_00 * green00
        elseif(nm == 3)then
            c1_00 = sol(1)
            c1_10 = sol(2)
            c1_01 = sol(3)
            s_ab(1) = c1_00 * green00 + c1_10*green10 + c1_01*green01
            s_ab(2) = xshift *s_ab(1) + xratio*(c1_00*green10 + c1_10*green20 + c1_01*green11)
            s_ab(3) = yshift *s_ab(1) + yratio*(c1_00*green01 + c1_10*green11 + c1_01*green02)
        elseif(nm==6)then
            c1_20 = sol(4)
            c1_11 = sol(5)
            c1_02 = sol(6)
            c1_00 = sol(1) - ( c1_20  +c1_02 )/12.
            c1_10 = sol(2)
            c1_01 = sol(3)

            s_ab(1) =  c1_00 * green00 + c1_10*green10 + c1_01*green01 &
                + c1_20 * green20 + c1_11*green11 + c1_02*green02

            temp2 = c1_00 * green10 + c1_10*green20 + c1_01*green11 &
                + c1_20 * green30 + c1_11*green21 + c1_02*green12
            s_ab(2) =  xshift* s_ab(1)  + xratio*temp2

            temp3 = c1_00 * green01 + c1_10*green11 + c1_01*green02 &
                + c1_20 * green21 + c1_11*green12 + c1_02*green03
            s_ab(3) = yshift* s_ab(1)  + yratio*temp3

            s_ab(4) =  (xshift*xshift-1./12.)* s_ab(1)+ temp2*2.*xshift*xratio &
                +  xratio**2*(c1_00 * green20 + c1_10*green30 + c1_01*green21 &
                + c1_20 * green40 + c1_11*green31 + c1_02*green22)

            s_ab(5) = xshift*yshift*s_ab(1)+ xratio*yshift * temp2 + xshift*yratio * temp3 &
                + xratio*yratio*(c1_00 * green11 + c1_10*green21 + c1_01*green12 &
                + c1_20 * green31 + c1_11*green22 + c1_02*green13)

            s_ab(6) = (yshift*yshift-1./12.)* s_ab(1)+ 2. * yshift*yratio * temp3 &
                + yratio**2*( c1_00 * green02 + c1_10*green12 + c1_01*green03 &
                + c1_20 * green22 + c1_11*green13 + c1_02*green04 )
        endif

    endif

    end subroutine segment_int