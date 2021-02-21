    subroutine green_p2qc_gauss3_outer(v1,v2,x_base,y_base,idx,idy,s_ab,v11,v22,v33,sol)
    implicit none
    real,intent(in) :: v1(1:2),v2(1:2),x_base,y_base,v11(1:2),v22(1:2),v33(1:2)
    integer,intent(in) :: idx,idy
    real,intent(out) :: s_ab(1:n_moment)
    real,intent(in) :: sol(1:n_moment)
    real :: x1,y1,x2,y2
    real :: xshift,yshift

  
    real :: green_00,green_10,green_01,green_20,green_11,green_02,green_30,green_21,green_12,green_03
    real :: green_40,green_31,green_22,green_13,green_04
    real :: c1_00,c1_10,c1_01,c1_20,c1_11,c1_02
 

    real :: temp2,temp3
    !*****
    real :: v_1(1:2),v_2(1:2),v_3(1:2)
    real :: a_trans,b_trans,c_trans,d_trans
    real :: e_x2,e_y2
    real :: ac1,bc1,cc1
    real :: ac2,bc2,cc2
    real :: v(1:2)
    real :: exx1,exx2,dx_xi
    real :: xig(1:3)
    real :: xi_gauss1,xi_gauss2,xi_gauss3
    real :: dxi_gauss1,dxi_gauss2,dxi_gauss3
    real :: yi_gauss1,yi_gauss2,yi_gauss3
    real :: dyi_gauss1,dyi_gauss2,dyi_gauss3

    x1 = (v1(1) -x(idx) )/dx ;
    y1 = (v1(2) -y(idy) )/dy ;
    x2 = (v2(1) -x(idx) )/dx ;
    y2 = (v2(2) -y(idy) )/dy ;
    xshift = ( x(idx) - x_base ) /dx
    yshift = ( y(idy) - y_base ) /dy

    c1_20 = sol(4)
    c1_11 = sol(5)
    c1_02 = sol(6)
    c1_00 = sol(1) -1./12.*( c1_20  +c1_02 )
    c1_10 = sol(2)
    c1_01 = sol(3)

    !**************************************
    v_1(1) = ( v11(1) -x(idx) )/dx ;
    v_1(2) = ( v11(2) -y(idy) )/dy ;
    v_2(1) = ( v22(1) -x(idx) )/dx ;
    v_2(2) = ( v22(2) -y(idy) )/dy ;
    v_3(1) = ( v33(1) -x(idx) )/dx ;
    v_3(2) = ( v33(2) -y(idy) )/dy ;
    call trans( v_1(1:2),v_3(1:2),a_trans,b_trans,c_trans,d_trans )
    e_x2 = trans_to_ex( a_trans ,b_trans , c_trans, v_2(1:2)   )
    e_y2 = trans_to_ey( a_trans ,b_trans , d_trans, v_2(1:2)  )

    ac1 = ( v_3(2)-v_1(2) )/2.*e_y2/( e_x2**2 -1. )
    bc1 = ( v_3(1)-v_1(1) )/2.
    cc1 = ( v_1(1)+v_3(1) )/2. - ( v_3(2)-v_1(2) )/2.*e_y2/( e_x2**2 -1. )

    ac2 = -( v_3(1)-v_1(1) )/2.*e_y2/( e_x2**2 -1. )
    bc2 = ( v_3(2) - v_1(2) )/2.
    cc2 = ( v_1(2)+v_3(2) )/2. + ( v_3(1)-v_1(1) )/2.*e_y2/( e_x2**2 -1. )

    v(1)=x1; v(2) =y1
    exx1 = trans_to_ex( a_trans ,b_trans , c_trans, v(1:2) )
    v(1)=x2; v(2) =y2
    exx2 = trans_to_ex( a_trans ,b_trans , c_trans, v(1:2) )
    dx_xi = exx2 -exx1

    xig(1) = (exx1+exx2)*0.5 + dx_xi * gau3(1,1)
    xig(2) = (exx1+exx2)*0.5 + dx_xi * gau3(2,1)
    xig(3) = (exx1+exx2)*0.5 + dx_xi * gau3(3,1)

    xi_gauss1 = ac1*xig(1)*xig(1) + bc1*xig(1) + cc1
    dxi_gauss1 = 2.*ac1*xig(1) + bc1
    xi_gauss2 = ac1*xig(2)*xig(2) + bc1*xig(2) + cc1
    dxi_gauss2 = 2.*ac1*xig(2) + bc1
    xi_gauss3 = ac1*xig(3)*xig(3) + bc1*xig(3) + cc1
    dxi_gauss3 = 2.*ac1*xig(3) + bc1

    yi_gauss1 = ac2*xig(1)*xig(1) + bc2*xig(1) + cc2
    dyi_gauss1 = 2.*ac2*xig(1) + bc2
    yi_gauss2 = ac2*xig(2)*xig(2) + bc2*xig(2) + cc2
    dyi_gauss2 = 2.*ac2*xig(2) + bc2
    yi_gauss3 = ac2*xig(3)*xig(3) + bc2*xig(3) + cc2
    dyi_gauss3 = 2.*ac2*xig(3) + bc2

    green_00 =  xi_gauss1 * dyi_gauss1*gau3(1,2) + xi_gauss2 * dyi_gauss2*gau3(2,2)+ xi_gauss3* dyi_gauss3*gau3(3,2)
    green_10 = ( xi_gauss1 * xi_gauss1 * dyi_gauss1*gau3(1,2) + xi_gauss2 *xi_gauss2 * dyi_gauss2*gau3(2,2)  &
        + xi_gauss3* xi_gauss3 * dyi_gauss3*gau3(3,2)    )*0.5
    green_01 = -( yi_gauss1*yi_gauss1*dxi_gauss1*gau3(1,2) + yi_gauss2*yi_gauss2*dxi_gauss2*gau3(2,2)        &
        +  yi_gauss3*yi_gauss3*dxi_gauss3*gau3(3,2)  )*0.5
    green_20 = ( xi_gauss1*xi_gauss1 * xi_gauss1 * dyi_gauss1*gau3(1,2)  &
        + xi_gauss2*xi_gauss2 *xi_gauss2 * dyi_gauss2*gau3(2,2) &
        + xi_gauss3*xi_gauss3 *xi_gauss3 * dyi_gauss3*gau3(3,2)   )/3.
    green_11 = (  xi_gauss1 * xi_gauss1*yi_gauss1 * dyi_gauss1 *gau3(1,2) &
        + xi_gauss2 *xi_gauss2*yi_gauss2 * dyi_gauss2*gau3(2,2)   &
        + xi_gauss3 *xi_gauss3*yi_gauss3 * dyi_gauss3*gau3(3,2)   )*0.5
    green_02 = -(  yi_gauss1*yi_gauss1*yi_gauss1*dxi_gauss1*gau3(1,2) &
        + yi_gauss2*yi_gauss2*yi_gauss2*dxi_gauss2*gau3(2,2)  &
        + yi_gauss3*yi_gauss3*yi_gauss3*dxi_gauss3*gau3(3,2)    )/3.
    !
    green_30 = ( xi_gauss1*xi_gauss1 * xi_gauss1 * xi_gauss1* dyi_gauss1 *gau3(1,2)  &
        + xi_gauss2*xi_gauss2 *xi_gauss2*xi_gauss2 * dyi_gauss2*gau3(2,2) &
        + xi_gauss3*xi_gauss3 *xi_gauss3*xi_gauss3 * dyi_gauss3*gau3(3,2)   )*0.25
    green_21 = ( xi_gauss1*xi_gauss1 * xi_gauss1*yi_gauss1 * dyi_gauss1*gau3(1,2)  &
        + xi_gauss2*xi_gauss2 *xi_gauss2*yi_gauss2 * dyi_gauss2*gau3(2,2)   &
        + xi_gauss3*xi_gauss3 *xi_gauss3*yi_gauss3* dyi_gauss3*gau3(3,2)   )/3.
    green_12 = -( xi_gauss1*yi_gauss1*yi_gauss1*yi_gauss1*dxi_gauss1*gau3(1,2) &
        + xi_gauss2*yi_gauss2*yi_gauss2*yi_gauss2*dxi_gauss2*gau3(2,2)   &
        + xi_gauss3*yi_gauss3*yi_gauss3*yi_gauss3*dxi_gauss3*gau3(3,2)    )/3.
    !
    green_03 =  -(  yi_gauss1*yi_gauss1*yi_gauss1*yi_gauss1*dxi_gauss1*gau3(1,2) &
        + yi_gauss2*yi_gauss2*yi_gauss2*yi_gauss2*dxi_gauss2*gau3(2,2)     &
        + yi_gauss3*yi_gauss3*yi_gauss3*yi_gauss3*dxi_gauss3*gau3(3,2)  )*0.25
    !
    green_40 = ( xi_gauss1*xi_gauss1 * xi_gauss1 * xi_gauss1* xi_gauss1* dyi_gauss1*gau3(1,2)  &
        + xi_gauss2*xi_gauss2 *xi_gauss2*xi_gauss2*xi_gauss2 * dyi_gauss2*gau3(2,2)  &
        + xi_gauss3*xi_gauss3 *xi_gauss3*xi_gauss3*xi_gauss3 * dyi_gauss3*gau3(3,2)   )*0.2
    green_31 = ( xi_gauss1*xi_gauss1 * xi_gauss1 * xi_gauss1*yi_gauss1* dyi_gauss1*gau3(1,2)  &
        + xi_gauss2*xi_gauss2 *xi_gauss2*xi_gauss2*yi_gauss2* dyi_gauss2*gau3(2,2) &
        + xi_gauss3*xi_gauss3 *xi_gauss3*xi_gauss3*yi_gauss3* dyi_gauss3*gau3(3,2)   )*0.25
    green_22 = ( xi_gauss1*xi_gauss1 * xi_gauss1*yi_gauss1*yi_gauss1* dyi_gauss1*gau3(1,2)  &
        + xi_gauss2*xi_gauss2 *xi_gauss2*yi_gauss2*yi_gauss2* dyi_gauss2*gau3(2,2)    &
        + xi_gauss3*xi_gauss3 *xi_gauss3*yi_gauss3*yi_gauss3* dyi_gauss3*gau3(3,2)  )/3.
    !
    green_13 = -( yi_gauss1*yi_gauss1*yi_gauss1*yi_gauss1*xi_gauss1*dxi_gauss1*gau3(1,2) &
        + yi_gauss2*yi_gauss2*yi_gauss2*yi_gauss2*xi_gauss2*dxi_gauss2*gau3(2,2)  &
        + yi_gauss3*yi_gauss3*yi_gauss3*yi_gauss3*xi_gauss3*dxi_gauss3*gau3(3,2) )*0.25
    !
    green_04 = -( yi_gauss1*yi_gauss1*yi_gauss1*yi_gauss1*yi_gauss1*dxi_gauss1*gau3(1,2) &
        + yi_gauss2*yi_gauss2*yi_gauss2*yi_gauss2*yi_gauss2*dxi_gauss2*gau3(2,2)        &
        + yi_gauss3*yi_gauss3*yi_gauss3*yi_gauss3*yi_gauss3*dxi_gauss3*gau3(3,2)    )*0.2

    s_ab(1) =  (c1_00 * green_00 + c1_10*green_10 + c1_01*green_01 &
        + c1_20 * green_20 + c1_11*green_11 + c1_02*green_02)*dx_xi
 
    temp2 = (c1_00 * green_10 + c1_10*green_20 + c1_01*green_11 &
        + c1_20 * green_30 + c1_11*green_21 + c1_02*green_12)*dx_xi
    s_ab(2) =  xshift* s_ab(1) + temp2

    temp3 = (c1_00 * green_01 + c1_10*green_11 + c1_01*green_02 &
        + c1_20 * green_21 + c1_11*green_12 + c1_02*green_03 )*dx_xi
    s_ab(3) = yshift* s_ab(1)   +temp3

    s_ab(4) =  (xshift*xshift-1./12. )* s_ab(1) &
        + temp2*2.*xshift &
        +  (c1_00 * green_20 + c1_10*green_30 + c1_01*green_21 &
        + c1_20 * green_40 + c1_11*green_31 + c1_02*green_22)*dx_xi

    s_ab(5) = xshift*yshift*s_ab(1) &
        + yshift * temp2 + xshift * temp3 &
        + (c1_00 * green_11 + c1_10*green_21 + c1_01*green_12 &
        + c1_20 * green_31 + c1_11*green_22 + c1_02*green_13)*dx_xi

    s_ab(6) = (yshift*yshift-1./12.)* s_ab(1) &
        + 2. * yshift * temp3 &
        +  (c1_00 * green_02 + c1_10*green_12 + c1_01*green_03 &
        + c1_20 * green_22 + c1_11*green_13 + c1_02*green_04)*dx_xi
    !**************************************

    end subroutine green_p2qc_gauss3_outer