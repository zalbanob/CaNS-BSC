! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_initgrid
  use mod_param, only:pi,is_gridpoint_natural_channel
  use mod_types
  implicit none
  private
  public initgrid
  contains
  !---------------- CUSTOM_GRID ----------------
  subroutine read_zf_from_file(n, lz, zf)
    integer, intent(in) :: n
    real(rp), intent(in) :: lz
    real(rp), intent(inout) :: zf(0:n+1)
    character(len=*), parameter :: fname = 'zgrid.txt'
    integer :: u, ios, count, idx
    real(rp) :: y, dyplus
    open(newunit=u, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'CUSTOM_FILE: cannot open zgrid.txt'
    count = 0
    do
      read(u,*,iostat=ios) idx, y, dyplus
      if (ios /= 0) exit
      if (count == 0) then
        zf(0) = y
      elseif (count <= n) then
        zf(count) = y
      end if
      count = count + 1
    end do
    close(u)
    if (count /= n+1) stop 'CUSTOM_FILE: need exactly n+1 lines of (i,y,dy+)'
  end subroutine read_zf_from_file
  !-------------- END CUSTOM_GRID --------------
  subroutine initgrid(gtype,n,gr,lz,dzc,dzf,zc,zf)
    !
    ! initializes the non-uniform grid along z
    !
    implicit none
    integer, parameter :: CLUSTER_TWO_END              = 1, &
                          CLUSTER_ONE_END              = 2, &
                          CLUSTER_ONE_END_R            = 3, &
                          CLUSTER_MIDDLE               = 4, &
                          HYBRID                       = 5, &
                          CUSTOM_FILE                  = 6
    integer , intent(in ) :: gtype,n
    real(rp), intent(in ) :: gr,lz
    real(rp), intent(out), dimension(0:n+1) :: dzc,dzf,zc,zf
    real(rp) :: z0
    integer :: k
    procedure (), pointer :: gridpoint => null()
    !
    ! step 1) determine coordinates of cell faces zf
    !
    zf(0) = 0.
    if (gtype == CUSTOM_FILE) then
      call read_zf_from_file(n, lz, zf)
    else
      select case(gtype)
      case(CLUSTER_TWO_END)
        gridpoint => gridpoint_cluster_two_end
      case(CLUSTER_ONE_END)
        gridpoint => gridpoint_cluster_one_end
      case(CLUSTER_ONE_END_R)
        gridpoint => gridpoint_cluster_one_end_r
      case(CLUSTER_MIDDLE)
        gridpoint => gridpoint_cluster_middle
      case(HYBRID)
        gridpoint => gridpoint_hybrid
      case default
        gridpoint => gridpoint_cluster_two_end
      end select
      if(.not.is_gridpoint_natural_channel) then
        do k=1,n
          z0  = (k-0.)/(1.*n)
          call gridpoint(gr,z0,zf(k))
        end do
      else
        do k=1,n
          call gridpoint_natural(k,n,zf(k))
        end do
      end if
      zf(1:n) = zf(1:n)*lz
    end if
    !
    ! step 2) determine grid spacing between faces dzf
    !
    do k=1,n
      dzf(k) = zf(k)-zf(k-1)
    end do
    dzf(0  ) = dzf(1)
    dzf(n+1) = dzf(n)
    !
    ! step 3) determine grid spacing between centers dzc
    !
    do k=0,n
      dzc(k) = .5*(dzf(k)+dzf(k+1))
    end do
    dzc(n+1) = dzc(n)
    !
    ! step 4) compute coordinates of cell centers zc and faces zf
    !
    zc(0)    = -dzc(0)/2.
    zf(0)    = 0.
    do k=1,n+1
      zc(k) = zc(k-1) + dzc(k-1)
      zf(k) = zf(k-1) + dzf(k)
    end do
  end subroutine initgrid
  !
  ! grid stretching functions
  ! see e.g., Fluid Flow Phenomena -- A Numerical Toolkit, by P. Orlandi
  !           Pirozzoli et al. JFM 788, 614–639 (commented)
  !
  subroutine gridpoint_cluster_two_end(alpha,z0,z)
    !
    ! clustered at the two sides
    !
    implicit none
    real(rp), intent(in) :: alpha,z0
    real(rp), intent(out) :: z
    if(alpha > epsilon(0._rp)) then
      z = 0.5*(1.+tanh((z0-0.5)*alpha)/tanh(alpha/2.))
      !z = 0.5*(1.+erf( (z0-0.5)*alpha)/erf( alpha/2.))
    else
      z = z0
    end if
  end subroutine gridpoint_cluster_two_end
  subroutine gridpoint_cluster_one_end(alpha,z0,z)
    !
    ! clustered at the lower side
    !
    implicit none
    real(rp), intent(in) :: alpha,z0
    real(rp), intent(out) :: z
    if(alpha > epsilon(0._rp)) then
      z = 1.0*(1.+tanh((z0-1.0)*alpha)/tanh(alpha/1.))
      !z = 1.0*(1.+erf( (z0-1.0)*alpha)/erf( alpha/1.))
    else
      z = z0
    end if
  end subroutine gridpoint_cluster_one_end
  !---------------- HYBRID ----------------
  subroutine gridpoint_hybrid(alpha,z0,z)
    !
    ! constant grid size over obstacle height followed by two-side stretching of the grid
    !
    implicit none
    real(rp), intent(in) :: alpha,z0
    real(rp), intent(out) :: z
    integer  :: n1, n2, n
    real(rp) :: z1, z2, l0
    n1 = 31   ! number of grid points in the constant grid size section
    n2 = 161  ! number of grid points in the variable grid size section
    n  = n1 + n2
    l0 = 0.1/1.05 ! height of the obstacle normalized by the height of the domain
    if(alpha > epsilon(0._rp)) then
      z1 = z0*n/(1.*n1)*l0
      if (z1 <= l0) then
        z = z1
      else
        z2 = (z0*n - n1)/(1.*(n-n1))
        z  = l0 + 0.5*(1.+tanh((z2-0.5)*alpha)/tanh(alpha/2.))*(1-l0)
      end if
    else
      z = z0
    end if
  end subroutine gridpoint_hybrid
  !-------------- END HYBRID --------------
  subroutine gridpoint_cluster_one_end_r(alpha,r0,r)
    !
    ! clustered at the upper side
    !
    implicit none
    real(rp), intent(in ) :: alpha,r0
    real(rp), intent(out) :: r
    if(alpha > epsilon(0._rp)) then
      r = 1._rp-1.0_rp*(1._rp+tanh((1._rp-r0-1.0_rp)*alpha)/tanh(alpha/1._rp))
      !r = 1._rp-1.0_rp*(1._rp+erf( (1._rp-r0-1.0_rp)*alpha)/erf( alpha/1._rp))
    else
      r = r0
    end if
  end subroutine gridpoint_cluster_one_end_r
  subroutine gridpoint_cluster_middle(alpha,z0,z)
    !
    ! clustered in the middle
    !
    implicit none
    real(rp), intent(in) :: alpha,z0
    real(rp), intent(out) :: z
    if(alpha > epsilon(0._rp)) then
      if(     z0 <= 0.5) then
        z = 0.5*(1.-1.+tanh(2.*alpha*(z0-0.))/tanh(alpha))
        !z = 0.5*(1.-1.+erf( 2.*alpha*(z0-0.))/erf( alpha))
      else if(z0  > 0.5) then
        z = 0.5*(1.+1.+tanh(2.*alpha*(z0-1.))/tanh(alpha))
        !z = 0.5*(1.+1.+erf( 2.*alpha*(z0-1.))/erf( alpha))
      end if
    else
      z = z0
    end if
  end subroutine gridpoint_cluster_middle
  subroutine gridpoint_natural(kg,nzg,z,kb_a,alpha_a,c_eta_a,dyp_a)
    !
    ! a physics-based, 'natural' grid stretching function for wall-bounded turbulence
    ! see Pirozzoli & Orlandi, JCP 439 - 110408 (2021)
    !
    ! clustered at the two sides
    !
    implicit none
    real(rp), parameter :: kb_p     = 32._rp,    &
                           alpha_p  = pi/1.5_rp, &
                           c_eta_p  = 0.8_rp,    &
                           dyp_p    = 0.05_rp
    integer , intent(in ) :: kg,nzg
    real(rp), intent(out) :: z
    real(rp), intent(in ), optional :: kb_a,alpha_a,c_eta_a,dyp_a
    real(rp)                        :: kb  ,alpha  ,c_eta  ,dyp
    real(rp) :: retau,n,k
    !
    ! handle input parameters
    !
    kb    = kb_p   ; if(present(kb_a   )) kb    = kb_a
    alpha = alpha_p; if(present(alpha_a)) alpha = alpha_a
    c_eta = c_eta_p; if(present(c_eta_a)) c_eta = c_eta_a
    dyp   = dyp_p  ; if(present(dyp_a  )) dyp   = dyp_a
    !
    ! determine retau
    !
    n = nzg/2._rp
    retau = 1._rp/(1._rp+(n/kb)**2)*(dyp*n+(3._rp/4._rp*alpha*c_eta*n)**(4._rp/3._rp)*(n/kb)**2)
#if 0
    if(kg==1) print*,'Grid targeting Retau = ',retau
#endif
    k = 1._rp*min(kg,(nzg-kg))
    !
    ! dermine z/(2h)
    !
    z = 1._rp/(1._rp+(k/kb)**2)*(dyp*k+(3._rp/4._rp*alpha*c_eta*k)**(4._rp/3._rp)*(k/kb)**2)/(2._rp*retau)
    if( kg > nzg-kg ) z = 1._rp-z
  end subroutine gridpoint_natural
end module mod_initgrid
