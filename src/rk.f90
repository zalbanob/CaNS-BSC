! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_rk
  use mod_mom  , only: momx_a,momy_a,momz_a, &
                       momx_d,momy_d,momz_d, &
                       momx_p,momy_p,momz_p, &
                       cmpt_wallshear, &
                       momx_d_xy,momy_d_xy,momz_d_xy, &
                       momx_d_z ,momy_d_z ,momz_d_z, &
                       mom_xyz_ad
  use mpi
  use mod_common_mpi, only: myid,ierr
  use mod_param, only: is_impdiff,is_impdiff_1d,is_boussinesq_buoyancy,is_fast_mom_kernels, &
                       at,datadir,restart,is_fringe,fringe_frac
  use mod_scal , only: scal,cmpt_scalflux,scalar
  use mod_utils, only: bulk_mean
  use mod_types
  implicit none
  public rk,rk_scal
  contains
  !---------------- TRIPPING ----------------
  subroutine read_trip_coeffs(fname,phi,phi_old,alpha,alpha_old)
    character(len=*), intent(in) :: fname
    real(rp), intent(inout)      :: phi(:),phi_old(:),alpha(:),alpha_old(:)
    integer :: iunit, ios, n, i
    open(newunit=iunit, file=fname, status="old", action="read", iostat=ios)
    if (ios /= 0) stop "Erreur ouverture trip.dat"
    read(iunit, *) n
    read(iunit, *) (phi(i),       i = 1, n)
    read(iunit, *) (phi_old(i),   i = 1, n)
    read(iunit, *) (alpha(i),     i = 1, n)
    read(iunit, *) (alpha_old(i), i = 1, n)
    close(iunit)
  end subroutine read_trip_coeffs

  subroutine write_trip_coeffs(fname,phi,phi_old,alpha,alpha_old)
    character(len=*), intent(in) :: fname
    real(rp), intent(in)          :: phi(:),phi_old(:),alpha(:),alpha_old(:)
    integer                       :: iunit, ios, n, i
    n = size(phi)
    open(newunit=iunit, file=fname, status="replace", action="write", iostat=ios)
    if (ios /= 0) stop "Erreur ouverture trip.dat"
    write(iunit,*) n
    write(iunit,*) (phi(i), i=1,n)
    write(iunit,*) (phi_old(i), i=1,n)
    write(iunit,*) (alpha(i), i=1,n)
    write(iunit,*) (alpha_old(i), i=1,n)
    close(iunit)
  end subroutine write_trip_coeffs
  !-------------- END TRIPPING --------------
  subroutine rk(rkpar,n,dli,dzci,dzfi,grid_vol_ratio_c,grid_vol_ratio_f,visc,dt,p, &
                is_forced,velf,bforce,gacc,beta,scalars,dudtrko,dvdtrko,dwdtrko,u,v,w,f, &
                istep,time,lo,ng,Ly,zc)
#if defined(_OPENACC)
    use mod_common_cudecomp, only: dudtrk_t => work, &
                                   dvdtrk_t => solver_buf_0, &
                                   dwdtrk_t => solver_buf_1
#endif
    !
    ! low-storage 3rd-order Runge-Kutta scheme
    ! for time integration of the momentum equations.
    !
    implicit none
    real(rp), intent(in   ), dimension(2)        :: rkpar
    integer , intent(in   ), dimension(3)        :: n
    real(rp), intent(in   )                      :: visc,dt
    real(rp), intent(in   ), dimension(3)        :: dli
    real(rp), intent(in   ), dimension(0:)       :: dzci,dzfi
    real(rp), intent(in   ), dimension(0:)       :: grid_vol_ratio_c,grid_vol_ratio_f
    real(rp), intent(in   ), dimension(0:,0:,0:) :: p
    logical , intent(in   ), dimension(3)        :: is_forced
    real(rp), intent(in   ), dimension(3)        :: velf,bforce
    real(rp), intent(in   ), dimension(3)        :: gacc
    real(rp), intent(in   )                      :: beta
    type(scalar), intent(in   ), target, dimension(:), optional :: scalars
    real(rp), intent(inout), dimension(1:,1:,1:) :: dudtrko,dvdtrko,dwdtrko
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out  ), dimension(3)        :: f
    !---------------- TRIPPING ----------------
    integer,  intent(in), optional :: istep
    real(rp), intent(in), optional :: time, Ly
    integer,  intent(in), optional, dimension(3) :: lo, ng
    real(rp), intent(in), optional, dimension(0:) :: zc
    !-------------- END TRIPPING --------------
    real(rp), pointer, contiguous, dimension(:,:,:) :: s
#if !defined(_OPENACC)
    real(rp), target       , allocatable, dimension(:,:,:), save :: dudtrk_t,dvdtrk_t,dwdtrk_t
#endif
    real(rp), pointer      , contiguous , dimension(:,:,:), save :: dudtrk  ,dvdtrk  ,dwdtrk
    real(rp),                allocatable, dimension(:,:,:), save :: dudtrkd ,dvdtrkd ,dwdtrkd
    logical, save :: is_first = .true.
    real(rp) :: factor1,factor2,factor12
    logical :: is_buoyancy
    integer  :: i,j,k
    !---------------- TRIPPING ----------------
    real(rp), parameter :: trip_delta0 = 1._RP
    real(rp), parameter :: trip_x0 = 5._RP
    real(rp), parameter :: trip_z0 = 0.175_rp*trip_delta0
    real(rp), parameter :: trip_lx = 1.4_rp * trip_delta0
    real(rp), parameter :: trip_lz = 0.35_rp * trip_delta0
    real(rp), parameter :: trip_pi = 3.141592_rp
    real(rp), parameter :: trip_ys = 0.6_rp*trip_delta0
    real(rp), parameter :: trip_ts = 1.4_rp*trip_delta0
    real(rp), allocatable, save :: trip_phi(:), trip_phi_old(:), trip_alpha(:), trip_alpha_old(:)
    real(rp), allocatable, save :: trip_hi(:), trip_hip1(:)
    real(rp), allocatable, save :: trip_exp_factor(:,:)
    integer,  save :: trip_iL = 0, trip_iR = -1, trip_zmin = 1, trip_zmax = 0
    integer,  save :: trip_icurrent = 0
    logical,  save :: trip_is_init = .false.
    logical,  save :: trip_restart_loaded = .false.
    !-------------- END TRIPPING --------------
    !---------------- FRINGE ----------------
    !real(rp), parameter :: fringe_xs_frac =
    real(rp), parameter :: fringe_beta = 2._RP
    real(rp), allocatable, save :: fringe_lambda(:)
    real(rp), allocatable, save :: fringe_u_tilde(:,:), fringe_u_tilde_tmp(:,:)
    integer,  save :: fringe_iL = 0, fringe_iR = -1, fringe_i_fringe = -1
    logical,  save :: fringe_is_init = .false.
    !-------------- END FRINGE --------------
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    !
    ! is_boussinesq_buoyancy = T will always imply present(scalars), as we do not allow for
    ! nscal = 0 with buoyancy on (see `param.f90`)
    !
    is_buoyancy = present(scalars).and.is_boussinesq_buoyancy
    if(is_buoyancy) then
      s => scalars(1)%val
    end if
    !
    ! initialization
    !
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      is_first = .false.
#if !defined(_OPENACC)
      allocate(dudtrk_t(n(1),n(2),n(3)),dvdtrk_t(n(1),n(2),n(3)),dwdtrk_t(n(1),n(2),n(3)))
#endif
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            dudtrko(i,j,k) = 0._rp
            dvdtrko(i,j,k) = 0._rp
            dwdtrko(i,j,k) = 0._rp
          end do
        end do
      end do
      if(is_impdiff) then
        allocate(dudtrkd(n(1),n(2),n(3)),dvdtrkd(n(1),n(2),n(3)),dwdtrkd(n(1),n(2),n(3)))
        !$acc enter data create(dudtrkd,dvdtrkd,dwdtrkd) async(1)
      end if
#if defined(_OPENACC)
      dudtrk(1:n(1),1:n(2),1:n(3)) => dudtrk_t(1:product(n(:)))
      dvdtrk(1:n(1),1:n(2),1:n(3)) => dvdtrk_t(1:product(n(:)))
      dwdtrk(1:n(1),1:n(2),1:n(3)) => dwdtrk_t(1:product(n(:)))
#else
      dudtrk  => dudtrk_t
      dvdtrk  => dvdtrk_t
      dwdtrk  => dwdtrk_t
#endif
    end if
    !
    if(is_fast_mom_kernels) then
      call mom_xyz_ad(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,u,v,w,dudtrk,dvdtrk,dwdtrk,dudtrkd,dvdtrkd,dwdtrkd)
    else
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            dudtrk(i,j,k) = 0._rp
            dvdtrk(i,j,k) = 0._rp
            dwdtrk(i,j,k) = 0._rp
          end do
        end do
      end do
      if(.not.is_impdiff) then
        call momx_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,u,dudtrk)
        call momy_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,v,dvdtrk)
        call momz_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,w,dwdtrk)
      else
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP parallel do   collapse(3) DEFAULT(shared)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              dudtrkd(i,j,k) = 0._rp
              dvdtrkd(i,j,k) = 0._rp
              dwdtrkd(i,j,k) = 0._rp
            end do
          end do
        end do
        if(.not.is_impdiff_1d) then
          call momx_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,u,dudtrkd)
          call momy_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,v,dvdtrkd)
          call momz_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,w,dwdtrkd)
        else
          call momx_d_xy(n(1),n(2),n(3),dli(1),dli(2),visc,u,dudtrk )
          call momy_d_xy(n(1),n(2),n(3),dli(1),dli(2),visc,v,dvdtrk )
          call momz_d_xy(n(1),n(2),n(3),dli(1),dli(2),visc,w,dwdtrk )
          call momx_d_z( n(1),n(2),n(3),dzci  ,dzfi  ,visc,u,dudtrkd)
          call momy_d_z( n(1),n(2),n(3),dzci  ,dzfi  ,visc,v,dvdtrkd)
          call momz_d_z( n(1),n(2),n(3),dzci  ,dzfi  ,visc,w,dwdtrkd)
        end if
        call momx_a(n(1),n(2),n(3),dli(1),dli(2),dzfi,u,v,w,dudtrk)
        call momy_a(n(1),n(2),n(3),dli(1),dli(2),dzfi,u,v,w,dvdtrk)
        call momz_a(n(1),n(2),n(3),dli(1),dli(2),dzci,u,v,w,dwdtrk)
      end if
    end if
    !
    !---------------- TRIPPING ----------------
    if(present(istep).and.present(time).and.present(lo).and.present(ng).and.present(Ly).and.present(zc)) then
      if(at /= 0._rp) then
        block
          integer :: Nf, iL_trip, iR_trip, imin, imax, itrip, rfn, zmin, zmax
          logical :: is_trip, need_realloc, trip_coeffs_changed
          real(rp) :: Af, gj, pt, b, xcoord, ycoord, zcoord
          Nf = int(Ly/trip_ys)
          if(.not.trip_is_init) then
            allocate(trip_phi(Nf), trip_phi_old(Nf), trip_alpha(Nf), trip_alpha_old(Nf))
            allocate(trip_hi(n(2)), trip_hip1(n(2)))
            trip_hi = 0._RP; trip_hip1 = 0._RP; trip_icurrent = 0
            trip_phi_old = 0._RP; trip_alpha_old = 0._RP
            trip_restart_loaded = .false.
            trip_is_init = .true.
          else if(size(trip_phi) /= Nf .or. size(trip_hi) /= n(2)) then
            deallocate(trip_phi, trip_phi_old, trip_alpha, trip_alpha_old, trip_hi, trip_hip1)
            trip_is_init = .false.
            allocate(trip_phi(Nf), trip_phi_old(Nf), trip_alpha(Nf), trip_alpha_old(Nf))
            allocate(trip_hi(n(2)), trip_hip1(n(2)))
            trip_hi = 0._RP; trip_hip1 = 0._RP; trip_icurrent = 0
            trip_phi_old = 0._RP; trip_alpha_old = 0._RP
            trip_restart_loaded = .false.
            trip_is_init = .true.
          end if
          imin    = max(1, int((trip_x0-2._rp*trip_lx)*dli(1)))
          imax    = min(ng(1), int((trip_x0+2._rp*trip_lx)*dli(1)))
          iL_trip = max(1,  imin - (lo(1)-1))
          iR_trip = min(n(1), imax - (lo(1)-1))
          zmin = -1
          zmax = -1
          do k = 1, n(3)
            if (zmin < 0 .and. zc(k) >= trip_z0 - trip_lz) zmin = k
            if (zc(k) <= trip_z0 + trip_lz)                zmax = k
          end do
          if (zmin < 0 .or. zmax < 0 .or. zmin > zmax) then
            zmin = 1
            zmax = 0
          end if
          is_trip = (iL_trip <= iR_trip) .and. (zmin <= zmax)
          Af      = at*(0.25_RP/trip_ts)
          itrip = int(time/trip_ts)
          pt    = time/trip_ts - real(itrip, rp)
          b     = 3.0_rp*pt**2 - 2.0_rp*pt**3
          trip_coeffs_changed = .false.
          if (restart .and. .not.trip_restart_loaded) then
            if (myid == 0) then
              call read_trip_coeffs(trim(datadir)//'trip_save.dat', trip_phi, trip_phi_old, trip_alpha, trip_alpha_old)
            end if
            call MPI_Bcast(trip_phi      , Nf, MPI_REAL_RP, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(trip_phi_old  , Nf, MPI_REAL_RP, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(trip_alpha    , Nf, MPI_REAL_RP, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(trip_alpha_old, Nf, MPI_REAL_RP, 0, MPI_COMM_WORLD, ierr)
            trip_hi = 0._RP
            trip_hip1 = 0._RP
            do rfn = 1, Nf
              do j = 1, n(2)
                ycoord = real(j+lo(2)-1, rp) / dli(2)
                trip_hip1(j) = trip_hip1(j) + trip_alpha     (rfn)*cos(2._RP*trip_pi*real(rfn,rp)*ycoord/Ly + trip_phi     (rfn))
                trip_hi  (j) = trip_hi  (j) + trip_alpha_old(rfn)*cos(2._RP*trip_pi*real(rfn,rp)*ycoord/Ly + trip_phi_old (rfn))
              end do
            end do
            trip_hip1 = trip_hip1/sqrt(real(Nf,rp))
            trip_hi   = trip_hi  /sqrt(real(Nf,rp))
            trip_icurrent = itrip
            trip_restart_loaded = .true.
            trip_coeffs_changed = .true.
          end if
          if (itrip > trip_icurrent) then
            trip_icurrent = itrip
            trip_hi  = trip_hip1
            trip_hip1 = 0._RP
            if (myid == 0) then
              do rfn = 1, Nf
                call random_number(trip_phi(rfn));   trip_phi(rfn)   = 2._RP*trip_pi*trip_phi(rfn)
                call random_number(trip_alpha(rfn))
              end do
            end if
            call MPI_Bcast(trip_phi  , Nf, MPI_REAL_RP, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(trip_alpha, Nf, MPI_REAL_RP, 0, MPI_COMM_WORLD, ierr)
            do rfn = 1, Nf
              do j = 1, n(2)
                ycoord = real(j+lo(2)-1, rp) / dli(2)
                trip_hip1(j) = trip_hip1(j) + trip_alpha(rfn)*cos(2._RP*trip_pi*real(rfn,rp)*ycoord/Ly + trip_phi(rfn))
              end do
            end do
            if (myid == 0) then
              call write_trip_coeffs(trim(datadir)//'trip.dat', trip_phi, trip_phi_old, trip_alpha, trip_alpha_old)
            end if
            trip_phi_old = trip_phi
            trip_alpha_old = trip_alpha
            trip_hip1 = trip_hip1 / sqrt(real(Nf, rp))
            trip_coeffs_changed = .true.
          end if
          if(is_trip) then
            need_realloc = .not.allocated(trip_exp_factor)
            need_realloc = need_realloc .or. (trip_iL /= iL_trip) .or. (trip_iR /= iR_trip) .or. (trip_zmin /= zmin) .or. (trip_zmax /= zmax)
            if(need_realloc) then
              if(allocated(trip_exp_factor)) deallocate(trip_exp_factor)
              allocate(trip_exp_factor(iL_trip:iR_trip, zmin:zmax))
              do k = zmin, zmax
                zcoord = zc(k)
                do i = iL_trip, iR_trip
                  xcoord = real(i+lo(1)-1, rp) / dli(1)
                  trip_exp_factor(i,k) = exp( -((xcoord-trip_x0)/trip_lx)**2 - ((zcoord-trip_z0)/trip_lz)**2 )
                end do
              end do
              !$acc enter data copyin(trip_exp_factor(iL_trip:iR_trip, zmin:zmax))
              trip_iL = iL_trip; trip_iR = iR_trip; trip_zmin = zmin; trip_zmax = zmax
              !$acc enter data copyin(trip_hi, trip_hip1)
            end if
            if (trip_coeffs_changed) then
              !$acc update device(trip_hi, trip_hip1)
            end if
            !$acc parallel loop gang collapse(2) default(present) private(gj) async(1)
            do k = zmin, zmax
              do j = 1, n(2)
                gj = Af*((1._rp - b)*trip_hi(j) + b*trip_hip1(j))
                !$acc loop vector
                do i = iL_trip, iR_trip
                  dwdtrk(i,j,k) = dwdtrk(i,j,k) + gj * trip_exp_factor(i,k)
                end do
              end do
            end do
          end if
        end block
      end if
    end if
    !-------------- END TRIPPING --------------

    !---------------- FRINGE ----------------
    if(present(istep).and.present(lo).and.present(ng)) then
      if(is_fringe) then
        block
          integer :: i_fringe, iL_fringe, iR_fringe
          logical :: is_fringe_local
          real(rp) :: xi, s_
          i_fringe  = max(1, int(fringe_frac*ng(1)))
          iL_fringe = max(1, i_fringe - (lo(1)-1))
          iR_fringe = n(1)
          is_fringe_local = (lo(1)+n(1)-1 >= i_fringe) .and. (iL_fringe <= iR_fringe)
          if (is_fringe_local) then
            if(.not.fringe_is_init .or. fringe_iL /= iL_fringe .or. fringe_iR /= iR_fringe .or. fringe_i_fringe /= i_fringe) then
              if(allocated(fringe_u_tilde)) deallocate(fringe_u_tilde, fringe_u_tilde_tmp, fringe_lambda)
              allocate(fringe_u_tilde(iL_fringe:iR_fringe, n(3)))
              allocate(fringe_u_tilde_tmp(iL_fringe:iR_fringe, n(3)))
              allocate(fringe_lambda(iL_fringe:iR_fringe))
              fringe_u_tilde = 0._rp
              fringe_u_tilde_tmp = 0._rp
              !$acc enter data create(fringe_u_tilde, fringe_u_tilde_tmp)
              !$acc parallel loop gang collapse(2) default(present) present(u,fringe_u_tilde) private(j)
              do k = 1, n(3)
                do i = iL_fringe, iR_fringe
                  s_ = 0._rp
                  !$acc loop reduction(+:s_)
                  do j = 1, n(2)
                    s_ = s_ + u(i,j,k)
                  end do
                  fringe_u_tilde(i,k) = s_ / real(n(2), rp)
                end do
              end do
              !$acc end parallel loop
              do i = iL_fringe, iR_fringe
                xi = real((i+lo(1)-1) - i_fringe, rp) / real(ng(1) - i_fringe, rp) - 1._RP
                fringe_lambda(i) = 0.5_rp * (1._rp + tanh(fringe_beta*xi))
              end do
              !$acc enter data copyin(fringe_lambda)
              fringe_iL = iL_fringe; fringe_iR = iR_fringe; fringe_i_fringe = i_fringe
              fringe_is_init = .true.
            end if
            if (mod(istep,100) == 0) then
              !$acc parallel loop gang collapse(2) default(present) present(u,fringe_u_tilde_tmp) private(j)
              do k = 1, n(3)
                do i = iL_fringe, iR_fringe
                  s_ = 0._rp
                  !$acc loop reduction(+:s_)
                  do j = 1, n(2)
                    s_ = s_ + u(i,j,k)
                  end do
                  fringe_u_tilde_tmp(i,k) = s_ / real(n(2), rp)
                end do
              end do
              !$acc end parallel loop
              !$acc parallel loop gang collapse(2) default(present) present(fringe_u_tilde,fringe_u_tilde_tmp)
              do k = 1, n(3)
                do i = iL_fringe, iR_fringe
                  fringe_u_tilde(i,k) = 0.9_rp*fringe_u_tilde(i,k) + 0.1_rp*fringe_u_tilde_tmp(i,k)
                end do
              end do
              !$acc end parallel loop
            end if
            !$acc parallel loop collapse(3) default(present) async(1)
            !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(shared)
            do k = 1, n(3)
              do j = 1, n(2)
                do i = iL_fringe, iR_fringe
                  dudtrk(i,j,k) = dudtrk(i,j,k) + fringe_lambda(i) * (fringe_u_tilde(i,k) - u(i,j,k))
                end do
              end do
            end do
          end if
        end block
      end if
    end if
    !-------------- END FRINGE --------------

#if !defined(_LOOP_UNSWITCHING)
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k) + &
                                factor12*(bforce(1) - dli(1)*( p(i+1,j,k)-p(i,j,k)))
          if(is_buoyancy) then
            u(i,j,k) = u(i,j,k) - factor12*gacc(1)*beta*0.5*(s(i+1,j,k)+s(i,j,k))
          end if
          !
          v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k) + &
                                factor12*(bforce(2) - dli(2)*( p(i,j+1,k)-p(i,j,k)))
          if(is_buoyancy) then
            v(i,j,k) = v(i,j,k) - factor12*gacc(2)*beta*0.5*(s(i,j+1,k)+s(i,j,k))
          end if
          !
          w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k) + &
                                factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k)))
          if(is_buoyancy) then
            w(i,j,k) = w(i,j,k) - factor12*gacc(3)*beta*0.5*(s(i,j,k+1)+s(i,j,k))
          end if
          !
          if(is_impdiff) then
            u(i,j,k) = u(i,j,k) + factor12*dudtrkd(i,j,k)
            v(i,j,k) = v(i,j,k) + factor12*dvdtrkd(i,j,k)
            w(i,j,k) = w(i,j,k) + factor12*dwdtrkd(i,j,k)
          end if
        end do
      end do
    end do
#else
    if(.not.is_impdiff .and. .not.is_buoyancy) then
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k) + &
                                  factor12*(bforce(1) - dli(1)*(p(i+1,j,k)-p(i,j,k)))
            v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k) + &
                                  factor12*(bforce(2) - dli(2)*(p(i,j+1,k)-p(i,j,k)))
            w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k) + &
                                  factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k)))
          end do
        end do
      end do
    else if(is_impdiff .and. .not.is_buoyancy) then
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k) + &
                                  factor12*(bforce(1) - dli(1)*(p(i+1,j,k)-p(i,j,k)) + &
                                            dudtrkd(i,j,k))
            v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k) + &
                                  factor12*(bforce(2) - dli(2)*(p(i,j+1,k)-p(i,j,k)) + &
                                            dvdtrkd(i,j,k))
            w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k) + &
                                  factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k)) + &
                                            dwdtrkd(i,j,k))
          end do
        end do
      end do
    else if(.not.is_impdiff .and. is_buoyancy) then
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k) + &
                                  factor12*(bforce(1) - dli(1)*(p(i+1,j,k)-p(i,j,k))) &
                                          - gacc(1)*beta*0.5*(s(i+1,j,k)+s(i,j,k))
            v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k) + &
                                  factor12*(bforce(2) - dli(2)*(p(i,j+1,k)-p(i,j,k))) &
                                          - gacc(2)*beta*0.5*(s(i,j+1,k)+s(i,j,k))
            w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k) + &
                                  factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k))) &
                                          - gacc(3)*beta*0.5*(s(i,j,k+1)+s(i,j,k))
          end do
        end do
      end do
    else
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k) + &
                                  factor12*(bforce(1) - dli(1)*(p(i+1,j,k)-p(i,j,k)) &
                                          - gacc(1)*beta*0.5*(s(i+1,j,k)+s(i,j,k)) + &
                                            dudtrkd(i,j,k))
            v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k) + &
                                  factor12*(bforce(2) - dli(2)*(p(i,j+1,k)-p(i,j,k)) &
                                          - gacc(2)*beta*0.5*(s(i,j+1,k)+s(i,j,k)) + &
                                            dvdtrkd(i,j,k))
            w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k) + &
                                  factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k)) &
                                          - gacc(3)*beta*0.5*(s(i,j,k+1)+s(i,j,k)) + &
                                            dwdtrkd(i,j,k))
          end do
        end do
      end do
    end if
#endif
    !
    ! replaced previous pointer swap to save memory on GPUs by using already allocated
    ! buffers
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP parallel do   collapse(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dudtrko(i,j,k) = dudtrk(i,j,k)
          dvdtrko(i,j,k) = dvdtrk(i,j,k)
          dwdtrko(i,j,k) = dwdtrk(i,j,k)
        end do
      end do
    end do
!#if 0 /*pressure gradient term treated explicitly above */
!    !$acc parallel loop collapse(3) default(present) async(1)
!    !$OMP parallel do   collapse(3) DEFAULT(shared)
!    do k=1,n(3)
!      do j=1,n(2)
!        do i=1,n(1)
!          dudtrk(i,j,k) = 0._rp
!          dvdtrk(i,j,k) = 0._rp
!          dwdtrk(i,j,k) = 0._rp
!        end do
!      end do
!    end do
!    call momx_p(n(1),n(2),n(3),dli(1),bforce(1),p,dudtrk)
!    call momy_p(n(1),n(2),n(3),dli(2),bforce(2),p,dvdtrk)
!    call momz_p(n(1),n(2),n(3),dzci  ,bforce(3),p,dwdtrk)
!    !$acc parallel loop collapse(3)
!    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
!    do k=1,n(3)
!      do j=1,n(2)
!        do i=1,n(1)
!          u(i,j,k) = u(i,j,k) + factor12*dudtrk(i,j,k)
!          v(i,j,k) = v(i,j,k) + factor12*dvdtrk(i,j,k)
!          w(i,j,k) = w(i,j,k) + factor12*dwdtrk(i,j,k)
!        end do
!      end do
!    end do
!#endif
    !
    ! compute bulk velocity forcing
    !
    call cmpt_bulk_forcing(n,is_forced,velf,grid_vol_ratio_c,grid_vol_ratio_f,u,v,w,f)
    !
    if(is_impdiff) then
      !
      ! compute rhs of Helmholtz equation
      !
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            u(i,j,k) = u(i,j,k) - .5_rp*factor12*dudtrkd(i,j,k)
            v(i,j,k) = v(i,j,k) - .5_rp*factor12*dvdtrkd(i,j,k)
            w(i,j,k) = w(i,j,k) - .5_rp*factor12*dwdtrkd(i,j,k)
          end do
        end do
      end do
    end if
  end subroutine rk
  !
  subroutine rk_scal(rkpar,n,dli,l,dzci,dzfi,grid_vol_ratio_f,alpha,dt,is_bound,u,v,w, &
                     is_forced,scalf,ssource,fluxo,dsdtrko,s,f)
#if defined(_OPENACC)
    use mod_common_cudecomp, only: dsdtrk_t => work
#endif
    !
    ! low-storage 3rd-order Runge-Kutta scheme
    ! for time integration of the scalar field.
    !
    implicit none
    logical , parameter :: is_cmpt_wallflux = .false.
    real(rp), intent(in   ), dimension(2) :: rkpar
    integer , intent(in   ), dimension(3) :: n
    real(rp), intent(in   ), dimension(3) :: dli,l
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in   ), dimension(:) :: grid_vol_ratio_f
    real(rp), intent(in   ) :: alpha,dt
    logical , intent(in   ), dimension(0:1,3)    :: is_bound
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w
    logical , intent(in   ) :: is_forced
    real(rp), intent(in   ) :: scalf,ssource
    real(rp), intent(inout), dimension(0:1,3) :: fluxo
    real(rp), intent(inout), dimension(1:,1:,1:) :: dsdtrko
    real(rp), intent(inout), dimension(0:,0:,0:) :: s
    real(rp), intent(out  ) :: f
    !
#if !defined(_OPENACC)
    real(rp), target       , allocatable, dimension(:,:,:), save :: dsdtrk_t
#endif
    real(rp), pointer      , contiguous , dimension(:,:,:), save :: dsdtrk
    real(rp), target       , allocatable, dimension(:,:,:), save :: dsdtrkd
    logical, save :: is_first = .true.
    !
    real(rp) :: factor1,factor2,factor12
    real(rp), dimension(0:1,3) :: flux
    integer :: i,j,k
    real(rp) :: mean
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      is_first = .false.
#if !defined(_OPENACC)
      allocate(dsdtrk_t(1:n(1),1:n(2),1:n(3)))
#endif
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            dsdtrko(i,j,k) = 0._rp
          end do
        end do
      end do
      if(is_impdiff) then
        allocate(dsdtrkd(n(1),n(2),n(3)))
        !$acc enter data create(dsdtrkd) async(1)
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP parallel do   collapse(3) DEFAULT(shared)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              dsdtrkd(i,j,k) = 0._rp
            end do
          end do
        end do
      end if
    end if
#if defined(_OPENACC)
    dsdtrk(1:n(1),1:n(2),1:n(3)) => dsdtrk_t(1:product(n(:)))
#else
    dsdtrk => dsdtrk_t
#endif
    !
    call scal(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,alpha,u,v,w,s,dsdtrk,dsdtrkd)
#if !defined(_LOOP_UNSWITCHING)
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          s(i,j,k) = s(i,j,k) + factor1*dsdtrk(i,j,k) + factor2*dsdtrko(i,j,k) + factor12*ssource
          if(is_impdiff) then
            s(i,j,k) = s(i,j,k) + factor12*dsdtrkd(i,j,k)
          end if
        end do
      end do
    end do
#else
    if(.not.is_impdiff) then
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            s(i,j,k) = s(i,j,k) + factor1*dsdtrk(i,j,k) + factor2*dsdtrko(i,j,k) + factor12*ssource
          end do
        end do
      end do
    else
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            s(i,j,k) = s(i,j,k) + factor1*dsdtrk(i,j,k) + factor2*dsdtrko(i,j,k) + &
                                  factor12*(ssource + dsdtrkd(i,j,k))
          end do
        end do
      end do
    end if
#endif
    !
    ! compute wall scalar flux
    !
    if(is_cmpt_wallflux) then
      call cmpt_scalflux(n,is_bound,l,dli,dzci,dzfi,alpha,s(:,:,:),flux)
      f = (factor1*sum((flux( 0,:)+flux( 1,:))/l(:)) + &
           factor2*sum((fluxo(0,:)+fluxo(1,:))/l(:)))
      fluxo(:,:) = flux(:,:)
    end if
    !
    ! bulk scalar forcing
    !
    if(is_forced) then
      call bulk_mean(n,grid_vol_ratio_f,s(:,:,:),mean)
      f = scalf - mean
    end if
    if(is_impdiff) then
      !
      ! compute rhs of Helmholtz equation
      !
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            s(i,j,k) = s(i,j,k) - .5_rp*factor12*dsdtrkd(i,j,k)
          end do
        end do
      end do
    end if
    !
    ! replaced previous pointer swap to save memory on GPUs by using already allocated
    ! buffers
    !
    !$acc parallel loop collapse(3) default(present) async(1) ! not really necessary
    !$OMP parallel do   collapse(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dsdtrko(i,j,k) = dsdtrk(i,j,k)
        end do
      end do
    end do
  end subroutine rk_scal
  !
  subroutine cmpt_bulk_forcing(n,is_forced,velf,grid_vol_ratio_c,grid_vol_ratio_f,u,v,w,f)
    implicit none
    integer , intent(in   ), dimension(3) :: n
    logical , intent(in   ), dimension(3) :: is_forced
    real(rp), intent(in   ), dimension(3) :: velf
    real(rp), intent(in   ), dimension(0:) :: grid_vol_ratio_c,grid_vol_ratio_f
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out  ), dimension(3) :: f
    real(rp) :: mean
    !
    ! bulk velocity forcing
    !
    f(:) = 0.
    if(is_forced(1)) then
      call bulk_mean(n,grid_vol_ratio_f,u,mean)
      f(1) = velf(1) - mean
    end if
    if(is_forced(2)) then
      call bulk_mean(n,grid_vol_ratio_f,v,mean)
      f(2) = velf(2) - mean
    end if
    if(is_forced(3)) then
      call bulk_mean(n,grid_vol_ratio_c,w,mean)
      f(3) = velf(3) - mean
    end if
  end subroutine cmpt_bulk_forcing
  !
  subroutine cmpt_bulk_forcing_alternative(rkpar,n,dli,l,dzci,dzfi,visc,dt,is_bound,is_forced,u,v,w,tauxo,tauyo,tauzo,f,is_first)
    !
    ! computes the pressure gradient to be added to the flow that perfectly balances the wall shear stresses
    ! this effectively prescribes zero net acceleration, which allows to sustain a constant mass flux
    !
    implicit none
    real(rp), intent(in), dimension(2) :: rkpar
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in) :: visc,dt
    real(rp), intent(in   ), dimension(3) :: dli,l
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    logical , intent(in   ), dimension(0:1,3)    :: is_bound
    logical , intent(in   ), dimension(3) :: is_forced
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(inout), dimension(0:1,3) :: tauxo,tauyo,tauzo
    real(rp), intent(inout), dimension(3) :: f
    real(rp), dimension(3) :: f_aux
    logical , intent(in   ) :: is_first
    real(rp), dimension(0:1,3) :: taux,tauy,tauz
    real(rp), dimension(3) :: taux_tot,tauy_tot,tauz_tot,tauxo_tot,tauyo_tot,tauzo_tot
    real(rp) :: factor1,factor2,factor12
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = (factor1 + factor2)/2.
    !
    call cmpt_wallshear(n,is_forced,is_bound,l,dli,dzci,dzfi,visc,u,v,w,taux,tauy,tauz)
    taux_tot(:) = sum(taux(0:1,:),1); tauxo_tot(:) = sum(tauxo(0:1,:),1)
    tauy_tot(:) = sum(tauy(0:1,:),1); tauyo_tot(:) = sum(tauyo(0:1,:),1)
    tauz_tot(:) = sum(tauz(0:1,:),1); tauzo_tot(:) = sum(tauzo(0:1,:),1)
    if(.not.is_impdiff) then
      if(is_first) then
        f(1) = (factor1*sum(taux_tot(:)/l(:)) + factor2*sum(tauxo_tot(:)/l(:)))
        f(2) = (factor1*sum(tauy_tot(:)/l(:)) + factor2*sum(tauyo_tot(:)/l(:)))
        f(3) = (factor1*sum(tauz_tot(:)/l(:)) + factor2*sum(tauzo_tot(:)/l(:)))
        tauxo(:,:) = taux(:,:)
        tauyo(:,:) = tauy(:,:)
        tauzo(:,:) = tauz(:,:)
      end if
    else
      if(is_impdiff_1d) then
        f_aux(1) = factor12*taux_tot(3)/l(3)
        f_aux(2) = factor12*tauy_tot(3)/l(3)
        if(is_first) then
          f(1) = factor1*taux_tot(2)/l(2) + factor2*tauxo_tot(2)/l(2) + f_aux(1)
          f(2) = factor1*tauy_tot(1)/l(1) + factor2*tauyo_tot(1)/l(1) + f_aux(2)
          f(3) = factor1*sum(tauz_tot(1:2)/l(1:2)) + factor2*sum(tauzo_tot(1:2)/l(1:2))
          tauxo(:,1:2) = taux(:,1:2)
          tauyo(:,1:2) = tauy(:,1:2)
          tauzo(:,1:2) = tauz(:,1:2)
        else
          f(1) = f(1) + f_aux(1)
          f(2) = f(2) + f_aux(2)
        end if
      else
        f_aux(:) = factor12*[sum(taux_tot(:)/l(:)), &
                             sum(tauy_tot(:)/l(:)), &
                             sum(tauz_tot(:)/l(:))]
        if(is_first) then
           f(:) = f_aux(:)
        else
           f(:) = f(:) + f_aux(:)
        end if
      end if
    end if
  end subroutine cmpt_bulk_forcing_alternative
end module mod_rk
