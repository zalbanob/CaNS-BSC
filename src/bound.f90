 
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_bound
  use iso_fortran_env, only: real32
  use mpi
  use mod_common_mpi, only: ierr,halo,myid
  use mod_param     , only: ipencil_axis
  use mod_types
  implicit none
  private
  public boundp,bounduvw,updt_rhs_b
  !---------------- INLET_PLAN ----------------
  public inflow_register,compute_means
  real(real32), pointer :: u2d_in(:,:), v2d_in(:,:), w2d_in(:,:), p2d_in(:,:) => null()
  !$acc declare create(u2d_in, v2d_in, w2d_in, p2d_in)
  !-------------- END INLET_PLAN --------------
  contains
  !---------------- INLET_PLAN ----------------
  subroutine inflow_register(u2d, v2d, w2d, p2d)
    real(real32), target, intent(in) :: u2d(:,:), v2d(:,:), w2d(:,:), p2d(:,:)
    u2d_in => u2d;  v2d_in => v2d;  w2d_in => w2d;  p2d_in => p2d
    !$acc enter data attach(u2d_in, v2d_in, w2d_in, p2d_in)
  end subroutine inflow_register

  subroutine compute_means(u2d, v2d, w2d)
    real(real32), intent(in)  :: u2d(:,:,:), v2d(:,:,:), w2d(:,:,:)
    integer :: s, ny, nz, ns
    ny = size(u2d,1); nz = size(u2d,2); ns = size(u2d,3)
  end subroutine compute_means

  subroutine inlet_replay_apply_uvw(istep, u, v, w, impose_norm_bc)
    integer , intent(in)    :: istep
    real(rp), intent(inout) :: u(:,:,:), v(:,:,:), w(:,:,:)
    logical , intent(in)    :: impose_norm_bc
    integer :: i0, ig, j0, j1, k0, k1, ny2d, nz2d
    integer :: j, k, jj, kk
    real(rp) :: gamma
    gamma = min(1._RP, real(istep,rp)/200._RP)
    ig = lbound(u,1)
    i0 = ig + 1
    j0 = lbound(u,2) + 1
    j1 = ubound(u,2) - 1
    k0 = lbound(u,3) + 1
    k1 = ubound(u,3) - 1
    ny2d = size(u2d_in,1)
    nz2d = size(u2d_in,2)
    if (impose_norm_bc) then
      !$acc parallel loop collapse(2) present(u,u2d_in)
      do k = k0, k1
        do j = j0, j1
          kk = k - k0 + 1
          jj = j - j0 + 1
          !u(ig, j, k) = ubar(jj, kk) + gamma*(u2d_in(jj, kk)-ubar(jj, kk))
          u(ig,j,k)=u2d_in(jj,kk)
        end do
      end do
    end if
      !$acc parallel loop collapse(2) present(v,v2d_in)
      do k = k0, k1
        do j = j0, j1
          kk = k - k0 + 1
          jj = j - j0 + 1
       !   !v(ig, j, k) = -v(i0, j, k) + 2.0_rp * ( vbar(jj, kk) + gamma*(v2d_in(jj, kk)-vbar(jj, kk)) )
          v(ig,j,k)=-v(i0,j,k) + 2.0_rp*v2d_in(jj,kk)
        end do
      end do
     !$acc parallel loop collapse(2) present(w,w2d_in)
      do k = k0, k1
        do j = j0, j1
          kk = k - k0 + 1
          jj = j - j0 + 1
          w(ig, j, k) = -w(i0,j,k) + 2.0_rp*w2d_in(jj,kk)
        end do
      end do
    !$acc wait
  end subroutine inlet_replay_apply_uvw

  !-------------- END INLET_PLAN --------------

  !---------------- BLASIUS ----------------
  subroutine blasius_profile_boundary(n,dl,visc,p,lo,zc)
    implicit none
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    real(rp), intent(in), dimension(0:) :: zc
    real(rp), intent(in) :: visc
    real(rp), intent(in) :: dl
    integer, intent(in) :: n,lo
    integer :: k
    real(rp) :: z, eta, cal_temp
    real(rp), parameter :: x0 = 43.3_RP
    !$acc parallel loop default(present)
    do k=1,n
      z = zc(k)
      eta = z*sqrt(1._RP/(visc*x0))
      cal_temp = 1._RP   &
                -1.4564_RP*exp(-1._RP*eta)  &
                + 1.2956_RP*exp(-1._RP*eta)*(1._RP-eta) &
                -0.8392_RP*exp(-2._RP*eta)
      p(0,:,k)= cal_temp
    end do
  end subroutine blasius_profile_boundary
  !-------------- END BLASIUS --------------

  !---------------- OUTFLOW_ADV ----------------
  subroutine advective_outflow(n,ng,lo,is_bound,istep,rkpar,dt,dl,u)
    implicit none
    integer,  intent(in), dimension(3)        :: n, ng, lo
    logical , intent(in), dimension(0:1,3  )  :: is_bound
    integer,  intent(in)                      :: istep
    real(rp), intent(in),    dimension(2)     :: rkpar
    real(rp), intent(in)                      :: dt
    real(rp), intent(in), dimension(3)        :: dl
    real(rp), intent(inout), dimension(0:,0:,0:) :: u
    integer :: i,j,k,i0_g,iL,iR,kg,ng3,cnt_i
    real(rp) :: factor1,factor2,invdl1
    real(rp), allocatable, save :: u_adv(:)
    real(rp), allocatable, save :: sum_loc(:)
    integer , allocatable, save :: cnt_loc(:), cnt_glob(:)
    real(rp), allocatable, save :: dadvtrk(:,:), dadvtrko(:,:)
    logical,  save :: is_first = .true.
    logical,  save :: u_adv_computed = .false.
    integer :: ii, jj
    logical :: first_u_adv_compute
    real(rp) :: s, meank, tmp
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    invdl1  = 1._rp/dl(1)
    ng3     = ng(3)
    i0_g = max(1, int(0.95_rp*ng(1)))
    iL   = max(1, i0_g - (lo(1)-1))
    iR   = min(n(1), ng(1) - (lo(1)-1))
    if (.not. allocated(u_adv))     allocate(u_adv(ng3))
    if (.not. allocated(sum_loc))   allocate(sum_loc(ng3))
    if (.not. allocated(cnt_loc))   allocate(cnt_loc(ng3))
    if (.not. allocated(cnt_glob))  allocate(cnt_glob(ng3))
    if (.not. allocated(dadvtrk))   allocate(dadvtrk(n(2),n(3)))
    if (.not. allocated(dadvtrko))  allocate(dadvtrko(n(2),n(3)))
    if (is_first) then
      dadvtrk   = 0._rp
      dadvtrko  = 0._rp
      u_adv     = 0._rp
      sum_loc   = 0._rp
      cnt_loc   = 0
      cnt_glob  = 0
      !$acc enter data copyin(u_adv, sum_loc, dadvtrk, dadvtrko)
      is_first = .false.
    end if
    !if (.not. u_adv_computed .or. mod(istep,600) == 0) then
    if (mod(istep,600) == 0) then
      if (lo(3) < 1 .or. lo(3) > ng3) return
      first_u_adv_compute = .not. u_adv_computed
      u_adv_computed = .true.
      cnt_i = max(0, iR - iL + 1)
      cnt_loc(:) = 0
      do k = 1, n(3)
        kg = lo(3) - 1 + k
        cnt_loc(kg) = cnt_i*n(2)
      end do
      call MPI_Allreduce(cnt_loc, cnt_glob, ng3, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      !$acc parallel loop present(sum_loc) default(present)
      do k = 1, ng3
        sum_loc(k) = 0._rp
      end do
      !$acc end parallel loop
      !$acc parallel loop gang present(u, sum_loc) default(present)
      do k = 1, n(3)
        s = 0._rp
        !$acc loop vector collapse(2) reduction(+:s)
        do jj = 1, n(2)
          do ii = iL, iR
            s = s + u(ii,jj,k)
          end do
        end do
        kg = lo(3) - 1 + k
        sum_loc(kg) = s
      end do
      !$acc end parallel loop
      !$acc update self(sum_loc)
      call MPI_Allreduce(MPI_IN_PLACE, sum_loc, ng3, MPI_REAL_RP, MPI_SUM, MPI_COMM_WORLD, ierr)
      do kg = 1, ng3
        if (cnt_glob(kg) > 0) then
          meank = sum_loc(kg)/real(cnt_glob(kg), rp)
          if (first_u_adv_compute) then
            tmp = meank
          else
            tmp = 0.9_rp*u_adv(kg) + 0.1_rp*meank
          end if
          u_adv(kg) = max(tmp, 0.30_rp)
        else
          u_adv(kg) = 0.30_rp
        end if
      end do
      !$acc update device(u_adv)
    end if
    if (is_bound(1,1)) then
      !$acc parallel loop collapse(2) present(u,u_adv,dadvtrk) default(present) private(kg)
      do k=1,n(3)
        do j=1,n(2)
          kg = lo(3) - 1 + k
          dadvtrk(j,k) = - 0.2_RP*u_adv(kg) * ( u(n(1),j,k) - u(n(1)-1,j,k) ) * invdl1
        end do
      end do
      !$acc end parallel loop
      !$acc parallel loop collapse(2) present(u,dadvtrk,dadvtrko) default(present)
      do k=1,n(3)
        do j=1,n(2)
          u(n(1),  j,k) = u(n(1),  j,k) + factor1*dadvtrk(j,k) + factor2*dadvtrko(j,k)
          u(n(1)+1,j,k) = u(n(1),  j,k)
        end do
      end do
      !$acc end parallel loop
      !$acc parallel loop collapse(2) present(dadvtrk,dadvtrko) default(present)
      do k=1,n(3)
        do j=1,n(2)
          dadvtrko(j,k) = dadvtrk(j,k)
        end do
      end do
      !$acc end parallel loop
    end if
  end subroutine advective_outflow
  !-------------- END OUTFLOW_ADV --------------
  subroutine bounduvw(cbc,n,bc,nb,is_bound,is_correc,dl,dzc,dzf,u,v,w, &
                      istep,rkpar,dt,visc,lo,ng,zc)
    !
    ! imposes velocity boundary conditions
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:1,3,3) :: bc
    integer , intent(in), dimension(0:1,3  ) :: nb
    logical , intent(in), dimension(0:1,3  ) :: is_bound
    logical , intent(in)                     :: is_correc
    real(rp), intent(in), dimension(3 ) :: dl
    real(rp), intent(in), dimension(0:) :: dzc,dzf
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    !---------------- INLET_PLAN ----------------
    integer , intent(in), optional :: istep
    real(rp), intent(in), optional, dimension(2) :: rkpar
    real(rp), intent(in), optional :: dt, visc
    integer , intent(in), optional, dimension(3) :: lo, ng
    real(rp), intent(in), optional, dimension(0:) :: zc
    !-------------- END INLET_PLAN --------------
    logical :: impose_norm_bc
    integer :: idir,nh
    !
    nh = 1
    !
#if !defined(_OPENACC)
    do idir = 1,3
      call updthalo(nh,halo(idir),nb(:,idir),idir,u)
      call updthalo(nh,halo(idir),nb(:,idir),idir,v)
      call updthalo(nh,halo(idir),nb(:,idir),idir,w)
    end do
#else
    call updthalo_gpu(nh,cbc(0,:,1)//cbc(1,:,1)==['PP','PP','PP'],u)
    call updthalo_gpu(nh,cbc(0,:,2)//cbc(1,:,2)==['PP','PP','PP'],v)
    call updthalo_gpu(nh,cbc(0,:,3)//cbc(1,:,3)==['PP','PP','PP'],w)
#endif
    !
    impose_norm_bc = (.not.is_correc).or.(cbc(0,1,1)//cbc(1,1,1) == 'PP')
    if(is_bound(0,1)) then
      if (cbc(0,1,1) == 'B') then
        if (impose_norm_bc) then
          if(.not.(present(zc).and.present(visc))) then
            if(myid == 0) print*, 'ERROR: BLASIUS BC requires present(zc,visc) in bounduvw().'
            error stop
          end if
          call blasius_profile_boundary(n(3),dl(3),visc,u,0,zc)
        end if
        call set_bc(cbc(0,1,2),0,1,nh,.true. ,bc(0,1,2),dl(1),v)
        call set_bc(cbc(0,1,3),0,1,nh,.true. ,bc(0,1,3),dl(1),w)
      else if (cbc(0,1,1) == 'I') then
        if(.not.present(istep)) then
          if(myid == 0) print*, 'ERROR: INLET_REPLAY BC requires present(istep) in bounduvw().'
          error stop
        end if
        call inlet_replay_apply_uvw(istep, u, v, w, impose_norm_bc)
      else
        if(impose_norm_bc) call set_bc(cbc(0,1,1),0,1,nh,.false.,bc(0,1,1),dl(1),u)
                           call set_bc(cbc(0,1,2),0,1,nh,.true. ,bc(0,1,2),dl(1),v)
                           call set_bc(cbc(0,1,3),0,1,nh,.true. ,bc(0,1,3),dl(1),w)
      end if
    end if
    if(impose_norm_bc .and. cbc(1,1,1) == 'A') then
      if(.not.(present(istep).and.present(rkpar).and.present(dt).and.present(lo).and.present(ng))) then
        if(myid == 0) print*, 'ERROR: OUTFLOW_ADV BC requires present(istep,rkpar,dt,lo,ng) in bounduvw().'
        error stop
      end if
      call advective_outflow(n,ng,lo,is_bound,istep,rkpar,dt,dl,u)
    end if
    if(is_bound(1,1)) then
      if (cbc(1,1,1) /= 'A') then
        if(impose_norm_bc) call set_bc(cbc(1,1,1),1,1,nh,.false.,bc(1,1,1),dl(1),u)
                           call set_bc(cbc(1,1,2),1,1,nh,.true. ,bc(1,1,2),dl(1),v)
                           call set_bc(cbc(1,1,3),1,1,nh,.true. ,bc(1,1,3),dl(1),w)
      else
                           call set_bc(cbc(1,1,2),1,1,nh,.true. ,bc(1,1,2),dl(1),v)
                           call set_bc(cbc(1,1,3),1,1,nh,.true. ,bc(1,1,3),dl(1),w)
      end if
    end if
    impose_norm_bc = (.not.is_correc).or.(cbc(0,2,2)//cbc(1,2,2) == 'PP')
    if(is_bound(0,2)) then
                         call set_bc(cbc(0,2,1),0,2,nh,.true. ,bc(0,2,1),dl(2),u)
      if(impose_norm_bc) call set_bc(cbc(0,2,2),0,2,nh,.false.,bc(0,2,2),dl(2),v)
                         call set_bc(cbc(0,2,3),0,2,nh,.true. ,bc(0,2,3),dl(2),w)
    end if
    if(is_bound(1,2)) then
                         call set_bc(cbc(1,2,1),1,2,nh,.true. ,bc(1,2,1),dl(2),u)
      if(impose_norm_bc) call set_bc(cbc(1,2,2),1,2,nh,.false.,bc(1,2,2),dl(2),v)
                         call set_bc(cbc(1,2,3),1,2,nh,.true. ,bc(1,2,3),dl(2),w)
    end if
    impose_norm_bc = (.not.is_correc).or.(cbc(0,3,3)//cbc(1,3,3) == 'PP')
    if(is_bound(0,3)) then
                         call set_bc(cbc(0,3,1),0,3,nh,.true. ,bc(0,3,1),dzc(0)   ,u)
                         call set_bc(cbc(0,3,2),0,3,nh,.true. ,bc(0,3,2),dzc(0)   ,v)
      if(impose_norm_bc) call set_bc(cbc(0,3,3),0,3,nh,.false.,bc(0,3,3),dzf(0)   ,w)
    end if
    if(is_bound(1,3)) then
                         call set_bc(cbc(1,3,1),1,3,nh,.true. ,bc(1,3,1),dzc(n(3)),u)
                         call set_bc(cbc(1,3,2),1,3,nh,.true. ,bc(1,3,2),dzc(n(3)),v)
      if(impose_norm_bc) call set_bc(cbc(1,3,3),1,3,nh,.false.,bc(1,3,3),dzf(n(3)),w)
    end if
  end subroutine bounduvw
  !
  subroutine boundp(cbc,n,bc,nb,is_bound,dl,dzc,p)
    !
    ! imposes pressure boundary conditions
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:1,3) :: bc
    integer , intent(in), dimension(0:1,3) :: nb
    logical , intent(in), dimension(0:1,3) :: is_bound
    real(rp), intent(in), dimension(3 ) :: dl
    real(rp), intent(in), dimension(0:) :: dzc
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    integer :: idir,nh
    !
    nh = 1
    !
#if !defined(_OPENACC)
    do idir = 1,3
      call updthalo(nh,halo(idir),nb(:,idir),idir,p)
    end do
#else
    call updthalo_gpu(nh,cbc(0,:)//cbc(1,:)==['PP','PP','PP'],p)
#endif
    !
    if(is_bound(0,1)) then
      call set_bc(cbc(0,1),0,1,nh,.true.,bc(0,1),dl(1),p)
    end if
    if(is_bound(1,1)) then
      call set_bc(cbc(1,1),1,1,nh,.true.,bc(1,1),dl(1),p)
    end if
    if(is_bound(0,2)) then
      call set_bc(cbc(0,2),0,2,nh,.true.,bc(0,2),dl(2),p)
     end if
    if(is_bound(1,2)) then
      call set_bc(cbc(1,2),1,2,nh,.true.,bc(1,2),dl(2),p)
    end if
    if(is_bound(0,3)) then
      call set_bc(cbc(0,3),0,3,nh,.true.,bc(0,3),dzc(0)   ,p)
    end if
    if(is_bound(1,3)) then
      call set_bc(cbc(1,3),1,3,nh,.true.,bc(1,3),dzc(n(3)),p)
    end if
  end subroutine boundp
  !
  subroutine set_bc(ctype,ibound,idir,nh,centered,rvalue,dr,p)
    implicit none
    character(len=1), intent(in) :: ctype
    integer , intent(in) :: ibound,idir,nh
    logical , intent(in) :: centered
    real(rp), intent(in) :: rvalue,dr
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp) :: factor,sgn
    integer  :: n,dh
    integer  :: i,j,k
    !
    n = size(p,idir) - 2*nh
    factor = rvalue
    if(ctype == 'D'.and.centered) then
      factor = 2.*factor
      sgn    = -1.
    end if
    if(ctype == 'N') then
      if(     ibound == 0) then
        factor = -dr*factor ! n.b.: only valid for nh /= 1 or factor /= 0
      else if(ibound == 1) then
        factor =  dr*factor ! n.b.: only valid for nh /= 1 or factor /= 0
      end if
      sgn    = 1.
    end if
    !
    do dh=0,nh-1
      select case(ctype)
      case('P')
        !
        ! n.b.: this periodic BC imposition assumes that the subroutine is only called for
        !       for non-decomposed directions, for which n is the domain length in index space;
        !       note that the is_bound(:,:) mask above (set under initmpi.f90) is only true along
        !       the (undecomposed) pencil direction;
        !       along decomposed directions, periodicity is naturally set via the halo exchange
        !
        select case(idir)
        case(1)
          !$acc parallel loop collapse(2) default(present) async(1)
          !$OMP parallel do   collapse(2) DEFAULT(shared)
         do k=1-nh,size(p,3)-nh
           do j=1-nh,size(p,2)-nh
              p(  0-dh,j,k) = p(n-dh,j,k)
              p(n+1+dh,j,k) = p(1+dh,j,k)
            end do
          end do
        case(2)
          !$acc parallel loop collapse(2) default(present) async(1)
          !$OMP parallel do   collapse(2) DEFAULT(shared)
          do k=1-nh,size(p,3)-nh
            do i=1-nh,size(p,1)-nh
              p(i,  0-dh,k) = p(i,n-dh,k)
              p(i,n+1+dh,k) = p(i,1+dh,k)
            end do
          end do
        case(3)
          !$acc parallel loop collapse(2) default(present) async(1)
          !$OMP parallel do   collapse(2) DEFAULT(shared)
          do j=1-nh,size(p,2)-nh
            do i=1-nh,size(p,1)-nh
              p(i,j,  0-dh) = p(i,j,n-dh)
              p(i,j,n+1+dh) = p(i,j,1+dh)
            end do
          end do
        end select
      case('D','N')
        if(centered) then
          select case(idir)
          case(1)
            if     (ibound == 0) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do j=1-nh,size(p,2)-nh
                  p(  0-dh,j,k) = factor+sgn*p(1+dh,j,k)
                end do
              end do
            else if(ibound == 1) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do j=1-nh,size(p,2)-nh
                  p(n+1+dh,j,k) = factor+sgn*p(n-dh,j,k)
                end do
              end do
            end if
          case(2)
            if     (ibound == 0) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do i=1-nh,size(p,1)-nh
                  p(i,  0-dh,k) = factor+sgn*p(i,1+dh,k)
                end do
              end do
            else if(ibound == 1) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do i=1-nh,size(p,1)-nh
                  p(i,n+1+dh,k) = factor+sgn*p(i,n-dh,k)
                end do
              end do
            end if
          case(3)
            if     (ibound == 0) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do j=1-nh,size(p,2)-nh
                do i=1-nh,size(p,1)-nh
                  p(i,j,  0-dh) = factor+sgn*p(i,j,1+dh)
                end do
              end do
            else if(ibound == 1) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do j=1-nh,size(p,2)-nh
                do i=1-nh,size(p,1)-nh
                  p(i,j,n+1+dh) = factor+sgn*p(i,j,n-dh)
                end do
              end do
            end if
          end select
        else if(.not.centered.and.ctype == 'D') then
          select case(idir)
          case(1)
            if     (ibound == 0) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do j=1-nh,size(p,2)-nh
                  p(0-dh,j,k) = factor
                end do
              end do
            else if(ibound == 1) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do j=1-nh,size(p,2)-nh
                  p(n+1,j,k) = p(n-1,j,k) ! unused
                  p(n+dh,j,k) = factor
                end do
              end do
            end if
          case(2)
            if     (ibound == 0) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do i=1-nh,size(p,1)-nh
                  p(i,0-dh,k) = factor
                end do
              end do
            else if(ibound == 1) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do i=1-nh,size(p,1)-nh
                  p(i,n+1,k) = p(i,n-1,k) ! unused
                  p(i,n+dh,k) = factor
                end do
              end do
            end if
          case(3)
            if     (ibound == 0) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do j=1-nh,size(p,2)-nh
                do i=1-nh,size(p,1)-nh
                  p(i,j,0-dh) = factor
                end do
              end do
            else if(ibound == 1) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do j=1-nh,size(p,2)-nh
                do i=1-nh,size(p,1)-nh
                  p(i,j,n+1) = p(i,j,n-1) ! unused
                  p(i,j,n+dh) = factor
                end do
              end do
            end if
          end select
        else if(.not.centered.and.ctype == 'N') then
          select case(idir)
          case(1)
            if     (ibound == 0) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do j=1-nh,size(p,2)-nh
                  !p(0-dh,j,k) = 1./3.*(-2.*factor+4.*p(1+dh,j,k)-p(2+dh,j,k))
                  p(0-dh,j,k) = 1.*factor + p(  1+dh,j,k)
                end do
              end do
            else if(ibound == 1) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do j=1-nh,size(p,2)-nh
                  !p(n+1,j,k) = 1./3.*(-2.*factor+4.*p(n-1,j,k)-p(n-2,j,k))
                  p(n+1,j,k) = p(n,j,k) ! unused
                  p(n+dh,j,k) = 1.*factor + p(n-1-dh,j,k)
                end do
              end do
            end if
          case(2)
            if     (ibound == 0) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do i=1-nh,size(p,1)-nh
                  !p(i,0-dh,k) = 1./3.*(-2.*factor+4.*p(i,1+dh,k)-p(i,2+dh,k))
                  p(i,0-dh,k) = 1.*factor + p(i,  1+dh,k)
                end do
              end do
            else if(ibound == 1) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do i=1-nh,size(p,1)-nh
                  !p(i,n+1,k) = 1./3.*(-2.*factor+4.*p(i,n-1,k)-p(i,n-2,k))
                  p(i,n+1,k) = p(i,n,k) ! unused
                  p(i,n+dh,k) = 1.*factor + p(i,n-1-dh,k)
                end do
              end do
            end if
          case(3)
            if     (ibound == 0) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do j=1-nh,size(p,2)-nh
                do i=1-nh,size(p,1)-nh
                  !p(i,j,0-dh) = 1./3.*(-2.*factor+4.*p(i,j,1+dh)-p(i,j,2+dh))
                  p(i,j,0-dh) = 1.*factor + p(i,j,  1+dh)
                end do
              end do
            else if(ibound == 1) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do j=1-nh,size(p,2)-nh
                do i=1-nh,size(p,1)-nh
                  !p(i,j,n+1) = 1./3.*(-2.*factor+4.*p(i,j,n-1)-p(i,j,n-2))
                  p(i,j,n+1) = p(i,j,n) ! unused
                  p(i,j,n+dh) = 1.*factor + p(i,j,n-1-dh)
                end do
              end do
            end if
          end select
        end if
      end select
    end do
  end subroutine set_bc
  !
  subroutine inflow(idir,is_bound,vel2d,u,v,w)
    implicit none
    integer , intent(in   )  :: idir
    logical , intent(in   ), dimension(0:1,3) :: is_bound
    real(rp), intent(in   ), dimension(0:,0:   ) :: vel2d
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    integer :: i,j,k
    integer, dimension(3) :: n
    !
    select case(idir)
      case(1) ! x direction
        if(is_bound(0,1)) then
          n(:) = shape(u) - 2*1
          i = 0
          !$acc parallel loop collapse(2) default(present) async(1)
          !$OMP parallel do   collapse(2) DEFAULT(shared)
          do k=1,n(3)
            do j=1,n(2)
              u(i,j,k) = vel2d(j,k)
            end do
          end do
        end if
      case(2) ! y direction
        if(is_bound(0,2)) then
          n(:) = shape(v) - 2*1
          j = 0
          !$acc parallel loop collapse(2) default(present) async(1)
          !$OMP parallel do   collapse(2) DEFAULT(shared)
          do k=1,n(3)
            do i=1,n(1)
              v(i,j,k) = vel2d(i,k)
            end do
          end do
        end if
      case(3) ! z direction
        if(is_bound(0,3)) then
          n(:) = shape(w) - 2*1
          k = 0
          !$acc parallel loop collapse(2) default(present) async(1)
          !$OMP parallel do   collapse(2) DEFAULT(shared)
          do j=1,n(2)
            do i=1,n(1)
              w(i,j,k) = vel2d(i,j)
            end do
          end do
        end if
    end select
  end subroutine inflow
  !
  subroutine updt_rhs_b(c_or_f,cbc,n,is_bound,rhsbx,rhsby,rhsbz,p,alpha)
    implicit none
    character(len=1), intent(in), dimension(3    ) :: c_or_f
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer , intent(in), dimension(3) :: n
    logical , intent(in), dimension(0:1,3) :: is_bound
    real(rp), intent(in), dimension(:,:,0:), optional :: rhsbx,rhsby,rhsbz
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    real(rp), intent(in), optional :: alpha
    integer , dimension(3) :: q
    integer :: idir
    integer :: nn
    integer :: i,j,k
    real(rp) :: norm
    q(:) = 0
    do idir = 1,3
      !---------------- BLASIUS ----------------
      if(c_or_f(idir) == 'f'.and.(cbc(1,idir) == 'D'.or.cbc(1,idir) == 'B')) q(idir) = 1
      !-------------- END BLASIUS --------------
    end do
    norm = 1.
    if(present(alpha)) norm = alpha
    !
    if(present(rhsbx)) then
      if(is_bound(0,1)) then
        !$acc parallel loop collapse(2) default(present) async(1)
        !$OMP parallel do   collapse(2) DEFAULT(shared)
        do k=1,n(3)
          do j=1,n(2)
            p(1 ,j,k) = p(1 ,j,k) + rhsbx(j,k,0)*norm
          end do
        end do
      end if
      if(is_bound(1,1)) then
        nn = n(1)-q(1)
        !$acc parallel loop collapse(2) default(present) async(1)
        !$OMP parallel do   collapse(2) DEFAULT(shared)
        do k=1,n(3)
          do j=1,n(2)
            p(nn,j,k) = p(nn,j,k) + rhsbx(j,k,1)*norm
          end do
        end do
      end if
    end if
    if(present(rhsby)) then
      if(is_bound(0,2)) then
        !$acc parallel loop collapse(2) default(present) async(1)
        !$OMP parallel do   collapse(2) DEFAULT(shared)
        do k=1,n(3)
          do i=1,n(1)
            p(i,1 ,k) = p(i,1 ,k) + rhsby(i,k,0)*norm
          end do
        end do
      end if
      if(is_bound(1,2)) then
        nn = n(2)-q(2)
        !$acc parallel loop collapse(2) default(present) async(1)
        !$OMP parallel do   collapse(2) DEFAULT(shared)
        do k=1,n(3)
          do i=1,n(1)
            p(i,nn,k) = p(i,nn,k) + rhsby(i,k,1)*norm
          end do
        end do
      end if
    end if
    if(present(rhsbz)) then
      if(is_bound(0,3)) then
        !$acc parallel loop collapse(2) default(present) async(1)
        !$OMP parallel do   collapse(2) DEFAULT(shared)
        do j=1,n(2)
          do i=1,n(1)
            p(i,j,1 ) = p(i,j,1 ) + rhsbz(i,j,0)*norm
          end do
        end do
      end if
      if(is_bound(1,3)) then
        nn = n(3)-q(3)
        !$acc parallel loop collapse(2) default(present) async(1)
        !$OMP parallel do   collapse(2) DEFAULT(shared)
        do j=1,n(2)
          do i=1,n(1)
            p(i,j,nn) = p(i,j,nn) + rhsbz(i,j,1)*norm
          end do
        end do
      end if
    end if
  end subroutine updt_rhs_b
  !
  subroutine updthalo(nh,halo,nb,idir,p)
    implicit none
    integer , intent(in) :: nh ! number of ghost points
    integer , intent(in) :: halo
    integer , intent(in), dimension(0:1) :: nb
    integer , intent(in) :: idir
    real(rp), dimension(1-nh:,1-nh:,1-nh:), intent(inout) :: p
    integer , dimension(3) :: lo,hi
#if defined(_ASYNC_HALO)
    integer :: requests(4)
#endif
    !
    !  this subroutine updates the halo that store info
    !  from the neighboring computational sub-domain
    !
    if(idir == ipencil_axis) return
    lo(:) = lbound(p)+nh
    hi(:) = ubound(p)-nh
    select case(idir)
    case(1) ! x direction
#if !defined(_ASYNC_HALO)
      call MPI_SENDRECV(p(lo(1)     ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),0, &
                        p(hi(1)+1   ,lo(2)-nh,lo(3)-nh),1,halo,nb(1),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      call MPI_SENDRECV(p(hi(1)-nh+1,lo(2)-nh,lo(3)-nh),1,halo,nb(1),0, &
                        p(lo(1)-nh  ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
#else
      call MPI_IRECV( p(hi(1)+1  ,lo(2)-nh,lo(3)-nh),1,halo,nb(1),0, &
                      MPI_COMM_WORLD,requests(1),ierr)
      call MPI_IRECV( p(lo(1)-nh ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),1, &
                      MPI_COMM_WORLD,requests(2),ierr)
      call MPI_ISEND(p(lo(1)     ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),0, &
                      MPI_COMM_WORLD,requests(3),ierr)
      call MPI_ISEND(p(hi(1)-nh+1,lo(2)-nh,lo(3)-nh),1,halo,nb(1),1, &
                      MPI_COMM_WORLD,requests(4),ierr)
      call MPI_WAITALL(4,requests,MPI_STATUSES_IGNORE,ierr)
#endif
    case(2) ! y direction
#if !defined(_ASYNC_HALO)
      call MPI_SENDRECV(p(lo(1)-nh,lo(2)     ,lo(3)-nh),1,halo,nb(0),0, &
                        p(lo(1)-nh,hi(2)+1   ,lo(3)-nh),1,halo,nb(1),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      call MPI_SENDRECV(p(lo(1)-nh,hi(2)-nh+1,lo(3)-nh),1,halo,nb(1),0, &
                        p(lo(1)-nh,lo(2)-nh  ,lo(3)-nh),1,halo,nb(0),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
#else
      call MPI_IRECV(p(lo(1)-nh,hi(2)+1   ,lo(3)-nh),1,halo,nb(1),0, &
                      MPI_COMM_WORLD,requests(1),ierr)
      call MPI_IRECV(p(lo(1)-nh,lo(2)-nh  ,lo(3)-nh),1,halo,nb(0),1, &
                      MPI_COMM_WORLD,requests(2),ierr)
      call MPI_ISEND(p(lo(1)-nh,lo(2)     ,lo(3)-nh),1,halo,nb(0),0, &
                      MPI_COMM_WORLD,requests(3),ierr)
      call MPI_ISEND(p(lo(1)-nh,hi(2)-nh+1,lo(3)-nh),1,halo,nb(1),1, &
                      MPI_COMM_WORLD,requests(4),ierr)
      call MPI_WAITALL(4,requests,MPI_STATUSES_IGNORE,ierr)
#endif
    case(3) ! z direction
#if !defined(_ASYNC_HALO)
      call MPI_SENDRECV(p(lo(1)-nh,lo(2)-nh,lo(3)     ),1,halo,nb(0),0, &
                        p(lo(1)-nh,lo(2)-nh,hi(3)+1   ),1,halo,nb(1),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      call MPI_SENDRECV(p(lo(1)-nh,lo(2)-nh,hi(3)-nh+1),1,halo,nb(1),0, &
                        p(lo(1)-nh,lo(2)-nh,lo(3)-nh  ),1,halo,nb(0),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
#else
      call MPI_IRECV(p(lo(1)-nh,lo(2)-nh,hi(3)+1   ),1,halo,nb(1),0, &
                      MPI_COMM_WORLD,requests(1),ierr)
      call MPI_IRECV(p(lo(1)-nh,lo(2)-nh,lo(3)-nh  ),1,halo,nb(0),1, &
                      MPI_COMM_WORLD,requests(2),ierr)
      call MPI_ISEND(p(lo(1)-nh,lo(2)-nh,lo(3)     ),1,halo,nb(0),0, &
                      MPI_COMM_WORLD,requests(3),ierr)
      call MPI_ISEND(p(lo(1)-nh,lo(2)-nh,hi(3)-nh+1),1,halo,nb(1),1, &
                      MPI_COMM_WORLD,requests(4),ierr)
      call MPI_WAITALL(4,requests,MPI_STATUSES_IGNORE,ierr)
#endif
    end select
  end subroutine updthalo
#if defined(_OPENACC)
  subroutine updthalo_gpu(nh,periods,p)
    use mod_types
#if !defined(_USE_DIEZDECOMP)
    use cudecomp
#else
    use diezdecomp
#endif
    use mod_common_cudecomp, only: work => work_halo, &
                                   ch => handle,gd => gd_halo, &
                                   dtype => cudecomp_real_rp, &
                                   istream => istream_acc_queue_1_comm_lib
    implicit none
    integer , intent(in) :: nh
    logical , intent(in) :: periods(3)
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    integer :: istat
#if !defined(_USE_DIEZDECOMP)
    !$acc host_data use_device(p,work)
#endif
    select case(ipencil_axis)
    case(1)
      istat = cudecompUpdateHalosX(ch,gd,p,work,dtype,[nh,nh,nh],periods,2,stream=istream)
      istat = cudecompUpdateHalosX(ch,gd,p,work,dtype,[nh,nh,nh],periods,3,stream=istream)
    case(2)
      istat = cudecompUpdateHalosY(ch,gd,p,work,dtype,[nh,nh,nh],periods,1,stream=istream)
      istat = cudecompUpdateHalosY(ch,gd,p,work,dtype,[nh,nh,nh],periods,3,stream=istream)
    case(3)
      istat = cudecompUpdateHalosZ(ch,gd,p,work,dtype,[nh,nh,nh],periods,1,stream=istream)
      istat = cudecompUpdateHalosZ(ch,gd,p,work,dtype,[nh,nh,nh],periods,2,stream=istream)
    end select
#if !defined(_USE_DIEZDECOMP)
    !$acc end host_data
#endif
  end subroutine updthalo_gpu
#endif
end module mod_bound
