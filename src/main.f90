! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
!
!        CCCCCCCCCCCCC                    NNNNNNNN        NNNNNNNN    SSSSSSSSSSSSSSS
!     CCC::::::::::::C                    N:::::::N       N::::::N  SS:::::::::::::::S
!   CC:::::::::::::::C                    N::::::::N      N::::::N S:::::SSSSSS::::::S
!  C:::::CCCCCCCC::::C                    N:::::::::N     N::::::N S:::::S     SSSSSSS
! C:::::C       CCCCCC   aaaaaaaaaaaaa    N::::::::::N    N::::::N S:::::S
!C:::::C                 a::::::::::::a   N:::::::::::N   N::::::N S:::::S
!C:::::C                 aaaaaaaaa:::::a  N:::::::N::::N  N::::::N  S::::SSSS
!C:::::C                          a::::a  N::::::N N::::N N::::::N   SS::::::SSSSS
!C:::::C                   aaaaaaa:::::a  N::::::N  N::::N:::::::N     SSS::::::::SS
!C:::::C                 aa::::::::::::a  N::::::N   N:::::::::::N        SSSSSS::::S
!C:::::C                a::::aaaa::::::a  N::::::N    N::::::::::N             S:::::S
! C:::::C       CCCCCC a::::a    a:::::a  N::::::N     N:::::::::N             S:::::S
!  C:::::CCCCCCCC::::C a::::a    a:::::a  N::::::N      N::::::::N SSSSSSS     S:::::S
!   CC:::::::::::::::C a:::::aaaa::::::a  N::::::N       N:::::::N S::::::SSSSSS:::::S
!     CCC::::::::::::C  a::::::::::aa:::a N::::::N        N::::::N S:::::::::::::::SS
!        CCCCCCCCCCCCC   aaaaaaaaaa  aaaa NNNNNNNN         NNNNNNN  SSSSSSSSSSSSSSS
!-------------------------------------------------------------------------------------
! CaNS -- Canonical Navier-Stokes Solver
!-------------------------------------------------------------------------------------
program cans
  use, intrinsic :: iso_fortran_env, only: compiler_version,compiler_options,real32
  use, intrinsic :: iso_c_binding  , only: C_PTR
  use, intrinsic :: ieee_arithmetic, only: is_nan => ieee_is_nan
  use mpi
  use decomp_2d
  use mod_bound          , only: boundp,bounduvw,updt_rhs_b, &
                                 inflow_register,compute_means
  use mod_chkdiv         , only: chkdiv
  use mod_chkdt          , only: chkdt
  use mod_common_mpi     , only: myid,ierr,dinfo_ptdma
  use mod_correc         , only: correc
  use mod_fft            , only: fftini,fftend
  use mod_fillps         , only: fillps
  use mod_initflow       , only: initflow,initscal,inflow_register_initflow
  use mod_initgrid       , only: initgrid
  use mod_initmpi        , only: initmpi
  use mod_initsolver     , only: initsolver
  use mod_load           , only: load_all
  use mod_mom            , only: bulk_forcing
  use mod_rk             , only: rk,rk_scal
  use mod_plan           , only: make_fname,read_plan_file_mpi,write_plan_file_mpi
  use mod_output         , only: out0d,gen_alias,out1d,out1d_chan,out2d,out3d,write_log_output,write_visu_2d,write_visu_3d,write_visu_3d_crop
  use mod_param          , only: ng,l,dl,dli, &
                                 gtype,gr, &
                                 at, &
                                 cfl,dtmax,dt_f, &
                                 visc,alpha_max, &
                                 inivel,is_wallturb, &
                                 nstep,time_max,tw_max,stop_type, &
                                 restart,is_overwrite_save,nsaves_max, &
                                 icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                                 cbcvel,bcvel,cbcpre,bcpre, &
                                 is_forced,bforce,velf, &
                                 gacc,nscal,beta, &
                                 dims, &
                                 nb,is_bound, &
                                 rkcoeff,small, &
                                 datadir, &
                                 read_input, &
                                 is_debug,is_debug_poisson, &
                                 is_timing, &
                                 is_impdiff,is_impdiff_1d, &
                                 is_poisson_pcr_tdma, &
                                 is_mask_divergence_check, &
                                 is_inlet_record,is_inlet_replay,ix0,save_freq,read_freq,istart, &
                                 iout3dcrop,imin,imax,jmin,jmax,kmin,kmax
  use mod_sanity         , only: test_sanity_input,test_sanity_solver
  use mod_scal           , only: scalar,initialize_scalars,bulk_forcing_s
  use mod_solve_helmholtz, only: solve_helmholtz,rhs_bound
#if !defined(_OPENACC)
  use mod_solver         , only: solver
#else
  use mod_solver_gpu     , only: solver => solver_gpu
  use mod_workspaces     , only: init_wspace_arrays,set_cufft_wspace,cudecomp_finalize
  use mod_common_cudecomp, only: istream_acc_queue_1,ap_z_ptdma
#endif
  use mod_timer          , only: timer_tic,timer_toc,timer_print
  use mod_updatep        , only: updatep
  use mod_utils          , only: bulk_mean,copy_file
#if defined(_OPENACC)
  use mod_utils          , only: device_memory_footprint
#endif
  use mod_types
  use omp_lib
  implicit none
  integer , dimension(3) :: lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z
  real(rp), target, allocatable, dimension(:,:,:) :: u,v,w,p,pp
  character(len=1), dimension(0:1,3) :: cbc_u, cbc_v, cbc_w
  real(rp), allocatable, dimension(:,:,:) :: dudtrko,dvdtrko,dwdtrko
  real(rp), dimension(0:1,3) :: tauxo,tauyo,tauzo
  real(rp), dimension(3) :: f
#if !defined(_OPENACC) || defined(_USE_HIP)
  type(C_PTR), dimension(2,2) :: arrplanp
#else
  integer    , dimension(2,2) :: arrplanp
#endif
  real(rp), allocatable, dimension(:,:) :: lambdaxyp
  real(rp), allocatable, dimension(:) :: ap,bp,cp
  integer , dimension(3) :: n_z_d
  real(rp), allocatable, dimension(:,:,:) :: ap_d,cp_d
  logical :: is_ptdma_update_p
  real(rp) :: normfftp
  type(rhs_bound) :: rhsbp
  real(rp) :: alpha
#if !defined(_OPENACC) || defined(_USE_HIP)
  type(C_PTR), dimension(2,2) :: arrplanu,arrplanv,arrplanw
#else
  integer    , dimension(2,2) :: arrplanu,arrplanv,arrplanw
#endif
  real(rp), allocatable, dimension(:,:) :: lambdaxyu,lambdaxyv,lambdaxyw
  real(rp), allocatable, dimension(:) :: au,av,aw,bu,bv,bw,cu,cv,cw
  real(rp) :: normfftu,normfftv,normfftw
  type(rhs_bound) :: rhsbu,rhsbv,rhsbw
  !
  real(rp) :: dt,dti,dt_cfl,time,dtrk,dtrki,divtot,divmax
  integer :: irk,istep
  !---------------- INLET_PLAN ----------------
  integer :: comm_inlet, color_inlet
  real(rp), allocatable :: u2d(:,:), u2dm(:,:), v2d(:,:), w2d(:,:), p2d(:,:)
  real(real32), allocatable :: u2d_pack(:,:,:), u2dm_pack(:,:,:), &
                               v2d_pack(:,:,:), w2d_pack(:,:,:), p2d_pack(:,:,:)
  real(real32), allocatable :: u2d_read_pack(:,:,:), v2d_read_pack(:,:,:), w2d_read_pack(:,:,:), p2d_read_pack(:,:,:)
  real(real32), allocatable, target :: u2d_read(:,:), v2d_read(:,:), w2d_read(:,:), p2d_read(:,:)
  character(len=512) :: fname_inlet
  logical :: is_inlet_window
  integer :: slot, itloc, cur_start, cur_end, ios_plan, ixloc
  !-------------- END INLET_PLAN --------------
  real(rp), allocatable, dimension(:) :: dzc  ,dzf  ,zc  ,zf  ,dzci  ,dzfi, &
                                         dzc_g,dzf_g,zc_g,zf_g,dzci_g,dzfi_g, &
                                         grid_vol_ratio_c,grid_vol_ratio_f
  real(rp) :: meanvelu,meanvelv,meanvelw
  real(rp), dimension(3) :: dpdl
  !
  type(scalar), target, allocatable, dimension(:) :: scalars
  type(scalar), pointer :: s
  real(rp) :: meanscal
  real(rp), allocatable, dimension(:) :: fs
  integer :: iscal
  !
  !real(rp), allocatable, dimension(:) :: var
  real(rp), dimension(42) :: var
  !
  real(rp) :: dt12,dt12av,dt12min,dt12max
  real(rp) :: twi,tw
  !
  integer  :: savecounter
  character(len=7  ) :: fldnum
  character(len=3  ) :: scalnum
  character(len=4  ) :: chkptnum
  character(len=100) :: filename,trip_filename
  integer :: k,kk,j
  logical :: is_done,kill
  integer :: gmin(3), gmax(3)
  integer :: lmin(3), lmax(3)
  integer :: comm_crop, ierr_crop, myid_crop
  integer :: color
  integer :: ng_crop(3), lo_crop(3), hi_crop(3)
  logical :: intersect,has_data

  !
  call MPI_INIT(ierr)
  call init_types_mpi(ierr)
  comm_inlet=MPI_COMM_NULL
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  !
  ! read parameter file
  !
  call read_input(myid)
  !
  ! initialize MPI/OpenMP
  !
  !$ call omp_set_num_threads(omp_get_max_threads())
  call initmpi(ng,dims,cbcpre,lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z,nb,is_bound)
  twi = MPI_WTIME()
  savecounter = 0
  !
  ! allocate variables
  !
  allocate(u( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           v( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           w( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           p( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           pp(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(dudtrko(n(1),n(2),n(3)), &
           dvdtrko(n(1),n(2),n(3)), &
           dwdtrko(n(1),n(2),n(3)))
  allocate(lambdaxyp(n_z(1),n_z(2)))
  allocate(ap(n_z(3)),bp(n_z(3)),cp(n_z(3)))
  if(is_poisson_pcr_tdma) then
#if defined(_OPENACC)
    n_z_d(:) = ap_z_ptdma%shape(:)
#else
    n_z_d(:) = dinfo_ptdma%zsz(:)
#endif
    allocate(ap_d(n_z_d(1),n_z_d(2),n_z_d(3)), &
             cp_d(n_z_d(1),n_z_d(2),n_z_d(3)))
  end if
  allocate(dzc( 0:n(3)+1), &
           dzf( 0:n(3)+1), &
           zc(  0:n(3)+1), &
           zf(  0:n(3)+1), &
           dzci(0:n(3)+1), &
           dzfi(0:n(3)+1))
  allocate(dzc_g( 0:ng(3)+1), &
           dzf_g( 0:ng(3)+1), &
           zc_g(  0:ng(3)+1), &
           zf_g(  0:ng(3)+1), &
           dzci_g(0:ng(3)+1), &
           dzfi_g(0:ng(3)+1))
  allocate(grid_vol_ratio_c,mold=dzc)
  allocate(grid_vol_ratio_f,mold=dzf)
  allocate(rhsbp%x(n(2),n(3),0:1), &
           rhsbp%y(n(1),n(3),0:1), &
           rhsbp%z(n(1),n(2),0:1))
  if(is_impdiff) then
    allocate(lambdaxyu(n_z(1),n_z(2)), &
             lambdaxyv(n_z(1),n_z(2)), &
             lambdaxyw(n_z(1),n_z(2)))
    allocate(au(n_z(3)),bu(n_z(3)),cu(n_z(3)), &
             av(n_z(3)),bv(n_z(3)),cv(n_z(3)), &
             aw(n_z(3)),bw(n_z(3)),cw(n_z(3)))
    allocate(rhsbu%x(n(2),n(3),0:1), &
             rhsbu%y(n(1),n(3),0:1), &
             rhsbu%z(n(1),n(2),0:1), &
             rhsbv%x(n(2),n(3),0:1), &
             rhsbv%y(n(1),n(3),0:1), &
             rhsbv%z(n(1),n(2),0:1), &
             rhsbw%x(n(2),n(3),0:1), &
             rhsbw%y(n(1),n(3),0:1), &
             rhsbw%z(n(1),n(2),0:1))
  end if
  !
  allocate(scalars(nscal))
  call initialize_scalars(scalars,nscal,n,n_z)
  allocate(fs(nscal))
  !$acc enter data copyin(scalars(:))
  !
  if(is_debug) then
    if(myid == 0) print*, 'This executable of CaNS was built with compiler: ', compiler_version()
    if(myid == 0) print*, 'Using the options: ', compiler_options()
    block
      character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: mpi_version
      integer :: ilen
      call MPI_GET_LIBRARY_VERSION(mpi_version,ilen,ierr)
      if(myid == 0) print*, 'MPI Version: ', trim(mpi_version)
    end block
    if(myid == 0) print*, ''
  end if
  if(myid == 0) print*, '*******************************'
  if(myid == 0) print*, '*** Beginning of simulation ***'
  if(myid == 0) print*, '*******************************'
  if(myid == 0) print*, ''
  call initgrid(gtype,ng(3),gr,l(3),dzc_g,dzf_g,zc_g,zf_g)
  if(myid == 0) then
    open(99,file=trim(datadir)//'grid.bin',action='write',form='unformatted',access='stream',status='replace')
    write(99) dzc_g(1:ng(3)),dzf_g(1:ng(3)),zc_g(1:ng(3)),zf_g(1:ng(3))
    close(99)
    open(99,file=trim(datadir)//'grid.out')
    do kk=0,ng(3)+1
      write(99,*) 0.,zf_g(kk),zc_g(kk),dzf_g(kk),dzc_g(kk)
    end do
    close(99)
    open(99,file=trim(datadir)//'geometry.out')
      write(99,*) ng(1),ng(2),ng(3)
      write(99,*) l(1),l(2),l(3)
    close(99)
  end if
  !$acc enter data copyin(lo,hi,n) async
  !$acc enter data copyin(bforce,dl,dli,l) async
  !$acc enter data copyin(gacc) async
  !$acc enter data copyin(zc_g,zf_g,dzc_g,dzf_g) async
  !$acc enter data create(zc,zf,dzc,dzf,dzci,dzfi,dzci_g,dzfi_g) async
  !
  !$acc parallel loop default(present) private(k) async
  do kk=lo(3)-1,hi(3)+1
    k = kk-(lo(3)-1)
    zc( k) = zc_g(kk)
    zf( k) = zf_g(kk)
    dzc(k) = dzc_g(kk)
    dzf(k) = dzf_g(kk)
    dzci(k) = dzc(k)**(-1)
    dzfi(k) = dzf(k)**(-1)
  end do
  !$acc data copy(ng) async
  !$acc parallel loop default(present) async
  do k=0,ng(3)+1
    dzci_g(k) = dzc_g(k)**(-1)
    dzfi_g(k) = dzf_g(k)**(-1)
  end do
  !$acc end data
  !$acc enter data create(grid_vol_ratio_c,grid_vol_ratio_f) async
  !$acc parallel loop default(present) async
  do k=0,n(3)+1
    grid_vol_ratio_c(k) = dl(1)*dl(2)*dzc(k)/(l(1)*l(2)*l(3))
    grid_vol_ratio_f(k) = dl(1)*dl(2)*dzf(k)/(l(1)*l(2)*l(3))
  end do
  !$acc update self(zc,zf,dzc,dzf,dzci,dzfi) async
  !$acc exit data copyout(zc_g,zf_g,dzc_g,dzf_g,dzci_g,dzfi_g) async ! not needed on the device
  !$acc wait
  !
  ! test input files before proceeding with the calculation
  !
  call test_sanity_input(ng,dims,stop_type,cbcvel,cbcpre,bcvel,bcpre,is_forced)
  !
  ! initialize Poisson solver
  !
  call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcpre,bcpre(:,:), &
                  lambdaxyp,['c','c','c'],ap,bp,cp,arrplanp,normfftp,rhsbp%x,rhsbp%y,rhsbp%z)
  !$acc enter data copyin(lambdaxyp,ap,bp,cp) async
  !$acc enter data copyin(rhsbp,rhsbp%x,rhsbp%y,rhsbp%z) async
  if(is_poisson_pcr_tdma) then
    !$acc enter data create(ap_d,cp_d) async
  end if
  !$acc wait
  if(is_impdiff) then
    ! Implicit diffusion solver only understands P/D/N. Keep these mapped
    ! boundary conditions consistent with the ones used at solve time.
    cbc_u = cbcvel(:,:,1)
    cbc_v = cbcvel(:,:,2)
    cbc_w = cbcvel(:,:,3)
    !if(cbc_u(0,1) == 'I' .or. cbc_u(0,1) == 'B') cbc_u(0,1) = 'D'
    !if(cbc_u(1,1) == 'A') cbc_u(1,1) = 'N'
    call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbc_u,bcvel(:,:,1), &
                    lambdaxyu,['f','c','c'],au,bu,cu,arrplanu,normfftu,rhsbu%x,rhsbu%y,rhsbu%z)
    call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbc_v,bcvel(:,:,2), &
                    lambdaxyv,['c','f','c'],av,bv,cv,arrplanv,normfftv,rhsbv%x,rhsbv%y,rhsbv%z)
    call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbc_w,bcvel(:,:,3), &
                    lambdaxyw,['c','c','f'],aw,bw,cw,arrplanw,normfftw,rhsbw%x,rhsbw%y,rhsbw%z)
    do iscal=1,nscal
      s => scalars(iscal)
      call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,s%cbc,s%bc, &
                      s%lambdaxy,['c','c','c'],s%a,s%b,s%c,s%arrplan,s%normfft, &
                      s%rhsb%x,s%rhsb%y,s%rhsb%z)
    end do
    if(is_impdiff_1d) then
      deallocate(lambdaxyu,lambdaxyv,lambdaxyw)
      call fftend(arrplanu)
      call fftend(arrplanv)
      call fftend(arrplanw)
      deallocate(rhsbu%x,rhsbu%y,rhsbv%x,rhsbv%y,rhsbw%x,rhsbw%y)
      do iscal=1,nscal
        s => scalars(iscal)
        deallocate(s%rhsb%x,s%rhsb%y)
        deallocate(s%lambdaxy)
        call fftend(s%arrplan)
      end do
    end if
    !$acc enter data copyin(au,bu,cu,av,bv,cv,aw,bw,cw) async
    if(.not.is_impdiff_1d) then
      !$acc enter data copyin(lambdaxyu,lambdaxyv,lambdaxyw) async
      !$acc enter data copyin(rhsbu,rhsbu%x,rhsbu%y,rhsbu%z) async
      !$acc enter data copyin(rhsbv,rhsbv%x,rhsbv%y,rhsbv%z) async
      !$acc enter data copyin(rhsbw,rhsbw%x,rhsbw%y,rhsbw%z) async
    else
      !$acc enter data copyin(rhsbu,rhsbu%z) async
      !$acc enter data copyin(rhsbv,rhsbv%z) async
      !$acc enter data copyin(rhsbw,rhsbw%z) async
    end if
    do iscal=1,nscal
      s => scalars(iscal)
      !$acc enter data copyin(s%a,s%b,s%c) async
      if(.not.is_impdiff_1d) then
        !$acc enter data copyin(s%lambdaxy) async
        !$acc enter data copyin(s%rhsb,s%rhsb%x,s%rhsb%y,s%rhsb%z) async
      else
        !$acc enter data copyin(s%rhsb,s%rhsb%z) async
      end if
    end do
    !$acc wait
  end if
#if defined(_OPENACC)
  !
  ! determine workspace sizes and allocate the memory
  !
  call init_wspace_arrays()
  call set_cufft_wspace(pack(arrplanp,.true.),istream_acc_queue_1)
  if(is_impdiff .and. .not.is_impdiff_1d) then
    call set_cufft_wspace(pack(arrplanu,.true.),istream_acc_queue_1)
    call set_cufft_wspace(pack(arrplanv,.true.),istream_acc_queue_1)
    call set_cufft_wspace(pack(arrplanw,.true.),istream_acc_queue_1)
  end if
  if(myid == 0) print*,'*** Device memory footprint (Gb): ', &
                  device_memory_footprint(n,n_z,nscal)/(1._sp*1024**3), ' ***'
#endif
  if(is_debug_poisson) then
    call test_sanity_solver(ng,lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z,dli,dzc,dzf,dzci,dzfi,dzci_g,dzfi_g, &
                            nb,is_bound,cbcvel,cbcpre,bcvel,bcpre)
  end if
  !
  is_ptdma_update_p = .true.
  !
  if(.not.restart) then
    istep = 0
    time = 0.
  else
    call load_all('r',trim(datadir)//'fld.bin', &
                  MPI_COMM_WORLD,ng,[1,1,1],lo,hi,nscal,u,v,w,p,scalars,time,istep)
    if(myid == 0) print*, '*** Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
  end if
  !---------------- INLET_PLAN ----------------
  if (is_inlet_record .or. is_inlet_replay .or. trim(inivel) == 'inl') then
    allocate(u2d(n(2),n(3)), u2dm(n(2),n(3)), v2d(n(2),n(3)), w2d(n(2),n(3)), p2d(n(2),n(3)))
    !$acc enter data create(u2d,u2dm,v2d,w2d,p2d)
  end if
  if (is_inlet_record) then
  if (save_freq <= 0) then
    if (myid == 0) print*, 'ERROR: is_inlet_record requires save_freq > 0.'
    error stop
  end if

  is_inlet_window = (ix0 >= lo(1)) .and. (ix0 <= hi(1))
  color_inlet = merge(1, MPI_UNDEFINED, is_inlet_window)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, color_inlet, myid, comm_inlet, ierr)

  if (is_inlet_window) then
    allocate(u2d_pack(n(2),n(3),save_freq), u2dm_pack(n(2),n(3),save_freq), &
             v2d_pack(n(2),n(3),save_freq), &
             w2d_pack(n(2),n(3),save_freq), p2d_pack(n(2),n(3),save_freq))
    u2d_pack = 0.0_real32
    u2dm_pack = 0.0_real32
    v2d_pack = 0.0_real32
    w2d_pack = 0.0_real32
    p2d_pack = 0.0_real32
    end if
  end if
  if (is_inlet_replay .or. trim(inivel) == 'inl') then
    if (read_freq <= 0) then
      if (myid == 0) print*, 'ERROR: is_inlet_replay/inivel=inl requires read_freq > 0.'
      error stop
    end if
    allocate(u2d_read_pack(n(2),n(3),read_freq), v2d_read_pack(n(2),n(3),read_freq), &
             w2d_read_pack(n(2),n(3),read_freq), p2d_read_pack(n(2),n(3),read_freq))
    allocate(u2d_read(n(2),n(3)), v2d_read(n(2),n(3)), w2d_read(n(2),n(3)), p2d_read(n(2),n(3)))
    cur_start = istart
    cur_end   = istart + read_freq - 1
    call make_fname(cur_start, read_freq, ix0, fname_inlet)
    call read_plan_file_mpi(fname_inlet, MPI_COMM_WORLD, read_freq, &
                        ng(2), ng(3), lo(2), lo(3), n(2), n(3), &
                        u2d_read_pack, v2d_read_pack, w2d_read_pack, p2d_read_pack, ios_plan)
    if (ios_plan /= 0) then
      if (myid == 0) print*, 'ERROR: cannot read inlet plan file: ', trim(fname_inlet)
      error stop
    end if
    u2d_read(:,:) = u2d_read_pack(:,:,1)
    v2d_read(:,:) = v2d_read_pack(:,:,1)
    w2d_read(:,:) = w2d_read_pack(:,:,1)
    p2d_read(:,:) = p2d_read_pack(:,:,1)
    !$acc enter data copyin(u2d_read,v2d_read,w2d_read,p2d_read)
    call inflow_register(u2d_read, v2d_read, w2d_read, p2d_read)
    !call compute_means(u2d_read_pack, v2d_read_pack, w2d_read_pack)
    if (.not.restart .and. trim(inivel) == 'inl') then
      call inflow_register_initflow(u2d_read, v2d_read, w2d_read, p2d_read)
    end if
  end if
  !-------------- END INLET_PLAN --------------
  if(.not.restart) then
    call initflow(inivel,bcvel,ng,lo,l,dl,zc,zf,dzc,dzf,visc,is_forced,velf,bforce,is_wallturb,u,v,w,p)
    do iscal=1,nscal
      s => scalars(iscal)
      call initscal(s%ini,s%bc,ng,lo,l,dl,zc,dzf,s%alpha,s%is_forced,s%scalf,s%val)
    end do
    if(myid == 0) print*, '*** Initial condition succesfully set ***'
  end if
  !$acc enter data copyin(u,v,w,p,dudtrko,dvdtrko,dwdtrko) create(pp)
  !
  !---------------- OUTFLOW_ADV ----------------
  ! compute initial dt after fields are present on device (chkdt uses PRESENT)
  call chkdt(n,dl,dzci,dzfi,visc,alpha_max,u,v,w,dt_cfl)
  dt = merge(dt_f,min(cfl*dt_cfl,dtmax),dt_f > 0.)
  if(myid == 0) print*, 'dt_cfl = ', dt_cfl, 'dt = ', dt
  dti = 1./dt
  !-------------- END OUTFLOW_ADV --------------
  call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,u,v,w,istep,rkcoeff(:,1),dt,visc,lo,ng,zc)
  call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,p)
  do iscal=1,nscal
    s => scalars(iscal)
    !$acc enter data copyin(s%val,s%dsdtrko) async(1)
    call boundp(s%cbc,n,s%bc,nb,is_bound,dl,dzc,s%val)
  end do
  !$acc wait
  !
  ! post-process and write initial condition
  !
  write(fldnum,'(i7.7)') istep
  !$acc wait ! not needed but to prevent possible future issues
  !$acc update self(u,v,w,p)
  do iscal=1,nscal
    !$acc update self(scalars(iscal)%val)
  end do
  if(iout1d > 0.and.mod(istep,max(iout1d,1)) == 0) then
    include 'out1d.h90'
  end if
  if(iout2d > 0.and.mod(istep,max(iout2d,1)) == 0) then
    include 'out2d.h90'
  end if
  if(iout3d > 0.and.mod(istep,max(iout3d,1)) == 0 ) then
    include 'out3d.h90'
  end if
  if(iout3dcrop > 0.and.mod(istep,max(iout3dcrop,1)) == 0 .and. istep>20) then
      !$acc wait
      !$acc update self(u,v,w,p)
      block
        use mpi
        real(rp), pointer, contiguous :: uc(:,:,:)

        include 'out3dcrop.h90'
      end block
    end if

  !
  kill = .false.
  !
  ! main loop
  !
  if(myid == 0) print*, '*** Calculation loop starts now ***'
  is_done = .false.
  do while(.not.is_done)
    if(is_timing) then
      !$acc wait(1)
      dt12 = MPI_WTIME()
    end if
    istep = istep + 1
    time = time + dt
    if(myid == 0) print*, 'Time step #', istep, 'Time = ', time
    tauxo(:,:) = 0.; tauyo(:,:) = 0.; tauzo(:,:) = 0.
    dpdl(:)     = 0.
    fs(1:nscal) = 0.
    !---------------- INLET_PLAN ----------------
    if (is_inlet_replay) then
      slot = 1 + mod(istep-1, read_freq)
      u2d_read(:,:) = u2d_read_pack(:,:,slot)
      v2d_read(:,:) = v2d_read_pack(:,:,slot)
      w2d_read(:,:) = w2d_read_pack(:,:,slot)
      p2d_read(:,:) = p2d_read_pack(:,:,slot)
      !$acc update device(u2d_read,v2d_read,w2d_read,p2d_read)
      if (mod(istep, read_freq) == 0 .and. istep/=nstep) then
        cur_start = cur_end + 1
        cur_end   = cur_start + read_freq - 1
        call make_fname(cur_start, read_freq, ix0, fname_inlet)
        call read_plan_file_mpi(fname_inlet, MPI_COMM_WORLD, read_freq, &
                        ng(2), ng(3), lo(2), lo(3), n(2), n(3), &
                        u2d_read_pack, v2d_read_pack, w2d_read_pack, p2d_read_pack, ios_plan)
        if (ios_plan /= 0) then
          if (myid == 0) print*, 'ERROR: cannot read inlet plan file: ', trim(fname_inlet)
          error stop
        end if
        !call compute_means(u2d_read_pack, v2d_read_pack, w2d_read_pack)
      end if
    end if
    !-------------- END INLET_PLAN --------------
    do irk=1,3
      dtrk = sum(rkcoeff(:,irk))*dt
      dtrki = dtrk**(-1)
      do iscal=1,nscal
        s => scalars(iscal)
        call rk_scal(rkcoeff(:,irk),n,dli,l,dzci,dzfi,grid_vol_ratio_f,s%alpha,dt,is_bound,u,v,w, &
                     s%is_forced,s%scalf,s%source,s%fluxo,s%dsdtrko,s%val,s%f)
        call bulk_forcing_s(n,s%is_forced,s%f,s%val)
        fs(iscal) = fs(iscal) + s%f
        if(is_impdiff) then
          if(is_impdiff_1d) then
            call solve_helmholtz(n,ng,hi,alpha=-0.5*s%alpha*dtrk, &
                                 a=s%a,b=s%b,c=s%c,rhsbz=s%rhsb%z,is_bound=is_bound, &
                                 cbc=s%cbc,c_or_f=['c','c','c'],p=s%val)
          else
            call solve_helmholtz(n,ng,hi,s%arrplan,s%normfft,-0.5*s%alpha*dtrk, &
                                 s%lambdaxy,s%a,s%b,s%c,s%rhsb%x,s%rhsb%y,s%rhsb%z,is_bound,s%cbc,['c','c','c'],s%val)
          end if
        end if
        call boundp(s%cbc,n,s%bc,nb,is_bound,dl,dzc,s%val)
      end do
      call rk(rkcoeff(:,irk),n,dli,dzci,dzfi,grid_vol_ratio_c,grid_vol_ratio_f,visc,dt,p, &
              is_forced,velf,bforce,gacc,beta,scalars,dudtrko,dvdtrko,dwdtrko,u,v,w,f, &
              istep,time,lo,ng,l(2),zc)
      call bulk_forcing(n,is_forced,f,u,v,w)
      dpdl(:) = dpdl(:) + f(:)
      if(is_impdiff) then
        !---------------- CUSTOM_BC ----------------
        ! Implicit diffusion solver only understands P/D/N. When using the custom BL BCs:
        ! - inflow 'I'/'B' (prescribed) is treated as Dirichlet 'D' for the Helmholtz solve
        ! - outflow 'A' (advective) is treated as Neumann 'N' for the Helmholtz solve
        !
        cbc_u = cbcvel(:,:,1)
        cbc_v = cbcvel(:,:,2)
        cbc_w = cbcvel(:,:,3)
        if(cbc_u(0,1) == 'I' .or. cbc_u(0,1) == 'B') cbc_u(0,1) = 'D'
        if(cbc_u(1,1) == 'A') cbc_u(1,1) = 'N'
        !-------------- END CUSTOM_BC --------------
        alpha = -.5*visc*dtrk
        if(is_impdiff_1d) then
          call solve_helmholtz(n,ng,hi,alpha=alpha,a=au,b=bu,c=cu,rhsbz=rhsbu%z, &
                               is_bound=is_bound,cbc=cbc_u,c_or_f=['f','c','c'],p=u)
          call solve_helmholtz(n,ng,hi,alpha=alpha,a=av,b=bv,c=cv,rhsbz=rhsbv%z, &
                               is_bound=is_bound,cbc=cbc_v,c_or_f=['c','f','c'],p=v)
          call solve_helmholtz(n,ng,hi,alpha=alpha,a=aw,b=bw,c=cw,rhsbz=rhsbw%z, &
                               is_bound=is_bound,cbc=cbc_w,c_or_f=['c','c','f'],p=w)
        else
          call solve_helmholtz(n,ng,hi,arrplanu,normfftu,alpha, &
                               lambdaxyu,au,bu,cu,rhsbu%x,rhsbu%y,rhsbu%z,is_bound,cbc_u,['f','c','c'],u)
          call solve_helmholtz(n,ng,hi,arrplanv,normfftv,alpha, &
                               lambdaxyv,av,bv,cv,rhsbv%x,rhsbv%y,rhsbv%z,is_bound,cbc_v,['c','f','c'],v)
          call solve_helmholtz(n,ng,hi,arrplanw,normfftw,alpha, &
                               lambdaxyw,aw,bw,cw,rhsbw%x,rhsbw%y,rhsbw%z,is_bound,cbc_w,['c','c','f'],w)
        end if
      end if
      call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,u,v,w,istep,rkcoeff(:,irk),dt,visc,lo,ng,zc)
      call fillps(n,dli,dzfi,dtrki,u,v,w,pp)
      call updt_rhs_b(['c','c','c'],cbcpre,n,is_bound,rhsbp%x,rhsbp%y,rhsbp%z,pp)
      call solver(n,ng,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre,['c','c','c'],pp,is_ptdma_update_p,ap_d,cp_d)
      call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,pp)
      call correc(n,dli,dzci,dtrk,pp,u,v,w)
      call bounduvw(cbcvel,n,bcvel,nb,is_bound,.true.,dl,dzc,dzf,u,v,w,istep,rkcoeff(:,irk),dt,visc,lo,ng,zc)
      call updatep(n,dli,dzci,dzfi,alpha,pp,p)
      call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,p)
    end do
    dpdl(:)     = -dpdl(:)*dti
    fs(1:nscal) = fs(1:nscal)*dti
    !---------------- INLET_PLAN ----------------
    if (is_inlet_record) then
      !if (myid==0) print*, "inlet_record"
      is_inlet_window = (ix0 >= lo(1)) .and. (ix0 <= hi(1))
      !if (is_inlet_window) print*,"myid= ", myid
      itloc = 1 + mod(istep-1, save_freq)
      if (is_inlet_window) then
        ixloc = ix0 - lo(1) + 1
        !$acc wait
        !$acc parallel loop collapse(2) default(present) present(u,v,w,p,u2d,u2dm,v2d,w2d,p2d) async(1)
        do k = 1, n(3)
          do j = 1, n(2)
            u2d(j,k) = u(ixloc,j,k)
            u2dm(j,k) = u(ixloc-1,j,k)
            v2d(j,k) = v(ixloc,j,k)
            w2d(j,k) = w(ixloc,j,k)
            p2d(j,k) = p(ixloc,j,k)
          end do
        end do
        !$acc update self(u2d,u2dm,v2d,w2d,p2d) async(1)
        !$acc wait(1)
        u2d_pack(:,:,itloc) = real(u2d(:,:), real32)
        u2dm_pack(:,:,itloc) = real(u2dm(:,:), real32)
        v2d_pack(:,:,itloc) = real(v2d(:,:), real32)
        w2d_pack(:,:,itloc) = real(w2d(:,:), real32)
        p2d_pack(:,:,itloc) = real(p2d(:,:), real32)
        if (mod(istep, save_freq) == 0) then
                call write_plan_file_mpi(comm_inlet, istep, save_freq, ix0, &
                         ng(2), ng(3), lo(2), lo(3), n(2), n(3), &
                         u2d_pack, u2dm_pack, v2d_pack, w2d_pack, p2d_pack, ios_plan)
        end if
      end if
    end if
    !-------------- END INLET_PLAN --------------
    !
    ! check simulation stopping criteria
    !
    if(stop_type(1)) then ! maximum number of time steps reached
      if(istep >= nstep   ) is_done = is_done.or..true.
    end if
    if(stop_type(2)) then ! maximum simulation time reached
      if(time  >= time_max) is_done = is_done.or..true.
    end if
    if(stop_type(3)) then ! maximum wall-clock time reached
      tw = (MPI_WTIME()-twi)/3600.
      if(tw    >= tw_max  ) is_done = is_done.or..true.
    end if
    if(icheck > 0.and.mod(istep,max(icheck,1)) == 0) then
      if(myid == 0) print*, 'Checking stability and divergence...'
      call chkdt(n,dl,dzci,dzfi,visc,alpha_max,u,v,w,dt_cfl)
      dt = merge(dt_f,min(cfl*dt_cfl,dtmax),dt_f > 0.)
      if(myid == 0) print*, 'dt_cfl = ', dt_cfl, 'dt = ', dt
      if(dt_cfl < small) then
        if(myid == 0) print*, 'ERROR: time step is too small.'
        if(myid == 0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      end if
      dti = 1./dt
      call chkdiv(lo,hi,dli,dzfi,u,v,w,divtot,divmax)
      if(myid == 0) print*, 'Total divergence = ', divtot, '| Maximum divergence = ', divmax
      if(.not.is_mask_divergence_check) then
        if(divmax > small.or.is_nan(divtot)) then
          if(myid == 0) print*, 'ERROR: maximum divergence is too large.'
          if(myid == 0) print*, 'Aborting...'
          is_done = .true.
          kill = .true.
        end if
      end if
    end if
    !
    ! output routines below
    !
    if(iout0d > 0.and.mod(istep,max(iout0d,1)) == 0) then
      !allocate(var(4))
      var(1) = 1.*istep
      var(2) = dt
      var(3) = time
      call out0d(trim(datadir)//'time.out',3,var)
      !
      if(any(is_forced(:)).or.any(abs(bforce(:)) > 0.)) then
        meanvelu = 0.
        meanvelv = 0.
        meanvelw = 0.
        if(is_forced(1).or.abs(bforce(1)) > 0.) then
          call bulk_mean(n,grid_vol_ratio_f,u,meanvelu)
        end if
        if(is_forced(2).or.abs(bforce(2)) > 0.) then
          call bulk_mean(n,grid_vol_ratio_f,v,meanvelv)
        end if
        if(is_forced(3).or.abs(bforce(3)) > 0.) then
          call bulk_mean(n,grid_vol_ratio_c,w,meanvelw)
        end if
        if(.not.any(is_forced(:))) dpdl(:) = -bforce(:) ! constant pressure gradient
        var(1)   = time
        var(2:4) = dpdl(1:3)
        var(5:7) = [meanvelu,meanvelv,meanvelw]
        call out0d(trim(datadir)//'forcing.out',7,var)
      end if
      !
      do iscal=1,nscal
        s => scalars(iscal)
        write(scalnum,'(i3.3)') iscal
        if(s%is_forced.or.abs(s%source) > 0.) then
          meanscal = 0.
          call bulk_mean(n,grid_vol_ratio_f,s%val,meanscal)
          if(.not.s%is_forced) fs(:) = s%source
          var(1:3) = [time,fs(iscal),meanscal]
          call out0d(trim(datadir)//'forcing_s_'//scalnum//'.out',3,var)
        end if
      end do
    end if
    write(fldnum,'(i7.7)') istep
    if(iout1d > 0.and.mod(istep,max(iout1d,1)) == 0) then
      !$acc wait
      !$acc update self(u,v,w,p)
      do iscal=1,nscal
        !$acc update self(scalars(iscal)%val)
      end do
      include 'out1d.h90'
    end if
    if(iout2d > 0.and.mod(istep,max(iout2d,1)) == 0) then
      !$acc wait
      !$acc update self(u,v,w,p)
      do iscal=1,nscal
        !$acc update self(scalars(iscal)%val)
      end do
      include 'out2d.h90'
    end if
    if(iout3d > 0.and.mod(istep,max(iout3d,1)) == 0) then
      if(myid == 0) print*, 'iout3d trigger at istep = ', istep
      !$acc wait
      !$acc update self(u,v,w,p)
      do iscal=1,nscal
        !$acc update self(scalars(iscal)%val)
      end do
      include 'out3d.h90'
    end if
    if(iout3dcrop > 0.and.mod(istep,max(iout3dcrop,1)) == 0 .and. istep>0) then
      !$acc wait
      !$acc update self(u,v,w,p)
      block
        use mpi
        real(rp), pointer, contiguous :: uc(:,:,:)

        include 'out3dcrop.h90'
      end block
    end if
    if(isave > 0.and.((mod(istep,max(isave,1)) == 0).or.(is_done.and..not.kill))) then
      if(is_overwrite_save) then
        filename = 'fld'
      else
        filename = 'fld_'//fldnum
        if(nsaves_max > 0) then
          if(savecounter >= nsaves_max) savecounter = 0
          savecounter = savecounter + 1
          write(chkptnum,'(i4.4)') savecounter
          filename = 'fld_'//chkptnum
          var(1) = 1.*istep
          var(2) = time
          var(3) = 1.*savecounter
          call out0d(trim(datadir)//'log_checkpoints.out',3,var)
        end if
      end if
      !$acc wait
      !$acc update self(u,v,w,p)
      do iscal=1,nscal
        !$acc update self(scalars(iscal)%val)
      end do
      call load_all('w',trim(datadir)//trim(filename)//'.bin', &
                    MPI_COMM_WORLD,ng,[1,1,1],lo,hi,nscal,u,v,w,p,scalars,time,istep)
      if(.not.is_overwrite_save) then
        !
        ! fld_*.bin -> last checkpoint file (symbolic link)
        !
        call gen_alias(myid,trim(datadir),trim(filename)//'.bin','fld.bin')
      end if
      !---------------- TRIPPING ----------------
      if (myid == 0 .and. at /= 0._rp) then
        if(is_overwrite_save) then
          call copy_file(trim(datadir)//'trip.dat', trim(datadir)//'trip_save.dat')
        else
          trip_filename = 'trip_'//fldnum
          if(nsaves_max > 0) trip_filename = 'trip_'//chkptnum
          call copy_file(trim(datadir)//'trip.dat', trim(datadir)//trim(trip_filename)//'.dat')
          call gen_alias(myid,trim(datadir),trim(trip_filename)//'.dat','trip_save.dat')
        end if
      end if
      !-------------- END TRIPPING --------------
      if(myid == 0) print*, '*** Checkpoint saved at time = ', time, 'time step = ', istep, '. ***'
    end if
    if(is_timing) then
      !$acc wait(1)
      dt12 = MPI_WTIME()-dt12
      call MPI_ALLREDUCE(dt12,dt12av ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
      if(myid == 0) print*, 'Avrg, min & max elapsed time: '
      if(myid == 0) print*, dt12av/(1.*product(dims)),dt12min,dt12max
    end if
  end do
  !
  ! clear ffts
  !
  call fftend(arrplanp)
  if(is_impdiff .and. .not.is_impdiff_1d) then
    call fftend(arrplanu)
    call fftend(arrplanv)
    call fftend(arrplanw)
  end if
  if(myid == 0.and.(.not.kill)) print*, '*** Fim ***'
  call decomp_2d_finalize
#if defined(_OPENACC)
  call cudecomp_finalize
#endif
  call MPI_FINALIZE(ierr)
end program cans
