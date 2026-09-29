! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
program mini_decomp_out3d
  use mpi
  use decomp_2d, only: decomp_2d_finalize
  use mod_common_mpi, only: myid, ierr
  use mod_initmpi, only: initmpi
  use mod_output, only: out3d, write_visu_3d_crop
  use mod_param, only: ng, dims, cbcpre, read_input, is_simple_writing, imin,imax,jmin,jmax,kmin,kmax,datadir
  use mod_types, only: rp, sp, i8, init_types_mpi
#if defined(_OPENACC)
  use mod_workspaces, only: cudecomp_finalize
#endif
  implicit none

  integer :: lo(3), hi(3), n(3), n_x_fft(3), n_y_fft(3), lo_z(3), hi_z(3), n_z(3)
  integer :: nb(0:1,3)
  logical :: is_bound(0:1,3)
  integer :: i, j, k
  integer :: ig, jg, kg
  integer :: bytes_per_value_write
  integer(i8) :: npts_global, bytes_global
  real(rp) :: t0_write, t1_write, dt_write
  real(rp), allocatable :: p(:,:,:)
    integer :: gmin(3), gmax(3)
  integer :: lmin(3), lmax(3)
  integer :: comm_crop, ierr_crop, myid_crop
  integer :: color
  integer :: ng_crop(3), lo_crop(3), hi_crop(3)
  logical :: intersect,has_data


  call MPI_INIT(ierr)
  call init_types_mpi(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

  call read_input(myid)
  call initmpi(ng, dims, cbcpre, lo, hi, n, n_x_fft, n_y_fft, lo_z, hi_z, n_z, nb, is_bound)

  allocate(p(n(1), n(2), n(3)))
  do k = 1, n(3)
    kg = lo(3) + k - 1
    do j = 1, n(2)
      jg = lo(2) + j - 1
      do i = 1, n(1)
        ig = lo(1) + i - 1
        p(i,j,k) = real(ig, rp) + 1.0e3_rp * real(jg, rp) + 1.0e6_rp * real(kg, rp)
      end do
    end do
  end do

  t0_write = MPI_WTIME()
  call out3d('mini_out3d.bin', [1,1,1], p)
  t1_write = MPI_WTIME()
  dt_write = t1_write - t0_write
  if (is_simple_writing) then
    bytes_per_value_write = storage_size(1.0_sp) / 8
  else
    bytes_per_value_write = storage_size(1.0_rp) / 8
  end if
  npts_global = int(ng(1), i8) * int(ng(2), i8) * int(ng(3), i8)
  bytes_global = npts_global * int(bytes_per_value_write, i8)

  if (myid == 0) then
    print '(A,3(I0,A))', 'mini_out3d global size (ng): ', ng(1), ' x ', ng(2), ' x ', ng(3)
    if (is_simple_writing) then
      print '(A)', 'mini_out3d write precision: float32 (simple)'
    else
      print '(A,I0,A)', 'mini_out3d write precision: ', storage_size(1.0_rp)/8, ' bytes/value'
    end if
    print '(A,I0,A)', 'mini_out3d expected file size: ', bytes_global, ' bytes'
    print '(A,F12.6,A)', 'mini_out3d write time: ', dt_write, ' s'
  end if

t0_write = MPI_WTIME()

  intersect = ( hi(1) >= imin .and. lo(1) <= imax .and. &
              hi(2) >= jmin .and. lo(2) <= jmax .and. &
              hi(3) >= kmin .and. lo(3) <= kmax )

if (intersect) then
  gmin = [ max(lo(1), imin), max(lo(2), jmin), max(lo(3), kmin) ]
  gmax = [ min(hi(1), imax), min(hi(2), jmax), min(hi(3), kmax) ]
  has_data = all(gmin <= gmax)
else
  has_data = .false.
end if

color = merge(1, MPI_UNDEFINED, has_data)
call MPI_Comm_split(MPI_COMM_WORLD, color, myid, comm_crop, ierr)

if (has_data) then
  lmin = gmin - lo + 1
  lmax = gmax - lo + 1

  ng_crop = [ imax-imin+1, jmax-jmin+1, kmax-kmin+1 ]
  lo_crop = gmin - [imin,jmin,kmin] + 1
  hi_crop = gmax - [imin,jmin,kmin] + 1
  
  call write_visu_3d_crop(comm_crop, datadir,  'vex_fld_crop_.bin', 'log_visu_3d_crop.out', &
                          [imin,jmin,kmin], [imax,jmax,kmax], 'Velocity_X', &
                          ng_crop, lo_crop, hi_crop, [1,1,1], 0._rp, 0, &
                          p(lmin(1):lmax(1), lmin(2):lmax(2), lmin(3):lmax(3)) )

  call MPI_Comm_free(comm_crop, ierr)
end if
  t1_write = MPI_WTIME()
  dt_write = t1_write - t0_write


  if (myid==0) then
    print '(A,F12.6,A)', 'crop write time: ', dt_write, ' s'
  end if
  deallocate(p)
  call decomp_2d_finalize
#if defined(_OPENACC)
  call cudecomp_finalize()
#endif
  call MPI_FINALIZE(ierr)
end program mini_decomp_out3d

