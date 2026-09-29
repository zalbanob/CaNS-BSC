! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_types
  use mpi, only: MPI_REAL, MPI_DOUBLE_PRECISION, MPI_DATATYPE_NULL, &
                 MPI_TYPECLASS_REAL, MPI_TYPE_MATCH_SIZE
  implicit none

  integer, parameter, public :: sp = selected_real_kind(6 , 37), &
                                dp = selected_real_kind(15,307), &
                                i8 = selected_int_kind(18)

#if defined(_SINGLE_PRECISION)
  integer, parameter, public :: rp = sp
  integer, parameter, public :: MPI_REAL_RP = MPI_REAL
#else
  integer, parameter, public :: rp = dp
  integer, parameter, public :: MPI_REAL_RP = MPI_DOUBLE_PRECISION
#endif

  integer, save, public :: MPI_REAL_SP = MPI_DATATYPE_NULL

contains

  subroutine init_types_mpi(ierr)
    integer, intent(out) :: ierr
    integer :: nbytes
    nbytes = storage_size(1.0_sp)/8
    call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, nbytes, MPI_REAL_SP, ierr)
  end subroutine init_types_mpi

end module mod_types

