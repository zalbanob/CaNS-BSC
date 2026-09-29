module mod_plan
  use, intrinsic :: iso_fortran_env, only: real32,int32
  use mpi
  implicit none
  private
  public :: make_fname,write_plan_file_mpi, read_plan_file_mpi

contains

  subroutine make_fname(start_it, save_freq, ix0, fname)
    integer, intent(in) :: start_it
    integer, intent(in) :: save_freq, ix0
    character(len=*), intent(out) :: fname
    integer :: end_it

    end_it = start_it + save_freq - 1
    write(fname,'("it",I0,"_",I0,"_",I0,".bin")') start_it, end_it, ix0
  end subroutine make_fname

  subroutine write_plan_file(istep_end, save_freq, ix0, ny, nz, &
                             ubuf, umbuf, vbuf, wbuf, pbuf, iostat_out)
    integer, intent(in) :: istep_end
    integer, intent(in) :: save_freq, ix0, ny, nz
    real(real32), intent(in) :: ubuf(ny,nz,save_freq), umbuf(ny,nz,save_freq), &
                                vbuf(ny,nz,save_freq), &
                                wbuf(ny,nz,save_freq), pbuf(ny,nz,save_freq)
    integer, intent(out) :: iostat_out

    character(len=256) :: fname
    integer :: istep_start
    integer :: iunit, ios

    iostat_out = 0
    istep_start = istep_end - save_freq + 1
    call make_fname(istep_start, save_freq, ix0, fname)

    open(newunit=iunit, file=trim(fname), access='stream', form='unformatted', &
         status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      iostat_out = ios
      return
    end if

    if (ios == 0) write(iunit, iostat=ios) ny, nz, save_freq, ix0
    if (ios == 0) write(iunit, iostat=ios) istep_start, istep_end
    if (ios /= 0) then
      iostat_out = ios
      close(iunit)
      return
    end if

    write(iunit, iostat=ios) ubuf
    if (ios == 0) write(iunit, iostat=ios) umbuf
    if (ios == 0) write(iunit, iostat=ios) vbuf
    if (ios == 0) write(iunit, iostat=ios) wbuf
    if (ios == 0) write(iunit, iostat=ios) pbuf
    if (ios /= 0) iostat_out = ios

    close(iunit)
  end subroutine write_plan_file
 
  subroutine write_plan_file_mpi(comm, istep_end, save_freq, ix0, ny_g, nz_g, lo_y, lo_z, ny_l, nz_l, &
                                 ubuf, umbuf, vbuf, wbuf, pbuf, iostat_out)
    integer, intent(in) :: comm
    integer, intent(in) :: istep_end
    integer, intent(in) :: save_freq, ix0, ny_g, nz_g, lo_y, lo_z, ny_l, nz_l
    real(real32), intent(in) :: ubuf(ny_l,nz_l,save_freq), umbuf(ny_l,nz_l,save_freq), &
                                vbuf(ny_l,nz_l,save_freq), &
                                wbuf(ny_l,nz_l,save_freq), pbuf(ny_l,nz_l,save_freq)
    integer, intent(out) :: iostat_out

    character(len=256) :: fname
    integer(int32) :: header(6)
    integer :: istep_start
    integer :: fh, myrank, ierr, ierr_all, bad, bad_all
    integer :: mpi_real32, mpi_int32
    integer(kind=MPI_OFFSET_KIND) :: header_bytes, field_bytes, total_bytes

    iostat_out = 0
    istep_start = istep_end - save_freq + 1
    call make_fname(istep_start, save_freq, ix0, fname)

    bad = 0
    if (save_freq <= 0 .or. ny_g <= 0 .or. nz_g <= 0 .or. ny_l <= 0 .or. nz_l <= 0) bad = 1
    if (lo_y < 1 .or. lo_z < 1) bad = 1
    if (lo_y + ny_l - 1 > ny_g .or. lo_z + nz_l - 1 > nz_g) bad = 1
    call MPI_ALLREDUCE(bad, bad_all, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      iostat_out = ierr
      return
    end if
    if (bad_all /= 0) then
      iostat_out = -3
      return
    end if

    call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, storage_size(1.0_real32)/8, mpi_real32, ierr)
    if (ierr == MPI_SUCCESS) call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, storage_size(header(1))/8, mpi_int32, ierr)
    if (ierr /= MPI_SUCCESS) then
      iostat_out = ierr
      return
    end if

    header_bytes = int(size(header), MPI_OFFSET_KIND)*int(storage_size(header(1))/8, MPI_OFFSET_KIND)
    field_bytes = int(ny_g, MPI_OFFSET_KIND)*int(nz_g, MPI_OFFSET_KIND)* &
                  int(save_freq, MPI_OFFSET_KIND)*int(storage_size(1.0_real32)/8, MPI_OFFSET_KIND)
    total_bytes = header_bytes + 5_MPI_OFFSET_KIND*field_bytes

    call MPI_FILE_OPEN(comm, trim(fname), MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
    if (ierr /= MPI_SUCCESS) then
      iostat_out = ierr
      return
    end if

    call MPI_FILE_SET_SIZE(fh, 0_MPI_OFFSET_KIND, ierr)
    iostat_out = ierr
    call MPI_COMM_RANK(comm, myrank, ierr)
    if (iostat_out == 0 .and. ierr /= MPI_SUCCESS) iostat_out = ierr
    ierr_all = iostat_out
    call MPI_ALLREDUCE(MPI_IN_PLACE, ierr_all, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
    if (ierr_all == 0 .and. ierr /= MPI_SUCCESS) ierr_all = ierr
    iostat_out = ierr_all

    if (iostat_out == 0 .and. myrank == 0) then
      header = [int(ny_g,int32), int(nz_g,int32), int(save_freq,int32), &
                int(ix0,int32), int(istep_start,int32), int(istep_end,int32)]
      call MPI_FILE_WRITE_AT(fh, 0_MPI_OFFSET_KIND, header, size(header), mpi_int32, MPI_STATUS_IGNORE, ierr)
      if (ierr /= MPI_SUCCESS) iostat_out = ierr
    end if
    call MPI_BCAST(iostat_out, 1, MPI_INTEGER, 0, comm, ierr)
    if (iostat_out == 0 .and. ierr /= MPI_SUCCESS) iostat_out = ierr
    call MPI_BARRIER(comm, ierr)
    if (iostat_out == 0 .and. ierr /= MPI_SUCCESS) iostat_out = ierr

    if (iostat_out == 0) call io_plan_field_mpi('w', fh, 0, header_bytes, ny_g, nz_g, save_freq, &
                                                lo_y, lo_z, ny_l, nz_l, ubuf, mpi_real32, iostat_out)
    if (iostat_out == 0) call io_plan_field_mpi('w', fh, 1, header_bytes, ny_g, nz_g, save_freq, &
                                                lo_y, lo_z, ny_l, nz_l, umbuf, mpi_real32, iostat_out)
    if (iostat_out == 0) call io_plan_field_mpi('w', fh, 2, header_bytes, ny_g, nz_g, save_freq, &
                                                lo_y, lo_z, ny_l, nz_l, vbuf, mpi_real32, iostat_out)
    if (iostat_out == 0) call io_plan_field_mpi('w', fh, 3, header_bytes, ny_g, nz_g, save_freq, &
                                                lo_y, lo_z, ny_l, nz_l, wbuf, mpi_real32, iostat_out)
    if (iostat_out == 0) call io_plan_field_mpi('w', fh, 4, header_bytes, ny_g, nz_g, save_freq, &
                                                lo_y, lo_z, ny_l, nz_l, pbuf, mpi_real32, iostat_out)

    ierr_all = iostat_out
    call MPI_ALLREDUCE(MPI_IN_PLACE, ierr_all, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
    if (ierr == MPI_SUCCESS .and. ierr_all == 0) call MPI_FILE_SET_SIZE(fh, total_bytes, ierr)
    if (ierr_all == 0 .and. ierr /= MPI_SUCCESS) ierr_all = ierr
    call MPI_FILE_CLOSE(fh, ierr)
    if (ierr_all == 0 .and. ierr /= MPI_SUCCESS) ierr_all = ierr
    iostat_out = ierr_all
  end subroutine write_plan_file_mpi
  subroutine read_plan_file(filename, save_freq, &
                            ubuf, vbuf, wbuf, pbuf, iostat_out)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: save_freq
    real(real32), intent(inout) :: ubuf(:,:,:), vbuf(:,:,:), wbuf(:,:,:), pbuf(:,:,:)
    integer, intent(out) :: iostat_out

    integer :: ny, nz, nrec, ix0_file, istep_start, istep_end
    integer :: iunit, ios
    integer :: nfields
    integer(kind=8) :: filesize, header_bytes, field_bytes, expected4, expected5
    real(real32), allocatable :: tmp(:,:,:)

    iostat_out = 0

    open(newunit=iunit, file=trim(filename), access='stream', form='unformatted', &
         status='old', action='read', iostat=ios)
    if (ios /= 0) then
      iostat_out = ios
      return
    end if


    if (ios == 0) then
      read(iunit, iostat=ios) ny, nz, nrec, ix0_file
      if (ios == 0) read(iunit, iostat=ios) istep_start, istep_end
      if (ios /= 0) then
        iostat_out = ios
        close(iunit)
        return
      end if
      if (ny /= size(ubuf,1) .or. nz /= size(ubuf,2) .or. nrec /= size(ubuf,3) .or. save_freq /= nrec) then
        iostat_out = -2
        close(iunit)
        return
      end if
    else
      rewind(iunit)
      nrec = size(ubuf,3)
      if (save_freq /= nrec) then
        iostat_out = -2
        close(iunit)
        return
      end if
    end if

    inquire(unit=iunit, size=filesize)
    header_bytes = 6_8*int(storage_size(ix0_file)/8, kind=8)
    field_bytes = int(nrec, kind=8)*int(ny, kind=8)*int(nz, kind=8)* &
                  int(storage_size(1.0_real32)/8, kind=8)
    expected4 = header_bytes + 4_8*field_bytes
    expected5 = header_bytes + 5_8*field_bytes
    if (filesize == expected5) then
      nfields = 5
    else if (filesize == expected4) then
      nfields = 4
    else
      iostat_out = -2
      close(iunit)
      return
    end if

    read(iunit, iostat=ios) ubuf
    if (nfields == 5 .and. ios == 0) then
      allocate(tmp(ny,nz,nrec))
      read(iunit, iostat=ios) tmp
      deallocate(tmp)
    end if
    if (ios == 0) read(iunit, iostat=ios) vbuf
    if (ios == 0) read(iunit, iostat=ios) wbuf
    if (ios == 0) read(iunit, iostat=ios) pbuf
    if (ios /= 0) iostat_out = ios

    close(iunit)
  end subroutine read_plan_file

  subroutine read_plan_file_mpi(filename, comm, save_freq, ny_g, nz_g, lo_y, lo_z, ny_l, nz_l, &
                                ubuf, vbuf, wbuf, pbuf, iostat_out)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: comm
    integer, intent(in) :: save_freq, ny_g, nz_g, lo_y, lo_z, ny_l, nz_l
    real(real32), intent(inout) :: ubuf(ny_l,nz_l,save_freq), vbuf(ny_l,nz_l,save_freq), &
                                   wbuf(ny_l,nz_l,save_freq), pbuf(ny_l,nz_l,save_freq)
    integer, intent(out) :: iostat_out

    integer(int32) :: header(6)
    integer :: fh, ierr, ierr_all, bad, bad_all
    integer :: mpi_real32, mpi_int32
    integer :: iv, iw, ip
    integer(kind=MPI_OFFSET_KIND) :: header_bytes, field_bytes, filesize, expected4, expected5

    iostat_out = 0

    bad = 0
    if (save_freq <= 0 .or. ny_g <= 0 .or. nz_g <= 0 .or. ny_l <= 0 .or. nz_l <= 0) bad = 1
    if (lo_y < 1 .or. lo_z < 1) bad = 1
    if (lo_y + ny_l - 1 > ny_g .or. lo_z + nz_l - 1 > nz_g) bad = 1
    call MPI_ALLREDUCE(bad, bad_all, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      iostat_out = ierr
      return
    end if
    if (bad_all /= 0) then
      iostat_out = -3
      return
    end if

    call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, storage_size(1.0_real32)/8, mpi_real32, ierr)
    if (ierr == MPI_SUCCESS) call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, storage_size(header(1))/8, mpi_int32, ierr)
    if (ierr /= MPI_SUCCESS) then
      iostat_out = ierr
      return
    end if

    header_bytes = int(size(header), MPI_OFFSET_KIND)*int(storage_size(header(1))/8, MPI_OFFSET_KIND)
    field_bytes = int(ny_g, MPI_OFFSET_KIND)*int(nz_g, MPI_OFFSET_KIND)* &
                  int(save_freq, MPI_OFFSET_KIND)*int(storage_size(1.0_real32)/8, MPI_OFFSET_KIND)

    call MPI_FILE_OPEN(comm, trim(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    if (ierr /= MPI_SUCCESS) then
      iostat_out = ierr
      return
    end if

    call MPI_FILE_READ_AT(fh, 0_MPI_OFFSET_KIND, header, size(header), mpi_int32, MPI_STATUS_IGNORE, ierr)
    if (ierr /= MPI_SUCCESS) iostat_out = ierr
    if (iostat_out == 0) then
      if (header(1) /= int(ny_g,int32) .or. header(2) /= int(nz_g,int32) .or. &
          header(3) /= int(save_freq,int32)) iostat_out = -2
    end if
    if (iostat_out == 0) then
      call MPI_FILE_GET_SIZE(fh, filesize, ierr)
      if (ierr /= MPI_SUCCESS) iostat_out = ierr
    end if
    if (iostat_out == 0) then
      expected4 = header_bytes + 4_MPI_OFFSET_KIND*field_bytes
      expected5 = header_bytes + 5_MPI_OFFSET_KIND*field_bytes
      if (filesize == expected5) then
        iv = 2
        iw = 3
        ip = 4
      else if (filesize == expected4) then
        iv = 1
        iw = 2
        ip = 3
      else
        iostat_out = -2
      end if
    end if

    if (iostat_out == 0) call io_plan_field_mpi('r', fh, 0, header_bytes, ny_g, nz_g, save_freq, &
                                                lo_y, lo_z, ny_l, nz_l, ubuf, mpi_real32, iostat_out)
    if (iostat_out == 0) call io_plan_field_mpi('r', fh, iv, header_bytes, ny_g, nz_g, save_freq, &
                                                lo_y, lo_z, ny_l, nz_l, vbuf, mpi_real32, iostat_out)
    if (iostat_out == 0) call io_plan_field_mpi('r', fh, iw, header_bytes, ny_g, nz_g, save_freq, &
                                                lo_y, lo_z, ny_l, nz_l, wbuf, mpi_real32, iostat_out)
    if (iostat_out == 0) call io_plan_field_mpi('r', fh, ip, header_bytes, ny_g, nz_g, save_freq, &
                                                lo_y, lo_z, ny_l, nz_l, pbuf, mpi_real32, iostat_out)

    ierr_all = iostat_out
    call MPI_ALLREDUCE(MPI_IN_PLACE, ierr_all, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
    if (ierr_all == 0 .and. ierr /= MPI_SUCCESS) ierr_all = ierr
    call MPI_FILE_CLOSE(fh, ierr)
    if (ierr_all == 0 .and. ierr /= MPI_SUCCESS) ierr_all = ierr
    iostat_out = ierr_all
  end subroutine read_plan_file_mpi

  subroutine io_plan_field_mpi(io, fh, ifield, header_bytes, ny_g, nz_g, nrec, lo_y, lo_z, ny_l, nz_l, &
                               buf, mpi_real32, iostat_out)
    character(len=1), intent(in) :: io
    integer, intent(in) :: fh, ifield, ny_g, nz_g, nrec, lo_y, lo_z, ny_l, nz_l, mpi_real32
    integer(kind=MPI_OFFSET_KIND), intent(in) :: header_bytes
    real(real32) :: buf(ny_l,nz_l,nrec)
    integer, intent(inout) :: iostat_out

    integer :: ierr, filetype
    integer :: sizes(3), subsizes(3), starts(3)
    integer(kind=MPI_OFFSET_KIND) :: field_bytes, disp

    sizes    = [ny_g, nz_g, nrec]
    subsizes = [ny_l, nz_l, nrec]
    starts   = [lo_y-1, lo_z-1, 0]

    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, mpi_real32, filetype, ierr)
    if (ierr /= MPI_SUCCESS) then
      iostat_out = ierr
      return
    end if
    call MPI_TYPE_COMMIT(filetype, ierr)
    if (ierr /= MPI_SUCCESS) then
      iostat_out = ierr
      call MPI_TYPE_FREE(filetype, ierr)
      return
    end if

    field_bytes = int(ny_g, MPI_OFFSET_KIND)*int(nz_g, MPI_OFFSET_KIND)* &
                  int(nrec, MPI_OFFSET_KIND)*int(storage_size(1.0_real32)/8, MPI_OFFSET_KIND)
    disp = header_bytes + int(ifield, MPI_OFFSET_KIND)*field_bytes

    call MPI_FILE_SET_VIEW(fh, disp, mpi_real32, filetype, 'native', MPI_INFO_NULL, ierr)
    if (ierr == MPI_SUCCESS) then
      select case(io)
      case('r')
        call MPI_FILE_READ_ALL(fh, buf, ny_l*nz_l*nrec, mpi_real32, MPI_STATUS_IGNORE, ierr)
      case('w')
        call MPI_FILE_WRITE_ALL(fh, buf, ny_l*nz_l*nrec, mpi_real32, MPI_STATUS_IGNORE, ierr)
      case default
        ierr = -4
      end select
    end if
    if (ierr /= MPI_SUCCESS) iostat_out = ierr

    call MPI_TYPE_FREE(filetype, ierr)
    if (iostat_out == 0 .and. ierr /= MPI_SUCCESS) iostat_out = ierr
  end subroutine io_plan_field_mpi

end module mod_plan
