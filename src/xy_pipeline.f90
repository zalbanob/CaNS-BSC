! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT

module mod_xy_pipeline
#if defined(_OPENACC) && !defined(_USE_HIP) && !defined(_USE_DIEZDECOMP)
  use mod_types
  use cudecomp
  use openacc
  !@cuf use cudafor
  use mod_common_cudecomp, only: handle,gd_poi,ap_x_poi,ap_y_poi,work, &
                                 cudecomp_real_rp,istream_acc_queue_1,cudecomp_is_t_in_place
  use mod_fft           , only: fftini_gpu,fftend,set_cufft_wspace, &
                                 signal_processing,fftf_gpu,fftb_gpu,copy
  use mod_utils         , only: f_sizeof
  implicit none
  private
  public xy_pipeline_init,xy_pipeline_set_workspace,xy_pipeline_finalize,xy_pipeline

  integer, parameter :: max_chunks = 4, max_descriptors = 3, transpose_queue = 3
  
  type(cudecompGridDesc) :: descriptors(max_descriptors)
  integer :: nchunks = 0, ndesc = 0
  integer :: width(max_chunks), first(max_chunks), descriptor(max_chunks)
  integer :: chunk_shape(3,2,max_descriptors), plans(2,2,max_descriptors)
  integer :: zpos(2)

  integer(acc_handle_kind)      :: transpose_stream
  type(cudaEvent)               :: ready(max_chunks),received(max_chunks)
  real(rp), allocatable, target :: stage_x(:,:),stage_y(:,:)
  real(rp), pointer, contiguous :: transpose_work(:),transpose_work_cuda(:)
  !@cuf attributes(device)      :: transpose_work_cuda

contains

  subroutine xy_pipeline_init(wsize_fft)
    use mod_param, only: ng
    implicit none
    integer(i8), intent(out)     :: wsize_fft
    type(cudecompGridDescConfig) :: conf
    type(cudecompPencilInfo)     :: pencil
    integer :: c,d, axis, base, rem, nlocal, nranks_z
    integer :: global_base, global_rem, chunk_base, chunk_rem, chunk_z
    integer :: descriptor_z(max_descriptors)
    integer(i8) :: wsize, max_transpose
    integer(i8) :: max_stage(2)

    wsize_fft = max(ap_x_poi%size,ap_y_poi%size)
    call check(cudecompGetGridDescConfig(handle,gd_poi,conf),'get parent descriptor')
    nranks_z    = conf%pdims(2)
    global_base = ng(3) / nranks_z
    global_rem  = mod(ng(3),nranks_z)

    zpos(1) = findloc(ap_x_poi%order,3,dim=1)
    zpos(2) = findloc(ap_y_poi%order,3,dim=1)
    nchunks = 1
    if(conf%pdims(1) > 1 .and. .not. cudecomp_is_t_in_place) nchunks = max(1,min(max_chunks,global_base))

    ! A local XY transpose or aliased pencil storage uses the existing full FFT
    ! plans and descriptor, are ordered on the compute queue
    if(nchunks == 1) return

    wsize_fft = 0
    nlocal = ap_x_poi%shape(zpos(1))
    base   = nlocal/nchunks
    rem    = mod(nlocal,nchunks)

    ! Split each existing Z partition, preserving ownership across X <-> Y.
    ! Compress these local chunks into a synthetic global Z extent. Splitting
    ! ng(3) first would repartition uneven Z slabs differently from the parent.
    chunk_base = global_base/nchunks
    chunk_rem  = mod(global_base,nchunks)
    ndesc = 0
    max_transpose = 1
    max_stage = 0

    do c=1,nchunks
      width(c) = base+merge(1,0,c <= rem)
      first(c) = 1+(c-1)*base+min(c-1,rem)
      chunk_z  = nranks_z*(chunk_base + merge(1 , 0, c <= chunk_rem)) + merge(global_rem, 0, c == chunk_rem+1)
      d = findloc(descriptor_z(:ndesc),chunk_z,dim=1)
      if(d == 0) then
        ndesc = ndesc+1
        d = ndesc
        descriptor_z(d) = chunk_z
        conf%gdims(3) = chunk_z
        conf%gdims_dist(3) = chunk_z

        ! Every rank creates the same descriptors in this global chunk order.
        call check(cudecompGridDescCreate(handle,descriptors(d),conf),'create chunk descriptor')
        call check(cudecompGetTransposeWorkspaceSize(handle,descriptors(d),wsize), 'get transpose workspace')
        max_transpose = max(max_transpose,wsize)

        do axis=1,2
          call check(cudecompGetPencilInfo(handle,descriptors(d),pencil,axis),'get chunk pencil')
          chunk_shape(:,axis,d) = pencil%shape
          max_stage(axis) = max(max_stage(axis),pencil%size)
          call fftini_gpu(ng(axis),pencil%shape,plans(:,axis,d),wsize)
          wsize_fft = max(wsize_fft,pencil%size,wsize)
        end do
      end if

      descriptor(c) = d
      call check(cudaEventCreateWithFlags(ready(c),cudaEventDisableTiming),'create ready event')
      call check(cudaEventCreateWithFlags(received(c),cudaEventDisableTiming),'create receive event')
    end do

    ! Two slots preserve the producer/consumer overlap when a Z chunk is not
    ! contiguous in the parent pencil. Each compact slot keeps its parent order.
    if(zpos(1) /= 3) then
      allocate(stage_x(max_stage(1),2))
      !$acc enter data create(stage_x)
    end if

    if(zpos(2) /= 3) then
      allocate(stage_y(max_stage(2),2))
      !$acc enter data create(stage_y)
    end if

    ! since transposition and FFT can run concurrently they need seperate workspaces
    allocate(transpose_work(max_transpose))
    call check(cudecompMalloc(handle,descriptors(1),transpose_work_cuda,max_transpose), 'allocate transpose workspace')
    call acc_map_data(transpose_work,transpose_work_cuda,max_transpose*f_sizeof(1._rp))
  end subroutine xy_pipeline_init

  subroutine xy_pipeline_set_workspace
    implicit none

    if(nchunks == 1) return
    !$acc wait(1) async(transpose_queue)
    transpose_stream = acc_get_cuda_stream(transpose_queue)

    if(transpose_stream == 0 .or. transpose_stream == istream_acc_queue_1) &
      error stop 'XY pipeline requires a separate transpose stream'
    
    call set_cufft_wspace(pack(plans(:,:,:ndesc),.true.),istream_acc_queue_1)
  end subroutine xy_pipeline_set_workspace

  subroutine xy_pipeline(direction,ng,bc,c_or_f,px,py,arrplan)
    implicit none
    character(len=1), intent(in) :: direction
    integer, intent(in) :: ng(3)
    integer, intent(in) :: arrplan(2,2)
    character(len=1), intent(in) :: bc(0:1,3),c_or_f(3)
    real(rp),  target, contiguous, intent(inout) :: px(:,:,:), py(:,:,:)
    real(rp), pointer, contiguous :: xchunk(:,:,:), ychunk(:,:,:)
    integer :: c, d, producer, consumer
    ! c is the chunk iterator
    ! d is the descriptor for the chunk

    producer = 1
    consumer = 2
    if(direction == 'B') then
      producer = 2
      consumer = 1
    end if

    call transform_chunk(1, producer)
    if(nchunks == 1) then
      !$acc host_data use_device(px,py,work)
      if(direction == 'F') then
        call check(cudecompTransposeXtoY(handle,gd_poi,px,py,work,cudecomp_real_rp,  stream=istream_acc_queue_1),'transpose X to Y')
      else
        call check(cudecompTransposeYtoX(handle,gd_poi,py,px,work,cudecomp_real_rp, stream=istream_acc_queue_1),'transpose Y to X')
      end if
      !$acc end host_data
      call transform_chunk(1,consumer)
      return
    end if

    call check(cudaEventRecord(ready(1), istream_acc_queue_1), 'record first producer')
    do c = 1, nchunks
      d = descriptor(c)
      call check(cudaStreamWaitEvent(transpose_stream,ready(c),0), 'wait for producer')
      if(c < nchunks) then
        call transform_chunk(c+1,producer)
        call check(cudaEventRecord(ready(c+1), istream_acc_queue_1), 'record next producer')
      end if

      call get_chunk(c, 1, px, xchunk)
      call get_chunk(c, 2, py, ychunk)
      !$acc host_data use_device(xchunk, ychunk, transpose_work)
      if(direction == 'F') then
        call check(cudecompTransposeXtoY(handle,descriptors(d), xchunk, ychunk, transpose_work, cudecomp_real_rp,stream=transpose_stream),'transpose X to Y')
      else
        call check(cudecompTransposeYtoX(handle,descriptors(d), ychunk, xchunk, transpose_work, cudecomp_real_rp,stream=transpose_stream),'transpose Y to X')
      end if
      !$acc end host_data

      ! unpack and the self-copy before consuming a complete chunk
      call check(cudaEventRecord(received(c),transpose_stream),'record received chunk')
      call check(cudaStreamWaitEvent(istream_acc_queue_1,received(c),0),'wait for received chunk')
      call transform_chunk(c,consumer)
    end do

  contains

    subroutine get_chunk(c,axis,parent,chunk)
      implicit none
      integer, intent(in) :: c,axis
      real(rp), target, contiguous, intent(inout) :: parent(:,:,:)
      real(rp), pointer, contiguous, intent(out)  :: chunk(:,:,:)
      integer :: k,slot,ns(3)

      if(nchunks == 1) then
        chunk => parent
      else if(zpos(axis) == 3) then
        k = first(c)
        chunk => parent(:,:,k:k+width(c)-1)
      else
        ns = chunk_shape(:,axis,descriptor(c))
        slot = 1 + mod(c - 1, 2)
        if(axis == 1) then
          chunk(1:ns(1),1:ns(2),1:ns(3)) => stage_x(1:product(ns),slot)
        else
          chunk(1:ns(1),1:ns(2),1:ns(3)) => stage_y(1:product(ns),slot)
        end if
      end if
    end subroutine get_chunk

    subroutine transform_chunk(c,axis)
      implicit none
      integer, intent(in) :: c,axis
      real(rp), pointer, contiguous :: parent(:,:,:),chunk(:,:,:)
      real(rp), pointer :: section(:,:,:)
      integer :: ns(3),lo(3),hi(3),idir,plan
      logical :: is_staged

      idir = 1
      if(direction == 'B') idir = 2
      if(axis == 1) then
        parent => px
      else
        parent => py
      end if

      call get_chunk(c,axis,parent,chunk)
      ns = shape(chunk)
      is_staged = nchunks > 1 .and. zpos(axis) /= 3
      if(nchunks == 1) then
        plan = arrplan(idir,axis)
      else
        plan = plans(idir,axis,descriptor(c))
      end if

      if(is_staged) then
        lo = 1
        hi = shape(parent)
        lo(zpos(axis)) = first(c)
        hi(zpos(axis)) = first(c) + width(c) - 1
        section => parent(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
        if(axis == producer) call copy(section,chunk)
      end if

      call signal_processing(0,direction, bc(0,axis)//bc(1,axis), c_or_f(axis), ng(axis), ns, 1, chunk)
      if(direction == 'F') then
        call fftf_gpu(plan,chunk)
      else
        call fftb_gpu(plan,chunk)
      end if
      call signal_processing(1,direction, bc(0,axis)//bc(1,axis), c_or_f(axis), ng(axis), ns, 1, chunk)
      if(is_staged .and. axis == consumer) call copy(chunk,section)
    end subroutine transform_chunk

  end subroutine xy_pipeline

  subroutine xy_pipeline_finalize
    implicit none
    integer :: c,d

    if(nchunks == 1) return

    !$acc wait
    do c = 1,nchunks
      call check(cudaEventDestroy(ready(c)),'destroy ready event')
      call check(cudaEventDestroy(received(c)),'destroy receive event')
    end do
    
    do d=1,ndesc
      call fftend(plans(:,:,d))
    end do

    call acc_unmap_data(transpose_work)
    call check(cudecompFree(handle,descriptors(1),transpose_work_cuda),'free transpose workspace')
    deallocate(transpose_work)

    if(allocated(stage_x)) then
      !$acc exit data delete(stage_x)
      deallocate(stage_x)
    end if

    if(allocated(stage_y)) then
      !$acc exit data delete(stage_y)
      deallocate(stage_y)
    end if

    do d=1,ndesc
      call check(cudecompGridDescDestroy(handle,descriptors(d)),'destroy chunk descriptor')
    end do

    ndesc = 0
    nchunks = 0
  end subroutine xy_pipeline_finalize

  subroutine check(status,operation)
    implicit none
    integer, intent(in) :: status
    character(len=*), intent(in) :: operation

    if(status /= 0) then
      write(*,*) 'XY pipeline: ',operation,', status ',status
      error stop
    end if
  end subroutine check
#endif
end module mod_xy_pipeline
