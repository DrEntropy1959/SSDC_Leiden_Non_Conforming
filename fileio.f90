module fileio
  use precision_vars
  implicit none

  include 'cgnslib_f.h'

  private

  public init_cgns_bc
  public cgnsReadUnstructuredGrid
  public cgnsPrelimUnstructuredGrid
  public cgnsParallelUnstructuredGrid
  public readSSDCstartup
  public init_code_bc_types
  public aflr3ReadUnstructuredGrid
  
contains

  subroutine init_cgns_bc()
    use variables, only: bctypes

    !   00 : BCTypeNull
    !   01 : BCTypeUserDefined
    !   02 : BCAxisymmetricWedge
    !   03 : BCDegenerateLine
    !   04 : BCDegeneratePoint
    !   05 : BCDirichlet
    !   06 : BCExtrapolate
    !   07 : BCFarfield
    !   08 : BCGeneral
    !   09 : BCInflow
    !   10 : BCInflowSubsonic
    !   11 : BCInflowSupersonic
    !   12 : BCNeumann
    !   13 : BCOutflow
    !   14 : BCOutflowSubsonic
    !   15 : BCOutflowSupersonic
    !   16 : BCSymmetryPlane
    !   17 : BCSymmetryPolar
    !   18 : BCTunnelInflow
    !   19 : BCTunnelOutflow
    !   20 : BCWall
    !   21 : BCWallInviscid
    !   22 : BCWallViscous
    !   23 : BCWallViscousHeatFlux
    !   24 : BCWallViscousIsothermal
    !   25 : FamilySpecified

    allocate(bctypes(-1:13))
    bctypes = Null
    bctypes(1) = BCGeneral
    bctypes(2) = UserDefined
    bctypes(3) = BCWallInviscid
    bctypes(4) = BCWallViscousIsothermal
    bctypes(5) = BCWallViscous
    bctypes(6) = BCOutflowSubsonic
    bctypes(7) = BCOutflowSupersonic
    bctypes(8) = BCInflowSubsonic
    bctypes(9) = BCInflowSupersonic
    ! bctypes(10) is reserved for periodic in the solver
    ! so FamilySpecified is only used for preprocessing
    bctypes(10) = FamilySpecified 
    bctypes(11) = BCFarfield
    bctypes(12) = BCDirichlet
    bctypes(13) = BCSymmetryPlane

  end subroutine init_cgns_bc

  subroutine error_check(iunit,ierr)
    integer,intent(in) :: ierr, iunit
    character(100) :: errmsg

    if(ierr /= 0) then
      call cg_get_error_f(errmsg)
      call cg_error_exit_f()
      !   call cg_error_print_f()
    end if

  end subroutine error_check

  subroutine cgnsReadUnstructuredGrid(filein)
    use referencevariables
    use variables, only: vx, bctypes, boundaryelems
    implicit none

    character(120), intent(in) :: filein

    integer :: iunit
    ! isize contains number of nodes and cells in each direction
    integer :: isize(3)
    ! indices to keep track of file, base, and zone number
    integer :: ibase, izone
    ! error flag
    integer :: ierr
    ! number of bases and zones
    integer :: nbases, nzones
    ! names of base and zone
    character(60) :: basename, zonename
    integer :: rangemin(3), rangemax(3)
    integer :: i
    ! number of boundary conditions
    integer :: nbocos
    ! boundary condition index
    integer :: iboco
    ! name of connections, donor, and boundary condition
    character(60) :: boconame, familyname
    character(60) :: gridcoordname, elementsectionname
    integer :: ngrids, igrid
    integer :: nsections, isection
    integer :: izonetype
    integer :: eltype, parentflag
    integer :: nbvpe
    ! boundary condition parameters (see cgns mid-level library doc)
    integer :: bocotype, ptsettype, boconpts, bocoNI, bocoNLF, bocoNDT, boconds
    integer, allocatable :: bocoNList(:), NList(:)
    integer :: nfam, ifam, nfambc, nfamgeo

    ! ensure that casefile exists
    call get_unit(iunit)
    open(unit=iunit,file=trim(filein)//'.cgns',status='old',iostat=ierr)
    if(ierr /= 0)then
      write(*,*) "cgnsGridRead:", trim(filein)," doesn't exist."
      stop
    end if
    close(iunit)

    ! Open CGNS File using cgns library
    call cg_open_f(trim(filein)//'.cgns',MODE_READ,iunit,ierr)
    call error_check(iunit,ierr)

    ! read the number of bases (should be 1)
    call cg_nbases_f(iunit,nbases,ierr)
    call error_check(iunit,ierr)
    if (nbases /= 1) then
      write(*,*) 'Error: More than 1 base in cgns file', trim(filein)
      write(*,*) 'number of bases in file: ', nbases
      stop
    end if
    ibase = 1
    ! read basename
    call cg_base_read_f(iunit, ibase, basename, ndim, nphysdim, ierr)
    call error_check(iunit,ierr)

    ! read number of zones (limited to 1 for now)
    call cg_nzones_f(iunit,ibase,nzones,ierr)
    call error_check(iunit,ierr)
    if (nzones > 1) then
      write(*,*) 'Error: currently only support 1 unstructured zone'
      write(*,*) 'number of zones in file: ', nzones
      stop
    end if
    izone = 1
    call cg_zone_type_f(iunit, ibase, izone, izonetype, ierr)
    call error_check(iunit, ierr)
    if (izonetype /= Unstructured) then
      write(*,*) 'Error: not an unstructured zone', izone
      stop
    end if
    call cg_zone_read_f(iunit, ibase, izone, zonename, isize, ierr)
    call error_check(iunit, ierr)

    ! set number of vertices
    nvertices = isize(1)
    ! allocate vertices
    allocate(vx(1:3,1:nvertices))

    ! read vertices
    call cg_ngrids_f(iunit, ibase, izone, ngrids, ierr)
    call error_check(iunit, ierr)
    if (ngrids > 1) then
      write(*,*) 'Error: only support one grid per zone'
      write(*,*) 'ngrids in zone: ', ngrids
      stop
    end if
    igrid = 1
    call cg_grid_read_f(iunit, ibase, izone, igrid, gridcoordname, ierr)
    call error_check(iunit, ierr)
    rangemin = 0
    rangemin(1) = 1
    rangemax = 0
    rangemax(1) = nvertices
    call cg_coord_read_f(iunit, ibase, izone, 'CoordinateX', RealDouble, &
      & rangemin, rangemax, vx(1,:), ierr)
    call error_check(iunit, ierr)
    call cg_coord_read_f(iunit, ibase, izone, 'CoordinateY', RealDouble, &
      & rangemin, rangemax, vx(2,:), ierr)
    call error_check(iunit, ierr)
    call cg_coord_read_f(iunit, ibase, izone, 'CoordinateZ', RealDouble, &
      & rangemin, rangemax, vx(3,:), ierr)
    call error_check(iunit, ierr)

    ! set number of elements
    nelems = isize(2)

    ! read elements
    call cg_nsections_f(iunit, ibase, izone, nsections, ierr)
    call error_check(iunit, ierr)
    nelemzones = nsections
    allocate( boundaryelems(nelemzones) )
    nvolumesections = 1
    allocate( isectionvolume(nvolumesections) )
    isectionvolume = 0
    nverticesperelem = 0
    do isection = 1, nsections
      call cg_section_read_f(iunit, ibase, izone, isection, elementsectionname, eltype, &
        & rangemin(1), rangemax(1), rangemin(2), parentflag, ierr)
      call error_check(iunit, ierr)
      boundaryelems(isection)%name = elementsectionname
      boundaryelems(isection)%eltype = eltype
      call cg_npe_f(eltype,nbvpe,ierr)
      call error_check(iunit, ierr)
      boundaryelems(isection)%nbvpe = nbvpe
      ! allocate boundary elements
      allocate( boundaryelems(isection)%belems(nbvpe,rangemin(1):rangemax(1)) )
      boundaryelems(isection)%ielstart = rangemin(1)
      boundaryelems(isection)%ielend = rangemax(1)
      boundaryelems(isection)%nbelems = rangemax(1)-rangemin(1) + 1
      call cg_elements_read_f(iunit, ibase, izone, isection, boundaryelems(isection)%belems, Null, ierr)
      call error_check(iunit, ierr)
      if (nbvpe > nverticesperelem) then
        nverticesperelem = nbvpe
        isectionvolume = isection
      end if
    end do

    ! read the number of boundary conditions
    call cg_nbocos_f(iunit,ibase,izone,nbocos,ierr)
    write(*,*) izone, 'number of boundary conditions: ',nbocos
    call error_check(iunit,ierr)
    bcloop:do iboco = 1,nbocos
      ! read boundary condition (see cgns documentation)
      call cg_boco_info_f(iunit,ibase,izone,iboco,boconame,bocotype,ptsettype, &
        & boconpts, bocoNI, bocoNLF, bocoNDT, boconds, ierr)
      write(*,*) iboco, 'bocotype: ', bocotype
      call error_check(iunit,ierr)
      allocate(bocoNList(boconpts*3))
      allocate(NList(boconpts))
      call cg_boco_read_f(iunit,ibase,izone,iboco,nlist(1:boconpts),bocoNList,ierr)
      call error_check(iunit,ierr)

      ! determine what type of bc
      if (bocotype==FamilySpecified) then
        call cg_goto_f(iunit, ibase, ierr, 'Zone_t',izone,'ZoneBC_t',1,"BC_t",iboco,'end')
        call error_check(iunit,ierr)
        call cg_famname_read_f(boconame, ierr)
        call error_check(iunit,ierr)
        write(*,*) 'boconame = ', boconame
        ! for FamilySpecified, we must read the family bc
        call cg_nfamilies_f(iunit, ibase, nfam, ierr)
        call error_check(iunit,ierr)
        do ifam = 1,nfam
          call cg_family_read_f(iunit,ibase,ifam,familyname,nfambc,nfamgeo,ierr)
          call error_check(iunit,ierr)
          if (familyname==boconame) exit
        end do
        call cg_fambc_read_f(iunit,ibase,ifam,1,boconame,bocotype,ierr)
        call error_check(iunit,ierr)
      end if
      do isection = 1,nsections
        if (boundaryelems(isection)%name == familyname) exit
      end do
      do i=3,size(bctypes(1:))
        if(bctypes(i)==bocotype)then
          boundaryelems(isection)%btype = i
          exit
        end if
      end do
      deallocate(bocoNList)
      deallocate(Nlist)
    end do bcloop

    call cg_close_f(iunit,ierr)
    call error_check(iunit,ierr)

  end subroutine cgnsReadUnstructuredGrid

  subroutine cgnsPrelimUnstructuredGrid(filein)
    !  CGNS Element types
    !  ElementType_t ElementTypeNull, ElementTypeUserDefined, NODE, BAR_2, BAR_3, TRI_3, TRI_6, QUAD_4, QUAD_8, QUAD_9, 
    !                                                         TETRA_4, TETRA_10, PYRA_5, PYRA_14, PENTA_6, PENTA_15, 
    !                                                         PENTA_18, HEXA_8, HEXA_20, HEXA_27, MIXED, PYRA_13, NGON_n, 
    !                                                         NFACE_n, BAR_4, TRI_9, TRI_10, QUAD_12, QUAD_16, TETRA_16, 
    !                                                         TETRA_20, PYRA_21, PYRA_29, PYRA_30, PENTA_24, PENTA_38, 
    !                                                         PENTA_40, HEXA_32, HEXA_56, HEXA_64
    !  CGNS BC Types
    !  BCType_t     BCTypeNull, BCTypeUserDefined, BCAxisymmetricWedge, BCDegenerateLine, BCExtrapolate, BCDegeneratePoint, 
    !               BCDirichlet, BCFarfield, BCNeumann, BCGeneral, BCInflow, BCOutflow, BCInflowSubsonic, BCOutflowSubsonic, 
    !               BCInflowSupersonic, BCOutflowSupersonic, BCSymmetryPlane, BCTunnelInflow, BCSymmetryPolar, BCTunnelOutflow, 
    !               BCWallViscous, BCWall, BCWallViscousHeatFlux, BCWallInviscid, BCWallViscousIsothermal, FamilySpecified
    !  CGNS BCData Types
    !  BCDataType_t   BCDataTypeNull, BCDataTypeUserDefined, Dirichlet, Neumann


    use referencevariables
    use variables, only: boundaryelems, code_bc_names
    implicit none

    character(120), intent(in) :: filein

    integer :: iunit
    ! isize contains number of nodes and cells in each direction
    integer :: isize(3)
    ! indices to keep track of file, base, and zone number
    integer :: ibase, izone
    ! error flag
    integer :: ierr
    ! number of bases and zones
    integer :: nbases, nzones
    ! names of base and zone
    character(60) :: basename, zonename
    integer :: rangemin(3), rangemax(3)
    integer :: i
    ! number of boundary conditions
    integer :: nbocos
    ! boundary condition index
    integer :: iboco
    ! name of connections, donor, and boundary condition
    character(60) :: boconame, familyname
    character(60) :: elementsectionname
    integer :: nsections, isection
    integer :: izonetype
    integer :: eltype, parentflag
    integer :: nbvpe
    ! boundary condition parameters (see cgns mid-level library doc)
    integer :: bocotype, ptsettype, boconpts, bocoNI, bocoNLF, bocoNDT, boconds
    integer, allocatable :: bocoNList(:), NList(:)
    integer :: nfam, ifam, nfambc, nfamgeo
    integer :: jvcount, jsection, ksection
    character(60) :: cgns_bc_descriptor_to_code_bc

    ! ensure that casefile exists
    call get_unit(iunit)
    open(unit=iunit,file=trim(filein)//'.cgns',status='old',iostat=ierr)
    if(ierr /= 0)then
      write(*,*) "cgnsGridRead:", trim(filein)," doesn't exist."
      stop
    end if
    close(iunit)

    ! Open CGNS File using cgns library
    call cg_open_f(trim(filein)//'.cgns',MODE_READ,iunit,ierr)
    call error_check(iunit,ierr)

    ! read the number of bases (should be 1)
    call cg_nbases_f(iunit,nbases,ierr)
    call error_check(iunit,ierr)
    if (nbases /= 1) then
      write(*,*) 'Error: More than 1 base in cgns file', trim(filein)
      write(*,*) 'number of bases in file: ', nbases
      stop
    end if
    ibase = 1

    ! read basename
    call cg_base_read_f(iunit, ibase, basename, ndim, nphysdim, ierr)
    call error_check(iunit,ierr)

    nvertices = 0
    nverticesperelem = 2**ndim
    nelems = 0
    nvolumesections = 0

    ! read number of zone
    call cg_nzones_f(iunit,ibase,nzones,ierr)
    call error_check(iunit,ierr)
    write(*,*) 'number of zones',nzones

!   allocate(znames(nzones))
!   nblks=nzones ! nblks is a global variable
    ! allocate number of blocks for each zone
!   call allocate_blocks(1,nblks)
!   allocate(blk(nblks))
!   ! fill in the names for each zone (must be done before reading connectivity)
!   do izone=1,nzones
!     call cg_zone_read_f(iunit, ibase, izone, zonename, isize, ierr)
!     call error_check(iunit ,ierr)
!     znames(izone) = zonename ! will correspond to donor name in connectivity
!     grdblk(izone)%nbid=-1 ! -1 indicates no neighbors
!     grdblk(izone)%gblkid = izone ! global block id
!     grdblk(izone)%btype = 0 ! 0 indicates no boundary condition
!     ! set up transform in default directions
!     do i = 1,3
!       grdblk(izone)%itransform(i,:,:) = i
!     end do
!   end do

    do izone = 1, nzones
      call cg_zone_type_f(iunit, ibase, izone, izonetype, ierr)
      call error_check(iunit, ierr)
      if (izonetype /= Unstructured) then
        write(*,*) 'Error: not an unstructured zone', izone
        write(*,*) 'convert all zones to unstructured and try again'
        stop
      end if
      call cg_zone_read_f(iunit, ibase, izone, zonename, isize, ierr)
      call error_check(iunit, ierr)

      ! set number of vertices
      nvertices = nvertices+isize(1)
      ! set number of elements
      nelems = nelems+isize(2)

      ! read elements
      call cg_nsections_f(iunit, ibase, izone, nsections, ierr)
      call error_check(iunit, ierr)
      nelemzones = nelemzones+nsections
      do isection = 1, nsections
        call cg_section_read_f(iunit, ibase, izone, isection, elementsectionname, eltype, &
          & rangemin(1), rangemax(1), rangemin(2), parentflag, ierr)
        call error_check(iunit, ierr)
        call cg_npe_f(eltype,nbvpe,ierr)    !   get number=nbvpe of elements of type=eltype
        call error_check(iunit, ierr)
        if (nbvpe == nverticesperelem) nvolumesections = nvolumesections+1
      end do
    end do

    !   write(*,*)'nbvpe    '      ,nbvpe
    !   write(*,*)'eltype'         ,eltype

    !   write(*,*)'nelemzones'     ,nelemzones
    !   write(*,*)'nvertices'      ,nvertices
    !   write(*,*)'nelems'         ,nelems
    !   write(*,*)'nvolumesections',nvolumesections

    allocate( boundaryelems(nelemzones) )
    allocate( isectionvolume(nvolumesections) )
    isectionvolume = 0
    jvcount = 0
    isection = 0
    do izone = 1,nzones
      ksection = isection
      call cg_nsections_f(iunit, ibase, izone, nsections, ierr)
      call error_check(iunit, ierr)
      do jsection = 1, nsections
        isection = isection+1
        call cg_section_read_f(iunit, ibase, izone, jsection, elementsectionname, eltype, &
          & rangemin(1), rangemax(1), rangemin(2), parentflag, ierr)
        call error_check(iunit, ierr)
        boundaryelems(isection)%name = elementsectionname
        boundaryelems(isection)%eltype = eltype
        call cg_npe_f(eltype,nbvpe,ierr)
        call error_check(iunit, ierr)
        boundaryelems(isection)%nbvpe = nbvpe
        ! allocate boundary elements
        !       write(*,*)'going through, section ',isection,elementsectionname, nbvpe, rangemin(1),rangemax(1)
        allocate( boundaryelems(isection)%belems(nbvpe,rangemin(1):rangemax(1)) )
        boundaryelems(isection)%ielstart = rangemin(1)
        boundaryelems(isection)%ielend = rangemax(1)
        boundaryelems(isection)%nbelems = rangemax(1)-rangemin(1) + 1
        call cg_elements_read_f(iunit, ibase, izone, jsection, boundaryelems(isection)%belems, Null, ierr)
        call error_check(iunit, ierr)
        if (nbvpe == nverticesperelem) then
          jvcount = jvcount + 1
          isectionvolume(jvcount) = isection
        end if
      end do
      ! read the number of boundary conditions
      call cg_nbocos_f(iunit,ibase,izone,nbocos,ierr)
      write(*,*) izone, 'number of boundary conditions: ',nbocos
      call error_check(iunit,ierr)
      bcloop:do iboco = 1,nbocos
        ! read boundary condition (see cgns documentation)
        call cg_boco_info_f(iunit,ibase,izone,iboco,boconame,bocotype,ptsettype, &
          & boconpts, bocoNI, bocoNLF, bocoNDT, boconds, ierr)
        write(*,*) iboco, 'bocotype: ', bocotype
        call error_check(iunit,ierr)
        allocate(bocoNList(boconpts*3))
        allocate(NList(boconpts))
        call cg_boco_read_f(iunit,ibase,izone,iboco,nlist(1:boconpts),bocoNList,ierr)
        call error_check(iunit,ierr)

        ! determine what type of bc
        if (bocotype==FamilySpecified) then
          call cg_goto_f(iunit, ibase, ierr, 'Zone_t',izone,'ZoneBC_t',1,"BC_t",iboco,'end')
          call error_check(iunit,ierr)
          call cg_famname_read_f(boconame, ierr)
          call error_check(iunit,ierr)
          write(*,*) 'boconame = ', boconame
          cgns_bc_descriptor_to_code_bc = boconame
          ! for FamilySpecified, we must read the family bc
          call cg_nfamilies_f(iunit, ibase, nfam, ierr)
          call error_check(iunit,ierr)
          do ifam = 1,nfam
            call cg_family_read_f(iunit,ibase,ifam,familyname,nfambc,nfamgeo,ierr)
            call error_check(iunit,ierr)
            if (familyname==boconame) exit
          end do
          call cg_fambc_read_f(iunit,ibase,ifam,1,boconame,bocotype,ierr)
          call error_check(iunit,ierr)
        end if
        isection = ksection
        do jsection = 1,nsections
          isection = isection+1
          if (boundaryelems(isection)%name == familyname) exit
        end do
!        do i = 3,size(bctypes(1:))
!          if(bctypes(i)==bocotype)then
!            boundaryelems(isection)%btype = i
!            exit
!          end if
!        end do
        do i = 1,size(code_bc_names(1:))
          !write(*,*) cgns_bc_descriptor_to_code_bc
          if(cgns_bc_descriptor_to_code_bc == code_bc_names(i)) then
            boundaryelems(isection)%btype = i
            !write(*,*) boundaryelems(isection)%btype
            exit
          end if
        end do

        deallocate(bocoNList)
        deallocate(Nlist)
      end do bcloop
    end do

    call cg_close_f(iunit,ierr)
    call error_check(iunit,ierr)

  end subroutine cgnsPrelimUnstructuredGrid

  subroutine cgnsParallelUnstructuredGrid(filein)
    use referencevariables
    use variables, only: vx, e2v, boundaryelems, &
      jelems
    implicit none

    character(120), intent(in) :: filein

    integer :: iunit
    ! isize contains number of nodes and cells in each direction
    integer :: isize(3)
    ! indices to keep track of file, base, and zone number
    integer :: ibase, izone
    ! error flag
    integer :: ierr
    ! number of bases and zones
    integer :: nbases, nzones
    ! names of base and zone
    character(60) :: basename, zonename
    integer :: rangemin(3), rangemax(3)
    integer :: i,j,k
    ! name of connections, donor, and boundary condition
    character(60) :: gridcoordname, elementsectionname
    integer :: ngrids, igrid
    integer :: nsections, isection
    integer :: izonetype
    integer :: eltype, parentflag
    integer :: nbvpe
    integer :: jvcount, jsection, ksection
    integer :: ielem
    integer, allocatable :: vlist(:), vlisttmp(:)

    ! ensure that casefile exists
    call get_unit(iunit)
    open(unit=iunit,file=trim(filein)//'.cgns',status='old',iostat=ierr)
    if(ierr /= 0)then
      write(*,*) "cgnsGridRead:", trim(filein)," doesn't exist."
      stop
    end if
    close(iunit)

    ! Open CGNS File using cgns library
    call cg_open_f(trim(filein)//'.cgns',MODE_READ,iunit,ierr)
    call error_check(iunit,ierr)

    ! read the number of bases (should be 1)
    call cg_nbases_f(iunit,nbases,ierr)
    call error_check(iunit,ierr)
    if (nbases /= 1) then
      write(*,*) 'Error: More than 1 base in cgns file', trim(filein)
      write(*,*) 'number of bases in file: ', nbases
      stop
    end if
    ibase = 1
    ! read basename
    call cg_base_read_f(iunit, ibase, basename, ndim, nphysdim, ierr)
    call error_check(iunit,ierr)

    nvertices = 0
    nverticesperelem = 2**ndim
    nelems = 0
    nvolumesections = 0
    ! read number of zones (limited to 1 for now)
    call cg_nzones_f(iunit,ibase,nzones,ierr)
    call error_check(iunit,ierr)
    do izone = 1, nzones
      call cg_zone_type_f(iunit, ibase, izone, izonetype, ierr)
      call error_check(iunit, ierr)
      if (izonetype /= Unstructured) then
        write(*,*) 'Error: not an unstructured zone', izone
        write(*,*) 'convert all zones to unstructured and try again'
        stop
      end if
      call cg_zone_read_f(iunit, ibase, izone, zonename, isize, ierr)
      call error_check(iunit, ierr)

      ! set number of vertices
      nvertices = nvertices+isize(1)
      ! set number of elements
      nelems = nelems+isize(2)

      rangemin(3) = minval(jelems)
      rangemin(3) = maxval(jelems)

      ! read elements
      call cg_nsections_f(iunit, ibase, izone, nsections, ierr)
      call error_check(iunit, ierr)
      do isection = 1, nsections
        call cg_section_read_f(iunit, ibase, izone, isection, elementsectionname, eltype, &
          & rangemin(1), rangemax(1), rangemin(2), parentflag, ierr)
        call error_check(iunit, ierr)
        if(rangemax(1) < rangemin(3) .or. rangemin(1) > rangemax(3))cycle
        nelemzones = nelemzones + 1
        call cg_npe_f(eltype,nbvpe,ierr)     
        call error_check(iunit, ierr)
        if (nbvpe == nverticesperelem) nvolumesections = nvolumesections+1
      end do
    end do
    allocate( boundaryelems(nelemzones) )
    allocate( isectionvolume(nvolumesections) )
    isectionvolume = 0
    jvcount = 0
    isection = 0
    do izone = 1,nzones
      ksection = isection
      call cg_nsections_f(iunit, ibase, izone, nsections, ierr)
      call error_check(iunit, ierr)
      rangemin(3) = minval(jelems)
      rangemin(3) = maxval(jelems)
      do jsection = 1, nsections
        call cg_section_read_f(iunit, ibase, izone, jsection, elementsectionname, eltype, &
          & rangemin(1), rangemax(1), rangemin(2), parentflag, ierr)
        call error_check(iunit, ierr)
        if(rangemax(1) < rangemin(3) .or. rangemin(1) > rangemax(3))cycle
        do ielem = ihelems(1), ihelems(2)
          call cg_elements_partial_read_f(iunit, ibase, izone, jsection, &
            jelems(ielem), jelems(ielem), e2v(:,ielem), Null, ierr)
        end do
        call cg_npe_f(eltype,nbvpe,ierr)
        call error_check(iunit, ierr)
        if (nbvpe == nverticesperelem) then
          jvcount = jvcount + 1
          isectionvolume(jvcount) = isection
        end if
      end do
    end do
    ! find the number of vertices
    allocate(vlisttmp(nvertices))
    vlisttmp = nvertices+100
    nvertices = 0
    jvcount = 1
    do ielem = ihelems(1), ihelems(2)
      vloop: do j = 1, nverticesperelem
        i = e2v(j,ielem)
        ! loop over known vertices
        do k = 1,nvertices
          if (i == vlisttmp(k)) then
            ! assign new local e2v entry
            e2v(j,ielem) = k
            cycle vloop
          end if
        end do
        ! add to vertex list
        vlisttmp(jvcount) = i
        ! assign new local e2v
        e2v(j,ielem) = jvcount
        jvcount = jvcount + 1
        nvertices = nvertices+1
      end do vloop
    end do

    allocate(vlist(nvertices))
    vlist(1:nvertices) = vlisttmp(1:nvertices)
    deallocate(vlisttmp)

    ! allocate vertices
    allocate(vx(1:3,1:nvertices))

    ! read vertices
    do izone = 1,nzones
      call cg_ngrids_f(iunit, ibase, izone, ngrids, ierr)
      call error_check(iunit, ierr)
      if (ngrids > 1) then
        write(*,*) 'Error: only support one grid per zone'
        write(*,*) 'ngrids in zone: ', ngrids
        stop
      end if
      igrid = 1
      call cg_grid_read_f(iunit, ibase, izone, igrid, gridcoordname, ierr)
      call error_check(iunit, ierr)
      rangemin = 0
      rangemin(1) = 1
      rangemax = 0
      rangemax(1) = nvertices
      do i = 1, nvertices
        rangemin = 0
        rangemin(1) = vlist(i)
        rangemax = 0
        rangemax(1) = vlist(i)
        call cg_coord_read_f(iunit, ibase, izone, 'CoordinateX', RealDouble, &
          & rangemin, rangemax, vx(1,i:i), ierr)
        call error_check(iunit, ierr)
        call cg_coord_read_f(iunit, ibase, izone, 'CoordinateY', RealDouble, &
          & rangemin, rangemax, vx(2,i:i), ierr)
        call error_check(iunit, ierr)
        call cg_coord_read_f(iunit, ibase, izone, 'CoordinateZ', RealDouble, &
          & rangemin, rangemax, vx(3,i:i), ierr)
        call error_check(iunit, ierr)
      end do
    end do

    nelems = ihelems(2) - ihelems(1) + 1

    deallocate(vlist)

    call cg_close_f(iunit,ierr)
    call error_check(iunit,ierr)

  end subroutine cgnsParallelUnstructuredGrid

  !============================================================================
  
  subroutine readSSDCstartup(file_name)

    ! Load module
    use referencevariables
    use controlvariables
    use nsereferencevariables
    use collocationvariables
    use SSWENOvariables

    ! Nothing is implicitly defined
    implicit none

    character(len=120), intent(in) :: file_name

    integer :: iunit, ierr


    namelist /PreProcParameters/ ndim, nequations, nprocs, npoly, npoly_max, npoly_DeltaF, &
      &casefile, grid_format, viscous, crossterms, RK_Method, IMEX_element, &
      &IMEX_penalty, physics, discretization, filter_solution, l01, l10, l00, &
      &periodic, periodic_distance, Entropy_Correction, variable_viscosity, &
      &heat_entropy_flow_wall_bc, flux_entropy_correction, alpha_ldg_flip_flop, &
      &Riemann_Diss, Riemann_Diss_BC, WENO_Bias, WENO_type, entropy_viscosity,  &
      &Grid_Topology, cylinder_x0, cylinder_x1, non_conforming, turbulent_viscosity, &
      &p_refine_strategy,entropy_flux,entropy_flux_BC

    namelist /GeneralParameters/ runcasename, solutionfile, outputfile, new, &
      & time_averaging, &
      & write_restart, write_restart_dir, write_restart_common_name, &
      & write_restart_formatted, write_restart_frequency, read_restart_dir, &
      & read_restart_common_name, read_restart_time, read_restart_formatted, &
      & read_restart_time_averaging, &
      & write_solution, write_solution_dir, write_solution_common_name, &
      & write_solution_formatted, write_solution_frequency, &
      & write_aero_coeffs, write_dke_dt, write_errors_wall, &
      & perturb_internal_vertices

    namelist /FlowParameters/ Re0, Re0inv, Pr0, Mach0, uniformFreeStreamAOA, &
      & membranelocation, referenceWaveSpeed, isothermalWallTemperature, &
      & InitialCondition

    namelist /AerodynamicParameters/ aero_coeffs_surface

    continue

    if (myprocid .eq. 0) then
      write(*,*) 'Reading the start-up file'
      write(*,*) '=========================================================================='
    end if


    ! Open the SSDCstartup file
    call get_unit(iunit) ! (precision_vars)
    open(unit = iunit, file = file_name, status = 'old', action = 'read', iostat = ierr)
    if(ierr /= 0)then
!      if (myprocid==0) then
        write(*,*)
        write(*,*) "Fatal Error: the following start-up file does not exist: ", file_name
        write(*,*) 'Exiting...'
!      endif
      stop
    end if

    ! Read the namelist for PreProcParameters
    read(iunit,PreProcParameters,iostat=ierr)
    if((Grid_Topology == 'cylinder') .and. (npoly == 1)) then
      Grid_Topology = 'linear'
      if(myprocid == 0) write(*,*) 'Inconsistancy between Grid_Topology and npoly. Switching Topology to linear'
    endif
    if(ierr /= 0) then
!      if (myprocid==0) then
        write(*,*)
        write(*,*) 'Error in startup file, cannot read PreProcParameters!', ierr
        write(*,*) 'Exiting...'
!      endif
      stop
    end if

    if (myprocid .eq. 0) then
      write(*,*) 'Polynomial order: ', npoly
      write(*,*) '==============================================================='
      write(*,*) 'Over-collocation: ', npoly_DeltaF
      write(*,*) '==============================================================='
    end if

    ! Read the namelist for GeneralParameters
    read(iunit,GeneralParameters,iostat=ierr)
    if(ierr /= 0) then
!      if (myprocid==0) then
        write(*,*)
        write(*,*) 'Error in startup file, cannot read GeneralParameters!', ierr
        write(*,*) 'Exiting...'
!      endif
      stop
    end if 

    ! Read the namelist for GeneralParameters
    read(iunit,FlowParameters,iostat=ierr)
    if(ierr /= 0) then
!      if (myprocid==0) then
        write(*,*)
        write(*,*) 'Error in startup file, cannot read FlowParameters!', ierr
        write(*,*) 'Exiting...'
!      endif
      stop
    end if

    ! Read the namelist for AerodynamicParameters
    read(iunit,AerodynamicParameters,iostat=ierr)
    if(ierr /= 0) then
!      if (myprocid==0) then
        write(*,*)
        write(*,*) 'Error in startup file, cannot read AerodynamicParameters!', ierr
        write(*,*) 'Exiting...'
!      endif
      stop
    end if 

    close(iunit)

    
    ! Check consistency of the input parameters
    ! -----------------------------------------
    ! -----------------------------------------

    ! Check Reynolds number for viscous flow
    if (viscous .and. Re0 < 1e-14) then
      if (myprocid==0) then
        write(*,*)
        write(*,*) 'Error in startup file: You have not specified the', &
          & 'reference Reynolds number or you have set it to zero!' 
        write(*,*) 'Exiting...'
      endif
      stop
    endif

    ! Delay the code so that we can think 2 more seconds if we really want to
    ! overwrite the restart files
    if (write_restart) then
      if (myprocid==0) then
        write(*,*)
        write(*,*) 'WARNING: In 2 seconds the execution of the code will', & 
          & ' continue. If you did not change the restart_dir name or the', & 
          & ' common_restart_name, all the files saved in the default', & 
          & ' directory or in the directory specified in the start-up file', &
          & ' will be overwritten.', &
          & ' Press Ctrl+c to stop the code execution. GOOD LUCK!'
        write(*,*)
      endif

      call sleep(2)
    
    endif
    
    ! Print warning when the user decides to restart a simulation
    if (.not. new) then
      if (myprocid==0) then
        write(*,*)
        write(*,*) 'WARNING: You have chosen to restart the simulation.', &
          & ' Please, verify the details of the restarting files, i.e.,', &
          & ' read_restart_dir, read_restart_common_name, read_restart_time', &
          & ' and read_restart_formatted in the SSDCstartup file.', &
          & ' Press Ctrl+c to stop the code execution.'
        write(*,*)
      endif

      call sleep(2)

    endif

    ! Call run time modifiable startup parameters subroutine
    call read_run_time_modifiable_startup_parameters(file_name)

    return
  end subroutine readSSDCstartup

  !============================================================================
  
  !============================================================================
  
  subroutine read_run_time_modifiable_startup_parameters(file_name)
    
    ! Load modules
    use controlvariables

    ! Nothing is implicitly defined
    implicit none

    character(len=120), intent(in) :: file_name

    integer :: iunit, ierr

    namelist /RunParameters/ timestep, timemaximum , verbose, variabletimestep, CFL, Dt_by_CFL

    continue

    ! Open the SSDCstartup file
    call get_unit(iunit) ! (precision_vars)
    open(unit=iunit,file=file_name,status='old',action='read',iostat=ierr)
    if(ierr /= 0)then
      write(*,*)
      write(*,*) "Fatal Error (modifiable startup parameters): the following start-up file does not exist: ", file_name
      write(*,*) 'Exiting...'
      stop
    end if


    ! read the namelist for RunParameters
    read(iunit,RunParameters,iostat=ierr)
    if(ierr /= 0) then
      write(*,*) 'Error in startup file, cannot read RunParameters!', ierr
      write(*,*) 'Exiting...'
      stop
    end if

    close(iunit)

    return
  end subroutine read_run_time_modifiable_startup_parameters

  !============================================================================
  
  !============================================================================
  ! aflr3ReadUnstructuredGrid - Reads the AFLR3 (.ugrid) grid.

  subroutine aflr3ReadUnstructuredGrid(filein)

    use referencevariables
    use variables, only: vx_master, nqface, if2nq, ifacetag, ic2nh, &
                       & aflr3_bc_ids, periodic_face_data_x1, &
                       & periodic_face_data_x2, periodic_face_data_x3, &
                       & wall_face_data 


    !  nnodesg    = number of nodes
    !  ntface     = number of boundary faces with triangles
    !  nqface     = number of boundary faces with quads
    !  ntet       = number of tetrahedra
    !  npyr       = number of pyramids
    !  nprz       = number of prisms
    !  nhex       = number of hexahedra
    !  x,y,z      = node coordinates
    !  if2nt      = face-to-node pointers for each TRI boundary face (3 nodes
    !               are pointed to by each triangle)
    !  if2nq      = face-to-node pointers for each QUAD boundary face (4 nodes
    !               are pointed to by each quadralateral)
    !  ifacetag   = tag number that groups boundary faces (e.g., all faces
    !               with the number "1" correspond to the same boundary
    !               entity (like wing) with a given boundary condition)
    !  ic2nt      = cell-to-node pointers for tets (4 nodes are pointed to by
    !               each tetrahedral cell)
    !  ic2np      = cell-to-node pointers for pyramids (5 nodes are pointed to by
    !               each pyramidal cell)
    !  ic2nz      = cell-to-node pointers for prisms (6 nodes are pointed to by
    !               each prismatic cell)
    !  ic2nh      = cell-to-node pointers for hexes (8 nodes are pointed to by
    !               each hexahedral cell)
    !  For my use:
    !  ibcnum     = number of surface elements associated with each
    !               number (e.g., ibcnum(1) = number of elements labeled
    !               as BC type=1, ibcnum(n) = number of elements labeled
    !               as BC type=n)
    !  itag       = element numbers of all the surface triangles
    !               (e.g., itag(1,n) = nth triangle element number in group
    !               of BC type = 1, itag(m,n) = nth triangle element number in
    !               group of BC type = m)
    !

!   integer :: nqface
!   integer , allocatable, dimension(:,:) :: if2nq
!   integer , allocatable, dimension(:)   :: ifacetag
!   integer , allocatable, dimension(:,:) :: ic2nh

    ! Nothing is implicitly defined
    implicit none

    character(120), intent(in) :: filein

    integer :: iunit
    integer :: ierr

    integer :: nnodesg, ntface, ntet, npyr, nprz, nhex

    integer , allocatable, dimension(:,:) :: if2nt
    integer , allocatable, dimension(:,:) :: ic2nt
    integer , allocatable, dimension(:,:) :: ic2np
    integer , allocatable, dimension(:,:) :: ic2nz

    integer :: i, j
    integer :: group_id
    integer :: i_face, i_group
    integer :: n_p_faces_x1_a, n_p_faces_x1_b, &
               n_p_faces_x2_a, n_p_faces_x2_b, &
               n_p_faces_x3_a, n_p_faces_x3_b
    integer :: n_p_faces_x1, n_p_faces_x2, n_p_faces_x3
    integer :: n_w_faces
    logical :: match

    continue

    ! Ensure that casefile exists
    call get_unit(iunit)
    open(unit=iunit,file=trim(filein)//'.b8.ugrid', &
!   open(unit=iunit,file=trim(filein)//'.ugrid',    &
      status='old',                                 &
      access='stream',                              &
      !convert='big_endian',                         &
      iostat=ierr)
    if(ierr /= 0)then
      write(*,*) "aflr3GridRead:", trim(filein)," doesn't exist."
      stop
    end if
    close(iunit)

    open(unit=iunit,file=trim(filein)//'.b8.ugrid', &
!   open(unit=iunit,file=trim(filein)//'.ugrid', &
      status='old',                                 &
      access='stream',                              &
      !convert='big_endian',                         &
      iostat=ierr)

    read(iunit) nnodesg, ntface, nqface, ntet, npyr, nprz, nhex
    !write(*,*) nnodesg, ntface, nqface, ntet, npyr, nprz, nhex


    ! Allocate vertices
    allocate(vx_master(1:3,1:nnodesg))

    ! Allocate faces
    if(ntface /= 0) allocate(if2nt(3,ntface))
    if(nqface /= 0) allocate(if2nq(4,nqface))

    if(ntface+nqface /= 0) allocate(ifacetag(ntface+nqface)) 

    ! Allocate elements
    if(ntet /= 0) allocate(ic2nt(4,ntet))
    if(npyr /= 0) allocate(ic2np(5,npyr))
    if(nprz /= 0) allocate(ic2nz(6,nprz))
    if(nhex /= 0) allocate(ic2nh(8,nhex))

    read(iunit) (vx_master(1,i),vx_master(2,i),vx_master(3,i),i=1,nnodesg), &
      ((if2nt(j,i),j=1,3),i=1,ntface),                                      &
      ((if2nq(j,i),j=1,4),i=1,nqface),                                      &
      (ifacetag(i),i=1,ntface+nqface),                                      &
      ((ic2nt(j,i),j=1,4),i=1,ntet),                                        &
      ((ic2np(j,i),j=1,5),i=1,npyr),                                        &
      ((ic2nz(j,i),j=1,6),i=1,nprz),                                        &
      ((ic2nh(j,i),j=1,8),i=1,nhex)

!      write(88,*) ic2nh
!
!    do i = 1,nnodesg
!      write(120,'(i6,1x,3(e15.8,1x))')i,vx(:,i)
!    enddo
!    do i = 1,nqface
!      write(120,'(i6,1x,5(i5,1x))')i,if2nq(:,i),ifacetag(i)
!    enddo

    ! Set number of nodes, vertices per element, and elements (serial)
    nvertices        = nnodesg
    nverticesperelem = 2**ndim
    nelems           = ntet + npyr + nprz + nhex

    ! Read BC information from the aflr3 (.b8.mapbc or .mapbc) grid file
    call aflr3_read_bc(filein)

    ! Re-set ifacetag to be the BC number. This allows to re-use the rest of the
    ! code without changing anything. In fact, ifacetag will be used then to set
    ! ef2e(2,i,j).
    ! =========================================================================
    do i_face = 1, size(ifacetag)

      ! Original face tag ID from AFLR3
      group_id = ifacetag(i_face)

      ! Initialize logical variable
      match = .false.

      do i_group = 1, size(aflr3_bc_ids(:,1))
        if (group_id == aflr3_bc_ids(i_group,1)) then
          ! Match found
          match = .true.

          ! Re-map ifacetage to the convention used in the code
          ifacetag(i_face) = aflr3_bc_ids(i_group,2)
        end if
      end do

      ! If the tag of the boundary face is not known throw an error
      if (match .eqv. .false.) then
        write(*,*) 'Unknown boundary face tag: ', group_id
        write(*,*) 'Please check the number of the boundary face tag.'
        write(*,*) 'Exting...'
        stop
      end if
    end do

    ! Count periodic and wall boundary faces
    ! =========================================================================
    ! Initialize number of periodic and wall boundary faces
    n_p_faces_x1_a = 0
    n_p_faces_x1_b = 0
    n_p_faces_x2_a = 0
    n_p_faces_x2_b = 0
    n_p_faces_x3_a = 0
    n_p_faces_x3_b = 0
    n_p_faces_x1 = 0
    n_p_faces_x2 = 0
    n_p_faces_x3 = 0
    n_w_faces = 0

    do i_face = 1, size(ifacetag)
      ! Count number of "periodic" faces in the x1 direction
      if (ifacetag(i_face) == 8) then 
        n_p_faces_x1_a = n_p_faces_x1_a + 1
        n_p_faces_x1 = n_p_faces_x1 + 1
      end if
      if (ifacetag(i_face) == 9) then  
        n_p_faces_x1_b = n_p_faces_x1_b + 1
        n_p_faces_x1 = n_p_faces_x1 + 1
      end if

      ! Count number of "periodic" faces in the x2 direction
      if (ifacetag(i_face) == 10) then 
        n_p_faces_x2_a = n_p_faces_x2_a + 1
        n_p_faces_x2 = n_p_faces_x2 + 1
      end if
      if (ifacetag(i_face) == 11) then  
        n_p_faces_x2_b = n_p_faces_x2_b + 1
        n_p_faces_x2 = n_p_faces_x2 + 1
      end if

      ! Count number of "periodic" faces in the x3 direction
      if (ifacetag(i_face) == 12) then 
        n_p_faces_x3_a = n_p_faces_x3_a + 1
        n_p_faces_x3 = n_p_faces_x3 + 1
      end if
      if (ifacetag(i_face) == 13) then  
        n_p_faces_x3_b = n_p_faces_x3_b + 1
        n_p_faces_x3 = n_p_faces_x3 + 1
      end if

      ! Count number of "wall" faces
      if (ifacetag(i_face) == 5 .or. ifacetag(i_face) == 6) then
        n_w_faces = n_w_faces + 1
      end if
    
    end do

    ! Check if the number of "periodic" faces in the x1 direction is correct
    if (n_p_faces_x1_a .ne. n_p_faces_x1_b) then
      write(*,*) 'Number of "periodic" faces in the x1_a plane does not match', &
        ' the number of "periodic" faces in the x1_b plane.'
      write(*,*) n_p_faces_x1_a, n_p_faces_x1_b
      write(*,*) 'Exiting...'
      stop
    end if
    
!    ! Print at screen the number of periodic faces in the x1 direction
!    write(*,*) 
!    write(*,*) 'Number of periodic faces in the x1 direction', n_p_faces_x1
!    write(*,*) 


    ! Check if the number of "periodic" faces in the x2 direction is correct
    if (n_p_faces_x2_a .ne. n_p_faces_x2_b) then
      write(*,*) 'Number of "periodic" faces in the x2_a plane does not match', &
        ' the number of "periodic" faces in the x2_b plane.'
      write(*,*) n_p_faces_x2_a, n_p_faces_x2_b
      write(*,*) 'Exiting...'
      stop
    end if

!    ! Print at screen the number of periodic faces in the x2 direction
!    write(*,*) 
!    write(*,*) 'Number of periodic faces in the x2 direction', n_p_faces_x2
!    write(*,*) 


    ! Check if the number of "periodic" faces in the x3 direction is correct
    if (n_p_faces_x3_a .ne. n_p_faces_x3_b) then
      write(*,*) 'Number of "periodic" faces in the x3_a plane does not match', &
        ' the number of "periodic" faces in the x3_b plane.'
      write(*,*) n_p_faces_x3_a, n_p_faces_x3_b
      write(*,*) 'Exiting...'
      stop
    end if

!    ! Print at screen the number of periodic faces in the x3 direction
!    write(*,*) 
!    write(*,*) 'Number of periodic faces in the x3 direction', n_p_faces_x3
!    write(*,*) 


    ! Allocate memory for periodic element face data
    allocate(periodic_face_data_x1(4+nverticesperface,n_p_faces_x1))
    periodic_face_data_x1 = 0

    allocate(periodic_face_data_x2(4+nverticesperface,n_p_faces_x2))
    periodic_face_data_x2 = 0

    allocate(periodic_face_data_x3(4+nverticesperface,n_p_faces_x3))
    periodic_face_data_x3 = 0


    ! Allocate memory for wall element face data
    allocate(wall_face_data(2,n_w_faces))
    wall_face_data = 0

!    write(*,*) 
!    write(*,*) 'Number of wall faces', n_w_faces
!    write(*,*)

    return
  end subroutine aflr3ReadUnstructuredGrid

  !============================================================================
  
  !============================================================================
  ! aflr3_read_bc - Reads the AFLR3 boundary conditions file and call the 
  ! subroutine that maps the AFLR3 boundary conditions numbering to the 
  ! numbering used in the code.

  subroutine aflr3_read_bc(file_in)
    
    ! Load modules
    use variables, only : aflr3_bc_ids 

    ! Nothing is implicitly defined
    implicit none

    character(120), intent(in) :: file_in
    integer :: i_unit
    integer :: i_err
    integer :: n_bc_groups
    integer :: i_group
    integer :: group_id, group_bc_id
    character(120) :: group_bc_name 
    character(120), allocatable, dimension(:) :: aflr3_bc_names

    continue

    ! Ensure that the  mapbc file exists
    ! =========================================================================
    call get_unit(i_unit)

    open(unit=i_unit,file=trim(file_in)//'.mapbc', status='old', &
      & action='READ', iostat=i_err)

    if(i_err /= 0)then
      write(*,*) "aflr3_read_bc:", trim(file_in)," doesn't exist."
      stop
    end if

    close(i_unit)

    ! Read mapbc file
    ! =========================================================================
    open(unit=i_unit,file=trim(file_in)//'.mapbc', status='old', &
      & action='READ', iostat=i_err)

    ! Read number of BC groups
    read(i_unit,*) n_bc_groups

    ! Allocate memory for AFLR3 BC info
    allocate(aflr3_bc_ids(n_bc_groups,2))
    allocate(aflr3_bc_names(n_bc_groups))

    do i_group = 1, n_bc_groups 
      ! Read group informations
      read(i_unit,*) group_id, group_bc_id, group_bc_name

      ! Assign information to the output array
      aflr3_bc_ids(i_group,1) = group_id
      aflr3_bc_ids(i_group,2) = group_bc_id
      aflr3_bc_names(i_group) = group_bc_name
    end do

    close(i_unit)

    ! Re-map the aflr3 group_bc_id to be the BC convention numbering used in
    ! this code
    ! =========================================================================
    call map_aflr3_bc_id(n_bc_groups,aflr3_bc_ids,aflr3_bc_names)

    return
  end subroutine aflr3_read_bc

  !============================================================================
  
  !============================================================================
  ! map_aflr3_bc_id - Maps the AFLR3 BC numbering to the convention used in the
  ! code.

  subroutine map_aflr3_bc_id(n_groups,bc_ids,bc_names)

    ! Load modules
    use variables, only : code_bc_names
    
    ! Nothing is implicitly defined
    implicit none

    integer, intent(in) :: n_groups
    integer, dimension(n_groups,2), intent(inout) :: bc_ids
    character(120), dimension(n_groups), intent(in) :: bc_names
    character(120) :: group_bc_name
    logical :: match
    integer :: i_group, i_code_bc

    continue

    do i_group = 1, n_groups
      ! Group BC name from AFLR3 grid
      group_bc_name = bc_names(i_group)

      ! Initialize logical variable to false
      match = .false.

      do i_code_bc = 1, size(code_bc_names)
        ! Check if the group BC name from AFLR3 matches one of the BC supported
        ! in the code
        if (group_bc_name == code_bc_names(i_code_bc)) then
          bc_ids(i_group,2) = i_code_bc
          match = .true.
        end if
      end do
      if (match .eqv. .false.) then
        write(*,*) 'Unknown BC type: ', group_bc_name
        write(*,*) 'Please check the name of the boundary condition or', &
          & ' implement this new type.'
        write(*,*) 'Exting...'
        stop
      end if
    end do

    return
  end subroutine map_aflr3_bc_id

  !============================================================================
  
  !============================================================================
  ! init_code_bc_types - Defines the IDs of the boundary conditions according to
  ! this code.

  subroutine init_code_bc_types()
    
    ! Load modules
    use variables, only : code_bc_names

    ! Nothing is implicitly defined
    implicit none

    integer, parameter :: n_code_bc_types = 15

    continue 

    ! Allocate memory for code BC names and IDs
    allocate(code_bc_names(n_code_bc_types))
    
    ! No BC
    code_bc_names(1) = 'null'

    ! Dirichlet BC. The actual BC subroutine will depends on the initial
    ! solution specific in the start-up file
    code_bc_names(2) = 'dirichlet'

    ! Uniform free-stream BC
    code_bc_names(3) = 'uniform_free_stream'

     ! Inviscid wall BC
    code_bc_names(4) = 'inviscid_wall'

    ! No slip adiabatic wall BC
    code_bc_names(5) = 'no_slip_adiabatic_wall'

    ! No slip wall entropy stable BC
    code_bc_names(6) = 'no_slip_wall_ss'

    ! Symmetry plane BC
    code_bc_names(7) = 'symmetry_plane'

    ! Periodic BC: 1a plane 
    code_bc_names(8) = 'periodic_x1_a'
    
    ! Periodic BC: 1b plane
    code_bc_names(9) = 'periodic_x1_b'

    ! Periodic BC: 2a plane 
    code_bc_names(10) = 'periodic_x2_a'
    
    ! Periodic BC: 2b plane
    code_bc_names(11) = 'periodic_x2_b'

    ! Periodic BC: 3a plane 
    code_bc_names(12) = 'periodic_x3_a'
    
    ! Periodic BC: 3b plane
    code_bc_names(13) = 'periodic_x3_b'

    ! Boundary Layer Inflow Profile
    code_bc_names(14) = 'BLayerProfile'

    ! Subsonic Pressure Boundary condition 
    code_bc_names(15) = 'SubsonicOutflow'

    return
  end subroutine init_code_bc_types
          
  !============================================================================
  
end module fileio
