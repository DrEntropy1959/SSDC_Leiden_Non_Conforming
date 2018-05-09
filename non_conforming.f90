module non_conforming

  use precision_vars
  use initcollocation

  !-- Nothing is implicitely defined
  implicit none

  !-- private subroutines functions etc
  private

  !-- public subroutines functions etc
  public Lagrange_interpolant_basis_1D
  public Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo
  public Vandermonde_1D_monomial
  public Rotate_xione_2_xitwo_and_back
  public Rotate_GL_2_G_and_back_I
  public h_refine
  public construct_h_refine_list
contains
!==================================================================================================
!
!
! Purpose: This subroutine 
! populates h_refine_list which is used to refine elements
!
!
!==================================================================================================
subroutine construct_h_refine_list
     use referencevariables, only : nelems, nfacesperelem
     use variables, only : ef2e, h_refine_list, nelems_to_refine

     !-- local variables
     integer :: ielem, iface, refine_method
     integer :: i_err

     !-- allocate list of which elements to refine
     allocate(h_refine_list(nelems))
     h_refine_list(1:nelems) = .false.

     refine_method = 2
     select case(refine_method)
     case(1)
       !-- loop over the elements identify which ones are boundary elements and
       !tag for refinment
       nelems_to_refine = 0
       e_loop_bc1 : do ielem = 1,nelems
         f_loop_bc1: do iface = 1,nfacesperelem
           !-- check if the element has a face with a boundary condition
           f_if_bc1: if (ef2e(1,iface,ielem) < 0) then
             h_refine_list(ielem) = .true.
             nelems_to_refine = nelems_to_refine+1
             !-- exit f_loop_bc2
             exit
           endif f_if_bc1
         end do f_loop_bc1
       enddo e_loop_bc1
     case(2)
       nelems_to_refine = 2
       h_refine_list(4) = .true.
       h_refine_list(9) = .true.
     case default
       write(*,*)'non_conforming: construct_h_refine_list: incorrect choice of refine_method = ',refine_method
       call PetscFinalize(i_err); stop
     end select

end subroutine construct_h_refine_list
!==================================================================================================
!
!
! Purpose: This subroutine h refines the mesh at the boundaries
! 
! Comments: Only setup for 3D. The connectivity in e_edge2e is setup assuming only h refinement. 
! Numbering for the splitting: we only consider a equal subvdivision
! of the hex in each direction. The nod numbering convention used here is given
! in the following figure:
!
!                        8-------15--------7
!                       /.                /|
!                      / .               / |
!                     /  .              /  |
!                    /                 /   |
!                   16   .            14   |
!                  /     20          /    19
!                 /      .          /      |
!                /       .         /       |
!               5---------13------6        |       
!               |        4.....11.|........3
!               |       .         |       /
!               |      .          |      /
!               17     .          18    /  zeta ^     
!               |    12           |    10       |    / eta
!               |   .             |   /         |   /
!               |  .              |  /          |  /
!               | .               | /           | /
!               |.                |/            |/
!               1-------9---------2             ----------> xi
!               5---------13------6    6---------14------7   7---------15------8     
!               |        |        |    |        |        |   |        |        |  
!               |        |        |    |        |        |   |        |        |  
!               |        |        |    |        |        |   |        |        |  
!               |        |        |    |        |        |   |        |        |
!               17       21       18   18       22       19  19       23       20
!               |        |        |    |        |        |   |        |        |
!               |        |        |    |        |        |   |        |        |
!               |        |        |    |        |        |   |        |        |
!               |        |        |    |        |        |   |        |        |
!               1-------9---------2    4-------10--------3   3-------11--------4   
!               8---------16------5    4---------11------3   8---------15------7     
!               |        |        |    |        |        |   |        |        |  
!               |        |        |    |        |        |   |        |        |  
!               |        |        |    |        |        |   |        |        |  
!               |        |        |    |        |        |   |        |        |
!               20       24       17   12       25       10  16       26       14
!               |        |        |    |        |        |   |        |        |
!               |        |        |    |        |        |   |        |        |
!               |        |        |    |        |        |   |        |        |
!               |        |        |    |        |        |   |        |        |
!               4-------12---------1   1-------9--------2    5-------13--------6   


!               the node in the middle of the volume is numbered 27
!
!               Each element is split into 8, in the figures below, in local numbering (1:8)
!               shows which sub-element touches which face
! 
!               face 1                 face 2                face 3
!               -------------------    -------------------   -------------------     
!               |        |        |    |        |        |   |        |        |  
!               |   4    |   3    |    |   5    |    6   |   |   6    |   7    |  
!               |        |        |    |        |        |   |        |        |  
!               |        |        |    |        |        |   |        |        |
!               -------------------    -------------------   -------------------
!               |        |        |    |        |        |   |        |        |
!               |   1    |    2   |    |    1   |    2   |   |    2   |    3   |
!               |        |        |    |        |        |   |        |        |
!               |        |        |    |        |        |   |        |        |
!               -------------------    -------------------   -------------------   
!
!               face 4                 face 5                face 6
!               -------------------    -------------------   -------------------     
!               |        |        |    |        |        |   |        |        |  
!               |   8    |   7    |    |   5    |    8   |   |   8    |   7    |  
!               |        |        |    |        |        |   |        |        |  
!               |        |        |    |        |        |   |        |        |
!               -------------------    -------------------   -------------------
!               |        |        |    |        |        |   |        |        |
!               |   4    |    3   |    |    1   |    4   |   |    5   |    6   |
!               |        |        |    |        |        |   |        |        |
!               |        |        |    |        |        |   |        |        |
!               -------------------    -------------------   -------------------   
! Outputs:
!         
!
!===================================================================================================
   subroutine h_refine()

     use referencevariables, only : nelems, nfacesperelem, ndim
     use variables, only : ef2e, e_edge2e, h_refine_list, nelems_to_refine
     use variables, only : vx_master, ic2nh

     !-- local variables
     integer :: ielem, iface, iedge, ipartner, max_partners, nvertex,&
                element_count, vertex_count, refine,ielem_old
     integer :: nelems_old, nvertex_old, nfaces_old
     integer, allocatable, dimension(:,:,:) :: ef2e_temp                                            !-- (7 ,nfaceperelem, nelements)
     integer, allocatable, dimension(:,:,:,:,:) :: e_edge2e_temp                                    !-- (3,number_of_edges_per_face,max_partners,nfaceperelem,nelems)
     integer, allocatable, dimension(:,:) :: e_old2e
     integer, allocatable, dimension(:,:) :: ic2nh_temp
     integer, allocatable, dimension(:) :: e2refine

     integer :: v1n, v2n, v3n, v4n, v5n, v6n, v7n, v8n,&
                v9n, v10n, v11n, v12n, v13n, v14n, v15n, v16n,&
                v17n, v18n, v19n, v20n, v21n, v22n, v23n, v24n,&
                v25n, v26n, v27n
     integer :: ielem1, ielem2, ielem3, ielem4, ielem5, ielem6, ielem7, ielem8,&
                kface, kelem1, kelem2, kelem3, kelem4, kelem_old,&
                kface1, kface2, kface3,kface4
     integer :: element_number_adjoining, faceID_adjoining, split_count
     integer :: face2elem(6,4), kelemsv(4), ielemsv(4), kfacesv(4)


     real(wp), dimension (3) :: xyz1, xyz2, xyz3, xyz4,& 
                                xyz5, xyz6, xyz7, xyz8,&
                                xyz9, xyz10, xyz11, xyz12,&
                                xyz13, xyz14, xyz15, xyz16,&
                                xyz17, xyz18, xyz19, xyz20,&
                                xyz21, xyz22, xyz23, xyz24,&
                                xyz25, xyz26, xyz27
     real(wp), allocatable, dimension(:,:) :: vx_master_temp                                         !-- (3,number of vertices)

     !-- update the number of faces per element (this is a maximum for one level refinment)
     nfaces_old = nfacesperelem
     nfacesperelem = 4*6
     
     !-- allocate the necessary temporary arrays
     nelems_old = nelems
     nelems = nelems-nelems_to_refine+8*nelems_to_refine
     nvertex_old = size(vx_master(1,:)) 
     nvertex = nvertex_old+19*nelems_to_refine
     max_partners = size(e_edge2e(1,1,:,1,1))
                                     
     allocate(vx_master_temp(3,nvertex))
     allocate(ef2e_temp(7,nfacesperelem,nelems))     
     allocate(e_edge2e_temp(3,2**(ndim-1),max_partners,nfacesperelem,nelems))
     allocate(ic2nh_temp(8,nelems))
     allocate(e_old2e(nelems_to_refine,8))
     allocate(e2refine(nelems_old))

     !-- store old information
     vx_master_temp(1:3,1:nvertex_old) = vx_master(1:3,1:nvertex_old)
     ic2nh_temp(1:8,1:nelems_old) = ic2nh(1:8,1:nelems_old)
     ef2e_temp = 0
     ef2e_temp(1:7,1:nfaces_old,1:nelems_old) = ef2e(1:7,1:nfaces_old,1:nelems_old)
     e_edge2e_temp(1:3,1:2**(ndim-1),1:max_partners,1:nfaces_old,1:nelems_old) = &
       e_edge2e(1:3,1:2**(ndim-1),1:max_partners,1:nfaces_old,1:nelems_old)

     !-- loop over the element and split
     element_count = nelems_old+1
     vertex_count = nvertex_old+1
     split_count = 1
     e_loop_bc3 : do ielem = 1,nelems_old

       !-- element to be split
       if(h_refine_list(ielem))then
       !-- vertex number of the original 8 vertices
       v1n = ic2nh(1,ielem);v2n = ic2nh(2,ielem); v3n = ic2nh(3,ielem); v4n = ic2nh(4,ielem);
       v5n = ic2nh(5,ielem);v6n = ic2nh(6,ielem); v7n = ic2nh(7,ielem); v8n = ic2nh(8,ielem);

       !-- spatial locations of the original 8 vertices
       xyz1 = Vx_master(1:3,v1n); xyz2 = Vx_master(1:3,v2n); xyz3 = Vx_master(1:3,v3n); xyz4 = Vx_master(1:3,v4n)
       xyz5 = Vx_master(1:3,v5n); xyz6 = Vx_master(1:3,v6n); xyz7 = Vx_master(1:3,v7n); xyz8 = Vx_master(1:3,v8n)

         !-- split the element


         !-- construct new vertices (see figure at top)
         xyz9 = 0.5_wp*(xyz1+xyz2); xyz10 = 0.5_wp*(xyz2+xyz3); xyz11 = 0.5_wp*(xyz3+xyz4); xyz12 = 0.5_wp*(xyz4+xyz1) 
         xyz13 = 0.5_wp*(xyz5+xyz6); xyz14 = 0.5_wp*(xyz6+xyz7); xyz15 = 0.5_wp*(xyz7+xyz8); xyz16 = 0.5_wp*(xyz8+xyz5)
         xyz17 = 0.5_wp*(xyz1+xyz5); xyz18 = 0.5_wp*(xyz2+xyz6); xyz19 = 0.5_wp*(xyz3+xyz7); xyz20 = 0.5_wp*(xyz4+xyz8)
         xyz21 = 0.5_wp*(xyz9+xyz13); xyz22 = 0.5_wp*(xyz10+xyz14); xyz23 = 0.5_wp*(xyz11+xyz15); xyz24 = 0.5_wp*(xyz12+xyz16)
         xyz25 = 0.5_wp*(xyz9+xyz11); xyz26 = 0.5_wp*(xyz13+xyz15);xyz27 = 0.5_wp*(xyz21+xyz23)

         !-- set the vertex numbers
         v9n = vertex_count; v10n = vertex_count+1; v11n = vertex_count+2; v12n = vertex_count+3;
         v13n = vertex_count+4; v14n = vertex_count+5; v15n = vertex_count+6; v16n = vertex_count+7;
         v17n = vertex_count+8; v18n = vertex_count+9; v19n = vertex_count+10; v20n = vertex_count+11;
         v21n = vertex_count+12; v22n = vertex_count+13; v23n = vertex_count+14; v24n = vertex_count+15;
         v25n = vertex_count+16; v26n = vertex_count+17; v27n = vertex_count+18;

         !-- store the vertices
         vx_master_temp(1:3,v1n) = xyz1;vx_master_temp(1:3,v2n) = xyz2;vx_master_temp(1:3,v3n) = xyz3;
         vx_master_temp(1:3,v4n) = xyz4;
         vx_master_temp(1:3,v5n) = xyz5;vx_master_temp(1:3,v6n) = xyz6;vx_master_temp(1:3,v7n) = xyz7;
         vx_master_temp(1:3,v8n) = xyz8; 
         vx_master_temp(1:3,v9n) = xyz9;vx_master_temp(1:3,v10n) = xyz10;vx_master_temp(1:3,v11n) = xyz11;
         vx_master_temp(1:3,v12n) = xyz12;  
         vx_master_temp(1:3,v13n) = xyz13;vx_master_temp(1:3,v14n) = xyz14;vx_master_temp(1:3,v15n) = xyz15;
         vx_master_temp(1:3,v16n) = xyz16; 
         vx_master_temp(1:3,v17n) = xyz17;vx_master_temp(1:3,v18n) = xyz18;vx_master_temp(1:3,v19n) = xyz19;
         vx_master_temp(1:3,v20n) = xyz20; 
         vx_master_temp(1:3,v21n) = xyz21;vx_master_temp(1:3,v22n) = xyz22;vx_master_temp(1:3,v23n) = xyz23;
         vx_master_temp(1:3,v24n) = xyz24;
         vx_master_temp(1:3,v25n) = xyz25;vx_master_temp(1:3,v26n) = xyz26;vx_master_temp(1:3,v27n) = xyz27;

         !-- store the vertex to element connectivity
         ic2nh_temp(1:8,ielem) = (/v1n,v9n,v25n,v12n,v17n,v21n,v27n,v24n/)
         ic2nh_temp(1:8,element_count) = (/V9n,V2n,V10n,V25n,v21n,v18n,v22n,v27n/)
         ic2nh_temp(1:8, element_count+1) = (/v25n,v10n,v3n,v11n,v27n,v22n,v19n,v23n/)
         ic2nh_temp(1:8, element_count+2) = (/v12n,v25n,v11n,v4n,v24n,v27n,v23n,v20n/)
         ic2nh_temp(1:8, element_count+3) = (/v17n,v21n,v27n,v24n,v5n,v13n,v26n,v16n/)
         ic2nh_temp(1:8, element_count+4) = (/v21n,v18n,v22n,v27n,v13n,v6n,v14n,v26n/)
         ic2nh_temp(1:8, element_count+5) = (/v27n,v22n,v19n,v23n,v26n,v14n,v7n,v15n/)
         ic2nh_temp(1:8, element_count+6) = (/v24n,v27n,v23n,v20n,v16n,v26n,v15n,v8n/)

         !-- update maping for split elements
         e_old2e(split_count,1:8) =&
           (/ielem,element_count,element_count+1,element_count+2,element_count+3,element_count+4,element_count+5,element_count+6/)
         !-- update map from old element numbering to split element map (e_old2e)
         e2refine(ielem) = split_count
        
         !-- update counters
         vertex_count = vertex_count+19
         element_count = element_count+7
         split_count = split_count+1
       endif
     enddo e_loop_bc3

     !-- construct face map which gives the local element numbers for a given face on a split element
     !   in counter-clockwise ordering (see the figures at the start)
     face2elem(1,1:4) = (/1,2,3,4/)
     face2elem(2,1:4) = (/1,2,6,5/)
     face2elem(3,1:4) = (/2,3,7,6/)
     face2elem(4,1:4) = (/4,3,7,8/)
     face2elem(5,1:4) = (/1,4,8,5/)
     face2elem(6,1:4) = (/5,6,7,8/)
!write(*,*)'vx_master_temp = [...'
!do ielem = 1,nvertex
!   write(*,*)vx_master_temp(1:3,ielem),';'
!enddo
!write(*,*)'];'
!
!write(*,*)'vx_master_old = [...'
!do ielem = 1,nvertex_old
!   write(*,*)vx_master(1:3,ielem),';'
!enddo
!write(*,*)'];'
!
!write(*,*)'ic2nh_temp = [...'
!do refine = 1,8
!   write(*,*)ic2nh_temp(refine,1:nelems),';'
!enddo
!write(*,*)'];'
!
!write(*,*)'ic2nh_old = [...'
!do refine = 1,8
!   write(*,*)ic2nh(refine,1:nelems_old),';'
!enddo
!write(*,*)'];'
!-- the edge numbering is counter-clockwise starting from the origin i.e.
!                 edge 3
!               ----------
!               |        |
!        edge4  |        | edge 2
!               |        |
!               |        |
!               ----------
!               edge 1

    !-- update ef2e
     ef2e_loop: do refine = 1,nelems_to_refine
          ielem_old = e_old2e(refine,1); ielem2 = e_old2e(refine,2); ielem3 = e_old2e(refine,3); ielem4 = e_old2e(refine,4)
          ielem5 = e_old2e(refine,5); ielem6 = e_old2e(refine,6); ielem7 = e_old2e(refine,7); ielem8 = e_old2e(refine,8)

          !-- internal connections 
          !-- element 1
          ef2e_temp(1,3,ielem_old) = 5; ef2e_temp(2,3,ielem_old) = ielem2
          ef2e_temp(1,4,ielem_old) = 2; ef2e_temp(2,4,ielem_old) = ielem4
          ef2e_temp(1,6,ielem_old) = 1; ef2e_temp(2,6,ielem_old) = ielem5
        
          !-- element 2
          ef2e_temp(1,4,ielem2) = 2; ef2e_temp(2,4,ielem2) = ielem3
          ef2e_temp(1,5,ielem2) = 3; ef2e_temp(2,5,ielem2) = ielem_old
          ef2e_temp(1,6,ielem2) = 1; ef2e_temp(2,6,ielem2) = ielem6

          !-- element 3
          ef2e_temp(1,2,ielem3) = 4; ef2e_temp(2,2,ielem3) = ielem2
          ef2e_temp(1,5,ielem3) = 3; ef2e_temp(2,5,ielem3) = ielem4
          ef2e_temp(1,6,ielem3) = 1; ef2e_temp(2,6,ielem3) = ielem7

          !-- element 4
          ef2e_temp(1,2,ielem4) = 4; ef2e_temp(2,2,ielem4) = ielem_old
          ef2e_temp(1,3,ielem4) = 5; ef2e_temp(2,3,ielem4) = ielem3
          ef2e_temp(1,6,ielem4) = 1; ef2e_temp(2,6,ielem4) = ielem8

          !-- element 5
          ef2e_temp(1,1,ielem5) = 6; ef2e_temp(2,1,ielem5) = ielem_old
          ef2e_temp(1,3,ielem5) = 5; ef2e_temp(2,3,ielem5) = ielem6
          ef2e_temp(1,4,ielem5) = 2; ef2e_temp(2,4,ielem5) = ielem8

          !-- element 6
          ef2e_temp(1,1,ielem6) = 6; ef2e_temp(2,1,ielem6) = ielem2
          ef2e_temp(1,4,ielem6) = 2; ef2e_temp(2,4,ielem6) = ielem7
          ef2e_temp(1,5,ielem6) = 3; ef2e_temp(2,5,ielem6) = ielem5

          !-- element 7
          ef2e_temp(1,1,ielem7) = 6; ef2e_temp(2,1,ielem7) = ielem3
          ef2e_temp(1,2,ielem7) = 4; ef2e_temp(2,2,ielem7) = ielem6
          ef2e_temp(1,5,ielem7) = 3; ef2e_temp(2,5,ielem7) = ielem8

          !-- element 8
          ef2e_temp(1,1,ielem8) = 6; ef2e_temp(2,1,ielem8) = ielem4
          ef2e_temp(1,2,ielem8) = 4; ef2e_temp(2,2,ielem8) = ielem5
          ef2e_temp(1,3,ielem8) = 5; ef2e_temp(2,3,ielem8) = ielem7

          !-- e_edge2e connections
          !-- new element 1 face 3 edge 2
          !e_edge2e_temp(1,2,1,3,ielem_old) = ielem2;e_edge2e_temp(1,2,2,3,ielem_old) = ielem3;e_edge2e_temp(1,2,3,3,ielem_old) = ielem4 
          !e_edge2e_temp(2,2,1:3,3,ielem_old) = ef2e(6,3,ielem_old)

          !-- new element 1 face 3 edge 3
          !e_edge2e_temp(1,3,1,3,ielem_old) = ielem2;e_edge2e_temp(1,3,2,3,ielem_old) = ielem5;e_edge2e_temp(1,3,3,3,ielem_old) = ielem6 
          !e_edge2e_temp(2,3,1:3,3,ielem_old) = ef2e(6,3,ielem_old)         

          !-- new element 1 face 4 edge 2
          !e_edge2e_temp(1,2,1,4,ielem_old) = ielem2;e_edge2e_temp(1,2,2,4,ielem_old) = ielem3;e_edge2e_temp(1,2,3,4,ielem_old) = ielem4 
          !e_edge2e_temp(2,2,1:3,4,ielem_old) = ef2e(6,4,ielem_old)

          !-- new element 1 face 4 edge 3
          !e_edge2e_temp(1,3,1,4,ielem_old) = ielem4;e_edge2e_temp(1,3,2,4,ielem_old) = ielem5;e_edge2e_temp(1,3,3,4,ielem_old) = ielem8 
          !e_edge2e_temp(2,3,1:3,4,ielem_old) = ef2e(6,4,ielem_old)       

          !-- new element 1 face 6 edge 2
          !e_edge2e_temp(1,2,1,6,ielem_old) = ielem2;e_edge2e_temp(1,2,2,6,ielem_old) = ielem5;e_edge2e_temp(1,2,3,6,ielem_old) = ielem6 
          !e_edge2e_temp(2,2,1:3,6,ielem_old) = ef2e(6,6,ielem_old)

          !-- new element 1 face 6 edge 3
          !e_edge2e_temp(1,3,1,6,ielem_old) = ielem4;e_edge2e_temp(1,3,2,6,ielem_old) = ielem5;e_edge2e_temp(1,3,3,6,ielem_old) = ielem8 
          !e_edge2e_temp(2,3,1:3,6,ielem_old) = ef2e(6,6,ielem_old)

          !-- loop over the original faces and determine the connectivity
          iface1: do iface = 1,nfaces_old
            bc_if: if(ef2e(1,iface,ielem_old)<0)then
              !-- this is a boundary face

              if(iface.EQ.1)then
                ef2e_temp(1,iface,ielem_old) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem_old) = ef2e(2,iface,ielem_old)
                ef2e_temp(1,iface,ielem2) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem2) = ef2e(2,iface,ielem_old)
                ef2e_temp(1,iface,ielem3) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem3) = ef2e(2,iface,ielem_old)
                ef2e_temp(1,iface,ielem4) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem4) = ef2e(2,iface,ielem_old)
              elseif(iface.EQ.2)then
                ef2e_temp(1,iface,ielem_old) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem_old) = ef2e(2,iface,ielem_old)
                ef2e_temp(1,iface,ielem2) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem2) = ef2e(2,iface,ielem_old)
                ef2e_temp(1,iface,ielem6) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem6) = ef2e(2,iface,ielem_old)
                ef2e_temp(1,iface,ielem5) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem5) = ef2e(2,iface,ielem_old)
              elseif(iface.EQ.3)then
                ef2e_temp(1,iface,ielem2) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem2) = ef2e(2,iface,ielem_old)
                ef2e_temp(1,iface,ielem3) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem3) = ef2e(2,iface,ielem_old)
                ef2e_temp(1,iface,ielem7) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem7) = ef2e(2,iface,ielem_old)
                ef2e_temp(1,iface,ielem6) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem6) = ef2e(2,iface,ielem_old)
              elseif(iface.EQ.4)then
                ef2e_temp(1,iface,ielem4) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem4) = ef2e(2,iface,ielem_old)
                ef2e_temp(1,iface,ielem3) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem3) = ef2e(2,iface,ielem_old)
                ef2e_temp(1,iface,ielem7) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem7) = ef2e(2,iface,ielem_old)
                ef2e_temp(1,iface,ielem8) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem8) = ef2e(2,iface,ielem_old)
              elseif(iface.EQ.5)then
                ef2e_temp(1,iface,ielem_old) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem_old) = ef2e(2,iface,ielem_old)
                ef2e_temp(1,iface,ielem4) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem4) = ef2e(2,iface,ielem_old)
                ef2e_temp(1,iface,ielem8) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem8) = ef2e(2,iface,ielem_old)
                ef2e_temp(1,iface,ielem5) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem5) = ef2e(2,iface,ielem_old)
              elseif(iface.EQ.6)then
                ef2e_temp(1,iface,ielem5) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem5) = ef2e(2,iface,ielem_old)
                ef2e_temp(1,iface,ielem6) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem6) = ef2e(2,iface,ielem_old)
                ef2e_temp(1,iface,ielem7) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem7) = ef2e(2,iface,ielem_old)
                ef2e_temp(1,iface,ielem8) = ef2e(1,iface,ielem_old); ef2e_temp(2,iface,ielem8) = ef2e(2,iface,ielem_old)
              endif
            elseif(h_refine_list(ef2e(2,iface,ielem_old)))then
              !-- both faces have been split

                kface = ef2e(1,iface,ielem_old)
                kelem_old = ef2e(2,iface,ielem_old) 
                kelem1 = e_old2e(e2refine(kelem_old),face2elem(kface,1))
                kelem2 = e_old2e(e2refine(kelem_old),face2elem(kface,2))
                kelem3 = e_old2e(e2refine(kelem_old),face2elem(kface,3))
                kelem4 = e_old2e(e2refine(kelem_old),face2elem(kface,4))

                ielem1 = e_old2e(refine,face2elem(iface,1)); ielem2 = e_old2e(refine,face2elem(iface,2))
                ielem3 = e_old2e(refine,face2elem(iface,3)); ielem4 = e_old2e(refine,face2elem(iface,4))
 
              if(ef2e(7,iface,ielem_old).EQ.0)then
                kelemsv = (/kelem1,kelem2,kelem3,kelem4/)
              elseif(ef2e(7,iface,ielem_old).EQ.1)then
                kelemsv =(/kelem2,kelem3,kelem4,kelem1/)
              elseif(ef2e(7,iface,ielem_old).EQ.2)then
                kelemsv = (/kelem3,kelem4,kelem1,kelem2/)
              elseif(ef2e(7,iface,ielem_old).EQ.3)then
                kelemsv = (/kelem4,kelem1,kelem2,kelem3/)
              elseif(ef2e(7,iface,ielem_old).EQ.4)then
                kelemsv = (/kelem4,kelem3,kelem2,kelem1/)
              elseif(ef2e(7,iface,ielem_old).EQ.5)then
                kelemsv = (/kelem1,kelem4,kelem3,kelem2/)
              elseif(ef2e(7,iface,ielem_old).EQ.6)then
                kelemsv = (/kelem2,kelem1,kelem4,kelem3/)
              elseif(ef2e(7,iface,ielem_old).EQ.7)then
                kelemsv = (/kelem3,kelem2,kelem1,kelem4/)
              endif
  
              !-- update connectivity             
              ef2e_temp(1,iface,ielem1) = kface; ef2e_temp(2,iface,ielem1) = kelemsv(1) 
              ef2e_temp(7,iface,ielem1) = ef2e(7,iface,ielem_old)
              ef2e_temp(1,iface,ielem2) = kface; ef2e_temp(2,iface,ielem2) = kelemsv(2) 
              ef2e_temp(7,iface,ielem2) = ef2e(7,iface,ielem_old)
              ef2e_temp(1,iface,ielem3) = kface; ef2e_temp(2,iface,ielem3) = kelemsv(3)
              ef2e_temp(7,iface,ielem3) = ef2e(7,iface,ielem_old)
              ef2e_temp(1,iface,ielem4) = kface; ef2e_temp(2,iface,ielem4) = kelemsv(4)
              ef2e_temp(7,iface,ielem4) = ef2e(7,iface,ielem_old)

            else
              !-- nonconforming face
              kface = ef2e(1,iface,ielem_old)
              kface1 = kface; kface2 = kface+6; kface3 = kface+12; kface4 = kface+18
              kelem_old = ef2e(2,iface,ielem_old) 

              ielem1 = e_old2e(refine,face2elem(iface,1)); ielem2 = e_old2e(refine,face2elem(iface,2))
              ielem3 = e_old2e(refine,face2elem(iface,3)); ielem4 = e_old2e(refine,face2elem(iface,4))


              if(ef2e(7,iface,ielem_old).EQ.0)then
                ielemsv = (/ielem1,ielem2,ielem3,ielem4/)
                kfacesv = (/kface1,kface2,kface3,kface4/)
              elseif(ef2e(7,iface,ielem_old).EQ.1)then
                ielemsv =(/ielem2,ielem3,ielem4,ielem1/)
                kfacesv =(/kface2,kface3,kface4,kface1/)
              elseif(ef2e(7,iface,ielem_old).EQ.2)then
                ielemsv = (/ielem3,ielem4,ielem1,ielem2/)
                kfacesv = (/kface3,kface4,kface1,kface2/)
              elseif(ef2e(7,iface,ielem_old).EQ.3)then
                ielemsv = (/ielem4,ielem1,ielem2,ielem3/)
                kfacesv = (/kface4,kface1,kface2,kface3/)
              elseif(ef2e(7,iface,ielem_old).EQ.4)then
                ielemsv = (/ielem4,ielem3,ielem2,ielem1/)
                kfacesv = (/kface4,kface3,kface2,kface1/)
              elseif(ef2e(7,iface,ielem_old).EQ.5)then
                ielemsv = (/ielem1,ielem4,ielem3,ielem2/)
                kfacesv = (/kface1,kface4,kface3,kface2/)
              elseif(ef2e(7,iface,ielem_old).EQ.6)then
                ielemsv = (/ielem2,ielem1,ielem4,ielem3/)
                kfacesv = (/kface2,kface1,kface4,kface3/)
              elseif(ef2e(7,iface,ielem_old).EQ.7)then
                ielemsv = (/ielem3,ielem2,ielem1,ielem4/)
                kfacesv = (/kface3,kface2,kface1,kface4/)
              endif

              !-- populate new elements
              ef2e_temp(1,iface,ielem1) = kfacesv(1); ef2e_temp(2,iface,ielem1) = kelem_old
              ef2e_temp(7,iface,ielem1) = ef2e(7,iface,ielem_old)
              ef2e_temp(1,iface,ielem2) = kfacesv(2); ef2e_temp(2,iface,ielem2) = kelem_old
              ef2e_temp(7,iface,ielem2) = ef2e(7,iface,ielem_old)
              ef2e_temp(1,iface,ielem3) = kfacesv(3); ef2e_temp(2,iface,ielem3) = kelem_old
              ef2e_temp(7,iface,ielem3) = ef2e(7,iface,ielem_old)
              ef2e_temp(1,iface,ielem4) = kfacesv(4); ef2e_temp(2,iface,ielem4) = kelem_old
              ef2e_temp(7,iface,ielem4) = ef2e(7,iface,ielem_old)

              !-- populate subfaces of kelem_old
              ef2e_temp(1,kface1,kelem_old) = iface; ef2e_temp(2,kface1,kelem_old) = ielemsv(1)
              ef2e_temp(7,kface1,kelem_old) = ef2e(7,kface,kelem_old)
              ef2e_temp(1,kface2,kelem_old) = iface; ef2e_temp(2,kface2,kelem_old) = ielemsv(2)
              ef2e_temp(7,kface2,kelem_old) = ef2e(7,kface,kelem_old)
              ef2e_temp(1,kface3,kelem_old) = iface; ef2e_temp(2,kface3,kelem_old) = ielemsv(3)
              ef2e_temp(7,kface3,kelem_old) = ef2e(7,kface,kelem_old)
              ef2e_temp(1,kface4,kelem_old) = iface; ef2e_temp(2,kface4,kelem_old) = ielemsv(4)
              ef2e_temp(7,kface4,kelem_old) = ef2e(7,kface,kelem_old)

            endif bc_if
          enddo iface1
      enddo ef2e_loop

     !-- assigne temp arrays to arrays used in main code
     deallocate(ef2e); allocate(ef2e(7,nfacesperelem,nelems))
     ef2e(:,:,:) = ef2e_temp(:,:,:)
     deallocate(e_edge2e);allocate(e_edge2e(3,2**(ndim-1),max_partners,nfacesperelem,nelems))
     e_edge2e(:,:,:,:,:) = e_edge2e_temp(:,:,:,:,:)
     deallocate(vx_master); allocate(vx_master(3,nvertex))
     vx_master(:,:) = vx_master_temp(:,:)
     deallocate(ic2nh);allocate(ic2nh(8,nelems))
     ic2nh(:,:) = ic2nh_temp(:,:)
do ielem = 1,nelems
write(*,*)'========================'
do iface = 1,nfacesperelem
write(*,*)'ielem = ',ielem,'iface = ',iface
write(*,*)'ef2e(1,iface,ielem) = ',ef2e(1,iface,ielem)
write(*,*)'ef2e(2,iface,ielem) = ',ef2e(2,iface,ielem)
write(*,*)'ef2e(7,iface,ielem) = ',ef2e(7,iface,ielem)
enddo
write(*,*)'========================'
enddo     
     !-- deallocate statements
     deallocate(ef2e_temp)
     deallocate(e_edge2e_temp)
     deallocate(vx_master_temp)
     deallocate(ic2nh_temp)
   end subroutine h_refine
!==================================================================================================
!
! Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo()
!
! Purpose: Constructs the Vandermonde matrix from the Lagrange basis functions on the nodal set 
!          xione evaluated at the nodes of xitwo
! 
! Comments: The computation of the interpolation basis functions is based on 
!           "Barycentric Lagrange Interpolation", SIAM REview, Vol 46, No 3 pp. 501-517 
!           (see docs/barycentric.pdf). The notation in this subroutine is consistent with that 
!           in the paper and specifically with the section 3 "An improved Lagrange Formula".
!
!           the barycentric weights, wj, are computed as wj = 1/(Pi_{k=0^{n-1,K!=j}(x_{j}-x_{k})}
!           the Lagrange basis functions are constructed as l_{j} =l*w_{j}/(x-x_{j}) where 
!           l = Pi_{k=0}^{n-1}(x-x_{j})
! Additional documentation: docs/barycentric.pdf
!
! Unit tests: unit_tests/Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo_program 
!             compile using sh compile_Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo
!
! Inputs: nxione (int): number of nodes in the nodal set xione
!         nxitwo (int): number of nodes in the nodal set xitwo
!         xione (real(wp)) size (nxioneX1): nodal set that is used to construct the Lagrange basis functions
!         xitwo (real(wp)) size (nxitwoX1): nodal set at which the Lagrange basis functions will be evaluated at
!
! Outputs:
!         Vandermonde (real(wp)) size (nxitwoXnxione): the above mentioned Vandermonde matrix
!
!===================================================================================================
   subroutine Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo(nxione,nxitwo,xione,xitwo,Vandermonde)

     !-- Nothing is implicitely defined
     implicit none

     !-- input variables
     integer, intent(in)                         :: nxione, nxitwo
     real(wp), intent(in)                        :: xione(nxione), xitwo(nxitwo)
     real(wp), intent(inout)                     :: Vandermonde(nxitwo,nxione)
     
     !-- local variables
     integer                                     :: i, j, k
     real(wp)                                    :: wj(nxione),l
 
     !-- construct the bary centric weights, wj,
     do j = 0,nxione-1 
       wj(j+1) = 1.0_wp
       do k = 0,nxione-1
         if (k.NE.j) then
           wj(j+1) = wj(j+1)*(xione(j+1)-xione(k+1))
         end if
       end do
       wj(j+1) = 1.0_wp/wj(j+1)
     end do
     !-- construct the Vandermonde matrix V(j,k) = l_k(xitwo(j)), where l_k is the Lagrange 
     !-- interpolation basis construted from xione
     do j = 1,nxione
       do k = 1,nxitwo
         !-- construct l
         l = 1.0_wp
         do i = 1, nxione
           l = l*(xitwo(k)-xione(i))
         end do
         !-- check to see if the nodal location on the second set of nodes matches the node 
         !-- location associated with the kth Lagrange basis function 
         if (abs(xitwo(k)-xione(j)) <=2.0_wp*epsilon(1.0_wp)) then
           Vandermonde(k,j) = 1.0_wp
         else
           Vandermonde(k,j)  = l*wj(j)/(xitwo(k)-xione(j))
         end if
       end do
     end do

   end subroutine Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo

!==================================================================================================
!
! Vandermonde_1D_monomial()
!
! Purpose: Constructs the Vandermonde matrix from the monomials on the node set xi 
! 
! Comments: 
!
! Additional documentation:
!
! Unit tests: unit_tests/Vandermonde_1D_monomial.f90 
!             compile using sh unit_tests/comiplie_Vandermonde_1D_monomial.sh
!
! Inputs: n (int): number of nodes in the nodal set xi
!         p (int): the highest degree monomial that will be evaluated
!         xi (real(wp)) size (nX1): nodal on which the monomials will be evaluated
!
! Outputs:
!         Vandermonde (real(wp)) size (nX(p+1)): the above mentioned Vandermonde matrix
!
!===================================================================================================
   subroutine Vandermonde_1D_monomial(n,p,xi,Vandermonde)

     !-- Nothing is implicitely defined
     implicit none

     !-- input variables
     integer, intent(in)                         :: n, p
     real(wp), intent(in)                        :: xi(n)
     real(wp), intent(inout)                     :: Vandermonde(n,(p+1))
     
     !-- local variables
     integer                                     :: i, j
 
     !-- construct the Vandermonde matrix
     do i = 1,n
       do j = 1,p+1
         Vandermonde(i,j) = xi(i)**(j-1)
       end do
     end do
   end subroutine Vandermonde_1D_monomial

!==================================================================================================
!
! Rotate_xione_2_xitwo_and_back()
!
! Purpose: This constructs the interpolants from xione to xitwo and back.
!
! Comments: This is done by constructing the Vandermonde matrix of the Lagrange basis functions, 
!           constructed from xione, evaluated on xitwo, Ixione2xitwo. The interpolant back 
!           is given by Pxione^-1*Ixione2xitwo^T*Pxitwo, where Pxione and Pxitwo are the 
!           diagonal norm matrices on xione and xitwo, respectively. 
!
! Additional documentation
!
! Inputs: nxione (int): number of nodes in the one-dimensional nodal distribution xione
!         nxitwo (int): number of nodes in the one-dimensional nodal distribution xitwo
!         xione (real(wp)) size (nxioneX1): first nodal distribution
!         xitwo (real(wp)) size (nxitwoX1): second nodal distribution
!         Bxione (real(wp) size (nxioneX1): quadrature weights on the first nodal distribution
!         Bxitwo (real(wp)) size (nxitwoX1): quadrature weights on the second nodal distribution
!
! Outputs:
!         xione2xitwo (real(wp)) size(nxioneXnxitwo): interpolation matrix from xione to xitwo
!         xitwo2xione (real(wp)) size(nxitwoXnxione): interpolation matrix from xitwo to xione
!
!===================================================================================================
   subroutine Rotate_xione_2_xitwo_and_back(nxione,nxitwo,xione,xitwo,Bxione,Bxitwo,xione2xitwo,xitwo2xione)

     !-- Nothing is implicitely defined
     implicit none

     !-- input variables
     integer, intent(in)                         :: nxione,nxitwo
     real(wp), intent(in)                        :: xione(nxione),xitwo(nxitwo),Bxione(nxione),Bxitwo(nxitwo)
     real(wp), intent(inout)                     :: xione2xitwo(nxitwo,nxione), xitwo2xione(nxione,nxitwo)

     !-- local variables
     integer                                     :: i,j
     !-- construct the Vandermonde matrix from the Lagrange basis funcitons on xione evaluated at
     !-- the nodal locations of xitwo

     call Vandermonde_1D_Lagrange_on_XIone_eval_at_XItwo(nxione,nxitwo,xione,xitwo,xione2xitwo)
    
     !-- construct the interpolant from xitwo to xione
     do i = 1, nxione
       do j = 1,nxitwo
         xitwo2xione(i,j) = 1.0_wp/Bxione(i)*(xione2xitwo(j,i)*Bxitwo(j))
       end do 
     end do 
   end subroutine Rotate_xione_2_xitwo_and_back

!==================================================================================================
!
! Rotate_GL_2_G_and_back_I()
!
! Purpose: This gives back the interpolant from GL to G and back
!
! Comments: The interpolants were solved for in the Maple script maple/Construct_1D_interpolants_Identity.mw
!
! Additional documentation
!
! Inputs: nGL (int): number of nodes in the one-dimensional Gauss Lobatto (GL) nodal distribution
!         nG (int): number of nodes in the one-dimensional Gauss (G) nodal distribution
!         xGL (real(wp)) size (nxioneX1): GL nodal distribution
!         xG (real(wp)) size (nxitwoX1): G nodal distribution
!         BxiGL (real(wp) size (nxioneX1): quadrature weights on the GL nodal distribution
!         BxiG (real(wp)) size (nxitwoX1): quadrature weights on the G nodal distribution
!
! Outputs:
!         IGL2G (real(wp)) size(nGXnGL): interpolation matrix from the GL nodes to the G nodes
!         IG2GL (real(wp)) size(nGLXnG): interpolation matrix from the G nodes to GL nodes
!
!===================================================================================================
   subroutine Rotate_GL_2_G_and_back_I(nGL,nG,xiGL,xiG,BxiGL,BxiG,IGL2G,IG2GL)

     !-- Nothing is implicitely defined
     implicit none

     !-- input variables
     integer, intent(in)                         :: nGL,nG
     real(wp), intent(in)                        :: xiGL(nGL),xiG(nG),BxiGL(nGL),BxiG(nG)
     real(wp), intent(inout)                     :: IGL2G(nG,nGL), IG2GL(nGL,nG)

     !-- local variables
     

    if (nGL==2.AND.nG==2) then
      IGL2G(1,1) = 1.000000000000000000000000000000_wp
      IGL2G(1,2) = 0.000000000000000000000000000000_wp
      IGL2G(2,1) = 0.000000000000000000000000000000_wp
      IGL2G(2,2) = 1.000000000000000000000000000000_wp
      
      IG2GL(1,1) = 1.000000000000000000000000000000_wp
      IG2GL(1,2) = 0.000000000000000000000000000000_wp
      IG2GL(2,1) = 0.000000000000000000000000000000_wp
      IG2GL(2,2) = 1.000000000000000000000000000000_wp
    else if (nGL==3.AND.nG==3) then
      IGL2G(1,1) = 0.343146490609516399718000303683_wp
      IGL2G(1,2) = 1.088303688022450577599852472590_wp
      IGL2G(1,3) = -0.431450178631966977317852776273_wp
      IGL2G(2,1) = 0.430189805014031610999907795369_wp
      IGL2G(2,2) = 0.139620389971936778000184409261_wp
      IGL2G(2,3) = 0.430189805014031610999907795369_wp
      IGL2G(3,1) = -0.431450178631966977317852776273_wp
      IGL2G(3,2) = 1.088303688022450577599852472590_wp
      IGL2G(3,3) = 0.343146490609516399718000303683_wp
      
      IG2GL(1,1) = 0.571910817682527332863333839473_wp
      IG2GL(1,2) = 1.147172813370750962666420787650_wp
      IG2GL(1,3) = -0.719083631053278295529754627123_wp
      IG2GL(2,1) = 0.453459870009354407333271863580_wp
      IG2GL(2,2) = 0.093080259981291185333456272841_wp
      IG2GL(2,3) = 0.453459870009354407333271863580_wp
      IG2GL(3,1) = -0.719083631053278295529754627123_wp
      IG2GL(3,2) = 1.147172813370750962666420787650_wp
      IG2GL(3,3) = 0.571910817682527332863333839473_wp      
    else if (nGL==4.AND.nG==4) then 
      IGL2G(1,1) = 0.670133597061877309375642466273_wp
      IGL2G(1,2) = 0.382690211445752349043404473625_wp
      IGL2G(1,3) = -0.059634895378013858655525747055_wp
      IGL2G(1,4) = 0.006811086870384200236478807157_wp
      IGL2G(2,1) = -0.124994022375115610291702237717_wp
      IGL2G(2,2) = 1.094392949182378686886803453950_wp
      IGL2G(2,3) = 0.011123163321311394153889248044_wp
      IGL2G(2,4) = 0.019477909871425529251009535715_wp
      IGL2G(3,1) = 0.019477909871425529251009535715_wp
      IGL2G(3,2) = 0.011123163321311394153889248052_wp
      IGL2G(3,3) = 1.094392949182378686886803453960_wp
      IGL2G(3,4) = -0.124994022375115610291702237717_wp
      IGL2G(4,1) = 0.006811086870384200236478807157_wp
      IGL2G(4,2) = -0.059634895378013858655525747055_wp
      IGL2G(4,3) = 0.382690211445752349043404473625_wp
      IGL2G(4,4) = 0.670133597061877309375642466273_wp
      
      IG2GL(1,1) = 1.398655311764185408399936481650_wp
      IG2GL(1,2) = -0.489085476472273963755402275747_wp
      IG2GL(1,3) = 0.076214547296997107993785971074_wp
      IG2GL(1,4) = 0.014215617411091447361679823021_wp
      IG2GL(2,1) = 0.159744773085697986032413787639_wp
      IG2GL(2,2) = 0.856443671190025127617787284033_wp
      IG2GL(2,3) = 0.008704700480085614686859237835_wp
      IG2GL(2,4) = -0.024893144755808728337060309506_wp
      IG2GL(3,1) = -0.024893144755808728337060309506_wp
      IG2GL(3,2) = 0.008704700480085614686859237829_wp
      IG2GL(3,3) = 0.856443671190025127617787284040_wp
      IG2GL(3,4) = 0.159744773085697986032413787639_wp
      IG2GL(4,1) = 0.014215617411091447361679823021_wp
      IG2GL(4,2) = 0.076214547296997107993785971074_wp
      IG2GL(4,3) = -0.489085476472273963755402275747_wp
      IG2GL(4,4) = 1.398655311764185408399936481650_wp
    else if (nGL==5.AND.nG==5) then
      IGL2G(1,1) = 0.470502966472490507460316236133_wp
      IGL2G(1,2) = 0.803126567593681665261854299464_wp
      IGL2G(1,3) = -0.491470918189824019737687496142_wp
      IGL2G(1,4) = 0.369914191164023953735669915576_wp
      IGL2G(1,5) = -0.152072807040372106720152955027_wp
      IGL2G(2,1) = 0.071767876791755725608492888402_wp
      IGL2G(2,2) = 0.529448395544771508275461087673_wp
      IGL2G(2,3) = 0.689001782387354883935218360338_wp
      IGL2G(2,4) = -0.492612611092600584063108759494_wp
      IGL2G(2,5) = 0.202394556368718466243936423086_wp
      IGL2G(3,1) = -0.187499999999999999999999999998_wp
      IGL2G(3,2) = 0.437500000000000000000000000001_wp
      IGL2G(3,3) = 0.500000000000000000000000000006_wp
      IGL2G(3,4) = 0.437500000000000000000000000001_wp
      IGL2G(3,5) = -0.187499999999999999999999999998_wp
      IGL2G(4,1) = 0.202394556368718466243936423086_wp
      IGL2G(4,2) = -0.492612611092600584063108759494_wp
      IGL2G(4,3) = 0.689001782387354883935218360338_wp
      IGL2G(4,4) = 0.529448395544771508275461087673_wp
      IGL2G(4,5) = 0.071767876791755725608492888402_wp
      IGL2G(5,1) = -0.152072807040372106720152955027_wp
      IGL2G(5,2) = 0.369914191164023953735669915576_wp
      IGL2G(5,3) = -0.491470918189824019737687496142_wp
      IGL2G(5,4) = 0.803126567593681665261854299464_wp
      IGL2G(5,5) = 0.470502966472490507460316236133_wp
      
      IG2GL(1,1) = 1.114748022560237464834855958820_wp
      IG2GL(1,2) = 0.343501634534003810651995244455_wp
      IG2GL(1,3) = -1.066666666666666666666666666660_wp
      IG2GL(1,4) = 0.968718374310688038674062097476_wp
      IG2GL(1,5) = -0.360301364738262647494246634073_wp
      IG2GL(2,1) = 0.349498057896440617267939081403_wp
      IG2GL(2,2) = 0.465445435697663304036304550688_wp
      IG2GL(2,3) = 0.457142857142857142857142857140_wp
      IG2GL(2,4) = -0.433062586135970603458639415191_wp
      IG2GL(2,5) = 0.160976235399009539297252925968_wp
      IG2GL(3,1) = -0.163747509950278330491814509669_wp
      IG2GL(3,2) = 0.463747509950278330491814509666_wp
      IG2GL(3,3) = 0.400000000000000000000000000005_wp
      IG2GL(3,4) = 0.463747509950278330491814509666_wp
      IG2GL(3,5) = -0.163747509950278330491814509669_wp
      IG2GL(4,1) = 0.160976235399009539297252925968_wp
      IG2GL(4,2) = -0.433062586135970603458639415191_wp
      IG2GL(4,3) = 0.457142857142857142857142857140_wp
      IG2GL(4,4) = 0.465445435697663304036304550688_wp
      IG2GL(4,5) = 0.349498057896440617267939081403_wp
      IG2GL(5,1) = -0.360301364738262647494246634073_wp
      IG2GL(5,2) = 0.968718374310688038674062097476_wp
      IG2GL(5,3) = -1.066666666666666666666666666660_wp
      IG2GL(5,4) = 0.343501634534003810651995244455_wp
      IG2GL(5,5) = 1.114748022560237464834855958820_wp     
    else if (nGL==6.AND.nG==6) then

    else if (nGL==7.AND.nG==7) then
    else if (nGL==8.AND.nG==8) then
    else if (nGL==9.AND.nG==9) then
    else if (nGL==10.AND.nG==10) then
    else if (nGL==11.AND.nG==11) then
    else if (nGL==12.AND.nG==12) then
    else if (nGL==13.AND.nG==13) then
    else if (nGL==15.AND.nG==15) then
    else if (nGL==16.AND.nG==16) then
    else if (nGL==17.AND.nG==17) then
    end if 
   end subroutine Rotate_GL_2_G_and_back_I

!==================================================================================================
!
! Lagrange_interpolant_basis_1D()
!
! Purpose: Evaluates the jth Lagrange basis function, constructed from the nodes xivec at the 
!          point xi.
!
! Additional documentation
!
! Inputs: 
!       n (int): number of nodes in the one-dimensional nodal distribution xivec
!       xivec (real(wp) vector of size [nX1]): one-dimensional nodal distribution
!       xi_eval (real(wp)): location at which the jth Lagrange basis function is evaluated
!       j (int): which Lagrange basis (numbered 0 through n-1)
!
! Outputs:
! 
!        lj_eval (real(wp)): evaluation of the jth Lagrange basis at xi_eval
!
!===================================================================================================
   subroutine Lagrange_interpolant_basis_1D(n,xivec,xi_eval,j,lj_eval)

     !-- arguments
     integer, intent(in)                         :: n,j
     real(wp), intent(in)                        :: xivec(1:n)
     real(wp), intent(in)                        :: xi_eval
     real(wp), intent(inout)                     :: lj_eval
    
     !-- local variables
     integer                                     :: k
     real(wp)                                    :: wj,l,xij

     !-- construct wj
     wj = 1.0_wp
     xij = xivec(j+1)
     do k = 0,n-1
       if (k.NE.j) then
         wj = wj*(xij-xivec(k+1))
       end if
     end do
     wj = 1.0_wp/wj

     !-- construct l
     l = 0.0_wp
     do k = 0,n-1
       l = l*(xi_eval)
     end do
     lj_eval = 1.0_wp
   end subroutine Lagrange_interpolant_basis_1D

!==================================================================================================
!
! Derivative_Lagrange_interpolant_basis_1D()
!
! Purpose:
!
! Additional documentation
!
! Inputs:
!
! Outputs:
!
!===================================================================================================
   !subroutine Derivative_Lagrange_interpolant_basis_1D()
   !end subroutine Derivative_Lagrange_interpolant_basis_1D
end module non_conforming
