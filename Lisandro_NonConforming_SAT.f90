     subroutine ldgmortarp( &
     app,              &
     nqm,R_m,          &
     nqp,R_p,          &
     facetable_m,      &
     facetable_p,      &
     P_m,dRdX_m,       &
     P_p,dRdX_p,       &
     id_m,U_m,         &
     id_p,U_p,         &
     G)                &
     bind(C)
  type(AppType),    intent(in)    :: app
  integer(kind=ik), intent(in),   value :: nqm
  real(kind=rk),    intent(in)    :: R_m(nqm)
  integer(kind=ik), intent(in),   value :: nqp
  real(kind=rk),    intent(in)    :: R_p(nqp)
  integer(kind=ik), parameter     :: npr = max(1,dim-1) * 2**(dim-1)
  integer(kind=ik), intent(in)    :: facetable_m(nqm**(dim-1),npr,2,dim)
  integer(kind=ik), intent(in)    :: facetable_p(nqp**(dim-1),npr,2,dim)
  real(kind=rk),    intent(in)    :: P_m(nqm)
  real(kind=rk),    intent(in)    :: dRdX_m(dim,dim,nqm**dim)
  real(kind=rk),    intent(in)    :: P_p(nqp)
  real(kind=rk),    intent(in)    :: dRdX_p(dim,dim,nqp**dim)
  integer(kind=ik), intent(in)    :: id_m(3)
  real(kind=rk),    intent(in)    :: U_m(dof,nqm**dim)
  integer(kind=ik), intent(in)    :: id_p(3)
  real(kind=rk),    intent(in)    :: U_p(dof,nqp**dim)
  real(kind=rk),    intent(inout) :: G(dim,dof,nqm**dim)

  integer(kind=ik) :: i,m
  integer(kind=ik) :: maxis,mside,mperm
  integer(kind=ik) :: paxis,pside,pperm
  real(kind=rk)    :: Pinv,sn,scale,nrm(dim)
  real(kind=rk)    :: U_i(dof,nen**(dim-1))
  real(kind=rk)    :: ldg(dim,dof)

  procedure(conv_fn), pointer :: entropy
  procedure(conv_fn), pointer :: primary
  call c_f_procpointer(app%entropy,entropy)
  call c_f_procpointer(app%primary,primary)

  maxis = id_m(1)+1; mside = id_m(2)+1; mperm = id_m(3)+1;
  paxis = id_p(1)+1; pside = id_p(2)+1; pperm = id_p(3)+1;
  Pinv = 1/P_m(1); sn = real(sign(1_ik,mside-2),rk)
  scale = 0.5_rk * Pinv * sn

  call interpolate()

  do i=1,nen**(dim-1)
     m = facetable_m(i,mperm,mside,maxis)+1

     nrm = scale * dRdX_m(:,maxis,m)
     call jumpn(U_m(:,m),U_i(:,i),nrm,ldg)

     G(:,:,m) = G(:,:,m) - ldg
  end do

contains
#include "nrmjump.f90.in"
#include "pinterp.f90.in"
#if DIM > 1
  subroutine interpolate()
    real(kind=rk) :: B(nen,0:nqp-1)
    real(kind=rk) :: W_p(dof,nqp**(dim-1))
    real(kind=rk) :: W_i(dof,nen**(dim-1))
    integer(kind=ik) :: i,tb(nqp**(dim-1))
    call lagrange(nqp-1,R_p,R_m,B)
    tb = facetable_p(:,pperm,pside,paxis)+1
    if (associated(entropy)) then
       do i=1,nqp**(dim-1)
          call entropy(app%ctx,U_p(:,tb(i)),W_p(:,i))
       end do
       call tensor(nqp-1,B,W_p,W_i)
       do i=1,nen**(dim-1)
          call primary(app%ctx,W_i(:,i),U_i(:,i))
       end do
    else
       W_p = U_p(:,tb)
       call tensor(nqp-1,B,W_p,U_i)
    endif
  end subroutine interpolate
#else
  subroutine interpolate
    integer(kind=ik) :: tb(nqp**(dim-1))
    tb = facetable_p(:,pperm,pside,paxis)+1
    U_i = U_p(:,tb)
  end subroutine interpolate
#endif
end subroutine ldgmortarp


subroutine satmortarp(  &
     app,               &
     nqm,R_m,           &
     nqp,R_p,           &
     facetable_m,       &
     facetable_p,       &
     P_m,detJ_m,Jnrm_m, &
     P_p,detJ_p,Jnrm_p, &
     id_m,U_m,G_m,      &
     id_p,U_p,G_p,      &
     F)                 &
     bind(C)
  type(AppType),    intent(in)    :: app
  integer(kind=ik), intent(in),   value :: nqm
  real(kind=rk),    intent(in)    :: R_m(nqm)
  integer(kind=ik), intent(in),   value :: nqp
  real(kind=rk),    intent(in)    :: R_p(nqp)
  integer(kind=ik), parameter     :: npr = max(1,dim-1) * 2**(dim-1)
  integer(kind=ik), intent(in)    :: facetable_m(nqm**(dim-1),npr,2,dim)
  integer(kind=ik), intent(in)    :: facetable_p(nqp**(dim-1),npr,2,dim)
  real(kind=rk),    intent(in)    :: P_m(nqm)
  real(kind=rk),    intent(in)    :: detJ_m(nqm**dim)
  real(kind=rk),    intent(in)    :: Jnrm_m(dim,dim,nqm**dim)
  real(kind=rk),    intent(in)    :: P_p(nqp)
  real(kind=rk),    intent(in)    :: detJ_p(nqp**dim)
  real(kind=rk),    intent(in)    :: Jnrm_p(dim,dim,nqp**dim)
  integer(kind=ik), intent(in)    :: id_m(3)
  real(kind=rk),    intent(in)    :: U_m(dof,nqm**dim)
  real(kind=rk),    intent(in)    :: G_m(dim,dof,nqm**dim)
  integer(kind=ik), intent(in)    :: id_p(3)
  real(kind=rk),    intent(in)    :: U_p(dof,nqp**dim)
  real(kind=rk),    intent(in)    :: G_p(dim,dof,nqp**dim)
  real(kind=rk),    intent(inout) :: F(dof,nqm**dim)

  integer(kind=ik) :: i,j,m,p
  integer(kind=ik) :: maxis,mside,mperm
  integer(kind=ik) :: paxis,pside,pperm
  real(kind=rk)    :: sn_m,Jn_m(dim)
  real(kind=rk)    :: sn_p,Jn_p(dim)
  real(kind=rk)    :: Pinv,Jn(dim)
  real(kind=rk)    :: Um(dof),Up(dof)
  !real(kind=rk)    :: nrm(dim),dJ,dWn(dim,dof)
  real(kind=rk)    :: GWm(dim,dof),GWp(dim,dof)
  real(kind=rk)    :: fIm(dof),fIp(dof),fIs(dof)
  real(kind=rk)    :: fVm(dof),fVp(dof),fVs(dof)
  real(kind=rk)    :: sat(dof)

  integer(kind=ik) :: tb_m(nqm**(dim-1))
  integer(kind=ik) :: tb_p(nqp**(dim-1))
  real(kind=rk)    :: Q(nqm,nqp)

  procedure(conv_fn),    pointer :: entropy
  procedure(iflux_fn),   pointer :: iflux
  procedure(iflux1n_fn), pointer :: iflux1n
  procedure(iflux2n_fn), pointer :: iflux2n
  procedure(upwind_fn),  pointer :: upwind
  procedure(vflux_fn),   pointer :: vflux
  call c_f_procpointer(app%entropy,entropy)
  call c_f_procpointer(app%iflux,iflux)
  call c_f_procpointer(app%iflux1n,iflux1n)
  call c_f_procpointer(app%iflux2n,iflux2n)
  call c_f_procpointer(app%upwind,upwind)
  call c_f_procpointer(app%vflux,vflux)

  maxis = id_m(1)+1; mside = id_m(2)+1; mperm = id_m(3)+1;
  paxis = id_p(1)+1; pside = id_p(2)+1; pperm = id_p(3)+1;
  tb_m = facetable_m(:,mperm,mside,maxis)+1
  tb_p = facetable_p(:,pperm,pside,paxis)+1
  sn_m = +real(sign(1_ik,mside-2),rk)
  sn_p = -real(sign(1_ik,pside-2),rk)
  Pinv = 1/P_m(1)

  call transferop(nqm,R_m,P_m,nqp,R_p,P_p,Q)

  if (associated(iflux2n)) then
     do i=1,nqm**(dim-1)
        sat = 0
        m = tb_m(i)
        Um = U_m(:,m)
        Jn_m = sn_m * Jnrm_m(:,maxis,m)
        call iflux1n(app%ctx,Um,Jn_m,fIm)
        sat = sat + fIm
        do j=1,nqp**(dim-1)
           p = tb_p(j)
           Up = U_p(:,p)
           Jn_p = sn_p * Jnrm_p(:,paxis,p)
           Jn = 0.5_rk * (Jn_m + Jn_p)
           call iflux2n(app%ctx,Um,Up,Jn,fIs)
           sat = sat - P2M(i,j) * fIs
        end do
        F(:,m) = F(:,m) + Pinv * sat
     end do
  else if (associated(iflux)) then
     do i=1,nqm**(dim-1)
        sat = 0
        m = tb_m(i)
        Um = U_m(:,m)
        Jn_m = sn_m * Jnrm_m(:,maxis,m)
        call ifluxn(app%ctx,Um,Jn_m,fIm)
        sat = sat + fIm
        do j=1,nqp**(dim-1)
           p = tb_p(j)
           Up = U_p(:,p)
           Jn_p = sn_p * Jnrm_p(:,paxis,p)
           Jn = 0.5_rk * (Jn_m + Jn_p)
           call ifluxn(app%ctx,Um,Jn,fIm)
           call ifluxn(app%ctx,Up,Jn,fIp)
           fIs = 0.5_rk * (fIm + fIp)
           sat = sat - P2M(i,j) * fIs
        end do
        F(:,m) = F(:,m) + Pinv * sat
     end do
  end if

  if (associated(upwind)) then
     do i=1,nqm**(dim-1)
        sat = 0
        m = tb_m(i)
        Um = U_m(:,m)
        Jn_m = sn_m * Jnrm_m(:,maxis,m)
        do j=1,nqp**(dim-1)
           p = tb_p(j)
           Up = U_p(:,p)
           Jn_p = sn_p * Jnrm_p(:,paxis,p)
           Jn = 0.5_rk * (Jn_m + Jn_p)
           call upwind(app%ctx,Um,Up,Jn,fIs)
           sat = sat - P2M(i,j) * fIs
        end do
        F(:,m) = F(:,m) + Pinv * sat
     end do
  end if

  if (associated(vflux)) then
     do i=1,nqm**(dim-1)
        sat = 0
        m = tb_m(i)
        Um = U_m(:,m)
        GWm = G_m(:,:,m) ! XXX
        Jn_m = sn_m * Jnrm_m(:,maxis,m)
        call vfluxn(app%ctx,Um,GWm,Jn_m,fVm)
        sat = sat + fVm
        do j=1,nqp**(dim-1)
           p = tb_p(j)
           Up = U_p(:,p)
           GWp = G_p(:,:,p) ! XXX
           Jn_p = sn_p * Jnrm_p(:,paxis,p)
           Jn = 0.5_rk * (Jn_m + Jn_p)
           call vfluxn(app%ctx,Um,GWm,Jn,fVm)
           call vfluxn(app%ctx,Up,GWp,Jn,fVp)
           fVs = 0.5_rk * (fVm + fVp)
           sat = sat - P2M(i,j) * fIs
        end do
        F(:,m) = F(:,m) + Pinv * sat
     end do
  endif

  !do i=1,nen**(dim-1)
  !   m = facetable_m(i,mperm,mside,maxis)+1
  !
  !   Jn = sn * Jnrm_m(:,maxis,m)
  !   dJ = detJ_m(m)
  !   Um = U_m(:,m)
  !   Up = U_i(:,i)
  !
  !   sat = 0
  !
  !   if (associated(iflux2n)) then
  !      call iflux1n(app%ctx,Um,Jn,fIm)
  !      call iflux2n(app%ctx,Um,Up,Jn,fIs)
  !      sat = sat + (fIm - fIs)
  !   else if (associated(iflux)) then
  !      call ifluxn(app%ctx,Um,Jn,fIm)
  !      call ifluxn(app%ctx,Up,Jn,fIp)
  !      sat = sat + 0.5_rk * (fIm - fIp)
  !   end if
  !
  !   if (associated(upwind)) then
  !      call upwind(app%ctx,Um,Up,Jn,vec)
  !      sat = sat - vec
  !   end if
  !
  !   if (associated(vflux)) then
  !      nrm = 0.5_rk * Pinv/dJ * Jn
  !      call jumpn(Um,Up,nrm,dWn)
  !      GWm = G_m(:,:,m) + dWn
  !      GWp = G_i(:,:,i) - dWn
  !      call vfluxn(app%ctx,Um,GWm,Jn,fVm)
  !      call vfluxn(app%ctx,Up,GWp,Jn,fVp)
  !      sat = sat - 0.5_rk * (fVm - fVp)
  !   end if
  !
  !   F(:,m) = F(:,m) + Pinv * sat
  !end do

contains
#include "nrmflux.f90.in"

!dir$ attributes forceinline :: transferop
pure subroutine transferop(m,Rm,Pm, &
                           p,Rp,Pp, &
                           Q)
  integer(kind=ik), intent(in)  :: m,p
  real(kind=rk),    intent(in)  :: Rm(m),Pm(m)
  real(kind=rk),    intent(in)  :: Rp(p),Pp(p)
  real(kind=rk),    intent(out) :: Q(m*p)
  integer(kind=ik)  :: k,j,i
  real(kind=rk)     :: L

  if (m < p) then

     do k=1,p
        do j=1,m
           L = Pp(k)/Pm(j)
           do i=1,m
              if (i == j) cycle
              L = L * (Rp(k)-Rm(i))/(Rm(j)-Rm(i))
           end do
           Q(j+(k-1)*m) = L
        end do
     end do

  else

     do k=1,m
        do j=1,p
           L = 1
           do i=1,p
              if (i == j) cycle
              L = L * (Rm(k)-Rp(i))/(Rp(j)-Rp(i))
           end do
           Q(k+(j-1)*m) = L
        end do
     end do

  endif

end subroutine transferop

#if DIM == 3

pure real(kind=rk) function P2M(i,j)
  integer(kind=ik), intent(in) :: i,j
  integer(kind=ik) :: i1,i2,j1,j2
  i1 = mod(i-1,nqm)+1; i2 = (i-1)/nqm+1;
  j1 = mod(j-1,nqp)+1; j2 = (j-1)/nqp+1;
  P2M = Q(i1,j1) * Q(i2,j2)
end function P2M

#elif DIM ==2

pure real(kind=rk) function P2M(i,j)
  integer(kind=ik), intent(in) :: i,j
  P2M = Q(i,j)
end function P2M

#else

pure real(kind=rk) function P2M(i,j)
  integer(kind=ik), intent(in) :: i,j
  P2M = 1
end function P2M

#endif

end subroutine satmortarp
