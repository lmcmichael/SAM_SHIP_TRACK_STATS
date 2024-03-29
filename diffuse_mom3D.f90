
subroutine diffuse_mom3D
	
!        momentum tendency due to SGS diffusion

use vars
use sgs, only: tk, grdf_x, grdf_y, grdf_z
use params, only: docolumn, dowallx, dowally, doMomentumHyperviscosity, &
                  tau_MomentumHyperviscosity
implicit none

real rdx2,rdy2,rdz2,rdz,rdx25,rdy25
real rdx21,rdy21,rdx251,rdy251,rdz25
real dxy,dxz,dyx,dyz,dzx,dzy

integer i,j,k,ic,ib,jb,jc,kc,kcu
real tkx, tky, tkz, rhoi, iadzw, iadz
real fu(0:nx,0:ny,nz),fv(0:nx,0:ny,nz),fw(0:nx,0:ny,nz)

!hyperviscosity 
real rdx16, rdy16, khyp

rdx2=1./(dx*dx)
rdy2=1./(dy*dy)

rdx25=0.25*rdx2
rdy25=0.25*rdy2

dxy=dx/dy
dxz=dx/dz
dyx=dy/dx
dyz=dy/dz

!hyperviscosity
rdx16 = rdx25*rdx25
rdy16 = rdy25*rdy25

!bloss: Set default timescale for damping 2-delta waves in x to 60 seconds.
khyp = SQRT(dx*dy)**4 / tau_MomentumHyperviscosity ! Units m^4/s

!-----------------------------------------
if(dowallx) then

  if(mod(rank,nsubdomains_x).eq.0) then
    do k=1,nzm
     do j=1,ny
         v(0,j,k) = v(1,j,k)
         w(0,j,k) = w(1,j,k)
     end do
    end do
  end if
  if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
    do k=1,nzm
     do j=1,ny
         v(nx+1,j,k) = v(nx,j,k)
         w(nx+1,j,k) = w(nx,j,k)
     end do
    end do
  end if

end if

if(dowally) then

  if(rank.lt.nsubdomains_x) then
    do k=1,nzm
       do i=1,nx
         u(i,1-YES3D,k) = u(i,1,k)
         w(i,1-YES3D,k) = w(i,1,k)
       end do
    end do
  end if
  if(rank.gt.nsubdomains-nsubdomains_x-1) then
    do k=1,nzm
       do i=1,nx
         u(i,ny+YES3D,k) = u(i,ny,k)
         w(i,ny+YES3D,k) = w(i,ny,k)
       end do
    end do
  end if

end if


do k=1,nzm
 kc=k+1
 kcu=min(kc,nzm)
 dxz=dx/(dz*adzw(kc))
 dyz=dy/(dz*adzw(kc))
  rdx21=rdx2    * grdf_x(k)
  rdy21=rdy2    * grdf_y(k)
  rdx251=rdx25  * grdf_x(k)
  rdy251=rdy25  * grdf_y(k)
  do j=1,ny
   jb=j-1
   do i=0,nx
    ic=i+1
    tkx=rdx21*tk(i,j,k)
    fu(i,j,k)=-2.*tkx*(u(ic,j,k)-u(i,j,k))
    tkx=rdx251*(tk(i,j,k)+tk(i,jb,k)+tk(ic,j,k)+tk(ic,jb,k))
    fv(i,j,k)=-tkx*(v(ic,j,k)-v(i,j,k)+(u(ic,j,k)-u(ic,jb,k))*dxy)
    tkx=rdx251*(tk(i,j,k)+tk(ic,j,k)+tk(i,j,kcu)+tk(ic,j,kcu)) 	
    fw(i,j,k)=-tkx*(w(ic,j,kc)-w(i,j,kc)+(u(ic,j,kcu)-u(ic,j,k))*dxz)
   end do 
   do i=1,nx
    ib=i-1
    dudt(i,j,k,na)=dudt(i,j,k,na)-(fu(i,j,k)-fu(ib,j,k))
    dvdt(i,j,k,na)=dvdt(i,j,k,na)-(fv(i,j,k)-fv(ib,j,k))
    dwdt(i,j,kc,na)=dwdt(i,j,kc,na)-(fw(i,j,k)-fw(ib,j,k))
   end do  
  end do 

  do j=0,ny
   jc=j+1
   do i=1,nx
    ib=i-1
    tky=rdy21*tk(i,j,k)
    fv(i,j,k)=-2.*tky*(v(i,jc,k)-v(i,j,k))
    tky=rdy251*(tk(i,j,k)+tk(ib,j,k)+tk(i,jc,k)+tk(ib,jc,k))
    fu(i,j,k)=-tky*(u(i,jc,k)-u(i,j,k)+(v(i,jc,k)-v(ib,jc,k))*dyx)
    tky=rdy251*(tk(i,j,k)+tk(i,jc,k)+tk(i,j,kcu)+tk(i,jc,kcu)) 	
    fw(i,j,k)=-tky*(w(i,jc,kc)-w(i,j,kc)+(v(i,jc,kcu)-v(i,jc,k))*dyz)
   end do 
  end do 
  do j=1,ny
    jb=j-1
    do i=1,nx	    
     dudt(i,j,k,na)=dudt(i,j,k,na)-(fu(i,j,k)-fu(i,jb,k))
     dvdt(i,j,k,na)=dvdt(i,j,k,na)-(fv(i,j,k)-fv(i,jb,k))
     dwdt(i,j,kc,na)=dwdt(i,j,kc,na)-(fw(i,j,k)-fw(i,jb,k))
   end do 
  end do 

end do 
 
!-------------------------
rdz=1./dz
dzx=dz/dx
dzy=dz/dy

do k=1,nzm-1
 kc=k+1
 uwsb(kc)=0.
 vwsb(kc)=0.
 iadz = 1./adz(k)
 iadzw= 1./adzw(kc)
 rdz2 = rdz*rdz * grdf_z(k)
 rdz25 = 0.25*rdz2
  do j=1,ny
   jb=j-1
   do i=1,nx
    ib=i-1
    tkz=rdz2*tk(i,j,k)
    fw(i,j,kc)=-2.*tkz*(w(i,j,kc)-w(i,j,k))*rho(k)*iadz
    tkz=rdz25*(tk(i,j,k)+tk(ib,j,k)+tk(i,j,kc)+tk(ib,j,kc))
    fu(i,j,kc)=-tkz*( (u(i,j,kc)-u(i,j,k))*iadzw + &
                       (w(i,j,kc)-w(ib,j,kc))*dzx)*rhow(kc) 	
    tkz=rdz25*(tk(i,j,k)+tk(i,jb,k)+tk(i,j,kc)+tk(i,jb,kc))
    fv(i,j,kc)=-tkz*( (v(i,j,kc)-v(i,j,k))*iadzw + &
                       (w(i,j,kc)-w(i,jb,kc))*dzy)*rhow(kc)
    uwsb(kc)=uwsb(kc)+fu(i,j,kc)
    vwsb(kc)=vwsb(kc)+fv(i,j,kc)
  end do 
 end do
end do

uwsb(1) = 0.
vwsb(1) = 0.
	
do j=1,ny
 do i=1,nx
   tkz=rdz2*grdf_z(nzm)*tk(i,j,nzm)
   fw(i,j,nz)=-2.*tkz*(w(i,j,nz)-w(i,j,nzm))/adz(nzm)*rho(nzm)
   fu(i,j,1)=fluxbu(i,j) * rdz * rhow(1)
   fv(i,j,1)=fluxbv(i,j) * rdz * rhow(1)
   fu(i,j,nz)=fluxtu(i,j) * rdz * rhow(nz)
   fv(i,j,nz)=fluxtv(i,j) * rdz * rhow(nz)
   uwsb(1) = uwsb(1) + fu(i,j,1)
   vwsb(1) = vwsb(1) + fv(i,j,1)
  end do
 end do
 
 do k=1,nzm
  kc=k+1
  rhoi = 1./(rho(k)*adz(k))
  do j=1,ny	  
   do i=1,nx
    dudt(i,j,k,na)=dudt(i,j,k,na)-(fu(i,j,kc)-fu(i,j,k))*rhoi
    dvdt(i,j,k,na)=dvdt(i,j,k,na)-(fv(i,j,kc)-fv(i,j,k))*rhoi
   end do
  end do
 end do ! k

 do k=2,nzm
  rhoi = 1./(rhow(k)*adzw(k))
  do j=1,ny
   do i=1,nx	 
    dwdt(i,j,k,na)=dwdt(i,j,k,na)-(fw(i,j,k+1)-fw(i,j,k))*rhoi
   end do
  end do
 end do ! k


 if(doMomentumHyperviscosity) then
   do k=1,nzm
      do j=1,ny
         do i=1,nx
            dudt(i,j,k,na) =  dudt(i,j,k,na) - khyp * (rdx16 *  &
                 (u(i-2,j,k) - 4*u(i-1,j,k) + 6*u(i,j,k) - 4*u(i+1,j,k) + u(i+2,j,k)) + &
                 rdy16* &
                 (u(i,j-2,k) - 4*u(i,j-1,k) + 6*u(i,j,k) - 4*u(i,j+1,k) + u(i,j+2,k)))

            dvdt(i,j,k,na) =  dvdt(i,j,k,na) - khyp * (rdx16 *  &
                 (v(i-2,j,k) - 4*v(i-1,j,k) + 6*v(i,j,k) - 4*v(i+1,j,k) + v(i+2,j,k)) + &
                 rdy16* &
                 (v(i,j-2,k) - 4*v(i,j-1,k) + 6*v(i,j,k) - 4*v(i,j+1,k) + v(i,j+2,k)))
         end do
      end do
   end do

   do k=2,nzm
      do j=1,ny
         do i=1,nx
            dwdt(i,j,k,na) =  dwdt(i,j,k,na) - khyp * (rdx16 *  &
                 (w(i-2, j,k) - 4*w(i-1,j,k) + 6*w(i,j,k) - 4*w(i+1,j,k) + w(i+2,j,k)) + &
                 rdy16* &
                 (w(i, j-2,k) - 4*w(i,j-1,k) + 6*w(i,j,k) - 4*w(i,j+1,k) + w(i,j+2,k)))
         end do
      end do
   end do
 end if

end subroutine diffuse_mom3D
