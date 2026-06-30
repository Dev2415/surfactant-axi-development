c***********************************************************************     
      subroutine cgstab(resmax,maxit)
c***********************************************************************     
      include 'mac3d.inc'
      dimension res(0:maxnx,0:maxny),reso(0:maxnx,0:maxny)
      dimension pk(0:maxnx,0:maxny),rk(0:maxnx,0:maxny)
      dimension zk(0:maxnx,0:maxny),vk(0:maxnx,0:maxny)
      dimension sk(0:maxnx,0:maxny),qk(0:maxnx,0:maxny)
      dimension ae(0:maxnx,0:maxny),aw(0:maxnx,0:maxny)
      dimension as(0:maxnx,0:maxny),an(0:maxnx,0:maxny)
      dimension ap(0:maxnx,0:maxny),d(0:maxnx,0:maxny)
      dimension xk(0:maxnx,0:maxny),Src(0:maxnx,0:maxny)

      eps=1.0d-30

      res=0.0d0
      reso=0.0d0
      pk=0.0d0
      rk=0.0d0
      zk=0.0d0
      vk=0.0d0
      sk=0.0d0
      qk=0.0d0
      xk=0.0d0
      ae=0.0d0
      aw=0.0d0
      as=0.0d0
      an=0.0d0
      ap=0.0d0
      d=0.0d0

      do j=2,nyp1
      do i=2,nxp1
         if(i.eq.2)then
            aw(i,j)=0.0d0
         else
            aw(i,j)=two*hxi*hxi/(r(i,j)+r(i-1,j))*x(i-1)/xs(i)
         endif
         if(i.eq.nxp1)then
            ae(i,j)=0.0d0
         else
            ae(i,j)=two*hxi*hxi/(r(i+1,j)+r(i,j))*x(i)/xs(i)
         endif
         if(j.eq.2)then
            as(i,j)=0.0d0
         else
            as(i,j)=two*hyi*hyi/(r(i,j)+r(i,j-1))
         endif
         an(i,j)=two*hyi*hyi/(r(i,j+1)+r(i,j))
         ap(i,j)=-ae(i,j)-aw(i,j)-an(i,j)-as(i,j)

         Sij=(x(i)*ut(i,j)-x(i-1)*ut(i-1,j))*hxi/xs(i)
     &                        +(vt(i,j)-vt(i,j-1))*hyi

         res(i,j)=Sij/dt
     &                  -( ae(i,j)*ps(i+1,j)+aw(i,j)*ps(i-1,j)
     &                    +an(i,j)*ps(i,j+1)+as(i,j)*ps(i,j-1)
     &                    +ap(i,j)*ps(i,j) )          
         reso(i,j)=res(i,j)
      end do
      end do

      do i=2,nxp1
      do j=2,nyp1
         d(i,j)=1.0d0/(ap(i,j)
     &              -aw(i,j)*d(i-1,j)*(ae(i-1,j)+0.995d0*an(i-1,j)) 
     &              -as(i,j)*d(i,j-1)*(an(i,j-1)+0.995d0*ae(i,j-1)))
      end do 
      end do

      do l=1,maxit

         rho_o=rho
         rho=0.0d0
         do i=2,nxp1
         do j=2,nyp1
            rho=rho+(res(i,j)+eps)*(reso(i,j)+eps)
         end do
         end do

         if(l.eq.1)then
            do i=2,nxp1
            do j=2,nyp1
                pk(i,j)=res(i,j)
            end do
            end do
         else
            beta=(rho+eps)/(rho_o+eps)*(alpha+eps)/(omega+eps)
            do i=2,nxp1
            do j=2,nyp1
                pk(i,j)=res(i,j)
     &             +(beta+eps)*(pk(i,j)+eps-(omega+eps)*(vk(i,j)+eps))
            end do
            end do
         endif

         do i=2,nxp1
         do j=2,nyp1
            zk(i,j)=( pk(i,j)-aw(i,j)*zk(i-1,j) 
     &                       -as(i,j)*zk(i,j-1)  )*d(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,nyp1
            zk(i,j)=(zk(i,j)+eps)/(d(i,j)+eps)
         end do
         end do

         do i=nxp1,2,-1
         do j=nyp1,2,-1
            zk(i,j)=( zk(i,j)-ae(i,j)*zk(i+1,j)
     &                       -an(i,j)*zk(i,j+1)  )*d(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,nyp1
            vk(i,j)=( ae(i,j)*zk(i+1,j)+aw(i,j)*zk(i-1,j)
     &               +an(i,j)*zk(i,j+1)+as(i,j)*zk(i,j-1)
     &               +ap(i,j)*zk(i,j) )          
         end do
         end do
 
         rtvk=zero
         do i=2,nxp1
         do j=2,nyp1
            rtvk=rtvk+(reso(i,j)+eps)*(vk(i,j)+eps)
         end do
         end do
         alpha=(rho+eps)/(rtvk+eps)

         do i=2,nxp1
         do j=2,nyp1
            sk(i,j)=res(i,j)-(alpha+eps)*(vk(i,j)+eps)
         end do
         end do

         do i=2,nxp1
         do j=2,nyp1
            qk(i,j)= sk(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,nyp1
            qk(i,j)=( sk(i,j)-aw(i,j)*qk(i-1,j) 
     &                       -as(i,j)*qk(i,j-1)  )*d(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,nyp1
            qk(i,j)=(qk(i,j)+eps)/(d(i,j)+eps)
         end do
         end do

         do i=nxp1,2,-1
         do j=nyp1,2,-1
            qk(i,j)=( qk(i,j)-ae(i,j)*qk(i+1,j)
     &                       -an(i,j)*qk(i,j+1)  )*d(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,nyp1
            rk(i,j)=( ae(i,j)*qk(i+1,j)+aw(i,j)*qk(i-1,j)
     &               +an(i,j)*qk(i,j+1)+as(i,j)*qk(i,j-1)
     &               +ap(i,j)*qk(i,j) )          
         end do
         end do

         rksk=0.0d0
         rkrk=0.0d0
         do i=2,nxp1
         do j=2,nyp1
            rksk=rksk+(rk(i,j)+eps)*(sk(i,j)+eps)
            rkrk=rkrk+(rk(i,j)+eps)*(rk(i,j)+eps)
         end do
         end do
         omega=(rksk+eps)/(rkrk+eps)

         rnorm=0.0d0
         psavg=0.0d0
         do i=2,nxp1
         do j=2,nyp1
            ps(i,j)=ps(i,j)+(alpha+eps)*(zk(i,j)+eps)
     &                     +(omega+eps)*(qk(i,j)+eps)
            psavg=psavg+ps(i,j)
            res(i,j)=(sk(i,j)+eps)-(omega+eps)*(rk(i,j)+eps)
            rnorm=dmax1(rnorm,dabs(res(i,j)))
         end do
         end do
         psavg=psavg/dfloat(nx*ny)

         if(rnorm.lt.resmax)then

            print *,'ps iter: ',l,rnorm

            ps(     0,0:nyp3)=ps(     2,0:nyp3)
            ps(     1,0:nyp3)=ps(     2,0:nyp3)
            ps(  nxp2,0:nyp3)=ps(  nxp1,0:nyp3)
            ps(  nxp3,0:nyp3)=ps(  nxp1,0:nyp3)
            ps(0:nxp3,     0)=ps(0:nxp3,     2)
            ps(0:nxp3,     1)=ps(0:nxp3,     2)
            ps(0:nxp3,  nyp2)=ps(0:nxp3,  nyp1)
            ps(0:nxp3,  nyp3)=ps(0:nxp3,  nyp1)
         return
         endif

      end do

      print *,'max ps iter: ',l,rnorm

      ps(     0,0:nyp3)=ps(     2,0:nyp3)
      ps(     1,0:nyp3)=ps(     2,0:nyp3)
      ps(  nxp2,0:nyp3)=ps(  nxp1,0:nyp3)
      ps(  nxp3,0:nyp3)=ps(  nxp1,0:nyp3)
      ps(0:nxp3,     0)=ps(0:nxp3,     2)
      ps(0:nxp3,     1)=ps(0:nxp3,     2)
      ps(0:nxp3,  nyp2)=ps(0:nxp3,  nyp1)
      ps(0:nxp3,  nyp3)=ps(0:nxp3,  nyp1)

      return
      end


c***********************************************************************     
      subroutine imp_u_bicgstab2D(resmax,maxit)
c***********************************************************************     
      include 'mac3d.inc'
      dimension res(0:maxnx,0:maxny),reso(0:maxnx,0:maxny)
      dimension pk(0:maxnx,0:maxny),rk(0:maxnx,0:maxny)
      dimension zk(0:maxnx,0:maxny),vk(0:maxnx,0:maxny)
      dimension sk(0:maxnx,0:maxny),qk(0:maxnx,0:maxny)
      dimension ae(0:maxnx,0:maxny),aw(0:maxnx,0:maxny)
      dimension as(0:maxnx,0:maxny),an(0:maxnx,0:maxny)
      dimension ap(0:maxnx,0:maxny),d(0:maxnx,0:maxny)
      dimension xk(0:maxnx,0:maxny)

      eps=1.0d-30

      res=0.0d0
      reso=0.0d0
      pk=0.0d0
      rk=0.0d0
      zk=0.0d0
      vk=0.0d0
      sk=0.0d0
      qk=0.0d0
      xk=0.0d0
      ae=0.0d0
      aw=0.0d0
      as=0.0d0
      an=0.0d0
      ap=0.0d0
      d=0.0d0
      utt=0.0d0
     
      do j=2,nyp1
      do i=2,nx
         ae(i,j)=-4.0d0*dt*hxi*hxi*vs(i+1,j)
     &                                  /(r(i+1,j)+r(i,j))*xs(i+1)/x(i)
         aw(i,j)=-4.0d0*dt*hxi*hxi*vs(i,j)/(r(i+1,j)+r(i,j))*xs(i)/x(i)
         an(i,j)=-0.5d0*dt*hyi*hyi*(vs(i,j)+vs(i+1,j)
     &               +vs(i,j+1)+vs(i+1,j+1))/(r(i+1,j)+r(i,j))
         as(i,j)=-0.5d0*dt*hyi*hyi*(vs(i,j-1)+vs(i+1,j-1)+
     &                    vs(i,j)+vs(i+1,j))/(r(i+1,j)+r(i,j))

         ap(i,j)=1.0d0-ae(i,j)-aw(i,j)-an(i,j)-as(i,j)

         res(i,j)=ut(i,j)-( ae(i,j)*utt(i+1,j)+aw(i,j)*utt(i-1,j)
     &                     +an(i,j)*utt(i,j+1)+as(i,j)*utt(i,j-1)
     &                     +ap(i,j)*utt(i,j) )          
         reso(i,j)=res(i,j)
      end do
      end do

      do j=2,nyp1
      do i=2,nx
         d(i,j)=1.0d0/(ap(i,j)
     &              -aw(i,j)*d(i-1,j)*(ae(i-1,j)+0.995d0*an(i-1,j)) 
     &              -as(i,j)*d(i,j-1)*(an(i,j-1)+0.995d0*ae(i,j-1)))
      end do 
      end do

      do l=1,maxit

         rho_o=rho
         rho=0.0d0
         do j=2,nyp1
         do i=2,nx
            rho=rho+res(i,j)*reso(i,j)
         end do
         end do

         if(l.eq.1)then
            do j=2,nyp1
            do i=2,nx
                pk(i,j)=res(i,j)
            end do
            end do
         else
            beta=(rho/rho_o)*(alpha/omega)
            do j=2,nyp1
            do i=2,nx
                pk(i,j)=res(i,j)+beta*(pk(i,j)-omega*vk(i,j))
            end do
            end do
         endif

         do j=2,nyp1
         do i=2,nx
            zk(i,j)=( pk(i,j)-aw(i,j)*zk(i-1,j) 
     &                       -as(i,j)*zk(i,j-1)  )*d(i,j)
         end do
         end do

         do j=2,nyp1
         do i=2,nx
            zk(i,j)=zk(i,j)/(d(i,j)+eps)
         end do
         end do

         do j=nyp1,2,-1
         do i=nx,2,-1
            zk(i,j)=( zk(i,j)-ae(i,j)*zk(i+1,j)
     &                       -an(i,j)*zk(i,j+1)  )*d(i,j)
         end do
         end do

         do j=2,nyp1
         do i=2,nx
            vk(i,j)=( ae(i,j)*zk(i+1,j)+aw(i,j)*zk(i-1,j)
     &               +an(i,j)*zk(i,j+1)+as(i,j)*zk(i,j-1)
     &               +ap(i,j)*zk(i,j) )          
         end do
         end do
 
         rtvk=zero
         do j=2,nyp1
         do i=2,nx
            rtvk=rtvk+reso(i,j)*vk(i,j)
         end do
         end do
         alpha=rho/(rtvk+eps)

         do j=2,nyp1
         do i=2,nx
            sk(i,j)=res(i,j)-alpha*vk(i,j)
         end do
         end do

         do j=2,nyp1
         do i=2,nx
            qk(i,j)= sk(i,j)
         end do
         end do

         do j=2,nyp1
         do i=2,nx
            qk(i,j)=( sk(i,j)-aw(i,j)*qk(i-1,j) 
     &                       -as(i,j)*qk(i,j-1)  )*d(i,j)
         end do
         end do

         do j=2,nyp1
         do i=2,nx
            qk(i,j)=qk(i,j)/(d(i,j)+eps)
         end do
         end do

         do j=nyp1,2,-1
         do i=nx,2,-1
            qk(i,j)=( qk(i,j)-ae(i,j)*qk(i+1,j)
     &                       -an(i,j)*qk(i,j+1)  )*d(i,j)
         end do
         end do

         do j=2,nyp1
         do i=2,nx
            rk(i,j)=( ae(i,j)*qk(i+1,j)+aw(i,j)*qk(i-1,j)
     &               +an(i,j)*qk(i,j+1)+as(i,j)*qk(i,j-1)
     &               +ap(i,j)*qk(i,j) )          
         end do
         end do

         rksk=0.0d0
         rkrk=0.0d0
         do j=2,nyp1
         do i=2,nx
            rksk=rksk+rk(i,j)*sk(i,j)
            rkrk=rkrk+rk(i,j)*rk(i,j)
         end do
         end do
         omega=rksk/(rkrk+eps)

         rnorm=0.0d0
         do j=2,nyp1
         do i=2,nx
            utt(i,j)=utt(i,j)+alpha*zk(i,j)+omega*qk(i,j)
            res(i,j)=sk(i,j)-omega*rk(i,j)
            rnorm=dmax1(rnorm,dabs(res(i,j)))
         end do
         end do

         if(rnorm.lt.resmax)then

            print *,'u iter: ',l,rnorm

            utt(     1,0:nyp3)= zero
            utt(  nxp1,0:nyp3)= zero
            utt(0:nxp3,     1)= zero
            utt(0:nxp3,  nyp2)= zero
         return
         endif

      end do

      return
      end
      

c***********************************************************************     
      subroutine imp_v_bicgstab2D(resmax,maxit)
c***********************************************************************     
      include 'mac3d.inc'
      dimension res(0:maxnx,0:maxny),reso(0:maxnx,0:maxny)
      dimension pk(0:maxnx,0:maxny),rk(0:maxnx,0:maxny)
      dimension zk(0:maxnx,0:maxny),vk(0:maxnx,0:maxny)
      dimension sk(0:maxnx,0:maxny),qk(0:maxnx,0:maxny)
      dimension ae(0:maxnx,0:maxny),aw(0:maxnx,0:maxny)
      dimension as(0:maxnx,0:maxny),an(0:maxnx,0:maxny)
      dimension ap(0:maxnx,0:maxny),d(0:maxnx,0:maxny)
      dimension xk(0:maxnx,0:maxny)

      eps=1.0d-30

      res=0.0d0
      reso=0.0d0
      pk=0.0d0
      rk=0.0d0
      zk=0.0d0
      vk=0.0d0
      sk=0.0d0
      qk=0.0d0
      xk=0.0d0
      ae=0.0d0
      aw=0.0d0
      as=0.0d0
      an=0.0d0
      ap=0.0d0
      d=0.0d0
      vtt=0.0d0

      do i=2,nxp1
      do j=2,ny
         if(i.eq.2)then
            aw(i,j)=0.0d0
         else
            aw(i,j)=-0.5d0*dt*hxi*hxi*(vs(i-1,j)+vs(i,j)
     &            +vs(i-1,j+1)+vs(i,j+1))/(r(i,j+1)+r(i,j))*x(i-1)/xs(i)
         endif
         ae(i,j)=-0.5d0*dt*hxi*hxi*(vs(i,j)+vs(i+1,j)
     &              +vs(i,j+1)+vs(i+1,j+1))/(r(i,j+1)+r(i,j))*x(i)/xs(i)
         an(i,j)=-4.0d0*dt*hyi*hyi*vs(i,j+1)/(r(i,j+1)+r(i,j))
         as(i,j)=-4.0d0*dt*hyi*hyi*vs(i,j)/(r(i,j+1)+r(i,j))
         ap(i,j)=1.0d0-ae(i,j)-aw(i,j)-an(i,j)-as(i,j)

         res(i,j)=vt(i,j)-( ae(i,j)*vtt(i+1,j)+aw(i,j)*vtt(i-1,j)
     &                     +an(i,j)*vtt(i,j+1)+as(i,j)*vtt(i,j-1)
     &                     +ap(i,j)*vtt(i,j) )          
         reso(i,j)=res(i,j)
      end do
      end do

      do i=2,nxp1
      do j=2,ny
         d(i,j)=1.0d0/(ap(i,j)
     &              -aw(i,j)*d(i-1,j)*(ae(i-1,j)+0.995d0*an(i-1,j)) 
     &              -as(i,j)*d(i,j-1)*(an(i,j-1)+0.995d0*ae(i,j-1)))
      end do 
      end do

      do l=1,maxit

         rho_o=rho
         rho=0.0d0
         do i=2,nxp1
         do j=2,ny
            rho=rho+res(i,j)*reso(i,j)
         end do
         end do

         if(l.eq.1)then
            do i=2,nxp1
            do j=2,ny
                pk(i,j)=res(i,j)
            end do
            end do
         else
            beta=(rho/rho_o)*(alpha/omega)
            do i=2,nxp1
            do j=2,ny
                pk(i,j)=res(i,j)+beta*(pk(i,j)-omega*vk(i,j))
            end do
            end do
         endif

         do i=2,nxp1
         do j=2,ny
            zk(i,j)=( pk(i,j)-aw(i,j)*zk(i-1,j) 
     &                       -as(i,j)*zk(i,j-1)  )*d(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,ny
            zk(i,j)=zk(i,j)/(d(i,j)+eps)
         end do
         end do

         do i=nxp1,2,-1
         do j=ny,2,-1
            zk(i,j)=( zk(i,j)-ae(i,j)*zk(i+1,j)
     &                       -an(i,j)*zk(i,j+1)  )*d(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,ny
            vk(i,j)=( ae(i,j)*zk(i+1,j)+aw(i,j)*zk(i-1,j)
     &               +an(i,j)*zk(i,j+1)+as(i,j)*zk(i,j-1)
     &               +ap(i,j)*zk(i,j) )          
         end do
         end do
 
         rtvk=zero
         do i=2,nxp1
         do j=2,ny
            rtvk=rtvk+reso(i,j)*vk(i,j)
         end do
         end do
         alpha=rho/(rtvk+eps)

         do i=2,nxp1
         do j=2,ny
            sk(i,j)=res(i,j)-alpha*vk(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,ny
            qk(i,j)= sk(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,ny
            qk(i,j)=( sk(i,j)-aw(i,j)*qk(i-1,j) 
     &                       -as(i,j)*qk(i,j-1)  )*d(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,ny
            qk(i,j)=qk(i,j)/(d(i,j)+eps)
         end do
         end do

         do i=nxp1,2,-1
         do j=ny,2,-1
            qk(i,j)=( qk(i,j)-ae(i,j)*qk(i+1,j)
     &                       -an(i,j)*qk(i,j+1)  )*d(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,ny
            rk(i,j)=( ae(i,j)*qk(i+1,j)+aw(i,j)*qk(i-1,j)
     &               +an(i,j)*qk(i,j+1)+as(i,j)*qk(i,j-1)
     &               +ap(i,j)*qk(i,j) )          
         end do
         end do

         rksk=0.0d0
         rkrk=0.0d0
         do i=2,nxp1
         do j=2,ny
            rksk=rksk+rk(i,j)*sk(i,j)
            rkrk=rkrk+rk(i,j)*rk(i,j)
         end do
         end do
         omega=rksk/(rkrk+eps)

         rnorm=0.0d0
         do i=2,nxp1
         do j=2,ny
            vtt(i,j)=vtt(i,j)+alpha*zk(i,j)+omega*qk(i,j)
            res(i,j)=sk(i,j)-omega*rk(i,j)
            rnorm=dmax1(rnorm,dabs(res(i,j)))
         end do
         end do

         if(rnorm.lt.resmax)then

            print *,'v iter: ',l,rnorm

            vtt(     1,0:nyp3)= vtt(     2,0:nyp3)
            vtt(  nxp2,0:nyp3)= zero
            vtt(0:nxp3,     1)= zero
            vtt(0:nxp3,  nyp1)= zero
         return
         endif

      end do

      return
      end
           
      
c***********************************************************************     
      subroutine imp_C_bicgstab2D(resmax,maxit)
c***********************************************************************     
      include 'mac3d.inc'
      dimension res(0:maxnx,0:maxny),reso(0:maxnx,0:maxny)
      dimension pk(0:maxnx,0:maxny),rk(0:maxnx,0:maxny)
      dimension zk(0:maxnx,0:maxny),vk(0:maxnx,0:maxny)
      dimension sk(0:maxnx,0:maxny),qk(0:maxnx,0:maxny)
      dimension ae(0:maxnx,0:maxny),aw(0:maxnx,0:maxny)
      dimension as(0:maxnx,0:maxny),an(0:maxnx,0:maxny)
      dimension ap(0:maxnx,0:maxny),d(0:maxnx,0:maxny)
      dimension xk(0:maxnx,0:maxny)

      eps=1.0d-30

      res=0.0d0
      reso=0.0d0
      pk=0.0d0
      rk=0.0d0
      zk=0.0d0
      vk=0.0d0
      sk=0.0d0
      qk=0.0d0
      xk=0.0d0
      ae=0.0d0
      aw=0.0d0
      as=0.0d0
      an=0.0d0
      ap=0.0d0
      d=0.0d0
      Cpp=0.0d0

      do i=2,nxp1
      do j=2,nyp1
         if(i.eq.2)then
            aw(i,j)=0.0d0
         else
            aw(i,j)=-Dc2*x(i-1)/xs(i)*hxi*hxi*dt
         endif
         ae(i,j)=-Dc2*x(i)/xs(i)*hxi*hxi*dt
         an(i,j)=-Dc2*hyi*hyi*dt
         as(i,j)=-Dc2*hyi*hyi*dt
         ap(i,j)=1.0d0-ae(i,j)-aw(i,j)-an(i,j)-as(i,j)

         res(i,j)=Cp(i,j)-( ae(i,j)*Cpp(i+1,j)+aw(i,j)*Cpp(i-1,j)
     &                     +an(i,j)*Cpp(i,j+1)+as(i,j)*Cpp(i,j-1)
     &                     +ap(i,j)*Cpp(i,j) )          
         reso(i,j)=res(i,j)
      end do
      end do

      do i=2,nxp1
      do j=2,nyp1
         d(i,j)=1.0d0/(ap(i,j)
     &              -aw(i,j)*d(i-1,j)*(ae(i-1,j)+0.995d0*an(i-1,j)) 
     &              -as(i,j)*d(i,j-1)*(an(i,j-1)+0.995d0*ae(i,j-1)))
      end do 
      end do

      do l=1,maxit

         rho_o=rho
         rho=0.0d0
         do i=2,nxp1
         do j=2,nyp1
            rho=rho+res(i,j)*reso(i,j)
         end do
         end do

         if(l.eq.1)then
            do i=2,nxp1
            do j=2,nyp1
                pk(i,j)=res(i,j)
            end do
            end do
         else
            beta=(rho/rho_o)*(alpha/omega)
            do i=2,nxp1
            do j=2,nyp1
                pk(i,j)=res(i,j)+beta*(pk(i,j)-omega*vk(i,j))
            end do
            end do
         endif

         do i=2,nxp1
         do j=2,nyp1
            zk(i,j)=( pk(i,j)-aw(i,j)*zk(i-1,j) 
     &                       -as(i,j)*zk(i,j-1)  )*d(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,nyp1
            zk(i,j)=zk(i,j)/(d(i,j)+eps)
         end do
         end do

         do i=nxp1,2,-1
         do j=nyp1,2,-1
            zk(i,j)=( zk(i,j)-ae(i,j)*zk(i+1,j)
     &                       -an(i,j)*zk(i,j+1)  )*d(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,nyp1
            vk(i,j)=( ae(i,j)*zk(i+1,j)+aw(i,j)*zk(i-1,j)
     &               +an(i,j)*zk(i,j+1)+as(i,j)*zk(i,j-1)
     &               +ap(i,j)*zk(i,j) )          
         end do
         end do
 
         rtvk=zero
         do i=2,nxp1
         do j=2,nyp1
            rtvk=rtvk+reso(i,j)*vk(i,j)
         end do
         end do
         alpha=rho/(rtvk+eps)

         do i=2,nxp1
         do j=2,nyp1
            sk(i,j)=res(i,j)-alpha*vk(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,nyp1
            qk(i,j)= sk(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,nyp1
            qk(i,j)=( sk(i,j)-aw(i,j)*qk(i-1,j) 
     &                       -as(i,j)*qk(i,j-1)  )*d(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,nyp1
            qk(i,j)=qk(i,j)/(d(i,j)+eps)
         end do
         end do

         do i=nxp1,2,-1
         do j=nyp1,2,-1
            qk(i,j)=( qk(i,j)-ae(i,j)*qk(i+1,j)
     &                       -an(i,j)*qk(i,j+1)  )*d(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,nyp1
            rk(i,j)=( ae(i,j)*qk(i+1,j)+aw(i,j)*qk(i-1,j)
     &               +an(i,j)*qk(i,j+1)+as(i,j)*qk(i,j-1)
     &               +ap(i,j)*qk(i,j) )          
         end do
         end do

         rksk=0.0d0
         rkrk=0.0d0
         do i=2,nxp1
         do j=2,nyp1
            rksk=rksk+rk(i,j)*sk(i,j)
            rkrk=rkrk+rk(i,j)*rk(i,j)
         end do
         end do
         omega=rksk/(rkrk+eps)

         rnorm=0.0d0
         do i=2,nxp1
         do j=2,nyp1
            Cpp(i,j)=Cpp(i,j)+alpha*zk(i,j)+omega*qk(i,j)
            res(i,j)=sk(i,j)-omega*rk(i,j)
            rnorm=dmax1(rnorm,dabs(res(i,j)))
         end do
         end do

         if(rnorm.lt.resmax)then

!            print *,'v iter: ',l,rnorm

            Cpp(     1,0:nyp3)= Cpp(     2,0:nyp3)
            Cpp(  nxp2,0:nyp3)= zero
            Cpp(0:nxp3,     1)= zero
            Cpp(0:nxp3,  nyp2)= zero
	      return
         endif

      end do

      return
      end
                        

c***********************************************************************     
      subroutine imp_rM_bicgstab2D(resmax,maxit)
c***********************************************************************     
      include 'mac3d.inc'
      dimension res(0:maxnx,0:maxny),reso(0:maxnx,0:maxny)
      dimension pk(0:maxnx,0:maxny),rk(0:maxnx,0:maxny)
      dimension zk(0:maxnx,0:maxny),vk(0:maxnx,0:maxny)
      dimension sk(0:maxnx,0:maxny),qk(0:maxnx,0:maxny)
      dimension ae(0:maxnx,0:maxny),aw(0:maxnx,0:maxny)
      dimension as(0:maxnx,0:maxny),an(0:maxnx,0:maxny)
      dimension ap(0:maxnx,0:maxny),d(0:maxnx,0:maxny)
      dimension xk(0:maxnx,0:maxny)

      eps=1.0d-30

      res=0.0d0
      reso=0.0d0
      pk=0.0d0
      rk=0.0d0
      zk=0.0d0
      vk=0.0d0
      sk=0.0d0
      qk=0.0d0
      xk=0.0d0
      ae=0.0d0
      aw=0.0d0
      as=0.0d0
      an=0.0d0
      ap=0.0d0
      d=0.0d0
      rMpp=0.0d0

      do i=2,nxp1
      do j=2,nyp1
         if(i.eq.2)then
            aw(i,j)=0.0d0
         else
            aw(i,j)=-Dm2*x(i-1)/xs(i)*hxi*hxi*dt
         endif
         ae(i,j)=-Dm2*x(i)/xs(i)*hxi*hxi*dt
         an(i,j)=-Dm2*hyi*hyi*dt
         as(i,j)=-Dm2*hyi*hyi*dt
         ap(i,j)=1.0d0-ae(i,j)-aw(i,j)-an(i,j)-as(i,j)

         res(i,j)=rMp(i,j)-( ae(i,j)*rMpp(i+1,j)+aw(i,j)*rMpp(i-1,j)
     &                      +an(i,j)*rMpp(i,j+1)+as(i,j)*rMpp(i,j-1)
     &                      +ap(i,j)*rMpp(i,j) )          
         reso(i,j)=res(i,j)
      end do
      end do

      do i=2,nxp1
      do j=2,nyp1
         d(i,j)=1.0d0/(ap(i,j)
     &              -aw(i,j)*d(i-1,j)*(ae(i-1,j)+0.995d0*an(i-1,j)) 
     &              -as(i,j)*d(i,j-1)*(an(i,j-1)+0.995d0*ae(i,j-1)))
      end do 
      end do

      do l=1,maxit

         rho_o=rho
         rho=0.0d0
         do i=2,nxp1
         do j=2,nyp1
            rho=rho+res(i,j)*reso(i,j)
         end do
         end do

         if(l.eq.1)then
            do i=2,nxp1
            do j=2,nyp1
                pk(i,j)=res(i,j)
            end do
            end do
         else
            beta=(rho/rho_o)*(alpha/omega)
            do i=2,nxp1
            do j=2,nyp1
                pk(i,j)=res(i,j)+beta*(pk(i,j)-omega*vk(i,j))
            end do
            end do
         endif

         do i=2,nxp1
         do j=2,nyp1
            zk(i,j)=( pk(i,j)-aw(i,j)*zk(i-1,j) 
     &                       -as(i,j)*zk(i,j-1)  )*d(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,nyp1
            zk(i,j)=zk(i,j)/(d(i,j)+eps)
         end do
         end do

         do i=nxp1,2,-1
         do j=nyp1,2,-1
            zk(i,j)=( zk(i,j)-ae(i,j)*zk(i+1,j)
     &                       -an(i,j)*zk(i,j+1)  )*d(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,nyp1
            vk(i,j)=( ae(i,j)*zk(i+1,j)+aw(i,j)*zk(i-1,j)
     &               +an(i,j)*zk(i,j+1)+as(i,j)*zk(i,j-1)
     &               +ap(i,j)*zk(i,j) )          
         end do
         end do
 
         rtvk=zero
         do i=2,nxp1
         do j=2,nyp1
            rtvk=rtvk+reso(i,j)*vk(i,j)
         end do
         end do
         alpha=rho/(rtvk+eps)

         do i=2,nxp1
         do j=2,nyp1
            sk(i,j)=res(i,j)-alpha*vk(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,nyp1
            qk(i,j)= sk(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,nyp1
            qk(i,j)=( sk(i,j)-aw(i,j)*qk(i-1,j) 
     &                       -as(i,j)*qk(i,j-1)  )*d(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,nyp1
            qk(i,j)=qk(i,j)/(d(i,j)+eps)
         end do
         end do

         do i=nxp1,2,-1
         do j=nyp1,2,-1
            qk(i,j)=( qk(i,j)-ae(i,j)*qk(i+1,j)
     &                       -an(i,j)*qk(i,j+1)  )*d(i,j)
         end do
         end do

         do i=2,nxp1
         do j=2,nyp1
            rk(i,j)=( ae(i,j)*qk(i+1,j)+aw(i,j)*qk(i-1,j)
     &               +an(i,j)*qk(i,j+1)+as(i,j)*qk(i,j-1)
     &               +ap(i,j)*qk(i,j) )          
         end do
         end do

         rksk=0.0d0
         rkrk=0.0d0
         do i=2,nxp1
         do j=2,nyp1
            rksk=rksk+rk(i,j)*sk(i,j)
            rkrk=rkrk+rk(i,j)*rk(i,j)
         end do
         end do
         omega=rksk/(rkrk+eps)

         rnorm=0.0d0
         do i=2,nxp1
         do j=2,nyp1
            rMpp(i,j)=rMpp(i,j)+alpha*zk(i,j)+omega*qk(i,j)
            res(i,j)=sk(i,j)-omega*rk(i,j)
            rnorm=dmax1(rnorm,dabs(res(i,j)))
         end do
         end do

         if(rnorm.lt.resmax)then

!            print *,'v iter: ',l,rnorm

            rMpp(     1,0:nyp3)= rMpp(     2,0:nyp3)
            rMpp(  nxp2,0:nyp3)= zero
            rMpp(0:nxp3,     1)= zero
            rMpp(0:nxp3,  nyp2)= zero
	      return
         endif

      end do

      return
      end
                        