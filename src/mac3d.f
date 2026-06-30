c***********************************************************************     
c        LEVEL CONTOUR RECONSTRUCTION MEHOTD WITH CONTACT MODEL
c
c                             Written by seungwon shin
c***********************************************************************     
   

c***********************************************************************     
      program mac3d
c***********************************************************************     
      include 'mac3d.inc'

      call inread 

      call intro 

      if (irestart.eq.0) then
         call init
       else
         call restartin
      endif

      open(15,file='history.dat',status='unknown')
      call solution
      close(15)

      stop
      end program
        
      
c***********************************************************************     
      subroutine inread
c***********************************************************************     
      include 'mac3d.inc'

      open(10,file='mac3d.in',status='unknown')

      read(10,*) dt
      read(10,*) nx
      read(10,*) ny
      read(10,*) xl
      read(10,*) yl
      read(10,*) r1	
      read(10,*) r2	
      read(10,*) vs1	
      read(10,*) vs2	
      read(10,*) surfT	
      read(10,*) the_re
      read(10,*) the_ad
      read(10,*) slip
      read(10,*) isurfactant_on
      read(10,*) Dc1
      read(10,*) Dc2
      read(10,*) rk_adsorpt
      read(10,*) rk_desorpt
      read(10,*) Dgm
      read(10,*) Gamma_inf
      read(10,*) Gamma_ini
      read(10,*) C_ini
      read(10,*) beta_s
      read(10,*) surf_eps
      read(10,*) imicelle_on
      read(10,*) Dm1
      read(10,*) Dm2
      read(10,*) rk_m_form
      read(10,*) rk_m_break
      read(10,*) rk_m_num
      read(10,*) rM_ini
      read(10,*) isubstrate_on
      read(10,*) Ds
      read(10,*) Cs_inf
      read(10,*) Cs_ini
      read(10,*) Cs_adsorpt
      read(10,*) Cs_desorpt
      read(10,*) Gm_Cs_adsorpt
      read(10,*) Gm_Cs_desorpt
      read(10,*) surf_slip_corr
      read(10,*) gry
      read(10,*) V_impact
      read(10,*) Radius
      read(10,*) pres
      read(10,*) irecon
      read(10,*) nfout 
      read(10,*) restep
      read(10,*) irestart
	  
      close(10)

      return
      end
     

c***********************************************************************     
      subroutine intro
c***********************************************************************     
      include 'mac3d.inc'
      
      ! DP change:
      CFL = 0.05d0
      zero=0.0d0
      half=0.5d0
      one=1.0d0
      two=2.0d0
      three=3.0d0
      four=4.0d0
      pi=four*datan(one)
      pd2=pi/two

      hx=xl/dfloat(nx)
      hy=yl/dfloat(ny)
      hxy=dmin1(hx,hy)
      hxh=hx/two
      hyh=hy/two
      hxi=dfloat(nx)/xl
      hyi=dfloat(ny)/yl
      nxp1=nx+1
      nyp1=ny+1
      nxp2=nx+2
      nyp2=ny+2
      nxp3=nx+3
      nyp3=ny+3

      eps=1.0d-5
      dsmin=1.0d0*hxy
      nf=0
      nrecon=1
      nrestart=1
      time=zero
      the_re_tmp=pi/180.0d0*the_re
      the_ad_tmp=pi/180.0d0*the_ad
      the_re=pi-the_ad_tmp
      the_ad=pi-the_re_tmp

      do i=-npad,nxp3+npad
         x(i)=hx*dfloat(i-1)
         xs(i)=x(i)-hxh
      enddo
      do j=-npad,nyp3+npad
         y(j)=hy*dfloat(j-1)
         ys(j)=y(j)-hyh
      enddo

      ps=zero
      u=zero
      v=zero
      ut=zero
      vt=zero
      utt=zero
      vtt=zero
      
      return
      end


c***********************************************************************     
      subroutine init
c***********************************************************************     
      include 'mac3d.inc'

      ! ycenter=Radius+hy
      ycenter=-Radius*dcos(pi/3.0d0) ! 60 degree initial angle 
      do j=0,nyp3
      do i=0,nxp3
         xcur=xs(i) 
         ycur=ys(j)      
         phi(i,j)=Radius-dsqrt(xcur*xcur+(ycur-ycenter)*(ycur-ycenter))
      enddo
      enddo

      ptrec=1.0d0
      call reconstruct

      Vol_ex=zero
      do j=2,nyp1
      do i=2,nxp1
         Vol_ex=Vol_ex+Heavi(Phi(i,j)/hxy)*hx*hy*2.0d0*pi*xs(i)
      enddo
      enddo

      do j=1,nyp3                                      
      do i=1,nxp3
         Hv1=Heavi(Phi(i,j)/hxy)
         Hv2=Heavi(Phi(i,j+1)/hxy)
         v(i,j)=V_impact*(Hv1+Hv2)*0.5d0
      enddo                                           
      enddo                                           

      if(isurfactant_on.eq.1)then
         do j=0,nyp3
         do i=0,nxp3
            C(i,j)=C_ini*Heavi(Phi(i,j)/hxy)
         enddo
         enddo
         Gmm_Xf=Gamma_ini  
         Gmm_ftn=Gamma_ini 

         if(imicelle_on.eq.1)then
            do j=0,nyp3
            do i=0,nxp3
               rM(i,j)=rM_ini*Heavi(Phi(i,j)/hxy)
            enddo
            enddo
         endif
         
         if(isubstrate_on.eq.1) Cs=Cs_ini
      endif
            
      call filenumber
      call print1
      write(*,*) 'printout # ',nf
      nf=nf+1

      return
      end


c***********************************************************************     
      subroutine solution
c***********************************************************************     
      include 'mac3d.inc'
      dimension Fu(0:maxnx,0:maxny),Fv(0:maxnx,0:maxny)
      dimension Phi_old(0:maxnx,0:maxny)
	  
      dt_ini = dt ! DP Change
      done=.false.
      do while (.not.done)

         call seekptn

         call contact_line

         call Find_distance

         do j=0,nyp3
         do i=0,nxp3
            r(i,j)=r1+(r2-r1)*Heavi(Phi(i,j)/hxy)
            vs(i,j)=vs1+(vs2-vs1)*Heavi(Phi(i,j)/hxy)
         enddo
         enddo

         if(isurfactant_on.eq.1)then
            if(isubstrate_on.eq.1) call substrate
            call surfactant
         endif

         call FMsource(Fu,Fv)

         if(isurfactant_on.eq.1)then
            call Concentration
            if(imicelle_on.eq.1) call Micelle
         endif

         call projection(Fu,Fv)

         Vol_cur=zero
         C_total=0.0d0
         rM_total=0.0d0         
         do j=2,nyp1
         do i=2,nxp1
            Vol_cur=Vol_cur
     &                 +Heavi(phi(i,j)/hxy)*hx*hy*2.0d0*pi*xs(i)
            C_total=C_total
     &          +C(i,j)*Heavi(phi(i,j)/hxy)*hx*hy*2.0d0*pi*xs(i)
            rM_total=rM_total
     &         +rM(i,j)*Heavi(phi(i,j)/hxy)*hx*hy*2.0d0*pi*xs(i)
         enddo
         enddo
         
         Cs_total=0.0d0         
         do i=2,nxp1
            Cs_total=Cs_total
     &         +Cs(i)*Heavi(phi(i,2)/hxy)*hx*2.0d0*pi*xs(i)
         enddo

         Gmass=0.0d0
         do m=1,nelemento
            ptxx=0.5d0*(ptxo(m,1)+ptxo(m,2))
            areaT=dsqrt( (ptxo(m,1)-ptxo(m,2))**two
     &                  +(ptyo(m,1)-ptyo(m,2))**two)
            Gmass=Gmass+Gmm_Xf(m)*areaT*2.0d0*pi*ptxx
         enddo

         if(time.ge.(pres-1.0d-10)*dfloat(nf)) then
            call filenumber
            call print1
            write(*,*) 'printout # ',nf
            if (nf.ge.nfout) then
               write(*,*) 'final output-time:',time
               done=.true.
            endif
            nf=nf+1
         endif    

         if(mod(int(time/dt),irecon).eq.0)then
            print *,'reconstructing surface...'
            Phi_old=Phi
            call Get_Phi
            call recon_fit
            nrecon=nrecon+1
            Phi=Phi_old

            nrecon=nrecon+1
         endif     

         ! DP Change!
         if(CFL.gt.zero.and.V_max.gt.zero)then
            !print *, V_max
            !pause
            dt_opt=CFL*hx/V_max
            ! dt_opt, dt_ini, CFL to include
            dt=min(dt_ini,dt_opt)
         endif       
         
         time=time+dt
         write(*,900)'t :',time,'#:',nelemento,ptx_max,pty_max,Gmass
         if(mod(int(time/dt),10).eq.0)then
            write(15,901)time,Vol_cur,ptx_max,pty_max,
     &                                Gmass,C_total,rM_total,Cs_total
         endif

         if(time.ge.(restep-1.0d-10)*dfloat(nrestart)) then
            write(6,*) 'writing restart file : ',nrestart
            call restartout
            nrestart=nrestart+1
         endif              

      end do
     
900   format(1x,a7,2x,f20.10,2x,a7,2x,i7,5(1x,e20.10))
901   format(20(e20.10,2x))

      return
      end   


c***********************************************************************     
      subroutine seekptn
c***********************************************************************     
      include 'mac3d.inc'

      uf=u
      vf=v

      if(isurfactant_on.eq.1)then
         Gmm_bt=0.0d0 
         do m=1,nelemento
	      x1=ptxo(m,1)
            y1=ptyo(m,1)
            x2=ptxo(m,2)
            y2=ptyo(m,2)
            if(iindx_pt(m,1).eq.1)then
               isf=floor((x1+hxh)*hxi)+1
               jsf=floor((y1+hyh)*hyi)+1
               pp=(x1-xs(isf))*hxi
               qq=(y1-ys(jsf))*hyi
               Gmm_bt=(1.0d0-pp)*(1.0d0-qq)*Gmm_ftn(isf  ,jsf  )
     &                      +pp *(1.0d0-qq)*Gmm_ftn(isf+1,jsf  )
     &               +(1.0d0-pp)*       qq *Gmm_ftn(isf  ,jsf+1)
     &                      +pp *       qq *Gmm_ftn(isf+1,jsf+1)
            endif
            if(iindx_pt(m,2).eq.1)then
               isf=floor((x2+hxh)*hxi)+1
               jsf=floor((y2+hyh)*hyi)+1
               pp=(x2-xs(isf))*hxi
               qq=(y2-ys(jsf))*hyi
               Gmm_bt=(1.0d0-pp)*(1.0d0-qq)*Gmm_ftn(isf  ,jsf  )
     &                      +pp *(1.0d0-qq)*Gmm_ftn(isf+1,jsf  )
     &               +(1.0d0-pp)*       qq *Gmm_ftn(isf  ,jsf+1)
     &                      +pp *       qq *Gmm_ftn(isf+1,jsf+1)
            endif
         enddo
         surf_coeff=dmax1(surf_eps,
     &       1.0d0+beta_s*log(1.0d0-dmin1(Gamma_inf,Gmm_bt)/Gamma_inf))
      else
         surf_coeff=1.0d0
      endif
      slip_factor=slip*dexp(-surf_slip_corr*(surf_coeff-1.0d0))
      
      uf(0:nxp3,1)=(4.0d0*slip_factor-1.0d0)*uf(0:nxp3,2)
      uf(0:nxp3,0)=uf(0:nxp3,1)

      ptx_max=0.0d0      
      ptx_min=xl      
      pty_max=0.0d0      
      pty_min=yl      

      Gmm_Ar=1.0d0
      do m=1,nelemento
         x1=ptxo(m,1)
         y1=ptyo(m,1)
         x2=ptxo(m,2)
         y2=ptyo(m,2)

         rnxx=y2-y1
         rnyy=-(x2-x1)
         rsum=dsqrt(rnxx**two+rnyy**two)
         rnxx=-rnxx/rsum
         rnyy=-rnyy/rsum

         ptxxp=0.5d0*(ptxo(m,1)+ptxo(m,2))
         areaTp=ptxxp*dsqrt( (ptxo(m,1)-ptxo(m,2))**two
     &                      +(ptyo(m,1)-ptyo(m,2))**two)
         do mm=1,2
            call seekvn(ptxo(m,mm),ptyo(m,mm),VVnx,VVny)
            if(iindx_pt(m,mm).eq.0)then
               rk1x=VVnx*dt
               rk1y=VVny*dt
               call seekvn(ptxo(m,mm)+rk1x,ptyo(m,mm)+rk1y,
     &                                                  VVnx,VVny)
               rk2x=VVnx*dt
               rk2y=VVny*dt
               ptxo(m,mm)=ptxo(m,mm)+(rk1x+rk2x)*0.5d0
               ptyo(m,mm)=ptyo(m,mm)+(rk1y+rk2y)*0.5d0

               ! DP Change
               V_test=sqrt(VVnx**two+VVny**two)
               if(V_test.gt.V_max)then
                  V_max=V_test
                  ! V_test to include
               endif

            else
               rnxs=0.0d0
               rnys=1.0d0
               ang_varc=acos(rnxx*rnxs+rnyy*rnys)
               if(ang_varc.gt.the_ad.or.ang_varc.lt.the_re)then
                  ptxo(m,mm)=ptxo(m,mm)+VVnx*dt
                  ptyo(m,mm)=ptyo(m,mm)+VVny*dt
               endif
            endif

            pty_max=dmax1(pty_max,ptyo(m,mm))
            pty_min=dmin1(pty_min,ptyo(m,mm))
            ptx_max=dmax1(ptx_max,ptxo(m,mm))
            ptx_min=dmin1(ptx_min,ptxo(m,mm))
		 enddo

         ptxxn=0.5d0*(ptxo(m,1)+ptxo(m,2))
         areaTn=ptxxn*dsqrt( (ptxo(m,1)-ptxo(m,2))**two
     &                      +(ptyo(m,1)-ptyo(m,2))**two)
         if ( areaTn.lt.1.0d-15 )then
            Gmm_Ar(m)=1.0d0
         else		  
            Gmm_Ar(m)=areaTp/areaTn
         end if		
      enddo

      ! DP Change
      do j=0,nyp3
      do i=0,nxp3
          V_test=sqrt(u(i,j)**two+v(i,j)**two)
          if(V_test.gt.V_max)then
             V_max=V_test
          endif
      enddo
      enddo

      return
      end


c***********************************************************************     
      subroutine contact_line
c***********************************************************************     
      include 'mac3d.inc'
	  
      nelement=nelemento
      ptx=ptxo
      pty=ptyo
	  
      Cs_con=0.0d0
      do m=1,nelemento
		 x1=ptxo(m,1)
         y1=ptyo(m,1)
         x2=ptxo(m,2)
         y2=ptyo(m,2)
		 
         areaT=dsqrt((x1-x2)**two+(y1-y2)**two)

         rnxx=y2-y1
         rnyy=-(x2-x1)
         rsum=dsqrt(rnxx**two+rnyy**two)
         rnxx=-rnxx/rsum
         rnyy=-rnyy/rsum
         
         if(iindx_pt(m,1).eq.1)then
            rnxs=0.0d0
            rnys=1.0d0
            rotc=dsign(1.0d0,rnxx*rnys-rnyy*rnxs)
            
            rntx=rnys*rotc
            rnty=-rnxs*rotc

            con_angc=acos(rnxx*rnxs+rnyy*rnys)
            if(con_angc.gt.the_ad)con_angc=the_ad
            if(con_angc.lt.the_re)con_angc=the_re

            rnbx=cos(con_angc*rotc)*rntx+sin(con_angc*rotc)*rnty
            rnby=-sin(con_angc*rotc)*rntx+cos(con_angc*rotc)*rnty
			
            nelement=nelement+1
            Gmm_Ar(nelement)=Gmm_Ar(m)
            Gmm_Xf(nelement)=Gmm_Xf(m)

            ptx(nelement,2)=x1
            pty(nelement,2)=y1
            iindx_pt(nelement,2)=1
          
            ptx(nelement,1)=x1+rnbx*hxy*1.0d0
            pty(nelement,1)=y1+rnby*hxy*1.0d0
            iindx_pt(nelement,1)=2

            isf=floor((x1+hxh)*hxi)+1
            jsf=floor((y1+hyh)*hyi)+1
            pp=(x1-xs(isf))*hxi
            qq=(y1-ys(jsf))*hyi
            Gm_contact=(1.0d0-pp)*(1.0d0-qq)*Gmm_ftn(isf  ,jsf  )
     &                       +pp *(1.0d0-qq)*Gmm_ftn(isf+1,jsf  )
     &                +(1.0d0-pp)*       qq *Gmm_ftn(isf  ,jsf+1)
     &                       +pp *       qq *Gmm_ftn(isf+1,jsf+1)
            Cs_contact=(1.0d0-pp)*Cs(isf)+pp*Cs(isf+1)

            Cs_src=Gm_Cs_adsorpt*(Gamma_inf-Gm_contact)
     &                               -Gm_Cs_desorpt*(Cs_inf-Cs_contact)
            do i=1,4
               iisf=isf-2+i
               xsf=(x1-xs(iisf))*hxi
               Cs_con(iisf)=Cs_con(iisf)+Cs_src*pes(xsf)*hxi
            enddo
            
         endif
 
         if(iindx_pt(m,2).eq.1)then
            rnxs=0.0d0
            rnys=1.0d0
            rotc=dsign(1.0d0,rnxx*rnys-rnyy*rnxs)
            
            rntx=rnys*rotc
            rnty=-rnxs*rotc

            con_angc=acos(rnxx*rnxs+rnyy*rnys)
            if(con_angc.gt.the_ad)con_angc=the_ad
            if(con_angc.lt.the_re)con_angc=the_re

            rnbx=cos(con_angc*rotc)*rntx+sin(con_angc*rotc)*rnty
            rnby=-sin(con_angc*rotc)*rntx+cos(con_angc*rotc)*rnty
			
            nelement=nelement+1
            Gmm_Ar(nelement)=Gmm_Ar(m)
            Gmm_Xf(nelement)=Gmm_Xf(m)

            ptx(nelement,1)=x2
            pty(nelement,1)=y2
            iindx_pt(nelement,1)=1
         
            ptx(nelement,2)=x2+rnbx*hxy*1.0d0
            pty(nelement,2)=y2+rnby*hxy*1.0d0
            iindx_pt(nelement,2)=2

            isf=floor((x2+hxh)*hxi)+1
            jsf=floor((y2+hyh)*hyi)+1
            pp=(x2-xs(isf))*hxi
            qq=(y2-ys(jsf))*hyi
            Gm_contact=(1.0d0-pp)*(1.0d0-qq)*Gmm_ftn(isf  ,jsf  )
     &                       +pp *(1.0d0-qq)*Gmm_ftn(isf+1,jsf  )
     &                +(1.0d0-pp)*       qq *Gmm_ftn(isf  ,jsf+1)
     &                       +pp *       qq *Gmm_ftn(isf+1,jsf+1)
            Cs_contact=(1.0d0-pp)*Cs(isf)+pp*Cs(isf+1)

            Cs_src=Gm_Cs_adsorpt*(Gamma_inf-Gm_contact)
     &                               -Gm_Cs_desorpt*(Cs_inf-Cs_contact)
            do i=1,4
               iisf=isf-2+i
               xsf=(x2-xs(iisf))*hxi
               Cs_con(iisf)=Cs_con(iisf)+Cs_src*pes(xsf)*hxi
            enddo
         endif
		 
	  enddo

      return
      end


c***********************************************************************     
      subroutine seekvn(ptxx_in,ptyy_in,uuf,vvf)
c***********************************************************************     
      include 'mac3d.inc'
 
      if(ptxx_in.lt.zero)then
         ptxx=-ptxx_in
      else
         ptxx=ptxx_in
      endif

      if(ptyy_in.lt.zero)then
         ptyy=0.0d0
      else
         ptyy=ptyy_in
      endif

      inf=floor(ptxx*hxi)+1
      jnf=floor(ptyy*hyi)+1
      isf=floor((ptxx+hxh)*hxi)+1
      jsf=floor((ptyy+hyh)*hyi)+1

      pp=(ptxx-x(inf))*hxi
      qq=(ptyy-ys(jsf))*hyi
      uuf=   (one-pp)*(one-qq)*uf(inf,jsf)
     &      +    pp*(one-qq)*uf(inf+1,jsf)
     &      +    (one-pp)*qq*uf(inf,jsf+1)
     &      +        pp*qq*uf(inf+1,jsf+1)
              
      pp=(ptxx-xs(isf))*hxi
      qq=(ptyy-y(jnf))*hyi
      vvf=   (one-pp)*(one-qq)*vf(isf,jnf)
     &      +    pp*(one-qq)*vf(isf+1,jnf)
     &      +    (one-pp)*qq*vf(isf,jnf+1)
     &      +        pp*qq*vf(isf+1,jnf+1)	  

      if(ptxx_in.lt.zero)then
         uuf=-uuf
      endif

      return
      end

  
c***********************************************************************     
      subroutine Concentration
c***********************************************************************     
      include 'mac3d.inc'

      do j=0,nyp3
      do i=0,nxp3
         if( phi(i,j).lt.0.0d0) then
            if( ptij(i,j).gt.0.5d0 )then 
               xbt=ptx_surf(i,j)
               ybt=pty_surf(i,j)
               rnxs=xbt-xs(i)
               rnys=ybt-ys(j)
               rsum=dsqrt(rnxs*rnxs+rnys*rnys)
               rnxs=rnxs/rsum
               rnys=rnys/rsum
               dis=dsqrt((xbt-xs(i))**2.0d0+(ybt-ys(j))**2.0d0)

               x_xf=xs(i)+(dis+1.5d0*hxy)*rnxs
               y_xf=ys(j)+(dis+1.5d0*hxy)*rnys

               isf=floor((x_xf+hxh)*hxi)+1
               jsf=floor((y_xf+hyh)*hyi)+1
               pp=(x_xf-xs(isf))*hxi
               qq=(y_xf-ys(jsf))*hyi
               C_xf=(1.0d0-pp)*(1.0d0-qq)*C(isf  ,jsf  )
     &                    +pp *(1.0d0-qq)*C(isf+1,jsf  )
     &             +(1.0d0-pp)*       qq *C(isf  ,jsf+1)
     &                    +pp *       qq *C(isf+1,jsf+1)

               S_xf=(1.0d0-pp)*(1.0d0-qq)*Gmm_src(isf  ,jsf  )
     &                    +pp *(1.0d0-qq)*Gmm_src(isf+1,jsf  )
     &             +(1.0d0-pp)*       qq *Gmm_src(isf  ,jsf+1)
     &                    +pp *       qq *Gmm_src(isf+1,jsf+1)

               C(i,j)=C_xf-S_xf/Dc2*(dis+1.5d0*hxy)
             else
               C(i,j)=0.0d0
             end if
         end if
      enddo	
      enddo	
      C(0,0:nyp3) = C(3,0:nyp3)
      C(1,0:nyp3) = C(2,0:nyp3)

      do j=2,nyp1
      do i=2,nxp1
         aR=u(i,j)
         if(aR.ne.0.0d0)then
            k1R=i+max(0,-int(dsign(1.0d0,aR)))
            aa=(C(k1R+1,j)-C(k1R,j))/2.0d0
            bb=(C(k1R,j)-C(k1R-1,j))/2.0d0
            if(dabs(aa).le.dabs(bb))then
               cc=aa
            else
               cc=bb
            endif
            CRx=C(k1R,j)+cc*dsign(1.0d0,aR)
         else
            CRx=(C(i,j)+C(i+1,j))*0.5d0
         endif
         aL=u(i-1,j)
         if(aL.ne.0.0d0)then
            k1L=i-1+max(0,-int(dsign(1.0d0,aL)))
            aa=(C(k1L+1,j)-C(k1L,j))/2.0d0
            bb=(C(k1L,j)-C(k1L-1,j))/2.0d0
            if(dabs(aa).le.dabs(bb))then
               cc=aa
            else
               cc=bb
            endif
            CLx=C(k1L,j)+cc*dsign(1.0d0,aL)
         else
            CLx=(C(i,j)+C(i-1,j))*0.5d0
         endif
         bR=v(i,j)
         if(bR.ne.0.0d0)then
            m1R=j+max(0,-int(dsign(1.0d0,bR)))
            aa=(C(i,m1R+1)-C(i,m1R))/2.0d0
            bb=(C(i,m1R)-C(i,m1R-1))/2.0d0
            if(dabs(aa).le.dabs(bb))then
               cc=aa
            else
               cc=bb
            endif
            CRy=C(i,m1R)+cc*dsign(1.0d0,bR)
         else
            CRy=(C(i,j)+C(i,j+1))*0.5d0
         endif
         bL=v(i,j-1)
         if(bL.ne.0.0d0)then
            m1L=j-1+max(0,-int(dsign(1.0d0,bL)))
            aa=(C(i,m1L+1)-C(i,m1L))/2.0d0
            bb=(C(i,m1L)-C(i,m1L-1))/2.0d0
            if(dabs(aa).le.dabs(bb))then
               cc=aa
            else
               cc=bb
            endif
            CLy=C(i,m1L)+cc*dsign(1.0d0,bL)
         else
            CLy=(C(i,j)+C(i,j-1))*0.5d0
         endif

         C_conv= -(x(i)*aR*CRx-x(i-1)*aL*CLx)*hxi/xs(i)
     &           -(bR*CRy-bL*CLy)*hyi

         C_RHS=  Dc2* (   (x(i)*(C(i+1,j)-C(i,j))
     &                    -x(i-1)*(C(i,j)-C(i-1,j)) )*hxi*hxi/xs(i)
     &                    +( C(i,j+1)+C(i,j-1)-two*C(i,j) )*hyi*hyi  )

         if(imicelle_on.eq.1)then
            C_src=rk_m_form*(C(i,j))**rk_m_num-rk_m_break*rM(i,j)
         else
            C_src=0.0d0
         endif
         	 
         Cp(i,j)=dt*(C_conv+C_RHS-rk_m_num*C_src)
      enddo
      enddo

      call imp_C_bicgstab2D(1.0d-10,10000)

      C(2:nxp1,2:nyp1)=C(2:nxp1,2:nyp1)+Cpp(2:nxp1,2:nyp1)

      C(     0,0:nyp3) = C(     3,0:nyp3)
      C(     1,0:nyp3) = C(     2,0:nyp3)
      C(  nxp2,0:nyp3) = C(  nxp1,0:nyp3)
      C(  nxp3,0:nyp3) = C(  nxp1,0:nyp3)
      C(0:nxp3,     1) = C(0:nxp3,     2)
      C(0:nxp3,     0) = C(0:nxp3,     2)
      C(0:nxp3,  nyp2) = C(0:nxp3,  nyp1)
      C(0:nxp3,  nyp3) = C(0:nxp3,  nyp1)

      do j=0,nyp3
      do i=0,nxp3
         if( phi(i,j).lt.0.0d0) then
            C(i,j)=0.0d0
         end if
      enddo	
      enddo	
      
      return
      end


c***********************************************************************     
      subroutine Micelle
c***********************************************************************     
      include 'mac3d.inc'

      do j=0,nyp3
      do i=0,nxp3
         if( phi(i,j).lt.0.0d0) then
            if( ptij(i,j).gt.0.5d0 )then 
               xbt=ptx_surf(i,j)
               ybt=pty_surf(i,j)
               rnxs=xbt-xs(i)
               rnys=ybt-ys(j)
               rsum=dsqrt(rnxs*rnxs+rnys*rnys)
               rnxs=rnxs/rsum
               rnys=rnys/rsum
               dis=dsqrt((xbt-xs(i))**2.0d0+(ybt-ys(j))**2.0d0)

               x_xf=xs(i)+(dis+1.5d0*hxy)*rnxs
               y_xf=ys(j)+(dis+1.5d0*hxy)*rnys

               isf=floor((x_xf+hxh)*hxi)+1
               jsf=floor((y_xf+hyh)*hyi)+1
               pp=(x_xf-xs(isf))*hxi
               qq=(y_xf-ys(jsf))*hyi
               rM_xf=(1.0d0-pp)*(1.0d0-qq)*rM(isf  ,jsf  )
     &                     +pp *(1.0d0-qq)*rM(isf+1,jsf  )
     &              +(1.0d0-pp)*       qq *rM(isf  ,jsf+1)
     &                     +pp *       qq *rM(isf+1,jsf+1)

               rM(i,j)=rM_xf
             else
               rM(i,j)=0.0d0
             end if
         end if
      enddo	
      enddo	
      rM(0,0:nyp3) = rM(3,0:nyp3)
      rM(1,0:nyp3) = rM(2,0:nyp3)

      do j=2,nyp1
      do i=2,nxp1
         aR=u(i,j)
         if(aR.ne.0.0d0)then
            k1R=i+max(0,-int(dsign(1.0d0,aR)))
            aa=(rM(k1R+1,j)-rM(k1R,j))/2.0d0
            bb=(rM(k1R,j)-rM(k1R-1,j))/2.0d0
            if(dabs(aa).le.dabs(bb))then
               cc=aa
            else
               cc=bb
            endif
            rMRx=rM(k1R,j)+cc*dsign(1.0d0,aR)
         else
            rMRx=(rM(i,j)+rM(i+1,j))*0.5d0
         endif
         aL=u(i-1,j)
         if(aL.ne.0.0d0)then
            k1L=i-1+max(0,-int(dsign(1.0d0,aL)))
            aa=(rM(k1L+1,j)-rM(k1L,j))/2.0d0
            bb=(rM(k1L,j)-rM(k1L-1,j))/2.0d0
            if(dabs(aa).le.dabs(bb))then
               cc=aa
            else
               cc=bb
            endif
            rMLx=rM(k1L,j)+cc*dsign(1.0d0,aL)
         else
            rMLx=(rM(i,j)+rM(i-1,j))*0.5d0
         endif
         bR=v(i,j)
         if(bR.ne.0.0d0)then
            m1R=j+max(0,-int(dsign(1.0d0,bR)))
            aa=(rM(i,m1R+1)-rM(i,m1R))/2.0d0
            bb=(rM(i,m1R)-rM(i,m1R-1))/2.0d0
            if(dabs(aa).le.dabs(bb))then
               cc=aa
            else
               cc=bb
            endif
            rMRy=rM(i,m1R)+cc*dsign(1.0d0,bR)
         else
            rMRy=(rM(i,j)+rM(i,j+1))*0.5d0
         endif
         bL=v(i,j-1)
         if(bL.ne.0.0d0)then
            m1L=j-1+max(0,-int(dsign(1.0d0,bL)))
            aa=(rM(i,m1L+1)-rM(i,m1L))/2.0d0
            bb=(rM(i,m1L)-rM(i,m1L-1))/2.0d0
            if(dabs(aa).le.dabs(bb))then
               cc=aa
            else
               cc=bb
            endif
            rMLy=rM(i,m1L)+cc*dsign(1.0d0,bL)
         else
            rMLy=(rM(i,j)+rM(i,j-1))*0.5d0
         endif

         rM_conv= -(x(i)*aR*rMRx-x(i-1)*aL*rMLx)*hxi/xs(i)
     &           -(bR*rMRy-bL*rMLy)*hyi

         rM_RHS= Dm2* (   (x(i)*(rM(i+1,j)-rM(i,j))
     &                  -x(i-1)*(rM(i,j)-rM(i-1,j)) )*hxi*hxi/xs(i)
     &                +( rM(i,j+1)+rM(i,j-1)-two*rM(i,j) )*hyi*hyi )

         rM_src=rk_m_form*(C(i,j))**rk_m_num-rk_m_break*rM(i,j)
	 
         rMp(i,j)=dt*(rM_conv+rM_RHS+rM_src)
      enddo
      enddo

      call imp_rM_bicgstab2D(1.0d-10,10000)

      rM(2:nxp1,2:nyp1)=rM(2:nxp1,2:nyp1)+rMpp(2:nxp1,2:nyp1)

      rM(     0,0:nyp3) = rM(     3,0:nyp3)
      rM(     1,0:nyp3) = rM(     2,0:nyp3)
      rM(  nxp2,0:nyp3) = rM(  nxp1,0:nyp3)
      rM(  nxp3,0:nyp3) = rM(  nxp1,0:nyp3)
      rM(0:nxp3,     1) = rM(0:nxp3,     2)
      rM(0:nxp3,     0) = rM(0:nxp3,     2)
      rM(0:nxp3,  nyp2) = rM(0:nxp3,  nyp1)
      rM(0:nxp3,  nyp3) = rM(0:nxp3,  nyp1)

      do j=0,nyp3
      do i=0,nxp3
         if( phi(i,j).lt.0.0d0) then
            rM(i,j)=0.0d0
         end if
      enddo	
      enddo	
      
      return
      end


c***********************************************************************     
      subroutine Substrate
c***********************************************************************     
      include 'mac3d.inc'

      do i=2,nxp1
         CXf=C(i,2)
         Cs_src=Cs_adsorpt*CXf*(Cs_inf-Cs(i))-Cs_desorpt*Cs(i)
         if(phi(i,2).gt.0.0d0) then
            if(phi(i+1,2).lt.0.0d0)then
               dCs_R=0.0d0
            else
               dCs_R=(Cs(i+1)-Cs(i))*hxi
            endif
            dCs_L=(Cs(i)-Cs(i-1))*hxi
            Cs_RHS= Ds*( (x(i)*dCs_R-x(i-1)*dCs_L)*hxi/xs(i) )
            Csn(i)=Cs(i)+dt*(Cs_RHS+Cs_src-Cs_con(i))

            C(i,1)=C(i,2)-Cs_src/Dc2*hy
         else
            Csn(i)=0.0d0
         endif
      enddo

      Cs=Csn
      
      Cs(     0) = Cs(     3)
      Cs(     1) = Cs(     2)
      Cs(  nxp2) = Cs(  nxp1)
      Cs(  nxp3) = Cs(  nxp1)

      return
      end
      

c***********************************************************************     
      subroutine Surfactant
c***********************************************************************     
      include 'mac3d.inc'
      dimension aux1(0:maxnx,0:maxny),aux2(0:maxnx,0:maxny)
      dimension aux3(0:maxnx,0:maxny)

      do m=1,nelemento
         x1=ptxo(m,1)
         y1=ptyo(m,1)
         x2=ptxo(m,2)
         y2=ptyo(m,2)
         
         ptxx=(x1+x2)*half
         ptyy=(y1+y2)*half

         rnxx=y2-y1
         rnyy=-(x2-x1)
         rsum=dsqrt(rnxx**two+rnyy**two)
         rnxx=-rnxx/rsum
         rnyy=-rnyy/rsum

         areaT=dsqrt((x1-x2)**two+(y1-y2)**two)

         t12x=(x2-x1)/areaT
         t12y=(y2-y1)/areaT
         
         call Search_Xf(x2,y2,t12x,t12y,xeval,yeval)
         io=floor((xeval+hxh)*hxi)+1
         jo=floor((yeval+hyh)*hyi)+1
         pp=(xeval-xs(io))*hxi
         qq=(yeval-ys(jo))*hyi
         Gmm_out=  (one-pp)*(one-qq)*Gmm_ftn(io,jo)
     &            +    pp*(one-qq)*Gmm_ftn(io+1,jo)
     &            +    (one-pp)*qq*Gmm_ftn(io,jo+1)
     &            +        pp*qq*Gmm_ftn(io+1,jo+1)
         call Search_Xf(x2,y2,-t12x,-t12y,xeval,yeval)
         io=floor((xeval+hxh)*hxi)+1
         jo=floor((yeval+hyh)*hyi)+1
         pp=(xeval-xs(io))*hxi
         qq=(yeval-ys(jo))*hyi
         Gmm_in=  (one-pp)*(one-qq)*Gmm_ftn(io,jo)
     &           +    pp*(one-qq)*Gmm_ftn(io+1,jo)
     &           +    (one-pp)*qq*Gmm_ftn(io,jo+1)
     &           +        pp*qq*Gmm_ftn(io+1,jo+1)
         Dels_Gm_2=(Gmm_out-Gmm_in)/hxy

         call Search_Xf(x1,y1,t12x,t12y,xeval,yeval)
         io=floor((xeval+hxh)*hxi)+1
         jo=floor((yeval+hyh)*hyi)+1
         pp=(xeval-xs(io))*hxi
         qq=(yeval-ys(jo))*hyi
         Gmm_out=  (one-pp)*(one-qq)*Gmm_ftn(io,jo)
     &            +    pp*(one-qq)*Gmm_ftn(io+1,jo)
     &            +    (one-pp)*qq*Gmm_ftn(io,jo+1)
     &            +        pp*qq*Gmm_ftn(io+1,jo+1)
         call Search_Xf(x1,y1,-t12x,-t12y,xeval,yeval)
         io=floor((xeval+hxh)*hxi)+1
         jo=floor((yeval+hyh)*hyi)+1
         pp=(xeval-xs(io))*hxi
         qq=(yeval-ys(jo))*hyi
         Gmm_in=  (one-pp)*(one-qq)*Gmm_ftn(io,jo)
     &           +    pp*(one-qq)*Gmm_ftn(io+1,jo)
     &           +    (one-pp)*qq*Gmm_ftn(io,jo+1)
     &           +        pp*qq*Gmm_ftn(io+1,jo+1)
         Dels_Gm_1=(Gmm_out-Gmm_in)/hxy

         io=floor((ptxx+hxh)*hxi)+1
         jo=floor((ptyy+hyh)*hyi)+1
         pp=(ptxx-xs(io))*hxi
         qq=(ptyy-ys(jo))*hyi
         CXf=  (one-pp)*(one-qq)*C(io,jo)
     &        +    pp*(one-qq)*C(io+1,jo)
     &        +    (one-pp)*qq*C(io,jo+1)
     &        +        pp*qq*C(io+1,jo+1)

         if(iindx_pt(m,1).eq.1)then
            isf=floor((x1+hxh)*hxi)+1
            pp=(x1-xs(isf))*hxi
            Cs_Src=(1.0d0-pp)*Cs_con(isf)+pp*Cs_con(isf+1)
         else if(iindx_pt(m,2).eq.1)then
            isf=floor((x2+hxh)*hxi)+1
            pp=(x2-xs(isf))*hxi
            Cs_Src=(1.0d0-pp)*Cs_con(isf)+pp*Cs_con(isf+1)
         else
            Cs_Src=0.0d0
         endif
         if(isubstrate_on.eq.0)Cs_Src=0.0d0

         Diff_s=Dgm*(x2*Dels_Gm_2-x1*Dels_Gm_1)/areaT/ptxx
         Source_s=rk_adsorpt*CXf*(Gamma_inf-Gmm_Xf(m))
     &                                -rk_desorpt*Gmm_Xf(m)
         Gmm_Sc(m)=Source_s
         Gmm_Xf(m)=Gmm_Ar(m)*Gmm_Xf(m)+(Diff_s+Source_s+Cs_Src)*dt
         Gmm_Xf(m)=dmax1(0.0d0,Gmm_Xf(m))
         Gmm_Xf(m)=dmin1(Gamma_inf,Gmm_Xf(m))
      enddo

      aux1=0.0d0
      aux2=0.0d0
      aux3=0.0d0
      do m=1,nelemento
		 x1=ptxo(m,1)
		 y1=ptyo(m,1)
		 x2=ptxo(m,2)
		 y2=ptyo(m,2)
	
         areaT=dsqrt((x1-x2)**two+(y1-y2)**two)

         ptxx=(x1+x2)*half
         ptyy=(y1+y2)*half

         isf=floor((ptxx+hxh)*hxi)+1
         jsf=floor((ptyy+hyh)*hyi)+1
         do j=1,4
         do i=1,4
            iisf=isf-2+i
            jjsf=jsf-2+j
            xsf=(ptxx-xs(iisf))*hxi
            ysf=(ptyy-ys(jjsf))*hyi
            aux1(iisf,jjsf)=aux1(iisf,jjsf)+pes(xsf)*pes(ysf)*areaT
            aux2(iisf,jjsf)=aux2(iisf,jjsf)+Gmm_Xf(m)
     &                                     *pes(xsf)*pes(ysf)*areaT
            aux3(iisf,jjsf)=aux3(iisf,jjsf)+Gmm_Sc(m)
     &                                     *pes(xsf)*pes(ysf)*areaT
         enddo
         enddo
      enddo

      Gmm_ftn=0.0d0
      Gmm_Src=0.0d0
      do j=0,nyp3
      do i=0,nxp3
         if(aux1(i,j).ne.0.0d0)then 
            Gmm_ftn(i,j)=aux2(i,j)/aux1(i,j)
            Gmm_Src(i,j)=aux3(i,j)/aux1(i,j)
         endif	 
      enddo
      enddo
      Gmm_ftn(0,0:nyp3) = Gmm_ftn(3,0:nyp3)
      Gmm_ftn(1,0:nyp3) = Gmm_ftn(2,0:nyp3)
      Gmm_Src(0,0:nyp3) = Gmm_Src(3,0:nyp3)
      Gmm_Src(1,0:nyp3) = Gmm_Src(2,0:nyp3)

      aux1=0.0d0
      aux2=0.0d0
      aux3=0.0d0
      do m=1,nelemento
		 x1=ptxo(m,1)
		 y1=ptyo(m,1)
		 x2=ptxo(m,2)
		 y2=ptyo(m,2)
	
         areaT=dsqrt((x1-x2)**two+(y1-y2)**two)

         ptxx=(x1+x2)*half
         ptyy=(y1+y2)*half

         io=floor((ptxx+hxh)*hxi)+1
         jo=floor((ptyy+hyh)*hyi)+1
         pp=(ptxx-xs(io))*hxi
         qq=(ptyy-ys(jo))*hyi
         Gmm_c=  (one-pp)*(one-qq)*Gmm_ftn(io,jo)
     &          +    pp*(one-qq)*Gmm_ftn(io+1,jo)
     &          +    (one-pp)*qq*Gmm_ftn(io,jo+1)
     &          +        pp*qq*Gmm_ftn(io+1,jo+1)

         isf=floor((ptxx+hxh)*hxi)+1
         jsf=floor((ptyy+hyh)*hyi)+1
         do j=1,4
         do i=1,4
            iisf=isf-2+i
            jjsf=jsf-2+j
            xsf=(ptxx-xs(iisf))*hxi
            ysf=(ptyy-ys(jjsf))*hyi
            aux1(iisf,jjsf)=aux1(iisf,jjsf)+pes(xsf)*pes(ysf)*areaT
            aux2(iisf,jjsf)=aux2(iisf,jjsf)+(Gmm_Xf(m)-Gmm_c)
     &                                     *pes(xsf)*pes(ysf)*areaT
         enddo
         enddo
      enddo

      do j=0,nyp3
      do i=0,nxp3
         if(aux1(i,j).ne.0.0d0)then 
            aux3(i,j)=dmax1(Gmm_ftn(i,j)+aux2(i,j)/aux1(i,j),0.0d0)
            aux3(i,j)=dmin1(aux3(i,j),Gamma_inf)
         endif	 
      enddo
      enddo
      Gmm_ftn=aux3

      Gmm_ftn(0,0:nyp3) = Gmm_ftn(3,0:nyp3)
      Gmm_ftn(1,0:nyp3) = Gmm_ftn(2,0:nyp3)
      Gmm_ftn(0:nxp3,0) = Gmm_ftn(0:nxp3,2)
      Gmm_ftn(0:nxp3,1) = Gmm_ftn(0:nxp3,2)

      return
      end


c***********************************************************************     
      subroutine FMsource(Fu,Fv)
c***********************************************************************     
      include 'mac3d.inc'
      dimension Fu(0:maxnx,0:maxny),Fv(0:maxnx,0:maxny)
      dimension divfx(0:maxnx,0:maxny),divfy(0:maxnx,0:maxny)
      dimension Fuc(0:maxnx,0:maxny),Fvc(0:maxnx,0:maxny)
      dimension divfxc(0:maxnx,0:maxny),divfyc(0:maxnx,0:maxny)
      dimension cv1(0:maxnx,0:maxny),cv2(0:maxnx,0:maxny)
      dimension Fut(0:maxnx,0:maxny),Fvt(0:maxnx,0:maxny)
      dimension Fun(0:maxnx,0:maxny),Fvn(0:maxnx,0:maxny)

      Fu=zero
      Fv=zero
      Fuc=zero
      Fvc=zero
      divfx=zero
      divfy=zero
      divfxc=zero
      divfyc=zero
      Fut=zero
      Fvt=zero      
      do m=1,nelemento
         x1=ptxo(m,1)
         y1=ptyo(m,1)
         x2=ptxo(m,2)
         y2=ptyo(m,2)

         ptxx=(x1+x2)*half
         ptyy=(y1+y2)*half

         areaT=dsqrt((x1-x2)**two+(y1-y2)**two)

         rnxx=y2-y1
         rnyy=-(x2-x1)
         rsum=dsqrt(rnxx**two+rnyy**two)
         rnxx=-rnxx/rsum
         rnyy=-rnyy/rsum

         t12x=(x2-x1)/areaT
         t12y=(y2-y1)/areaT
         t21x=(x1-x2)/areaT
         t21y=(y1-y2)/areaT

         call Search_Xf(ptxx,ptyy,t12x,t12y,xeval,yeval)
         io=floor((xeval+hxh)*hxi)+1
         jo=floor((yeval+hyh)*hyi)+1
         pp=(xeval-xs(io))*hxi
         qq=(yeval-ys(jo))*hyi
         Gmm_out=  (one-pp)*(one-qq)*Gmm_ftn(io,jo)
     &            +    pp*(one-qq)*Gmm_ftn(io+1,jo)
     &            +    (one-pp)*qq*Gmm_ftn(io,jo+1)
     &            +        pp*qq*Gmm_ftn(io+1,jo+1)
         call Search_Xf(x2,y2,-t12x,-t12y,xeval,yeval)
         io=floor((xeval+hxh)*hxi)+1
         jo=floor((yeval+hyh)*hyi)+1
         pp=(xeval-xs(io))*hxi
         qq=(yeval-ys(jo))*hyi
         Gmm_in=  (one-pp)*(one-qq)*Gmm_ftn(io,jo)
     &           +    pp*(one-qq)*Gmm_ftn(io+1,jo)
     &           +    (one-pp)*qq*Gmm_ftn(io,jo+1)
     &           +        pp*qq*Gmm_ftn(io+1,jo+1)

         surf_out=dmax1(surf_eps,
     &       1.0d0+beta_s*log(1.0d0-dmin1(Gamma_inf,Gmm_out)/Gamma_inf))
         surf_in=dmax1(surf_eps,
     &        1.0d0+beta_s*log(1.0d0-dmin1(Gamma_inf,Gmm_in)/Gamma_inf))

         Dels_sig_12=(surf_out-surf_in)/hxy

         surf_coef=dmax1(surf_eps,
     &     1.0d0+beta_s*log(1.0d0-dmin1(Gamma_inf,Gmm_Xf(m))/Gamma_inf))
        
         if(iindx_pt(m,1).ne.2)then
         is1=floor((x1+hxh)*hxi)+1
         js1=floor((y1+hyh)*hyi)+1
         do 101 j=1,4
         do 101 i=1,4
            iis1=is1-2+i
            jjs1=js1-2+j
            xsf1=(x1-xs(iis1))*hxi
            ysf1=(y1-ys(jjs1))*hyi
            Fu(iis1,jjs1)=Fu(iis1,jjs1)
     &                       +surf_coef*t12x*pes(xsf1)*pes(ysf1)*hxi*hyi
            Fv(iis1,jjs1)=Fv(iis1,jjs1)
     &                       +surf_coef*t12y*pes(xsf1)*pes(ysf1)*hxi*hyi
            divfx(iis1,jjs1)=divfx(iis1,jjs1)
     &                      +half*rnxx*areaT*pes(xsf1)*pes(ysf1)*hxi*hyi
            divfy(iis1,jjs1)=divfy(iis1,jjs1)
     &                      +half*rnyy*areaT*pes(xsf1)*pes(ysf1)*hxi*hyi
101      continue
         endif

         if(iindx_pt(m,2).ne.2)then
         is2=floor((x2+hxh)*hxi)+1
         js2=floor((y2+hyh)*hyi)+1
         do 102 j=1,4
         do 102 i=1,4
            iis2=is2-2+i
            jjs2=js2-2+j
            xsf2=(x2-xs(iis2))*hxi
            ysf2=(y2-ys(jjs2))*hyi
            Fu(iis2,jjs2)=Fu(iis2,jjs2)
     &                       +surf_coef*t21x*pes(xsf2)*pes(ysf2)*hxi*hyi
            Fv(iis2,jjs2)=Fv(iis2,jjs2)
     &                       +surf_coef*t21y*pes(xsf2)*pes(ysf2)*hxi*hyi
            divfx(iis2,jjs2)=divfx(iis2,jjs2)
     &                      +half*rnxx*areaT*pes(xsf2)*pes(ysf2)*hxi*hyi
            divfy(iis2,jjs2)=divfy(iis2,jjs2)
     &                      +half*rnyy*areaT*pes(xsf2)*pes(ysf2)*hxi*hyi
102      continue
         endif

         in3=floor(ptxx*hxi)+1
         jn3=floor(ptyy*hyi)+1
         is3=floor((ptxx+hxh)*hxi)+1
         js3=floor((ptyy+hyh)*hyi)+1
         do 103 j=1,4
         do 103 i=1,4
            iin3=in3-2+i
            jjn3=jn3-2+j
            iis3=is3-2+i
            jjs3=js3-2+j
            xnf3=(ptxx-x(iin3))*hxi
            ynf3=(ptyy-y(jjn3))*hyi
            xsf3=(ptxx-xs(iis3))*hxi
            ysf3=(ptyy-ys(jjs3))*hyi
            Fuc(iis3,jjs3)=Fuc(iis3,jjs3)
     &            -surf_coef*rnxx*rnxx*pes(xsf3)*pes(ysf3)*hxi*hyi*areaT
            Fvc(iis3,jjs3)=Fvc(iis3,jjs3)
     &            -surf_coef*rnxx*rnyy*pes(xsf3)*pes(ysf3)*hxi*hyi*areaT
            divfxc(iis3,jjs3)=divfxc(iis3,jjs3)
     &                      +rnxx*ptxx*pes(xsf3)*pes(ysf3)*hxi*hyi*areaT
            divfyc(iis3,jjs3)=divfyc(iis3,jjs3)
     &                      +rnyy*ptxx*pes(xsf3)*pes(ysf3)*hxi*hyi*areaT
            Fut(iin3,jjs3)=Fut(iin3,jjs3)
     &               +Dels_sig_12*t12x*pes(xnf3)*pes(ysf3)*hxi*hyi*areaT
            Fvt(iis3,jjn3)=Fvt(iis3,jjn3)
     &               +Dels_sig_12*t12y*pes(xsf3)*pes(ynf3)*hxi*hyi*areaT
103      continue
      enddo

      if(isurfactant_on.eq.0)then
         Fut=0.0d0
         Fvt=0.0d0
      endif
      
      cv1=zero
      cv2=zero
      do 300 j=0,nyp3
      do 300 i=0,nxp3
         if(ptij(i,j).eq.1.0d0)then 
            xbt=ptx_surf(i,j)
            ybt=pty_surf(i,j)
            call Phi_evalc(xbt,ybt,Fus,Fu)
            call Phi_evalc(xbt,ybt,Fvs,Fv)
            call Phi_evalc(xbt,ybt,divsx,divfx)
            call Phi_evalc(xbt,ybt,divsy,divfy)
            Fuv=Fus*divsx+Fvs*divsy
            divs=divsx**two+divsy**two
            cv1(i,j)=Fuv/divs

            call Phi_evalc(xbt,ybt,Fusc,Fuc)
            call Phi_evalc(xbt,ybt,Fvsc,Fvc)
            call Phi_evalc(xbt,ybt,divsxc,divfxc)
            call Phi_evalc(xbt,ybt,divsyc,divfyc)
            Fuvc=Fusc*divsxc+Fvsc*divsyc
            divsc=divsxc**two+divsyc**two
            cv2(i,j)=Fuvc/divsc
         else
            cv1(i,j)=zero
            cv2(i,j)=zero
         endif
300   continue
      cv1(0,0:nyp3)=cv1(3,0:nyp3)
      cv1(1,0:nyp3)=cv1(2,0:nyp3)
      cv2(0,0:nyp3)=cv2(3,0:nyp3)
      cv2(1,0:nyp3)=cv2(2,0:nyp3)

      Fu=zero
      Fv=zero
      do 400 j=0,nyp3
      do 400 i=0,nxp3
          cvx_deno=(ptij(i,j)+ptij(i+1,j))
          cvy_deno=(ptij(i,j)+ptij(i,j+1))
          cvx1=(cv1(i,j)*ptij(i,j)+cv1(i+1,j)*ptij(i+1,j))/cvx_deno
          if(cvx_deno.eq.zero)cvx1=zero
          cvy1=(cv1(i,j)*ptij(i,j)+cv1(i,j+1)*ptij(i,j+1))/cvy_deno
          if(cvy_deno.eq.zero)cvy1=zero
          cvx2=(cv2(i,j)*ptij(i,j)+cv2(i+1,j)*ptij(i+1,j))/cvx_deno
          if(cvx_deno.eq.zero)cvx2=zero
          cvy2=(cv2(i,j)*ptij(i,j)+cv2(i,j+1)*ptij(i,j+1))/cvy_deno
          if(cvy_deno.eq.zero)cvy2=zero

          Fu(i,j)=surfT*((cvx1+cvx2)*(Heavi(Phi(i+1,j)/hxy)
     &                            -Heavi(Phi(i,j)/hxy))*hxi+Fut(i,j))
          Fv(i,j)=surfT*((cvy1+cvy2)*(Heavi(Phi(i,j+1)/hxy)
     &                            -Heavi(Phi(i,j)/hxy))*hyi+Fvt(i,j))
400   continue 

      return
      end


c***********************************************************************     
      subroutine projection(Fu,Fv)
c***********************************************************************     
      include 'mac3d.inc'
      dimension Fu(0:maxnx,0:maxny),Fv(0:maxnx,0:maxny)

      do j=2,nyp1
      do i=2,nx
         aR=(u(i+1,j)+u(i,j))/2.0d0
         if(aR.ne.zero)then
            k1R=i+max(0,-int(dsign(1.0d0,aR)))
            aa=(u(k1R+1,j)-u(k1R,j))/2.0d0
            bb=(u(k1R,j)-u(k1R-1,j))/2.0d0
            if(dabs(aa).le.dabs(bb))then
               cc=aa
            else
               cc=bb
            endif
            uRx=u(k1R,j)+cc*dsign(1.0d0,aR)
         else
            uRx=(u(i,j)+u(i+1,j))*half
         endif
         aL=(u(i,j)+u(i-1,j))/2.0d0
         if(aL.ne.zero)then
            k1L=i-1+max(0,-int(dsign(1.0d0,aL)))
            aa=(u(k1L+1,j)-u(k1L,j))/2.0d0
            bb=(u(k1L,j)-u(k1L-1,j))/2.0d0
            if(dabs(aa).le.dabs(bb))then
               cc=aa
            else
               cc=bb
            endif
            uLx=u(k1L,j)+cc*dsign(1.0d0,aL)
         else
            uLx=(u(i,j)+u(i-1,j))*half
         endif

         bR=(v(i,j)+v(i+1,j))/2.0d0
         if(bR.ne.zero)then
            m1R=j+max(0,-int(dsign(1.0d0,bR)))
            aa=(u(i,m1R+1)-u(i,m1R))/2.0d0
            bb=(u(i,m1R)-u(i,m1R-1))/2.0d0
            if(dabs(aa).le.dabs(bb))then
               cc=aa
            else
               cc=bb
            endif
            uRy=u(i,m1R)+cc*dsign(1.0d0,bR)
         else
            uRy=(u(i,j)+u(i,j+1))*half
         endif

         bL=(v(i,j-1)+v(i+1,j-1))/2.0d0
         if(bL.ne.zero)then
            m1L=j-1+max(0,-int(dsign(1.0d0,bL)))
            aa=(u(i,m1L+1)-u(i,m1L))/2.0d0
            bb=(u(i,m1L)-u(i,m1L-1))/2.0d0
            if(dabs(aa).le.dabs(bb))then
               cc=aa
            else
               cc=bb
            endif
            uLy=u(i,m1L)+cc*dsign(1.0d0,bL)
         else
            uLy=(u(i,j)+u(i,j-1))*half
         endif

         Au_adv= -(xs(i+1)*aR*uRx-xs(i)*aL*uLx)*hxi/x(i)
     &           -(bR*uRy-bL*uLy)*hyi

         Au_vis=
     &     +(                        vs(i+1,j)*xs(i+1)*(u(i+1,j)-u(i,j))
     &                                  -vs(i,j)*xs(i)*(u(i,j)-u(i-1,j))
     &                                                )*two*hxi*hxi/x(i)
     &    +( (vs(i,j)+vs(i+1,j)+vs(i,j+1)+vs(i+1,j+1))*(u(i,j+1)-u(i,j))
     &     -(vs(i,j)+vs(i+1,j)+vs(i,j-1)+vs(i+1,j-1))*(u(i,j)-u(i,j-1)))
     &                                                     *hyi*hyi/four
     &    +( (vs(i,j)+vs(i+1,j)+vs(i,j+1)+vs(i+1,j+1))*(v(i+1,j)-v(i,j))
     &     -(vs(i,j)+vs(i+1,j)+vs(i,j-1)+vs(i+1,j-1))
     &                              *(v(i+1,j-1)-v(i,j-1)))*hxi*hyi/four
     &    -(vs(i,j)+vs(i+1,j))*u(i,j)/x(i)/x(i)
     
		 ut(i,j)=( Au_adv+grx
     &                      +(Fu(i,j)+Au_vis)*two/(r(i,j)+r(i+1,j)) )*dt
      enddo
      enddo

      call imp_u_bicgstab2D(1.0d-10,1000)

      do j=2,ny
      do i=2,nxp1
         aR=(u(i,j)+u(i,j+1))/2.0d0
         if(aR.ne.zero)then
            k1R=i+max(0,-int(dsign(1.0d0,aR)))
            aa=(v(k1R+1,j)-v(k1R,j))/2.0d0
            bb=(v(k1R,j)-v(k1R-1,j))/2.0d0
            if(dabs(aa).le.dabs(bb))then
               cc=aa
            else
               cc=bb
            endif
            vRx=v(k1R,j)+cc*dsign(1.0d0,aR)
         else
            vRx=(v(i,j)+v(i+1,j))/2.0d0
         endif
         aL=(u(i-1,j)+u(i-1,j+1))/2.0d0
         if(aL.ne.zero)then
            k1L=i-1+max(0,-int(dsign(1.0d0,aL)))
            aa=(v(k1L+1,j)-v(k1L,j))/2.0d0
            bb=(v(k1L,j)-v(k1L-1,j))/2.0d0
            if(dabs(aa).le.dabs(bb))then
               cc=aa
            else
               cc=bb
            endif
            vLx=v(k1L,j)+cc*dsign(1.0d0,aL)
         else
            vLx=(v(i,j)+v(i-1,j))/2.0d0
         endif

         bR=(v(i,j)+v(i,j+1))/2.0d0
         if(bR.ne.zero)then
            m1R=j+max(0,-int(dsign(1.0d0,bR)))
            aa=(v(i,m1R+1)-v(i,m1R))/2.0d0
            bb=(v(i,m1R)-v(i,m1R-1))/2.0d0
            if(dabs(aa).le.dabs(bb))then
               cc=aa
            else
               cc=bb
            endif
            vRy=v(i,m1R)+cc*dsign(1.0d0,bR)
         else
            vRy=(v(i,j)+v(i,j+1))/2.0d0
         endif
         bL=(v(i,j-1)+v(i,j))/2.0d0
         if(bL.ne.zero)then
            m1L=j-1+max(0,-int(dsign(1.0d0,bL)))
            aa=(v(i,m1L+1)-v(i,m1L))/2.0d0
            bb=(v(i,m1L)-v(i,m1L-1))/2.0d0
            if(dabs(aa).le.dabs(bb))then
               cc=aa
            else
               cc=bb
            endif
            vLy=v(i,m1L)+cc*dsign(1.0d0,bL)
         else
            vLy=(v(i,j-1)+v(i,j))/2.0d0
         endif

         Av_adv=-(x(i)*aR*vRx-x(i-1)*aL*vLx)*hxi/xs(i)
     &           -(bR*vRy-bL*vLy)*hyi

         Av_vis=                                            two*hyi*hyi*
     &           (vs(i,j+1)*(v(i,j+1)-v(i,j))-vs(i,j)*(v(i,j)-v(i,j-1)))
     &     +(  (vs(i,j)+vs(i+1,j)+vs(i,j+1)+vs(i+1,j+1))
     &                          *(v(i+1,j)-v(i,j))*x(i)
     &        -(vs(i,j)+vs(i-1,j)+vs(i,j+1)+vs(i-1,j+1))
     &                    *(v(i,j)-v(i-1,j))*x(i-1) )*hxi*hxi/four/xs(i)
     &     +(  (vs(i,j)+vs(i+1,j)+vs(i,j+1)+vs(i+1,j+1))
     &           *(u(i,j+1)-u(i,j))*x(i)
     &        -(vs(i,j)+vs(i-1,j)+vs(i,j+1)+vs(i-1,j+1))
     &              *(u(i-1,j+1)-u(i-1,j))*x(i-1)   )*hxi*hyi/four/xs(i)
     
          vt(i,j)=( Av_adv+gry
     &                      +(Fv(i,j)+Av_vis)*two/(r(i,j)+r(i,j+1)) )*dt
      enddo
      enddo

      call imp_v_bicgstab2D(1.0d-10,1000)

      ut(2:nx,2:nyp1)=u(2:nx,2:nyp1)+utt(2:nx,2:nyp1)
      vt(2:nxp1,2:ny)=v(2:nxp1,2:ny)+vtt(2:nxp1,2:ny)

      ut(     0,0:nyp3)= zero
      ut(     1,0:nyp3)= zero
      ut(  nxp1,0:nyp3)= zero
      ut(  nxp2,0:nyp3)= zero
      ut(  nxp3,0:nyp3)= zero
      ut(0:nxp3,     0)= -ut(0:nxp3,     3) ! Dp Change 
      ut(0:nxp3,     1)= -ut(0:nxp3,     2) ! Dp change
      ut(0:nxp3,  nyp2)= zero
      ut(0:nxp3,  nyp3)= zero
      
      vt(     0,0:nyp3)= vt(     3,0:nyp3)
      vt(     1,0:nyp3)= vt(     2,0:nyp3)
      vt(  nxp2,0:nyp3)= zero !vt(  nxp1,0:nyp3)
      vt(  nxp3,0:nyp3)= zero !vt(  nxp1,0:nyp3)
      vt(0:nxp3,     0)= zero
      vt(0:nxp3,     1)= zero
      vt(0:nxp3,  nyp1)= vt(0:nxp3,    ny)   !zero
      vt(0:nxp3,  nyp2)= vt(0:nxp3,    ny)   !zero
      vt(0:nxp3,  nyp3)= vt(0:nxp3,    ny)   !zero

      call cgstab(1.0d-10,100000)
      
      u(2:nx,2:nyp1)=ut(2:nx,2:nyp1)
     &       -dt*hxi*(ps(3:nxp1,2:nyp1)-ps(2:nx,2:nyp1))
     &            *two/(r(3:nxp1,2:nyp1)+r(2:nx,2:nyp1))

      v(2:nxp1,2:ny)=vt(2:nxp1,2:ny)
     &       -dt*hyi*(ps(2:nxp1,3:nyp1)-ps(2:nxp1,2:ny))
     &            *two/(r(2:nxp1,3:nyp1)+r(2:nxp1,2:ny))

      u(     0,0:nyp3)= zero
      u(     1,0:nyp3)= zero
      u(  nxp1,0:nyp3)= zero
      u(  nxp2,0:nyp3)= zero
      u(  nxp3,0:nyp3)= zero
      u(0:nxp3,     0)= -u(0:nxp3,     3) !zero
      u(0:nxp3,     1)= -u(0:nxp3,     2) !zero
      u(0:nxp3,  nyp2)= zero
      u(0:nxp3,  nyp3)= zero
      
      v(     0,0:nyp3)= v(     3,0:nyp3)
      v(     1,0:nyp3)= v(     2,0:nyp3)
      v(  nxp2,0:nyp3)= zero !v(  nxp1,0:nyp3)
      v(  nxp3,0:nyp3)= zero !v(  nxp1,0:nyp3)
      v(0:nxp3,     0)= zero
      v(0:nxp3,     1)= zero
      v(0:nxp3,  nyp1)= v(0:nxp3,    ny)!zero
      v(0:nxp3,  nyp2)= v(0:nxp3,    ny)!zero
      v(0:nxp3,  nyp3)= v(0:nxp3,    ny)!zero

      return
      end


c***********************************************************************     
      subroutine Get_Phi
c***********************************************************************     
      include 'mac3d.inc'
      dimension delIp(maxpt,2),fIp(maxpt,2)
      dimension Dist(0:maxnx,0:maxny)
      dimension Distp(0:maxnx,0:maxny)
      dimension weight(0:maxnx,0:maxny),coeff(maxpt,2)

      Weight=zero
      do m=1,nelemento
         x1=ptxo(m,1)
         y1=ptyo(m,1)
         x2=ptxo(m,2)
         y2=ptyo(m,2)

         call Phi_evalc(x1,y1,Dist1,Phi)
         call Phi_evalc(x2,y2,Dist2,Phi)

         delIp(m,1)=Dist1
         delIp(m,2)=Dist2
         
         isf=floor((x1+hxh)*hxi)+1
         jsf=floor((y1+hyh)*hyi)+1
         do j=1,4
         do i=1,4
            iisf=isf-2+i
            jjsf=jsf-2+j
            xsf=(x1-xs(iisf))*hxi
            ysf=(y1-ys(jjsf))*hyi
            if(        (iisf.ge.0.and.iisf.le.nxp3)
     &            .and.(jjsf.ge.0.and.jjsf.le.nyp3) )then
               Weight(iisf,jjsf)=Weight(iisf,jjsf)+fM3(xsf)*fM3(ysf)
            endif
         enddo
         enddo

         isf=floor((x2+hxh)*hxi)+1
         jsf=floor((y2+hyh)*hyi)+1
         do j=1,4
         do i=1,4
            iisf=isf-2+i
            jjsf=jsf-2+j
            xsf=(x2-xs(iisf))*hxi
            ysf=(y2-ys(jjsf))*hyi
            if(        (iisf.ge.0.and.iisf.le.nxp3)
     &            .and.(jjsf.ge.0.and.jjsf.le.nyp3) )then
               Weight(iisf,jjsf)=Weight(iisf,jjsf)+fM3(xsf)*fM3(ysf)
            endif
         enddo
         enddo         
      enddo

      do m=1,nelemento
         x1=ptxo(m,1)
         y1=ptyo(m,1)
         x2=ptxo(m,2)
         y2=ptyo(m,2)

         call Phi_evalc(x1,y1,Wei1,Weight)
         call Phi_evalc(x2,y2,Wei2,Weight)

         coeff(m,1)=Wei1
         coeff(m,2)=Wei2
      enddo         

      fIp=zero
      Dist=zero
      rmax_Dist=1.0d9
      niter=1
      do while(niter.lt.25)

         Distp=zero
         do m=1,nelemento
            x1=ptxo(m,1)
            y1=ptyo(m,1)
            x2=ptxo(m,2)
            y2=ptyo(m,2)

            areaT=dsqrt((x1-x2)**two+(y1-y2)**two)

            ptxx=(x1+x2)*half
            ptyy=(y1+y2)*half

            Dist1=fIp(m,1)
            Dist2=fIp(m,2)
            Distc=(Dist1+Dist2)*half

            isf=floor((x1+hxh)*hxi)+1
            jsf=floor((y1+hyh)*hyi)+1
            do j=1,4
            do i=1,4
               iisf=isf-2+i
               jjsf=jsf-2+j
               xsf=(x1-xs(iisf))*hxi
               ysf=(y1-ys(jjsf))*hyi
               if(        (iisf.ge.0.and.iisf.le.nxp3)
     &               .and.(jjsf.ge.0.and.jjsf.le.nyp3) )then
                  Distp(iisf,jjsf)=Distp(iisf,jjsf)
     &                            +Dist1*fM3(xsf)*fM3(ysf)/coeff(m,1)
               endif
            enddo
            enddo

            isf=floor((x2+hxh)*hxi)+1
            jsf=floor((y2+hyh)*hyi)+1
            do j=1,4
            do i=1,4
               iisf=isf-2+i
               jjsf=jsf-2+j
               xsf=(x2-xs(iisf))*hxi
               ysf=(y2-ys(jjsf))*hyi
               if(        (iisf.ge.0.and.iisf.le.nxp3)
     &               .and.(jjsf.ge.0.and.jjsf.le.nyp3) )then
                  Distp(iisf,jjsf)=Distp(iisf,jjsf)
     &                            +Dist2*fM3(xsf)*fM3(ysf)/coeff(m,2)
               endif
            enddo
            enddo
         enddo

         Distp(0,0:nyp3)=Distp(3,0:nyp3)
         Distp(1,0:nyp3)=Distp(2,0:nyp3)
         Dist=Dist+Distp
 
         rmax_Dist=-1.0d9
         do m=1,nelemento
            x1=ptxo(m,1)
            y1=ptyo(m,1)
            x2=ptxo(m,2)
            y2=ptyo(m,2)
            call Phi_evalc(x1,y1,Dist1,Dist)
            call Phi_evalc(x2,y2,Dist2,Dist)
            fIp_old1=fIp(m,1)
            fIp_old2=fIp(m,2)
            fIp(m,1)=delIp(m,1)-Dist1
            fIp(m,2)=delIp(m,2)-Dist2
            rmax_Dist=dmax1(rmax_Dist,dabs(fIp(m,1)),dabs(fIp(m,2)))
         enddo
         niter=niter+1
      enddo

      do j=0,nyp3
      do i=0,nxp3
         Phi(i,j)=Phi(i,j)-Dist(i,j)
      enddo
      enddo

      return
      end


c***********************************************************************     
      subroutine recon_fit
c***********************************************************************     
      include 'mac3d.inc'
      dimension Phi_org(0:maxnx,0:maxny)
      dimension aux1(0:maxnx,0:maxny),aux2(0:maxnx,0:maxny)
      dimension aux3(0:maxnx,0:maxny)

      if(isurfactant_on.eq.1)then
         aux1=0.0d0
         aux2=0.0d0
         aux3=0.0d0
         do m=1,nelemento
            x1=ptxo(m,1)
            y1=ptyo(m,1)
            x2=ptxo(m,2)
            y2=ptyo(m,2)
	
            areaT=dsqrt((x1-x2)**two+(y1-y2)**two)

            ptxx=(x1+x2)*half
            ptyy=(y1+y2)*half

            io=floor((ptxx+hxh)*hxi)+1
            jo=floor((ptyy+hyh)*hyi)+1
            pp=(ptxx-xs(io))*hxi
            qq=(ptyy-ys(jo))*hyi
            Gmm_c=  (one-pp)*(one-qq)*Gmm_ftn(io,jo)
     &             +    pp*(one-qq)*Gmm_ftn(io+1,jo)
     &             +    (one-pp)*qq*Gmm_ftn(io,jo+1)
     &             +        pp*qq*Gmm_ftn(io+1,jo+1)

            isf=floor((ptxx+hxh)*hxi)+1
            jsf=floor((ptyy+hyh)*hyi)+1
            do j=1,4
            do i=1,4
               iisf=isf-2+i
               jjsf=jsf-2+j
               xsf=(ptxx-xs(iisf))*hxi
               ysf=(ptyy-ys(jjsf))*hyi
               aux1(iisf,jjsf)=aux1(iisf,jjsf)+pes(xsf)*pes(ysf)*areaT
               aux2(iisf,jjsf)=aux2(iisf,jjsf)+(Gmm_Xf(m)-Gmm_c)
     &                                        *pes(xsf)*pes(ysf)*areaT
            enddo
            enddo
         enddo

         do j=0,nyp3
         do i=0,nxp3
            if(aux1(i,j).ne.0.0d0)then 
               aux3(i,j)=dmax1(Gmm_ftn(i,j)+aux2(i,j)/aux1(i,j),0.0d0)
               aux3(i,j)=dmin1(aux3(i,j),Gamma_inf)
            endif	 
         enddo
         enddo
         Gmm_ftn=aux3

         Gmm_ftn(0,0:nyp3) = Gmm_ftn(3,0:nyp3)
         Gmm_ftn(1,0:nyp3) = Gmm_ftn(2,0:nyp3)
         Gmm_ftn(0:nxp3,0) = Gmm_ftn(0:nxp3,2)
         Gmm_ftn(0:nxp3,1) = Gmm_ftn(0:nxp3,2)
      endif
      
      Phi_org=Phi

      rI_val1=-0.3d0*hxy
      Phi=Phi_org+rI_val1
      call Getvolume(Vol_cal1)

      rI_val2=0.29d0*hxy
      Phi=Phi_org+rI_val2
      call Getvolume(Vol_cal2)

      if((Vol_cal1-Vol_ex)*(Vol_cal2-Vol_ex).lt.0.0d0)then
         Vol_cal3=1.0d0
         niter=0
         do while( (dabs((Vol_cal3-Vol_ex)/Vol_ex).gt.1.0d-5 ) 
     &                                  .and. (niter.lt.100) )
            niter=niter+1
            rI_val3=rI_val1+(Vol_ex-vol_cal1)*(rI_val1-rI_val2)
     &                                     /(vol_cal1-vol_cal2)
            Phi=Phi_org+rI_val3
            call Getvolume(Vol_cal3)
 
            if((Vol_cal3-Vol_ex).lt.0.0d0)then
               rI_val1=rI_val3
               vol_cal1=vol_cal3
            else if((Vol_cal3-Vol_ex).gt.0.0d0)then
               rI_val2=rI_val3
               vol_cal2=vol_cal3
            end if
         enddo      

         print *,niter,rI_val3,dabs((Vol_cal3-Vol_ex)/Vol_ex)

      else
         Phi=Phi_org
         call Getvolume(Vol_cal0)
      endif

      if(isurfactant_on.eq.1)then
         Gmm_Xf=0.0d0
         do m=1,nelemento
            x1=ptxo(m,1)
            y1=ptyo(m,1)
            x2=ptxo(m,2)
            y2=ptyo(m,2)
	
            ptxx=(x1+x2)*half
            ptyy=(y1+y2)*half

            io=floor((ptxx+hxh)*hxi)+1
            jo=floor((ptyy+hyh)*hyi)+1
            pp=(ptxx-xs(io))*hxi
            qq=(ptyy-ys(jo))*hyi
            Gmm_Xf(m)=  (one-pp)*(one-qq)*Gmm_ftn(io,jo)
     &                 +    pp*(one-qq)*Gmm_ftn(io+1,jo)
     &                 +    (one-pp)*qq*Gmm_ftn(io,jo+1)
     &                 +        pp*qq*Gmm_ftn(io+1,jo+1)

            Gmm_Xf(m)=dmax1(0.0d0,Gmm_Xf(m))
            Gmm_Xf(m)=dmin1(Gamma_inf,Gmm_Xf(m))
         enddo
      endif
      
      return
      end


c***********************************************************************     
      subroutine Getvolume(Volc)
c***********************************************************************     
      include 'mac3d.inc'

      call reconstruct

      call Find_distance

      Volc=zero
      do j=2,nyp1
      do i=2,nxp1
         Volc=Volc+Heavi(Phi(i,j)/hxy)*hx*hy*2.0d0*pi*xs(i)
      enddo
      enddo

      return
      end


c***********************************************************************     
      subroutine reconstruct
c***********************************************************************     
      include 'mac3d.inc'
      dimension ptx_org(3),pty_org(3)

      nelemento=0
      j_cur=0
      do 100 j=1,ny
         j_cur=j_cur+1
         i_cur=0
         do 100 i=1,nx
            i_cur=i_cur+1

            icr=i
            if(i.lt.0)icr=3-i

            if(ptrec(icr,j).gt.0.5d0)then

            ptx_org(1)=x(i+mod(i_cur+1,2))    
            pty_org(1)=y(j+mod(j_cur,2))
            ptx_org(2)=x(i+mod(mod(i_cur+1,2)+1,2))  
            pty_org(2)=y(j+mod(j_cur,2))
            ptx_org(3)=x(i+mod(i_cur+1,2))    
            pty_org(3)=y(j+mod(mod(j_cur,2)+1,2))
            call Getsurface(ptx_org,pty_org)

            ptx_org(1)=x(i+mod(i_cur,2))
            pty_org(1)=y(j+mod(j_cur+1,2))
            ptx_org(2)=x(i+mod(mod(i_cur,2)+1,2))
            pty_org(2)=y(j+mod(j_cur+1,2))
            ptx_org(3)=x(i+mod(i_cur,2))
            pty_org(3)=y(j+mod(mod(j_cur+1,2)+1,2))
            call Getsurface(ptx_org,pty_org)

            endif

100   continue

      iindx_pt=0
      do m=1,nelemento
         if(ptyo(m,1).eq.0.0d0)then
            iindx_pt(m,1)=1
         endif
         if(ptyo(m,2).eq.0.0d0)then
            iindx_pt(m,2)=1
         endif
      enddo

      return
      end


c***********************************************************************     
      subroutine Getsurface(ptx_org,pty_org)
c***********************************************************************     
      include 'mac3d.inc'
      dimension pta_loc(2),ptb_loc(2)
      dimension ptx_org(3),pty_org(3)

      iline_chk=0
      call surface(ptx_org,pty_org,iline_chk,pta_loc,ptb_loc)
      if(iline_chk.eq.1)then
         x1=pta_loc(1)
         y1=ptb_loc(1)
         x2=pta_loc(2)
         y2=ptb_loc(2)

         rnxx=y2-y1
         rnyy=-(x2-x1)
         rsum=dsqrt(rnxx**two+rnyy**two)
         rnxx=-rnxx/rsum
         rnyy=-rnyy/rsum

         ptxx=(x1+x2)/two
         ptyy=(y1+y2)/two

         rmax_DD=0.0d0
         do m=1,3
            rcur_DD=Valinf(ptx_org(m),pty_org(m))
            if(dabs(rcur_DD).gt.dabs(rmax_DD))then
               mmax_DD=m
               rmax_DD=rcur_DD
            endif
         enddo
         x0=ptx_org(mmax_DD)
         y0=pty_org(mmax_DD)
         DD_sign=sign(1.0d0,rmax_DD)
         rnxc=x0-ptxx
         rnyc=y0-ptyy

         if(DD_sign*(rnxx*rnxc+rnyy*rnyc).lt.zero)then
            nelemento=nelemento+1
            ptxo(nelemento,1)=x2
            ptyo(nelemento,1)=y2
            ptxo(nelemento,2)=x1
            ptyo(nelemento,2)=y1
         else
            nelemento=nelemento+1
            ptxo(nelemento,1)=x1
            ptyo(nelemento,1)=y1
            ptxo(nelemento,2)=x2
            ptyo(nelemento,2)=y2
         endif

      endif

      return
      end


c***********************************************************************     
      subroutine surface(ptx_loc,pty_loc,i_chk,pta_loc,ptb_loc)
c***********************************************************************     
      include 'mac3d.inc'
      dimension ptx_loc(3),pty_loc(3)
      dimension pta_loc(2),ptb_loc(2)

      D11=Valinf(ptx_loc(1),pty_loc(1))
      D22=Valinf(ptx_loc(2),pty_loc(2))
      D33=Valinf(ptx_loc(3),pty_loc(3))

      npoint=0
      if(D11.eq.zero)then
         npoint=npoint+1
         pta_loc(npoint)=ptx_loc(1)
         ptb_loc(npoint)=pty_loc(1)
      endif
      if(D11*D22.lt.zero)then
         npoint=npoint+1
         call Searchpt(ptx_loc(1),pty_loc(1),
     &                 ptx_loc(2),pty_loc(2),xsol,ysol)
         pta_loc(npoint)=xsol
         ptb_loc(npoint)=ysol
      endif
      if(D22.eq.zero)then
         npoint=npoint+1
         pta_loc(npoint)=ptx_loc(2)
         ptb_loc(npoint)=pty_loc(2)
      endif
      if(D22*D33.lt.zero)then
         npoint=npoint+1
         call Searchpt(ptx_loc(2),pty_loc(2),
     &                 ptx_loc(3),pty_loc(3),xsol,ysol)
         pta_loc(npoint)=xsol
         ptb_loc(npoint)=ysol
      endif
      if(D33.eq.zero)then
        npoint=npoint+1
         pta_loc(npoint)=ptx_loc(3)
         ptb_loc(npoint)=pty_loc(3)
      endif
      if(D33*D11.lt.zero)then
         npoint=npoint+1
         call Searchpt(ptx_loc(3),pty_loc(3),
     &                 ptx_loc(1),pty_loc(1),xsol,ysol)
         pta_loc(npoint)=xsol
         ptb_loc(npoint)=ysol
      endif

      if(npoint.eq.2)then
         if( dsqrt( (pta_loc(1)-pta_loc(2))**two
     &             +(ptb_loc(1)-ptb_loc(2))**two ).gt.1.0d-7*hxy)then
            i_chk=1
         endif
      endif

      return
      end


c***********************************************************************     
      subroutine Searchpt(xini,yini,xfin,yfin,xsol,ysol)
c***********************************************************************     
      include 'mac3d.inc'

      x1=xini
      y1=yini
      FimeI1=Valinf(x1,y1)
 
      x2=xfin
      y2=yfin
      FimeI2=Valinf(x2,y2)

100   x3=(x1+x2)*half
      y3=(y1+y2)*half
      FimeI3=Valinf(x3,y3)

      if(dabs(FimeI3).gt.1.0d-7*hxy)then
         if(FimeI3*FimeI2.lt.0.0d0)then
            x1=x3
            y1=y3
            goto 100
         else if(FimeI3*FimeI1.lt.0.0d0)then
            x2=x3
            y2=y3
            goto 100
         endif
      endif

      xsol=x3
      ysol=y3

      return
      end


c***********************************************************************     
      subroutine Search_Xf(xpos,ypos,rtx,rty,xsol,ysol)
c***********************************************************************     
      include 'mac3d.inc'

      ang1=-0.45d0*pi
      rotx=dcos(ang1)*rtx-dsin(ang1)*rty      
      roty=dsin(ang1)*rtx+dcos(ang1)*rty      
      x1=xpos+0.5d0*hxy*rotx
      y1=ypos+0.5d0*hxy*roty
      FimeI1=Valinf(x1,y1)

      ang2=0.45d0*pi
      rotx=dcos(ang2)*rtx-dsin(ang2)*rty      
      roty=dsin(ang2)*rtx+dcos(ang2)*rty      
      x2=xpos+0.5d0*hxy*rotx
      y2=ypos+0.5d0*hxy*roty
      FimeI2=Valinf(x2,y2)

100   ang3=0.5d0*(ang1+ang2)
      rotx=dcos(ang3)*rtx-dsin(ang3)*rty      
      roty=dsin(ang3)*rtx+dcos(ang3)*rty      
      x3=xpos+0.5d0*hxy*rotx
      y3=ypos+0.5d0*hxy*roty
      FimeI3=Valinf(x3,y3)

      if(dabs(FimeI3).gt.1.0d-7*hxy)then
         if(FimeI3*FimeI2.lt.0.0d0)then
            ang1=ang3
            goto 100
         else if(FimeI3*FimeI1.lt.0.0d0)then
            ang2=ang3
            goto 100
         endif
      endif

      xsol=x3
      ysol=y3

      return
      end


c***********************************************************************     
      subroutine Find_distance
c***********************************************************************     
      include 'mac3d.inc'
      dimension Dist(0:maxnx,0:maxny)
      dimension Dists(0:maxnx,0:maxny)
      dimension schk(0:maxnx,0:maxny),ichk(0:maxnx,0:maxny)

      ichk(0:nxp3,0:nyp3)=0
      schk(0:nxp3,0:nyp3)=0.0d0
      Dist(0:nxp3,0:nyp3)=10.0d0*hxy
      Dists(0:nxp3,0:nyp3)=10.0d0*hxy
      ptij(0:nxp3,0:nyp3)=0.0d0
      do 100 m=1,nelemento
         x1=ptxo(m,1)
         y1=ptyo(m,1)
         x2=ptxo(m,2)
         y2=ptyo(m,2)

         ptxxc=(x1+x2)/two
         ptyyc=(y1+y2)/two

         areaT=dsqrt((x1-x2)**two+(y1-y2)**two)

         rnxx=y2-y1
         rnyy=-(x2-x1)
         rsum=dsqrt(rnxx**two+rnyy**two)
         rnxx=-rnxx/rsum
         rnyy=-rnyy/rsum

         rnx12=x2-x1
         rny12=y2-y1
         rsum=dsqrt(rnx12**two+rny12**two)
         rnx12=rnx12/rsum
         rny12=rny12/rsum

         rnx21=x1-x2
         rny21=y1-y2
         rsum=dsqrt(rnx21**two+rny21**two)
         rnx21=rnx21/rsum
         rny21=rny21/rsum

         isf=floor((ptxxc+hxh)*hxi)+1
         jsf=floor((ptyyc+hyh)*hyi)+1
         do 150 j=jsf-1,jsf+2
         do 150 i=isf-1,isf+2
            if(i.lt.0.or.i.gt.nxp3)goto 100
            if(j.lt.0.or.j.gt.nyp3)goto 100
            Dist_p=dabs(Dists(i,j))
            ichk(i,j)=1

            ptrec(i,j)=1.0d0

            x0=xs(i)
            y0=ys(j)

            Dist_q=dsqrt((ptxxc-x0)**two+(ptyyc-y0)**two)
            rnxc=x0-ptxxc
            rnyc=y0-ptyyc
            schk_c=rnxx*rnxc+rnyy*rnyc
            if(Dist_q.lt.Dist_p)then
               Dists(i,j)=Dist_q
               schk(i,j)=schk_c
            endif            
            
150      continue
      
         do 200 j=jsf-2,jsf+3
         do 200 i=isf-2,isf+3
            if(i.lt.0.or.i.gt.nxp3)goto 200
            if(j.lt.0.or.j.gt.nyp3)goto 200
            Dist_o=dabs(Dist(i,j))

            ptij(i,j)=1.0d0

            x0=xs(i)
            y0=ys(j)
            tt=(x0-x1)*rnx12+(y0-y1)*rny12
            xb=x1+rnx12*tt
            yb=y1+rny12*tt
            rnx1b=xb-x1
            rny1b=yb-y1
            rnx2b=xb-x2
            rny2b=yb-y2
            chk12=rnx1b*rnx12+rny1b*rny12
            chk21=rnx2b*rnx21+rny2b*rny21
            if(chk12.ge.zero.and.chk21.ge.zero)then
               Dist_c=dsqrt((x0-xb)**two+(y0-yb)**two)
               if(Dist_c.le.Dist_o)then
                  Dist(i,j)=Dist_c
                  ptx_surf(i,j)=xb
                  pty_surf(i,j)=yb
               endif      
            else
               if(chk12.lt.zero)then
                  xchk=x1
                  ychk=y1 
               else 
                  xchk=x2
                  ychk=y2 
               endif
               Dist_c=dsqrt((x0-xchk)**two+(y0-ychk)**two)
               if(Dist_c.le.Dist_o)then
                  Dist(i,j)=Dist_c
                  ptx_surf(i,j)=xchk
                  pty_surf(i,j)=ychk
               endif
            endif
200      continue
100   continue

      do m=1,600
         do j=1,nyp2
         do i=1,nxp2
            if(ichk(i,j).eq.0)then
               schk_c=0.0d0
               rsum=0.0d0
               do ii=-1,1
               do jj=-1,1
                  if((ii.ne.0).and.(jj.ne.0))then
                     schk_c=schk_c+schk(i+ii,j+jj)
                     rsum=rsum+1.0d0
                  end if
               enddo            
               enddo            
               schk(i,j)=schk_c/rsum
            end if            
         enddo            
         enddo            
         schk(     0,0:nyp3) = schk(     3,0:nyp3)
         schk(     1,0:nyp3) = schk(     2,0:nyp3)
         schk(  nxp2,0:nyp3) =-1.0d0
         schk(  nxp3,0:nyp3) =-1.0d0
         schk(0:nxp3,     0) = schk(0:nxp3,     1)
         schk(0:nxp3,  nyp3) =-1.0d0
      enddo            

      do 600 j=0,nyp3
      do 600 i=0,nxp3
         phi(i,j)=Dist(i,j)*dsign(1.0d0,schk(i,j))
600   continue
      phi(   0,0:nyp3) = phi(   3,0:nyp3)
      phi(   1,0:nyp3) = phi(   2,0:nyp3)
      phi(0:nxp3,nyp3) = phi(0:nxp3,nyp2)

      return
      end


c***********************************************************************     
      subroutine print1
c***********************************************************************     
      include 'mac3d.inc'

      open(11,file='xf'//number//'.m',status='unknown')
      write(11,*) 'mptp=['
      do m=1,nelemento
         x1=ptxo(m,1)
         y1=ptyo(m,1)
         x2=ptxo(m,2)
         y2=ptyo(m,2)
         write(11,905)x1,y1,x2,y2
         write(11,*) ';'
      enddo
      write(11,*) '];'
      write(11,*) 'hold on'            
      write(11,*) 'for i=1:',nelemento,','
      write(11,*) 'x(1)=mptp(i,1);'
      write(11,*) 'y(1)=mptp(i,2);'
      write(11,*) 'x(2)=mptp(i,3);'
      write(11,*) 'y(2)=mptp(i,4);'
      write(11,*) "plot(x,y,'b-');"
      write(11,*) 'end'
      write(11,*) "axis('equal')"
      write(11,680) 'axis([',zero,xl,zero,yl,'])'
      write(11,*) 'box on'
      close(11)
	
	  
      open(8,file='2D_Para'//number//'.vtk',status='unknown')
      write(8,'(a26)')'# vtk DataFile Version 3.1'
      write(8,'(a22)')'uv p data for paraview'
      write(8,'(a5)')'ASCII'
      write(8,'(a24)')'DATASET RECTILINEAR_GRID'
      write(8,'(a10,3(2x,i4))')'DIMENSIONS',nxp1,nyp1,1
      write(8,'(a13,2x,i4,2x,a5)')'X_COORDINATES',nxp1,'float'
      write(8,907) (x(i),i=1,nxp1)
      write(8,'(a13,2x,i4,2x,a5)')'Y_COORDINATES',nyp1,'float'
      write(8,907) (y(j),j=1,nyp1)
      write(8,'(a21)')'Z_COORDINATES 1 float'
      write(8,'(a15)')'0'
      write(8,'(a10,2x,i6)')'POINT_DATA',nxp1*nyp1
      write(8,'(a30)')'SCALARS distance FLOAT'
      write(8,'(a30)')'LOOKUP_TABLE default'
      do j=1,nyp1
         write(8,907) (0.25d0*(phi(i,j)+phi(i+1,j)
     &                         +phi(i,j+1)+phi(i+1,j+1)),i=1,nxp1)
      enddo
      write(8,'(a30)')'SCALARS pressure FLOAT'
      write(8,'(a30)')'LOOKUP_TABLE default'
      do j=1,nyp1
         write(8,907) (0.25d0*(ps(i,j)+ps(i+1,j)
     &                        +ps(i,j+1)+ps(i+1,j+1)),i=1,nxp1)
      enddo
      write(8,'(a30)')'SCALARS concentration FLOAT'
      write(8,'(a30)')'LOOKUP_TABLE default'
      do j=1,nyp1
         write(8,907) (0.25d0*(C(i,j)+C(i+1,j)
     &                         +C(i,j+1)+C(i+1,j+1)),i=1,nxp1)
      enddo
      write(8,'(a30)')'SCALARS micelle FLOAT'
      write(8,'(a30)')'LOOKUP_TABLE default'
      do j=1,nyp1
         write(8,907) (0.25d0*(rM(i,j)+rM(i+1,j)
     &                 +rM(i,j+1)+rM(i+1,j+1)),i=1,nxp1)
      enddo
      write(8,'(a30)')'SCALARS Gamma FLOAT'
      write(8,'(a30)')'LOOKUP_TABLE default'
      do j=1,nyp1
         write(8,907) (0.25d0*(Gmm_ftn(i,j)+Gmm_ftn(i+1,j)
     &                 +Gmm_ftn(i,j+1)+Gmm_ftn(i+1,j+1)),i=1,nxp1)
      enddo
      write(8,'(a22)')'VECTORS velocity FLOAT'
      do j=1,nyp1
      do i=1,nxp1
       write(8,907)0.5d0*(u(i,j)+u(i,j+1)),0.5d0*(v(i,j)+v(i+1,j)),0.0d0
      enddo
      enddo
      close(8) 

      call print_npz(number)

901   format(40(e20.10,2x))
902   format(8(e20.10,2x))
905   format(9(e20.10,2x))
941   format(800(e20.10,a1))
680   format(a7,4(2x,f20.10),2x,3a)
907   format(5000(f20.10,1x))

      return
      end

! DP Change	
c***********************************************************************     
      subroutine print_npz(number_in)
c***********************************************************************     
      include 'mac3d.inc'
      character*3 number_in

      open(99, file='snap'//number_in//'.bin',
     &     status='unknown', form='unformatted', access='stream')

c     --- header: grid sizes (4 integers) ---
      write(99) nx, ny, nxp1, nyp1

c     --- scalars ---
      write(99) time, xl, yl, hx, hy

c     --- 1-D grid arrays (1:nxp1) ---
      write(99) (x(i),  i=1,nxp1)
      write(99) (y(j),  j=1,nyp1)
      write(99) (xs(i), i=1,nxp1)
      write(99) (ys(j), j=1,nyp1)

c     --- 2-D fields (1:nxp1, 1:nyp1) ---
      write(99) ((phi(i,j),    i=1,nxp1), j=1,nyp1)
      write(99) ((ps(i,j),     i=1,nxp1), j=1,nyp1)
      write(99) ((u(i,j),      i=1,nxp1), j=1,nyp1)
      write(99) ((v(i,j),      i=1,nxp1), j=1,nyp1)
      write(99) ((C(i,j),      i=1,nxp1), j=1,nyp1)
      write(99) ((rM(i,j),     i=1,nxp1), j=1,nyp1)
      write(99) ((Gmm_ftn(i,j),i=1,nxp1), j=1,nyp1)

c     --- surface elements ---
      write(99) nelemento
      write(99) ((ptxo(m,k), k=1,2), m=1,nelemento)
      write(99) ((ptyo(m,k), k=1,2), m=1,nelemento)

      close(99)

      return
      end
! DP Change

c***********************************************************************     
      subroutine filenumber
c***********************************************************************     
      include 'mac3d.inc'

	i1=mod(nf,10)
	i2=mod(nf-i1,100)/10
	i3=(nf-i1-i2)/100

      number=char(i3+48)//char(i2+48)//char(i1+48)
	
	return
	end


c***********************************************************************     
      double precision function pes(zr)
c***********************************************************************     
      implicit double precision (a-h,o-z)

      zero=0.0d0
      half=0.5d0
      one=1.0d0
      two=2.0d0
      if (abs(zr).ge.two) then
         pes=zero
         return
      endif
      if (abs(zr).le.one) then
         pes=pes2(zr)
      else
         pes=half-pes2(two-abs(zr))
      endif   

      return
      end


c***********************************************************************     
      double precision function pes2(zr)
c***********************************************************************     
      implicit double precision (a-h,o-z)

      zero=0.0d0
      half=0.5d0
      one=1.0d0
      two=2.0d0
      three=3.0d0
      four=4.0d0
      zr=abs(zr)
      pes2=(three-two*zr+sqrt(one+four*zr-four*zr*zr))/8.0d0

      return
      end


c***********************************************************************     
      double precision function Heavi(zr)
c***********************************************************************     
      implicit double precision (a-h,o-z)

      zero=0.0d0
      half=0.5d0
      one=1.0d0
      two=2.0d0
      three=3.0d0
      four=4.0d0
      pi=four*datan(one)

      if (zr.gt.two) then
         Heavi=one
         return
      else if (zr.lt.-two) then
         Heavi=zero
         return
      else
         Heavi=(zr+two)/four+dsin(pi*zr/two)/two/pi
         return
      endif   

      end


c***********************************************************************     
      subroutine Teval(xeval_in,yeval_in,Tval)
c***********************************************************************     
      include 'mac3d.inc'

      if(xeval_in.lt.zero)then
         xeval=-xeval_in
      else
         xeval=xeval_in
      endif

      if(yeval_in.lt.zero)then
         yeval=0.0d0
      else
         yeval=yeval_in
      endif
      	  
      io=floor((xeval+hxh)*hxi)+1
      jo=floor((yeval+hyh)*hyi)+1

      x1=xs(io)
      y1=ys(jo)
      pp=(xeval-x1)*hxi
      qq=(yeval-y1)*hyi

      Tval=  (one-pp)*(one-qq)*C(io,jo)
     &      +    pp*(one-qq)*C(io+1,jo)
     &      +    (one-pp)*qq*C(io,jo+1)
     &      +        pp*qq*C(io+1,jo+1)
     
      return
      end


c***********************************************************************     
      subroutine Phi_evalc(xeval_in,yeval_in,Deval,Phic)
c***********************************************************************     
      include 'mac3d.inc'
      dimension Phic(0:maxnx,0:maxny)

      if(xeval_in.lt.zero)then
         xeval=-xeval_in
      elseif(xeval_in.gt.xl)then
         xeval=xl
      else
         xeval=xeval_in
      endif

      if(yeval_in.lt.zero)then
         yeval=0.0d0
      else
         yeval=yeval_in
      endif

      isf=floor((xeval+hxh)*hxi)+1
      jsf=floor((yeval+hyh)*hyi)+1

      Deval=zero
      do 100 j=1,4
      do 100 i=1,4
         iisf=isf-2+i
         jjsf=jsf-2+j
         xsf=(xeval-xs(iisf))*hxi
         ysf=(yeval-ys(jjsf))*hyi

         Deval=Deval+fM3(xsf)*fM3(ysf)*Phic(iisf,jjsf)
100   continue

      return
      end


c***********************************************************************     
      subroutine normal(ptxx_in,ptyy_in,rn1,rn2)
c***********************************************************************     
      include 'mac3d.inc'

      if(ptxx_in.lt.zero)then
         ptxx=-ptxx_in
      else
         ptxx=ptxx_in
      endif

      if(ptyy_in.lt.zero)then
         ptyy=0.0d0
      else
         ptyy=ptyy_in
      endif

      rn1=zero
      rn2=zero
      is=int((ptxx+hxh)*hxi)+1
      js=int((ptyy+hyh)*hyi)+1
      do 10 j=1,4
      do 10 i=1,4
         iis=is-2+i
         jjs=js-2+j
         xsf=(ptxx-xs(iis))*hxi
         ysf=(ptyy-ys(jjs))*hyi

         rn1=rn1+d1fM3(xsf)/hx*fM3(ysf)*Phi(iis,jjs)
         rn2=rn2+d1fM3(ysf)/hy*fM3(xsf)*Phi(iis,jjs)
10    continue

      sum=dsqrt(rn1*rn1+rn2*rn2)
      if(ptxx_in.lt.zero)then
         rn1=-rn1/sum
         rn2=rn2/sum
      else
         rn1=rn1/sum
         rn2=rn2/sum
      endif

      return
      end


c***********************************************************************     
      double precision function Valinf(ptxx_in,ptyy_in)
c***********************************************************************     
      include 'mac3d.inc'

      if(ptxx_in.lt.zero)then
         ptxx=-ptxx_in
      else
         ptxx=ptxx_in
      endif

      if(ptyy_in.lt.zero)then
         ptyy=0.0d0
      else
         ptyy=ptyy_in
      endif

      isf=floor((ptxx+hxh)*hxi)+1
      jsf=floor((ptyy+hyh)*hyi)+1

      Valinf=zero
      do 1 j=1,4
      do 1 i=1,4
         iisf=isf-2+i
         jjsf=jsf-2+j
         xsf=(ptxx-xs(iisf))*hxi
         ysf=(ptyy-ys(jjsf))*hyi
  
         Valinf=Valinf+fM3(xsf)*fM3(ysf)*Phi(iisf,jjsf)
1     continue

      return
      end


c***********************************************************************     
      double precision function fM3(zr)
c***********************************************************************     
      implicit double precision (a-h,o-z)

      zero=0.0d0
      half=0.5d0
      one=1.0d0
      two=2.0d0
      x_h=dabs(zr)
      if (x_h.ge.two) then
         fM3=zero
         return
      endif
      if (x_h.le.one) then
         fM3=2.0d0/3.0d0-x_h**2.0d0+x_h**3.0d0/2.0d0
      else
         fM3=(2.0d0-x_h)**3.0d0/6.0d0
      endif   

      return
      end


c***********************************************************************     
      double precision function d1fM3(zr)
c***********************************************************************     
      implicit double precision (a-h,o-z)

      zero=0.0d0
      half=0.5d0
      one=1.0d0
      two=2.0d0
      three=3.0d0

      if(zr.gt.two)then
         d1fM3=0.0d0
      else if(zr.gt.one.and.zr.le.two)then
         d1fM3=-(2.0d0-zr)**2.0d0/2.0d0
      else if(zr.gt.zero.and.zr.le.one)then
         d1fM3=-2.0d0*zr+3.0d0*zr**2.0d0/2.0d0
      else if(zr.gt.-one.and.zr.le.zero)then
         d1fM3=-2.0d0*zr-3.0d0*zr**2.0d0/2.0d0
      else if(zr.gt.-two.and.zr.le.-one)then
         d1fM3=(2.0d0+zr)**2.0d0/2.0d0
      else if(zr.le.-two)then
         d1fM3=0.0d0
      endif   
  
      return
      end


c***********************************************************************     
      subroutine restartin
c***********************************************************************     
      include 'mac3d.inc'

      write(6,*) 'restartin in'
      
      open(21,file='restart1.out',status='old')
      read(21,904) time,Vol_ex,nf,nrestart
      read(21,905) ((ps(i,j),i=0,nxp3),j=0,nyp3)
      read(21,905) ((u(i,j),i=0,nxp3),j=0,nyp3)
      read(21,905) ((v(i,j),i=0,nxp3),j=0,nyp3)
      read(21,905) ((Phi(i,j),i=0,nxp3),j=0,nyp3)
      read(21,905) ((C(i,j),i=0,nxp3),j=0,nyp3)
      read(21,905) ((rM(i,j),i=0,nxp3),j=0,nyp3)
      read(21,905) ((Gmm_ftn(i,j),i=0,nxp3),j=0,nyp3)
      read(21,905) (Cs(i),i=0,nxp3)
      close(21)

      open(19,file='restart2.out',status='old')
      read(19,907)nelemento
      do m=1,nelemento
         read(19,906)ptxo(m,1),ptyo(m,1),
     &               ptxo(m,2),ptyo(m,2),
     &               Gmm_Xf(m),iindx_pt(m,1),iindx_pt(m,2)
      enddo
      close(19) 

      nrestart=nrestart+1
      nrecon=0

904   format(e22.11,e22.11,e22.11,3i10)
905   format(900000e22.11)
906   format(5(e22.11,2x),2(i3,2x))
907   format(i20)

      write(6,*) 'restartin out'

      return
      end


c***********************************************************************     
      subroutine restartout
c***********************************************************************     
      include 'mac3d.inc'

	i1=mod(nrestart,10)
	i2=mod(nrestart-i1,100)/10
	i3=(nrestart-i1-i2)/100

      number=char(i3+48)//char(i2+48)//char(i1+48)

      open(22,file='restart1_'//number//'.out',status='unknown')
      write(22,904) time,Vol_ex,nf,nrestart
      write(22,905) ((ps(i,j),i=0,nxp3),j=0,nyp3)
      write(22,905) ((u(i,j),i=0,nxp3),j=0,nyp3)
      write(22,905) ((v(i,j),i=0,nxp3),j=0,nyp3)
      write(22,905) ((Phi(i,j),i=0,nxp3),j=0,nyp3)
      write(22,905) ((C(i,j),i=0,nxp3),j=0,nyp3)
      write(22,905) ((rM(i,j),i=0,nxp3),j=0,nyp3)
      write(22,905) ((Gmm_ftn(i,j),i=0,nxp3),j=0,nyp3)
      write(22,905) (Cs(i),i=0,nxp3)
      close(22)

      open(24,file='restart2_'//number//'.out',status='unknown')
      write(24,907)nelemento
      do m=1,nelemento
         write(24,906)ptxo(m,1),ptyo(m,1),
     &                ptxo(m,2),ptyo(m,2),
     &                Gmm_Xf(m),iindx_pt(m,1),iindx_pt(m,2)
      enddo
      close(24) 

904   format(e22.11,e22.11,e22.11,3i10)
905   format(900000e22.11)
906   format(5(e22.11,2x),2(i3,2x))
907   format(i20)

      return
      end          