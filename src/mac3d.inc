      implicit double precision (a-h,o-z)
      parameter (maxpt=1000)
      parameter (maxnx=1000)
      parameter (maxny=1000)
      parameter (npad=10)
      logical done
      
      character*3 number

      common/alllog/done

      common/allchar/number

      common/allint/      nf,nx,ny,nxp1,nyp1,nxp2,nyp2,nxp3,nyp3,
     &           nelement,nfout,nrecon,nelemento,irecon,irestart,
     &         nrestart,isurfactant_on,imicelle_on,isubstrate_on

      common/allreal/            xl,yl,hx,hy,hxi,hyi,hxh,hyh,hxy,
     &     pi,pd2,zero,half,one,two,three,four,grx,gry,eps,dsmin,
     &     r1,r2,vs1,vs2,surfT,the_ad,the_re,slip,surf_slip_corr,
     &                 time,dt,pres,Radius,Vol_cur,Vol_ex,restep,
     &                           ptx_max,pty_max,ptx_min,pty_min,
     &                    Dc1,Dc2,Dm1,Dm2,surf_eps,Ds,beta_s,Dgm,
     &       rk_adsorpt,rk_desorpt,rk_m_form,rk_m_break,rk_m_num,
     &         Cs_adsorpt,Cs_desorpt,Gm_Cs_adsorpt,Gm_Cs_desorpt,
     &   C_ini,gamma_inf,gamma_ini,Cs_inf,Cs_ini,V_impact,rM_ini,
     &                           V_test,V_max,CFL,dt_ini, dt_opt 

      common/grid/
     &                   x(-npad:maxnx+npad),y(-npad:maxny+npad),
     &                 xs(-npad:maxnx+npad),ys(-npad:maxny+npad)

      common/prime_variables/                ps(0:maxnx,0:maxny),
     &                    u(0:maxnx,0:maxny),ut(0:maxnx,0:maxny),
     &                    v(0:maxnx,0:maxny),vt(0:maxnx,0:maxny),
     &                   uf(0:maxnx,0:maxny),vf(0:maxnx,0:maxny),
     &                 utt(0:maxnx,0:maxny),vtt(0:maxnx,0:maxny),
     &                    C(0:maxnx,0:maxny),Cp(0:maxnx,0:maxny),
     &                Cpp(0:maxnx,0:maxny),rMpp(0:maxnx,0:maxny),
     &                  rM(0:maxnx,0:maxny),rMp(0:maxnx,0:maxny),
     &                  Cs(0:maxnx),Csn(0:maxnx),Cs_con(0:maxnx)
     
      common/properties/  r(0:maxnx,0:maxny),vs(0:maxnx,0:maxny) 
    
      common/Indicator/                     phi(0:maxnx,0:maxny)

      common/surf_points/            ptxo(maxpt,2),ptyo(maxpt,2),
     &           ptij(0:maxnx,0:maxny),ptx(maxpt,2),pty(maxpt,2),
     &       ptx_surf(0:maxnx,0:maxny),pty_surf(0:maxnx,0:maxny),
     &                 Gmm_Xf(maxpt),Gmm_Ar(maxpt),Gmm_Sc(maxpt),
     &         Gmm_ftn(0:maxnx,0:maxny),Gmm_Src(0:maxnx,0:maxny),
     &                                    ptrec(0:maxnx,0:maxny)

      common/contact_indx/                    iindx_pt(maxpt,2)      
