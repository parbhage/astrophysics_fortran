            program TOV
c
c THIS PROGRAM SOLVES OPENHEIMER-VOLKOFF EQN. FOR NEUTRON STAR
c
         implicit real*8(a-h,o-z)
         dimension prod(9000),rnb(9000),chi(9000),ener(9000),pres(9000)
              open(10,file='TOV.data',status='unknown')
              open(12,file='TOV.out',status='unknown')
c
c ****   READ THE FOLLOWING  ****
c
c    prodi = starting density 9normalised) where mass & radii will be calculated
c  prodf = final density(normalised)  up to which the mass & radii be calculated
c    ngap = interval of points(density) where you want to 
c           calculate the mass & radi
c  ngap = 1 means it will calculate at every density starting from prodi,
c         2 means it will calculate after skipping one density and so on.
c  imaxw = 1/0 with/without maxwellian construction for pressure if dp/drho<0
c  iprof = 1/0 YES/NOT calculate the mass-radius profile for a given density
c** BE AWARE **  put prodi and prodf only from the ndep.dat file
c
c
          read(10,*)prodi,prodf,ngap,imaxw
          read(10,*)index,iprof
ccc
             if(index .eq. 1)then
             open(15,file='MRk120.dat',status='unknown')
c           open(17,file='pe2.dat',status='unknown')
             endif
            if(iprof .eq. 1)then
           open(16,file='k2n1pr.dat',status='unknown')
            endif
ccc
             pie = dacos(-1.0d0)
             dr = 0.001d0          !! km
             dkm = 1.3234d-06      !! conversion of MeV/fm^3 to km^-2
             xmsun = 1.4766d0      !! mass of sun (in km)
c
c  read energy and pressure in MeV/fm^3
c  check the format 300
          open(unit=14,file='eosk120.dat',status='unknown')
             nstart = 0
           do i = 1,9000
            read(14,300)pd,rb,ch,enr,prs
cc          read(14,900)pd,enr,prs
            prod(i) = pd
            rnb(i) = rb
            chi(i) = ch
            ener(i) = enr*dkm      !! energy in km^-2
            pres(i) = prs*dkm      !! press in km^-2
            if(prod(i) .le. prodi)nstart = nstart + 1
            if(prod(i) .ge. prodf)go to 5
           end do
 5           ntot=i
            niter = nstart - ngap
c
        if(imaxw .eq. 0)go to 15
c
c***   Maxwell construction
             nst = 2
             np = nst-1
           do 8 jj = nst,ntot
            if( pres(jj) .ge. pres(jj-1) )then
              np = np + 1
               if(jj .gt. np)then
                pcon = pres(jj-1) + (pdi-pres(jj-1))/2.
                 do k = np-2,1,-1
                  if(pres(k) .le. pcon)go to 6
                  pres(k)=pcon
                 end do
  6               do k = np-1,jj-1
                   pres(k) = pcon
                  end do
                 do kk = jj,ntot
                  if(pres(kk) .ge. pcon)go to 7
                  pres(kk)=pcon
                 end do
 7              nst = kk+1
                np = nst-1
                go to 8
             endif
              pdi = pres(jj)
            endif
 8        continue

c****
c
c  start to calculate the mass and radius of a star at a given ener and prs
c
 15      continue
           niter = niter + ngap
           xmas = 0.0d0            !! mass of the star (initial value)
           rr =  1.0d-05           !!  radius of the star in km (initial value)
           sener = ener(niter)
           spres = pres(niter)
         write(12,*)'ntot,niter,prod(niter),spres,sener'
         write(12,*)ntot,niter,prod(niter),spres/dkm,sener/dkm
 20     continue
           dm = 4.*pie*sener*rr*rr*dr
           dp = - (sener + spres)*(xmas + 4.*pie*rr**3*spres)*dr/rr
     &             /(rr-2.*xmas)
           spres = spres + dp
         if( spres .le. pres(1) )go to 50  !! small press cond. satisfied
c
c finding energy for the new pressure (spres) from P-E curve by interpolation
             do i = niter-1,1,-1
               if(pres(i) .le. spres)go to 40
             end do
 40      sener = ener(i) + (ener(i+1)-ener(i))*(spres-pres(i))
     &                          /(pres(i+1)-pres(i))
           xmas = xmas + dm  !! mass in km
           rr = rr + dr      !! radius in km
c
          write(16,250)sener/dkm,spres/dkm,xmas/xmsun,rr
c
c          write(12,*)'dm,xmas,rr'
c          write(12,*)dm,xmas,rr
c          write(12,*)'dp,spres/dkm,sener/dkm'
c          write(12,*)dp,spres/dkm,sener/dkm
c
             go to 20
 50       xmr = xmas/xmsun
cc        write(15,350)prod(niter),rnb(niter),ener(niter)/dkm,xmr,rr
          write(15,350)prod(niter),ener(niter)/dkm,xmr,rr
          write(12,350)prod(niter),rnb(niter),ener(niter)/dkm,xmr,rr
          write(17,250)prod(niter),rnb(niter),ener(niter)/dkm
     &               ,pres(niter)/dkm
            if( prod(niter+ngap) .ne. 0.0d0 ) go to 15  
c
            write(*,*)'  PROGRAM COMPLETED'
 300    format(5e16.7)
 350    format(5e15.6)
 250    format(4e15.6)
 900     format(7(2x,f9.4))
           stop
            end

