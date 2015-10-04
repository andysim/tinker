c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine induce  --  evaluate induced dipole moments  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "induce" computes the induced dipole moments at polarizable
c     sites due to direct or mutual polarization
c
c     assumes multipole components have already been rotated into
c     the global coordinate frame; computes induced dipoles based
c     on full system, use of active or inactive atoms is ignored
c
c
      subroutine induce
      use sizes
      use inform
      use iounit
      use limits
      use mpole
      use polar
      use potent
      use solute
      use units
      use uprior
      implicit none
      integer i,j,k
      real*8 norm
      logical header
c
c
c     choose the method for computation of induced dipoles
c
      if (solvtyp(1:2) .eq. 'PB') then
         call induce0e
      else if (solvtyp(1:2) .eq. 'GK') then
         call induce0d
      else
         call induce0a
      end if
c
c     update the lists of previous induced dipole values
c
      if (use_pred) then
         nualt = min(nualt+1,maxualt)
         do i = 1, npole
            do j = 1, 3
               do k = nualt, 2, -1
                  udalt(k,j,i) = udalt(k-1,j,i)
                  upalt(k,j,i) = upalt(k-1,j,i)
               end do
               udalt(1,j,i) = uind(j,i)
               upalt(1,j,i) = uinp(j,i)
               if (use_solv) then
                  do k = nualt, 2, -1
                     usalt(k,j,i) = usalt(k-1,j,i)
                     upsalt(k,j,i) = upsalt(k-1,j,i)
                  end do
                  usalt(1,j,i) = uinds(j,i)
                  upsalt(1,j,i) = uinps(j,i)
               end if
            end do
         end do
      end if
c
c     print out a list of the final induced dipole moments
c
      if (debug) then
         header = .true.
         do i = 1, npole
            if (polarity(i) .ne. 0.0d0) then
               if (header) then
                  header = .false.
                  if (solvtyp.eq.'GK' .or. solvtyp.eq.'PB') then
                     write (iout,10)
   10                format (/,' Vacuum Induced Dipole Moments',
     &                          ' (Debyes) :')
                  else
                     write (iout,20)
   20                format (/,' Induced Dipole Moments (Debyes) :')
                  end if
                  if (digits .ge. 8) then
                     write (iout,30)
   30                format (/,4x,'Atom',14x,'X',15x,'Y',15x,'Z',
     &                          15x,'Total',/)
                  else if (digits .ge. 6) then
                     write (iout,40)
   40                format (/,4x,'Atom',14x,'X',13x,'Y',13x,'Z',
     &                          12x,'Total',/)
                  else
                     write (iout,50)
   50                format (/,4x,'Atom',14x,'X',11x,'Y',11x,'Z',
     &                          9x,'Total',/)
                  end if
               end if
               k = ipole(i)
               norm = sqrt(uind(1,i)**2+uind(2,i)**2+uind(3,i)**2)
               if (digits .ge. 8) then
                  write (iout,60)  k,(debye*uind(j,i),j=1,3),debye*norm
   60             format (i8,3x,4f16.8)
               else if (digits .ge. 6) then
                  write (iout,70)  k,(debye*uind(j,i),j=1,3),debye*norm
   70             format (i8,4x,4f14.6)
               else
                  write (iout,80)  k,(debye*uind(j,i),j=1,3),debye*norm
   80             format (i8,5x,4f12.4)
               end if
            end if
         end do
         header = .true.
         if (solvtyp.eq.'GK' .or. solvtyp.eq.'PB') then
            do i = 1, npole
               if (polarity(i) .ne. 0.0d0) then
                  if (header) then
                     header = .false.
                     write (iout,90)
   90                format (/,' SCRF Induced Dipole Moments',
     &                          ' (Debyes) :')
                     if (digits .ge. 8) then
                        write (iout,100)
  100                   format (/,4x,'Atom',14x,'X',15x,'Y',15x,'Z',
     &                             15x,'Total',/)
                     else if (digits .ge. 6) then
                        write (iout,110)
  110                   format (/,4x,'Atom',14x,'X',13x,'Y',13x,'Z',
     &                             12x,'Total',/)
                     else
                        write (iout,120)
  120                   format (/,4x,'Atom',14x,'X',11x,'Y',11x,'Z',
     &                             9x,'Total',/)
                     end if
                  end if
                  k = ipole(i)
                  norm = sqrt(uinds(1,i)**2+uinds(2,i)**2+uinds(3,i)**2)
                  if (digits .ge. 8) then
                     write (iout,130)  k,(debye*uinds(j,i),j=1,3),
     &                                 debye*norm
  130                format (i8,3x,4f16.8)
                  else if (digits .ge. 6) then
                     write (iout,140)  k,(debye*uinds(j,i),j=1,3),
     &                                 debye*norm
  140                format (i8,4x,4f14.6)
                  else
                     write (iout,150)  k,(debye*uinds(j,i),j=1,3),
     &                                 debye*norm
  150                format (i8,5x,4f12.4)
                  end if
               end if
            end do
         end if
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine induce0a  --  conjugate gradient dipole solver  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "induce0a" computes the induced dipole moments at polarizable
c     sites using a preconditioned conjugate gradient solver
c
c
      subroutine induce0a
      use sizes
      use atoms
      use inform
      use iounit
      use limits
      use mpole
      use polar
      use polpot
      use potent
      use units
      use uprior
      use files
      implicit none
      integer i,j,k,iter
c OPT IMPLEMENTATION
      integer ptord
c OPT IMPLEMENTATION
      integer maxiter
      real*8 polmin
      real*8 eps,epsold
      real*8 epsd,epsp
      real*8 udsum,upsum
      real*8 a,ap,b,bp
      real*8 sum,sump
      real*8, allocatable :: poli(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8, allocatable :: udir(:,:)
      real*8, allocatable :: udirp(:,:)
      real*8, allocatable :: rsd(:,:)
      real*8, allocatable :: rsdp(:,:)
      real*8, allocatable :: zrsd(:,:)
      real*8, allocatable :: zrsdp(:,:)
      real*8, allocatable :: conj(:,:)
      real*8, allocatable :: conjp(:,:)
      real*8, allocatable :: vec(:,:)
      real*8, allocatable :: vecp(:,:)
c OPT IMPLEMENTATION
      integer ioff, addr, imin, fitord, imaxd, imaxp
      real*8 csum, erd(3), erp(3), sumd, resmin, grdmin
      real*8 maxd, maxp, error
      real*8, allocatable :: ptd(:,:), ptp(:,:), errors(:)
      integer ipdb
      integer freeunit
      character*120 pdbfile
c OPT IMPLEMENTATION
      logical done
      character*6 mode
c OPT IMPLEMENTATION
      external savestate
      external computegrad
      external computehess
c OPT IMPLEMENTATION
c
c
c     zero out the induced dipoles at each site
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
         end do
      end do
      if (.not. use_polar)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (field(3,npole))
      allocate (fieldp(3,npole))
      allocate (udir(3,npole))
      allocate (udirp(3,npole))

c
c     get the electrostatic field due to permanent multipoles
c
      if (use_ewald) then
         call dfield0c (field,fieldp)
      else if (use_mlist) then
         call dfield0b (field,fieldp)
      else
         call dfield0a (field,fieldp)
      end if
c
c     set induced dipoles to polarizability times direct field
c
      do i = 1, npole
         do j = 1, 3
            udir(j,i) = polarity(i) * field(j,i)
            udirp(j,i) = polarity(i) * fieldp(j,i)
            uind(j,i) = udir(j,i)
            uinp(j,i) = udirp(j,i)
         end do
      end do

c OPT IMPLEMENTATION
      if (poltyp(1:3) .eq. 'OPT') then
         ptuind(:,:,0) = uind
         ptuinp(:,:,0) = uinp

         do ptord = 1, ptmaxord
           uind = ptuind(:,:,ptord-1)
           uinp = ptuinp(:,:,ptord-1)
           ptpointer = ptord-1
           if (use_ewald) then
              call ufield0c (field,fieldp)
           else if (use_mlist) then
              call ufield0b (field,fieldp)
           else
              call ufield0a (field,fieldp)
           end if
           do i = 1, npole
              ptuind(:,i,ptord) = polarity(i)*field(:,i)
              ptuinp(:,i,ptord) = polarity(i)*fieldp(:,i)
           enddo
         enddo

         uind = 0d0
         uinp = 0d0
         do ptord = 0, ptmaxord
           uind = uind + ptcoefs(ptord) * ptuind(:,:,ptord)
           uinp = uinp + ptcoefs(ptord) * ptuinp(:,:,ptord)
           !write(*,*) "PT(", ptord, ")"
           !write(*,*) "   UIND"
           !do i = 1, npole
           !  write(*,'(3F16.10)') ptuind(:,i,ptord)
           !enddo
           !write(*,*) "   UINP"
           !do i = 1, npole
           !  write(*,'(3F16.10)') ptuinp(:,i,ptord)
           !enddo
         enddo
         !write(*,*) "UIND"
         !do i = 1,npole
         !   write(*,'(3F16.10)') uind(:,i)
         !enddo
         !write(*,*) "UINP"
         !do i = 1,npole
         !   write(*,'(3F16.10)') uinp(:,i)
         !enddo
      endif
c
c     set tolerances for computation of mutual induced dipoles
c
c OPT IMPLEMENTATION
      if (poltyp .eq. 'MUTUAL' .or. dofit) then
c OPT IMPLEMENTATION
         ptpointer = -1
c OPT IMPLEMENTATION
         done = .false.
         maxiter = 500
         iter = 0
         polmin = 0.00000001d0
         eps = 100.0d0
c
c     estimated induced dipoles from polynomial predictor
c
         if (use_pred .and. nualt.eq.maxualt) then
            call ulspred
            do i = 1, npole
               do j = 1, 3
                  udsum = 0.0d0
                  upsum = 0.0d0
                  do k = 1, nualt-1
                     udsum = udsum + bpred(k)*udalt(k,j,i)
                     upsum = upsum + bpredp(k)*upalt(k,j,i)
                  end do
                  uind(j,i) = udsum
                  uinp(j,i) = upsum
               end do
            end do
         end if
c
c     perform dynamic allocation of some local arrays
c
         allocate (poli(npole))
         allocate (rsd(3,npole))
         allocate (rsdp(3,npole))
         allocate (zrsd(3,npole))
         allocate (zrsdp(3,npole))
         allocate (conj(3,npole))
         allocate (conjp(3,npole))
         allocate (vec(3,npole))
         allocate (vecp(3,npole))
c
c     get the electrostatic field due to induced dipoles
c
         if (use_ewald) then
            call ufield0c (field,fieldp)
         else if (use_mlist) then
            call ufield0b (field,fieldp)
         else
            call ufield0a (field,fieldp)
         end if
c
c     set initial conjugate gradient residual and conjugate vector
c
         do i = 1, npole
            poli(i) = max(polmin,polarity(i))
            do j = 1, 3
               rsd(j,i) = (udir(j,i)-uind(j,i))/poli(i)
     &                       + field(j,i)
               rsdp(j,i) = (udirp(j,i)-uinp(j,i))/poli(i)
     &                       + fieldp(j,i)
            end do
         end do
         mode = 'BUILD'
         if (use_mlist) then
            call uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
            mode = 'APPLY'
            call uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
         else
            call uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
            mode = 'APPLY'
            call uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
         end if
         do i = 1, npole
            do j = 1, 3
               conj(j,i) = zrsd(j,i)
               conjp(j,i) = zrsdp(j,i)
            end do
         end do
c
c     conjugate gradient iteration of the mutual induced dipoles
c
         do while (.not. done)
            iter = iter + 1
            do i = 1, npole
               do j = 1, 3
                  vec(j,i) = uind(j,i)
                  vecp(j,i) = uinp(j,i)
                  uind(j,i) = conj(j,i)
                  uinp(j,i) = conjp(j,i)
               end do
            end do
            if (use_ewald) then
               call ufield0c (field,fieldp)
            else if (use_mlist) then
               call ufield0b (field,fieldp)
            else
               call ufield0a (field,fieldp)
            end if
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = vec(j,i)
                  uinp(j,i) = vecp(j,i)
                  vec(j,i) = conj(j,i)/poli(i) - field(j,i)
                  vecp(j,i) = conjp(j,i)/poli(i) - fieldp(j,i)
               end do
            end do
            a = 0.0d0
            ap = 0.0d0
            sum = 0.0d0
            sump = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  a = a + conj(j,i)*vec(j,i)
                  ap = ap + conjp(j,i)*vecp(j,i)
                  sum = sum + rsd(j,i)*zrsd(j,i)
                  sump = sump + rsdp(j,i)*zrsdp(j,i)
               end do
            end do
            if (a .ne. 0.0d0)  a = sum / a
            if (ap .ne. 0.0d0)  ap = sump / ap
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = uind(j,i) + a*conj(j,i)
                  uinp(j,i) = uinp(j,i) + ap*conjp(j,i)
                  rsd(j,i) = rsd(j,i) - a*vec(j,i)
                  rsdp(j,i) = rsdp(j,i) - ap*vecp(j,i)
               end do
            end do
            if (use_mlist) then
               call uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
            else
               call uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
            end if
            b = 0.0d0
            bp = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  b = b + rsd(j,i)*zrsd(j,i)
                  bp = bp + rsdp(j,i)*zrsdp(j,i)
               end do
            end do
            if (sum .ne. 0.0d0)  b = b / sum
            if (sump .ne. 0.0d0)  bp = bp / sump
            epsd = 0.0d0
            epsp = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  conj(j,i) = zrsd(j,i) + b*conj(j,i)
                  conjp(j,i) = zrsdp(j,i) + bp*conjp(j,i)
                  epsd = epsd + rsd(j,i)*rsd(j,i)
                  epsp = epsp + rsdp(j,i)*rsdp(j,i)
               end do
            end do
c
c     check the convergence of the mutual induced dipoles
c
            epsold = eps
            eps = max(epsd,epsp)
            eps = debye * sqrt(eps/dble(npolar))
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. politer)  done = .true.
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (poli)
         deallocate (rsd)
         deallocate (rsdp)
         deallocate (zrsd)
         deallocate (zrsdp)
         deallocate (conj)
         deallocate (conjp)
         deallocate (vec)
         deallocate (vecp)
c
c     print the results from the conjugate gradient iteration
c
         if (debug) then
            write (iout,30)  iter,eps
   30       format (/,' Induced Dipoles :',6x,'Iterations',i5,
     &                 6x,'RMS Change',f15.10)
         end if
c
c     terminate the calculation if dipoles failed to converge
c
         if (iter.ge.maxiter .or. eps.gt.epsold) then
            write (iout,40)
   40       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
            call prterr
            call fatal
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (field)
      deallocate (fieldp)
      deallocate (udir)
      deallocate (udirp)
c OPT IMPLEMENTATION
      if (dofit) then
         allocate(ptd(3,npole))
         allocate(ptp(3,npole))
         allocate(errors(npole))
         do fitord = 2,6,2
           grdmin = (1d-10)**fitord
           fitpow = fitord
           write(*,*) "Results for residual raised to power of",fitpow
           call tncg('DTNCG ', 'NONE  ', ptmaxord+1, ptcoefs, resmin, 
           !call tncg('TNCG  ', 'DIAG  ', ptmaxord+1, ptcoefs, resmin, 
     &               grdmin, computegrad, computehess, savestate)

           write(*,*)
           write(*,*) "Quality of PT induced dipole fit:-"
           write(*,*)
           write(*,'(A5,A21,A20,A20)') " Atom", "UIND", "", "UINP"
           write(*,'(A8,A8,A10,A14,A8,A8,A10,A14)')
     &         " ", "Exact", "PT", "Error"," ", "Exact", "PT", "Error"
           ptd = 0d0
           ptp = 0d0
           do j = 0,ptmaxord
             do i = 1,npole
               ptd(:,i) = ptd(:,i) + ptcoefs(j)*ptuind(:,i,j)
               ptp(:,i) = ptp(:,i) + ptcoefs(j)*ptuinp(:,i,j)
             enddo
           enddo
           sumd = 0d0
           sump = 0d0
           maxd = 0.d0
           maxp = 0.d0
           do i = 1,npole
              erd = ptd(:,i)-uind(:,i)
              error = erd(1)**2 + erd(2)**2 + erd(3)**2
              errors(i) = 100.0*sqrt(error)*debye
              sumd = sumd + error
              if(error .gt. maxd) then
                imaxd = i
                maxd = error
              endif
              erp = ptp(:,i)-uinp(:,i)
              error = erp(1)**2 + erp(2)**2 + erp(3)**2
              sump = sump + error
              if(error .gt. maxp) then
                imaxp = i
                maxp = error
              endif
              write(*,'(I5,A1,3F12.8,A4,3F12.8)')
     &                i,"X",uind(1,i)*debye,ptd(1,i)*debye,erd(1)*debye,
     &                   "",uinp(1,i)*debye,ptp(1,i)*debye,erp(1)*debye
              write(*,'(I5,A1,3F12.8,A4,3F12.8)')
     &                i,"Y",uind(2,i)*debye,ptd(2,i)*debye,erd(2)*debye,
     &                   "",uinp(2,i)*debye,ptp(2,i)*debye,erp(2)*debye
              write(*,'(I5,A1,3F12.8,A4,3F12.8)')
     &                i,"Z",uind(3,i)*debye,ptd(3,i)*debye,erd(3)*debye,
     &                   "",uinp(3,i)*debye,ptp(3,i)*debye,erp(3)*debye
           enddo
           maxd = sqrt(maxd)*debye
           maxp = sqrt(maxp)*debye
           write(*,*) 
           write(*,'(A16,F12.10,A16,F12.10,A5,I5)') "Uind RMS error: ",
     &       sqrt(sumd/(3d0*npole))*debye, "Max: ", maxd, " for ", imaxd
           write(*,'(A16,F12.10,A16,F12.10,A5,I5)') "Uinp RMS error: ",
     &       sqrt(sump/(3d0*npole))*debye, "Max: ", maxp, " for ", imaxp
           write(*,*) 
           write(*,'(A23,20F8.4)')  "Optimized coefficients: ",
     &                                    ptcoefsf(0:ptmaxord)
           write(*,*)  
           csum = 0d0
           do i = 0,ptmaxord
             csum = csum + ptcoefsf(i)
           enddo
           write(*,'(A20,F8.4)')  "Sum of coefficients:", csum
           write(*,*) 
c
c          open the Protein Data Bank file to be used for output
c
           ipdb = freeunit ()
           pdbfile = filename(1:leng)//'.pdb'
           call version (pdbfile,'new')
           open (unit=ipdb,file=pdbfile,status='new')
           call makepdb
           call prtpdb2 (ipdb,errors)
           close (unit=ipdb)
         enddo
         deallocate(errors)
         deallocate(ptd)
         deallocate(ptp)
         stop
      endif
c OPT IMPLEMENTATION
      return
      end
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prtpdb  --  output of Protein Data Bank file  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtpdb" writes out a set of Protein Data Bank coordinates
c     to an external disk file
c
c
      subroutine prtpdb2 (ipdb,bfacs)
      use sizes
      use files
      use pdb
      use sequen
      use titles
      implicit none
      integer i,k,ipdb
      integer start,stop
      integer resmax,resnumb
      integer, allocatable :: resid(:)
      real*8 crdmin,crdmax
      real*8 bfacs(*)
      logical opened
      logical rename
      logical reformat
      character*1 chnname
      character*1, allocatable :: chain(:)
      character*2 atmc,resc
      character*3 resname
      character*6 crdc
      character*38 fstr
      character*120 pdbfile
c
c
c     set flags for residue naming and large value formatting
c
      rename = .false.
      reformat = .true.
c
c     open the output unit if not already done
c
      inquire (unit=ipdb,opened=opened)
      if (.not. opened) then
         pdbfile = filename(1:leng)//'.pdb'
         call version (pdbfile,'new')
         open (unit=ipdb,file=pdbfile,status='new')
      end if
c
c     write out the header lines and the title
c
      if (ltitle .eq. 0) then
         fstr = '(''HEADER'',/,''COMPND'',/,''SOURCE'')'
         write (ipdb,fstr(1:32))
      else
         fstr = '(''HEADER'',4x,a,/,''COMPND'',/,''SOURCE'')'
         write (ipdb,fstr(1:37))  title(1:ltitle)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (resid(maxres))
      allocate (chain(maxres))
c
c     find the chain name and chain position for each residue
c
      do i = 1, nchain
         start = ichain(1,i)
         stop = ichain(2,i)
         do k = start, stop
            resid(k) = k - start + 1
            chain(k) = chnnam(i)
         end do
      end do
c
c     change some TINKER residue names to match PDB standards
c
      if (rename) then
         do i = 1, npdb
            if (pdbres(i) .eq. 'CYX')  pdbres(i) = 'CYS'
            if (pdbres(i) .eq. 'CYD')  pdbres(i) = 'CYS'
            if (pdbres(i) .eq. 'TYD')  pdbres(i) = 'TYR'
            if (pdbres(i) .eq. 'HID')  pdbres(i) = 'HIS'
            if (pdbres(i) .eq. 'HIE')  pdbres(i) = 'HIS'
            if (pdbres(i) .eq. 'HIP')  pdbres(i) = 'HIS'
            if (pdbres(i) .eq. 'ASH')  pdbres(i) = 'ASP'
            if (pdbres(i) .eq. 'GLH')  pdbres(i) = 'GLU'
            if (pdbres(i) .eq. 'LYD')  pdbres(i) = 'LYS'
         end do
      end if
c
c     set formatting to match the PDB fixed format standard
c
      atmc = 'i5'
      resc = 'i4'
      crdc = '3f8.3'
c
c     check for large values requiring extended formatting
c
      if (reformat) then
         resmax = 0
         crdmin = 0.0d0
         crdmax = 0.0d0
         do i = 1, npdb
            if (pdbtyp(i) .eq. 'ATOM  ') then
               resmax = max(resmax,resid(resnum(i)))
            else
               resmax = max(resmax,resnum(i))
            end if
            crdmin = min(crdmin,xpdb(i),ypdb(i),zpdb(i))
            crdmax = max(crdmax,xpdb(i),ypdb(i),zpdb(i))
         end do
         if (npdb .ge. 100000)  atmc = 'i6'
         if (resmax .ge. 10000)  resc = 'i5'
         if (resmax .ge. 100000)  resc = 'i6'
         if (crdmin .le. -100.0d0)  crdc = '3f9.3 '
         if (crdmax .ge. 1000.0d0)  crdc = '3f9.3 '
         if (crdmin .le. -1000.0d0)  crdc = '3f10.3'
         if (crdmax .ge. 10000.0d0)  crdc = '3f10.3'
      end if
c
c     write info and coordinates for each PDB atom
c
      fstr = '(a6,'//atmc//',1x,a4,1x,a3,1x,a1,'//resc//
     &          ',4x,'//crdc//',f6.2)'
      do i = 1, npdb
         resname = pdbres(i)
         if (resname(2:3) .eq. '  ')  resname = '  '//resname(1:1)
         if (resname(3:3) .eq. ' ')  resname = ' '//resname(1:2)
         if (pdbtyp(i) .eq. 'ATOM  ') then
            resnumb = resid(resnum(i))
            chnname = chain(resnum(i))
         else
            resnumb = resnum(i)
            chnname = ' '
         end if
         write (ipdb,'(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3,f6.2)')
     &                      pdbtyp(i),i,pdbatm(i),resname,chnname,
     &                      resnumb,xpdb(i),ypdb(i),zpdb(i),bfacs(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (resid)
      deallocate (chain)
c
c     check for large values requiring extended formatting
c
      if (reformat) then
         if (npdb .ge. 100000)  atmc = 'i7'
         if (npdb .ge. 10000)  atmc = 'i6'
      end if
c
c     write any connectivity records for PDB atoms
c
      fstr = '(''CONECT'',9'//atmc//')'
      do i = 1, npdb
         if (npdb12(i) .ne. 0) then
            write (ipdb,fstr(1:14))  i,(ipdb12(k,i),k=1,npdb12(i))
         end if
      end do
      fstr = '(''END'')'
      write (ipdb,fstr(1:7))
c
c     close the output unit if opened by this routine
c
c     if (.not. opened)  close (unit=ipdb)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine makepdb  --  build PDB from Cartesian coords  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "makepdb" cconstructs a Protein Data Bank file from a set
c     of Cartesian coordinates with special handling for systems
c     consisting of biopolymer chains, ligands and water molecules
c
c
      subroutine makepdb
      use sizes
      use atomid
      use atoms
      use couple
      use files
      use molcul
      use pdb
      use resdue
      use sequen
      implicit none
      integer i,j,k,m,ii
      integer kp,ka,kn
      integer iseq,freeunit
      integer start,stop
      integer pdbnum,atmnum
      integer justify,cbi
      integer noxy,nhydro
      integer, allocatable :: ni(:)
      integer, allocatable :: cai(:)
      integer, allocatable :: ci(:)
      integer, allocatable :: oi(:)
      integer, allocatable :: poi(:)
      integer, allocatable :: op1(:)
      integer, allocatable :: op2(:)
      integer, allocatable :: op3(:)
      integer, allocatable :: c5i(:)
      integer, allocatable :: o5i(:)
      integer, allocatable :: c4i(:)
      integer, allocatable :: o4i(:)
      integer, allocatable :: c3i(:)
      integer, allocatable :: o3i(:)
      integer, allocatable :: c2i(:)
      integer, allocatable :: o2i(:)
      integer, allocatable :: c1i(:)
      logical exist,generic
      logical cbone,nbone,obone
      logical first
      logical, allocatable :: water(:)
      logical, allocatable :: hetmol(:)
      character*3 resname
      character*4 atmname
      character*7, allocatable :: restyp(:)
      character*120 seqfile
      save first
      data first  / .true. /
c
c
c     perform dynamic allocation of some global arrays
c
      if (first) then
         first = .false.
         if (.not. allocated(resnum))  allocate (resnum(maxatm))
         if (.not. allocated(resatm))  allocate (resatm(2,maxatm))
         if (.not. allocated(npdb12))  allocate (npdb12(maxatm))
         if (.not. allocated(ipdb12))  allocate (ipdb12(maxval,maxatm))
         if (.not. allocated(pdblist))  allocate (pdblist(maxatm))
         if (.not. allocated(xpdb))  allocate (xpdb(maxatm))
         if (.not. allocated(ypdb))  allocate (ypdb(maxatm))
         if (.not. allocated(zpdb))  allocate (zpdb(maxatm))
         if (.not. allocated(pdbres))  allocate (pdbres(maxatm))
         if (.not. allocated(pdbatm))  allocate (pdbatm(maxatm))
         if (.not. allocated(pdbtyp))  allocate (pdbtyp(maxatm))
      end if
c
c     initialize number of PDB atoms and atom mapping
c
      npdb = 0
      do i = 1, n
         pdblist(i) = 0
      end do
c
c     read the biopolymer sequence file if one exists
c
      iseq = freeunit ()
      seqfile = filename(1:leng)//'.seq'
      call version (seqfile,'old')
      inquire (file=seqfile,exist=exist)
      if (exist) then
         open (unit=iseq,file=seqfile,status='old')
         rewind (unit=iseq)
         call readseq (iseq)
         close (unit=iseq)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (ni(maxres))
      allocate (cai(maxres))
      allocate (ci(maxres))
      allocate (oi(maxres))
      allocate (poi(maxres))
      allocate (op1(maxres))
      allocate (op2(maxres))
      allocate (op3(maxres))
      allocate (c5i(maxres))
      allocate (o5i(maxres))
      allocate (c4i(maxres))
      allocate (o4i(maxres))
      allocate (c3i(maxres))
      allocate (o3i(maxres))
      allocate (c2i(maxres))
      allocate (o2i(maxres))
      allocate (c1i(maxres))
      allocate (restyp(maxres))
c
c     zero out the backbone atoms for biopolymer sequences
c
      do i = 1, nseq
         ni(i) = 0
         cai(i) = 0
         ci(i) = 0
         oi(i) = 0
         poi(i) = 0
         op1(i) = 0
         op2(i) = 0
         op3(i) = 0
         c5i(i) = 0
         o5i(i) = 0
         c4i(i) = 0
         o4i(i) = 0
         c3i(i) = 0
         o3i(i) = 0
         c2i(i) = 0
         o2i(i) = 0
         c1i(i) = 0
      end do
c
c     set the molecule type for each residue via chain type
c
      generic = .true.
      do i = 1, nchain
         do j = ichain(1,i), ichain(2,i)
            restyp(j) = 'GENERIC'
            if (chntyp(i) .eq. 'PEPTIDE')  restyp(j) = 'PEPTIDE'
            if (chntyp(i) .eq. 'NUCLEIC')  restyp(j) = 'NUCLEIC'
         end do
         if (restyp(j)  .ne. 'GENERIC')  generic = .false.
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (water(nmol))
c
c     check each molecule to see if it is a water molecule
c
      do i = 1, nmol
         water(i) = .false.
         if (imol(2,i)-imol(1,i) .eq. 2) then
            noxy = 0
            nhydro = 0
            do j = imol(1,i), imol(2,i)
               k = kmol(j)
               if (atomic(k) .eq. 8)  noxy = noxy + 1
               if (atomic(k) .eq. 1)  nhydro = nhydro + 1
            end do
            if (noxy.eq.1 .and. nhydro.eq.2)  water(i) = .true.
         end if
      end do
c
c     for general structures, transfer each atom to PDB format
c
      if (generic) then
         do i = 1, nmol
            do j = imol(1,i), imol(2,i)
               k = kmol(j)
               atmname = ' '//name(k)
               if (water(i)) then
                  resname = 'HOH'
               else
                  justify = 0
                  call numeral (type(k),resname,justify)
               end if
               pdbnum = i
               call pdbatom (atmname,resname,pdbnum,k)
               pdbtyp(npdb) = 'HETATM'
            end do
         end do
         do i = 1, nmol
            do j = imol(1,i), imol(2,i)
               k = kmol(j)
               kp = pdblist(k)
               npdb12(kp) = n12(k)
               do m = 1, n12(k)
                  ipdb12(m,kp) = pdblist(i12(m,k))
               end do
            end do
         end do
      end if
c
c     find the amide nitrogens and other peptide backbone atoms
c
      m = 1
      do i = 1, n
         if (restyp(m) .eq. 'PEPTIDE') then
            resname = amino(seqtyp(m))
            if (resname .eq. 'H2N') then
               m = m + 1
            else if (resname .eq. 'FOR') then
               if (atomic(i) .eq. 6) then
                  nbone = .false.
                  obone = .false.
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k) .eq. 7) then
                        nbone = .true.
                     else if (atomic(k) .eq. 8) then
                        obone = .true.
                     end if
                  end do
                  if (nbone .and. obone) then
                     cai(m) = i
                     ci(m) = i
                     oi(m) = i + 1
                     m = m + 1
                  end if
               end if
            else if (resname .eq. 'ACE') then
               if (atomic(i) .eq. 6) then
                  nbone = .false.
                  obone = .false.
                  do j = 1, n13(i)
                     k = i13(j,i)
                     if (atomic(k) .eq. 7) then
                        nbone = .true.
                     else if (atomic(k) .eq. 8) then
                        obone = .true.
                     end if
                  end do
                  if (nbone .and. obone) then
                     cai(m) = i
                     ci(m) = i + 1
                     oi(m) = i + 2
                     m = m + 1
                  end if
               end if
            else if (resname .eq. 'COH') then
               if (n12(i) .gt. 1) then
                  if (atomic(i) .eq. 8) then
                     nbone = .false.
                     obone = .false.
                     do j = 1, n13(i)
                        k = i13(j,i)
                        if (atomic(k) .eq. 8) then
                           obone = .true.
                        end if
                     end do
                     do j = 1, n14(i)
                        k = i14(j,i)
                        if (atomic(k) .eq. 7) then
                           nbone = .true.
                        end if
                     end do
                     if (nbone .and. obone) then
                        ni(m) = i
                        m = m + 1
                     end if
                  end if
               end if
            else if (resname .eq. 'NH2') then
               if (atomic(i) .eq. 7) then
                  nbone = .false.
                  obone = .false.
                  do j = 1, n13(i)
                     k = i13(j,i)
                     if (atomic(k) .eq. 8) then
                        obone = .true.
                     end if
                  end do
                  do j = 1, n14(i)
                     k = i14(j,i)
                     if (atomic(k) .eq. 7) then
                        nbone = .true.
                     end if
                  end do
                  if (nbone .and. obone) then
                     ni(m) = i
                     m = m + 1
                  end if
               end if
            else if (resname .eq. 'NME') then
               if (atomic(i) .eq. 7) then
                  nbone = .false.
                  obone = .false.
                  do j = 1, n13(i)
                     k = i13(j,i)
                     if (atomic(k) .eq. 8) then
                        obone = .true.
                     end if
                  end do
                  do j = 1, n14(i)
                     k = i14(j,i)
                     if (atomic(k) .eq. 7) then
                        nbone = .true.
                     end if
                  end do
                  if (nbone .and. obone) then
                     ni(m) = i
                     cai(m) = i + 1
                     m = m + 1
                  end if
               end if
            else
               if (atomic(i) .eq. 7) then
                  obone = .false.
                  do j = 1, n14(i)
                     k = i14(j,i)
                     if (atomic(k) .eq. 8) then
                        obone = .true.
                     end if
                  end do
                  if (obone) then
                     ni(m) = i
                     cai(m) = i + 1
                     ci(m) = i + 2
                     oi(m) = i + 3
                     m = m + 1
                  end if
               end if
            end if
c
c     find the phosphates and sugar C1 nucleotide backbone atoms
c
         else if (restyp(m) .eq. 'NUCLEIC') then
            resname = nuclz(seqtyp(m))
            if (resname .eq. 'MP ') then
               if (atomic(i) .eq. 15) then
                  poi(m) = i
                  m = m + 1
               end if
            end if
            if (atomic(i).eq.6 .and. n12(i).eq.4) then
               cbone = .false.
               nbone = .false.
               obone = .false.
               do j = 1, n12(i)
                  k = i12(j,i)
                  ka = atomic(k)
                  kn = n12(k)
                  if (ka .eq. 6)  cbone = .true.
                  if (ka.eq.7 .and. kn.eq.3)  nbone = .true.
                  if (ka.eq.8 .and. kn.eq.2)  obone = .true.
               end do
               if (cbone .and. nbone .and. obone) then
                  c1i(m) = i
                  m = m + 1
               end if
            end if
         end if
         if (m .gt. nseq)  goto 10
      end do
   10 continue
c
c     find the remainder of the nucleotide backbone atoms
c
      do ii = 1, nchain
         if (chntyp(ii) .eq. 'NUCLEIC') then
            start = ichain(1,ii)
            stop = ichain(2,ii)
            do i = start, stop
               m = c1i(i)
               if (m .ne. 0) then
                  do j = 1, n12(m)
                     k = i12(j,m)
                     ka = atomic(k)
                     if (ka .eq. 6)  c2i(i) = k
                     if (ka .eq. 7)  ni(i) = k
                     if (ka .eq. 8)  o4i(i) = k
                  end do
                  m = o4i(i)
                  do j = 1, n12(m)
                     k = i12(j,m)
                     ka = atomic(k)
                     if (ka.eq.6 .and. k.ne.c1i(i))  c4i(i) = k
                  end do
                  m = c2i(i)
                  do j = 1, n12(m)
                     k = i12(j,m)
                     ka = atomic(k)
                     if (ka .eq. 8)  o2i(i) = k
                     if (ka.eq.6 .and. k.ne.c1i(i))  c3i(i) = k
                  end do
                  m = c3i(i)
                  do j = 1, n12(m)
                     k = i12(j,m)
                     ka = atomic(k)
                     if (ka .eq. 8)  o3i(i) = k
                  end do
                  m = c4i(i)
                  do j = 1, n12(m)
                     k = i12(j,m)
                     ka = atomic(k)
                     if (ka.eq.6 .and. k.ne.c3i(i))  c5i(i) = k
                  end do
                  m = c5i(i)
                  do j = 1, n12(m)
                     k = i12(j,m)
                     ka = atomic(k)
                     if (ka .eq. 8)  o5i(i) = k
                  end do
                  m = o5i(i)
                  do j = 1, n12(m)
                     k = i12(j,m)
                     ka = atomic(k)
                     if (ka .eq. 15)  poi(i) = k
                  end do
               end if
               if (i .gt. 1) then
                  resname = nuclz(seqtyp(i-1))
                  if (resname .eq. 'MP ')  poi(i) = 0
                  if (resname .eq. 'DP ')  poi(i) = 0
                  if (resname .eq. 'TP ')  poi(i) = 0
               end if
               m = poi(i)
               if (m .ne. 0) then
                  do j = 1, n12(m)
                     k = i12(j,m)
                     ka = atomic(k)
                     if (ka.eq.8 .and. n12(k).eq.1) then
                        if (op1(i) .eq. 0) then
                           op1(i) = k
                        else if (op2(i) .eq. 0) then
                           op2(i) = k
                        else
                           op3(i) = k
                        end if
                     end if
                  end do
               end if
            end do
         end if
      end do
c
c     copy the atoms of each biopolymer residue into PDB format
c
      do m = 1, nchain
         start = ichain(1,m)
         stop = ichain(2,m)
         if (chntyp(m) .eq. 'PEPTIDE') then
            do i = start, stop
               resname = amino(seqtyp(i))
               if (resname .eq. 'H2N') then
                  continue
               else if (resname .eq. 'FOR') then
                  call pdbatom (' C  ',resname,i,ci(i))
                  call pdbatom (' O  ',resname,i,oi(i))
               else if (resname .eq. 'ACE') then
                  call pdbatom (' CH3',resname,i,cai(i))
                  call pdbatom (' C  ',resname,i,ci(i))
                  call pdbatom (' O  ',resname,i,oi(i))
               else if (resname .eq. 'COH') then
                  call pdbatom (' OH ',resname,i,ni(i))
               else if (resname .eq. 'NH2') then
                  call pdbatom (' N  ',resname,i,ni(i))
               else if (resname .eq. 'NME') then
                  call pdbatom (' N  ',resname,i,ni(i))
                  call pdbatom (' CH3',resname,i,cai(i))
               else
                  call pdbatom (' N  ',resname,i,ni(i))
                  call pdbatom (' CA ',resname,i,cai(i))
                  call pdbatom (' C  ',resname,i,ci(i))
                  call pdbatom (' O  ',resname,i,oi(i))
               end if
               call getside (resname,i,ci(i),cai(i),cbi)
               if ((resname.eq.'CYS'.or.resname.eq.'CYX')
     &               .and. cbi.ne.0) then
                  resname = 'CYS'
                  do j = 1, n13(cbi)
                     if (atomic(i13(j,cbi)) .eq. 16)  resname = 'CYX'
                  end do
               end if
               if (i.eq.stop .and. ci(i).ne.0) then
                  do j = 1, n12(ci(i))
                     k = i12(j,ci(i))
                     if (atomic(k).eq.8 .and. k.ne.oi(i)) then
                        call pdbatom (' OXT',resname,i,k)
                        goto 20
                     end if
                  end do
   20             continue
               end if
               call getproh (resname,i,m,ni(i),cai(i),cbi)
            end do
         else if (chntyp(m) .eq. 'NUCLEIC') then
            do i = start, stop
               resname = nuclz(seqtyp(i))
               if (resname .eq. 'MP ') then
                  call pdbatom (' P  ',resname,i,poi(i))
                  call pdbatom (' OP1',resname,i,op1(i))
                  call pdbatom (' OP2',resname,i,op2(i))
                  call pdbatom (' OP3',resname,i,op3(i))
               else if (resname .eq. 'DP ') then
               else if (resname .eq. 'TP ') then
               else
                  call pdbatom (' P  ',resname,i,poi(i))
                  call pdbatom (' OP1',resname,i,op1(i))
                  call pdbatom (' OP2',resname,i,op2(i))
                  call pdbatom (' O5''',resname,i,o5i(i))
                  call pdbatom (' C5''',resname,i,c5i(i))
                  call pdbatom (' C4''',resname,i,c4i(i))
                  call pdbatom (' O4''',resname,i,o4i(i))
                  call pdbatom (' C3''',resname,i,c3i(i))
                  call pdbatom (' O3''',resname,i,o3i(i))
                  call pdbatom (' C2''',resname,i,c2i(i))
                  call pdbatom (' O2''',resname,i,o2i(i))
                  call pdbatom (' C1''',resname,i,c1i(i))
                  call getbase (resname,i,ni(i))
                  call getnuch (resname,i,ni(i),c1i(i),c2i(i),o2i(i),
     &                          c3i(i),o3i(i),c4i(i),c5i(i),o5i(i))
               end if
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (ni)
      deallocate (cai)
      deallocate (ci)
      deallocate (oi)
      deallocate (poi)
      deallocate (op1)
      deallocate (op2)
      deallocate (op3)
      deallocate (c5i)
      deallocate (o5i)
      deallocate (c4i)
      deallocate (o4i)
      deallocate (c3i)
      deallocate (o3i)
      deallocate (c2i)
      deallocate (o2i)
      deallocate (c1i)
      deallocate (restyp)
c
c     perform dynamic allocation of some local arrays
c
      allocate (hetmol(nmol))
c
c     copy any water, ions or ligands following biopolymer chains
c
      if (.not. generic) then
         do i = 1, nmol
            hetmol(i) = .true.
         end do
         do i = 1, n
            if (pdblist(i) .ne. 0)  hetmol(molcule(i)) = .false.
         end do
         do i = 1, nmol
            if (hetmol(i)) then
               do j = imol(1,i), imol(2,i)
                  k = kmol(j)
                  atmnum = atomic(k)
                  atmname = ' '//name(k)
                  justify = 0
                  call numeral (type(k),resname,justify)
                  if (water(i)) then
                     if (atmnum .eq. 1)  atmname = ' H  '
                     if (atmnum .eq. 8)  atmname = ' O  '
                     resname = 'HOH'
                  else if (atmnum .eq. 11) then
                     atmname = 'NA  '
                     resname = ' NA'
                  else if (atmnum .eq. 12) then
                     atmname = 'MG  '
                     resname = ' MG'
                  else if (atmnum .eq. 17) then
                     atmname = 'CL  '
                     resname = ' CL'
                  else if (atmnum .eq. 19) then
                     atmname = ' K  '
                     resname = '  K'
                  else if (atmnum .eq. 20) then
                     atmname = 'CA  '
                     resname = ' CA'
                  else if (atmnum .eq. 35) then
                     atmname = 'BR  '
                     resname = ' BR'
                  else if (atmnum .eq. 53) then
                     atmname = ' I  '
                     resname = '  I'
                  end if
                  pdbnum = nseq + i - 1
                  call pdbatom (atmname,resname,pdbnum,k)
                  pdbtyp(npdb) = 'HETATM'
               end do
            end if
         end do
         do i = 1, nmol
            if (hetmol(i)) then
               do j = imol(1,i), imol(2,i)
                  k = kmol(j)
                  kp = pdblist(k)
                  npdb12(kp) = n12(k)
                  do m = 1, n12(k)
                     ipdb12(m,kp) = pdblist(i12(m,k))
                  end do
               end do
            end if
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (water)
      deallocate (hetmol)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine pdbatom  --  add a single atom to PDB file  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "pdbatom" adds an atom to the Protein Data Bank file
c
c
      subroutine pdbatom (atmname,resname,ires,icoord)
      use sizes
      use atoms
      use pdb
      implicit none
      integer ires,icoord
      character*3 resname
      character*4 atmname
c
c
c     for each atom set the sequential number, record type, atom
c     name, residue name, residue number and atomic coordinates
c
      if (icoord .ne. 0) then
         npdb = npdb + 1
         pdbtyp(npdb) = 'ATOM  '
         pdbatm(npdb) = atmname
         pdbres(npdb) = resname
         resnum(npdb) = ires
         xpdb(npdb) = x(icoord)
         ypdb(npdb) = y(icoord)
         zpdb(npdb) = z(icoord)
         npdb12(npdb) = 0
         pdblist(icoord) = npdb
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine getside  --  extract the amino acid side chains  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "getside" finds the side chain heavy atoms for a single amino
c     acid residue and copies the names and coordinates to the Protein
c     Data Bank file
c
c
      subroutine getside (resname,ires,ci,cai,cbi)
      use sizes
      use atomid
      use atoms
      use couple
      implicit none
      integer i,j,ires
      integer ci,cai,cbi
      character*3 resname
c
c
c     if residue is a terminal cap, there is no side chain
c
      cbi = 0
      if (resname .eq. 'H2N')  return
      if (resname .eq. 'FOR')  return
      if (resname .eq. 'ACE')  return
      if (resname .eq. 'COH')  return
      if (resname .eq. 'NH2')  return
      if (resname .eq. 'NME')  return
c
c     find the beta carbon atom for the current residue
c
      do i = 1, n
         if (i.ne.ci .and. atomic(i).eq.6) then
            do j = 1, 4
               if (i12(j,i) .eq. cai) then
                  cbi = i
                  if (resname .ne. 'AIB') then
                     call pdbatom (' CB ',resname,ires,cbi)
                  else
                     call pdbatom (' CB1',resname,ires,cbi)
                  end if
                  goto 10
               end if
            end do
         end if
      end do
   10 continue
      if (cbi .eq. 0)  return
c
c     glycine residue  (GLY)
c
      if (resname .eq. 'GLY') then
         continue
c
c     alanine residue  (ALA)
c
      else if (resname .eq. 'ALA') then
         continue
c
c     valine residue  (VAL)
c
      else if (resname .eq. 'VAL') then
         call pdbatom (' CG1',resname,ires,cbi+1)
         call pdbatom (' CG2',resname,ires,cbi+2)
c
c     leucine residue  (LEU)
c
      else if (resname .eq. 'LEU') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD1',resname,ires,cbi+2)
         call pdbatom (' CD2',resname,ires,cbi+3)
c
c     isoleucine residue  (ILE)
c
      else if (resname .eq. 'ILE') then
         call pdbatom (' CG1',resname,ires,cbi+1)
         call pdbatom (' CG2',resname,ires,cbi+2)
         call pdbatom (' CD1',resname,ires,cbi+3)
c
c     serine residue  (SER)
c
      else if (resname .eq. 'SER') then
         call pdbatom (' OG ',resname,ires,cbi+1)
c
c     threonine residue  (THR)
c
      else if (resname .eq. 'THR') then
         call pdbatom (' OG1',resname,ires,cbi+1)
         call pdbatom (' CG2',resname,ires,cbi+2)
c
c     cysteine residue  (CYS)
c
      else if (resname .eq. 'CYS') then
         call pdbatom (' SG ',resname,ires,cbi+1)
c
c     cysteine residue  (CYX)
c
      else if (resname .eq. 'CYX') then
         call pdbatom (' SG ',resname,ires,cbi+1)
c
c     deprotonated cysteine residue  (CYD)
c
      else if (resname .eq. 'CYD') then
         call pdbatom (' SG ',resname,ires,cbi+1)
c
c     proline residue  (PRO)
c
      else if (resname .eq. 'PRO') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD ',resname,ires,cbi+2)
c
c     phenylalanine residue  (PHE)
c
      else if (resname .eq. 'PHE') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD1',resname,ires,cbi+2)
         call pdbatom (' CD2',resname,ires,cbi+3)
         call pdbatom (' CE1',resname,ires,cbi+4)
         call pdbatom (' CE2',resname,ires,cbi+5)
         call pdbatom (' CZ ',resname,ires,cbi+6)
c
c     tyrosine residue  (TYR)
c
      else if (resname .eq. 'TYR') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD1',resname,ires,cbi+2)
         call pdbatom (' CD2',resname,ires,cbi+3)
         call pdbatom (' CE1',resname,ires,cbi+4)
         call pdbatom (' CE2',resname,ires,cbi+5)
         call pdbatom (' CZ ',resname,ires,cbi+6)
         call pdbatom (' OH ',resname,ires,cbi+7)
c
c     deprotonated tyrosine residue  (TYD)
c
      else if (resname .eq. 'TYD') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD1',resname,ires,cbi+2)
         call pdbatom (' CD2',resname,ires,cbi+3)
         call pdbatom (' CE1',resname,ires,cbi+4)
         call pdbatom (' CE2',resname,ires,cbi+5)
         call pdbatom (' CZ ',resname,ires,cbi+6)
         call pdbatom (' OH ',resname,ires,cbi+7)
c
c     tryptophan residue  (TRP)
c
      else if (resname .eq. 'TRP') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD1',resname,ires,cbi+2)
         call pdbatom (' CD2',resname,ires,cbi+3)
         call pdbatom (' NE1',resname,ires,cbi+4)
         call pdbatom (' CE2',resname,ires,cbi+5)
         call pdbatom (' CE3',resname,ires,cbi+6)
         call pdbatom (' CZ2',resname,ires,cbi+7)
         call pdbatom (' CZ3',resname,ires,cbi+8)
         call pdbatom (' CH2',resname,ires,cbi+9)
c
c     histidine (HD and HE) residue  (HIS)
c
      else if (resname .eq. 'HIS') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' ND1',resname,ires,cbi+2)
         call pdbatom (' CD2',resname,ires,cbi+3)
         call pdbatom (' CE1',resname,ires,cbi+4)
         call pdbatom (' NE2',resname,ires,cbi+5)
c
c     histidine (HD only) residue  (HID)
c
      else if (resname .eq. 'HID') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' ND1',resname,ires,cbi+2)
         call pdbatom (' CD2',resname,ires,cbi+3)
         call pdbatom (' CE1',resname,ires,cbi+4)
         call pdbatom (' NE2',resname,ires,cbi+5)
c
c     histidine (HE only) residue  (HIE)
c
      else if (resname .eq. 'HIE') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' ND1',resname,ires,cbi+2)
         call pdbatom (' CD2',resname,ires,cbi+3)
         call pdbatom (' CE1',resname,ires,cbi+4)
         call pdbatom (' NE2',resname,ires,cbi+5)
c
c     aspartic acid residue  (ASP)
c
      else if (resname .eq. 'ASP') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' OD1',resname,ires,cbi+2)
         call pdbatom (' OD2',resname,ires,cbi+3)
c
c     protonated aspartic acid residue  (ASH)
c
      else if (resname .eq. 'ASH') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' OD1',resname,ires,cbi+2)
         call pdbatom (' OD2',resname,ires,cbi+3)
c
c     asparagine residue  (ASN)
c
      else if (resname .eq. 'ASN') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' OD1',resname,ires,cbi+2)
         call pdbatom (' ND2',resname,ires,cbi+3)
c
c     glutamic acid residue  (GLU)
c
      else if (resname .eq. 'GLU') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD ',resname,ires,cbi+2)
         call pdbatom (' OE1',resname,ires,cbi+3)
         call pdbatom (' OE2',resname,ires,cbi+4)
c
c     protonated glutamic acid residue  (GLH)
c
      else if (resname .eq. 'GLH') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD ',resname,ires,cbi+2)
         call pdbatom (' OE1',resname,ires,cbi+3)
         call pdbatom (' OE2',resname,ires,cbi+4)
c
c     glutamine residue  (GLN)
c
      else if (resname .eq. 'GLN') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD ',resname,ires,cbi+2)
         call pdbatom (' OE1',resname,ires,cbi+3)
         call pdbatom (' NE2',resname,ires,cbi+4)
c
c     methionine residue  (MET)
c
      else if (resname .eq. 'MET') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' SD ',resname,ires,cbi+2)
         call pdbatom (' CE ',resname,ires,cbi+3)
c
c     lysine residue  (LYS)
c
      else if (resname .eq. 'LYS') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD ',resname,ires,cbi+2)
         call pdbatom (' CE ',resname,ires,cbi+3)
         call pdbatom (' NZ ',resname,ires,cbi+4)
c
c     deprotonated lysine residue  (LYD)
c
      else if (resname .eq. 'LYD') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD ',resname,ires,cbi+2)
         call pdbatom (' CE ',resname,ires,cbi+3)
         call pdbatom (' NZ ',resname,ires,cbi+4)
c
c     arginine residue  (ARG)
c
      else if (resname .eq. 'ARG') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD ',resname,ires,cbi+2)
         call pdbatom (' NE ',resname,ires,cbi+3)
         call pdbatom (' CZ ',resname,ires,cbi+4)
         call pdbatom (' NH1',resname,ires,cbi+5)
         call pdbatom (' NH2',resname,ires,cbi+6)
c
c     ornithine residue  (ORN)
c
      else if (resname .eq. 'ORN') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD ',resname,ires,cbi+2)
         call pdbatom (' NE ',resname,ires,cbi+3)
c
c     methylalanine residue  (AIB)
c
      else if (resname .eq. 'AIB') then
         call pdbatom (' CB2',resname,ires,cbi+1)
c
c     pyroglutamic acid residue  (PCA)
c
      else if (resname .eq. 'PCA') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD ',resname,ires,cbi+2)
         call pdbatom (' OE ',resname,ires,cbi+3)
c
c     unknown residue  (UNK)
c
      else if (resname .eq. 'UNK') then
         continue
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine getproh  --  extract the amino acid hydrogens  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "getproh" finds the hydrogen atoms for a single amino acid
c     residue and copies the names and coordinates to the Protein
c     Data Bank file
c
c
      subroutine getproh (resname,ires,jchain,ni,cai,cbi)
      use sizes
      use atomid
      use atoms
      use couple
      use fields
      use sequen
      implicit none
      integer i,nh,hca
      integer ires,jchain
      integer ni,cai,cbi
      logical allatom
      character*3 resname
      character*4 atmname
c
c
c     get any amide hydrogen atoms for non-N-terminal residues
c
      if (ires.ne.ichain(1,jchain) .or. n12(ni).ne.4) then
         if (resname .ne. 'PRO') then
            do i = 1, n
               if (atomic(i).eq.1 .and. i12(1,i).eq.ni) then
                  if (resname .eq. 'COH') then
                     call pdbatom (' HO ',resname,ires,i)
                  else if (resname .eq. 'NH2') then
                     call pdbatom (' H1 ',resname,ires,i)
                     call pdbatom (' H2 ',resname,ires,i+1)
                  else
                     call pdbatom (' H  ',resname,ires,i)
                  end if
                  goto 10
               end if
            end do
         end if
c
c     get any amide hydrogen atoms for N-terminal residues
c
      else
         if (resname .eq. 'PRO') then
            nh = 0
            do i = 1, n
               if (atomic(i).eq.1 .and. i12(1,i).eq.ni) then
                  nh = nh + 1
                  if (nh .eq. 1) then
                     atmname = ' H1 '
                  else if (nh .eq. 2) then
                     atmname = ' H2 '
                  end if
                  call pdbatom (atmname,resname,ires,i)
                  if (nh .eq. 2)  goto 10
               end if
            end do
         else if (resname .eq. 'PCA') then
            do i = 1, n
               if (atomic(i).eq.1 .and. i12(1,i).eq.ni) then
                  atmname = ' H  '
                  call pdbatom (atmname,resname,ires,i)
                  goto 10
               end if
            end do
         else
            nh = 0
            do i = 1, n
               if (atomic(i).eq.1 .and. i12(1,i).eq.ni) then
                  nh = nh + 1
                  if (nh .eq. 1) then
                     atmname = ' H1 '
                  else if (nh .eq. 2) then
                     atmname = ' H2 '
                  else if (nh .eq. 3) then
                     atmname = ' H3 '
                  end if
                  call pdbatom (atmname,resname,ires,i)
                  if (nh .eq. 3)  goto 10
               end if
            end do
         end if
      end if
   10 continue
c
c     get the alpha hydrogen atom for the current residue
c
      hca = 0
      do i = 1, n
         if (atomic(i).eq.1 .and. i12(1,i).eq.cai) then
            hca = i
            if (resname .eq. 'GLY') then
               atmname = ' HA2'
            else if (resname .eq. 'FOR') then
               atmname = ' H  '
            else if (resname .eq. 'ACE') then
               atmname = ' H1 '
            else if (resname .eq. 'NME') then
               atmname = ' H1 '
            else
               atmname = ' HA '
            end if
            call pdbatom (atmname,resname,ires,i)
            goto 20
         end if
      end do
   20 continue
c
c     backbone only if no alpha hydrogen or beta carbon
c
      if (hca.eq.0 .and. cbi.eq.0)  return
c
c     if no alpha hydrogen, then united atom force field
c
      if (hca .ne. 0) then
         allatom = .true.
      else if (resname .eq. 'AIB') then
         if (n12(cbi) .eq. 1) then
            allatom = .false.
         else
            allatom = .true.
         end if
      else
         allatom = .false.
      end if
c
c     glycine residue  (GLY)
c
      if (resname .eq. 'GLY') then
         if (allatom) then
            call pdbatom (' HA3',resname,ires,hca+1)
         end if
c
c     alanine residue  (ALA)
c
      else if (resname .eq. 'ALA') then
         if (allatom) then
            call pdbatom (' HB1',resname,ires,hca+2)
            call pdbatom (' HB2',resname,ires,hca+3)
            call pdbatom (' HB3',resname,ires,hca+4)
         end if
c
c     valine residue  (VAL)
c
      else if (resname .eq. 'VAL') then
         if (allatom) then
            call pdbatom (' HB ',resname,ires,hca+4)
            call pdbatom ('HG11',resname,ires,hca+5)
            call pdbatom ('HG12',resname,ires,hca+6)
            call pdbatom ('HG13',resname,ires,hca+7)
            call pdbatom ('HG21',resname,ires,hca+8)
            call pdbatom ('HG22',resname,ires,hca+9)
            call pdbatom ('HG23',resname,ires,hca+10)
         end if
c
c     leucine residue  (LEU)
c
      else if (resname .eq. 'LEU') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+5)
            call pdbatom (' HB3',resname,ires,hca+6)
            call pdbatom (' HG ',resname,ires,hca+7)
            call pdbatom ('HD11',resname,ires,hca+8)
            call pdbatom ('HD12',resname,ires,hca+9)
            call pdbatom ('HD13',resname,ires,hca+10)
            call pdbatom ('HD21',resname,ires,hca+11)
            call pdbatom ('HD22',resname,ires,hca+12)
            call pdbatom ('HD23',resname,ires,hca+13)
         end if
c
c     isoleucine residue  (ILE)
c
      else if (resname .eq. 'ILE') then
         if (allatom) then
            call pdbatom (' HB ',resname,ires,hca+5)
            call pdbatom ('HG12',resname,ires,hca+6)
            call pdbatom ('HG13',resname,ires,hca+7)
            call pdbatom ('HG21',resname,ires,hca+8)
            call pdbatom ('HG22',resname,ires,hca+9)
            call pdbatom ('HG23',resname,ires,hca+10)
            call pdbatom ('HD11',resname,ires,hca+11)
            call pdbatom ('HD12',resname,ires,hca+12)
            call pdbatom ('HD13',resname,ires,hca+13)
         end if
c
c     serine residue  (SER)
c
      else if (resname .eq. 'SER') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+3)
            call pdbatom (' HB3',resname,ires,hca+4)
            call pdbatom (' HG ',resname,ires,hca+5)
         else
            call pdbatom (' HG ',resname,ires,cbi+2)
         end if
c
c     threonine residue  (THR)
c
      else if (resname .eq. 'THR') then
         if (allatom) then
            call pdbatom (' HB ',resname,ires,hca+4)
            call pdbatom (' HG1',resname,ires,hca+5)
            call pdbatom ('HG21',resname,ires,hca+6)
            call pdbatom ('HG22',resname,ires,hca+7)
            call pdbatom ('HG23',resname,ires,hca+8)
         else
            call pdbatom (' HG1',resname,ires,cbi+3)
         end if
c
c     cysteine residue  (CYS)
c
      else if (resname .eq. 'CYS') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+3)
            call pdbatom (' HB3',resname,ires,hca+4)
            call pdbatom (' HG ',resname,ires,hca+5)
         else if (biotyp(86) .ne. 0) then
            call pdbatom (' HG ',resname,ires,cbi+2)
         end if
c
c     cystine residue  (CYX)
c
      else if (resname .eq. 'CYX') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+3)
            call pdbatom (' HB3',resname,ires,hca+4)
         end if
c
c     deprotonated cysteine residue  (CYD)
c
      else if (resname .eq. 'CYD') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+3)
            call pdbatom (' HB3',resname,ires,hca+4)
         end if
c
c     proline residue  (PRO)
c
      else if (resname .eq. 'PRO') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+4)
            call pdbatom (' HB3',resname,ires,hca+5)
            call pdbatom (' HG2',resname,ires,hca+6)
            call pdbatom (' HG3',resname,ires,hca+7)
            call pdbatom (' HD2',resname,ires,hca+8)
            call pdbatom (' HD3',resname,ires,hca+9)
         end if
c
c     phenylalanine residue  (PHE)
c
      else if (resname .eq. 'PHE') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+8)
            call pdbatom (' HB3',resname,ires,hca+9)
            call pdbatom (' HD1',resname,ires,hca+10)
            call pdbatom (' HD2',resname,ires,hca+11)
            call pdbatom (' HE1',resname,ires,hca+12)
            call pdbatom (' HE2',resname,ires,hca+13)
            call pdbatom (' HZ ',resname,ires,hca+14)
         else if (biotyp(126) .ne. 0) then
            call pdbatom (' HD1',resname,ires,cbi+7)
            call pdbatom (' HD2',resname,ires,cbi+8)
            call pdbatom (' HE1',resname,ires,cbi+9)
            call pdbatom (' HE2',resname,ires,cbi+10)
            call pdbatom (' HZ ',resname,ires,cbi+11)
         end if
c
c     tyrosine residue  (TYR)
c
      else if (resname .eq. 'TYR') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+9)
            call pdbatom (' HB3',resname,ires,hca+10)
            call pdbatom (' HD1',resname,ires,hca+11)
            call pdbatom (' HD2',resname,ires,hca+12)
            call pdbatom (' HE1',resname,ires,hca+13)
            call pdbatom (' HE2',resname,ires,hca+14)
            call pdbatom (' HH ',resname,ires,hca+15)
         else if (biotyp(141) .ne. 0) then
            call pdbatom (' HD1',resname,ires,cbi+8)
            call pdbatom (' HD2',resname,ires,cbi+9)
            call pdbatom (' HE1',resname,ires,cbi+10)
            call pdbatom (' HE2',resname,ires,cbi+11)
            call pdbatom (' HH ',resname,ires,cbi+12)
         else
            call pdbatom (' HH ',resname,ires,cbi+8)
         end if
c
c     deprotonated tyrosine residue  (TYD)
c
      else if (resname .eq. 'TYD') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+9)
            call pdbatom (' HB3',resname,ires,hca+10)
            call pdbatom (' HD1',resname,ires,hca+11)
            call pdbatom (' HD2',resname,ires,hca+12)
            call pdbatom (' HE1',resname,ires,hca+13)
            call pdbatom (' HE2',resname,ires,hca+14)
         else if (biotyp(141) .ne. 0) then
            call pdbatom (' HD1',resname,ires,cbi+8)
            call pdbatom (' HD2',resname,ires,cbi+9)
            call pdbatom (' HE1',resname,ires,cbi+10)
            call pdbatom (' HE2',resname,ires,cbi+11)
         end if
c
c     tryptophan residue  (TRP)
c
      else if (resname .eq. 'TRP') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+11)
            call pdbatom (' HB3',resname,ires,hca+12)
            call pdbatom (' HD1',resname,ires,hca+13)
            call pdbatom (' HE1',resname,ires,hca+14)
            call pdbatom (' HE3',resname,ires,hca+15)
            call pdbatom (' HZ2',resname,ires,hca+16)
            call pdbatom (' HZ3',resname,ires,hca+17)
            call pdbatom (' HH2',resname,ires,hca+18)
         else if (biotyp(172) .ne. 0) then
            call pdbatom (' HD1',resname,ires,cbi+10)
            call pdbatom (' HE1',resname,ires,cbi+11)
            call pdbatom (' HE3',resname,ires,cbi+12)
            call pdbatom (' HZ2',resname,ires,cbi+13)
            call pdbatom (' HZ3',resname,ires,cbi+14)
            call pdbatom (' HH2',resname,ires,cbi+15)
         else
            call pdbatom (' HE1',resname,ires,cbi+10)
         end if
c
c     histidine (HD and HE) residue  (HIS)
c
      else if (resname .eq. 'HIS') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+7)
            call pdbatom (' HB3',resname,ires,hca+8)
            call pdbatom (' HD1',resname,ires,hca+9)
            call pdbatom (' HD2',resname,ires,hca+10)
            call pdbatom (' HE1',resname,ires,hca+11)
            call pdbatom (' HE2',resname,ires,hca+12)
         else if (biotyp(197) .ne. 0) then
            call pdbatom (' HD1',resname,ires,cbi+6)
            call pdbatom (' HD2',resname,ires,cbi+7)
            call pdbatom (' HE1',resname,ires,cbi+8)
            call pdbatom (' HE2',resname,ires,cbi+9)
         else
            call pdbatom (' HD1',resname,ires,cbi+6)
            call pdbatom (' HE2',resname,ires,cbi+7)
         end if
c
c     histidine (HD only) residue  (HID)
c
      else if (resname .eq. 'HID') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+7)
            call pdbatom (' HB3',resname,ires,hca+8)
            call pdbatom (' HD1',resname,ires,hca+9)
            call pdbatom (' HD2',resname,ires,hca+10)
            call pdbatom (' HE1',resname,ires,hca+11)
         else if (biotyp(214) .ne. 0) then
            call pdbatom (' HD1',resname,ires,cbi+6)
            call pdbatom (' HD2',resname,ires,cbi+7)
            call pdbatom (' HE1',resname,ires,cbi+8)
         else
            call pdbatom (' HD1',resname,ires,cbi+6)
         end if
c
c     histidine (HE only) residue  (HIE)
c
      else if (resname .eq. 'HIE') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+7)
            call pdbatom (' HB3',resname,ires,hca+8)
            call pdbatom (' HD2',resname,ires,hca+9)
            call pdbatom (' HE1',resname,ires,hca+10)
            call pdbatom (' HE2',resname,ires,hca+11)
         else if (biotyp(229) .ne. 0) then
            call pdbatom (' HD2',resname,ires,cbi+6)
            call pdbatom (' HE1',resname,ires,cbi+7)
            call pdbatom (' HE2',resname,ires,cbi+8)
         else
            call pdbatom (' HE2',resname,ires,cbi+6)
         end if
c
c     aspartic acid residue  (ASP)
c
      else if (resname .eq. 'ASP') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+5)
            call pdbatom (' HB3',resname,ires,hca+6)
         end if
c
c     protonated aspartic acid residue  (ASH)
c
      else if (resname .eq. 'ASH') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+5)
            call pdbatom (' HB3',resname,ires,hca+6)
            call pdbatom (' HD2',resname,ires,hca+7)
         else
            call pdbatom (' HD2',resname,ires,cbi+4)
         end if
c
c     asparagine residue  (ASN)
c
      else if (resname .eq. 'ASN') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+5)
            call pdbatom (' HB3',resname,ires,hca+6)
            call pdbatom ('HD21',resname,ires,hca+7)
            call pdbatom ('HD22',resname,ires,hca+8)
         else
            call pdbatom ('HD21',resname,ires,cbi+4)
            call pdbatom ('HD22',resname,ires,cbi+5)
         end if
c
c     glutamic acid residue  (GLU)
c
      else if (resname .eq. 'GLU') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+6)
            call pdbatom (' HB3',resname,ires,hca+7)
            call pdbatom (' HG2',resname,ires,hca+8)
            call pdbatom (' HG3',resname,ires,hca+9)
         end if
c
c     protonated glutamic acid residue  (GLH)
c
      else if (resname .eq. 'GLH') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+6)
            call pdbatom (' HB3',resname,ires,hca+7)
            call pdbatom (' HG2',resname,ires,hca+8)
            call pdbatom (' HG3',resname,ires,hca+9)
            call pdbatom (' HE2',resname,ires,hca+10)
         else
            call pdbatom (' HE2',resname,ires,cbi+5)
         end if
c
c     glutamine residue  (GLN)
c
      else if (resname .eq. 'GLN') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+6)
            call pdbatom (' HB3',resname,ires,hca+7)
            call pdbatom (' HG2',resname,ires,hca+8)
            call pdbatom (' HG3',resname,ires,hca+9)
            call pdbatom ('HE21',resname,ires,hca+10)
            call pdbatom ('HE22',resname,ires,hca+11)
         else
            call pdbatom ('HE21',resname,ires,cbi+5)
            call pdbatom ('HE22',resname,ires,cbi+6)
         end if
c
c     methionine residue  (MET)
c
      else if (resname .eq. 'MET') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+5)
            call pdbatom (' HB3',resname,ires,hca+6)
            call pdbatom (' HG2',resname,ires,hca+7)
            call pdbatom (' HG3',resname,ires,hca+8)
            call pdbatom (' HE1',resname,ires,hca+9)
            call pdbatom (' HE2',resname,ires,hca+10)
            call pdbatom (' HE3',resname,ires,hca+11)
         end if
c
c     lysine residue  (LYS)
c
      else if (resname .eq. 'LYS') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+6)
            call pdbatom (' HB3',resname,ires,hca+7)
            call pdbatom (' HG2',resname,ires,hca+8)
            call pdbatom (' HG3',resname,ires,hca+9)
            call pdbatom (' HD2',resname,ires,hca+10)
            call pdbatom (' HD3',resname,ires,hca+11)
            call pdbatom (' HE2',resname,ires,hca+12)
            call pdbatom (' HE3',resname,ires,hca+13)
            call pdbatom (' HZ1',resname,ires,hca+14)
            call pdbatom (' HZ2',resname,ires,hca+15)
            call pdbatom (' HZ3',resname,ires,hca+16)
         else
            call pdbatom (' HZ1',resname,ires,cbi+5)
            call pdbatom (' HZ2',resname,ires,cbi+6)
            call pdbatom (' HZ3',resname,ires,cbi+7)
         end if
c
c     deprotonated lysine residue  (LYD)
c
      else if (resname .eq. 'LYD') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+6)
            call pdbatom (' HB3',resname,ires,hca+7)
            call pdbatom (' HG2',resname,ires,hca+8)
            call pdbatom (' HG3',resname,ires,hca+9)
            call pdbatom (' HD2',resname,ires,hca+10)
            call pdbatom (' HD3',resname,ires,hca+11)
            call pdbatom (' HE2',resname,ires,hca+12)
            call pdbatom (' HE3',resname,ires,hca+13)
            call pdbatom (' HZ1',resname,ires,hca+14)
            call pdbatom (' HZ2',resname,ires,hca+15)
         else
            call pdbatom (' HZ1',resname,ires,cbi+5)
            call pdbatom (' HZ2',resname,ires,cbi+6)
         end if
c
c     arginine residue  (ARG)
c
      else if (resname .eq. 'ARG') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+8)
            call pdbatom (' HB3',resname,ires,hca+9)
            call pdbatom (' HG2',resname,ires,hca+10)
            call pdbatom (' HG3',resname,ires,hca+11)
            call pdbatom (' HD2',resname,ires,hca+12)
            call pdbatom (' HD3',resname,ires,hca+13)
            call pdbatom (' HE ',resname,ires,hca+14)
            call pdbatom ('HH11',resname,ires,hca+15)
            call pdbatom ('HH12',resname,ires,hca+16)
            call pdbatom ('HH21',resname,ires,hca+17)
            call pdbatom ('HH22',resname,ires,hca+18)
         else
            call pdbatom (' HE ',resname,ires,cbi+7)
            call pdbatom ('HH11',resname,ires,cbi+8)
            call pdbatom ('HH12',resname,ires,cbi+9)
            call pdbatom ('HH21',resname,ires,cbi+10)
            call pdbatom ('HH22',resname,ires,cbi+11)
         end if
c
c     ornithine residue  (ORN)
c
      else if (resname .eq. 'ORN') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+5)
            call pdbatom (' HB3',resname,ires,hca+6)
            call pdbatom (' HG2',resname,ires,hca+7)
            call pdbatom (' HG3',resname,ires,hca+8)
            call pdbatom (' HD2',resname,ires,hca+9)
            call pdbatom (' HD3',resname,ires,hca+10)
            call pdbatom (' HE1',resname,ires,hca+11)
            call pdbatom (' HE2',resname,ires,hca+12)
            call pdbatom (' HE3',resname,ires,hca+13)
         else
            call pdbatom (' HE1',resname,ires,cbi+4)
            call pdbatom (' HE2',resname,ires,cbi+5)
            call pdbatom (' HE3',resname,ires,cbi+6)
         end if
c
c     methylalanine residue  (AIB)
c
      else if (resname .eq. 'AIB') then
         if (allatom) then
            call pdbatom ('HB11',resname,ires,cbi+2)
            call pdbatom ('HB12',resname,ires,cbi+3)
            call pdbatom ('HB13',resname,ires,cbi+4)
            call pdbatom ('HB21',resname,ires,cbi+5)
            call pdbatom ('HB22',resname,ires,cbi+6)
            call pdbatom ('HB23',resname,ires,cbi+7)
         end if
c
c     pyroglutamic acid residue  (PCA)
c
      else if (resname .eq. 'PCA') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+5)
            call pdbatom (' HB3',resname,ires,hca+6)
            call pdbatom (' HG2',resname,ires,hca+7)
            call pdbatom (' HG3',resname,ires,hca+8)
         end if
c
c     unknown residue  (UNK)
c
      else if (resname .eq. 'UNK') then
         if (allatom) then
            call pdbatom (' HA3',resname,ires,hca+1)
         end if
c
c     N-terminal deprotonated residue  (H2N)
c
      else if (resname .eq. 'H2N') then
         continue
c
c     N-terminal formyl residue  (FOR)
c
      else if (resname .eq. 'FOR') then
         continue
c
c     N-terminal acetyl residue  (ACE)
c
      else if (resname .eq. 'ACE') then
         if (allatom) then
            call pdbatom (' H2 ',resname,ires,hca+1)
            call pdbatom (' H3 ',resname,ires,hca+2)
         end if
c
c     C-terminal protonated residue (COH)
c
      else if (resname .eq. 'COH') then
         continue
c
c     C-terminal amide residue  (NH2)
c
      else if (resname .eq. 'NH2') then
         continue
c
c     C-terminal N-methylamide residue  (NME)
c
      else if (resname .eq. 'NME') then
         if (allatom) then
            call pdbatom (' H2 ',resname,ires,hca+1)
            call pdbatom (' H3 ',resname,ires,hca+2)
         end if
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine getbase  --  extract the nucleotide side chains  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "getbase" finds the base heavy atoms for a single nucleotide
c     residue and copies the names and coordinates to the Protein
c     Data Bank file
c
c
      subroutine getbase (resname,ires,ni)
      implicit none
      integer ires,ni
      character*3 resname
c
c
c     adenine in adenosine residue  (A)
c
      if (resname .eq. 'A  ') then
         call pdbatom (' N9 ',resname,ires,ni)
         call pdbatom (' C8 ',resname,ires,ni+1)
         call pdbatom (' N7 ',resname,ires,ni+2)
         call pdbatom (' C5 ',resname,ires,ni+3)
         call pdbatom (' C6 ',resname,ires,ni+4)
         call pdbatom (' N6 ',resname,ires,ni+5)
         call pdbatom (' N1 ',resname,ires,ni+6)
         call pdbatom (' C2 ',resname,ires,ni+7)
         call pdbatom (' N3 ',resname,ires,ni+8)
         call pdbatom (' C4 ',resname,ires,ni+9)
c
c     guanine in guanosine residue  (G)
c
      else if (resname .eq. 'G  ') then
         call pdbatom (' N9 ',resname,ires,ni)
         call pdbatom (' C8 ',resname,ires,ni+1)
         call pdbatom (' N7 ',resname,ires,ni+2)
         call pdbatom (' C5 ',resname,ires,ni+3)
         call pdbatom (' C6 ',resname,ires,ni+4)
         call pdbatom (' O6 ',resname,ires,ni+5)
         call pdbatom (' N1 ',resname,ires,ni+6)
         call pdbatom (' C2 ',resname,ires,ni+7)
         call pdbatom (' N2 ',resname,ires,ni+8)
         call pdbatom (' N3 ',resname,ires,ni+9)
         call pdbatom (' C4 ',resname,ires,ni+10)
c
c     cytosine in cytidine residue  (C)
c
      else if (resname .eq. 'C  ') then
         call pdbatom (' N1 ',resname,ires,ni)
         call pdbatom (' C2 ',resname,ires,ni+1)
         call pdbatom (' O2 ',resname,ires,ni+2)
         call pdbatom (' N3 ',resname,ires,ni+3)
         call pdbatom (' C4 ',resname,ires,ni+4)
         call pdbatom (' N4 ',resname,ires,ni+5)
         call pdbatom (' C5 ',resname,ires,ni+6)
         call pdbatom (' C6 ',resname,ires,ni+7)
c
c     uracil in uridine residue  (U)
c
      else if (resname .eq. 'U  ') then
         call pdbatom (' N1 ',resname,ires,ni)
         call pdbatom (' C2 ',resname,ires,ni+1)
         call pdbatom (' O2 ',resname,ires,ni+2)
         call pdbatom (' N3 ',resname,ires,ni+3)
         call pdbatom (' C4 ',resname,ires,ni+4)
         call pdbatom (' O4 ',resname,ires,ni+5)
         call pdbatom (' C5 ',resname,ires,ni+6)
         call pdbatom (' C6 ',resname,ires,ni+7)
c
c     adenine in deoxyadenosine residue  (DA)
c
      else if (resname .eq. 'DA ') then
         call pdbatom (' N9 ',resname,ires,ni)
         call pdbatom (' C8 ',resname,ires,ni+1)
         call pdbatom (' N7 ',resname,ires,ni+2)
         call pdbatom (' C5 ',resname,ires,ni+3)
         call pdbatom (' C6 ',resname,ires,ni+4)
         call pdbatom (' N6 ',resname,ires,ni+5)
         call pdbatom (' N1 ',resname,ires,ni+6)
         call pdbatom (' C2 ',resname,ires,ni+7)
         call pdbatom (' N3 ',resname,ires,ni+8)
         call pdbatom (' C4 ',resname,ires,ni+9)
c
c     guanine in deoxyguanosine residue  (DG)
c
      else if (resname .eq. 'DG ') then
         call pdbatom (' N9 ',resname,ires,ni)
         call pdbatom (' C8 ',resname,ires,ni+1)
         call pdbatom (' N7 ',resname,ires,ni+2)
         call pdbatom (' C5 ',resname,ires,ni+3)
         call pdbatom (' C6 ',resname,ires,ni+4)
         call pdbatom (' O6 ',resname,ires,ni+5)
         call pdbatom (' N1 ',resname,ires,ni+6)
         call pdbatom (' C2 ',resname,ires,ni+7)
         call pdbatom (' N2 ',resname,ires,ni+8)
         call pdbatom (' N3 ',resname,ires,ni+9)
         call pdbatom (' C4 ',resname,ires,ni+10)
c
c     cytosine in deoxycytidine residue  (DC)
c
      else if (resname .eq. 'DC ') then
         call pdbatom (' N1 ',resname,ires,ni)
         call pdbatom (' C2 ',resname,ires,ni+1)
         call pdbatom (' O2 ',resname,ires,ni+2)
         call pdbatom (' N3 ',resname,ires,ni+3)
         call pdbatom (' C4 ',resname,ires,ni+4)
         call pdbatom (' N4 ',resname,ires,ni+5)
         call pdbatom (' C5 ',resname,ires,ni+6)
         call pdbatom (' C6 ',resname,ires,ni+7)
c
c     thymine in deoxythymidine residue  (DT)
c
      else if (resname .eq. 'DT ') then
         call pdbatom (' N1 ',resname,ires,ni)
         call pdbatom (' C2 ',resname,ires,ni+1)
         call pdbatom (' O2 ',resname,ires,ni+2)
         call pdbatom (' N3 ',resname,ires,ni+3)
         call pdbatom (' C4 ',resname,ires,ni+4)
         call pdbatom (' O4 ',resname,ires,ni+5)
         call pdbatom (' C5 ',resname,ires,ni+6)
         call pdbatom (' C7 ',resname,ires,ni+7)
         call pdbatom (' C6 ',resname,ires,ni+8)
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine getnuch  --  extract the nucleotide hydrogens  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "getnuch" finds the nucleotide hydrogen atoms for a single
c     residue and copies the names and coordinates to the Protein
c     Data Bank file
c
c
      subroutine getnuch (resname,ires,ni,c1i,c2i,
     &                    o2i,c3i,o3i,c4i,c5i,o5i)
      use sizes
      use atomid
      use couple
      implicit none
      integer i,k
      integer ires,ni
      integer c1i,c4i
      integer c2i,o2i
      integer c3i,o3i
      integer c5i,o5i
      logical allatom,done
      character*3 resname
c
c
c     if no ribose C1 hydrogen, then united atom force field
c
      allatom = .true.
      if (n12(c1i) .ne. 4)  allatom = .false.
c
c     get sugar ring hydrogen atoms for the current residue
c
      done = .false.
      do i = 1, n12(c5i)
         k = i12(i,c5i)
         if (atomic(k).eq.1) then
            if (.not. done) then
               call pdbatom (' H5''',resname,ires,k)
               done = .true.
            else
               call pdbatom ('H5''''',resname,ires,k)
            end if
         end if
      end do
      do i = 1, n12(c4i)
         k = i12(i,c4i)
         if (atomic(k) .eq. 1)  call pdbatom (' H4''',resname,ires,k)
      end do
      do i = 1, n12(c3i)
         k = i12(i,c3i)
         if (atomic(k) .eq. 1)  call pdbatom (' H3''',resname,ires,k)
      end do
      done = .false.
      do i = 1, n12(c2i)
         k = i12(i,c2i)
         if (atomic(k) .eq. 1) then
            if (.not. done) then
               call pdbatom (' H2''',resname,ires,k)
               done = .true.
            else
               call pdbatom ('H2''''',resname,ires,k)
            end if
         end if
      end do
      if (o2i .ne. 0) then
         do i = 1, n12(o2i)
            k = i12(i,o2i)
            if (atomic(k) .eq. 1)  call pdbatom ('HO2''',resname,ires,k)
         end do
      end if
      do i = 1, n12(c1i)
         k = i12(i,c1i)
         if (atomic(k) .eq. 1)  call pdbatom (' H1''',resname,ires,k)
      end do
c
c     adenine in adenosine residue  (A)
c
      if (resname .eq. 'A  ') then
         if (allatom) then
            call pdbatom (' H8 ',resname,ires,ni+10)
            call pdbatom (' H61',resname,ires,ni+11)
            call pdbatom (' H62',resname,ires,ni+12)
            call pdbatom (' H2 ',resname,ires,ni+13)
         else
            call pdbatom (' H61',resname,ires,ni+10)
            call pdbatom (' H62',resname,ires,ni+11)
         end if
c
c     guanine in guanosine residue  (G)
c
      else if (resname .eq. 'G  ') then
         if (allatom) then
            call pdbatom (' H8 ',resname,ires,ni+11)
            call pdbatom (' H1 ',resname,ires,ni+12)
            call pdbatom (' H21',resname,ires,ni+13)
            call pdbatom (' H22',resname,ires,ni+14)
         else
            call pdbatom (' H1 ',resname,ires,ni+11)
            call pdbatom (' H21',resname,ires,ni+12)
            call pdbatom (' H22',resname,ires,ni+13)
         end if
c
c     cytosine in cytidine residue  (C)
c
      else if (resname .eq. 'C  ') then
         if (allatom) then
            call pdbatom (' H41',resname,ires,ni+8)
            call pdbatom (' H42',resname,ires,ni+9)
            call pdbatom (' H5 ',resname,ires,ni+10)
            call pdbatom (' H6 ',resname,ires,ni+11)
         else
            call pdbatom (' H41',resname,ires,ni+8)
            call pdbatom (' H42',resname,ires,ni+9)
         end if
c
c     uracil in uridine residue  (U)
c
      else if (resname .eq. 'U  ') then
         if (allatom) then
            call pdbatom (' H3 ',resname,ires,ni+8)
            call pdbatom (' H5 ',resname,ires,ni+9)
            call pdbatom (' H6 ',resname,ires,ni+10)
         else
            call pdbatom (' H3 ',resname,ires,ni+8)
         end if
c
c     adenine in deoxyadenosine residue  (DA)
c
      else if (resname .eq. 'DA ') then
         if (allatom) then
            call pdbatom (' H8 ',resname,ires,ni+10)
            call pdbatom (' H61',resname,ires,ni+11)
            call pdbatom (' H62',resname,ires,ni+12)
            call pdbatom (' H2 ',resname,ires,ni+13)
         else
            call pdbatom (' H61',resname,ires,ni+10)
            call pdbatom (' H62',resname,ires,ni+11)
         end if
c
c     guanine in deoxyguanosine residue  (DG)
c
      else if (resname .eq. 'DG ') then
         if (allatom) then
            call pdbatom (' H8 ',resname,ires,ni+11)
            call pdbatom (' H1 ',resname,ires,ni+12)
            call pdbatom (' H21',resname,ires,ni+13)
            call pdbatom (' H22',resname,ires,ni+14)
         else
            call pdbatom (' H1 ',resname,ires,ni+11)
            call pdbatom (' H21',resname,ires,ni+12)
            call pdbatom (' H22',resname,ires,ni+13)
         end if
c
c     cytosine in deoxycytidine residue  (DC)
c
      else if (resname .eq. 'DC ') then
         if (allatom) then
            call pdbatom (' H41',resname,ires,ni+8)
            call pdbatom (' H42',resname,ires,ni+9)
            call pdbatom (' H5 ',resname,ires,ni+10)
            call pdbatom (' H6 ',resname,ires,ni+11)
         else
            call pdbatom (' H41',resname,ires,ni+8)
            call pdbatom (' H42',resname,ires,ni+9)
         end if
c
c     thymine in deoxythymidine residue  (DT)
c
      else if (resname .eq. 'DT ') then
         if (allatom) then
            call pdbatom (' H3 ',resname,ires,ni+9)
            call pdbatom (' H71',resname,ires,ni+10)
            call pdbatom (' H72',resname,ires,ni+11)
            call pdbatom (' H73',resname,ires,ni+12)
            call pdbatom (' H6 ',resname,ires,ni+13)
         else
            call pdbatom (' H3 ',resname,ires,ni+9)
         end if
      end if
c
c     get any capping hydrogen atoms for the current residue
c
      do i = 1, n12(o5i)
         k = i12(i,o5i)
         if (atomic(k) .eq. 1)  call pdbatom (' H5T',resname,ires,k)
      end do
      do i = 1, n12(o3i)
         k = i12(i,o3i)
         if (atomic(k) .eq. 1)  call pdbatom (' H3T',resname,ires,k)
      end do
      return
      end

      function computegrad (xx,g)
      use polar
      use mpole
      implicit none
      integer i,j,jj,nvar,ncart,limit
      real*8 computegrad
      real*8 xx(0:ptmaxord), g(0:ptmaxord)
      real*8 c
      real*8, allocatable :: difd(:,:), difp(:,:)

      nvar = ptmaxord+1
      ncart = 3*npole
      allocate(difd(3,npole))
      if (includeuinp .ne. 0) allocate(difp(3,npole))

c     difference is just uind - ptuind
      difd = uind
      if (includeuinp .ne. 0) difp = uinp

      do i = 0,ptmaxord
        do j = 1,npole
          difd(:,j) = difd(:,j) - xx(i)*ptuind(:,j,i)
        enddo
        if (includeuinp .ne. 0) then
           do j = 1,npole
             difp(:,j) = difp(:,j) - xx(i)*ptuinp(:,j,i)
           enddo
        endif
      enddo

c     the residual is simply |dif|**fitpow
      computegrad = 0d0
      do i = 1,npole
        computegrad = computegrad + 
     &     (difd(1,i)**fitpow + difd(2,i)**fitpow + difd(3,i)**fitpow)
      enddo
      if(includeuinp .ne. 0) then
        do i = 1,npole
          computegrad = computegrad + 
     &       (difp(1,i)**fitpow + difp(2,i)**fitpow + difp(3,i)**fitpow)
        enddo
      endif

      do i = 0,ptmaxord
        g(i) = 0d0
        do j = 1,npole
          g(i) = g(i) - DBLE(fitpow)*(
     &             ptuind(1,j,i)*difd(1,j)**(fitpow-1)
     &            +ptuind(2,j,i)*difd(2,j)**(fitpow-1)
     &            +ptuind(3,j,i)*difd(3,j)**(fitpow-1))
        enddo
        if (includeuinp .ne. 0) then
          do j = 1,npole
            g(i) = g(i) - DBLE(fitpow)*(
     &                  ptuinp(1,j,i)*difp(1,j)**(fitpow-1)
     &                 +ptuinp(2,j,i)*difp(2,j)**(fitpow-1)
     &                 +ptuinp(3,j,i)*difp(3,j)**(fitpow-1))
          enddo
        endif
      enddo

      deallocate(difd)
      if (includeuinp .ne. 0) deallocate(difp)
      return
      end


      function computehess (mode,xx,h,hinit,hstop,hindex,hdiag)
      use polar
      use mpole
      implicit none
      integer hinit(*)
      integer hstop(*)
      integer hindex(*)
      integer i,j,k,ncart,addr,nvar
      real*8, allocatable :: difd(:,:), difp(:,:)
      real*8 computehess
      real*8 xx(*)
      real*8 hdiag(*)
      real*8 h(*)
      character*4 mode

      if (mode .eq. 'NONE')  return

      nvar = ptmaxord+1

      allocate(difd(3,npole))
      if (includeuinp .ne. 0) allocate(difp(3,npole))

c     difference is just uind - ptuind
      difd = uind
      if (includeuinp .ne. 0) difp = uinp

      do i = 0,ptmaxord
        do j = 1,npole
          difd(:,j) = difd(:,j) - xx(i)*ptuind(:,j,i)
        enddo
        if (includeuinp .ne. 0) then
           do j = 1,npole
             difp(:,j) = difp(:,j) - xx(i)*ptuinp(:,j,i)
           enddo
        endif
      enddo

c     compute the diagonals of the hessian matrix
      do i = 1,nvar
        hdiag(i) = 0d0
        do k = 1,npole
          hdiag(i) = hdiag(i)-DBLE(fitpow*(fitpow-1))*ptuind(1,k,i-1)
     &                        *ptuind(1,k,i-1)*difd(1,k)**(fitpow-2)
     &                       -DBLE(fitpow*(fitpow-1))*ptuind(2,k,i-1)
     &                        *ptuind(2,k,i-1)*difd(2,k)**(fitpow-2)
     &                       -DBLE(fitpow*(fitpow-1))*ptuind(3,k,i-1)
     &                        *ptuind(3,k,i-1)*difd(3,k)**(fitpow-2)
        enddo
        if (includeuinp .ne. 0) then
          do k = 1,npole
            hdiag(i) = hdiag(i)-DBLE(fitpow*(fitpow-1))*ptuinp(1,k,i-1)
     &                          *ptuinp(1,k,i-1)*difp(1,k)**(fitpow-2)
     &                         -DBLE(fitpow*(fitpow-1))*ptuinp(2,k,i-1)
     &                          *ptuinp(2,k,i-1)*difp(2,k)**(fitpow-2)
     &                         -DBLE(fitpow*(fitpow-1))*ptuinp(3,k,i-1)
     &                          *ptuinp(3,k,i-1)*difp(3,k)**(fitpow-2)
          enddo
        endif
      enddo


      if (mode .eq. 'DIAG')  return

c     now compute the off-diagonals of the hessian matrix
      addr = 1
      do i = 1,nvar
        hinit(i) = addr
        do j = i+1,nvar
          h(addr) = 0d0
          do k = 1,npole
            h(addr) = h(addr)-DBLE(fitpow*(fitpow-1))*ptuind(1,k,i-1)
     &                        *ptuind(1,k,j-1)*difd(1,k)**(fitpow-2)
     &                       -DBLE(fitpow*(fitpow-1))*ptuind(2,k,i-1)
     &                        *ptuind(2,k,j-1)*difd(2,k)**(fitpow-2)
     &                       -DBLE(fitpow*(fitpow-1))*ptuind(3,k,i-1)
     &                        *ptuind(3,k,j-1)*difd(3,k)**(fitpow-2)
          enddo
          if (includeuinp .ne. 0) then
            do k = 1,npole
              h(addr) = h(addr)-DBLE(fitpow*(fitpow-1))*ptuinp(1,k,i-1)
     &                          *ptuinp(1,k,j-1)*difp(1,k)**(fitpow-2)
     &                         -DBLE(fitpow*(fitpow-1))*ptuinp(2,k,i-1)
     &                          *ptuinp(2,k,j-1)*difp(2,k)**(fitpow-2)
     &                         -DBLE(fitpow*(fitpow-1))*ptuinp(3,k,i-1)
     &                          *ptuinp(3,k,j-1)*difp(3,k)**(fitpow-2)
            enddo
          endif
          hindex(addr) = j
          addr = addr + 1
        enddo
        hstop(i) = addr-1
      enddo

      deallocate(difd)
      if (includeuinp .ne. 0) deallocate(difp)
      return

      end

      subroutine savestate (ncycle,f,xx)
      use polar
      implicit none
      integer i,j,iopt,iend
      integer ncycle,nvar
      integer lext,freeunit
      real*8 f,xx(0:ptmaxord+1)
      logical exist
      character*7 ext
      character*120 optfile
      character*120 endfile
c      
c     map the partial coefficients back to full coefficients
c
      do i = ptmaxord,0,-1
         ptcoefsf(i) = xx(i)
         do j = i+1,ptmaxord
            ptcoefsf(i) = ptcoefsf(i) - ptcoefsf(j)
         enddo
      enddo
c
c     test for requested termination of the optimization
c
CCCC      endfile = 'tinker.end'
CCCC      inquire (file=endfile,exist=exist)
CCCC      if (.not. exist) then
CCCC         endfile = filename(1:leng)//'.end'
CCCC         inquire (file=endfile,exist=exist)
CCCC         if (exist) then
CCCC            iend = freeunit ()
CCCC            open (unit=iend,file=endfile,status='old')
CCCC            close (unit=iend,status='delete')
CCCC         end if
CCCC      end if
CCCC      if (exist) then
CCCC         write (iout,10)
CCCC   10    format (/,' OPTSAVE  --  Optimization Calculation Ending',
CCCC     &              ' due to User Request')
CCCC         call fatal
CCCC      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine dfield0a  --  direct induction via double loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dfield0a" computes the direct electrostatic field due to
c     permanent multipole moments via a double loop
c
c
      subroutine dfield0a (field,fieldp)
      use sizes
      use atoms
      use bound
      use cell
      use couple
      use group
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr3,rr5,rr7
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      logical proceed
      character*6 mode
c
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
c
c     find the electrostatic field due to permanent multipoles
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = i+1, npole
            dscale(ipole(j)) = 1.0d0
            pscale(ipole(j)) = 1.0d0
         end do
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = 1.0d0 - expdamp
                        scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                        scale7 = 1.0d0 - expdamp
     &                              *(1.0d0-damp+0.6d0*damp**2)
                     end if
                  end if
                  rr3 = scale3 / (r*r2)
                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
                  rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
                  dir = dix*xr + diy*yr + diz*zr
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qir = qix*xr + qiy*yr + qiz*zr
                  dkr = dkx*xr + dky*yr + dkz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
                  qkr = qkx*xr + qky*yr + qkz*zr
                  fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkx + 2.0d0*rr5*qkx
                  fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dky + 2.0d0*rr5*qky
                  fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkz + 2.0d0*rr5*qkz
                  fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*dix - 2.0d0*rr5*qix
                  fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diy - 2.0d0*rr5*qiy
                  fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diz - 2.0d0*rr5*qiz
                  do j = 1, 3
                     field(j,i) = field(j,i) + fid(j)*dscale(kk)
                     field(j,k) = field(j,k) + fkd(j)*dscale(kk)
                     fieldp(j,i) = fieldp(j,i) + fid(j)*pscale(kk)
                     fieldp(j,k) = fieldp(j,k) + fkd(j)*pscale(kk)
                  end do
               end if
            end if
         end do
      end do
c
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
         do i = 1, npole
            ii = ipole(i)
            pdi = pdamp(i)
            pti = thole(i)
            ci = rpole(1,i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            qixx = rpole(5,i)
            qixy = rpole(6,i)
            qixz = rpole(7,i)
            qiyy = rpole(9,i)
            qiyz = rpole(10,i)
            qizz = rpole(13,i)
            do j = i, npole
               dscale(ipole(j)) = 1.0d0
               pscale(ipole(j)) = 1.0d0
            end do
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = p2scale
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = p3scale
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = p4scale
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = p5scale
            end do
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = d1scale
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = d2scale
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = d3scale
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = d4scale
            end do
            do k = i, npole
               kk = ipole(k)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
               do m = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call imager (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     scale7 = 1.0d0
                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(k))
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           expdamp = exp(damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                           scale7 = 1.0d0 - expdamp
     &                                 *(1.0d0-damp+0.6d0*damp**2)
                        end if
                     end if
                     rr3 = scale3 / (r*r2)
                     rr5 = 3.0d0 * scale5 / (r*r2*r2)
                     rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
                     dir = dix*xr + diy*yr + diz*zr
                     qix = qixx*xr + qixy*yr + qixz*zr
                     qiy = qixy*xr + qiyy*yr + qiyz*zr
                     qiz = qixz*xr + qiyz*yr + qizz*zr
                     qir = qix*xr + qiy*yr + qiz*zr
                     dkr = dkx*xr + dky*yr + dkz*zr
                     qkx = qkxx*xr + qkxy*yr + qkxz*zr
                     qky = qkxy*xr + qkyy*yr + qkyz*zr
                     qkz = qkxz*xr + qkyz*yr + qkzz*zr
                     qkr = qkx*xr + qky*yr + qkz*zr
                     fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkx + 2.0d0*rr5*qkx
                     fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dky + 2.0d0*rr5*qky
                     fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkz + 2.0d0*rr5*qkz
                     fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*dix - 2.0d0*rr5*qix
                     fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diy - 2.0d0*rr5*qiy
                     fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diz - 2.0d0*rr5*qiz
                     do j = 1, 3
                        fip(j) = fid(j)
                        fkp(j) = fkd(j)
                     end do
                     if (use_polymer .and. r2 .le. polycut2) then
                        do j = 1, 3
                           fid(j) = fid(j) * dscale(kk)
                           fip(j) = fip(j) * pscale(kk)
                           fkd(j) = fkd(j) * dscale(kk)
                           fkp(j) = fkp(j) * pscale(kk)
                        end do
                     end if
                     do j = 1, 3
                        field(j,i) = field(j,i) + fid(j)
                        fieldp(j,i) = fieldp(j,i) + fip(j)
                        if (ii .ne. kk) then
                           field(j,k) = field(j,k) + fkd(j)
                           fieldp(j,k) = fieldp(j,k) + fkp(j)
                        end if
                     end do
                  end if
               end do
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ufield0a  --  mutual induction via double loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ufield0a" computes the mutual electrostatic field due to
c     induced dipole moments via a double loop
c
c
      subroutine ufield0a (field,fieldp)
      use sizes
      use atoms
      use bound
      use cell
      use group
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr3,rr5
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 duir,puir
      real*8 dukr,pukr
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8, allocatable :: dscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      logical proceed
      character*6 mode
c
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
c
c     find the electrostatic field due to mutual induced dipoles
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         duix = uind(1,i)
         duiy = uind(2,i)
         duiz = uind(3,i)
         puix = uinp(1,i)
         puiy = uinp(2,i)
         puiz = uinp(3,i)
         do j = i+1, npole
            dscale(ipole(j)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = u4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  dukx = uind(1,k)
                  duky = uind(2,k)
                  dukz = uind(3,k)
                  pukx = uinp(1,k)
                  puky = uinp(2,k)
                  pukz = uinp(3,k)
                  scale3 = dscale(kk)
                  scale5 = dscale(kk)
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = scale3 * (1.0d0-expdamp)
                        scale5 = scale5 * (1.0d0-expdamp*(1.0d0-damp))
                     end if
                  end if
                  rr3 = scale3 / (r*r2)
                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
                  duir = xr*duix + yr*duiy + zr*duiz
                  dukr = xr*dukx + yr*duky + zr*dukz
                  puir = xr*puix + yr*puiy + zr*puiz
                  pukr = xr*pukx + yr*puky + zr*pukz
                  fid(1) = -rr3*dukx + rr5*dukr*xr
                  fid(2) = -rr3*duky + rr5*dukr*yr
                  fid(3) = -rr3*dukz + rr5*dukr*zr
                  fkd(1) = -rr3*duix + rr5*duir*xr
                  fkd(2) = -rr3*duiy + rr5*duir*yr
                  fkd(3) = -rr3*duiz + rr5*duir*zr
                  fip(1) = -rr3*pukx + rr5*pukr*xr
                  fip(2) = -rr3*puky + rr5*pukr*yr
                  fip(3) = -rr3*pukz + rr5*pukr*zr
                  fkp(1) = -rr3*puix + rr5*puir*xr
                  fkp(2) = -rr3*puiy + rr5*puir*yr
                  fkp(3) = -rr3*puiz + rr5*puir*zr
                  do j = 1, 3
                     field(j,i) = field(j,i) + fid(j)
                     field(j,k) = field(j,k) + fkd(j)
                     fieldp(j,i) = fieldp(j,i) + fip(j)
                     fieldp(j,k) = fieldp(j,k) + fkp(j)
                  end do
               end if
            end if
         end do
      end do
c
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
         do i = 1, npole
            ii = ipole(i)
            pdi = pdamp(i)
            pti = thole(i)
            duix = uind(1,i)
            duiy = uind(2,i)
            duiz = uind(3,i)
            puix = uinp(1,i)
            puiy = uinp(2,i)
            puiz = uinp(3,i)
            do j = i, npole
               dscale(ipole(j)) = 1.0d0
            end do
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = u1scale
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = u2scale
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = u3scale
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = u4scale
            end do
            do k = i, npole
               kk = ipole(k)
               dukx = uind(1,k)
               duky = uind(2,k)
               dukz = uind(3,k)
               pukx = uinp(1,k)
               puky = uinp(2,k)
               pukz = uinp(3,k)
               proceed = .true.
               do m = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call imager (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(k))
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           expdamp = exp(damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                        end if
                     end if
                     rr3 = scale3 / (r*r2)
                     rr5 = 3.0d0 * scale5 / (r*r2*r2)
                     duir = xr*duix + yr*duiy + zr*duiz
                     dukr = xr*dukx + yr*duky + zr*dukz
                     puir = xr*puix + yr*puiy + zr*puiz
                     pukr = xr*pukx + yr*puky + zr*pukz
                     fid(1) = -rr3*dukx + rr5*dukr*xr
                     fid(2) = -rr3*duky + rr5*dukr*yr
                     fid(3) = -rr3*dukz + rr5*dukr*zr
                     fkd(1) = -rr3*duix + rr5*duir*xr
                     fkd(2) = -rr3*duiy + rr5*duir*yr
                     fkd(3) = -rr3*duiz + rr5*duir*zr
                     fip(1) = -rr3*pukx + rr5*pukr*xr
                     fip(2) = -rr3*puky + rr5*pukr*yr
                     fip(3) = -rr3*pukz + rr5*pukr*zr
                     fkp(1) = -rr3*puix + rr5*puir*xr
                     fkp(2) = -rr3*puiy + rr5*puir*yr
                     fkp(3) = -rr3*puiz + rr5*puir*zr
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           do j = 1, 3
                              fid(j) = fid(j) * dscale(kk)
                              fkd(j) = fkd(j) * dscale(kk)
                              fip(j) = fip(j) * dscale(kk)
                              fkp(j) = fkp(j) * dscale(kk)
                           end do
                        end if
                     end if
                     do j = 1, 3
                        field(j,i) = field(j,i) + fid(j)
                        fieldp(j,i) = fieldp(j,i) + fip(j)
                        if (ii .ne. kk) then
                           field(j,k) = field(j,k) + fkd(j)
                           fieldp(j,k) = fieldp(j,k) + fkp(j)
                        end if
                     end do
                  end if
               end do
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine dfield0b  --  direct induction via pair list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "dfield0b" computes the mutual electrostatic field due to
c     permanent multipole moments via a pair list
c
c
      subroutine dfield0b (field,fieldp)
      use sizes
      use atoms
      use bound
      use couple
      use group
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr3,rr5,rr7
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      logical proceed
      character*6 mode
c
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
c
c     find the electrostatic field due to permanent multipoles
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = 1, nelst(i)
            dscale(ipole(elst(j,i))) = 1.0d0
            pscale(ipole(elst(j,i))) = 1.0d0
         end do
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = 1.0d0 - expdamp
                        scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                        scale7 = 1.0d0 - expdamp
     &                              *(1.0d0-damp+0.6d0*damp**2)
                     end if
                  end if
                  rr3 = scale3 / (r*r2)
                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
                  rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
                  dir = dix*xr + diy*yr + diz*zr
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qir = qix*xr + qiy*yr + qiz*zr
                  dkr = dkx*xr + dky*yr + dkz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
                  qkr = qkx*xr + qky*yr + qkz*zr
                  fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkx + 2.0d0*rr5*qkx
                  fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dky + 2.0d0*rr5*qky
                  fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkz + 2.0d0*rr5*qkz
                  fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*dix - 2.0d0*rr5*qix
                  fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diy - 2.0d0*rr5*qiy
                  fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diz - 2.0d0*rr5*qiz
                  do j = 1, 3
                     field(j,i) = field(j,i) + fid(j)*dscale(kk)
                     field(j,k) = field(j,k) + fkd(j)*dscale(kk)
                     fieldp(j,i) = fieldp(j,i) + fid(j)*pscale(kk)
                     fieldp(j,k) = fieldp(j,k) + fkd(j)*pscale(kk)
                  end do
               end if
            end if
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ufield0b  --  mutual induction via pair list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ufield0b" computes the mutual electrostatic field due to
c     induced dipole moments via a pair list
c
c
      subroutine ufield0b (field,fieldp)
      use sizes
      use atoms
      use bound
      use group
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr3,rr5
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 duir,puir
      real*8 dukr,pukr
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8, allocatable :: dscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      logical proceed
      character*6 mode
c
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
c
c     find the electrostatic field due to mutual induced dipoles
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         duix = uind(1,i)
         duiy = uind(2,i)
         duiz = uind(3,i)
         puix = uinp(1,i)
         puiy = uinp(2,i)
         puiz = uinp(3,i)
         do j = 1, nelst(i)
            dscale(ipole(elst(j,i))) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = u4scale
         end do
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  dukx = uind(1,k)
                  duky = uind(2,k)
                  dukz = uind(3,k)
                  pukx = uinp(1,k)
                  puky = uinp(2,k)
                  pukz = uinp(3,k)
                  scale3 = dscale(kk)
                  scale5 = dscale(kk)
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = scale3 * (1.0d0-expdamp)
                        scale5 = scale5 * (1.0d0-(1.0d0-damp)*expdamp)
                     end if
                  end if
                  rr3 = scale3 / (r*r2)
                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
                  duir = xr*duix + yr*duiy + zr*duiz
                  dukr = xr*dukx + yr*duky + zr*dukz
                  puir = xr*puix + yr*puiy + zr*puiz
                  pukr = xr*pukx + yr*puky + zr*pukz
                  fid(1) = -rr3*dukx + rr5*dukr*xr
                  fid(2) = -rr3*duky + rr5*dukr*yr
                  fid(3) = -rr3*dukz + rr5*dukr*zr
                  fkd(1) = -rr3*duix + rr5*duir*xr
                  fkd(2) = -rr3*duiy + rr5*duir*yr
                  fkd(3) = -rr3*duiz + rr5*duir*zr
                  fip(1) = -rr3*pukx + rr5*pukr*xr
                  fip(2) = -rr3*puky + rr5*pukr*yr
                  fip(3) = -rr3*pukz + rr5*pukr*zr
                  fkp(1) = -rr3*puix + rr5*puir*xr
                  fkp(2) = -rr3*puiy + rr5*puir*yr
                  fkp(3) = -rr3*puiz + rr5*puir*zr
                  do j = 1, 3
                     field(j,i) = field(j,i) + fid(j)
                     field(j,k) = field(j,k) + fkd(j)
                     fieldp(j,i) = fieldp(j,i) + fip(j)
                     fieldp(j,k) = fieldp(j,k) + fkp(j)
                  end do
               end if
            end if
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine dfield0c  --  direct induction via Ewald sum  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "dfield0c" computes the mutual electrostatic field due to
c     permanent multipole moments via Ewald summation
c
c
      subroutine dfield0c (field,fieldp)
      use sizes
      use atoms
      use boxes
      use ewald
      use limits
      use math
      use mpole
      use polar
      implicit none
      integer i,j,ii
      real*8 term
      real*8 ucell(3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
c
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do
c
c     get the reciprocal space part of the electrostatic field
c
      call udirect1 (field)
      do i = 1, npole
         do j = 1, 3
            fieldp(j,i) = field(j,i)
         end do
      end do
c
c     get the real space portion of the electrostatic field
c
      if (use_mlist) then
         call udirect2b (field,fieldp)
      else
         call udirect2a (field,fieldp)
      end if
c
c     get the self-energy portion of the electrostatic field
c
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do i = 1, npole
         do j = 1, 3
            field(j,i) = field(j,i) + term*rpole(j+1,i)
            fieldp(j,i) = fieldp(j,i) + term*rpole(j+1,i)
         end do
      end do
c
c     compute the cell dipole boundary correction to field
c
      if (boundary .eq. 'VACUUM') then
         do i = 1, 3
            ucell(i) = 0.0d0
         end do
         do i = 1, npole
            ii = ipole(i)
            ucell(1) = ucell(1) + rpole(2,i) + rpole(1,i)*x(ii)
            ucell(2) = ucell(2) + rpole(3,i) + rpole(1,i)*y(ii)
            ucell(3) = ucell(3) + rpole(4,i) + rpole(1,i)*z(ii)
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do i = 1, npole
            do j = 1, 3
               field(j,i) = field(j,i) - term*ucell(j)
               fieldp(j,i) = fieldp(j,i) - term*ucell(j)
            end do
         end do
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ufield0c  --  mutual induction via Ewald sum  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ufield0c" computes the mutual electrostatic field due to
c     induced dipole moments via Ewald summation
c
c
      subroutine ufield0c (field,fieldp)
      use sizes
      use atoms
      use boxes
      use ewald
      use limits
      use math
      use mpole
      use polar
      implicit none
      integer i,j
      real*8 term
      real*8 ucell(3)
      real*8 ucellp(3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
c
c
c     zero out the electrostatic field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do
c
c     get the reciprocal space part of the electrostatic field
c
      call umutual1 (field,fieldp)
c
c     get the real space portion of the electrostatic field
c
      if (use_mlist) then
         call umutual2b (field,fieldp)
      else
         call umutual2a (field,fieldp)
      end if
c
c     get the self-energy portion of the electrostatic field
c
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do i = 1, npole
         do j = 1, 3
            field(j,i) = field(j,i) + term*uind(j,i)
            fieldp(j,i) = fieldp(j,i) + term*uinp(j,i)
         end do
      end do
c
c     compute the cell dipole boundary correction to the field
c
      if (boundary .eq. 'VACUUM') then
         do i = 1, 3
            ucell(i) = 0.0d0
            ucellp(i) = 0.0d0
         end do
         do i = 1, npole
            do j = 1, 3
               ucell(j) = ucell(j) + uind(j,i)
               ucellp(j) = ucellp(j) + uinp(j,i)
            end do
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do i = 1, npole
            do j = 1, 3
               field(j,i) = field(j,i) - term*ucell(j)
               fieldp(j,i) = fieldp(j,i) - term*ucellp(j)
            end do
         end do
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine udirect1  --  Ewald recip direct induced field  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "udirect1" computes the reciprocal space contribution of the
c     permanent atomic multipole moments to the field
c
c
      subroutine udirect1 (field)
      use sizes
      use bound
      use boxes
      use ewald
      use math
      use mpole
      use pme
c OPT IMPLEMENTATION
      use polar
c OPT IMPLEMENTATION
      implicit none
      integer i,j,k,ntot
      integer k1,k2,k3
      integer m1,m2,m3
      integer nff,nf1,nf2,nf3
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 field(3,*)
      real*8, allocatable :: cmp(:,:)
      real*8, allocatable :: fmp(:,:)
      real*8, allocatable :: cphi(:,:)
      real*8, allocatable :: fphi(:,:)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cmp(10,npole))
      allocate (fmp(10,npole))
      allocate (cphi(10,npole))
      allocate (fphi(20,npole))
c
c     copy multipole moments and coordinates to local storage
c
      do i = 1, npole
         cmp(1,i) = rpole(1,i)
         cmp(2,i) = rpole(2,i)
         cmp(3,i) = rpole(3,i)
         cmp(4,i) = rpole(4,i)
         cmp(5,i) = rpole(5,i)
         cmp(6,i) = rpole(9,i)
         cmp(7,i) = rpole(13,i)
         cmp(8,i) = 2.0d0 * rpole(6,i)
         cmp(9,i) = 2.0d0 * rpole(7,i)
         cmp(10,i) = 2.0d0 * rpole(10,i)
      end do
c
c     compute B-spline coefficients and spatial decomposition
c
      call bspline_fill
      call table_fill
c
c     convert Cartesian multipoles to fractional coordinates
c
      call cmp_to_fmp (cmp,fmp)
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_mpole (fmp)
      call fftfront
c OPT IMPLEMENTATION
      permgridf = qgrid
c OPT IMPLEMENTATION
c
c     make the scalar summation over reciprocal lattice
c
      qfac(1,1,1) = 0.0d0
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      ntot = nfft1 * nfft2 * nfft3
      do i = 1, ntot-1
         k3 = i/nff + 1
         j = i - (k3-1)*nff
         k2 = j/nfft1 + 1
         k1 = j - (k2-1)*nfft1 + 1
         m1 = k1 - 1
         m2 = k2 - 1
         m3 = k3 - 1
         if (k1 .gt. nf1)  m1 = m1 - nfft1
         if (k2 .gt. nf2)  m2 = m2 - nfft2
         if (k3 .gt. nf3)  m3 = m3 - nfft3
         r1 = dble(m1)
         r2 = dble(m2)
         r3 = dble(m3)
         h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
         h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
         h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
         hsq = h1*h1 + h2*h2 + h3*h3
         term = -pterm * hsq
         expterm = 0.0d0
         if (term .gt. -50.0d0) then
            denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
            expterm = exp(term) / denom
            if (.not. use_bounds) then
               expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
            else if (octahedron) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
            end if
         end if
         qfac(k1,k2,k3) = expterm
      end do
c
c     account for the zeroth grid point for a finite system
c
      qfac(1,1,1) = 0.0d0
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         qfac(1,1,1) = expterm
      end if
c
c     complete the transformation of the PME grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               term = qfac(i,j,k)
               qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
               qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
            end do
         end do
      end do
c
c     perform 3-D FFT backward transform and get field
c
      call fftback
      call fphi_mpole (fphi)
c OPT IMPLEMENTATION
      fphiperm = fphi
c OPT IMPLEMENTATION
c
c     convert the field from fractional to Cartesian
c
      call fphi_to_cphi (fphi,cphi)
c
c     increment the field at each multipole site
c
      do i = 1, npole
         field(1,i) = field(1,i) - cphi(2,i)
         field(2,i) = field(2,i) - cphi(3,i)
         field(3,i) = field(3,i) - cphi(4,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cmp)
      deallocate (fmp)
      deallocate (cphi)
      deallocate (fphi)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine udirect2a  --  Ewald real direct field via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "udirect2a" computes the real space contribution of the permanent
c     atomic multipole moments to the field via a double loop
c
c
      subroutine udirect2a (field,fieldp)
      use sizes
      use atoms
      use boxes
      use bound
      use cell
      use couple
      use ewald
      use math
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      use units
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xr,yr,zr,r,r2
      real*8 rr1,rr2,rr3
      real*8 rr5,rr7
      real*8 erfc,bfac,exp2a
      real*8 ci,dix,diy,diz
      real*8 qixx,qiyy,qizz
      real*8 qixy,qixz,qiyz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkyy,qkzz
      real*8 qkxy,qkxz,qkyz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 ralpha,aefac
      real*8 aesq2,aesq2n
      real*8 pdi,pti,pgamma
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 bn(0:3),bcn(3)
      real*8 fimd(3),fkmd(3)
      real*8 fimp(3),fkmp(3)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      character*6 mode
      external erfc
c
c
c     check for multipoles and set cutoff coefficients
c
      if (npole .eq. 0)  return
      mode = 'EWALD'
      call switch (mode)
      aesq2 = 2.0 * aewald * aewald
      aesq2n = 0.0d0
      if (aewald .gt. 0.0d0)  aesq2n = 1.0d0 / (sqrtpi*aewald)
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
      end do
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. cut2) then
               r = sqrt(r2)
               rr1 = 1.0d0 / r
               rr2 = rr1 * rr1
               rr3 = rr2 * rr1
               rr5 = rr2 * rr3
               rr7 = rr2 * rr5
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
c
c     calculate the error function damping factors
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) * rr1
               exp2a = exp(-ralpha**2)
               aefac = aesq2n
               do j = 1, 3
                  bfac = dble(j+j-1)
                  aefac = aesq2 * aefac
                  bn(j) = (bfac*bn(j-1)+aefac*exp2a) * rr2
               end do
c
c     compute the polarization damping scale factors
c
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                     scale7 = 1.0d0 - expdamp
     &                           *(1.0d0-damp+0.6d0*damp**2)
                  end if
               end if
c
c     find the field terms for the current interaction
c
               dir = dix*xr + diy*yr + diz*zr
               qix = qixx*xr + qixy*yr + qixz*zr
               qiy = qixy*xr + qiyy*yr + qiyz*zr
               qiz = qixz*xr + qiyz*yr + qizz*zr
               qir = qix*xr + qiy*yr + qiz*zr
               dkr = dkx*xr + dky*yr + dkz*zr
               qkx = qkxx*xr + qkxy*yr + qkxz*zr
               qky = qkxy*xr + qkyy*yr + qkyz*zr
               qkz = qkxz*xr + qkyz*yr + qkzz*zr
               qkr = qkx*xr + qky*yr + qkz*zr
               bcn(1) = bn(1) - (1.0d0-scale3*dscale(kk))*rr3
               bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*dscale(kk))*rr5
               bcn(3) = bn(3) - 15.0d0*(1.0d0-scale7*dscale(kk))*rr7
               fimd(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
               fimd(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dky + 2.0d0*bcn(2)*qky
               fimd(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
               fkmd(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*dix - 2.0d0*bcn(2)*qix
               fkmd(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diy - 2.0d0*bcn(2)*qiy
               fkmd(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diz - 2.0d0*bcn(2)*qiz
               bcn(1) = bn(1) - (1.0d0-scale3*pscale(kk))*rr3
               bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*pscale(kk))*rr5
               bcn(3) = bn(3) - 15.0d0*(1.0d0-scale7*pscale(kk))*rr7
               fimp(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
               fimp(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dky + 2.0d0*bcn(2)*qky
               fimp(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
               fkmp(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*dix - 2.0d0*bcn(2)*qix
               fkmp(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diy - 2.0d0*bcn(2)*qiy
               fkmp(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diz - 2.0d0*bcn(2)*qiz
c
c     increment the field at each site due to this interaction
c
               do j = 1, 3
                  field(j,i) = field(j,i) + fimd(j)
                  field(j,k) = field(j,k) + fkmd(j)
                  fieldp(j,i) = fieldp(j,i) + fimp(j)
                  fieldp(j,k) = fieldp(j,k) + fkmp(j)
               end do
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
         do i = 1, npole
            ii = ipole(i)
            pdi = pdamp(i)
            pti = thole(i)
            ci = rpole(1,i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            qixx = rpole(5,i)
            qixy = rpole(6,i)
            qixz = rpole(7,i)
            qiyy = rpole(9,i)
            qiyz = rpole(10,i)
            qizz = rpole(13,i)
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = p2scale
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = p3scale
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = p4scale
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = p5scale
            end do
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = d1scale
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = d2scale
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = d3scale
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = d4scale
            end do
            do k = i, npole
               kk = ipole(k)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
               do m = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call imager (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
c
c     calculate the error function damping factors
c
                  if (r2 .le. cut2) then
                     r = sqrt(r2)
                     rr1 = 1.0d0 / r
                     rr2 = rr1 * rr1
                     rr3 = rr2 * rr1
                     rr5 = rr2 * rr3
                     rr7 = rr2 * rr5
                     ralpha = aewald * r
                     bn(0) = erfc(ralpha) * rr1
                     exp2a = exp(-ralpha**2)
                     aefac = aesq2n
                     do j = 1, 3
                        bfac = dble(j+j-1)
                        aefac = aesq2 * aefac
                        bn(j) = (bfac*bn(j-1)+aefac*exp2a) * rr2
                     end do
c
c     compute the polarization damping scale factors
c
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     scale7 = 1.0d0
                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(k))
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           expdamp = exp(damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                           scale7 = 1.0d0 - expdamp
     &                                 *(1.0d0-damp+0.6d0*damp**2)
                        end if
                     end if
                     dsc3 = scale3
                     dsc5 = scale5
                     dsc7 = scale7
                     psc3 = scale3
                     psc5 = scale5
                     psc7 = scale7
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           dsc3 = scale3 * dscale(kk)
                           dsc5 = scale5 * dscale(kk)
                           dsc7 = scale7 * dscale(kk)
                           psc3 = scale3 * pscale(kk)
                           psc5 = scale5 * pscale(kk)
                           psc7 = scale7 * pscale(kk)
                        end if
                     end if
c
c     find the field terms for the current interaction
c
                     dir = dix*xr + diy*yr + diz*zr
                     qix = qixx*xr + qixy*yr + qixz*zr
                     qiy = qixy*xr + qiyy*yr + qiyz*zr
                     qiz = qixz*xr + qiyz*yr + qizz*zr
                     qir = qix*xr + qiy*yr + qiz*zr
                     dkr = dkx*xr + dky*yr + dkz*zr
                     qkx = qkxx*xr + qkxy*yr + qkxz*zr
                     qky = qkxy*xr + qkyy*yr + qkyz*zr
                     qkz = qkxz*xr + qkyz*yr + qkzz*zr
                     qkr = qkx*xr + qky*yr + qkz*zr
                     bcn(1) = bn(1) - (1.0d0-dsc3)*rr3
                     bcn(2) = bn(2) - 3.0d0*(1.0d0-dsc5)*rr5
                     bcn(3) = bn(3) - 15.0d0*(1.0d0-dsc7)*rr7
                     fimd(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                           - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
                     fimd(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                           - bcn(1)*dky + 2.0d0*bcn(2)*qky
                     fimd(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                           - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
                     fkmd(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                           - bcn(1)*dix - 2.0d0*bcn(2)*qix
                     fkmd(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                           - bcn(1)*diy - 2.0d0*bcn(2)*qiy
                     fkmd(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                           - bcn(1)*diz - 2.0d0*bcn(2)*qiz
                     bcn(1) = bn(1) - (1.0d0-psc3)*rr3
                     bcn(2) = bn(2) - 3.0d0*(1.0d0-psc5)*rr5
                     bcn(3) = bn(3) - 15.0d0*(1.0d0-psc7)*rr7
                     fimp(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                           - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
                     fimp(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                           - bcn(1)*dky + 2.0d0*bcn(2)*qky
                     fimp(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                           - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
                     fkmp(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                           - bcn(1)*dix - 2.0d0*bcn(2)*qix
                     fkmp(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                           - bcn(1)*diy - 2.0d0*bcn(2)*qiy
                     fkmp(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                           - bcn(1)*diz - 2.0d0*bcn(2)*qiz
c
c     increment the field at each site due to this interaction
c
                     do j = 1, 3
                        field(j,i) = field(j,i) + fimd(j)
                        fieldp(j,i) = fieldp(j,i) + fimd(j)
                        if (ii .ne. kk) then
                           field(j,k) = field(j,k) + fkmp(j)
                           fieldp(j,k) = fieldp(j,k) + fkmp(j)
                        end if
                     end do
                  end if
               end do
            end do
c
c     reset interaction scaling coefficients for connected atoms
c
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = 1.0d0
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = 1.0d0
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = 1.0d0
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = 1.0d0
            end do
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = 1.0d0
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = 1.0d0
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = 1.0d0
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine udirect2b  --  Ewald real direct field via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "udirect2b" computes the real space contribution of the permanent
c     atomic multipole moments to the field via a neighbor list
c
c
      subroutine udirect2b (field,fieldp)
      use sizes
      use atoms
      use boxes
      use bound
      use couple
      use ewald
      use math
      use mpole
      use neigh
      use openmp
      use polar
      use polgrp
      use polpot
      use shunt
      use tarray
      use units
      implicit none
      integer i,j,k,m
      integer ii,kk,kkk
      integer nlocal,maxlocal
      integer tid,toffset0
!$    integer omp_get_thread_num
      integer, allocatable :: toffset(:)
      integer, allocatable :: ilocal(:,:)
      real*8 xr,yr,zr,r,r2
      real*8 rr1,rr2,rr3
      real*8 rr5,rr7
      real*8 erfc,bfac,exp2a
      real*8 ci,dix,diy,diz
      real*8 qixx,qiyy,qizz
      real*8 qixy,qixz,qiyz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkyy,qkzz
      real*8 qkxy,qkxz,qkyz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 ralpha,aefac
      real*8 aesq2,aesq2n
      real*8 pdi,pti,pgamma
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 bn(0:3),bcn(3)
      real*8 fimd(3),fkmd(3)
      real*8 fimp(3),fkmp(3)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8, allocatable :: fieldt(:,:)
      real*8, allocatable :: fieldtp(:,:)
      real*8, allocatable :: dlocal(:,:)
      character*6 mode
      external erfc
c
c
c     check for multipoles and set cutoff coefficients
c
      if (npole .eq. 0)  return
      mode = 'EWALD'
      call switch (mode)
      aesq2 = 2.0 * aewald * aewald
      aesq2n = 0.0d0
      if (aewald .gt. 0.0d0)  aesq2n = 1.0d0 / (sqrtpi*aewald)
      nlocal = 0
      toffset0 = 0
      maxlocal = int(dble(npole)*dble(maxelst)/dble(nthread))
c
c     perform dynamic allocation of some local arrays
c
      allocate (toffset(0:nthread-1))
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
      allocate (fieldt(3,npole))
      allocate (fieldtp(3,npole))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
      end do
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(n,npole,ipole,x,y,z,pdamp,thole,
!$OMP& rpole,p2scale,p3scale,p4scale,p41scale,p5scale,d1scale,d2scale,
!$OMP& d3scale,d4scale,u1scale,u2scale,u3scale,u4scale,n12,i12,n13,i13,
!$OMP& n14,i14,n15,i15,np11,ip11,np12,ip12,np13,ip13,np14,ip14,nelst,
!$OMP& elst,cut2,aewald,aesq2,aesq2n,poltyp,ntpair,tindex,tdipdip,
!$OMP& toffset,toffset0,field,fieldp,fieldt,fieldtp,maxlocal)
!$OMP& firstprivate(pscale,dscale,uscale,nlocal)
c
c     perform dynamic allocation of some local arrays
c
c OPT IMPLEMENTATION
c      if (poltyp .eq. 'MUTUAL') then
       if (poltyp .eq. 'MUTUAL' .or. poltyp(1:3) .eq. 'OPT') then
c OPT IMPLEMENTATION
         allocate (ilocal(2,maxlocal))
         allocate (dlocal(6,maxlocal))
      end if
c
c     initialize local variables for OpenMP calculation
c
!$OMP DO collapse(2)
      do i = 1, npole
         do j = 1, 3
            fieldt(j,i) = 0.0d0
            fieldtp(j,i) = 0.0d0
         end do
      end do
!$OMP END DO
c
c     compute the real space portion of the Ewald summation
c
!$OMP DO reduction(+:fieldt,fieldtp) schedule(guided)
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
         end do
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. cut2) then
               r = sqrt(r2)
               rr1 = 1.0d0 / r
               rr2 = rr1 * rr1
               rr3 = rr2 * rr1
               rr5 = rr2 * rr3
               rr7 = rr2 * rr5
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
c
c     calculate the error function damping factors
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) * rr1
               exp2a = exp(-ralpha**2)
               aefac = aesq2n
               do j = 1, 3
                  bfac = dble(j+j-1)
                  aefac = aesq2 * aefac
                  bn(j) = (bfac*bn(j-1)+aefac*exp2a) * rr2
               end do
c
c     compute the polarization damping scale factors
c
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                     scale7 = 1.0d0 - expdamp
     &                           *(1.0d0-damp+0.6d0*damp**2)
                  end if
               end if
c
c     find the field terms for the current interaction
c
               dir = dix*xr + diy*yr + diz*zr
               qix = qixx*xr + qixy*yr + qixz*zr
               qiy = qixy*xr + qiyy*yr + qiyz*zr
               qiz = qixz*xr + qiyz*yr + qizz*zr
               qir = qix*xr + qiy*yr + qiz*zr
               dkr = dkx*xr + dky*yr + dkz*zr
               qkx = qkxx*xr + qkxy*yr + qkxz*zr
               qky = qkxy*xr + qkyy*yr + qkyz*zr
               qkz = qkxz*xr + qkyz*yr + qkzz*zr
               qkr = qkx*xr + qky*yr + qkz*zr
               bcn(1) = bn(1) - (1.0d0-scale3*dscale(kk))*rr3
               bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*dscale(kk))*rr5
               bcn(3) = bn(3) - 15.0d0*(1.0d0-scale7*dscale(kk))*rr7
               fimd(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
               fimd(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dky + 2.0d0*bcn(2)*qky
               fimd(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
               fkmd(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*dix - 2.0d0*bcn(2)*qix
               fkmd(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diy - 2.0d0*bcn(2)*qiy
               fkmd(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diz - 2.0d0*bcn(2)*qiz
               bcn(1) = bn(1) - (1.0d0-scale3*pscale(kk))*rr3
               bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*pscale(kk))*rr5
               bcn(3) = bn(3) - 15.0d0*(1.0d0-scale7*pscale(kk))*rr7
               fimp(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
               fimp(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dky + 2.0d0*bcn(2)*qky
               fimp(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
               fkmp(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*dix - 2.0d0*bcn(2)*qix
               fkmp(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diy - 2.0d0*bcn(2)*qiy
               fkmp(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diz - 2.0d0*bcn(2)*qiz
c
c     find terms needed later to compute mutual polarization
c
c OPT IMPLEMENTATION
c              if (poltyp .eq. 'MUTUAL') then
               if (poltyp .eq. 'MUTUAL' .or. poltyp(1:3).eq.'OPT') then
c OPT IMPLEMENTATION
                  bcn(1) = bn(1) - (1.0d0-scale3*uscale(kk))*rr3
                  bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*uscale(kk))*rr5
                  nlocal = nlocal + 1
                  ilocal(1,nlocal) = i
                  ilocal(2,nlocal) = k
                  dlocal(1,nlocal) = -bcn(1) + bcn(2)*xr*xr
                  dlocal(2,nlocal) = bcn(2)*xr*yr
                  dlocal(3,nlocal) = bcn(2)*xr*zr
                  dlocal(4,nlocal) = -bcn(1) + bcn(2)*yr*yr
                  dlocal(5,nlocal) = bcn(2)*yr*zr
                  dlocal(6,nlocal) = -bcn(1) + bcn(2)*zr*zr
               end if
c
c     increment the field at each site due to this interaction
c
               do j = 1, 3
                  fieldt(j,i) = fieldt(j,i) + fimd(j)
                  fieldt(j,k) = fieldt(j,k) + fkmd(j)
                  fieldtp(j,i) = fieldtp(j,i) + fimp(j)
                  fieldtp(j,k) = fieldtp(j,k) + fkmp(j)
               end do
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            uscale(ip11(j,ii)) = 1.0d0
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            uscale(ip12(j,ii)) = 1.0d0
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            uscale(ip13(j,ii)) = 1.0d0
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            uscale(ip14(j,ii)) = 1.0d0
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
!$OMP END DO
c
c     transfer the results from local to global arrays
c
!$OMP DO
      do i = 1, npole
         do j = 1, 3
            field(j,i) = fieldt(j,i) + field(j,i)
            fieldp(j,i) = fieldtp(j,i) + fieldp(j,i)
         end do
      end do
!$OMP END DO
c
c     store terms needed later to compute mutual polarization
c
!$OMP CRITICAL
      tid = 0
!$    tid = omp_get_thread_num ()
      toffset(tid) = toffset0
      toffset0 = toffset0 + nlocal
      ntpair = toffset0
!$OMP END CRITICAL
c OPT IMPLEMENTATION
c     if (poltyp .eq. 'MUTUAL') then
      if (poltyp .eq. 'MUTUAL' .or. poltyp(1:3).eq.'OPT') then
c OPT IMPLEMENTATION
         k = toffset(tid)
         do i = 1, nlocal
            m = k + i
            tindex(1,m) = ilocal(1,i)
            tindex(2,m) = ilocal(2,i)
            do j = 1, 6
               tdipdip(j,m) = dlocal(j,i)
            end do
         end do
         deallocate (ilocal)
         deallocate (dlocal)
      end if
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (toffset)
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      deallocate (fieldt)
      deallocate (fieldtp)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine umutual1  --  Ewald recip mutual induced field  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "umutual1" computes the reciprocal space contribution of the
c     induced atomic dipole moments to the field
c
c
      subroutine umutual1 (field,fieldp)
      use sizes
      use boxes
      use ewald
      use math
      use mpole
      use pme
      use polar
      implicit none
      integer i,j,k
      real*8 term
      real*8 a(3,3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8, allocatable :: fuind(:,:)
      real*8, allocatable :: fuinp(:,:)
      real*8, allocatable :: fdip_phi1(:,:)
      real*8, allocatable :: fdip_phi2(:,:)
      real*8, allocatable :: fdip_sum_phi(:,:)
      real*8, allocatable :: dipfield1(:,:)
      real*8, allocatable :: dipfield2(:,:)
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (fuind(3,npole))
      allocate (fuinp(3,npole))
c OPT IMPLEMENTATION
      allocate (fdip_phi1(10,npole))
      allocate (fdip_phi2(10,npole))
c OPT IMPLEMENTATION
      allocate (fdip_sum_phi(20,npole))
      allocate (dipfield1(3,npole))
      allocate (dipfield2(3,npole))
c
c     convert Cartesian dipoles to fractional coordinates
c
      do i = 1, 3
         a(1,i) = dble(nfft1) * recip(i,1)
         a(2,i) = dble(nfft2) * recip(i,2)
         a(3,i) = dble(nfft3) * recip(i,3)
      end do
      do i = 1, npole
         do k = 1, 3
            fuind(k,i) = a(k,1)*uind(1,i) + a(k,2)*uind(2,i)
     &                      + a(k,3)*uind(3,i)
            fuinp(k,i) = a(k,1)*uinp(1,i) + a(k,2)*uinp(2,i)
     &                      + a(k,3)*uinp(3,i)
         end do
      end do
c
c     assign PME grid and perform 3-D FFT forward transform
c
c OPT IMPLEMENTATION
      if (ptpointer .ge. 0) then
c       Save the p and d grids separately.  This could use real->complex FFT
        dipfield1 = 0d0
        call grid_uind (fuind,dipfield1)
        call fftfront
        uindgridf(:,:,:,:,ptpointer) = qgrid
        call grid_uind (dipfield1,fuinp)
        call fftfront
        do k = 1, nfft3
           do j = 1, nfft2
              do i = 1, nfft1
                 uinpgridf(1,i,j,k,ptpointer) = qgrid(2,i,j,k)
                 uinpgridf(2,i,j,k,ptpointer) = -qgrid(1,i,j,k)
              end do
           end do
        end do
        qgrid = qgrid + uindgridf(:,:,:,:,ptpointer)
      else
        call grid_uind (fuind,fuinp)
        call fftfront
      endif
c OPT IMPLEMENTATION
c
c     complete the transformation of the PME grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               term = qfac(i,j,k)
               qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
               qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
            end do
         end do
      end do
c
c     perform 3-D FFT backward transform and get field
c
      call fftback
      call fphi_uind (fdip_phi1,fdip_phi2,fdip_sum_phi)
c OPT IMPLEMENTATION
      if (ptpointer .ge. 0) then
        ptfphid(:,:,ptpointer) = fdip_phi1
        ptfphip(:,:,ptpointer) = fdip_phi2
      endif
c OPT IMPLEMENTATION
c
c     convert the dipole fields from fractional to Cartesian
c
      do i = 1, 3
         a(i,1) = dble(nfft1) * recip(i,1)
         a(i,2) = dble(nfft2) * recip(i,2)
         a(i,3) = dble(nfft3) * recip(i,3)
      end do
      do i = 1, npole
         do k = 1, 3
            dipfield1(k,i) = a(k,1)*fdip_phi1(2,i)
     &                          + a(k,2)*fdip_phi1(3,i)
     &                          + a(k,3)*fdip_phi1(4,i)
            dipfield2(k,i) = a(k,1)*fdip_phi2(2,i)
     &                          + a(k,2)*fdip_phi2(3,i)
     &                          + a(k,3)*fdip_phi2(4,i)
         end do
      end do
c
c     increment the field at each multipole site
c
      do i = 1, npole
         do k = 1, 3
            field(k,i) = field(k,i) - dipfield1(k,i)
            fieldp(k,i) = fieldp(k,i) - dipfield2(k,i)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (fuind)
      deallocate (fuinp)
      deallocate (fdip_phi1)
      deallocate (fdip_phi2)
      deallocate (fdip_sum_phi)
      deallocate (dipfield1)
      deallocate (dipfield2)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine umutual2a  --  Ewald real mutual field via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "umutual2a" computes the real space contribution of the induced
c     atomic dipole moments to the field via a double loop
c
c
      subroutine umutual2a (field,fieldp)
      use sizes
      use atoms
      use boxes
      use bound
      use cell
      use couple
      use ewald
      use math
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      use units
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xr,yr,zr,r,r2
      real*8 rr1,rr2,rr3,rr5
      real*8 erfc,bfac,exp2a
      real*8 duir,dukr
      real*8 puir,pukr
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 ralpha,aefac
      real*8 aesq2,aesq2n
      real*8 pdi,pti,pgamma
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 bn(0:2)
      real*8 fimd(3),fkmd(3)
      real*8 fimp(3),fkmp(3)
      real*8, allocatable :: uscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      character*6 mode
      external erfc
c
c
c     check for multipoles and set cutoff coefficients
c
      if (npole .eq. 0)  return
      mode = 'EWALD'
      call switch (mode)
      aesq2 = 2.0 * aewald * aewald
      aesq2n = 0.0d0
      if (aewald .gt. 0.0d0)  aesq2n = 1.0d0 / (sqrtpi*aewald)
c
c     perform dynamic allocation of some local arrays
c
      allocate (uscale(n))
c
c     set array needed to scale connected atom interactions
c
      do i = 1, n
         uscale(i) = 1.0d0
      end do
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         duix = uind(1,i)
         duiy = uind(2,i)
         duiz = uind(3,i)
         puix = uinp(1,i)
         puiy = uinp(2,i)
         puiz = uinp(3,i)
         do j = 1, np11(ii)
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            uscale(ip14(j,ii)) = u4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. cut2) then
               r = sqrt(r2)
               rr1 = 1.0d0 / r
               rr2 = rr1 * rr1
               rr3 = rr2 * rr1
               rr5 = rr2 * rr3
               dukx = uind(1,k)
               duky = uind(2,k)
               dukz = uind(3,k)
               pukx = uinp(1,k)
               puky = uinp(2,k)
               pukz = uinp(3,k)
c
c     calculate the error function damping factors
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) * rr1
               exp2a = exp(-ralpha**2)
               aefac = aesq2n
               do j = 1, 2
                  bfac = dble(j+j-1)
                  aefac = aesq2 * aefac
                  bn(j) = (bfac*bn(j-1)+aefac*exp2a) * rr2
               end do
c
c     compute the polarization damping scale factors
c
               scale3 = uscale(kk)
               scale5 = uscale(kk)
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = scale3 * (1.0d0-expdamp)
                     scale5 = scale5 * (1.0d0-(1.0d0-damp)*expdamp)
                  end if
               end if
               bn(1) = -bn(1) + (1.0d0-scale3)*rr3
               bn(2) = bn(2) - 3.0d0*(1.0d0-scale5)*rr5
c
c     find the field terms for the current interaction
c
               duir = xr*duix + yr*duiy + zr*duiz
               dukr = xr*dukx + yr*duky + zr*dukz
               puir = xr*puix + yr*puiy + zr*puiz
               pukr = xr*pukx + yr*puky + zr*pukz
               fimd(1) = bn(1)*dukx + bn(2)*dukr*xr
               fimd(2) = bn(1)*duky + bn(2)*dukr*yr
               fimd(3) = bn(1)*dukz + bn(2)*dukr*zr
               fkmd(1) = bn(1)*duix + bn(2)*duir*xr
               fkmd(2) = bn(1)*duiy + bn(2)*duir*yr
               fkmd(3) = bn(1)*duiz + bn(2)*duir*zr
               fimp(1) = bn(1)*pukx + bn(2)*pukr*xr
               fimp(2) = bn(1)*puky + bn(2)*pukr*yr
               fimp(3) = bn(1)*pukz + bn(2)*pukr*zr
               fkmp(1) = bn(1)*puix + bn(2)*puir*xr
               fkmp(2) = bn(1)*puiy + bn(2)*puir*yr
               fkmp(3) = bn(1)*puiz + bn(2)*puir*zr
c
c     increment the field at each site due to this interaction
c
               do j = 1, 3
                  field(j,i) = field(j,i) + fimd(j)
                  field(j,k) = field(j,k) + fkmd(j)
                  fieldp(j,i) = fieldp(j,i) + fimp(j)
                  fieldp(j,k) = fieldp(j,k) + fkmp(j)
               end do
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, np11(ii)
            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            uscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
         do i = 1, npole
            ii = ipole(i)
            pdi = pdamp(i)
            pti = thole(i)
            duix = uind(1,i)
            duiy = uind(2,i)
            duiz = uind(3,i)
            puix = uinp(1,i)
            puiy = uinp(2,i)
            puiz = uinp(3,i)
            do j = 1, np11(ii)
               uscale(ip11(j,ii)) = u1scale
            end do
            do j = 1, np12(ii)
               uscale(ip12(j,ii)) = u2scale
            end do
            do j = 1, np13(ii)
               uscale(ip13(j,ii)) = u3scale
            end do
            do j = 1, np14(ii)
               uscale(ip14(j,ii)) = u4scale
            end do
            do k = i, npole
               kk = ipole(k)
               dukx = uind(1,k)
               duky = uind(2,k)
               dukz = uind(3,k)
               pukx = uinp(1,k)
               puky = uinp(2,k)
               pukz = uinp(3,k)
               do m = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call imager (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
c
c     calculate the error function damping factors
c
                  if (r2 .le. cut2) then
                     r = sqrt(r2)
                     rr1 = 1.0d0 / r
                     rr2 = rr1 * rr1
                     rr3 = rr2 * rr1
                     rr5 = rr2 * rr3
                     ralpha = aewald * r
                     bn(0) = erfc(ralpha) * rr1
                     exp2a = exp(-ralpha**2)
                     aefac = aesq2n
                     do j = 1, 2
                        bfac = dble(j+j-1)
                        aefac = aesq2 * aefac
                        bn(j) = (bfac*bn(j-1)+aefac*exp2a) * rr2
                     end do
c
c     compute the polarization damping scale factors
c
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(k))
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           expdamp = exp(damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                        end if
                     end if
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           scale3 = scale3 * uscale(kk)
                           scale5 = scale5 * uscale(kk)
                        end if
                     end if
                     bn(1) = -bn(1) + (1.0d0-scale3)*rr3
                     bn(2) = bn(2) - 3.0d0*(1.0d0-scale5)*rr5
c
c     find the field terms for the current interaction
c
                     duir = xr*duix + yr*duiy + zr*duiz
                     dukr = xr*dukx + yr*duky + zr*dukz
                     puir = xr*puix + yr*puiy + zr*puiz
                     pukr = xr*pukx + yr*puky + zr*pukz
                     fimd(1) = bn(1)*dukx + bn(2)*dukr*xr
                     fimd(2) = bn(1)*duky + bn(2)*dukr*yr
                     fimd(3) = bn(1)*dukz + bn(2)*dukr*zr
                     fkmd(1) = bn(1)*duix + bn(2)*duir*xr
                     fkmd(2) = bn(1)*duiy + bn(2)*duir*yr
                     fkmd(3) = bn(1)*duiz + bn(2)*duir*zr
                     fimp(1) = bn(1)*pukx + bn(2)*pukr*xr
                     fimp(2) = bn(1)*puky + bn(2)*pukr*yr
                     fimp(3) = bn(1)*pukz + bn(2)*pukr*zr
                     fkmp(1) = bn(1)*puix + bn(2)*puir*xr
                     fkmp(2) = bn(1)*puiy + bn(2)*puir*yr
                     fkmp(3) = bn(1)*puiz + bn(2)*puir*zr
c
c     increment the field at each site due to this interaction
c
                     do j = 1, 3
                        field(j,i) = field(j,i) + fimd(j)
                        fieldp(j,i) = fieldp(j,i) + fimp(j)
                        if (ii .ne. kk) then
                           field(j,k) = field(j,k) + fkmd(j)
                           fieldp(j,k) = fieldp(j,k) + fkmp(j)
                        end if
                     end do
                  end if
               end do
            end do
c
c     reset interaction scaling coefficients for connected atoms
c
            do j = 1, np11(ii)
               uscale(ip11(j,ii)) = 1.0d0
            end do
            do j = 1, np12(ii)
               uscale(ip12(j,ii)) = 1.0d0
            end do
            do j = 1, np13(ii)
               uscale(ip13(j,ii)) = 1.0d0
            end do
            do j = 1, np14(ii)
               uscale(ip14(j,ii)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (uscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine umutual2b  --  Ewald real mutual field via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "umutual2b" computes the real space contribution of the induced
c     atomic dipole moments to the field via a neighbor list
c
c
      subroutine umutual2b (field,fieldp)
      use sizes
      use mpole
      use polar
      use tarray
      implicit none
      integer i,j,k,m
      real*8 fimd(3),fkmd(3)
      real*8 fimp(3),fkmp(3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8, allocatable :: fieldt(:,:)
      real*8, allocatable :: fieldtp(:,:)
c
c
c     check for multipoles and set cutoff coefficients
c
      if (npole .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (fieldt(3,npole))
      allocate (fieldtp(3,npole))
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,uind,uinp,ntpair,tindex,
!$OMP& tdipdip,field,fieldp,fieldt,fieldtp)
c
c     initialize local variables for OpenMP calculation
c
!$OMP DO collapse(2)
      do i = 1, npole
         do j = 1, 3
            fieldt(j,i) = 0.0d0
            fieldtp(j,i) = 0.0d0
         end do
      end do
!$OMP END DO
c
c     find the field terms for each pairwise interaction
c
!$OMP DO reduction(+:fieldt,fieldtp) schedule(guided)
      do m = 1, ntpair
         i = tindex(1,m)
         k = tindex(2,m)
         fimd(1) = tdipdip(1,m)*uind(1,k) + tdipdip(2,m)*uind(2,k)
     &                + tdipdip(3,m)*uind(3,k)
         fimd(2) = tdipdip(2,m)*uind(1,k) + tdipdip(4,m)*uind(2,k)
     &                + tdipdip(5,m)*uind(3,k)
         fimd(3) = tdipdip(3,m)*uind(1,k) + tdipdip(5,m)*uind(2,k)
     &                + tdipdip(6,m)*uind(3,k)
         fkmd(1) = tdipdip(1,m)*uind(1,i) + tdipdip(2,m)*uind(2,i)
     &                + tdipdip(3,m)*uind(3,i)
         fkmd(2) = tdipdip(2,m)*uind(1,i) + tdipdip(4,m)*uind(2,i)
     &                + tdipdip(5,m)*uind(3,i)
         fkmd(3) = tdipdip(3,m)*uind(1,i) + tdipdip(5,m)*uind(2,i)
     &                + tdipdip(6,m)*uind(3,i)
         fimp(1) = tdipdip(1,m)*uinp(1,k) + tdipdip(2,m)*uinp(2,k)
     &                + tdipdip(3,m)*uinp(3,k)
         fimp(2) = tdipdip(2,m)*uinp(1,k) + tdipdip(4,m)*uinp(2,k)
     &                + tdipdip(5,m)*uinp(3,k)
         fimp(3) = tdipdip(3,m)*uinp(1,k) + tdipdip(5,m)*uinp(2,k)
     &                + tdipdip(6,m)*uinp(3,k)
         fkmp(1) = tdipdip(1,m)*uinp(1,i) + tdipdip(2,m)*uinp(2,i)
     &                + tdipdip(3,m)*uinp(3,i)
         fkmp(2) = tdipdip(2,m)*uinp(1,i) + tdipdip(4,m)*uinp(2,i)
     &                + tdipdip(5,m)*uinp(3,i)
         fkmp(3) = tdipdip(3,m)*uinp(1,i) + tdipdip(5,m)*uinp(2,i)
     &                + tdipdip(6,m)*uinp(3,i)
c
c     increment the field at each site due to this interaction
c
         do j = 1, 3
            fieldt(j,i) = fieldt(j,i) + fimd(j)
            fieldt(j,k) = fieldt(j,k) + fkmd(j)
            fieldtp(j,i) = fieldtp(j,i) + fimp(j)
            fieldtp(j,k) = fieldtp(j,k) + fkmp(j)
         end do
      end do
!$OMP END DO
c
c     end OpenMP directives for the major loop structure
c
!$OMP DO
      do i = 1, npole
         do j = 1, 3
            field(j,i) = fieldt(j,i) + field(j,i)
            fieldp(j,i) = fieldtp(j,i) + fieldp(j,i)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (fieldt)
      deallocate (fieldtp)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine induce0d  --  Kirkwood SCRF induced dipoles  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "induce0d" computes the induced dipole moments at polarizable
c     sites for generalized Kirkwood SCRF and vacuum environments
c
c
      subroutine induce0d
      use sizes
      use atoms
      use inform
      use iounit
      use mpole
      use polar
      use polpot
      use potent
      use units
      use uprior
      implicit none
      integer i,j,k,iter
      integer maxiter
      real*8 polmin
      real*8 eps,epsold
      real*8 epsd,epsp
      real*8 epsds,epsps
      real*8 udsum,upsum
      real*8 ussum,upssum
      real*8 a,ap,as,aps
      real*8 b,bp,bs,bps
      real*8 sum,sump
      real*8 sums,sumps
      real*8, allocatable :: poli(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8, allocatable :: fields(:,:)
      real*8, allocatable :: fieldps(:,:)
      real*8, allocatable :: udir(:,:)
      real*8, allocatable :: udirp(:,:)
      real*8, allocatable :: udirs(:,:)
      real*8, allocatable :: udirps(:,:)
      real*8, allocatable :: rsd(:,:)
      real*8, allocatable :: rsdp(:,:)
      real*8, allocatable :: rsds(:,:)
      real*8, allocatable :: rsdps(:,:)
      real*8, allocatable :: zrsd(:,:)
      real*8, allocatable :: zrsdp(:,:)
      real*8, allocatable :: zrsds(:,:)
      real*8, allocatable :: zrsdps(:,:)
      real*8, allocatable :: conj(:,:)
      real*8, allocatable :: conjp(:,:)
      real*8, allocatable :: conjs(:,:)
      real*8, allocatable :: conjps(:,:)
      real*8, allocatable :: vec(:,:)
      real*8, allocatable :: vecp(:,:)
      real*8, allocatable :: vecs(:,:)
      real*8, allocatable :: vecps(:,:)
      logical done
      character*6 mode
c
c
c     zero out the induced dipoles at each site; uind and uinp are
c     vacuum dipoles, uinds and uinps are SCRF dipoles
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
            uinds(j,i) = 0.0d0
            uinps(j,i) = 0.0d0
         end do
      end do
      if (.not.use_polar .and. .not.use_solv)  return
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (field(3,npole))
      allocate (fieldp(3,npole))
      allocate (fields(3,npole))
      allocate (fieldps(3,npole))
      allocate (udir(3,npole))
      allocate (udirp(3,npole))
      allocate (udirs(3,npole))
      allocate (udirps(3,npole))
c
c     compute the direct induced dipole moment at each atom, and
c     another set that also includes RF due to permanent multipoles
c
      call dfield0d (field,fieldp,fields,fieldps)
c
c     set vacuum induced dipoles to polarizability times direct field;
c     set SCRF induced dipoles to polarizability times direct field
c     plus the GK reaction field due to permanent multipoles
c
      do i = 1, npole
         do j = 1, 3
            udir(j,i) = polarity(i) * field(j,i)
            udirp(j,i) = polarity(i) * fieldp(j,i)
            udirs(j,i) = polarity(i) * fields(j,i)
            udirps(j,i) = polarity(i) * fieldps(j,i)
            uind(j,i) = udir(j,i)
            uinp(j,i) = udirp(j,i)
            uinds(j,i) = udirs(j,i)
            uinps(j,i) = udirps(j,i)
         end do
      end do
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 500
         iter = 0
         polmin = 0.00000001d0
         eps = 100.0d0
c
c     estimated induced dipoles from polynomial predictor
c
         if (use_pred .and. nualt.eq.maxualt) then
            do i = 1, npole
               do j = 1, 3
                  udsum = 0.0d0
                  upsum = 0.0d0
                  ussum = 0.0d0
                  upssum = 0.0d0
                  do k = 1, nualt-1
                     udsum = udsum + bpred(k)*udalt(k,j,i)
                     upsum = upsum + bpredp(k)*upalt(k,j,i)
                     ussum = ussum + bpreds(k)*usalt(k,j,i)
                     upssum = upssum + bpredps(k)*upsalt(k,j,i)
                  end do
                  uind(j,i) = udsum
                  uinp(j,i) = upsum
                  uinds(j,i) = ussum
                  uinps(j,i) = upssum
               end do
            end do
         end if
c
c     perform dynamic allocation of some local arrays
c
         allocate (poli(npole))
         allocate (rsd(3,npole))
         allocate (rsdp(3,npole))
         allocate (rsds(3,npole))
         allocate (rsdps(3,npole))
         allocate (zrsd(3,npole))
         allocate (zrsdp(3,npole))
         allocate (zrsds(3,npole))
         allocate (zrsdps(3,npole))
         allocate (conj(3,npole))
         allocate (conjp(3,npole))
         allocate (conjs(3,npole))
         allocate (conjps(3,npole))
         allocate (vec(3,npole))
         allocate (vecp(3,npole))
         allocate (vecs(3,npole))
         allocate (vecps(3,npole))
c
c     set initial conjugate gradient residual and conjugate vector
c
         call ufield0d (field,fieldp,fields,fieldps)
         do i = 1, npole
            poli(i) = max(polmin,polarity(i))
            do j = 1, 3
               rsd(j,i) = (udir(j,i)-uind(j,i))/poli(i)
     &                       + field(j,i)
               rsdp(j,i) = (udirp(j,i)-uinp(j,i))/poli(i)
     &                        + fieldp(j,i)
               rsds(j,i) = (udirs(j,i)-uinds(j,i))/poli(i)
     &                        + fields(j,i)
               rsdps(j,i) = (udirps(j,i)-uinps(j,i))/poli(i)
     &                         + fieldps(j,i)
               zrsd(j,i) = rsd(j,i) * poli(i)
               zrsdp(j,i) = rsdp(j,i) * poli(i)
               zrsds(j,i) = rsds(j,i) * poli(i)
               zrsdps(j,i) = rsdps(j,i) * poli(i)
               conj(j,i) = zrsd(j,i)
               conjp(j,i) = zrsdp(j,i)
               conjs(j,i) = zrsds(j,i)
               conjps(j,i) = zrsdps(j,i)
            end do
         end do
c
c     conjugate gradient iteration of the mutual induced dipoles
c
         do while (.not. done)
            iter = iter + 1
            do i = 1, npole
               do j = 1, 3
                  vec(j,i) = uind(j,i)
                  vecp(j,i) = uinp(j,i)
                  vecs(j,i) = uinds(j,i)
                  vecps(j,i) = uinps(j,i)
                  uind(j,i) = conj(j,i)
                  uinp(j,i) = conjp(j,i)
                  uinds(j,i) = conjs(j,i)
                  uinps(j,i) = conjps(j,i)
               end do
            end do
            call ufield0d (field,fieldp,fields,fieldps)
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = vec(j,i)
                  uinp(j,i) = vecp(j,i)
                  uinds(j,i) = vecs(j,i)
                  uinps(j,i) = vecps(j,i)
                  vec(j,i) = conj(j,i)/poli(i) - field(j,i)
                  vecp(j,i) = conjp(j,i)/poli(i) - fieldp(j,i)
                  vecs(j,i) = conjs(j,i)/poli(i) - fields(j,i)
                  vecps(j,i) = conjps(j,i)/poli(i) - fieldps(j,i)
               end do
            end do
            a = 0.0d0
            ap = 0.0d0
            as = 0.0d0
            aps = 0.0d0
            sum = 0.0d0
            sump = 0.0d0
            sums = 0.0d0
            sumps = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  a = a + conj(j,i)*vec(j,i)
                  ap = ap + conjp(j,i)*vecp(j,i)
                  as = as + conjs(j,i)*vecs(j,i)
                  aps = aps + conjps(j,i)*vecps(j,i)
                  sum = sum + rsd(j,i)*zrsd(j,i)
                  sump = sump + rsdp(j,i)*zrsdp(j,i)
                  sums = sums + rsds(j,i)*zrsds(j,i)
                  sumps = sumps + rsdps(j,i)*zrsdps(j,i)
               end do
            end do
            if (a .ne. 0.0d0)  a = sum / a
            if (ap .ne. 0.0d0)  ap = sump / ap
            if (as .ne. 0.0d0)  as = sums / as
            if (aps .ne. 0.0d0)  aps = sumps / aps
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = uind(j,i) + a*conj(j,i)
                  uinp(j,i) = uinp(j,i) + ap*conjp(j,i)
                  uinds(j,i) = uinds(j,i) + as*conjs(j,i)
                  uinps(j,i) = uinps(j,i) + aps*conjps(j,i)
                  rsd(j,i) = rsd(j,i) - a*vec(j,i)
                  rsdp(j,i) = rsdp(j,i) - ap*vecp(j,i)
                  rsds(j,i) = rsds(j,i) - as*vecs(j,i)
                  rsdps(j,i) = rsdps(j,i) - aps*vecps(j,i)
               end do
            end do
            b = 0.0d0
            bp = 0.0d0
            bs = 0.0d0
            bps = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  zrsd(j,i) = rsd(j,i) * poli(i)
                  zrsdp(j,i) = rsdp(j,i) * poli(i)
                  zrsds(j,i) = rsds(j,i) * poli(i)
                  zrsdps(j,i) = rsdps(j,i) * poli(i)
                  b = b + rsd(j,i)*zrsd(j,i)
                  bp = bp + rsdp(j,i)*zrsdp(j,i)
                  bs = bs + rsds(j,i)*zrsds(j,i)
                  bps = bps + rsdps(j,i)*zrsdps(j,i)
               end do
            end do
            if (sum .ne. 0.0d0)  b = b / sum
            if (sump .ne. 0.0d0)  bp = bp / sump
            if (sums .ne. 0.0d0)  bs = bs / sums
            if (sumps .ne. 0.0d0)  bps = bps / sumps
            epsd = 0.0d0
            epsp = 0.0d0
            epsds = 0.0d0
            epsps = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  conj(j,i) = zrsd(j,i) + b*conj(j,i)
                  conjp(j,i) = zrsdp(j,i) + bp*conjp(j,i)
                  conjs(j,i) = zrsds(j,i) + bs*conjs(j,i)
                  conjps(j,i) = zrsdps(j,i) + bps*conjps(j,i)
                  epsd = epsd + rsd(j,i)*rsd(j,i)
                  epsp = epsp + rsdp(j,i)*rsdp(j,i)
                  epsds = epsds + rsds(j,i)*rsds(j,i)
                  epsps = epsps + rsdps(j,i)*rsdps(j,i)
               end do
            end do
c
c     check the convergence of the mutual induced dipoles
c
            epsold = eps
            eps = max(epsd,epsp,epsds,epsps)
            eps = debye * sqrt(eps/dble(npolar))
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. politer)  done = .true.
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (poli)
         deallocate (rsd)
         deallocate (rsdp)
         deallocate (rsds)
         deallocate (rsdps)
         deallocate (zrsd)
         deallocate (zrsdp)
         deallocate (zrsds)
         deallocate (zrsdps)
         deallocate (conj)
         deallocate (conjp)
         deallocate (conjs)
         deallocate (conjps)
         deallocate (vec)
         deallocate (vecp)
         deallocate (vecs)
         deallocate (vecps)
c
c     print the results from the conjugate gradient iteration
c
         if (debug) then
            write (iout,30)  iter,eps
   30       format (/,' Induced Dipoles :',6x,'Iterations',i5,
     &                 6x,'RMS Change',f15.10)
         end if
c
c     terminate the calculation if dipoles failed to converge
c
         if (iter.ge.maxiter .or. eps.gt.epsold) then
            write (iout,40)
   40       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
            call prterr
            call fatal
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (field)
      deallocate (fieldp)
      deallocate (fields)
      deallocate (fieldps)
      deallocate (udir)
      deallocate (udirp)
      deallocate (udirs)
      deallocate (udirps)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine dfield0d  --  generalized Kirkwood direct field  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "dfield0d" computes the direct electrostatic field due to
c     permanent multipole moments for use with with generalized
c     Kirkwood implicit solvation
c
c
      subroutine dfield0d (field,fieldp,fields,fieldps)
      use sizes
      use atoms
      use couple
      use gkstuf
      use group
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      use solute
      implicit none
      integer i,j,k,ii,kk
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 fgrp,r,r2
      real*8 rr3,rr5,rr7
      real*8 ci,uxi,uyi,uzi
      real*8 qxxi,qxyi,qxzi
      real*8 qyyi,qyzi,qzzi
      real*8 ck,uxk,uyk,uzk
      real*8 qxxk,qxyk,qxzk
      real*8 qyyk,qyzk,qzzk
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 pdi,pti,pgamma
      real*8 rb2,rbi,rbk
      real*8 dwater,fc,fd,fq
      real*8 gf,gf2,gf3,gf5,gf7
      real*8 expterm,expc,expc1
      real*8 dexpc,expcdexpc
      real*8 a(0:3,0:2)
      real*8 gc(4),gux(10)
      real*8 guy(10),guz(10)
      real*8 gqxx(4),gqxy(4)
      real*8 gqxz(4),gqyy(4)
      real*8 gqyz(4),gqzz(4)
      real*8 fid(3),fkd(3)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8 fields(3,*)
      real*8 fieldps(3,*)
      real*8, allocatable :: fieldt(:,:)
      real*8, allocatable :: fieldtp(:,:)
      real*8, allocatable :: fieldts(:,:)
      real*8, allocatable :: fieldtps(:,:)
      logical proceed
c
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
            fields(j,i) = 0.0d0
            fieldps(j,i) = 0.0d0
         end do
      end do
c
c     set dielectric constant and scaling factors for water
c
      dwater = 78.3d0
      fc = 1.0d0 * (1.0d0-dwater) / (1.0d0*dwater)
      fd = 2.0d0 * (1.0d0-dwater) / (1.0d0+2.0d0*dwater)
      fq = 3.0d0 * (1.0d0-dwater) / (2.0d0+3.0d0*dwater)
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         dscale(i) = 1.0d0
         pscale(i) = 1.0d0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (fieldt(3,npole))
      allocate (fieldtp(3,npole))
      allocate (fieldts(3,npole))
      allocate (fieldtps(3,npole))
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,pdamp,thole,rborn,
!$OMP& rpole,n12,n13,n14,n15,np11,np12,np13,np14,i12,i13,i14,i15,
!$OMP% ip11,ip12,ip13,ip14,p2scale,p3scale,p4scale,p41scale,p5scale,
!$OMP& d1scale,d2scale,d3scale,d4scale,use_intra,x,y,z,off2,fc,fd,fq,
!$OMP& gkc,field,fieldp,fields,fieldps)
!$OMP& firstprivate(dscale,pscale)
!$OMP% shared(fieldt,fieldtp,fieldts,fieldtps)
c
c     initialize local variables for OpenMP calculation
c
!$OMP DO collapse(2)
      do i = 1, npole
         do j = 1, 3
            fieldt(j,i) = 0.0d0
            fieldtp(j,i) = 0.0d0
            fieldts(j,i) = 0.0d0
            fieldtps(j,i) = 0.0d0
         end do
      end do
!$OMP END DO
c
c     find the field terms for each pairwise interaction
c
!$OMP DO reduction(+:fieldt,fieldtp,fieldts,fieldtps)
!$OMP& schedule(guided)
c
c     compute the direct induced dipole moment at each atom, and
c     another set that also includes RF due to permanent multipoles
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         rbi = rborn(ii)
         ci = rpole(1,i)
         uxi = rpole(2,i)
         uyi = rpole(3,i)
         uzi = rpole(4,i)
         qxxi = rpole(5,i)
         qxyi = rpole(6,i)
         qxzi = rpole(7,i)
         qyyi = rpole(9,i)
         qyzi = rpole(10,i)
         qzzi = rpole(13,i)
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
         do k = i, npole
            kk = ipole(k)
            rbk = rborn(kk)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               xr2 = xr * xr
               yr2 = yr * yr
               zr2 = zr * zr
               r2 = xr2 + yr2 + zr2
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  uxk = rpole(2,k)
                  uyk = rpole(3,k)
                  uzk = rpole(4,k)
                  qxxk = rpole(5,k)
                  qxyk = rpole(6,k)
                  qxzk = rpole(7,k)
                  qyyk = rpole(9,k)
                  qyzk = rpole(10,k)
                  qzzk = rpole(13,k)
c
c     self-interactions for the solute field are skipped
c
                  if (i .ne. k) then
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     scale7 = 1.0d0
                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(k))
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           expdamp = exp(damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                           scale7 = 1.0d0 - expdamp
     &                                 *(1.0d0-damp+0.6d0*damp**2)
                        end if
                     end if
                     rr3 = scale3 / (r*r2)
                     rr5 = 3.0d0 * scale5 / (r*r2*r2)
                     rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
                     dir = uxi*xr + uyi*yr + uzi*zr
                     qix = qxxi*xr + qxyi*yr + qxzi*zr
                     qiy = qxyi*xr + qyyi*yr + qyzi*zr
                     qiz = qxzi*xr + qyzi*yr + qzzi*zr
                     qir = qix*xr + qiy*yr + qiz*zr
                     dkr = uxk*xr + uyk*yr + uzk*zr
                     qkx = qxxk*xr + qxyk*yr + qxzk*zr
                     qky = qxyk*xr + qyyk*yr + qyzk*zr
                     qkz = qxzk*xr + qyzk*yr + qzzk*zr
                     qkr = qkx*xr + qky*yr + qkz*zr
                     fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                           - rr3*uxk + 2.0d0*rr5*qkx
                     fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                           - rr3*uyk + 2.0d0*rr5*qky
                     fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                           - rr3*uzk + 2.0d0*rr5*qkz
                     fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                           - rr3*uxi - 2.0d0*rr5*qix
                     fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                           - rr3*uyi - 2.0d0*rr5*qiy
                     fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                           - rr3*uzi - 2.0d0*rr5*qiz
                     do j = 1, 3
                        fieldt(j,i) = fieldt(j,i) + fid(j)*dscale(kk)
                        fieldt(j,k) = fieldt(j,k) + fkd(j)*dscale(kk)
                        fieldtp(j,i) = fieldtp(j,i) + fid(j)*pscale(kk)
                        fieldtp(j,k) = fieldtp(j,k) + fkd(j)*pscale(kk)
                     end do
                  end if
c
c     set the reaction potential auxiliary terms
c
                  rb2 = rbi * rbk
                  expterm = exp(-r2/(gkc*rb2))
                  expc = expterm / gkc
                  dexpc = -2.0d0 / (gkc*rb2)
                  gf2 = 1.0d0 / (r2+rb2*expterm)
                  gf = sqrt(gf2)
                  gf3 = gf2 * gf
                  gf5 = gf3 * gf2
                  gf7 = gf5 * gf2
                  a(0,0) = gf
                  a(1,0) = -gf3
                  a(2,0) = 3.0d0 * gf5
                  a(3,0) = -15.0d0 * gf7
c
c     set the reaction potential gradient auxiliary terms
c
                  expc1 = 1.0d0 - expc
                  a(0,1) = expc1 * a(1,0)
                  a(1,1) = expc1 * a(2,0)
                  a(2,1) = expc1 * a(3,0)
c
c     dipole second reaction potential gradient auxiliary term
c
                  expcdexpc = -expc * dexpc
                  a(1,2) = expc1*a(2,1) + expcdexpc*a(2,0)
c
c     multiply the auxiliary terms by dielectric functions
c
                  a(0,1) = fc * a(0,1)
                  a(1,0) = fd * a(1,0)
                  a(1,1) = fd * a(1,1)
                  a(1,2) = fd * a(1,2)
                  a(2,0) = fq * a(2,0)
                  a(2,1) = fq * a(2,1)
c
c     unweighted dipole reaction potential tensor
c
                  gux(1) = xr * a(1,0)
                  guy(1) = yr * a(1,0)
                  guz(1) = zr * a(1,0)
c
c     unweighted reaction potential gradient tensor
c
                  gc(2) = xr * a(0,1)
                  gc(3) = yr * a(0,1)
                  gc(4) = zr * a(0,1)
                  gux(2) = a(1,0) + xr2*a(1,1)
                  gux(3) = xr * yr * a(1,1)
                  gux(4) = xr * zr * a(1,1)
                  guy(2) = gux(3)
                  guy(3) = a(1,0) + yr2*a(1,1)
                  guy(4) = yr * zr * a(1,1)
                  guz(2) = gux(4)
                  guz(3) = guy(4)
                  guz(4) = a(1,0) + zr2*a(1,1)
                  gqxx(2) = xr * (2.0d0*a(2,0)+xr2*a(2,1))
                  gqxx(3) = yr * xr2*a(2,1)
                  gqxx(4) = zr * xr2*a(2,1)
                  gqyy(2) = xr * yr2*a(2,1)
                  gqyy(3) = yr * (2.0d0*a(2,0)+yr2*a(2,1))
                  gqyy(4) = zr * yr2 * a(2,1)
                  gqzz(2) = xr * zr2 * a(2,1)
                  gqzz(3) = yr * zr2 * a(2,1)
                  gqzz(4) = zr * (2.0d0*a(2,0)+zr2*a(2,1))
                  gqxy(2) = yr * (a(2,0)+xr2*a(2,1))
                  gqxy(3) = xr * (a(2,0)+yr2*a(2,1))
                  gqxy(4) = zr * xr * yr * a(2,1)
                  gqxz(2) = zr * (a(2,0)+xr2*a(2,1))
                  gqxz(3) = gqxy(4)
                  gqxz(4) = xr * (a(2,0)+zr2*a(2,1))
                  gqyz(2) = gqxy(4)
                  gqyz(3) = zr * (a(2,0)+yr2*a(2,1))
                  gqyz(4) = yr * (a(2,0)+zr2*a(2,1))
c
c     unweighted dipole second reaction potential gradient tensor
c
                  gux(5) = xr * (3.0d0*a(1,1)+xr2*a(1,2))
                  gux(6) = yr * (a(1,1)+xr2*a(1,2))
                  gux(7) = zr * (a(1,1)+xr2*a(1,2))
                  gux(8) = xr * (a(1,1)+yr2*a(1,2))
                  gux(9) = zr * xr * yr * a(1,2)
                  gux(10) = xr * (a(1,1)+zr2*a(1,2))
                  guy(5) = yr * (a(1,1)+xr2*a(1,2))
                  guy(6) = xr * (a(1,1)+yr2*a(1,2))
                  guy(7) = gux(9)
                  guy(8) = yr * (3.0d0*a(1,1)+yr2*a(1,2))
                  guy(9) = zr * (a(1,1)+yr2*a(1,2))
                  guy(10) = yr * (a(1,1)+zr2*a(1,2))
                  guz(5) = zr * (a(1,1)+xr2*a(1,2))
                  guz(6) = gux(9)
                  guz(7) = xr * (a(1,1)+zr2*a(1,2))
                  guz(8) = zr * (a(1,1)+yr2*a(1,2))
                  guz(9) = yr * (a(1,1)+zr2*a(1,2))
                  guz(10) = zr * (3.0d0*a(1,1)+zr2*a(1,2))
c
c     generalized Kirkwood permanent reaction field
c
                  fid(1) = uxk*gux(2) + uyk*gux(3) + uzk*gux(4)
     &                        + 0.5d0 * (ck*gux(1) + qxxk*gux(5)
     &                            + qyyk*gux(8) + qzzk*gux(10)
     &                            + 2.0d0*(qxyk*gux(6)+qxzk*gux(7)
     &                                         +qyzk*gux(9)))
     &                        + 0.5d0 * (ck*gc(2) + qxxk*gqxx(2)
     &                            + qyyk*gqyy(2) + qzzk*gqzz(2)
     &                            + 2.0d0*(qxyk*gqxy(2)+qxzk*gqxz(2)
     &                                         +qyzk*gqyz(2)))
                  fid(2) = uxk*guy(2) + uyk*guy(3) + uzk*guy(4)
     &                        + 0.5d0 * (ck*guy(1) + qxxk*guy(5)
     &                            + qyyk*guy(8) + qzzk*guy(10)
     &                            + 2.0d0*(qxyk*guy(6)+qxzk*guy(7)
     &                                         +qyzk*guy(9)))
     &                        + 0.5d0 * (ck*gc(3) + qxxk*gqxx(3)
     &                            + qyyk*gqyy(3) + qzzk*gqzz(3)
     &                            + 2.0d0*(qxyk*gqxy(3)+qxzk*gqxz(3)
     &                                         +qyzk*gqyz(3)))
                  fid(3) = uxk*guz(2) + uyk*guz(3) + uzk*guz(4)
     &                        + 0.5d0 * (ck*guz(1) + qxxk*guz(5)
     &                            + qyyk*guz(8) + qzzk*guz(10)
     &                            + 2.0d0*(qxyk*guz(6)+qxzk*guz(7)
     &                                         +qyzk*guz(9)))
     &                        + 0.5d0 * (ck*gc(4) + qxxk*gqxx(4)
     &                            + qyyk*gqyy(4) + qzzk*gqzz(4)
     &                            + 2.0d0*(qxyk*gqxy(4)+qxzk*gqxz(4)
     &                                         +qyzk*gqyz(4)))
                  fkd(1) = uxi*gux(2) + uyi*gux(3) + uzi*gux(4)
     &                        - 0.5d0 * (ci*gux(1) + qxxi*gux(5)
     &                            + qyyi*gux(8) + qzzi*gux(10)
     &                            + 2.0d0*(qxyi*gux(6)+qxzi*gux(7)
     &                                         +qyzi*gux(9)))
     &                        - 0.5d0 * (ci*gc(2) + qxxi*gqxx(2)
     &                            + qyyi*gqyy(2) + qzzi*gqzz(2)
     &                            + 2.0d0*(qxyi*gqxy(2)+qxzi*gqxz(2)
     &                                         +qyzi*gqyz(2)))
                  fkd(2) = uxi*guy(2) + uyi*guy(3) + uzi*guy(4)
     &                        - 0.5d0 * (ci*guy(1) + qxxi*guy(5)
     &                            + qyyi*guy(8) + qzzi*guy(10)
     &                            + 2.0d0*(qxyi*guy(6)+qxzi*guy(7)
     &                                         +qyzi*guy(9)))
     &                        - 0.5d0 * (ci*gc(3) + qxxi*gqxx(3)
     &                            + qyyi*gqyy(3) + qzzi*gqzz(3)
     &                            + 2.0d0*(qxyi*gqxy(3)+qxzi*gqxz(3)
     &                                         +qyzi*gqyz(3)))
                  fkd(3) = uxi*guz(2) + uyi*guz(3) + uzi*guz(4)
     &                        - 0.5d0 * (ci*guz(1) + qxxi*guz(5)
     &                            + qyyi*guz(8) + qzzi*guz(10)
     &                            + 2.0d0*(qxyi*guz(6)+qxzi*guz(7)
     &                                         +qyzi*guz(9)))
     &                        - 0.5d0 * (ci*gc(4) + qxxi*gqxx(4)
     &                            + qyyi*gqyy(4) + qzzi*gqzz(4)
     &                            + 2.0d0*(qxyi*gqxy(4)+qxzi*gqxz(4)
     &                                         +qyzi*gqyz(4)))
c
c     scale the self-field by half, such that it sums to one below
c
                  if (i .eq. k) then
                     do j = 1, 3
                        fid(j) = 0.5d0 * fid(j)
                        fkd(j) = 0.5d0 * fkd(j)
                     end do
                  end if
                  do j = 1, 3
                     fieldts(j,i) = fieldts(j,i) + fid(j)
                     fieldts(j,k) = fieldts(j,k) + fkd(j)
                     fieldtps(j,i) = fieldtps(j,i) + fid(j)
                     fieldtps(j,k) = fieldtps(j,k) + fkd(j)
                  end do
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
!$OMP END DO
c
c     add local copies to global variables for OpenMP calculation
c
!$OMP DO
      do i = 1, npole
         do j = 1, 3
            field(j,i) = field(j,i) + fieldt(j,i)
            fieldp(j,i) = fieldp(j,i) + fieldtp(j,i)
            fields(j,i) = fields(j,i) + fieldts(j,i)
            fieldps(j,i) = fieldps(j,i) + fieldtps(j,i)
         end do
      end do
!$OMP END DO
c
c     combine permanent multipole field and GK reaction field
c
!$OMP DO
      do i = 1, npole
         do j = 1, 3
            fields(j,i) = field(j,i) + fields(j,i)
            fieldps(j,i) = fieldp(j,i) + fieldps(j,i)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
      deallocate (fieldt)
      deallocate (fieldtp)
      deallocate (fieldts)
      deallocate (fieldtps)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ufield0d  --  generalized Kirkwood mutual field  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ufield0d" computes the mutual electrostatic field due to
c     induced dipole moments for use with with generalized Kirkwood
c     implicit solvation
c
c
      subroutine ufield0d (field,fieldp,fields,fieldps)
      use sizes
      use atoms
      use gkstuf
      use group
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      use solute
      implicit none
      integer i,j,k,ii,kk
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 fgrp,r,r2
      real*8 rr3,rr5
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 duir,dukr
      real*8 puir,pukr
      real*8 duixs,duiys,duizs
      real*8 puixs,puiys,puizs
      real*8 dukxs,dukys,dukzs
      real*8 pukxs,pukys,pukzs
      real*8 duirs,puirs
      real*8 dukrs,pukrs
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 pdi,pti,pgamma
      real*8 rb2,rbi,rbk
      real*8 dwater,fd
      real*8 gf,gf2,gf3,gf5
      real*8 expterm,expc
      real*8 expc1,dexpc
      real*8 a(0:3,0:2)
      real*8 gux(10),guy(10)
      real*8 guz(10)
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 fids(3),fkds(3)
      real*8 fips(3),fkps(3)
      real*8, allocatable :: dscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8 fields(3,*)
      real*8 fieldps(3,*)
      real*8, allocatable :: fieldt(:,:)
      real*8, allocatable :: fieldtp(:,:)
      real*8, allocatable :: fieldts(:,:)
      real*8, allocatable :: fieldtps(:,:)
      logical proceed
c
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
            fields(j,i) = 0.0d0
            fieldps(j,i) = 0.0d0
         end do
      end do
c
c     set dielectric constant and scaling factor for water
c
      dwater = 78.3d0
      fd = 2.0d0 * (1.0d0-dwater) / (1.0d0+2.0d0*dwater)
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
c
c     set array needed to scale connected atom interactions
c
      do i = 1, n
         dscale(i) = 1.0d0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (fieldt(3,npole))
      allocate (fieldtp(3,npole))
      allocate (fieldts(3,npole))
      allocate (fieldtps(3,npole))
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,pdamp,thole,rborn,
!$OMP& uind,uinp,uinds,uinps,np11,np12,np13,np14,ip11,ip12,ip13,ip14,
!$OMP& u1scale,u2scale,u3scale,u4scale,use_intra,x,y,z,off2,fd,gkc,
!$OMP& field,fieldp,fields,fieldps)
!$OMP& firstprivate(dscale) shared(fieldt,fieldtp,fieldts,fieldtps)
c
c     initialize local variables for OpenMP calculation
c
!$OMP DO collapse(2)
      do i = 1, npole
         do j = 1, 3
            fieldt(j,i) = 0.0d0
            fieldtp(j,i) = 0.0d0
            fieldts(j,i) = 0.0d0
            fieldtps(j,i) = 0.0d0
         end do
      end do
!$OMP END DO
c
c     find the field terms for each pairwise interaction
c
!$OMP DO reduction(+:fieldt,fieldtp,fieldts,fieldtps)
!$OMP& schedule(guided)
c
c     compute the mutual electrostatic field at each atom,
c     and another field including RF due to induced dipoles
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         rbi = rborn(ii)
         duix = uind(1,i)
         duiy = uind(2,i)
         duiz = uind(3,i)
         puix = uinp(1,i)
         puiy = uinp(2,i)
         puiz = uinp(3,i)
         duixs = uinds(1,i)
         duiys = uinds(2,i)
         duizs = uinds(3,i)
         puixs = uinps(1,i)
         puiys = uinps(2,i)
         puizs = uinps(3,i)
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = u4scale
         end do
         do k = i, npole
            kk = ipole(k)
            rbk = rborn(kk)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               xr2 = xr * xr
               yr2 = yr * yr
               zr2 = zr * zr
               r2 = xr2 + yr2 + zr2
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  dukx = uind(1,k)
                  duky = uind(2,k)
                  dukz = uind(3,k)
                  pukx = uinp(1,k)
                  puky = uinp(2,k)
                  pukz = uinp(3,k)
                  dukxs = uinds(1,k)
                  dukys = uinds(2,k)
                  dukzs = uinds(3,k)
                  pukxs = uinps(1,k)
                  pukys = uinps(2,k)
                  pukzs = uinps(3,k)
                  if (i .ne. k) then
                     scale3 = dscale(kk)
                     scale5 = dscale(kk)
                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(k))
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           expdamp = exp(damp)
                           scale3 = scale3 * (1.0d0-expdamp)
                           scale5 = scale5 * (1.0d0-(1.0d0-damp)
     &                                           *expdamp)
                        end if
                     end if
                     rr3 = scale3 / (r*r2)
                     rr5 = 3.0d0 * scale5 / (r*r2*r2)
                     duir = xr*duix + yr*duiy + zr*duiz
                     dukr = xr*dukx + yr*duky + zr*dukz
                     puir = xr*puix + yr*puiy + zr*puiz
                     pukr = xr*pukx + yr*puky + zr*pukz
                     duirs = xr*duixs + yr*duiys + zr*duizs
                     dukrs = xr*dukxs + yr*dukys + zr*dukzs
                     puirs = xr*puixs + yr*puiys + zr*puizs
                     pukrs = xr*pukxs + yr*pukys + zr*pukzs
                     fid(1) = -rr3*dukx + rr5*dukr*xr
                     fid(2) = -rr3*duky + rr5*dukr*yr
                     fid(3) = -rr3*dukz + rr5*dukr*zr
                     fkd(1) = -rr3*duix + rr5*duir*xr
                     fkd(2) = -rr3*duiy + rr5*duir*yr
                     fkd(3) = -rr3*duiz + rr5*duir*zr
                     fip(1) = -rr3*pukx + rr5*pukr*xr
                     fip(2) = -rr3*puky + rr5*pukr*yr
                     fip(3) = -rr3*pukz + rr5*pukr*zr
                     fkp(1) = -rr3*puix + rr5*puir*xr
                     fkp(2) = -rr3*puiy + rr5*puir*yr
                     fkp(3) = -rr3*puiz + rr5*puir*zr
                     fids(1) = -rr3*dukxs + rr5*dukrs*xr
                     fids(2) = -rr3*dukys + rr5*dukrs*yr
                     fids(3) = -rr3*dukzs + rr5*dukrs*zr
                     fkds(1) = -rr3*duixs + rr5*duirs*xr
                     fkds(2) = -rr3*duiys + rr5*duirs*yr
                     fkds(3) = -rr3*duizs + rr5*duirs*zr
                     fips(1) = -rr3*pukxs + rr5*pukrs*xr
                     fips(2) = -rr3*pukys + rr5*pukrs*yr
                     fips(3) = -rr3*pukzs + rr5*pukrs*zr
                     fkps(1) = -rr3*puixs + rr5*puirs*xr
                     fkps(2) = -rr3*puiys + rr5*puirs*yr
                     fkps(3) = -rr3*puizs + rr5*puirs*zr
                     do j = 1, 3
                        fieldt(j,i) = fieldt(j,i) + fid(j)
                        fieldt(j,k) = fieldt(j,k) + fkd(j)
                        fieldtp(j,i) = fieldtp(j,i) + fip(j)
                        fieldtp(j,k) = fieldtp(j,k) + fkp(j)
                        fieldts(j,i) = fieldts(j,i) + fids(j)
                        fieldts(j,k) = fieldts(j,k) + fkds(j)
                        fieldtps(j,i) = fieldtps(j,i) + fips(j)
                        fieldtps(j,k) = fieldtps(j,k) + fkps(j)
                     end do
                  end if
c
c     unweighted dipole reaction potential gradient tensor
c
                  rb2 = rbi * rbk
                  expterm = exp(-r2/(gkc*rb2))
                  expc = expterm / gkc
                  dexpc = -2.0d0 / (gkc*rbi*rbk)
                  gf2 = 1.0d0 / (r2+rb2*expterm)
                  gf = sqrt(gf2)
                  gf3 = gf2 * gf
                  gf5 = gf3 * gf2
                  a(1,0) = -gf3
                  a(2,0) = 3.0d0 * gf5
                  expc1 = 1.0d0 - expc
                  a(1,1) = expc1 * a(2,0)
                  gux(2) = fd * (a(1,0) + xr2*a(1,1))
                  gux(3) = fd * xr*yr*a(1,1)
                  gux(4) = fd * xr*zr*a(1,1)
                  guy(2) = gux(3)
                  guy(3) = fd * (a(1,0) + yr2*a(1,1))
                  guy(4) = fd * yr*zr*a(1,1)
                  guz(2) = gux(4)
                  guz(3) = guy(4)
                  guz(4) = fd * (a(1,0) + zr2*a(1,1))
                  fids(1) = dukxs*gux(2) + dukys*guy(2) + dukzs*guz(2)
                  fids(2) = dukxs*gux(3) + dukys*guy(3) + dukzs*guz(3)
                  fids(3) = dukxs*gux(4) + dukys*guy(4) + dukzs*guz(4)
                  fkds(1) = duixs*gux(2) + duiys*guy(2) + duizs*guz(2)
                  fkds(2) = duixs*gux(3) + duiys*guy(3) + duizs*guz(3)
                  fkds(3) = duixs*gux(4) + duiys*guy(4) + duizs*guz(4)
                  fips(1) = pukxs*gux(2) + pukys*guy(2) + pukzs*guz(2)
                  fips(2) = pukxs*gux(3) + pukys*guy(3) + pukzs*guz(3)
                  fips(3) = pukxs*gux(4) + pukys*guy(4) + pukzs*guz(4)
                  fkps(1) = puixs*gux(2) + puiys*guy(2) + puizs*guz(2)
                  fkps(2) = puixs*gux(3) + puiys*guy(3) + puizs*guz(3)
                  fkps(3) = puixs*gux(4) + puiys*guy(4) + puizs*guz(4)
                  if (i .eq. k) then
                     do j = 1, 3
                        fids(j) = 0.5d0 * fids(j)
                        fkds(j) = 0.5d0 * fkds(j)
                        fips(j) = 0.5d0 * fips(j)
                        fkps(j) = 0.5d0 * fkps(j)
                     end do
                  end if
                  do j = 1, 3
                     fieldts(j,i) = fieldts(j,i) + fids(j)
                     fieldts(j,k) = fieldts(j,k) + fkds(j)
                     fieldtps(j,i) = fieldtps(j,i) + fips(j)
                     fieldtps(j,k) = fieldtps(j,k) + fkps(j)
                  end do
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
!$OMP END DO
c
c     add local copies to global variables for OpenMP calculation
c
!$OMP DO
      do i = 1, npole
         do j = 1, 3
            field(j,i) = field(j,i) + fieldt(j,i)
            fieldp(j,i) = fieldp(j,i) + fieldtp(j,i)
            fields(j,i) = fields(j,i) + fieldts(j,i)
            fieldps(j,i) = fieldps(j,i) + fieldtps(j,i)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (fieldt)
      deallocate (fieldtp)
      deallocate (fieldts)
      deallocate (fieldtps)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine induce0e  --  Poisson-Boltzmann induced dipoles  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "induce0e" computes the induced dipole moments at polarizable
c     sites for Poisson-Boltzmann SCRF and vacuum environments
c
c
      subroutine induce0e
      use sizes
      use atoms
      use inform
      use iounit
      use mpole
      use polar
      use polpot
      use potent
      use units
      use uprior
      implicit none
      integer i,j,k,iter
      integer maxiter
      real*8 polmin
      real*8 eps,epsold
      real*8 epsd,epsp
      real*8 epsds,epsps
      real*8 udsum,upsum
      real*8 ussum,upssum
      real*8 a,ap,as,aps
      real*8 b,bp,bs,bps
      real*8 sum,sump
      real*8 sums,sumps
      real*8, allocatable :: poli(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8, allocatable :: fields(:,:)
      real*8, allocatable :: fieldps(:,:)
      real*8, allocatable :: udir(:,:)
      real*8, allocatable :: udirp(:,:)
      real*8, allocatable :: udirs(:,:)
      real*8, allocatable :: udirps(:,:)
      real*8, allocatable :: rsd(:,:)
      real*8, allocatable :: rsdp(:,:)
      real*8, allocatable :: rsds(:,:)
      real*8, allocatable :: rsdps(:,:)
      real*8, allocatable :: zrsd(:,:)
      real*8, allocatable :: zrsdp(:,:)
      real*8, allocatable :: zrsds(:,:)
      real*8, allocatable :: zrsdps(:,:)
      real*8, allocatable :: conj(:,:)
      real*8, allocatable :: conjp(:,:)
      real*8, allocatable :: conjs(:,:)
      real*8, allocatable :: conjps(:,:)
      real*8, allocatable :: vec(:,:)
      real*8, allocatable :: vecp(:,:)
      real*8, allocatable :: vecs(:,:)
      real*8, allocatable :: vecps(:,:)
      logical done
      character*6 mode
c
c
c     zero out the induced dipoles; uind and uinp are vacuum dipoles,
c     uinds and uinps are Poisson-Boltzmann SCRF dipoles
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
            uinds(j,i) = 0.0d0
            uinps(j,i) = 0.0d0
         end do
      end do
      if (.not.use_polar .or. .not.use_solv)  return
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (field(3,npole))
      allocate (fieldp(3,npole))
      allocate (fields(3,npole))
      allocate (fieldps(3,npole))
      allocate (udir(3,npole))
      allocate (udirp(3,npole))
      allocate (udirs(3,npole))
      allocate (udirps(3,npole))
c
c     compute the direct induced dipole moment at each atom, and
c     another set that also includes RF due to permanent multipoles
c
      call dfield0e (field,fieldp,fields,fieldps)
c
c     set vacuum induced dipoles to polarizability times direct field;
c     SCRF induced dipoles are polarizability times direct field
c     plus the reaction field due to permanent multipoles
c
      do i = 1, npole
         do j = 1, 3
            udir(j,i) = polarity(i) * field(j,i)
            udirp(j,i) = polarity(i) * fieldp(j,i)
            udirs(j,i) = polarity(i) * fields(j,i)
            udirps(j,i) = polarity(i) * fieldps(j,i)
            uind(j,i) = udir(j,i)
            uinp(j,i) = udirp(j,i)
            uinds(j,i) = udirs(j,i)
            uinps(j,i) = udirps(j,i)
         end do
      end do
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 500
         iter = 0
         polmin = 0.00000001d0
         eps = 100.0d0
c
c     estimated induced dipoles from polynomial predictor
c
         if (use_pred .and. nualt.eq.maxualt) then
            do i = 1, npole
               do j = 1, 3
                  udsum = 0.0d0
                  upsum = 0.0d0
                  ussum = 0.0d0
                  upssum = 0.0d0
                  do k = 1, nualt-1
                     udsum = udsum + bpred(k)*udalt(k,j,i)
                     upsum = upsum + bpredp(k)*upalt(k,j,i)
                     ussum = ussum + bpreds(k)*usalt(k,j,i)
                     upssum = upssum + bpredps(k)*upsalt(k,j,i)
                  end do
                  uind(j,i) = udsum
                  uinp(j,i) = upsum
                  uinds(j,i) = ussum
                  uinps(j,i) = upssum
               end do
            end do
         end if
c
c     perform dynamic allocation of some local arrays
c
         allocate (poli(npole))
         allocate (rsd(3,npole))
         allocate (rsdp(3,npole))
         allocate (rsds(3,npole))
         allocate (rsdps(3,npole))
         allocate (zrsd(3,npole))
         allocate (zrsdp(3,npole))
         allocate (zrsds(3,npole))
         allocate (zrsdps(3,npole))
         allocate (conj(3,npole))
         allocate (conjp(3,npole))
         allocate (conjs(3,npole))
         allocate (conjps(3,npole))
         allocate (vec(3,npole))
         allocate (vecp(3,npole))
         allocate (vecs(3,npole))
         allocate (vecps(3,npole))
c
c     set initial conjugate gradient residual and conjugate vector
c
         call ufield0e (field,fieldp,fields,fieldps)
         do i = 1, npole
            poli(i) = max(polmin,polarity(i))
            do j = 1, 3
               rsd(j,i) = (udir(j,i)-uind(j,i))/poli(i)
     &                       + field(j,i)
               rsdp(j,i) = (udirp(j,i)-uinp(j,i))/poli(i)
     &                        + fieldp(j,i)
               rsds(j,i) = (udirs(j,i)-uinds(j,i))/poli(i)
     &                        + fields(j,i)
               rsdps(j,i) = (udirps(j,i)-uinps(j,i))/poli(i)
     &                         + fieldps(j,i)
               zrsd(j,i) = rsd(j,i) * poli(i)
               zrsdp(j,i) = rsdp(j,i) * poli(i)
               zrsds(j,i) = rsds(j,i) * poli(i)
               zrsdps(j,i) = rsdps(j,i) * poli(i)
               conj(j,i) = zrsd(j,i)
               conjp(j,i) = zrsdp(j,i)
               conjs(j,i) = zrsds(j,i)
               conjps(j,i) = zrsdps(j,i)
            end do
         end do
c
c     conjugate gradient iteration of the mutual induced dipoles
c
         do while (.not. done)
            iter = iter + 1
            do i = 1, npole
               do j = 1, 3
                  vec(j,i) = uind(j,i)
                  vecp(j,i) = uinp(j,i)
                  vecs(j,i) = uinds(j,i)
                  vecps(j,i) = uinps(j,i)
                  uind(j,i) = conj(j,i)
                  uinp(j,i) = conjp(j,i)
                  uinds(j,i) = conjs(j,i)
                  uinps(j,i) = conjps(j,i)
               end do
            end do
            call ufield0e (field,fieldp,fields,fieldps)
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = vec(j,i)
                  uinp(j,i) = vecp(j,i)
                  uinds(j,i) = vecs(j,i)
                  uinps(j,i) = vecps(j,i)
                  vec(j,i) = conj(j,i)/poli(i) - field(j,i)
                  vecp(j,i) = conjp(j,i)/poli(i) - fieldp(j,i)
                  vecs(j,i) = conjs(j,i)/poli(i) - fields(j,i)
                  vecps(j,i) = conjps(j,i)/poli(i) - fieldps(j,i)
               end do
            end do
            a = 0.0d0
            ap = 0.0d0
            as = 0.0d0
            aps = 0.0d0
            sum = 0.0d0
            sump = 0.0d0
            sums = 0.0d0
            sumps = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  a = a + conj(j,i)*vec(j,i)
                  ap = ap + conjp(j,i)*vecp(j,i)
                  as = as + conjs(j,i)*vecs(j,i)
                  aps = aps + conjps(j,i)*vecps(j,i)
                  sum = sum + rsd(j,i)*zrsd(j,i)
                  sump = sump + rsdp(j,i)*zrsdp(j,i)
                  sums = sums + rsds(j,i)*zrsds(j,i)
                  sumps = sumps + rsdps(j,i)*zrsdps(j,i)
               end do
            end do
            if (a .ne. 0.0d0)  a = sum / a
            if (ap .ne. 0.0d0)  ap = sump / ap
            if (as .ne. 0.0d0)  as = sums / as
            if (aps .ne. 0.0d0)  aps = sumps / aps
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = uind(j,i) + a*conj(j,i)
                  uinp(j,i) = uinp(j,i) + ap*conjp(j,i)
                  uinds(j,i) = uinds(j,i) + as*conjs(j,i)
                  uinps(j,i) = uinps(j,i) + aps*conjps(j,i)
                  rsd(j,i) = rsd(j,i) - a*vec(j,i)
                  rsdp(j,i) = rsdp(j,i) - ap*vecp(j,i)
                  rsds(j,i) = rsds(j,i) - as*vecs(j,i)
                  rsdps(j,i) = rsdps(j,i) - aps*vecps(j,i)
               end do
            end do
            b = 0.0d0
            bp = 0.0d0
            bs = 0.0d0
            bps = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  zrsd(j,i) = rsd(j,i) * poli(i)
                  zrsdp(j,i) = rsdp(j,i) * poli(i)
                  zrsds(j,i) = rsds(j,i) * poli(i)
                  zrsdps(j,i) = rsdps(j,i) * poli(i)
                  b = b + rsd(j,i)*zrsd(j,i)
                  bp = bp + rsdp(j,i)*zrsdp(j,i)
                  bs = bs + rsds(j,i)*zrsds(j,i)
                  bps = bps + rsdps(j,i)*zrsdps(j,i)
               end do
            end do
            if (sum .ne. 0.0d0)  b = b / sum
            if (sump .ne. 0.0d0)  bp = bp / sump
            if (sums .ne. 0.0d0)  bs = bs / sums
            if (sumps .ne. 0.0d0)  bps = bps / sumps
            epsd = 0.0d0
            epsp = 0.0d0
            epsds = 0.0d0
            epsps = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  conj(j,i) = zrsd(j,i) + b*conj(j,i)
                  conjp(j,i) = zrsdp(j,i) + bp*conjp(j,i)
                  conjs(j,i) = zrsds(j,i) + bs*conjs(j,i)
                  conjps(j,i) = zrsdps(j,i) + bps*conjps(j,i)
                  epsd = epsd + rsd(j,i)*rsd(j,i)
                  epsp = epsp + rsdp(j,i)*rsdp(j,i)
                  epsds = epsds + rsds(j,i)*rsds(j,i)
                  epsps = epsps + rsdps(j,i)*rsdps(j,i)
               end do
            end do
c
c     check the convergence of the mutual induced dipoles
c
            epsold = eps
            eps = max(epsd,epsp,epsds,epsps)
            eps = debye * sqrt(eps/dble(npolar))
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. politer)  done = .true.
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (poli)
         deallocate (rsd)
         deallocate (rsdp)
         deallocate (rsds)
         deallocate (rsdps)
         deallocate (zrsd)
         deallocate (zrsdp)
         deallocate (zrsds)
         deallocate (zrsdps)
         deallocate (conj)
         deallocate (conjp)
         deallocate (conjs)
         deallocate (conjps)
         deallocate (vec)
         deallocate (vecp)
         deallocate (vecs)
         deallocate (vecps)
c
c     print the results from the conjugate gradient iteration
c
         if (debug) then
            write (iout,30)  iter,eps
   30       format (/,' Induced Dipoles :',6x,'Iterations',i5,
     &                 6x,'RMS Change',f15.10)
         end if
c
c     terminate the calculation if dipoles failed to converge
c
         if (iter.ge.maxiter .or. eps.gt.epsold) then
            write (iout,40)
   40       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
            call prterr
            call fatal
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (field)
      deallocate (fieldp)
      deallocate (fields)
      deallocate (fieldps)
      deallocate (udir)
      deallocate (udirp)
      deallocate (udirs)
      deallocate (udirps)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine dfield0e  --  Poisson-Boltzmann direct field  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "dfield0e" computes the direct electrostatic field due to
c     permanent multipole moments for use with in Poisson-Boltzmann
c
c
      subroutine dfield0e (field,fieldp,fields,fieldps)
      use sizes
      use atoms
      use couple
      use group
      use mpole
      use pbstuf
      use polar
      use polgrp
      use polpot
      use shunt
      use solute
      implicit none
      integer i,j,k,ii,kk
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 fgrp,r,r2
      real*8 rr3,rr5,rr7
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8 fields(3,*)
      real*8 fieldps(3,*)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      logical proceed
c
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
      end do
c
c     compute the direct electrostatic field at each atom, and
c     another field including RF due to permanent multipoles;
c     note self-interactions for the solute field are skipped
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               xr2 = xr * xr
               yr2 = yr * yr
               zr2 = zr * zr
               r2 = xr2 + yr2 + zr2
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = 1.0d0 - expdamp
                        scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                        scale7 = 1.0d0 - expdamp
     &                              *(1.0d0-damp+0.6d0*damp**2)
                     end if
                  end if
                  rr3 = scale3 / (r*r2)
                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
                  rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
                  dir = dix*xr + diy*yr + diz*zr
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qir = qix*xr + qiy*yr + qiz*zr
                  dkr = dkx*xr + dky*yr + dkz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
                  qkr = qkx*xr + qky*yr + qkz*zr
                  fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkx + 2.0d0*rr5*qkx
                  fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dky + 2.0d0*rr5*qky
                  fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkz + 2.0d0*rr5*qkz
                  fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*dix - 2.0d0*rr5*qix
                  fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diy - 2.0d0*rr5*qiy
                  fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diz - 2.0d0*rr5*qiz
                  do j = 1, 3
                     field(j,i) = field(j,i) + fid(j)*dscale(kk)
                     field(j,k) = field(j,k) + fkd(j)*dscale(kk)
                     fieldp(j,i) = fieldp(j,i) + fid(j)*pscale(kk)
                     fieldp(j,k) = fieldp(j,k) + fkd(j)*pscale(kk)
                  end do
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
c
c     find the Poisson-Boltzmann reaction field at each site
c
      call pbempole
c
c     combine permanent multipole field and PB reaction field
c
      do i = 1, npole
         ii = ipole(i)
         do j = 1, 3
            fields(j,i) = field(j,i) + pbep(j,ii)
            fieldps(j,i) = fieldp(j,i) + pbep(j,ii)
         end do
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ufield0e  --  Poisson-Boltzmann mutual field  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ufield0e" computes the mutual electrostatic field due to
c     induced dipole moments via a Poisson-Boltzmann solver
c
c
      subroutine ufield0e (field,fieldp,fields,fieldps)
      use sizes
      use atoms
      use group
      use mpole
      use pbstuf
      use polar
      use polgrp
      use polpot
      use shunt
      use solute
      implicit none
      integer i,j,k,ii,kk
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 fgrp,r,r2
      real*8 rr3,rr5
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 duir,puir
      real*8 dukr,pukr
      real*8 duixs,duiys,duizs
      real*8 puixs,puiys,puizs
      real*8 dukxs,dukys,dukzs
      real*8 pukxs,pukys,pukzs
      real*8 duirs,puirs
      real*8 dukrs,pukrs
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 fids(3),fkds(3)
      real*8 fips(3),fkps(3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8 fields(3,*)
      real*8 fieldps(3,*)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: indpole(:,:)
      real*8, allocatable :: inppole(:,:)
      logical proceed
c
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
            fields(j,i) = 0.0d0
            fieldps(j,i) = 0.0d0
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
c
c     set array needed to scale connected atom interactions
c
      do i = 1, n
         dscale(i) = 1.0d0
      end do
c
c     compute the mutual electrostatic field at each atom,
c     and another field including RF due to induced dipoles
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         duix = uind(1,i)
         duiy = uind(2,i)
         duiz = uind(3,i)
         puix = uinp(1,i)
         puiy = uinp(2,i)
         puiz = uinp(3,i)
         duixs = uinds(1,i)
         duiys = uinds(2,i)
         duizs = uinds(3,i)
         puixs = uinps(1,i)
         puiys = uinps(2,i)
         puizs = uinps(3,i)
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = u4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               xr2 = xr * xr
               yr2 = yr * yr
               zr2 = zr * zr
               r2 = xr2 + yr2 + zr2
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  dukx = uind(1,k)
                  duky = uind(2,k)
                  dukz = uind(3,k)
                  pukx = uinp(1,k)
                  puky = uinp(2,k)
                  pukz = uinp(3,k)
                  dukxs = uinds(1,k)
                  dukys = uinds(2,k)
                  dukzs = uinds(3,k)
                  pukxs = uinps(1,k)
                  pukys = uinps(2,k)
                  pukzs = uinps(3,k)
                  scale3 = dscale(kk)
                  scale5 = dscale(kk)
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = scale3 * (1.0d0-expdamp)
                        scale5 = scale5 * (1.0d0-(1.0d0-damp)*expdamp)
                     end if
                  end if
                  rr3 = scale3 / (r*r2)
                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
                  duir = xr*duix + yr*duiy + zr*duiz
                  dukr = xr*dukx + yr*duky + zr*dukz
                  puir = xr*puix + yr*puiy + zr*puiz
                  pukr = xr*pukx + yr*puky + zr*pukz
                  duirs = xr*duixs + yr*duiys + zr*duizs
                  dukrs = xr*dukxs + yr*dukys + zr*dukzs
                  puirs = xr*puixs + yr*puiys + zr*puizs
                  pukrs = xr*pukxs + yr*pukys + zr*pukzs
                  fid(1) = -rr3*dukx + rr5*dukr*xr
                  fid(2) = -rr3*duky + rr5*dukr*yr
                  fid(3) = -rr3*dukz + rr5*dukr*zr
                  fkd(1) = -rr3*duix + rr5*duir*xr
                  fkd(2) = -rr3*duiy + rr5*duir*yr
                  fkd(3) = -rr3*duiz + rr5*duir*zr
                  fip(1) = -rr3*pukx + rr5*pukr*xr
                  fip(2) = -rr3*puky + rr5*pukr*yr
                  fip(3) = -rr3*pukz + rr5*pukr*zr
                  fkp(1) = -rr3*puix + rr5*puir*xr
                  fkp(2) = -rr3*puiy + rr5*puir*yr
                  fkp(3) = -rr3*puiz + rr5*puir*zr
                  fids(1) = -rr3*dukxs + rr5*dukrs*xr
                  fids(2) = -rr3*dukys + rr5*dukrs*yr
                  fids(3) = -rr3*dukzs + rr5*dukrs*zr
                  fkds(1) = -rr3*duixs + rr5*duirs*xr
                  fkds(2) = -rr3*duiys + rr5*duirs*yr
                  fkds(3) = -rr3*duizs + rr5*duirs*zr
                  fips(1) = -rr3*pukxs + rr5*pukrs*xr
                  fips(2) = -rr3*pukys + rr5*pukrs*yr
                  fips(3) = -rr3*pukzs + rr5*pukrs*zr
                  fkps(1) = -rr3*puixs + rr5*puirs*xr
                  fkps(2) = -rr3*puiys + rr5*puirs*yr
                  fkps(3) = -rr3*puizs + rr5*puirs*zr
                  do j = 1, 3
                     field(j,i) = field(j,i) + fid(j)
                     field(j,k) = field(j,k) + fkd(j)
                     fieldp(j,i) = fieldp(j,i) + fip(j)
                     fieldp(j,k) = fieldp(j,k) + fkp(j)
                     fields(j,i) = fields(j,i) + fids(j)
                     fields(j,k) = fields(j,k) + fkds(j)
                     fieldps(j,i) = fieldps(j,i) + fips(j)
                     fieldps(j,k) = fieldps(j,k) + fkps(j)
                  end do
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(pbeuind))  allocate (pbeuind(3,n))
      if (.not. allocated(pbeuinp))  allocate (pbeuinp(3,n))
c
c     perform dynamic allocation of some local arrays
c
      allocate (indpole(3,n))
      allocate (inppole(3,n))
c
c     zero out the PB reaction field at each atomic site
c
      do i = 1, n
         do j = 1, 3
            indpole(j,i) = 0.0d0
            inppole(j,i) = 0.0d0
            pbeuind(j,i) = 0.0d0
            pbeuinp(j,i) = 0.0d0
         end do
      end do
c
c     find the Poisson-Boltzmann reaction field at each site
c
      do i = 1, npole
         ii = ipole(i)
         do j = 1, 3
            indpole(j,ii) = uinds(j,i)
            inppole(j,ii) = uinps(j,i)
         end do
      end do
      call apbsinduce (indpole,pbeuind)
      call apbsnlinduce (inppole,pbeuinp)
c
c     perform deallocation of some local arrays
c
      deallocate (indpole)
      deallocate (inppole)
c
c     combine mutual induced dipole field and PB reaction field
c
      do i = 1, npole
         ii = ipole(i)
         do j = 1, 3
            fields(j,i) = fields(j,i) + pbeuind(j,ii)
            fieldps(j,i) = fieldps(j,i) + pbeuinp(j,ii)
         end do
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ulspred  --  induced dipole prediction coeffs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ulspred" uses standard extrapolation or a least squares fit
c     to set coefficients of an induced dipole predictor polynomial
c
c     literature references:
c
c     J. Kolafa, "Time-Reversible Always Stable Predictor-Corrector
c     Method for Molecular Dynamics of Polarizable Molecules", Journal
c     of Computational Chemistry, 25, 335-342 (2004)
c
c     W. Wang and R. D. Skeel, "Fast Evaluation of Polarizable Forces",
c     Journal of Chemical Physics, 123, 164107 (2005)
c
c
      subroutine ulspred
      use sizes
      use mpole
      use uprior
      implicit none
      integer i,j,k,m
      real*8 coeff,udk,upk
      real*8 amax,apmax
      real*8 b(maxualt)
      real*8 bp(maxualt)
      real*8 a(maxualt*(maxualt+1)/2)
      real*8 ap(maxualt*(maxualt+1)/2)
      real*8 c(maxualt,maxualt)
      real*8 cp(maxualt,maxualt)
c
c
c     set the Gear predictor binomial coefficients
c
      if (polpred .eq. 'GEAR') then
         do i = 1, nualt
            coeff = gear(i)
            bpred(i) = coeff
            bpredp(i) = coeff
            bpreds(i) = coeff
            bpredps(i) = coeff
         end do
c
c     set always stable predictor-corrector (ASPC) coefficients
c
      else if (polpred .eq. 'ASPC') then
         do i = 1, nualt
            coeff = aspc(i)
            bpred(i) = coeff
            bpredp(i) = coeff
            bpreds(i) = coeff
            bpredps(i) = coeff
         end do
c
c     derive normal equations corresponding to least squares fit
c
      else
         do k = 1, nualt
            b(k) = 0.0d0
            bp(k) = 0.0d0
            do m = k, nualt
               c(k,m) = 0.0d0
               cp(k,m) = 0.0d0
            end do
         end do
         do i = 1, npole
            do j = 1, 3
               do k = 1, nualt
                  udk = udalt(k,j,i)
                  upk = upalt(k,j,i)
                  do m = k, nualt
                     c(k,m) = c(k,m) + udk*udalt(m,j,i)
                     cp(k,m) = cp(k,m) + upk*upalt(m,j,i)
                  end do
               end do
            end do
         end do
         i = 0
         do k = 2, nualt
            b(k-1) = c(1,k)
            bp(k-1) = cp(1,k)
            do m = k, nualt
               i = i + 1
               a(i) = c(k,m)
               ap(i) = cp(k,m)
            end do
         end do
c
c     check for nonzero coefficients and solve normal equations
c
         k = nualt - 1
         amax = 0.0d0
         apmax = 0.0d0
         do i = 1, k*(k+1)/2
            amax = max(amax,a(i))
            apmax = max(apmax,ap(i))
         end do
         if (amax .ne. 0.0d0)  call cholesky (k,a,b)
         if (apmax .ne. 0.0d0)  call cholesky (k,ap,bp)
c
c     transfer the final solution to the coefficient vector
c
         do k = 1, nualt-1
            bpred(k) = b(k)
            bpredp(k) = bp(k)
            bpreds(k) = b(k)
            bpredps(k) = bp(k)
         end do
         bpred(nualt) = 0.0d0
         bpredp(nualt) = 0.0d0
         bpreds(nualt) = 0.0d0
         bpredps(nualt) = 0.0d0
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine uscale0a  --  dipole preconditioner via loop  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "uscale0a" builds and applies a preconditioner for the conjugate
c     gradient induced dipole solver using a double loop
c
c
      subroutine uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
      use sizes
      use atoms
      use limits
      use mpole
      use polar
      use polgrp
      use polpot
      use usolve
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5
      real*8 pdi,pti
      real*8 polmin
      real*8 poli,polik
      real*8 damp,expdamp
      real*8 pgamma,off2
      real*8 scale3,scale5
      real*8 m1,m2,m3
      real*8 m4,m5,m6
      real*8, allocatable :: dscale(:)
      real*8 rsd(3,*)
      real*8 rsdp(3,*)
      real*8 zrsd(3,*)
      real*8 zrsdp(3,*)
      character*6 mode
c
c
c     apply the preconditioning matrix to the current residual
c
      if (mode .eq. 'APPLY') then
c
c     use diagonal preconditioner elements as first approximation
c
         polmin = 0.00000001d0
         do i = 1, npole
            poli = udiag * max(polmin,polarity(i))
            do j = 1, 3
               zrsd(j,i) = poli * rsd(j,i)
               zrsdp(j,i) = poli * rsdp(j,i)
            end do
         end do
c
c     use the off-diagonal preconditioner elements in second phase
c
         off2 = usolvcut * usolvcut
         j = 0
         do i = 1, npole-1
            ii = ipole(i)
            do k = i+1, npole
               kk = ipole(k)
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  m1 = minv(j+1)
                  m2 = minv(j+2)
                  m3 = minv(j+3)
                  m4 = minv(j+4)
                  m5 = minv(j+5)
                  m6 = minv(j+6)
                  j = j + 6
                  zrsd(1,i) = zrsd(1,i) + m1*rsd(1,k) + m2*rsd(2,k)
     &                           + m3*rsd(3,k)
                  zrsd(2,i) = zrsd(2,i) + m2*rsd(1,k) + m4*rsd(2,k)
     &                           + m5*rsd(3,k)
                  zrsd(3,i) = zrsd(3,i) + m3*rsd(1,k) + m5*rsd(2,k)
     &                           + m6*rsd(3,k)
                  zrsd(1,k) = zrsd(1,k) + m1*rsd(1,i) + m2*rsd(2,i)
     &                           + m3*rsd(3,i)
                  zrsd(2,k) = zrsd(2,k) + m2*rsd(1,i) + m4*rsd(2,i)
     &                           + m5*rsd(3,i)
                  zrsd(3,k) = zrsd(3,k) + m3*rsd(1,i) + m5*rsd(2,i)
     &                           + m6*rsd(3,i)
                  zrsdp(1,i) = zrsdp(1,i) + m1*rsdp(1,k) + m2*rsdp(2,k)
     &                            + m3*rsdp(3,k)
                  zrsdp(2,i) = zrsdp(2,i) + m2*rsdp(1,k) + m4*rsdp(2,k)
     &                            + m5*rsdp(3,k)
                  zrsdp(3,i) = zrsdp(3,i) + m3*rsdp(1,k) + m5*rsdp(2,k)
     &                            + m6*rsdp(3,k)
                  zrsdp(1,k) = zrsdp(1,k) + m1*rsdp(1,i) + m2*rsdp(2,i)
     &                            + m3*rsdp(3,i)
                  zrsdp(2,k) = zrsdp(2,k) + m2*rsdp(1,i) + m4*rsdp(2,i)
     &                            + m5*rsdp(3,i)
                  zrsdp(3,k) = zrsdp(3,k) + m3*rsdp(1,i) + m5*rsdp(2,i)
     &                            + m6*rsdp(3,i)
               end if
            end do
         end do
c
c     construct off-diagonal elements of preconditioning matrix
c
      else if (mode .eq. 'BUILD') then
c
c     perform dynamic allocation of some local arrays
c
         allocate (dscale(n))
c
c     set array needed to scale connected atom interactions
c
         do i = 1, n
            dscale(i) = 1.0d0
         end do
c
c     determine the off-diagonal elements of the preconditioner
c
         off2 = usolvcut * usolvcut
         m = 0
         do i = 1, npole-1
            ii = ipole(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
            pdi = pdamp(i)
            pti = thole(i)
            poli = polarity(i)
            do j = i+1, npole
               dscale(ipole(j)) = 1.0d0
            end do
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = u1scale
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = u2scale
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = u3scale
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = u4scale
            end do
            do k = i+1, npole
               kk = ipole(k)
               xr = x(kk) - xi
               yr = y(kk) - yi
               zr = z(kk) - zi
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  scale3 = dscale(kk)
                  scale5 = dscale(kk)
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = scale3 * (1.0d0-expdamp)
                        scale5 = scale5 * (1.0d0-expdamp*(1.0d0-damp))
                     end if
                  end if
                  polik = poli * polarity(k)
                  rr3 = scale3 * polik / (r*r2)
                  rr5 = 3.0d0 * scale5 * polik / (r*r2*r2)
                  minv(m+1) = rr5*xr*xr - rr3
                  minv(m+2) = rr5*xr*yr
                  minv(m+3) = rr5*xr*zr
                  minv(m+4) = rr5*yr*yr - rr3
                  minv(m+5) = rr5*yr*zr
                  minv(m+6) = rr5*zr*zr - rr3
                  m = m + 6
               end if
            end do
c
c     reset interaction scaling coefficients for connected atoms
c
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = 1.0d0
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = 1.0d0
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = 1.0d0
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = 1.0d0
            end do
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (dscale)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine uscale0b  --  dipole preconditioner via list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "uscale0b" builds and applies a preconditioner for the conjugate
c     gradient induced dipole solver using a neighbor pair list
c
c
      subroutine uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
      use sizes
      use atoms
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use usolve
      implicit none
      integer i,j,k,m
      integer ii,kk,kkk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5
      real*8 pdi,pti
      real*8 polmin
      real*8 poli,polik
      real*8 damp,expdamp
      real*8 pgamma
      real*8 scale3,scale5
      real*8 m1,m2,m3
      real*8 m4,m5,m6
      real*8, allocatable :: dscale(:)
      real*8 rsd(3,*)
      real*8 rsdp(3,*)
      real*8 zrsd(3,*)
      real*8 zrsdp(3,*)
      real*8, allocatable :: zrsdt(:,:)
      real*8, allocatable :: zrsdtp(:,:)
      character*6 mode
c
c
c     apply the preconditioning matrix to the current residual
c
      if (mode .eq. 'APPLY') then
c
c     perform dynamic allocation of some local arrays
c
         allocate (zrsdt(3,npole))
         allocate (zrsdtp(3,npole))
c
c     use diagonal preconditioner elements as first approximation
c
         polmin = 0.00000001d0
         do i = 1, npole
            poli = udiag * max(polmin,polarity(i))
            do j = 1, 3
               zrsd(j,i) = 0.0d0
               zrsdp(j,i) = 0.0d0
               zrsdt(j,i) = poli * rsd(j,i)
               zrsdtp(j,i) = poli * rsdp(j,i)
            end do
         end do
c
c     use the off-diagonal preconditioner elements in second phase
c
!$OMP PARALLEL default(private) shared(npole,mindex,minv,nulst,ulst,
!$OMP& rsd,rsdp,zrsd,zrsdp,zrsdt,zrsdtp)
!$OMP DO reduction(+:zrsdt,zrsdtp) schedule(guided)
         do i = 1, npole
            m = mindex(i)
            do kk = 1, nulst(i)
               k = ulst(kk,i)
               m1 = minv(m+1)
               m2 = minv(m+2)
               m3 = minv(m+3)
               m4 = minv(m+4)
               m5 = minv(m+5)
               m6 = minv(m+6)
               m = m + 6
               zrsdt(1,i) = zrsdt(1,i) + m1*rsd(1,k) + m2*rsd(2,k)
     &                        + m3*rsd(3,k)
               zrsdt(2,i) = zrsdt(2,i) + m2*rsd(1,k) + m4*rsd(2,k)
     &                        + m5*rsd(3,k)
               zrsdt(3,i) = zrsdt(3,i) + m3*rsd(1,k) + m5*rsd(2,k)
     &                        + m6*rsd(3,k)
               zrsdt(1,k) = zrsdt(1,k) + m1*rsd(1,i) + m2*rsd(2,i)
     &                        + m3*rsd(3,i)
               zrsdt(2,k) = zrsdt(2,k) + m2*rsd(1,i) + m4*rsd(2,i)
     &                        + m5*rsd(3,i)
               zrsdt(3,k) = zrsdt(3,k) + m3*rsd(1,i) + m5*rsd(2,i)
     &                        + m6*rsd(3,i)
               zrsdtp(1,i) = zrsdtp(1,i) + m1*rsdp(1,k) + m2*rsdp(2,k)
     &                         + m3*rsdp(3,k)
               zrsdtp(2,i) = zrsdtp(2,i) + m2*rsdp(1,k) + m4*rsdp(2,k)
     &                         + m5*rsdp(3,k)
               zrsdtp(3,i) = zrsdtp(3,i) + m3*rsdp(1,k) + m5*rsdp(2,k)
     &                         + m6*rsdp(3,k)
               zrsdtp(1,k) = zrsdtp(1,k) + m1*rsdp(1,i) + m2*rsdp(2,i)
     &                         + m3*rsdp(3,i)
               zrsdtp(2,k) = zrsdtp(2,k) + m2*rsdp(1,i) + m4*rsdp(2,i)
     &                         + m5*rsdp(3,i)
               zrsdtp(3,k) = zrsdtp(3,k) + m3*rsdp(1,i) + m5*rsdp(2,i)
     &                         + m6*rsdp(3,i)
            end do
         end do
!$OMP END DO
c
c     transfer the results from local to global arrays
c
!$OMP DO
         do i = 1, npole
            do j = 1, 3
               zrsd(j,i) = zrsdt(j,i) + zrsd(j,i)
               zrsdp(j,i) = zrsdtp(j,i) + zrsdp(j,i)
            end do
         end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
         deallocate (zrsdt)
         deallocate (zrsdtp)
c
c     build the off-diagonal elements of preconditioning matrix
c
      else if (mode .eq. 'BUILD') then
         m = 0
         do i = 1, npole
            mindex(i) = m
            m = m + 6*nulst(i)
         end do
c
c     perform dynamic allocation of some local arrays
c
         allocate (dscale(n))
c
c     set array needed to scale connected atom interactions
c
         do i = 1, n
            dscale(i) = 1.0d0
         end do
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(n,npole,ipole,x,y,z,pdamp,
!$OMP& thole,polarity,u1scale,u2scale,u3scale,u4scale,np11,ip11,
!$OMP& np12,ip12,np13,ip13,np14,ip14,nulst,ulst,mindex,minv)
!$OMP& firstprivate (dscale)
c
c     determine the off-diagonal elements of the preconditioner
c
!$OMP DO schedule(guided)

         do i = 1, npole
            ii = ipole(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
            pdi = pdamp(i)
            pti = thole(i)
            poli = polarity(i)
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = u1scale
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = u2scale
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = u3scale
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = u4scale
            end do
            m = mindex(i)
            do kkk = 1, nulst(i)
               k = ulst(kkk,i)
               kk = ipole(k)
               xr = x(kk) - xi
               yr = y(kk) - yi
               zr = z(kk) - zi
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               r = sqrt(r2)
               scale3 = dscale(kk)
               scale5 = dscale(kk)
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = scale3 * (1.0d0-expdamp)
                     scale5 = scale5 * (1.0d0-expdamp*(1.0d0-damp))
                  end if
               end if
               polik = poli * polarity(k)
               rr3 = scale3 * polik / (r*r2)
               rr5 = 3.0d0 * scale5 * polik / (r*r2*r2)
               minv(m+1) = rr5*xr*xr - rr3
               minv(m+2) = rr5*xr*yr
               minv(m+3) = rr5*xr*zr
               minv(m+4) = rr5*yr*yr - rr3
               minv(m+5) = rr5*yr*zr
               minv(m+6) = rr5*zr*zr - rr3
               m = m + 6
            end do
c
c     reset interaction scaling coefficients for connected atoms
c
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = 1.0d0
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = 1.0d0
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = 1.0d0
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = 1.0d0
            end do
         end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
         deallocate (dscale)
      end if
      return
      end
