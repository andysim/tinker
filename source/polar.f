c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module polar  --  induced dipole moments & polarizability  ##        
c     ##                                                             ##
c     #################################################################
c
c
c     npolar    total number of polarizable sites in the system
c     polarity  dipole polarizability for each multipole site (Ang**3)
c     thole     Thole polarizability damping value for each site
c     pdamp     value of polarizability scale factor for each site
c     uind      induced dipole components at each multipole site
c     uinp      induced dipoles in field used for energy interactions
c     uinds     GK or PB induced dipoles at each multipole site
c     uinps     induced dipoles in field used for GK or PB energy
c OPT IMPLEMENTATION
c     dofit     whether to fit the coeffiecients in the PT series to converged dipoles
c     maxpt     the maximum allowed order of perturbation theory
c     maxord    the maximum order or perturbation theory used in this calculation
c     ptcoefs   the coefficients of each (partial contribution) term in the PT expantion
c     ptcoefsf  the coefficients of each (full, i.e., summed) term in the PT expansion
c     uindgridf the fourier space representation of µd(0), µd(1), µd(2), etc.
c     uinpgridf the fourier space representation of µp(0), µp(1), µp(2), etc.
c     ptfphid   the fractional coordinate rec space potential due to µd(0), µd(1), µd(2), etc.
c     ptfphip   the fractional coordinate rec space potential due to µp(0), µp(1), µp(2), etc.
c     permgridf the fourier space representation of the permanent moments
c     permgridr the fractional coordinate rec space potential due to permanent moments
c     ptpointer a hack used to figure out where to store the current grids
c     includeuinp whether to include uinp when determining optimal PT coeffiecients
c OPT IMPLEMENTATION
c
c
      module polar
      implicit none
      integer npolar
      real*8, allocatable :: polarity(:)
      real*8, allocatable :: thole(:)
      real*8, allocatable :: pdamp(:)
      real*8, allocatable :: uind(:,:)
      real*8, allocatable :: uinp(:,:)
      real*8, allocatable :: uinds(:,:)
      real*8, allocatable :: uinps(:,:)
c OPT IMPLEMENTATION
      integer, parameter :: MAXPT = 8
      integer ptmaxord
      real*8 :: ptcoefs(0:MAXPT), ptcoefsf(0:MAXPT)
      logical :: dofit = .false.
      integer :: includeuinp = 1
      integer :: ptpointer
      real*8, allocatable :: ptuind(:,:,:), ptuinp(:,:,:)
      real*8, allocatable :: ptfphip(:,:,:)
      real*8, allocatable :: ptfphid(:,:,:)
      real*8, allocatable :: fphiperm(:,:)
      real*8, allocatable :: uinpgridf(:,:,:,:,:)
      real*8, allocatable :: uindgridf(:,:,:,:,:)
      real*8, allocatable :: permgridf(:,:,:,:)
c OPT IMPLEMENTATION
      save
      end
