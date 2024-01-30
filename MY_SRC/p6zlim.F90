MODULE p6zlim
   !!======================================================================
   !!                         ***  MODULE p6zlim  ***
   !! TOP :   PISCES-QUOTA : Computes the various nutrient limitation terms
   !!                        of phytoplankton includes explicit diazotroph
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-04  (O. Aumont, C. Ethe) Limitation for iron modelled in quota 
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
   !!----------------------------------------------------------------------
   !!   p6z_lim        :   Compute the nutrients limitation terms 
   !!   p6z_lim_init   :   Read the namelist 
   !!----------------------------------------------------------------------
   USE oce_trc         ! Shared ocean-passive tracers variables
   USE trc             ! Tracers defined
   USE p2zlim          ! Nutrient limitation
   USE p4zlim          ! Nutrient limitation 
   USE sms_pisces      ! PISCES variables
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC p6z_lim           ! called in p4zbio.F90  
   PUBLIC p6z_lim_init      ! called in trcsms_pisces.F90 
   PUBLIC p6z_lim_alloc     ! called in trcini_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::  concpno3    !:  NO3 half saturation for picophyto  
   REAL(wp), PUBLIC ::  concpnh4    !:  NH4 half saturation for picophyto
   REAL(wp), PUBLIC ::  concnpo4    !:  PO4 half saturation for nanophyto
   REAL(wp), PUBLIC ::  concppo4    !:  PO4 half saturation for picophyto
   REAL(wp), PUBLIC ::  concdpo4    !:  PO4 half saturation for diatoms
   REAL(wp), PUBLIC ::  concpfer    !:  Iron half saturation for picophyto
   REAL(wp), PUBLIC ::  concbpo4    !:  PO4 half saturation for bacteria
   REAL(wp), PUBLIC ::  xsizepic    !:  Minimum size criteria for picophyto
   REAL(wp), PUBLIC ::  xsizerp     !:  Size ratio for picophytoplankton
   REAL(wp), PUBLIC ::  qfnopt      !:  optimal Fe quota for nanophyto
   REAL(wp), PUBLIC ::  qfpopt      !:  optimal Fe quota for picophyto
   REAL(wp), PUBLIC ::  qfdopt      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  qnnmin      !:  minimum N  quota for nanophyto
   REAL(wp), PUBLIC ::  qnnmax      !:  maximum N quota for nanophyto
   REAL(wp), PUBLIC ::  qpnmin      !:  minimum P quota for nanophyto
   REAL(wp), PUBLIC ::  qpnmax      !:  maximum P quota for nanophyto
   REAL(wp), PUBLIC ::  qnpmin      !:  minimum N quota for nanophyto
   REAL(wp), PUBLIC ::  qnpmax      !:  maximum N quota for nanophyto
   REAL(wp), PUBLIC ::  qppmin      !:  minimum P quota for nanophyto
   REAL(wp), PUBLIC ::  qppmax      !:  maximum P quota for nanophyto
   REAL(wp), PUBLIC ::  qndmin      !:  minimum N quota for diatoms
   REAL(wp), PUBLIC ::  qndmax      !:  maximum N quota for diatoms
   REAL(wp), PUBLIC ::  qpdmin      !:  minimum P quota for diatoms
   REAL(wp), PUBLIC ::  qpdmax      !:  maximum P quota for diatoms
   REAL(wp), PUBLIC ::  qfnmax      !:  maximum Fe quota for nanophyto
   REAL(wp), PUBLIC ::  qfpmax      !:  maximum Fe quota for picophyto
   REAL(wp), PUBLIC ::  qfdmax      !:  maximum Fe quota for diatoms
   REAL(wp), PUBLIC ::  xpsinh4     !:  respiration cost of NH4 assimilation
   REAL(wp), PUBLIC ::  xpsino3     !:  respiration cost of NO3 assimilation
   REAL(wp), PUBLIC ::  xpsiuptk    !:  Mean respiration cost
   ! Diazotrophy
   REAL(wp), PUBLIC ::  concdzno3   !:  Nitrate half saturation for diazo
   REAL(wp), PUBLIC ::  concdznh4   !:  NH4 half saturation for diazo
   REAL(wp), PUBLIC ::  concdzpo4   !:  PO4 half saturation for diazo
   REAL(wp), PUBLIC ::  concdzfer   !:  Fe half saturation for diazo
   REAL(wp), PUBLIC ::  xsizedz     !:  Minimum size criteria for phyto
   REAL(wp), PUBLIC ::  xsizerdz    !:  size ratio of diazotrophs
   REAL(wp), PUBLIC ::  qfdzopt     !:  optimal Fe quota of diazotrophs 
   REAL(wp), PUBLIC ::  qndzmin     !:  Minimal N quota of diazotrophs
   REAL(wp), PUBLIC ::  qndzmax     !:  Maximal N quota of diazotrophs
   REAL(wp), PUBLIC ::  qpdzmin     !:  Minimal P quota of diazotrophs
   REAL(wp), PUBLIC ::  qpdzmax     !:  Maximal P quota of diazotrophs
   REAL(wp), PUBLIC ::  qfdzmax     !:  Maximal Fe quota of diazotrophs
   REAL(wp), PUBLIC ::  xpsinfix    !:  Cost of biosynthesis associated with Nfix
   REAL(wp), PUBLIC ::  xkdop       !:  half saturation of DOP uptake phytos
   REAL(wp), PUBLIC ::  xkdopdz     !:  half saturation of DOP uptake diazos
   REAL(wp), PUBLIC ::  Facul_lim   !:  Diazos Facultative sensitivity
   REAL(wp), PUBLIC ::  kustkaFe    !: Diazo Fe limitation based on kustka 2003
   REAL(wp), PUBLIC ::  maxFescale  !: MAximum Fe cost of nfix can be scaled
   REAL(wp), PUBLIC ::  maxPminscale!: Maximum Pmin scaling 
   LOGICAL , PUBLIC ::  ln_tiue     !: Boolean to activate temperature dependence on Fe cost of nfix (IUE)
   LOGICAL , PUBLIC ::  ln_tpue     !: Boolean to activate temperature dependence on QPmin for diazo (PUE)

   !!*  Allometric variations of the quotas
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqnnmin    !: Minimum N quota of nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqnnmax    !: Maximum N quota of nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpnmin    !: Minimum P quota of nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpnmax    !: Maximum P quota of picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqnpmin    !: Minimum N quota of picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqnpmax    !: Maximum N quota of picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqppmin    !: Minimum P quota of picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqppmax    !: Maximum P quota of picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqndmin    !: Minimum N quota of diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqndmax    !: Maximum N quota of diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpdmin    !: Minimum P quota of diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpdmax    !: Maximum P quota of diatoms
   ! Diazotrophy
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqndzmin   !: Minimum N quota of diazos
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqndzmax   !: Maximum N quota of diazos
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpdzmin   !: Minimum P quota of diazos
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpdzmax   !: Maximum P quota of diazos
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   qfemindz   !: QFe min of diazotrophs

   !!* Phytoplankton nutrient limitation terms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpicono3   !: Limitation of NO3 uptake by picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpiconh4   !: Limitation of NH4 uptake by picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpicopo4   !: Limitation of PO4 uptake by picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanodop   !: Limitation of DOP uptake by nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpicodop   !: Limitation of DOP uptake by picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatdop   !: Limitation of DOP uptake by diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpicofer   !: Limitation of Fe uptake by picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimpic    !: Limitation of picophyto PP by nutrients
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimpics   !: Limitation of picophyto PP by nutrients
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimphys   !: Limitation of nanophyto PP by nutrients
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdias   !: Limitation of diatoms PP by nutrients
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimpfe    !: Limitation of picophyto PP by Fe
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   fvnuptk    !: Maximum potential uptake rate of nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   fvpuptk    !: Maximum potential uptake rate of picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   fvduptk    !: Maximum potential uptake rate of diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xqfuncfecp !: 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimnpn, xlimnpp, xlimnpd
   !Diazotrophy
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiazno3   !: Limitation of NO3 uptake by diazos
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiaznh4   !: Limitation of NH4 uptake by diazos
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiazpo4   !: Limitation of PO4 uptake by diazos
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiazdop   !: Limitation of DOP uptake by diazos
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiazfer   !: Limitation of Fe uptake by diazos
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdiaz   !: Limitation of diazo PP by C and N (P by proxy)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdzfe   !: Limitation of diazo PP by Fe
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   fvdzuptk   !: 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdzp    !: Limitation of diazo Nfix by P
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdiazo  !: Limitation of diazo PP by nutrients
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   qfenfixdz  !: Fe cost of nitrogen fixation
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimnpdz    !: Limitation of Nfix by N
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   Fe_scale_nfix !: Temp dependant Fe cost scaling
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   Pmin_scale_nfix !: Temp dependent QPmin scaling
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xqfuncfecdz !: 
   ! Coefficient for iron limitation following Flynn and Hipkin (1999)
   REAL(wp) ::  xcoef1   = 0.00167  / 55.85
   REAL(wp) ::  xcoef2   = 1.21E-5 * 14. / 55.85 / 7.625 * 0.5 * 1.5
   REAL(wp) ::  xcoef3   = 1.15E-4 * 14. / 55.85 / 7.625 * 0.5 
   ! Diazotrophy
   REAL(wp) ::  xcoef4   = 7.5E-4  * 14. / 55.85 / 7.625 * 0.5

    LOGICAL  :: l_dia_nut_lim, l_dia_iron_lim, l_dia_fracal
    LOGICAL  :: l_dia_size_lim, l_dia_size_pro
   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p6zlim.F90 10070 2018-08-28 14:30:54Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p6z_lim( kt, knt, Kbb, Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p6z_lim  ***
      !!
      !! ** Purpose :   Compute the co-limitations by the various nutrients
      !!                for the various phytoplankton species. Quota based
      !!                approach. The quota model is derived from theoretical
      !!                models proposed by Pahlow and Oschlies (2009) and 
      !!                Flynn (2001). Various adaptations from several 
      !!                publications by these authors have been also adopted.
      !!                Explicit Diazotrophy Added 
      !!
      !! ** Method  : Quota based approach. The quota model is derived from 
      !!              theoretical models by Pahlow and Oschlies (2009) and 
      !!              Flynn (2001). Various adaptations from several publications
      !!              by these authors have been also adopted.
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in)  :: kt, knt
      INTEGER, INTENT(in)  :: Kbb, Kmm  ! time level indices
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zlim1, zlim2, zlim3, zlim4, zno3, zferlim
      REAL(wp) ::   z1_trndia, z1_trnpic, z1_trnphy, ztem1, ztem2, zetot1
      REAL(wp) ::   zratio, zration, zratiof, znutlim, zfalim, zxpsiuptk
      REAL(wp) ::   zconc1d, zconc1dnh4, zconc0n, zconc0nnh4, zconc0npo4, zconc0dpo4
      REAL(wp) ::   zconc0p, zconc0pnh4, zconc0ppo4, zconcpfe, zconcnfe, zconcdfe
      REAL(wp) ::   fanano, fananop, fananof, fadiat, fadiatp, fadiatf
      REAL(wp) ::   fapico, fapicop, fapicof, zlimpo4, zlimdop
      REAL(wp) ::   zrpho, zrass, zcoef, zfuptk, ztrn, ztrp
      REAL(wp) ::   zfvn, zfvp, zfvf, zsizen, zsizep, zsized, znanochl, zpicochl, zdiatchl
      REAL(wp) ::   zqfemn, zqfemp, zqfemd, zbiron
      REAL(wp) ::   znutlimtot, zlimno3, zlimnh4, zlim1f, zsizetmp
      !! Explicit diazotroph PFT
      REAL(wp) ::   z1_trndiaz, zconc0dz, zconc0dznh4, zconc0dzpo4, zconcdzfe
      REAL(wp) ::   fadiaz, fadiazp, fadiazf, zsizedz, zdiazchl, zqfemdz, zratiop, qfenfix
      REAL(wp) ::   facul, IUE_scale, PUE_scale
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p6z_lim')

      IF( kt == nittrc000 )  THEN
         l_dia_nut_lim  = iom_use( "LNnut"   ) .OR. iom_use( "LDnut" ) .OR. iom_use( "LPnut" ) .OR. iom_use( "LDZnut" )
         l_dia_iron_lim = iom_use( "LNFe"    ) .OR. iom_use( "LDFe"  ) .OR. iom_use( "LPFe"  ) .OR. iom_use( "LDZFe" )
         l_dia_size_lim = iom_use( "SIZEN"   ) .OR. iom_use( "SIZED" ) .OR. iom_use( "SIZEP" ) .OR. iom_use( "SIZEDZ" )
         l_dia_size_pro = iom_use( "RASSN"   ) .OR. iom_use( "RASSP" ) .OR. iom_use( "RASSP" ) .OR. iom_use( "RASSDZ" )
         l_dia_fracal   = iom_use( "xfracal" )
      ENDIF
      !
      sizena(:,:,:) = 0.0  ;  sizepa(:,:,:) = 0.0  ;  sizeda(:,:,:) = 0.0 ; sizedza(:,:,:) = 0.0
      !
      DO_3D( 0, 0, 0, 0, 1, jpkm1)
         ! Computation of the Chl/C ratio of each phytoplankton group
         ! -------------------------------------------------------
         z1_trnphy   = 1. / ( tr(ji,jj,jk,jpphy,Kbb) + rtrn )
         z1_trnpic   = 1. / ( tr(ji,jj,jk,jppic,Kbb) + rtrn )
         z1_trndia   = 1. / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )
         znanochl = tr(ji,jj,jk,jpnch,Kbb) * z1_trnphy
         zpicochl = tr(ji,jj,jk,jppch,Kbb) * z1_trnpic
         zdiatchl = tr(ji,jj,jk,jpdch,Kbb) * z1_trndia
         ! Diazotrophy
         z1_trndiaz   = 1. / ( tr(ji,jj,jk,jpcdz,Kbb) + rtrn )
         zdiazchl = tr(ji,jj,jk,jpchd,Kbb) * z1_trndiaz

         ! Computation of a variable Ks for the different phytoplankton
         ! group as a function of their relative size. Allometry
         ! from Edwards et al. (2012)
         !------------------------------------------------

         ! diatoms
         zsized            = sized(ji,jj,jk)**0.81
         zconcdfe          = concdfer * zsized
         zconc1d           = concdno3 * zsized
         zconc1dnh4        = concdnh4 * zsized
         zconc0dpo4        = concdpo4 * zsized

         ! picophytoplankton
         zsizep            = sizep(ji,jj,jk)**0.81
         zconcpfe          = concpfer * zsizep
         zconc0p           = concpno3 * zsizep
         zconc0pnh4        = concpnh4 * zsizep
         zconc0ppo4        = concppo4 * zsizep

         ! nanophytoplankton
         zsizen            = sizen(ji,jj,jk)**0.81
         zconcnfe          = concnfer * zsizen
         zconc0n           = concnno3 * zsizen
         zconc0nnh4        = concnnh4 * zsizen
         zconc0npo4        = concnpo4 * zsizen

         ! diazotrophs
         zsizedz           = sizedz(ji,jj,jk)**0.81
         zconcdzfe         = concdzfer * zsizedz
         zconc0dz          = concdzno3 * zsizedz
         zconc0dznh4       = concdznh4 * zsizedz
         zconc0dzpo4       = concdzpo4 * zsizedz

         ! Allometric variations of the minimum and maximum quotas
         ! From Talmy et al. (2014) and Maranon et al. (2013)
         ! -------------------------------------------------------
         xqnnmin(ji,jj,jk) = qnnmin * sizen(ji,jj,jk)**(-0.18)
         xqnnmax(ji,jj,jk) = qnnmax
         xqndmin(ji,jj,jk) = qndmin * sized(ji,jj,jk)**(-0.18)
         xqndmax(ji,jj,jk) = qndmax
         xqnpmin(ji,jj,jk) = qnpmin * sizep(ji,jj,jk)**(-0.18)
         xqnpmax(ji,jj,jk) = qnpmax
         ! diazotrophy
         xqndzmin(ji,jj,jk) = qndzmin * sizedz(ji,jj,jk)**(-0.18)
         xqndzmax(ji,jj,jk) = qndzmax
         !
         ! Michaelis-Menten Limitation term for nutrients Small flagellates
         ! -----------------------------------------------
         ztrn    = tr(ji,jj,jk,jpnh4,Kbb) + tr(ji,jj,jk,jpno3,Kbb)
         ztrp    = tr(ji,jj,jk,jppo4,Kbb) + tr(ji,jj,jk,jpdop,Kbb) / 200.0

         ! Computation of the optimal allocation parameters
         ! Based on the different papers by Pahlow et al., and Smith et al.
         ! -----------------------------------------------------------------
         zbiron = ( 75.0 * ( 1.0 - plig(ji,jj,jk) ) + plig(ji,jj,jk) ) * biron(ji,jj,jk)
               
         ! Nanophytoplankton
         znutlim = ztrn / zconc0n
         fanano = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
         znutlim = ztrp / zconc0npo4
         fananop = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
         znutlim = zbiron / zconcnfe
         fananof = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )

         ! Picophytoplankton
         znutlim = ztrn / zconc0p
         fapico = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
         znutlim = ztrp / zconc0npo4
         fapicop = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
         znutlim = zbiron / zconcpfe
         fapicof = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )

         ! Diatoms
         znutlim = ztrn / zconc1d
         fadiat = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
         znutlim = ztrp / zconc0dpo4
         fadiatp = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
         znutlim = zbiron / zconcdfe
         fadiatf = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )

         ! Diazotrophs
!         znutlim = ztrn / zconc0dz
!         fadiaz = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
!         znutlim = ztrp / zconc0dzpo4
!         fadiazp = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
!         znutlim = zbiron / zconcdzfe
!         fadiazf = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
         ! Nemo4 Diazocode
         znutlim = MAX( tr(ji,jj,jk,jpnh4,Kbb) / zconc0dznh4,    &
         &         tr(ji,jj,jk,jpno3,Kbb) / zconc0dz)
         fadiaz  = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
         znutlim = tr(ji,jj,jk,jppo4,Kbb) / zconc0dzpo4
         fadiazp = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
         znutlim = zbiron / zconcdzfe
         fadiazf = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )         

         !
         ! Michaelis-Menten Limitation term by nutrients of
         !  heterotrophic bacteria
         ! -------------------------------------------------------------
         zlim1   = ( tr(ji,jj,jk,jpnh4,Kbb) + tr(ji,jj,jk,jpno3,Kbb) )  &
             &      / ( concbno3 + tr(ji,jj,jk,jpnh4,Kbb) + tr(ji,jj,jk,jpno3,Kbb) )
         !
         zlim2    = tr(ji,jj,jk,jppo4,Kbb) / ( tr(ji,jj,jk,jppo4,Kbb) + concbpo4)
         zlim3    = biron(ji,jj,jk) / ( concbfe + biron(ji,jj,jk) )
         zlim4    = tr(ji,jj,jk,jpdoc,Kbb) / ( xkdoc   + tr(ji,jj,jk,jpdoc,Kbb) )

         ! Xlimbac is used for DOC solubilization whereas xlimbacl
         ! is used for all the other bacterial-dependent term
         ! -------------------------------------------------------
         xlimbacl(ji,jj,jk) = MIN( zlim1, zlim2, zlim3 )
         xlimbac (ji,jj,jk) = xlimbacl(ji,jj,jk) * zlim4
         !
         ! Limitation of N based nutrients uptake (NO3 and NH4)
         zfalim  = (1.-fanano) / fanano
         zlimnh4 = tr(ji,jj,jk,jpnh4,Kbb) / ( zconc0n + tr(ji,jj,jk,jpnh4,Kbb) )
         zlimno3 = tr(ji,jj,jk,jpno3,Kbb) / ( zconc0n + tr(ji,jj,jk,jpno3,Kbb) )
         znutlimtot = (1. - fanano) * ztrn  / ( zfalim * zconc0n + ztrn )
         xnanonh4(ji,jj,jk) = znutlimtot * 5.0 * zlimnh4 / ( zlimno3 + 5.0 * zlimnh4 + rtrn )
         xnanono3(ji,jj,jk) = znutlimtot * zlimno3 / ( zlimno3 + 5.0 * zlimnh4 + rtrn )
         !
         ! Limitation of P based nutrients (PO4 and DOP)
         zfalim  = (1.-fananop) / fananop
         zlimpo4 = tr(ji,jj,jk,jppo4,Kbb) / ( tr(ji,jj,jk,jppo4,Kbb) + zconc0npo4 )
         zlimdop = tr(ji,jj,jk,jpdop,Kbb) / ( tr(ji,jj,jk,jpdop,Kbb) + zconc0npo4 )
         znutlimtot = (1. - fananop) * ztrp / ( zfalim * zconc0npo4 + ztrp )
         xnanopo4(ji,jj,jk) = znutlimtot * 100.0 * zlimpo4 / ( zlimdop + 100.0 * zlimpo4 + rtrn )
         xnanodop(ji,jj,jk) = znutlimtot * zlimdop / ( zlimdop + 100.0 * zlimpo4 + rtrn )
         !
         ! Limitation of Fe uptake
         zfalim = (1.-fananof) / fananof
         xnanofer(ji,jj,jk) = (1. - fananof) * zbiron / ( zbiron + zfalim * zconcnfe )
         !
         ! The minimum iron quota depends on the size of PSU, respiration
         ! and the reduction of nitrate following the parameterization 
         ! proposed by Flynn and Hipkin (1999)
         zratiof   = tr(ji,jj,jk,jpnfe,Kbb) * z1_trnphy
         zqfemn = xcoef1 * znanochl + xcoef2 + xcoef3 * xnanono3(ji,jj,jk)
         xqfuncfecn(ji,jj,jk) = zqfemn + qfnopt
         !
         zration = tr(ji,jj,jk,jpnph,Kbb) * z1_trnphy
         zration = MIN(xqnnmax(ji,jj,jk), MAX( xqnnmin(ji,jj,jk), zration ))
         fvnuptk(ji,jj,jk) = 2.5 * xpsiuptk * xqnnmin(ji,jj,jk) / (zration + rtrn)  &
         &                   * MAX(0., (1. - ratchl * znanochl / 12. ) )
         !
         zlim1  = (zration - xqnnmin(ji,jj,jk) ) / (xqnnmax(ji,jj,jk) - xqnnmin(ji,jj,jk) )

         ! The value of the optimal quota in the formulation below
         ! has been found by solving a non linear equation
         zlim1f = ( 1.13 - xqnnmin(ji,jj,jk) ) / (xqnnmax(ji,jj,jk) - xqnnmin(ji,jj,jk) )
         zlim3  = MAX( 0.,( zratiof - zqfemn ) / qfnopt )
         ! computation of the various limitation terms of nanophyto
         ! growth and PP
         xlimnfe (ji,jj,jk) = MIN( 1., zlim3 )
         xlimphy (ji,jj,jk) = MIN( 1., zlim1, zlim3 )
         xlimphys(ji,jj,jk) = MIN( 1., zlim1/( zlim1f + rtrn ), zlim3 )
         xlimnpn (ji,jj,jk) = MIN( 1., zlim1)
         !
         ! Michaelis-Menten Limitation term for nutrients picophytoplankton
         ! ----------------------------------------------------------------
         ! Limitation of N based nutrients uptake (NO3 and NH4) 
         zfalim = (1.-fapico) / fapico 
         zlimnh4 = tr(ji,jj,jk,jpnh4,Kbb) / ( zconc0p + tr(ji,jj,jk,jpnh4,Kbb) )
         zlimno3 = tr(ji,jj,jk,jpno3,Kbb) / ( zconc0p + tr(ji,jj,jk,jpno3,Kbb) )
         znutlimtot = (1. - fapico) * ztrn / ( zfalim * zconc0p + ztrn )
         xpiconh4(ji,jj,jk) = znutlimtot * 5.0 * zlimnh4 / ( zlimno3 + 5.0 * zlimnh4 + rtrn )
         xpicono3(ji,jj,jk) = znutlimtot * zlimno3 / ( zlimno3 + 5.0 * zlimnh4 + rtrn )
         !
         ! Limitation of P based nutrients uptake (PO4 and DOP)
         zfalim = (1.-fapicop) / fapicop 
         zlimpo4 = tr(ji,jj,jk,jppo4,Kbb) / ( tr(ji,jj,jk,jppo4,Kbb) + zconc0ppo4 )
         zlimdop = tr(ji,jj,jk,jpdop,Kbb) / ( tr(ji,jj,jk,jpdop,Kbb) + zconc0ppo4 )
         znutlimtot = (1. - fapicop) * ztrp / ( zfalim * zconc0ppo4 + ztrp)
         xpicopo4(ji,jj,jk) = znutlimtot * 100.0 * zlimpo4 / ( zlimdop + 100.0 * zlimpo4 + rtrn )
         xpicodop(ji,jj,jk) = znutlimtot * zlimdop / ( zlimdop + 100.0 * zlimpo4 + rtrn )
         !
         zfalim = (1.-fapicof) / fapicof
         xpicofer(ji,jj,jk) = (1. - fapicof) * zbiron / ( zbiron + zfalim * zconcpfe )
         !
         ! The minimum iron quota depends on the size of PSU, respiration
         ! and the reduction of nitrate following the parameterization 
         ! proposed by Flynn and Hipkin (1999)
         zratiof = tr(ji,jj,jk,jppfe,Kbb) * z1_trnpic
         zqfemp = xcoef1 * zpicochl + xcoef2 + xcoef3 * xpicono3(ji,jj,jk)
         xqfuncfecp(ji,jj,jk) = zqfemp + qfpopt
         !
         zration   = tr(ji,jj,jk,jpnpi,Kbb) * z1_trnpic
         zration = MIN(xqnpmax(ji,jj,jk), MAX( xqnpmin(ji,jj,jk), zration ))
         fvpuptk(ji,jj,jk) = 2.5 * xpsiuptk * xqnpmin(ji,jj,jk) / (zration + rtrn)  &
         &                   * MAX(0., (1. - ratchl * zpicochl / 12. ) ) 
         !
         zlim1    = (zration - xqnpmin(ji,jj,jk) ) / (xqnpmax(ji,jj,jk) - xqnpmin(ji,jj,jk) )

         ! The value of the optimal quota in the formulation below
         ! has been found by solving a non linear equation
         zlim1f   = (1.13 - xqnpmin(ji,jj,jk) ) / (xqnpmax(ji,jj,jk) - xqnpmin(ji,jj,jk) )
         zlim3    = MAX( 0.,( zratiof - zqfemp ) / qfpopt )

         ! computation of the various limitation terms of picophyto
         ! growth and PP
         xlimpfe (ji,jj,jk) = MIN( 1., zlim3 )
         xlimpic (ji,jj,jk) = MIN( 1., zlim1, zlim3 )
         xlimnpp (ji,jj,jk) = MIN( 1., zlim1 )
         xlimpics(ji,jj,jk) = MIN( 1., zlim1/( zlim1f + rtrn ), zlim3 )
         !
         !   Michaelis-Menten Limitation term for nutrients Diatoms
         !   ------------------------------------------------------
         !
         ! Limitation of N based nutrients uptake (NO3 and NH4)
         zfalim = (1.-fadiat) / fadiat 
         zlimnh4 = tr(ji,jj,jk,jpnh4,Kbb) / ( zconc1d + tr(ji,jj,jk,jpnh4,Kbb) )
         zlimno3 = tr(ji,jj,jk,jpno3,Kbb) / ( zconc1d + tr(ji,jj,jk,jpno3,Kbb) )
         znutlimtot = (1.0 - fadiat) * ztrn / ( zfalim * zconc1d + ztrn )
         xdiatnh4(ji,jj,jk) = znutlimtot * 5.0 * zlimnh4 / ( zlimno3 + 5.0 * zlimnh4 + rtrn )
         xdiatno3(ji,jj,jk) = znutlimtot * zlimno3 / ( zlimno3 + 5.0 * zlimnh4 + rtrn )
         !
         ! Limitation of P based nutrients uptake (PO4 and DOP)
         zfalim = (1.-fadiatp) / fadiatp
         zlimpo4 = tr(ji,jj,jk,jppo4,Kbb) / ( tr(ji,jj,jk,jppo4,Kbb) + zconc0dpo4 )
         zlimdop = tr(ji,jj,jk,jpdop,Kbb) / ( tr(ji,jj,jk,jpdop,Kbb) + zconc0dpo4 )
         znutlimtot = (1. - fadiatp) * ztrp / ( zfalim * zconc0dpo4 + ztrp )
         xdiatpo4(ji,jj,jk) = znutlimtot * 100.0 * zlimpo4 / ( zlimdop + 100.0 * zlimpo4 + rtrn )
         xdiatdop(ji,jj,jk) = znutlimtot * zlimdop / ( zlimdop + 100.0 * zlimpo4 + rtrn )
         !
         ! Limitation of Fe uptake
         zfalim = (1.-fadiatf) / fadiatf
         xdiatfer(ji,jj,jk) = (1. - fadiatf) * zbiron / ( zbiron + zfalim * zconcdfe )
         !
         ! The minimum iron quota depends on the size of PSU, respiration
         ! and the reduction of nitrate following the parameterization 
         ! proposed by Flynn and Hipkin (1999)
         zratiof   = tr(ji,jj,jk,jpdfe,Kbb) * z1_trndia
         zqfemd = xcoef1 * zdiatchl + xcoef2 + xcoef3 * xdiatno3(ji,jj,jk)
         xqfuncfecd(ji,jj,jk) = zqfemd + qfdopt
         !
         zration   = tr(ji,jj,jk,jpndi,Kbb) * z1_trndia
         zration   = MIN(xqndmax(ji,jj,jk), MAX( xqndmin(ji,jj,jk), zration ))
         fvduptk(ji,jj,jk) = 2.5 * xpsiuptk * xqndmin(ji,jj,jk) / (zration + rtrn)   &
         &                   * MAX(0., (1. - ratchl * zdiatchl / 12. ) ) 
         !
         zlim1    = (zration - xqndmin(ji,jj,jk) ) / (xqndmax(ji,jj,jk) - xqndmin(ji,jj,jk) )
         ! The value of the optimal quota in the formulation below
         ! has been found by solving a non linear equation
         zlim1f   = (1.13 - xqndmin(ji,jj,jk) ) / (xqndmax(ji,jj,jk) - xqndmin(ji,jj,jk) )
         zlim3    = tr(ji,jj,jk,jpsil,Kbb) / ( tr(ji,jj,jk,jpsil,Kbb) + xksi(ji,jj) )
         zlim4    = MAX( 0., ( zratiof - zqfemd ) / qfdopt )
         ! computation of the various limitation terms of diatoms
         ! growth and PP
         xlimdfe(ji,jj,jk) = MIN( 1., zlim4 )
         xlimdia(ji,jj,jk) = MIN( 1., zlim1, zlim3, zlim4 )
         xlimdias(ji,jj,jk) = MIN (1.0, zlim1 / (zlim1f + rtrn ), zlim3, zlim4 )
         xlimsi(ji,jj,jk)  = MIN( zlim1, zlim4 )
         xlimnpd(ji,jj,jk) = MIN( 1., zlim1 )
         !
         ! Michaelis-Menten Limitation term for nutrients diazotrophs
         ! ----------------------------------------------------------------
         ! Limitation of N based nutrients uptake (NO3 and NH4) 
         zfalim  = (1.-fadiaz) / fadiaz
         zlimnh4 = tr(ji,jj,jk,jpnh4,Kbb) / ( zconc0dz + tr(ji,jj,jk,jpnh4,Kbb) )
         zlimno3 = tr(ji,jj,jk,jpno3,Kbb) / ( zconc0dz + tr(ji,jj,jk,jpno3,Kbb) )
         znutlimtot = (1. - fadiaz) * ztrn  / ( zfalim * zconc0dz + ztrn )
         xdiaznh4(ji,jj,jk) = znutlimtot * 5.0 * zlimnh4 / ( zlimno3 + 5.0 * zlimnh4 + rtrn )
         xdiazno3(ji,jj,jk) = znutlimtot * zlimno3 / ( zlimno3 + 5.0 * zlimnh4 + rtrn )
         !
         !
         ! Limitation of P based nutrients (PO4 and DOP)
         zfalim  = (1.-fadiazp) / fadiazp
         zlimpo4 = tr(ji,jj,jk,jppo4,Kbb) / ( tr(ji,jj,jk,jppo4,Kbb) + zconc0dzpo4 )
         zlimdop = tr(ji,jj,jk,jpdop,Kbb) / ( tr(ji,jj,jk,jpdop,Kbb) + zconc0dzpo4 )
         znutlimtot = (1. - fadiazp) * ztrp / ( zfalim * zconc0dzpo4 + ztrp )
         xdiazpo4(ji,jj,jk) = znutlimtot * 100.0 * zlimpo4 / ( zlimdop + 100.0 * zlimpo4 + rtrn )
         xdiazdop(ji,jj,jk) = znutlimtot * zlimdop / ( zlimdop + 100.0 * zlimpo4 + rtrn )
         !
        ! xdiazpo4(ji,jj,jk) = (1. - fadiazp) * tr(ji,jj,jk,jppo4,Kbb) / ( tr(ji,jj,jk,jppo4,Kbb) + zfalim * zconc0dzpo4 )
        ! xdiazdop(ji,jj,jk) = tr(ji,jj,jk,jpdop,Kbb) / ( tr(ji,jj,jk,jpdop,Kbb)+ xkdoc )   &
        ! &                    * ( 1.0 - xdiazpo4(ji,jj,jk) )

         !
         ! Limitation of Fe uptake
         zfalim = (1.-fadiazf) / fadiazf
         xdiazfer(ji,jj,jk) = (1. - fadiazf) * zbiron / ( zbiron + zfalim * zconcdzfe )
         !
         ! The minimum iron quota depends on the size of PSU, respiration
         ! and the reduction of nitrate following the parameterization 
         ! proposed by Flynn and Hipkin (1999)
         zratiof   = tr(ji,jj,jk,jpfed,Kbb) * z1_trndiaz

         ! Temperature dependence of diazotroph IUE of nitrogen fixation
         !
         IF ( ln_tiue ) THEN
         ! Trichodesmium values based upon Jiang et al. (2018)
           IF ( ln_tricho ) THEN
         !Tricho Temp driven Fe cost scaling set to 1 where Trico growth =0.1 d-1 ~20.5oC
              IUE_scale = (-(0.001392*(ts(ji,jj,jk,jp_tem,Kmm)**5)) + (0.1559*(ts(ji,jj,jk,jp_tem,Kmm)**4)) &
         &                -(6.7685*(ts(ji,jj,jk,jp_tem,Kmm)**3)) + (141.81*(ts(ji,jj,jk,jp_tem,Kmm)**2))   &
         &                -(1421.1*ts(ji,jj,jk,jp_tem,Kmm)) + 5388.1)/33.49
              Fe_scale_nfix(ji,jj,jk) = 1/(IUE_scale+rtrn)
              IF(ts(ji,jj,jk,jp_tem,Kmm) .LT. 17 .OR. ts(ji,jj,jk,jp_tem,Kmm).GT. 35 ) THEN
                Fe_scale_nfix(ji,jj,jk) = 0.
              ENDIF
              IF (Fe_scale_nfix(ji,jj,jk) .GT. maxFescale ) THEN
                Fe_scale_nfix(ji,jj,jk) = maxFescale
              ENDIF
         ! Crocosphaera values based upon Yang et al. (2021)
           ELSE
         !Croco Temp driven Fe cost scaling set to 1 where Croco growth =0.1 d-1 ~21.3oC
              IUE_scale = ((0.02092*(ts(ji,jj,jk,jp_tem,Kmm)**4)) - (2.302*(ts(ji,jj,jk,jp_tem,Kmm)**3)) &
         &                +(92.08*(ts(ji,jj,jk,jp_tem,Kmm)**2)) -(1582*ts(ji,jj,jk,jp_tem,Kmm))   &
         &                 + 9881)/20.64
              Fe_scale_nfix(ji,jj,jk) = 1/(IUE_scale+rtrn)
              IF(ts(ji,jj,jk,jp_tem,Kmm) .LT. 20 .OR. ts(ji,jj,jk,jp_tem,Kmm).GT. 35 ) THEN
                Fe_scale_nfix(ji,jj,jk) = 0.
              ENDIF
              IF ( Fe_scale_nfix(ji,jj,jk) .GT. maxFescale ) THEN
                Fe_scale_nfix(ji,jj,jk) = maxFescale
              ENDIF
           ENDIF

         ELSE
         !No temperature dependance on Fe cost
           Fe_scale_nfix(ji,jj,jk) = 1

         ENDIF

         IF ( ln_facul ) THEN
         facul = (1-(xdiazno3(ji,jj,jk)+xdiaznh4(ji,jj,jk)))  ! Diazotroph Facultative term (1 = all nfix, 0 = no nfix)
         ELSE

         facul = 1.
         xdiaznh4(ji,jj,jk) = 0
         xdiazno3(ji,jj,jk) = 0        

         ENDIF

         qfenfix = kustkaFe * facul * Fe_scale_nfix(ji,jj,jk) ! IUE scalesthe Fe cost of Nfix based upon Kustka et al. (2003)
                                                                    ! kustkaFe = 13 umol/mol
         zqfemdz = xcoef1 * zdiazchl + xcoef2 + xcoef3 * xdiazno3(ji,jj,jk) + qfenfix
        ! zqfemdz = xcoef1 * zdiazchl + xcoef2 + xcoef3 * xdiazno3(ji,jj,jk) + ((xcoef3 * 60) * facul * Fe_scale_nfix(ji,jj,jk))
         qfemindz(ji,jj,jk) = zqfemdz
         xqfuncfecdz(ji,jj,jk) = zqfemdz + qfdzopt
         !
         zration = tr(ji,jj,jk,jpndz,Kbb) * z1_trndiaz
         zration = MIN(xqndzmax(ji,jj,jk), MAX( xqndzmin(ji,jj,jk), zration ))
         fvdzuptk(ji,jj,jk) = 2.5 * xpsiuptk * xqndzmin(ji,jj,jk) / (zration + rtrn)  &
         &                   * MAX(0., (1. - ratchl * zdiazchl / 12. ) )
         !
         zlim1  = (zration - xqndzmin(ji,jj,jk) ) / (xqndzmax(ji,jj,jk) - xqndzmin(ji,jj,jk) )
         
         ! The value of the optimal quota in the formulation below
         ! has been found by solving a non linear equation
         zlim1f = ( 1.13 - xqndzmin(ji,jj,jk) ) / (xqndzmax(ji,jj,jk) - xqndzmin(ji,jj,jk) )
         !
         zlim3  = MAX( 0.,( zratiof - zqfemdz ) / qfdzopt )

         ! computation of the various limitation terms of nanophyto
         ! growth and PP
         xlimdzfe (ji,jj,jk) = MIN( 1., zlim3 )
         xlimdiaz (ji,jj,jk) = MIN( 1., zlim1, zlim3 )
         xlimdiazo(ji,jj,jk) = MIN( 1., zlim1/( zlim1f + rtrn ), zlim3 )
         xlimnpdz (ji,jj,jk) = MIN( 1., zlim1)

      END_3D

      !
      ! Compute the phosphorus quota values. It is based on Litchmann et al., 2004 and Daines et al, 2013.
      ! The relative contribution of three fonctional pools are computed: light harvesting apparatus, 
      ! nutrient uptake pool and assembly machinery. DNA is assumed to represent 1% of the dry mass of 
      ! phytoplankton (see Daines et al., 2013). 
      ! --------------------------------------------------------------------------------------------------
      DO_3D( 0, 0, 0, 0, 1, jpkm1)
         ztrp    = tr(ji,jj,jk,jppo4,Kbb) + tr(ji,jj,jk,jpdop,Kbb) / 200.0
         ! Size estimation of nanophytoplankton based on total biomass
         ! Assumes that larger biomass implies addition of larger cells
         ! ------------------------------------------------------------
         zcoef = tr(ji,jj,jk,jpphy,Kbb) - MIN(xsizephy, tr(ji,jj,jk,jpphy,Kbb) )
         sizena(ji,jj,jk) = 1. + ( xsizern -1.0 ) * zcoef / ( xsizephy + zcoef )
         ! N/P ratio of nanophytoplankton
         ! ------------------------------
         zfuptk = 0.2 + 0.12 / ( 3.0 * sizen(ji,jj,jk) + rtrn )
         ! Computed from Inomura et al. (2020) using Pavlova Lutheri
         zrpho  = 11.55 * tr(ji,jj,jk,jpnch,Kbb) / ( tr(ji,jj,jk,jpphy,Kbb) * 12. + rtrn )
         zrass = 0.62 * (0.15 + 0.85 * ( 1. - zrpho - zfuptk ) * xlimnpn(ji,jj,jk) )
         xqpnmin(ji,jj,jk) = ( 0.0078 + 0.62 * 0.15 * 0.0783 ) * 16.
         xqpnmax(ji,jj,jk) = ( zrpho * 0.0089 + zrass * 0.0783 ) * 16.
         xqpnmax(ji,jj,jk) = xqpnmax(ji,jj,jk) + ( 0.0078 + 0.022 ) * 16. + 3500 * ztrp
         xqpnmax(ji,jj,jk) = MIN( qpnmax, xqpnmax(ji,jj,jk) )

         ! Size estimation of picophytoplankton based on total biomass
         ! Assumes that larger biomass implies addition of larger cells
         ! ------------------------------------------------------------
         zcoef = tr(ji,jj,jk,jppic,Kbb) - MIN(xsizepic, tr(ji,jj,jk,jppic,Kbb) )
         sizepa(ji,jj,jk) = 1. + ( xsizerp -1.0 ) * zcoef / ( xsizepic + zcoef )

         ! N/P ratio of picophytoplankton
         ! ------------------------------
         zfuptk = 0.2 + 0.12 / ( 0.7 * sizep(ji,jj,jk) + rtrn )
         ! Computed from Inomura et al. (2020) using a synechococcus
         zrpho = 13.4 * tr(ji,jj,jk,jppch,Kbb) / ( tr(ji,jj,jk,jppic,Kbb) * 12. + rtrn )
         zrass = 0.4 * ( 0.15 + 0.85 * ( 1. - zrpho - zfuptk ) * xlimnpp(ji,jj,jk) )
         xqppmin(ji,jj,jk) = ( 0.0078 + 0.4/4. * 0.0517 ) * 16.
         xqppmax(ji,jj,jk) = ( zrpho * 0.0076 + zrass * 0.0517 ) * 16.
         xqppmax(ji,jj,jk) = xqppmax(ji,jj,jk) + ( 0.0078 + 0.022 ) * 16. + 1500 * ztrp
         xqppmax(ji,jj,jk) = MIN( qppmax, xqppmax(ji,jj,jk) )

         ! Size estimation of diatoms based on total biomass
         ! Assumes that larger biomass implies addition of larger cells
         ! ------------------------------------------------------------
         zcoef = tr(ji,jj,jk,jpdia,Kbb) - MIN(xsizedia, tr(ji,jj,jk,jpdia,Kbb) )
         sizeda(ji,jj,jk) = 1. + ( xsizerd - 1.0 ) * zcoef / ( xsizedia + zcoef )
         ! N/P ratio of diatoms
         ! --------------------
         zfuptk = 0.2 + 0.12 / ( 5.0 * sized(ji,jj,jk) + rtrn )
         ! Computed from Inomura et al. (2020) using a synechococcus
         zrpho = 8.08 * tr(ji,jj,jk,jpdch,Kbb) / ( tr(ji,jj,jk,jpndi,Kbb) * 12. + rtrn )
         zrass = 0.66 * ( 0.15 + 0.85 * ( 1. - zrpho - zfuptk ) * xlimnpd(ji,jj,jk) )
         xqpdmin(ji,jj,jk) = ( 0.0078 + 0.66/4. * 0.0783 ) * 16.
         xqpdmax(ji,jj,jk) = ( zrpho * 0.0135 + zrass * 0.0783 ) * 16.
         xqpdmax(ji,jj,jk) = xqpdmax(ji,jj,jk) + ( 0.0078 + 0.022 ) * 16. + 5000 * ztrp
         xqpdmax(ji,jj,jk) = MIN(qpdmax, xqpdmax(ji,jj,jk) )

         ! Size estimation of diazotrophys based on total biomass
         ! Assumes that larger biomass implies addition of larger cells
         ! ------------------------------------------------------------
         zcoef = tr(ji,jj,jk,jpcdz,Kbb) - MIN(xsizedz, tr(ji,jj,jk,jpcdz,Kbb) )
         sizedza(ji,jj,jk) = 1. + ( xsizerdz -1.0 ) * zcoef / ( xsizedz + zcoef )
         ! N/P ratio of nanophytoplankton
         ! ------------------------------
         zfuptk = 0.2 + 0.12 / ( 3.0 * sizedz(ji,jj,jk) + rtrn )
         ! Computed from Inomura et al. (2020) using Pavlova Lutheri
         zrpho  = 11.55 * tr(ji,jj,jk,jpchd,Kbb) / ( tr(ji,jj,jk,jpcdz,Kbb) * 12. + rtrn )
         zrass = 0.62 * (0.15 + 0.85 * ( 1. - zrpho - zfuptk ) * xlimnpdz(ji,jj,jk) )
         !xqpdzmin(ji,jj,jk) = ( 0.0078 + 0.62 * 0.15 * 0.0783 ) * 16.
         xqpdzmin(ji,jj,jk) = qpdzmin
         !
         !zrpho  = 1.54 * tr(ji,jj,jk,jpchd,Kbb) / ( tr(ji,jj,jk,jpndz,Kbb) *rno3 * 14. + rtrn )
         !zrass = MAX(0.62/4., ( 1. - zrpho - zfuptk ) *xlimnpdz(ji,jj,jk))
         !xqpdzmin(ji,jj,jk) = ( 0.0 + 0.0078 + 0.62/4. * 0.0783 *xqndzmin(ji,jj,jk) ) * 16.
         ! Temperature dependence of diazotroph PUE of nitrogen fixation
         IF ( ln_tpue ) THEN
         ! Trichodesmium PUE values based upon Jiang et al. (2018)
           IF (ln_tricho ) THEN
           ! Tricho Temp driven QPmin scaling set to 1 where Tricho growth  = 0.1 d-1 20.5oC
              PUE_scale = ((-0.000133*(ts(ji,jj,jk,jp_tem,Kmm)**4)) + (0.012452*(ts(ji,jj,jk,jp_tem,Kmm)**3)) &
              &            -(0.4294*(ts(ji,jj,jk,jp_tem,Kmm)**2)) + (6.538*ts(ji,jj,jk,jp_tem,Kmm))       &
              &            - 37.11)/0.25
              Pmin_scale_nfix(ji,jj,jk) = 1/(PUE_scale+rtrn)
              IF (ts(ji,jj,jk,jp_tem,Kmm) .LT. 17 .OR. ts(ji,jj,jk,jp_tem,Kmm).GT. 35 ) THEN
                 Pmin_scale_nfix(ji,jj,jk) = 0.
              ENDIF
              IF ( Pmin_scale_nfix(ji,jj,jk) .GT. maxPminscale ) THEN
                 Pmin_scale_nfix(ji,jj,jk) = maxPminscale
              ENDIF
           ! Crocosphaera PUE values based upon Yang et al. (2021)
           ELSE
           !Croco Temp driven QPmin scaling set to 1 where Croco growth =0.1 d-1 ~21.3oC       
              PUE_scale = ((-0.0004429*(ts(ji,jj,jk,jp_tem,Kmm)**4)) + (0.04684*(ts(ji,jj,jk,jp_tem,Kmm)**3)) &
              &            -(1.83905*(ts(ji,jj,jk,jp_tem,Kmm)**2)) + (31.829*ts(ji,jj,jk,jp_tem,Kmm))       &
              &            - 204.815)/0.2628
              Pmin_scale_nfix(ji,jj,jk) = 1/(PUE_scale+rtrn)
              IF (ts(ji,jj,jk,jp_tem,Kmm) .LT. 20 .OR. ts(ji,jj,jk,jp_tem,Kmm).GT. 35 ) THEN
                 Pmin_scale_nfix(ji,jj,jk) = 0.
              ENDIF
              IF ( Pmin_scale_nfix(ji,jj,jk) .GT. maxPminscale ) THEN
                 Pmin_scale_nfix(ji,jj,jk) = maxPminscale
              ENDIF
           ENDIF

         ELSE
         ! No temperature dependance on QPmin
           Pmin_scale_nfix(ji,jj,jk) = 1

         ENDIF
         xqpdzmin(ji,jj,jk) = xqpdzmin(ji,jj,jk) * Pmin_scale_nfix(ji,jj,jk)
         xqpdzmax(ji,jj,jk) = ( zrpho * 0.0089 + zrass * 0.0783 ) * 16.
         xqpdzmax(ji,jj,jk) = xqpdzmax(ji,jj,jk) + ( 0.0078 + 0.022 ) * 16. + 3500 * ztrp
         !
         !xqpdzmax(ji,jj,jk) = ( zrpho * 0.0128 + zrass * 0.0783 ) * 16.
         !xqpdzmax(ji,jj,jk) = xqpdzmax(ji,jj,jk) * tr(ji,jj,jk,jpndz,Kbb) / (tr(ji,jj,jk,jpcdz,Kbb) + rtrn )  &
         !&      + (0.033 + 0.0078 ) * 16.
         !
         xqpdzmax(ji,jj,jk) = MIN( qpdzmax, xqpdzmax(ji,jj,jk) )
      END_3D

      ! Compute the fraction of nanophytoplankton that is made of calcifiers
      ! This is a purely adhoc formulation described in Aumont et al. (2015)
      ! This fraction depends on nutrient limitation, light, temperature
      ! --------------------------------------------------------------------
      DO_3D( 0, 0, 0, 0, 1, jpkm1)
         ztem1  = MAX( 0., ts(ji,jj,jk,jp_tem,Kmm) + 1.8 )
         ztem2  = ts(ji,jj,jk,jp_tem,Kmm) - 10.
         zetot1 = MAX( 0., etot_ndcy(ji,jj,jk) - 1.) / ( 4. + etot_ndcy(ji,jj,jk) ) * 30. / ( 30. + etot_ndcy(ji,jj,jk) ) 

         xfracal(ji,jj,jk) = caco3r * xlimphy(ji,jj,jk) * ztem1 / ( 0.1 + ztem1 )     &
            &                * MAX( 1., tr(ji,jj,jk,jpphy,Kbb) / xsizephy )   &
            &                * ( 1. + EXP(-ztem2 * ztem2 / 25. ) )         &
            &                * zetot1 * MIN( 1., 50. / ( hmld(ji,jj) + rtrn ) )
         xfracal(ji,jj,jk) = MAX( 0.02, MIN( 0.8 , xfracal(ji,jj,jk) ) )
      END_3D
      !
      IF( lk_iomput .AND. knt == nrdttrc ) THEN        ! save output diagnostics
        !
        IF( l_dia_fracal ) THEN   ! fraction of calcifiers
          ALLOCATE( zw3d(A2D(0),jpk) )  ;  zw3d(A2D(0),jpk) = 0._wp
          zw3d(A2D(0),1:jpkm1) = xfracal(A2D(0),1:jpkm1) * tmask(A2D(0),1:jpkm1)
          CALL iom_put( "xfracal",  zw3d)
          DEALLOCATE( zw3d )
        ENDIF
        !
        IF( l_dia_nut_lim ) THEN   ! Nutrient limitation term
          ALLOCATE( zw3d(A2D(0),jpk) )  ;  zw3d(A2D(0),jpk) = 0._wp
          zw3d(A2D(0),1:jpkm1) = xlimphy(A2D(0),1:jpkm1) * tmask(A2D(0),1:jpkm1)
          CALL iom_put( "LNnut",  zw3d)
          zw3d(A2D(0),1:jpkm1) = xlimdia(A2D(0),1:jpkm1) * tmask(A2D(0),1:jpkm1)
          CALL iom_put( "LDnut",  zw3d)
          zw3d(A2D(0),1:jpkm1) = xlimpic(A2D(0),1:jpkm1) * tmask(A2D(0),1:jpkm1)
          CALL iom_put( "LPnut",  zw3d)
          zw3d(A2D(0),1:jpkm1) = xlimdiaz(A2D(0),1:jpkm1) * tmask(A2D(0),1:jpkm1)
          CALL iom_put( "LDZnut",  zw3d)
          zw3d(A2D(0),1:jpkm1) = Pmin_scale_nfix(A2D(0),1:jpkm1) * tmask(A2D(0),1:jpkm1)
          CALL iom_put( "Pmin_scale", zw3d)
          zw3d(A2D(0),1:jpkm1) = xlimnpdz(A2D(0),1:jpkm1) * tmask(A2D(0),1:jpkm1)
          CALL iom_put( "LDZNP", zw3d)
          DEALLOCATE( zw3d )
        ENDIF
        !
        IF( l_dia_iron_lim ) THEN   ! Iron limitation term
          ALLOCATE( zw3d(A2D(0),jpk) )  ;  zw3d(A2D(0),jpk) = 0._wp
          zw3d(A2D(0),1:jpkm1) = xlimnfe(A2D(0),1:jpkm1) * tmask(A2D(0),1:jpkm1)
          CALL iom_put( "LNFe",  zw3d)
          zw3d(A2D(0),1:jpkm1) = xlimdfe(A2D(0),1:jpkm1) * tmask(A2D(0),1:jpkm1)
          CALL iom_put( "LDFe",  zw3d)
          zw3d(A2D(0),1:jpkm1) = xlimpfe(A2D(0),1:jpkm1) * tmask(A2D(0),1:jpkm1)
          CALL iom_put( "LPFe",  zw3d)
          zw3d(A2D(0),1:jpkm1) = xlimdzfe(A2D(0),1:jpkm1) * tmask(A2D(0),1:jpkm1)
          CALL iom_put( "LDZFe",  zw3d)
          zw3d(A2D(0),1:jpkm1) = Fe_scale_nfix(A2D(0),1:jpkm1) * tmask(A2D(0),1:jpkm1)
          CALL iom_put( "Fecost_scale", zw3d)
          zw3d(A2D(0),1:jpkm1) = qfemindz(A2D(0),1:jpkm1) * tmask(A2D(0),1:jpkm1)
          CALL iom_put( "MinFeQuotaDz", zw3d)
          DEALLOCATE( zw3d )
        ENDIF
        !
        IF( l_dia_size_lim ) THEN   ! Size limitation term
          ALLOCATE( zw3d(A2D(0),jpk) )  ;  zw3d(A2D(0),jpk) = 0._wp
          zw3d(A2D(0),1:jpkm1) = sizen(A2D(0),1:jpkm1) * tmask(A2D(0),1:jpkm1)
          CALL iom_put( "SIZEN",  zw3d)
          zw3d(A2D(0),1:jpkm1) = sized(A2D(0),1:jpkm1) * tmask(A2D(0),1:jpkm1)
          CALL iom_put( "SIZED",  zw3d)
          zw3d(A2D(0),1:jpkm1) = sizep(A2D(0),1:jpkm1) * tmask(A2D(0),1:jpkm1)
          CALL iom_put( "SIZEP",  zw3d)
          zw3d(A2D(0),1:jpkm1) = sizedz(A2D(0),1:jpkm1) * tmask(A2D(0),1:jpkm1)
          CALL iom_put( "SIZEDZ",  zw3d)
          DEALLOCATE( zw3d )
        ENDIF
        !
        IF( l_dia_size_pro ) THEN   ! Size of the protein machinery
          ALLOCATE( zw3d(A2D(0),jpk) )  ;  zw3d(A2D(0),jpk) = 0._wp
          DO_3D( 0, 0, 0, 0, 1, jpkm1)
             zfuptk = 0.2 + 0.12 / ( 3.0 * sizen(ji,jj,jk) + rtrn )
             zrpho  = 11.55 * tr(ji,jj,jk,jpnch,Kbb) / ( tr(ji,jj,jk,jpphy,Kbb) * 12. + rtrn )
             zw3d(ji,jj,jk) = MAX(0.62/4., ( 1. - zrpho - zfuptk ) * xlimnpn(ji,jj,jk) ) * tmask(ji,jj,jk)
          END_3D
          CALL iom_put( "RASSN",  zw3d)
          DO_3D( 0, 0, 0, 0, 1, jpkm1)
            zfuptk = 0.2 + 0.12 / ( 3.0 * sizep(ji,jj,jk) + rtrn )
             zrpho  = 11.55 * tr(ji,jj,jk,jppch,Kbb) / ( tr(ji,jj,jk,jppic,Kbb) * 12. + rtrn )
             zw3d(ji,jj,jk) = MAX(0.62/4., ( 1. - zrpho - zfuptk ) * xlimnpp(ji,jj,jk) ) * tmask(ji,jj,jk)
          END_3D
          CALL iom_put( "RASSP",  zw3d)
          DO_3D( 0, 0, 0, 0, 1, jpkm1)
             zfuptk = 0.2 + 0.12 / ( 3.0 * sized(ji,jj,jk) + rtrn )
             zrpho  = 11.55 * tr(ji,jj,jk,jpdch,Kbb) / ( tr(ji,jj,jk,jpndi,Kbb) * 12. + rtrn )
             zw3d(ji,jj,jk) = MAX(0.62/4., ( 1. - zrpho - zfuptk ) * xlimnpd(ji,jj,jk) ) * tmask(ji,jj,jk)
          END_3D
          CALL iom_put( "RASSD",  zw3d)
          DO_3D( 0, 0, 0, 0, 1, jpkm1)
             zfuptk = 0.2 + 0.12 / ( 3.0 * sizedz(ji,jj,jk) + rtrn )
             zrpho  = 11.55 * tr(ji,jj,jk,jpchd,Kbb) / ( tr(ji,jj,jk,jpcdz,Kbb) * 12. + rtrn )
             zw3d(ji,jj,jk) = MAX(0.62/4., ( 1. - zrpho - zfuptk ) * xlimnpdz(ji,jj,jk) ) * tmask(ji,jj,jk)
          END_3D
          CALL iom_put( "RASSDZ",  zw3d)
          DEALLOCATE( zw3d )
        ENDIF
        !
        !IF( iom_use( "QFEMINDZ"  ) ) CALL iom_put( "QFEMINDZ"  , qfemindz(:,:,:) * tmask(:,:,:) ) ! min Qfe diazos
        !IF( iom_use( "Pmin_scale") )
        !  ALLOCATE( zw3d(A2D(0),jpk) )  ;  zw3d(A2D(0),jpk) = 0._wp
        !  DO_3D( 0, 0, 0, 0, 1, jpkm1)             zw3d(A2D(0).1:jpkm1)) = Fe_scale_nfix(ji,jj,jk) * tmask(ji,jj,jk)
        !     zw3d(A2D(0),1:jpkm1) = Fe_scale_nfix(A2D(0),1:jpkm1) * tmask(A2D(0),1:jpkm1)
        !  END_3D
        !  CALL iom_put( "Fecost_scale", zw3d)!Fe_scale_nfix(:,:,:) * tmask(:,:,:) ) !Fe cost of nfix Temp dependant scaling
        !  DEALLOCATE( zw3d )
        !ENDIF
        !IF( iom_use( "Pmin_scale") ) CALL iom_put( "Pmin_scale", Pmin_scale_nfix(:,:,:) * tmask(:,:,:) ) !QPmin temperature dependance scaling
      ENDIF
      !
      IF( ln_timing )  CALL timing_stop('p6z_lim')
      !
   END SUBROUTINE p6z_lim


   SUBROUTINE p6z_lim_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p6z_lim_init  ***
      !!
      !! ** Purpose :   Initialization of nutrient limitation parameters
      !!
      !! ** Method  :   Read the namp6zlim and nampisquota namelists and check
      !!      the parameters called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist namp6zlim
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!
      NAMELIST/namp6zlim/ concnno3, concpno3, concdno3, concnnh4, concpnh4, concdnh4,  &
         &                concnfer, concpfer, concdfer, concbfe, concnpo4, concppo4,   &
         &                concdpo4, concbno3, concbnh4, concbpo4, xsizedia, xsizepic,  &
         &                xsizephy, xsizern, xsizerp, xsizerd, xksi1, xksi2, xkdoc,    &
         &                caco3r, oxymin, ratchl, concdzno3, concdznh4, concdzpo4,concdzfer,  &
         &                xsizerdz, xsizedz, Facul_lim, kustkaFe, maxFescale, maxPminscale,   &
         &                ln_tiue, ln_tpue
         !
      NAMELIST/namp6zquota/ qnnmin, qnnmax, qpnmin, qpnmax, qnpmin, qnpmax, qppmin,      &
         &                  qppmax, qndmin, qndmax, qpdmin, qpdmax, qfnmax, qfpmax, qfdmax,  &
         &                  qfnopt, qfpopt, qfdopt, qfdzopt, qndzmin, qndzmax, qpdzmin, qpdzmax, &
         &                  qfdzmax
      !!----------------------------------------------------------------------
      !
      READ  ( numnatp_ref, namp6zlim, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namp6zlim in reference namelist' )
      !
      READ  ( numnatp_cfg, namp6zlim, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'namp6zlim in configuration namelist' )
      IF(lwm) WRITE ( numonp, namp6zlim )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for nutrient limitations, namp6zlim'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    mean rainratio                           caco3r    = ', caco3r
         WRITE(numout,*) '    C associated with Chlorophyll            ratchl    = ', ratchl
         WRITE(numout,*) '    NO3 half saturation of nanophyto         concnno3  = ', concnno3
         WRITE(numout,*) '    NO3 half saturation of picophyto         concpno3  = ', concpno3
         WRITE(numout,*) '    NO3 half saturation of diatoms           concdno3  = ', concdno3
         WRITE(numout,*) '    NH4 half saturation for phyto            concnnh4  = ', concnnh4
         WRITE(numout,*) '    NH4 half saturation for pico             concpnh4  = ', concpnh4
         WRITE(numout,*) '    NH4 half saturation for diatoms          concdnh4  = ', concdnh4
         WRITE(numout,*) '    PO4 half saturation for phyto            concnpo4  = ', concnpo4
         WRITE(numout,*) '    PO4 half saturation for pico             concppo4  = ', concppo4
         WRITE(numout,*) '    PO4 half saturation for diatoms          concdpo4  = ', concdpo4
         WRITE(numout,*) '    half saturation constant for Si uptake   xksi1     = ', xksi1
         WRITE(numout,*) '    half saturation constant for Si/C        xksi2     = ', xksi2
         WRITE(numout,*) '    half-sat. of DOC remineralization        xkdoc     = ', xkdoc
         WRITE(numout,*) '    Iron half saturation for nanophyto       concnfer  = ', concnfer
         WRITE(numout,*) '    Iron half saturation for picophyto       concpfer  = ', concpfer
         WRITE(numout,*) '    Iron half saturation for diatoms         concdfer  = ', concdfer
         WRITE(numout,*) '    size ratio for nanophytoplankton         xsizern   = ', xsizern
         WRITE(numout,*) '    size ratio for picophytoplankton         xsizerp   = ', xsizerp
         WRITE(numout,*) '    size ratio for diatoms                   xsizerd   = ', xsizerd
         WRITE(numout,*) '    NO3 half saturation of bacteria          concbno3  = ', concbno3
         WRITE(numout,*) '    NH4 half saturation for bacteria         concbnh4  = ', concbnh4
         WRITE(numout,*) '    Minimum size criteria for diatoms        xsizedia  = ', xsizedia
         WRITE(numout,*) '    Minimum size criteria for picophyto      xsizepic  = ', xsizepic
         WRITE(numout,*) '    Minimum size criteria for nanophyto      xsizephy  = ', xsizephy
         WRITE(numout,*) '    Fe half saturation for bacteria          concbfe   = ', concbfe
         WRITE(numout,*) '    halk saturation constant for anoxia       oxymin   =' , oxymin
         WRITE(numout,*) '    NO3 half saturation of diazotrophs       concdzno3 = ', concdzno3
         WRITE(numout,*) '    NH4 half saturation of diazotrophs       concdznh4 = ', concdznh4
         WRITE(numout,*) '    PO4 half saturation of diazotrophs       concdzpo4 = ', concdzpo4
         WRITE(numout,*) '    Fe half saturation of diazotrophs        concdzfer = ', concdzfer
         WRITE(numout,*) '    Size ratio for diazotrophs               xsizerdz  = ', xsizerdz
         WRITE(numout,*) '    Minimum size criteria for diazotrophs    xsizedz   = ', xsizedz
         WRITE(numout,*) '    Diazo Facultative limit                  Facul_lim = ', Facul_lim
         WRITE(numout,*) '    Kustka Fe limitation of diazos           kustkaFe  = ', kustkaFe
         WRITE(numout,*) '    Max Fe cost of nfix can be scaled      maxFescale  = ', maxFescale
         WRITE(numout,*) '    Max QPmin can be scaled              maxPminscale  = ', maxPminscale
         WRITE(numout,*) '    turn temp dependence of nfix fe cost     ln_tiue   = ', ln_tiue
         WRITE(numout,*) '    turn temp dependence of nfix QPmin       ln_tpue   = ', ln_tpue
      ENDIF

      READ  ( numnatp_ref, namp6zquota, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisquota in reference namelist' )
      !
      READ  ( numnatp_cfg, namp6zquota, IOSTAT = ios, ERR = 904 )
904   IF( ios >  0 ) CALL ctl_nam ( ios , 'nampisquota in configuration namelist' )
      IF(lwm) WRITE ( numonp, namp6zquota )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for nutrient limitations, namp6zquota'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    optimal Fe quota for nano.               qfnopt    = ', qfnopt
         WRITE(numout,*) '    optimal Fe quota for pico.               qfpopt    = ', qfpopt
         WRITE(numout,*) '    Optimal Fe quota for diatoms             qfdopt    = ', qfdopt
         WRITE(numout,*) '    Minimal N quota for nano                 qnnmin    = ', qnnmin
         WRITE(numout,*) '    Maximal N quota for nano                 qnnmax    = ', qnnmax
         WRITE(numout,*) '    Minimal P quota for nano                 qpnmin    = ', qpnmin
         WRITE(numout,*) '    Maximal P quota for nano                 qpnmax    = ', qpnmax
         WRITE(numout,*) '    Minimal N quota for pico                 qnpmin    = ', qnpmin
         WRITE(numout,*) '    Maximal N quota for pico                 qnpmax    = ', qnpmax
         WRITE(numout,*) '    Minimal P quota for pico                 qppmin    = ', qppmin
         WRITE(numout,*) '    Maximal P quota for pico                 qppmax    = ', qppmax
         WRITE(numout,*) '    Minimal N quota for diatoms              qndmin    = ', qndmin
         WRITE(numout,*) '    Maximal N quota for diatoms              qndmax    = ', qndmax
         WRITE(numout,*) '    Minimal P quota for diatoms              qpdmin    = ', qpdmin
         WRITE(numout,*) '    Maximal P quota for diatoms              qpdmax    = ', qpdmax
         WRITE(numout,*) '    Maximal Fe quota for nanophyto.          qfnmax    = ', qfnmax
         WRITE(numout,*) '    Maximal Fe quota for picophyto.          qfpmax    = ', qfpmax
         WRITE(numout,*) '    Maximal Fe quota for diatoms             qfdmax    = ', qfdmax
         WRITE(numout,*) '    optimal Fe quota for diazotrophs         qfdzopt   = ', qfdzopt
         WRITE(numout,*) '    Minimal N quota for diazotrophs          qndzmin   = ', qndzmin
         WRITE(numout,*) '    Maximal N quota for diazotrophs          qndzmax   = ', qndzmax
         WRITE(numout,*) '    Minimal P quota for diazotrophs          qpdzmin   = ', qpdzmin
         WRITE(numout,*) '    Maximal P quota for diazotrophs          qpdzmax   = ', qpdzmax
         WRITE(numout,*) '    Maximal Fe quota for diazotrophs         qfdzmax   = ', qfdzmax
      ENDIF
      !
      ! Metabolic cost of nitrate and ammonium utilisation
      xpsino3  = 2.3 * rno3
      xpsinh4  = 1.8 * rno3
      xpsiuptk = 1.0 / 6.625
      ! Diazotrophy (biosynthesis on N2)
      xpsinfix = (2.0 / 0.7) * 2.3 * rno3
      !
      xfracal(:,:,jpk) = 0._wp
      xlimphy(:,:,jpk) = 0._wp
      xlimpic(:,:,jpk) = 0._wp
      xlimdia(:,:,jpk) = 0._wp
      xlimnfe(:,:,jpk) = 0._wp
      xlimpfe(:,:,jpk) = 0._wp
      xlimdfe(:,:,jpk) = 0._wp
      sizen  (:,:,jpk) = 0._wp
      sizep  (:,:,jpk) = 0._wp
      sized  (:,:,jpk) = 0._wp
      ! Diazotrophy
      xlimdiaz(:,:,jpk) = 0._wp
      xlimdzfe(:,:,jpk) = 0._wp
      sizedz  (:,:,jpk) = 0._wp
      !
   END SUBROUTINE p6z_lim_init


   INTEGER FUNCTION p6z_lim_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p6z_lim_alloc  ***
      !!----------------------------------------------------------------------
      USE lib_mpp , ONLY: ctl_stop
      INTEGER ::   ierr(2)        ! Local variables
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      !*  Biological arrays for phytoplankton growth
      ALLOCATE( xpicono3(A2D(0),jpk), xpiconh4(A2D(0),jpk),       &
         &      xpicopo4(A2D(0),jpk), xpicodop(A2D(0),jpk),       &
         &      xnanodop(A2D(0),jpk), xdiatdop(A2D(0),jpk),       &
         &      xpicofer(A2D(0),jpk), xlimpfe (A2D(0),jpk),       &
         &      fvnuptk (A2D(0),jpk), fvduptk (A2D(0),jpk),       &
         &      xlimphys(A2D(0),jpk), xlimdias(A2D(0),jpk),       &
         &      xlimnpp (A2D(0),jpk), xlimnpn (A2D(0),jpk),       &
         &      xlimnpd (A2D(0),jpk),                             &
         &      xlimpics(A2D(0),jpk), xqfuncfecp(A2D(0),jpk),     &
         &      fvpuptk (A2D(0),jpk), xlimpic (A2D(0),jpk),       & 
         &      xdiazno3(A2D(0),jpk), xdiaznh4(A2D(0),jpk),       &
         &      xdiazpo4(A2D(0),jpk), xdiazdop(A2D(0),jpk),       &
         &      xdiazfer(A2D(0),jpk), xlimdiaz(A2D(0),jpk),       &
         &      xlimdzfe(A2D(0),jpk), fvdzuptk(A2D(0),jpk),       &
         &      xlimdzp(A2D(0),jpk),  xlimdiazo(A2D(0),jpk),      &
         &      qfemindz(A2D(0),jpk), qfenfixdz(A2D(0),jpk),      &
         &      xlimnpdz(A2D(0),jpk), Fe_scale_nfix(A2D(0),jpk),  &
         &      xqfuncfecdz(A2D(0),jpk),                          &
         &      Pmin_scale_nfix(A2D(0),jpk),       STAT=ierr(1) )
         !
      !*  Minimum/maximum quotas of phytoplankton
      ALLOCATE( xqnnmin (A2D(0),jpk), xqnnmax(A2D(0),jpk),       &
         &      xqpnmin (A2D(0),jpk), xqpnmax(A2D(0),jpk),       &
         &      xqnpmin (A2D(0),jpk), xqnpmax(A2D(0),jpk),       &
         &      xqppmin (A2D(0),jpk), xqppmax(A2D(0),jpk),       &
         &      xqndmin (A2D(0),jpk), xqndmax(A2D(0),jpk),       &
         &      xqpdmin (A2D(0),jpk), xqpdmax(A2D(0),jpk),       &
         &      xqndzmin (A2D(0),jpk), xqndzmax(A2D(0),jpk),     &
         &      xqpdzmin (A2D(0),jpk), xqpdzmax(A2D(0),jpk), STAT=ierr(2) )
         !
      p6z_lim_alloc = MAXVAL( ierr )
      !
      IF( p6z_lim_alloc /= 0 ) CALL ctl_stop( 'STOP', 'p6z_lim_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p6z_lim_alloc
   !!======================================================================
END MODULE p6zlim
