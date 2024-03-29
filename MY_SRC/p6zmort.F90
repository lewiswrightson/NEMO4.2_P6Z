MODULE p6zmort
   !!======================================================================
   !!                         ***  MODULE p6zmort  ***
   !! TOP :   PISCES-QUOTA Compute the mortality terms for phytoplankton
   !!         including explicit diazotrophy
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont)  Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
   !!----------------------------------------------------------------------
   !!   p6z_mort       :   Compute the mortality terms for phytoplankton
   !!   p6z_mort_init  :   Initialize the mortality params for phytoplankton
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p2zlim
   USE p4zlim          !  Phytoplankton limitation terms (p4z)
   USE p6zlim          !  Phytoplankton limitation terms (p6z)
   USE prtctl          !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p6z_mort           ! Called from p4zbio.F90 
   PUBLIC   p6z_mort_init      ! Called from trcini_pisces.F90 

   !! * Shared module variables
   REAL(wp), PUBLIC :: wchln   !! Quadratic mortality rate of nanophytoplankton
   REAL(wp), PUBLIC :: wchlp   !: Quadratic mortality rate of picophytoplankton
   REAL(wp), PUBLIC :: wchld   !: Quadratic mortality rate of diatoms
   REAL(wp), PUBLIC :: mpratn  !: Linear mortality rate of nanophytoplankton
   REAL(wp), PUBLIC :: mpratp  !: Linear mortality rate of picophytoplankton
   REAL(wp), PUBLIC :: mpratd  !: Linear mortality rate of diatoms
   REAL(wp), PUBLIC :: wchldz  !: Quadratic mortality rate of diazotrophs
   REAL(wp), PUBLIC :: mpratdz !: Linear mortality rate of diazotrophs

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p6zmort.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p6z_mort( kt, Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p6z_mort  ***
      !!
      !! ** Purpose :   Calls the different subroutine to compute
      !!                the different phytoplankton mortality terms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Krhs  ! time level indices
      !!---------------------------------------------------------------------

      CALL p6z_mort_nano( Kbb, Krhs )            ! nanophytoplankton
      CALL p6z_mort_pico( Kbb, Krhs )            ! picophytoplankton
      CALL p6z_mort_diat( Kbb, Krhs )            ! diatoms
      CALL p6z_mort_diazo( Kbb, Krhs )           ! diazotrophs

   END SUBROUTINE p6z_mort


   SUBROUTINE p6z_mort_nano( Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p6z_mort_nano  ***
      !!
      !! ** Purpose :   Compute the mortality terms for nanophytoplankton
      !!
      !! ** Method  : - Both quadratic and simili linear mortality terms
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kbb, Krhs  ! time level indices
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompaph, zlim1, zlim2
      REAL(wp) :: zfactfe, zfactch, zfactn, zfactp, zprcaca
      REAL(wp) :: ztortp , zrespp , zmortp
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p6z_mort_nano')
      !
      prodcal(:,:,:) = 0._wp   ! calcite production variable set to zero
      DO_3D( 0, 0, 0, 0, 1, jpkm1)
         zcompaph = MAX( ( tr(ji,jj,jk,jpphy,Kbb) - 1e-9 ), 0.e0 )

         ! Quadratic mortality of nano due to aggregation during
         ! blooms (Doney et al. 1996)
         ! -----------------------------------------------------
         zlim2   = xlimphy(ji,jj,jk) * xlimphy(ji,jj,jk)
         zlim1   = 0.25 * ( 1. - zlim2 ) / ( 0.25 + zlim2 ) * tr(ji,jj,jk,jpphy,Kbb)
         zrespp = wchln * 1.e6 * xstep * zlim1 * xdiss(ji,jj,jk) * zcompaph

         ! Phytoplankton linear mortality
         ! A michaelis-menten like term is introduced to avoid 
         ! extinction of nanophyto in highly limited areas
         ! ----------------------------------------------------
         ztortp = mpratn * xstep * zcompaph * tr(ji,jj,jk,jpphy,Kbb) / ( xkmort + tr(ji,jj,jk,jpphy,Kbb) )
         zmortp = zrespp + ztortp

         !   Update the arrays TRA which contains the biological sources and sinks
         zfactn  = tr(ji,jj,jk,jpnph,Kbb)/(tr(ji,jj,jk,jpphy,Kbb)+rtrn)
         zfactp  = tr(ji,jj,jk,jppph,Kbb)/(tr(ji,jj,jk,jpphy,Kbb)+rtrn)
         zfactfe = tr(ji,jj,jk,jpnfe,Kbb)/(tr(ji,jj,jk,jpphy,Kbb)+rtrn)
         zfactch = tr(ji,jj,jk,jpnch,Kbb)/(tr(ji,jj,jk,jpphy,Kbb)+rtrn)
         tr(ji,jj,jk,jpphy,Krhs) = tr(ji,jj,jk,jpphy,Krhs) - zmortp
         tr(ji,jj,jk,jpnph,Krhs) = tr(ji,jj,jk,jpnph,Krhs) - zmortp * zfactn
         tr(ji,jj,jk,jppph,Krhs) = tr(ji,jj,jk,jppph,Krhs) - zmortp * zfactp
         tr(ji,jj,jk,jpnch,Krhs) = tr(ji,jj,jk,jpnch,Krhs) - zmortp * zfactch
         tr(ji,jj,jk,jpnfe,Krhs) = tr(ji,jj,jk,jpnfe,Krhs) - zmortp * zfactfe

                       ! Production PIC particles due to mortality
         zprcaca = xfracal(ji,jj,jk) * zmortp
         prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
         !
         tr(ji,jj,jk,jpdic,Krhs) = tr(ji,jj,jk,jpdic,Krhs) - zprcaca
         tr(ji,jj,jk,jptal,Krhs) = tr(ji,jj,jk,jptal,Krhs) - 2. * zprcaca
         tr(ji,jj,jk,jpcal,Krhs) = tr(ji,jj,jk,jpcal,Krhs) + zprcaca
         tr(ji,jj,jk,jppoc,Krhs) = tr(ji,jj,jk,jppoc,Krhs) + zmortp
         tr(ji,jj,jk,jppon,Krhs) = tr(ji,jj,jk,jppon,Krhs) + zmortp * zfactn
         tr(ji,jj,jk,jppop,Krhs) = tr(ji,jj,jk,jppop,Krhs) + zmortp * zfactp
         prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + zmortp
         tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) + zmortp * zfactfe
      END_3D
      !
       IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('nano')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
         CALL prt_ctl(tab4d_1=tr(:,:,:,:,Krhs), mask1=tmask, clinfo=ctrcnm)
       ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p6z_mort_nano')
      !
   END SUBROUTINE p6z_mort_nano


   SUBROUTINE p6z_mort_pico( Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p6z_mort_pico  ***
      !!
      !! ** Purpose :   Compute the mortality terms for picophytoplankton
      !!
      !! ** Method  : - Both quadratic and semilininear terms are used
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kbb, Krhs  ! time level indices
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompaph, zlim1, zlim2
      REAL(wp) :: zfactfe, zfactch, zfactn, zfactp
      REAL(wp) :: ztortp , zrespp , zmortp 
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p6z_mort_pico')
      !
      DO_3D( 0, 0, 0, 0, 1, jpkm1)
         zcompaph = MAX( ( tr(ji,jj,jk,jppic,Kbb) - 1e-9 ), 0.e0 )

         ! Quadratic mortality of pico due to aggregation during
         ! blooms (Doney et al. 1996)
         ! -----------------------------------------------------
         zlim2   = xlimpic(ji,jj,jk) * xlimpic(ji,jj,jk)
         zlim1   = 0.25 * ( 1. - zlim2 ) / ( 0.25 + zlim2 ) * tr(ji,jj,jk,jppic,Kbb)
         zrespp = wchlp * 1.e6 * xstep * zlim1 * xdiss(ji,jj,jk) * zcompaph

         ! Phytoplankton linear mortality
         ! A michaelis-menten like term is introduced to avoid 
         ! extinction of picophyto in highly limited areas
         ! ----------------------------------------------------
         ztortp = mpratp * xstep  * zcompaph * tr(ji,jj,jk,jppic,Kbb) /  ( xkmort + tr(ji,jj,jk,jppic,Kbb) )
         zmortp = zrespp + ztortp

         !   Update the arrays TRA which contains the biological sources and sinks
         zfactn = tr(ji,jj,jk,jpnpi,Kbb)/(tr(ji,jj,jk,jppic,Kbb)+rtrn)
         zfactp = tr(ji,jj,jk,jpppi,Kbb)/(tr(ji,jj,jk,jppic,Kbb)+rtrn)
         zfactfe = tr(ji,jj,jk,jppfe,Kbb)/(tr(ji,jj,jk,jppic,Kbb)+rtrn)
         zfactch = tr(ji,jj,jk,jppch,Kbb)/(tr(ji,jj,jk,jppic,Kbb)+rtrn)
         tr(ji,jj,jk,jppic,Krhs) = tr(ji,jj,jk,jppic,Krhs) - zmortp
         tr(ji,jj,jk,jpnpi,Krhs) = tr(ji,jj,jk,jpnpi,Krhs) - zmortp * zfactn
         tr(ji,jj,jk,jpppi,Krhs) = tr(ji,jj,jk,jpppi,Krhs) - zmortp * zfactp
         tr(ji,jj,jk,jppch,Krhs) = tr(ji,jj,jk,jppch,Krhs) - zmortp * zfactch
         tr(ji,jj,jk,jppfe,Krhs) = tr(ji,jj,jk,jppfe,Krhs) - zmortp * zfactfe
         tr(ji,jj,jk,jppoc,Krhs) = tr(ji,jj,jk,jppoc,Krhs) + zmortp
         tr(ji,jj,jk,jppon,Krhs) = tr(ji,jj,jk,jppon,Krhs) + zmortp * zfactn
         tr(ji,jj,jk,jppop,Krhs) = tr(ji,jj,jk,jppop,Krhs) + zmortp * zfactp
         tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) + zmortp * zfactfe
         prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + zmortp
      END_3D
      !
       IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('pico')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
         CALL prt_ctl(tab4d_1=tr(:,:,:,:,Krhs), mask1=tmask, clinfo=ctrcnm)
       ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p6z_mort_pico')
      !
   END SUBROUTINE p6z_mort_pico


   SUBROUTINE p6z_mort_diat( Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p6z_mort_diat  ***
      !!
      !! ** Purpose :   Compute the mortality terms for diatoms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kbb, Krhs  ! time level indices
      INTEGER  ::  ji, jj, jk
      REAL(wp) ::  zfactfe,zfactsi,zfactch, zfactn, zfactp, zcompadi
      REAL(wp) ::  zrespp, ztortp, zmortp
      REAL(wp) ::  zlim2, zlim1
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p6z_mort_diat')
      !

      DO_3D( 0, 0, 0, 0, 1, jpkm1)

         zcompadi = MAX( ( tr(ji,jj,jk,jpdia,Kbb) - 1E-9), 0. )

         !   Aggregation term for diatoms is increased in case of nutrient
         !   stress as observed in reality. The stressed cells become more
         !   sticky and coagulate to sink quickly out of the euphotic zone
         !   -------------------------------------------------------------
         !  Phytoplankton squared mortality
         !  -------------------------------
         zlim2   = xlimdia(ji,jj,jk) * xlimdia(ji,jj,jk)
         zlim1   = 0.25 * ( 1. - zlim2 ) / ( 0.25 + zlim2 ) 
         zrespp  = 1.e6 * xstep * wchld * zlim1 * xdiss(ji,jj,jk) * zcompadi * tr(ji,jj,jk,jpdia,Kbb)

         ! Phytoplankton linear mortality
         ! A michaelis-menten like term is introduced to avoid 
         ! extinction of diatoms in highly limited areas
         !  ---------------------------------------------------
         ztortp  = mpratd * xstep  * zcompadi * tr(ji,jj,jk,jpdia,Kbb) /  ( xkmort + tr(ji,jj,jk,jpdia,Kbb) )
         zmortp  = zrespp + ztortp

         !   Update the arrays tr(:,:,:,:,Krhs) which contains the biological sources and sinks
         !   ---------------------------------------------------------------------
         zfactn  = tr(ji,jj,jk,jpndi,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )
         zfactp  = tr(ji,jj,jk,jppdi,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )
         zfactch = tr(ji,jj,jk,jpdch,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )
         zfactfe = tr(ji,jj,jk,jpdfe,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )
         zfactsi = tr(ji,jj,jk,jpdsi,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )
         tr(ji,jj,jk,jpdia,Krhs) = tr(ji,jj,jk,jpdia,Krhs) - zmortp 
         tr(ji,jj,jk,jpndi,Krhs) = tr(ji,jj,jk,jpndi,Krhs) - zmortp * zfactn
         tr(ji,jj,jk,jppdi,Krhs) = tr(ji,jj,jk,jppdi,Krhs) - zmortp * zfactp
         tr(ji,jj,jk,jpdch,Krhs) = tr(ji,jj,jk,jpdch,Krhs) - zmortp * zfactch
         tr(ji,jj,jk,jpdfe,Krhs) = tr(ji,jj,jk,jpdfe,Krhs) - zmortp * zfactfe
         tr(ji,jj,jk,jpdsi,Krhs) = tr(ji,jj,jk,jpdsi,Krhs) - zmortp * zfactsi
         tr(ji,jj,jk,jpgsi,Krhs) = tr(ji,jj,jk,jpgsi,Krhs) + zmortp * zfactsi
         tr(ji,jj,jk,jpgoc,Krhs) = tr(ji,jj,jk,jpgoc,Krhs) + zrespp 
         tr(ji,jj,jk,jpgon,Krhs) = tr(ji,jj,jk,jpgon,Krhs) + zrespp * zfactn
         tr(ji,jj,jk,jpgop,Krhs) = tr(ji,jj,jk,jpgop,Krhs) + zrespp * zfactp
         tr(ji,jj,jk,jpbfe,Krhs) = tr(ji,jj,jk,jpbfe,Krhs) + zrespp * zfactfe
         tr(ji,jj,jk,jppoc,Krhs) = tr(ji,jj,jk,jppoc,Krhs) + ztortp
         tr(ji,jj,jk,jppon,Krhs) = tr(ji,jj,jk,jppon,Krhs) + ztortp * zfactn
         tr(ji,jj,jk,jppop,Krhs) = tr(ji,jj,jk,jppop,Krhs) + ztortp * zfactp
         tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) + ztortp * zfactfe
         prodpoc(ji,jj,jk)   = prodpoc(ji,jj,jk) + ztortp
         prodgoc(ji,jj,jk)   = prodgoc(ji,jj,jk) + zrespp

      END_3D
      !
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('diat')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
         CALL prt_ctl(tab4d_1=tr(:,:,:,:,Krhs), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p6z_mort_diat')
      !
   END SUBROUTINE p6z_mort_diat

   SUBROUTINE p6z_mort_diazo( Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p6z_mort_diazo  ***
      !!
      !! ** Purpose :   Compute the mortality terms for diazotrophs
      !!
      !! ** Method  : - Both quadratic and simili linear mortality terms
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kbb, Krhs  ! time level indices
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompaph, zlim1, zlim2
      REAL(wp) :: zfactfe, zfactch, zfactn, zfactp, zprcaca
      REAL(wp) :: ztortp , zrespp , zmortp
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p6z_mort_diazo')
      !
      prodcal(:,:,:) = 0._wp   ! calcite production variable set to zero
      DO_3D( 0, 0, 0, 0, 1, jpkm1)
         zcompaph = MAX( ( tr(ji,jj,jk,jpcdz,Kbb) - 1e-9 ), 0.e0 )

         ! Quadratic mortality of diatom due to aggregation during
         ! blooms (Doney et al. 1996)
         ! -----------------------------------------------------
         zlim2   = xlimdiaz(ji,jj,jk) * xlimdiaz(ji,jj,jk)
         zlim1   = 0.25 * ( 1. - zlim2 ) / ( 0.25 + zlim2 ) * tr(ji,jj,jk,jpcdz,Kbb)
         zrespp = wchldz * 1.e6 * xstep * zlim1 * xdiss(ji,jj,jk) * zcompaph

         ! Phytoplankton linear mortality
         ! A michaelis-menten like term is introduced to avoid 
         ! extinction of nanophyto in highly limited areas
         ! ----------------------------------------------------
         ztortp = mpratdz * xstep * zcompaph * tr(ji,jj,jk,jpcdz,Kbb) / ( xkmort + tr(ji,jj,jk,jpcdz,Kbb) )
         zmortp = zrespp + ztortp

         !   Update the arrays TRA which contains the biological sources and
         !   sinks
         zfactn  = tr(ji,jj,jk,jpndz,Kbb)/(tr(ji,jj,jk,jpcdz,Kbb)+rtrn)
         zfactp  = tr(ji,jj,jk,jppdz,Kbb)/(tr(ji,jj,jk,jpcdz,Kbb)+rtrn)
         zfactfe = tr(ji,jj,jk,jpfed,Kbb)/(tr(ji,jj,jk,jpcdz,Kbb)+rtrn)
         zfactch = tr(ji,jj,jk,jpchd,Kbb)/(tr(ji,jj,jk,jpcdz,Kbb)+rtrn)
         tr(ji,jj,jk,jpcdz,Krhs) = tr(ji,jj,jk,jpcdz,Krhs) - zmortp
         tr(ji,jj,jk,jpndz,Krhs) = tr(ji,jj,jk,jpndz,Krhs) - zmortp * zfactn
         tr(ji,jj,jk,jppdz,Krhs) = tr(ji,jj,jk,jppdz,Krhs) - zmortp * zfactp
         tr(ji,jj,jk,jpchd,Krhs) = tr(ji,jj,jk,jpchd,Krhs) - zmortp * zfactch
         tr(ji,jj,jk,jpfed,Krhs) = tr(ji,jj,jk,jpfed,Krhs) - zmortp * zfactfe
         !
         tr(ji,jj,jk,jppoc,Krhs) = tr(ji,jj,jk,jppoc,Krhs) + ztortp
         tr(ji,jj,jk,jppon,Krhs) = tr(ji,jj,jk,jppon,Krhs) + ztortp * zfactn
         tr(ji,jj,jk,jppop,Krhs) = tr(ji,jj,jk,jppop,Krhs) + ztortp * zfactp
         !
         tr(ji,jj,jk,jpgoc,Krhs) = tr(ji,jj,jk,jpgoc,Krhs) + zrespp
         tr(ji,jj,jk,jpgon,Krhs) = tr(ji,jj,jk,jpgon,Krhs) + zrespp * zfactn
         tr(ji,jj,jk,jpgop,Krhs) = tr(ji,jj,jk,jpgop,Krhs) + zrespp * zfactp
         tr(ji,jj,jk,jpbfe,Krhs) = tr(ji,jj,jk,jpbfe,Krhs) + zrespp * zfactfe
         prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + zmortp
         prodgoc(ji,jj,jk)   = prodgoc(ji,jj,jk) + zrespp
         tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) + ztortp * zfactfe
      END_3D
      !
       IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('diazo')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
         CALL prt_ctl(tab4d_1=tr(:,:,:,:,Krhs), mask1=tmask, clinfo=ctrcnm)
       ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p6z_mort_diazo')
      !
   END SUBROUTINE p6z_mort_diazo

   SUBROUTINE p6z_mort_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p6z_mort_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton mortality parameters
      !!
      !! ** Method  :   Read the namp6zmort namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist namp6zmort
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ios   ! Local integer output status for namelist read
      !!
      NAMELIST/namp6zmort/ wchln, wchlp, wchld, mpratn, mpratp, mpratd, wchldz, mpratdz
      !!----------------------------------------------------------------------

      READ  ( numnatp_ref, namp6zmort, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namp6zmort in reference namelist' )

      READ  ( numnatp_cfg, namp6zmort, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'namp6zmort in configuration namelist' )
      IF(lwm) WRITE ( numonp, namp6zmort )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for phytoplankton mortality, namp6zmort'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    quadratic mortality of phytoplankton      wchln     =', wchln
         WRITE(numout,*) '    quadratic mortality of picophyto.         wchlp     =', wchlp
         WRITE(numout,*) '    quadratic mortality of diatoms            wchld     =', wchld
         WRITE(numout,*) '    nanophyto. mortality rate                 mpratn    =', mpratn
         WRITE(numout,*) '    picophyto. mortality rate                 mpratp    =', mpratp
         WRITE(numout,*) '    Diatoms mortality rate                    mpratd    =', mpratd
         WRITE(numout,*) '    quadratic mortality of diazotrophs        wchldz    =', wchldz
         WRITE(numout,*) '    diazotrophs mortality rate                mpratdz   =', mpratdz
      ENDIF

   END SUBROUTINE p6z_mort_init

   !!======================================================================
END MODULE p6zmort
