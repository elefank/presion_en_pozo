! ***********************************************************************************************************
!
!    This module provides the following functions and constants
!    (1) RadianToDegree()     - Convierte el argumento de radianes a grados
!    (2) DegreeToRadian()     - Convierte el argumento de grados a radianes
!    (3) MySIN()              - Calculael seno del argumento en grados
!
!   Bitácora de revisiones:
!       Fecha       Programador     Descripción del cambio
!     ==========    ===========     ======================
!     02/01/2016    Hernán Q. V.    Código Original

!************************************************************************************************************

MODULE Funciones
use Variables
use Aceite_PVT
use Gas_PVT
use rutina


   IMPLICIT   NONE

   REAL, PARAMETER :: PI        = 3.1415926       ! Declaramos constantes
   REAL, PARAMETER :: Degree180 = 180.0
   REAL, PARAMETER :: R_to_D    = Degree180/PI
   REAL, PARAMETER :: D_to_R    = PI/Degree180

CONTAINS

! --------------------------------------------------------------------
! FUNCION  DegreeToRadian():
!    Esta Funcion toma un argumento REAL en grados y lo convierte a
!    su equivalente en radianes.
! --------------------------------------------------------------------

   REAL FUNCTION  DegreeToRadian(Degree)
      IMPLICIT  NONE
      REAL, INTENT(IN) :: Degree

      DegreeToRadian = Degree * D_to_R
   END FUNCTION  DegreeToRadian

! --------------------------------------------------------------------
! FUNCION  MySIN():
!    Esta funcion toma un argumento REAL en grados y calcula el seno.
!    El cálculo se lleva a cabo en radianes utilizando la función
!    intrínseca de Fortran SIN()y la función definida anteriormente.
! --------------------------------------------------------------------

   REAL FUNCTION  MySIN(x)
      IMPLICIT  NONE
      REAL, INTENT(IN) :: x

      MySIN = SIN(DegreeToRadian(x))
   END FUNCTION  MySIN


!--------------------------------------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------------------------------------
    Subroutine PVT(Presion,Temperatura,Sgo,SGg,RGA,Bo,Bg,Rs,Rho_o,Rho_g,Mu_o,Mu_g,Sigma_og,ask3)
     Real, Intent(IN)  ::Presion,Temperatura,SGo,SGg,RGA
     Real, Intent(OUT) :: Bo,Bg,Rs,Rho_o,Rho_g,Mu_o,Mu_g,Sigma_og
     Character,intent(IN)::ask3
        !Vazquez(API,Temperatura,Presion,SGg,SGo,RsV,PbV,SGcien,Psep,Tsep)
     Call Vazquez(API,Temperatura,Presion,SGg,SGo,Rs,Pb,SGcien,Psep,Tsep,ask3)

        !VazquezBo(API,SGo,Temperatura,SGcien,RsV,BoV)
     Call VazquezBo(API,SGo,Temperatura,SGcien,Rs,Bo)

        !Standing_SGd(API,SGo,Rr,SGd)
     Call Standing_SGd(API,SGo,Rs,SGd)

        !DensityOilCY (SGo,SGd,Rs,Bo,Do_cy)
     Call DensityOilCY (SGo,SGd,Rs,Bo,Rho_o)

        !Oil_viscosity(Sgo,API,P,T,Pb,Rs,Mu_ob)
     Call Oil_viscosity(Sgo,API,Presion,Temperatura,Pb,Rs,Mu_o)

!------------------------------------------------------------------------------------------------------------
!                               Propiedades del GAS

        !SG_Free_Gas(RGa,Rs,SGd,SGt,SGf,Criterio)
     Call SG_Free_Gas(RGA,Rs,SGd,SGg,SGf,Criterio)

        !GasZ_DAK(Presion,Temperatura,SGg,Z)
     Call GasZ_DAK(Presion,Temperatura,SGg,Z)

        !GasBg(Z,P,T,Bg)
     Call GasBg(Z,Presion,Temperatura,Bg)

        !GasDensity(SGf,P,T,Z,Dgas_cy)
     Call GasDensity(SGf,Presion,Temperatura,Z,Rho_g)

        !Gasviscosity(SGf,T,Dgas_cy,Mu_g)
     Call Gasviscosity(SGf,Temperatura,Rho_g,Mu_g)

        !Oil_Gas_Interfacial_Tension(Rs,SGo,API,T,Sigma_og)
     Call Oil_Gas_Interfacial_Tension(Rs,SGo,API,Temperatura,Sigma_og)

End Subroutine PVT

END MODULE  Funciones
