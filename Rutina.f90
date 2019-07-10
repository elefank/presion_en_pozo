module Rutina
!***************************************************************************************************************
!   Propósito:
!    Alojar los procedimientos para llevar a cabo los cálculos de rutina en el cálculo del gradiente de presión
!
!   Bitácora de revisiones:
!       Fecha       Programador     Descripción del cambio
!     ==========    ===========     ======================
!     11/12/2015    Hernán Q. V.    Código Original
!     01/02/2016    ''              Cálulo de las propiedades sin resbalamiento
!***************************************************************************************************************

use variables

contains

    Subroutine Props(Qo_ce,RGA,Diam,Bo,Bg,Rs,Rho_o,Rho_g,Mu_o,Mu_g,VsL,VsG,Vm,Lambda,Mu_n)
    Implicit none
    !Objetivo : Calcular las popiedades de la mezcla necesarias para el cálculo
        Real, Intent(IN) ::Qo_ce,RGA,Diam
        Real, Intent(IN) :: Bo,Bg,Rs,Rho_o,Rho_g,Mu_o,Mu_g
        Real, Intent(OUT) :: VsL,VsG,Vm,Lambda,Mu_n
        Real :: Area,Qo_cf,Qg_cf,Pi
        Pi = 3.141569
    ! Cuerpo
            Area = Pi/4.0*(Diam/12.0)**2
    ! Liquido
            Qo_cf = Qo_ce*Bo*5.614/86400.0  ! Factores de conversión a [ft3/s].
            VsL   = Qo_cf/Area
    !Gas
            Qg_cf = Qo_ce*(RGA-Rs)*Bg*(1/86400.0) ! Factores de conversión a [ft3/s].
            VsG   = Qg_cf/Area

            Vm    = VsL + VsG
            Lambda= VsL/(VsL+vsG)! Qo_cf/(Qg_cf+Qo_cf)

            Mu_n  = Mu_o*Lambda  + Mu_g*(1-Lambda)
            Rho_n = Rho_o*Lambda + Rho_g*(1-Lambda)
    End Subroutine Props

End module Rutina
