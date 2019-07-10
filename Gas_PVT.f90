Module Gas_PVT
!************************************************************************************************************
!   Propósito:
!    Este módulo contiene las subrutinas que permiten calcular las porpiedades PVT para gas,
!    utilizando diversas correlaciones correlaciones.
!
!   Observaciones:
!     En esta parte del código solo se hacen llamadas a archivos para entrada y salida de datos y para llamar
!     subrutinas.
!
!   Bitácora de revisiones:
!       Fecha       Programador     Descripción del cambio
!     ==========    ===========     ======================
!     01/02/2016    Hernán Q. V.    Código Original
!
!************************************************************************************************************
use variables

implicit none

contains

!------------------------------------------------------------------------------------------------------------
!   Subrutina GasZ_DAK:
!       Calcula el factor de desviacion del gas por el metodo de Dranchuk-Abu-Kasem,
!       utilizando el metodo de Newton Rhapson.
!------------------------------------------------------------------------------------------------------------
        Subroutine GasZ_DAK(Presion,Temperatura,SGg,Z)
        Implicit None
        Real(kind=4),Intent(IN)::Presion,Temperatura,SGg
        Real(kind=4),Intent(OUT)::Z
        Real(kind=4)::Psc,Tsc,Psr,Tsr,dr,z2,Fz,dFz
        Real(kind=4)::A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11
        Real(kind=4)::C1,C2,C3,C4
        Integer::i

            A1 =  0.3265
            A2 = -1.07
            A3 = -0.5339
            A4 =  0.01569
            A5 = -0.05165
            A6 =  0.5475
            A7 = -0.7361
            A8 =  0.1844
            A9 =  0.1056
           A10 =  0.6134
           A11 =  0.7210

            Psc = 756.8-131*SGg-3.6*SGg**2
            Tsc = 169.2+349.5*SGg-74*SGg**2

            Psr = Presion/Psc
            Tsr = (Temperatura+460)/Tsc

            z2 = 0.1
            dr = 0.27*(Psr/(z2*Tsr))

            C1 = A1+A2/Tsr+A3/Tsr**3+A4/Tsr**4+A5/Tsr**5
            C2 = A6+A7/Tsr+A8/Tsr**2
            C3 = A9*(A7/Tsr+A8/Tsr**2)
            C4 = A10*(1+A11*dr**2)*(dr**2/Tsr**3)*exp(-A11*dr**2)

            Fz = z2-( 1 + C1*dr + C2*dr**2 - C3*dr**5 + C4 )
            i = 1

            Do While (ABS(Fz).GT.0.00001)
                dFz = 1 + C1*dr/z2 + 2*C2*dr**2/z2 - 5*C3*dr**5/z2 + 2*A10*dr**2/(Tsr**3*z2)*&
                (1 + A11*dr**2 - (A11*dr**2)**2)*exp(-A11*dr**2)

                z2 = z2-(Fz/dFz)
                dr = 0.27*(Psr/(z2*Tsr))
                Fz = z2-( 1 + C1*dr + C2*dr**2 - C3*dr**5 + C4 )
                i = i+1
            End Do

            Z = z2
   End Subroutine GasZ_DAK


!------------------------------------------------------------------------------------------------------------
!   Subrutina GasDensity:
!       Calcula la densidad del gas
!------------------------------------------------------------------------------------------------------------
        Subroutine GasDensity(SGf,P,T,Z,Dgas_cy)
        Implicit None
        Real(kind=4),INTENT(IN)::SGf,P,T,Z
        Real(kind=4),INTENT(OUT)::Dgas_cy

           Dgas_cy = (2.7*SGf*P)/(Z*(T+460))
!          Dgas_cy = (SGf*0.0764)*Bg
        End Subroutine Gasdensity


!------------------------------------------------------------------------------------------------------------
!   Subrutina GasBg:
!       Calcula el factor de volumen del gas
!------------------------------------------------------------------------------------------------------------
        Subroutine GasBg(Z,P,T,Bg)
        IMPLICIT NONE
        Real(kind=4),INTENT(IN)::Z,P,T
        Real(kind=4),INTENT(OUT)::Bg

            Bg = 0.0283*Z*(T+460)/P
        End Subroutine GasBg


!------------------------------------------------------------------------------------------------------------
!   Subrutina Gasviscosity:
!       Calcula la viscosidad del gas q condicionesde yacimiento.
!------------------------------------------------------------------------------------------------------------
        Subroutine Gasviscosity(SGf,T,Dgas_cy,Mu_g)
        IMPLICIT NONE
        Real(kind=4),INTENT(IN)::SGf,T,Dgas_cy
        Real(kind=4),INTENT(OUT)::Mu_g
        Real(kind=4)::Mg,X,Y,K

            Mg = 28.97*SGf
            !Write(*,*)Mg
            X = 3.5+(986/(T+460))+0.01*Mg
            !Write(*,*)X
            Y = 2.4-0.2*X
            !Write(*,*)Y
            K = ((9.4+0.02*Mg)*((T+460))**1.5)/(209+19*Mg+(T+460))
            !Write(*,*)K
            Mu_g = 10.0**(-4) * K *exp( X* (Dgas_cy/62.4)**Y)
            !Write(*,*)Mu_g
        End Subroutine Gasviscosity


!--------------------------------------------------------------------------------------------------
!       Calcula la densidad relativa del Gas libre
!--------------------------------------------------------------------------------------------------
        Subroutine SG_Free_Gas(RGa,Rs,SGd,SGt,SGf,Criterio)
        Implicit none
        Real(kind=4),intent(IN)::RGa,Rs,SGt
        Real(kind=4),intent(INOUT)::SGd
        Real(kind=4),intent(OUT)::SGf,Criterio

        SGf = (RGA*SGt-Rs*SGd)/(RGA-Rs)

        IF ( (SGt.LE.SGd).and.(SGd.GE.0.56).and.(Sgt.GE.SGf).and.(SGf.GE.0.56) ) then
!            write(9,6)
!            write(*,6)
!            6 Format('',T14,'La composición no es constante',/)
        ELse
!            write(9,5)
!            write(*,5)
!            5 Format('',T14,'Se asume composición constante',/)
           Criterio = 101.d0

        End IF

        IF (Criterio .EQ. 101.d0) then
            SGd = SGg
            SGf = SGg
        END If

        End subroutine SG_Free_Gas


End module Gas_PVT
