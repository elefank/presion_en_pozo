module Aceite_PVT
!************************************************************************************************************
!   Propósito:
!    Este módulo contiene las subrutinas que permiten calcular las porpiedades PVT para aceite
!    negro, utilizando diversas correlaciones correlaciones.
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
contains

!------------------------------------------------------------------------------------------------------------
!   Subrutina Vazquez:
!       Calcula la presión de burbuja y la relación de solubilidad, utilizando las correlaciones de Vazquez
!       pregunta al usuario por las condiciones de separación si estan disponibles, de lo contrario, utiliza
!       la gravedad específica del gas total producido.
!------------------------------------------------------------------------------------------------------------
        Subroutine Vazquez(API,Temperatura,Presion,SGg,SGo,RsV,PbV,SGcien,Psep,Tsep,ask3)
         IMPLICIT NONE
         Real, External    :: API
         Real,INTENT(IN)   :: Temperatura,Presion,SGg,SGo
         Real,Intent(IN):: Psep, Tsep
         Real,INTENT(OUT)  :: RsV,PbV,SGcien
         Real::c1,c2,c3
         Character,intent(IN)::ask3

            IF (ask3 =='y'.OR. ask3 == 'Y') Then

                SGcien = SGg *( 1.0 + 5.912E-5 * API(SGo)* Tsep * LOG10(Psep/114.7) )
            Else
                SGcien = SGg
            End IF


            IF (API(SGo).LE.30.0)Then
                c1=0.0362
                c2=1.0937
                c3=25.724
            Else
                c1=0.0178
                c2=1.1870
                c3=23.931
            End IF

            PbV= (Rga/(c1*SGcien*EXP((c3*API(SGo))/(Temperatura+460))))**(1.0/c2)
            RsV= c1*SGcien*Presion**c2*EXP(c3*API(SGo)/(Temperatura+460))

            End subroutine Vazquez


!------------------------------------------------------------------------------------------------------------
!   Subrutina VazquezBo:
!       Calcula el factor de volumen del aceitep oor el método de Vazquez.
!------------------------------------------------------------------------------------------------------------
        Subroutine VazquezBo(API,SGo,Temperatura,SGcien,RsV,BoV)
         IMPLICIT NONE
         Real,External   :: API
         Real,INTENT(IN) ::SGo,Temperatura,SGcien,RsV
         Real,INTENT(OUT)::BoV
         Real::d1,d2,d3

            IF (API(SGo).LE.30.0)Then
             d1 = 4.677E-4
             d2 = 1.751E-5
             d3 = -1.811E-8
            Else
             d1 = 4.670E-4
             d2 = 1.100E-5
             d3 = 1.337E-9
            END IF

            BoV = 1.0+d1*RsV+(Temperatura-60.0)*(API(SGo)/SGcien)*(d2+d3*RsV)
        End subroutine VazquezBo


!------------------------------------------------------------------------------------------------------------
!   Subrutina Standing:
!       Calcula la presion de saturación por el método de Standing
!------------------------------------------------------------------------------------------------------------
        Subroutine Standing(Rga,T,SGg,SGo,API,PbS)
         Implicit none
         Real(kind=4)::API
         Real(kind=4),intent(IN)::Rga,T,SGg,Sgo
         Real(kind=4),intent(OUT)::PbS

          PbS = 18.2*((Rga/SGg)**0.83*10**(0.00091*T-0.0125*API(SGo))-1.4)

        End Subroutine Standing


!------------------------------------------------------------------------------------------------------------
!   Subrutina Standing_Rs_SGd:
!       Se calcula la relacion de solubilidad y la gravedad especifica del gas disuelto
!       a partir de las ecuaciones de Katz y Standing (Garaicochea).
!------------------------------------------------------------------------------------------------------------
        Subroutine Standing_SGd(API,SGo,Rs,SGd)
         Implicit none
         Real(Kind=4),external::API
         Real(kind=4),intent(IN)::SGo,Rs
         Real(kind=4),intent(OUT)::SGd

         !write(*,*)'El valor de RS entra como ',Rs
          SGd = 0.25  +  0.02*API(SGo)  +  0.000001*Rs *(0.6874 - 3.5864*API(SGo))

!          RsS =    455.334473
!          SGd =   0.909999967
        End Subroutine Standing_SGd


!------------------------------------------------------------------------------------------------------------
!   Subrutina Standing_Bo:
!       Factor de Volumen del aceite por el metodo de Standing
!------------------------------------------------------------------------------------------------------------
        Subroutine Standing_Bo(Rs,SGg,SGo,T,BoS)
         Implicit none
         Real(kind=4),intent(in)::Rs,SGg,SGo,T
         Real(kind=4),intent(out)::BoS

            BoS = 0.9759+0.00012*(Rs*(SGg/SGo)**0.5+1.25*T)**1.2

        End subroutine Standing_Bo


!------------------------------------------------------------------------------------------------------------
!   Subrutina Glaso:
!       Calcula la presión de burbuja y la relación de solubilidad, utilizando las correlacion de Glaso.
!------------------------------------------------------------------------------------------------------------
        Subroutine Glaso(API,Rga,Temperatura,Presion,SGg,SGo,RsG,PbG)
        IMPLICIT NONE
        Real,External   :: API
        Real,INTENT(IN) :: Rga,Temperatura,Presion,SGg,SGo
        Real,INTENT(OUT):: RsG,PbG
        Real            ::PG,PGD,LoPb

            PG = 10**(2.8869-(14.1811-3.3093*LOG10(Presion))**0.5)
            RsG = SGg*(API(SGo)**0.989/(Temperatura**0.172)*PG)**1.2255

            IF (API(SGo).GT.40.0)Then
            PGD = (Rga/SGg)**0.816*Temperatura**0.130*API(SGo)**(-0.989)
            Else
            PGD = (Rga/SGg)**0.816*Temperatura**0.172*API(SGo)**(-0.989)
            End IF

            LoPb = 1.7669+(1.7447*LOG10(PGD)-0.30218*LOG10(PGD)**2)
            PbG = 10**(LoPb)
        End subroutine Glaso


!------------------------------------------------------------------------------------------------------------
!   Subrutina GlasoBo:
!       Calcula el factor de volumen del aceite, utilizando las correlacion de Glaso.
!------------------------------------------------------------------------------------------------------------
        Subroutine GlasoBo(RsG,SGo,SGg,Temperatura,BoG)
        IMPLICIT NONE
        Real,INTENT(IN)::RsG,SGo,SGg,Temperatura
        Real, INTENT(OUT)::BoG
        Real::A,Gbob

        Gbob = RsG*(SGg/SGo)**0.526+0.968*Temperatura
        A = -6.58511+2.91329*LOG10(Gbob )-0.27683*LOG10(Gbob)**2
        BoG = 1+10**A

        end subroutine GlasoBo


!--------------------------------------------------------------------------------------------------
!   Subrutina DensityOilCY:
!       Calcula la densidad del aceite para el caso de un yacimiento saturado
!--------------------------------------------------------------------------------------------------
        Subroutine DensityOilCY (SGo,SGd,Rs,Bo,Do_cy)
         implicit none
         real(kind=4),intent(in)::SGo,SGd,Rs,Bo
         real(kind=4),intent(out)::Do_cy

         Do_cy = ( 62.4*SGo + 0.0136*Rs*SGd)/Bo

        End subroutine DensityOilCY


!------------------------------------------------------------------------------------------------------------
!   Subrutina Oil_viscosity:
!       Calcula la viscosidad del aceite: a condiciones de yacimiento; para el calculo de la viscosidad
!       del aceite muerto, y en la presión de burbuja la viscosidad se calcula utilizando la correlación
!       de Beggs & Robinson, la viscosidad a condiciones de flujo se calcula con la corrección de Vazquez
!       & Beggs.
!------------------------------------------------------------------------------------------------------------
        Subroutine Oil_viscosity(Sgo,API,P,T,Pb,Rs,Mu_ob)
         implicit none
         Real,External :: API
         Real(kind=4)::X,b,m,a
         Real(kind=4),intent(in)::SGo,P,T,Pb,Rs
         Real(kind=4),intent(out)::Mu_ob
         Real :: Mu_od
        !------------------------------------------------------------
        !               Aceite muerto
            X = 10**(3.0324-0.02023*API(SGo))/T**1.163
            Mu_od = 10**X-1
        !-------------------------------------------------------------
        !               Aceite Saturado
            b = 5.44*(Rs+150)**(-0.338)
            Mu_ob = (10.715*(Rs+100)**(-0.515))*Mu_od**b
!        !-------------------------------------------------------------
!        !               Aceite Bajo saturado
!            a = -3.9*10**(-5)*P-5
!            m = 2.6*P**1.187*10**a
!            Mu_o = Mu_ob*(P/Pb)**m

        End subroutine Oil_viscosity


!------------------------------------------------------------------------------------------------------------
!   Subrutina Oil_Gas_Interfacial_Tension:
!       Calcula la tension interfacial Gas-Aceite con la correlación de Abul-Majael et al.
!------------------------------------------------------------------------------------------------------------
        Subroutine Oil_Gas_Interfacial_Tension(Rs,SGo,API,T,Sigma_og)
         implicit none
         Real,External::API
         real(kind=4)::k, factor
         real(kind=4),intent(in)::Rs,SGo,T
         real(kind=4),intent(out)::Sigma_og

            If (Rs .LE. 280.75) then
             k = (1+4.4183E-3*Rs**1.0157)**(-1)
            Else
             k = 227.786*Rs**(-1.1367)
            End if

            Sigma_og = k*(38.085-0.259*API(SGo))*(1.7101-1.6949*10**(-3)*T)

            Factor = 0.056379 + 0.94362*exp(-3.8491E-3*Rs)

            Sigma_og = Sigma_og*factor

        End subroutine Oil_Gas_Interfacial_Tension


!------------------------------------------------------------------------------------------------------------
!   Función API:
!       Convierte el argumento (gravedad específica del aceite), a su equivalente en Grados API
!       Agua = 10°API, Agua = 1 Gravedad específica.
!------------------------------------------------------------------------------------------------------------
        Function API(x)
         Real(kind=4),intent(in)::x
         Real(kind=4)::API

            API = (141.5/x)-131.5

        End function API


end module Aceite_PVT
