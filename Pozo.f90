module Pozo
!***************************************************************************************************************
!   Propósito:
!    Este módulo contiene la subrutina que permite encontrar la presión de fondo en un pozo de aceite, conociendo
!    las presion en la cabeza y las temperaturas en el fondo y la cabeza del pozo, además de el gasto de aceite
!    relación gas aceite producido y gravedades especificas del gas y aceite utilizando correlaciones de flujo
!    multifásico.
!
!   Observaciones:
!     En esta parte del código solo se hacen llamadas a archivos para entrada y salida de datos y para llamar
!     subrutinas.
!
!   Bitácora de revisiones:
!       Fecha       Programador     Descripción del cambio
!     ==========    ===========     ======================
!     02/01/2016    Hernán Q. V.    Código original utilizando la correlación de Beggs y Brill.
!     03/01/2016    ''              ...
!
!***************************************************************************************************************


    use Variables
    use BeggsyBrill
    use Funciones

    implicit none

contains

        Subroutine Calculo(Pwh,Twh,L,Tf,Qo_ce,RGA,Diam,Rugo,angulo,ask3,direccion,Pwf)
        Implicit none
        Real(kind=4),intent(in) :: Pwh,Twh,L,Tf,Qo_ce,RGA,Diam,Rugo,angulo
        Integer,     intent(in) :: direccion
        Character,   intent(in) :: ask3
        Real(kind=4),intent(out):: Pwf
        Real(kind=4)            :: Dz,dif,Pprom,Tprom
        Real                    :: P2,Pzdz,P1,particion
        integer                 :: i,part

        !Arreglos de Presión y Temperatura
        Real(kind=4), allocatable,dimension(:)::Presiones
        Real(kind=4), allocatable,dimension(:)::Temperaturas
        Real(kind=4), allocatable,dimension(:)::TempProm
!----------------------------------------------------------------------------------------------------------
        !estos paámetros determinan el tamaño del cálculo ambos deben ser iguales y enteros
        part = 10         !Sin punto decimal
        particion = 10.0  !Con punto decimal
!----------------------------------------------------------------------------------------------------------
        Dz = L/particion !divide la longitud total entre el numero de partes

        allocate(Presiones(part+1),Temperaturas(part+1),TempProm(part))

        Presiones(1)=Pwh
        !Presiones(part+1)=Pwf

        ! Guardamos el perfil de temperaturas en el arreglo
        Do i=1 , (part+1)
          Temperaturas(i)= Tf-((Tf-Twh)/L*(L-(Dz*i)))
        End do
        ! Guardamos el perfil de temperaturas PROMEDIO en el arreglo
        Do i=1, (part)
          If (i .eq. 1) Then
          TempProm(i)= ( Twf+Temperaturas(i))/2.0
          End IF
          TempProm(i)= ( Temperaturas(i-1)+Temperaturas(i) )/2.0
        End do

!        Resolviendo el sistema de la cabeza hacia abajo con referencia en el fondo

      Do i=1, (part)
         Tprom = TempProm(i) ! Temperatura promedio en la celda
         P2 = Presiones(i) ! Presion en la cabeza para el primer paso
         Pprom = ( P2 + Presiones(i) )/2.0  ! Se calcula el valor supuesto de la presión promedio
!                       ----P1----- supuesta

         Call PVT(Pprom,Tprom,Sgo,SGg,RGA,Bo,Bg,Rs,Rho_o,Rho_g,Mu_o,Mu_g,Sigma_og,ask3)
         Call Grad_Beggs_Brill (Qo_ce,RGA,Diam,Rugo,angulo,Bo,Bg,Rs,Rho_o,Rho_g,Mu_o,Mu_g,Sigma_og,Direccion &
                                 ,NFroude,NLV,Reynolds,Delta_P)

        Pzdz = (Delta_P + P2/Dz)*Dz
        dif = (Pzdz-P1)/P1

!        Mientras la diferencia sea mayor que 0.00001 se entra en el siguiente ciclo
        Do While((Abs(dif)).GT. 0.00001)
            P1 = Pzdz
            Pprom = (P2+P1)/2.0
!            Presion = Pprom

            !SUBRUTINA_PVT(IN=Tprom,Pprom :: OUT=Bo,Bg,Rs,Rho_o,Rho_g,Mu_o,Mu_g,Sigma_og)
            Call PVT(Pprom,Tprom,Sgo,SGg,RGA,Bo,Bg,Rs,Rho_o,Rho_g,Mu_o,Mu_g,Sigma_og,ask3)

            Call Grad_Beggs_Brill (Qo_ce,RGA,Diam,Rugo,angulo,Bo,Bg,Rs,Rho_o,Rho_g,Mu_o,Mu_g,Sigma_og,Direccion &
                                 ,NFroude,NLV,Reynolds,Delta_P) ! Solo arroja un Delta_P.

                Pzdz = (Delta_P + P2/Dz)*Dz   !! Despeja e valor en el límite de la celda (P1).
                dif = (Pzdz-P1)/P1
        End do
        !Almacenamos el valor de P1 en el arreglo para el siguiente valor de (i).
        Presiones(i+1)=P1

      !IMPRESIÓN DE CÁLCULOS
      Write(*,20)i
      Write(9,20)i
      20 Format (' Celda ('I3')' )


      Write(*,*)'Presion promedio     :',Pprom,'[PSI]'
      Write(9,*)'Presion promedio     :',Pprom,'[PSI]'

      Write(*,*)'Temperatura promedio :',Temperaturas(i),'[°F]'
      Write(9,*)'Temperatura promedio :',Temperaturas(i),'[°F]'

      Write(*,*)'Gradiente de presión :',Delta_P,'[PSI/ft]'
      Write(9,*)'Gradiente de presión :',Delta_P,'[PSI/ft]'

      write(*,*)''
      write(9,*)''

      End do

      Pwf = P1
      write(*,*)''
      write(*,*)'Presión en el fondo  :',Pwf,'[PSI]'
      write(9,*)'Presión en el fondo  :',Pwf,'[PSI]'
      write(*,*)''


    deallocate(Presiones)
    deallocate(Temperaturas)
    deallocate(TempProm)

    End subroutine Calculo


End module Pozo
