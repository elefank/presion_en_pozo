module BeggsyBrill
!***************************************************************************************************************
!   Propósito:
!    Llevar a cabo el calculo del gradiente de presión usando la correlación de Beggs & Brill
!
!   Observaciones:
!    Blah Blah Blah
!
!   Bitácora de revisiones:
!       Fecha       Programador     Descripción del cambio
!     ==========    ===========     ======================
!     11/12/2015    Hernán Q. V.    Identificacion del patron de flujo.
!     12/12/2015    Hernan Q. V.    Cálculo del colgamiento para flujo burbuja.
!     26/12/2015    Hernán Q. V.    Todos los cálculos necesarios para el gradiente de presión.

!***************************************************************************************************************

    use Variables   ! Hacemos accesibles las declaraciones de las variables contenidas en el módulo "Variables.f90"
    use Rutina      ! Hacemos accesibles los procedimientios contenidos en el módulo "Rutina.f90"
    use Funciones   ! Hacemos accesibles las subrutinas contenidas en el modulo "Funciones.f90"
    implicit none

    contains

!-----------------------------------------------------------------------------------------------------------------
!   Subrutina Grad_Beggs_Brill
!       Condensa las demás subrutinas de este módulo, para calcular el gradiente de presión en [PSI/ft],
!       a partir de los datos de producción y propiedades PVT.
!       Identifica el partón de flujo presente, y corrige el colgamiento por angulo de inclinación, para
!       para todos los patrones de flujo, en dirección ascendente y descendente.
!------------------------------------------------------------------------------------------------------------
    Subroutine Grad_Beggs_Brill (Qo_ce,RGA,Diam,Rugo,angulo,Bo,Bg,Rs,Rho_o,Rho_g,Mu_o,Mu_g,Sigma_og,Direccion &
                                 ,NFroude,NLV,Reynolds,Delta_P)

    Real,    Intent(IN)  :: Qo_ce,RGA,Diam,Rugo,angulo
    Real,    Intent(IN)  :: Bo,Bg,Rs,Rho_o,Rho_g,Mu_o,Mu_g,Sigma_og !De correlaciones PVT
    Integer, Intent(IN)  :: Direccion
    Real,    External    :: NFroude,NLV,Reynolds
    Real,    Intent(OUT) :: Delta_P

             !Props(Qo_ce,RGA,Diam,Bo,Bg,Rs,Rho_o,Rho_g,Mu_o,Mu_g,VsL,VsG,Vm,Lambda,Mu_n)
        Call Props(Qo_ce,RGA,Diam,Bo,Bg,Rs,Rho_o,Rho_g,Mu_o,Mu_g,VsL,VsG,Vm,Lambda,Mu_n)
            !Patron_flujo(Vm,Diam,Lambda,NFroude,Patron,AA)
        Call Patron_flujo(Vm,Diam,Lambda,NFroude,Patron,AA)
            !Colgamiento(Patron,Lambda,Vm,Diam,NFroude,HL,HL_seg,HL_int)
        Call Colgamiento(Patron,Lambda,Vm,Diam,NFroude,HL,HL_seg,HL_int)
            !Correc_angulo(Vm,Diam,Lambda,VsL,Rho_L,Sigma_GL,HL,Patron,Direccion,angulo,NLV,Nfroude,MySIN,HL_C)
        Call Correc_angulo(Vm,Diam,Lambda,VsL,Rho_o,Sigma_og,HL,Patron,Direccion,angulo,NLV,Nfroude,MySIN,HL_C)
        ! si el patrón de flujo es segregado
        IF (Patron .EQ. 2)Then
                !transicion(Rho_L,Sigma_GL,HL_seg,HL_int,AA,direccion,angulo,HL_C)
            Call transicion(Rho_o,Sigma_og,HL_seg,HL_int,AA,direccion,angulo,HL_C)
        End IF

            !Correc_Payne(HL_C,angulo,HL_P)
        Call Correc_Payne(HL_C,angulo,HL_P)
            !frictionBeggs(Reynolds,Rugo,Rho_n,Mu_n,Lambda,HL_C,F_Beggs)
        Call frictionBeggs(Reynolds,Rugo,Rho_n,Mu_n,Lambda,HL_C,F_Beggs)
            !Delta_Presion(MySIN,Presion,Rho_o,Rho_g,Rho_n,Vm,Diam,angulo,HL,Fr,Delta_P)
        Call Delta_Presion(MySIN,Presion,Rho_o,Rho_g,Rho_n,Vm,Diam,angulo,HL_P,F_Beggs,Delta_P)

    End Subroutine Grad_Beggs_Brill


!------------------------------------------------------------------------------------------------------------
!   Subrutina Patron_flujo
!       Calcula los parámetros característicos de la ecuación de Beggs & Brill, para identificar el patrón
!       de flujo presente en la tubería.
!------------------------------------------------------------------------------------------------------------
     Subroutine Patron_flujo(Vm,Diam,Lambda,NFroude,Patron,AA)
     implicit none
!       Objetivo: determinar el patrón de flujo

        Real,Intent(IN) :: Vm,Diam,Lambda   ! Definidas previamente en el módulo "Variables.f90"
        Integer, Intent(OUT) :: Patron      ! Variable que sirve para indicar el patrón de flujo presente
        Real,Intent(OUT) :: AA              ! Variable para ponderar en caso de que exista flujo en transición
        Real, External :: NFroude
        Real :: NFrut                       ! Número de Froude (definido en una función fuera de la subrutina)
        Real :: L1,L2,L3,L4                 ! Variables locales auxiliares

        ! Cuerpo de la subrutina
         L1 = 316.0*Lambda**0.302
         L2 = 0.000925*Lambda**(-2.468)
         L3 = 0.10*Lambda**(-1.452)
         L4 = 0.5*Lambda**(-6.738)

         AA = (L3-NFroude(Vm,Diam))/(L3-L2) ! esta variable será usada en caso de que exista flujo en transición

        ! NFrut = Vm**2.0 / ( 32.174*(Diam/12) )

          NFrut = NFroude(Vm,Diam)

        IF (Lambda.LT.0.01 .and. NFrut.LT.L1) Then
            Patron = 1
        Else IF (L1.GE.0.01 .and. Nfrut.LT.L2)Then
            Patron = 1
        Else IF (L1.GE.0.01 .and. L2.LE.Nfrut .and. Nfrut.LE.L3)Then
            Patron = 2
        Else IF (0.01.LE.Lambda .and. lambda.LT.0.4 .and. L3.LT.Nfrut .and. Nfrut.LE.L1)Then
            Patron = 3
        Else IF (Lambda.GE.0.4 .and. L3.LT.Nfrut .and. Nfrut.LE.L4)Then
            Patron = 3
        Else IF (Lambda.LT.0.4 .and. Nfrut.GE.L1)Then
            Patron = 4
        Else IF(Lambda.GE.0.4 .and. Nfrut.GT.L4)Then
            Patron = 4
        End IF
       ! Patron igual a
       ! 1 ! El patron de flujo es segregado
       ! 2 ! El patron de flujo es transición
       ! 3 ! El patron de flujo es intermitente
       ! 4 ! El patrón de flujo es disrtibuido
     End Subroutine Patron_flujo

!------------------------------------------------------------------------------------------------------------
!   Subrutina Colgamiento:
!       Calcula el colgamiento considerando flujo horizontal, dependiendo del patrón de flujo presente
!------------------------------------------------------------------------------------------------------------
     Subroutine Colgamiento(Patron,Lambda,Vm,Diam,NFroude,HL,HL_seg,HL_int)
     implicit none
     !Objetivo : Calcular el colgamiento
      Integer,Intent(IN) :: Patron
      Real, Intent(IN)   ::Lambda,Vm,Diam
      Real, Intent(OUT)  :: HL,HL_seg,HL_int         ! Colgamiento sin ángulo
      Real, External :: NFroude                     ! Función definida fuera de la subrutina
      Real :: a,b,c                              ! Constantes auxiliares

        Select Case(Patron)
            Case(1)    ! Segregado
                a = 0.980
                b = 0.4846
                c = 0.0868
                HL = a*Lambda**b/(NFroude(Vm,Diam)**c)
            Case(2)    ! Transición
!               Calcular Colgamientos en ambos patrones Segregado e Intermitente
                a = 0.980
                b = 0.4846
                c = 0.0868
                HL_seg = a*Lambda**b/(NFroude(Vm,Diam)**c)

                a = 0.845
                b = 0.5351
                c = 0.0173
                HL_int = a*Lambda**b/(NFroude(Vm,Diam)**c)

            Case(3)    ! Intermitente
                a = 0.845
                b = 0.5351
                c = 0.0173
                HL = a*Lambda**b/(NFroude(Vm,Diam)**c)
            Case(4)    ! Distribuido
                a = 1.065
                b = 0.5824
                c = 0.0609
                HL = a*Lambda**b/(NFroude(Vm,Diam)**c)
        End Select

    End subroutine Colgamiento

!-----------------------------------------------------------------------------------------------------------------
!   Subrutina Correc_agulo:
!       Aplica la correción al colgamiento por inclinación, dependiendo del patrón de flujo presente,
!       en el caso de flujo en transición no calcula nada.
!------------------------------------------------------------------------------------------------------------
    Subroutine Correc_angulo(Vm,Diam,Lambda,VsL,Rho_L,Sigma_GL,HL,Patron,Direccion,angulo,NLV,Nfroude,MySIN,HL_C)
    implicit none
    !Objetivo : Calcular el colgamiento corregido por el ángulo de inlcinación, según el patrón de flujo presente
       Real,External :: NLV,Nfroude,MySIN
       Real,Intent (IN) :: Lambda,Vm,Diam,VsL,Rho_L,Sigma_GL! Variables definidas anteriormente
       Integer,Intent(IN) :: Patron,Direccion
       Real,Intent (IN) :: HL,angulo                        ! Colgamiento en flujo horizontal
       Real,Intent (OUT) :: HL_C                            ! Colgamiento corregido por el ángulo
       Real :: Chi,C,theta                                  ! Variables Auxiliares
       Real :: e,f,g,h                                      ! Variables Auxiliares
           !Cuerpo
       IF (Direccion .eq. 1)Then   ! Cuando el flujo en la tubería es ascendente "uphill".
         Select Case(Patron)
            Case(1)    ! Segregado
                e = 0.011
                f = -3.7680
                g = 3.5390
                h = -1.6140
!            Case(2)    ! Transición
!                !HAY QUE PONDERAR
            Case(3)    ! Intermitente
                e = 2.960
                f = 0.3050
                g = -0.4473
                h = 0.0978
            Case(4)    ! Distribuido
                HL_C = HL  ! NO HAY CORRECCIÓN
                Return
        End Select
       Else IF(Direccion .eq. 2) Then  ! Cuando el flujo en la tubería es descendente "downhill".
        e = 4.7
        f = -0.3692
        g = 0.1244          ! Aplica para todos lo patrones de flujo
        h = -0.5056
       End IF

           C = (1.0-Lambda)*LOG( e*(Lambda**f) *(( NLV(VsL,Rho_L,Sigma_GL) )**g) *(( NFroude(Vm,Diam) )**h))
           theta = angulo
           IF (C .GE. 0) Then    ! Restricción del método
            Chi = 1.0 + C*(MySIN(1.8*theta)-0.333*MySIN(1.8*theta)**3)
           Else
            Chi = 1.0
           End IF

           HL_C = HL * Chi

    End Subroutine Correc_angulo


!-----------------------------------------------------------------------------------------------------------------
! Aqui va la corrección del colgamiento por angulo de inclinacion, para patrón de flujo en transición:
! Se obtiene después de las correcciones de ángulo de inclinación antes de la corrección de Payne,
! La subrutina llama recursivamente a la subrutina que corrige por el ángulo, cambiando la variable Patron
! y el colgamiento de entrada por los valores de HL_seg y HL_int con Patron=1 y Patron=3 respectivamente,
! El valor de colgamiento resultante se guarda en una variable interna, para luego ponderar y obtener el
! colgamiento ponderado para flujo en transicion y corregido por ángulo.
!-----------------------------------------------------------------------------------------------------------------
    Subroutine transicion(Rho_L,Sigma_GL,HL_seg,HL_int,AA,direccion,angulo,HL_C)
    implicit none
    !Objetivo : Calcular el colgamiento corregido para el ángulo de inclinación cuando el patrón de flujo presente es transición
        Real,Intent(IN) :: HL_seg,HL_int,angulo,AA
        Real,Intent(In) :: Rho_L,Sigma_GL
        Integer,Intent(IN) :: direccion
        Real,Intent(OUT) :: HL_C
!        Integer :: Patron
        Real :: Segregado,Intermitente

        ! Primero Para flujo segregado
        ! Patron = 1  *Segregado
           ! Correc_angulo(Vm,Diam,Lambda,VsL,Rho_L,Sigma_GL,HL,Patron,Direccion,angulo,NLV,Nfroude,HL_C)
        call Correc_angulo(Vm,Diam,Lambda,VsL,Rho_L,Sigma_GL,HL_seg,1,Direccion,angulo,NLV,Nfroude,MySIN,Segregado)

        ! Luego para flujo intermitente
        ! Patron = 3 *Intermitente
           ! Correc_angulo(Vm,Diam,Lambda,VsL,Rho_L,Sigma_GL,HL,Patron,Direccion,angulo,NLV,Nfroude,HL_C)
        call Correc_angulo(Vm,Diam,Lambda,VsL,Rho_L,Sigma_GL,HL_int,3,Direccion,angulo,NLV,Nfroude,MySIN,Intermitente)


                HL_C = AA*Segregado + (1-AA)*Intermitente

    End subroutine transicion

!----------------------------------------------------------------------------------------------------------------------
!   Subrutina Correc_Payne:
!       Aplica la correción de Payne al colgamiento, según el ángulo de inclinación.
!------------------------------------------------------------------------------------------------------------
    Subroutine Correc_Payne(HL_C,angulo,HL_P)
    implicit none
    !Objetivo : aplicar la corrección propuesta por Payne et. al.
        Real,Intent(IN) :: HL_C,angulo
        Real,Intent(OUT) :: HL_P

            !Cuerpo
            IF (angulo>0) Then
                HL_P = 0.924*HL_C
            Else IF(angulo<0) Then
                HL_P = 0.685*HL_C
            End IF
            !Restricción
            IF(HL_P<Lambda)Then
                HL_P = HL_C
            End IF

    End Subroutine Correc_Payne

!----------------------------------------------------------------------------------------------------------------
!   Subrutina frictionBeggs:
!       Calcula el factor de fricción característico de la mezcla utilizado en la correlación de
!       Beggs & Brill.
!------------------------------------------------------------------------------------------------------------
    Subroutine frictionBeggs(Reynolds,Rugo,Rho_n,Mu_n,Lambda,HL_C,F_Beggs)
    implicit none
    !Objetivo calcular el factor de fricción corregido.
     Real, External :: Reynolds
     Real,Intent(IN):: Rugo,Rho_n,Mu_n,Lambda,HL_C
     Real,Intent(OUT) :: F_Beggs
     Real :: y,s,ffn,Reynolds_m,friccion_m

        y = Lambda/(HL_C**2)

        If( 1.0.LT.y .and. y.LT.1.2  )Then
            s = LOG(2.2*y-1.2)
        Else
            s = LOG(y)/( -0.0523 + 3.182*LOG(y) - 0.8725*(LOG(y))**2 + 0.01853*(LOG(y))**4 )
        End If

        ffn = EXP(s)

        Reynolds_m = Reynolds(Rho_n,Vm,Diam,Mu_n)
!        Write(*,*)'Reynolds_m = ' ,Reynolds_m

        call Factor_de_friccion(Rugo,Diam,Reynolds_m,friccion_m)
!        Write(*,*)'Friccion = ',Friccion_m

        F_Beggs = friccion_m * ffn

    End subroutine frictionBeggs

!------------------------------------------------------------------------------------------------------------
!   Subrutina Delta_Presion:
!        Calcula el gradiente de presión usando la correlación de Beggs & Brill,
!        utilizando el término de pérdida convección.
!------------------------------------------------------------------------------------------------------------
    Subroutine Delta_Presion(MySIN,Presion,Rho_o,Rho_g,Rho_n,Vm,Diam,angulo,HL,Fr,Delta_P)
    Real, External     :: MySIN
    Real, Intent (IN)  :: Presion,Rho_o,Rho_g,Rho_n,Vm,Diam,angulo
    Real, Intent (IN)  :: HL,Fr
    Real, Intent (OUT) :: Delta_P
    Real               :: Rho_s,Ek

    Rho_s = Rho_o*HL + Rho_g*(1-HL)

    Ek = Vm*VsG*Rho_n/Presion

    Delta_P = ( Fr*Rho_n*(Vm**2.0) /( 2.0*32.174*(Diam/12.0) ) +  Rho_s*MySIN(angulo)   ) !/ (1-Ek)
    Delta_P = Delta_P / 144.0

    End subroutine Delta_Presion

!------------------------------------------------------------------------------------------------------------
! Subrutina Factor de Fricción:
!    Calcula el factor de friccion usando la ecuacion de Colebrook y White resolviendo iterativamante.
!------------------------------------------------------------------------------------------------------------
    Subroutine Factor_de_friccion(e,D,Re,f)
    implicit none
    ! Objetivo : Calcular el factor de fricción con la ecuacion de Colebrook.
     Real,Intent(IN) :: e,D,Re
     Real,Intent(OUT):: f
     Real            :: fs,fs2,df
        !e = 0.0006
        !D = 2.5
        !Re = 12915.736765
        f = 0.01
        fs = 1/sqrt(f)
        fs2= -0.869*log((e/D)/3.7+2.523/(Re*sqrt(f)))
        df =fs-fs2

        ! Cuerpo del cálculo
        Do while ((ABS(df)).GT. 0.0001)
            f = f+0.00000001
            fs = 1/sqrt(f)
            fs2 = -0.869*log((e/D)/3.7+2.523/(Re*sqrt(f)))
            df = fs-fs2
        end do
    End Subroutine Factor_de_friccion

!------------------------------------------------------------------------------------------------------------
!   Funcion NLV :
!       Calcula el número adimensional del líquido de Duns & Ross, parametro que sirve para
!       caracterizar el patrón e flujo en la correlación de Beggs &  Brill.
!------------------------------------------------------------------------------------------------------------
    Function NLV(VsL,Rho_L,Sigma_GL)          ! Número de velocidad del líquido
        implicit none
        Real,Intent(IN) ::VsL                 ! Velocidad superficial del líquido [ft/s]
        Real,Intent(IN) ::Rho_L               ! Densidad del líquido [Lbm/ft^3]
        Real,Intent(IN) ::Sigma_GL            ! Tensión superficial Gas-Liquido [dina/cm]
        Real :: NLV
            !Cuerpo de la función:
            NLV = 1.938*VsL*(Rho_L/Sigma_GL)**(1.0/4.0)
    End Function NLV

!------------------------------------------------------------------------------------------------------------
!   Función Nfroude :
!       Calcula el paramero adimensional del número de Froude, que sirve para identificar el patrón
!       de flujo en la correlaci+on de Beggs & Brill.
!------------------------------------------------------------------------------------------------------------
     Function NFroude(Vm,Diam)         ! Numero de Froude
        implicit none
        Real,Intent(IN) :: Vm          ! Velocidad superficial de la mezcla en en pies por segundo [ft/s]
        Real,Intent(IN) :: Diam        ! Diámetro en pulgadas [in]
        Real :: NFroude
         ! Cuerpo de la función:
         NFroude = Vm**2.0 / ( 32.174*(Diam/12.0) )
     End Function NFroude

!------------------------------------------------------------------------------------------------------------
!   Función Reynolds:
!       Calcula el parámetro adimensional del número de Reynolds, que caracteriza el régimen de flujo
!       que sirve para calcular el factor de fricción.
!------------------------------------------------------------------------------------------------------------
     Function Reynolds(Rho,V,Diam,Mu)      ! Numero de Reynolds
        implicit none
        Real,Intent(IN) :: Rho,V,Diam      ! Velocidad superficial de la mezcla en en pies por segundo [ft/s]
        Real,Intent(IN) :: Mu              ! Diámetro en pulgadas [in]
        Real :: Reynolds
         ! Cuerpo de la función:
         Reynolds = Rho*V*(Diam/12.0)*1488.0 / Mu
     End Function Reynolds

end module BeggsyBrill
