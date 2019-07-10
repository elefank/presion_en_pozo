module Variables
!************************************************************************************************************
!   Propósito:
!    Alojar todas las variables usadas por los demás modulos y subrutinas del programa;
!
!   Observaciones:
!    El uso de la sentencia SAVE indica que la variable no será modificada durante la
!    ejecución del programa.
!
!   Bitácora de revisiones:
!       Fecha       Programador     Descripción del cambio
!     ==========    ===========     ======================
!     11/12/2015    Hernán Q. V.    Código Original
!     02/01/2015    ''              Se agregaron las variables para el cálculo de la presión en un pozo

!************************************************************************************************************

!------------------------------------------------------------------------------------------------------
!                           Glosario de términos y unidades.
!------------------------------------------------------------------------------------------------------
!                           Variables del problema para un pozo

    Real :: Pwh             ! Presión en la cabeza del pozo     [psi].
    Real :: Twh             ! Temperatura en la cabeza del pozo [°F].
    Real :: Pwf             ! Presión en el fonodo del pozo     [psi].
    Real :: Twf             ! Temperatura en el fondo el pozo   [°F].
    Real :: Long            ! Longitud del pozo                 [ft]
!   Real :: Part            ! Particiones para el calculo       [entero]
!------------------------------------------------------------------------------------------------------
!                           Variables necesarias para calcular las propiedades de los fluidos.

    Real :: Presion         ! Presión del sistema [psi].
    Real :: Temperatura     ! Temperatura [°F].
    Real :: SGo             ! Gravedad específica del aceite [agua = 1].
    Real :: SGg             ! Gravedad específica del gas total producido [aire = 1].
    Real :: RGA             ! Relación gas aceite producido [SCF/STB].
    Real :: Qo_ce           ! Gasto de aceite a condiciones estándar [bbl/día].

!------------------------------------------------------------------------------------------------------
!                           Variables del pozo

    Real :: Diam            ! Diámetro de la tubería [in].
    Real :: Rugo            ! Rugosidad de la tubería [in].
!------------------------------------------------------------------------------------------------------
!                             Variables PVT
    Real :: Pb              ! Presión de saturación del aceite [psi]
    Real :: SGcien          ! Gravedad especifica del separador.
    Real :: Tsep            ! Temperatura del separador [psi]
    Real :: Psep            ! Presión del separador [psi]
    Real :: SGd             ! Gravedad específica del gas disuelto
    Real :: SGf             ! Gravedad específica del gas libre
    Real :: SGt             ! Gravedad específica del gas total
    Real :: Bg              ! Factor de volumen del gas
    Real :: Z               ! Factor de compresibilidad del gas
    Real :: Bo              ! Factor de volumen del aceite [bls/STB]
    Real :: Rs              ! Relación de solubilidad
    Real :: Rho_o           ! Densidad del aceite a condiiones de yacimiento [lb/ft3].
    Real :: Rho_g           ! Densidad del gas a condiiones de yacimiento [lb/ft3].
    Real :: Mu_o            ! Viscosidad del aceite a condiiones de yacimiento [cp].
    Real :: Mu_g            ! Viscosidad del gas a condiiones de yacimiento [cp].
    Real :: Sigma_og        ! Tesnión interfacial gas-aceite [dina/cm]
!------------------------------------------------------------------------------------------------------
!                           Variables calculadas a partir de las propiedades del fluido.

    Real :: VsL, VsG        ! Velocidades superficiales del liquido y gas [ft/s].
    Real :: Mu_n            ! Viscosidad de la mezcla sin resbalamiento [cp].
    Real :: Rho_n           ! Densidad de la mezcla sin resbalamiento [cp].
    Real :: Vm              ! Velocidad de la mezcla en pies sobre segundo [ft/s].
    Real :: Lambda          ! Colgamiento sin rebalamiento o fracción volumétrica de entrada.
!------------------------------------------------------------------------------------------------------
!                           Variables usadas en el método de Beggs & Brill
!
!   Real :: NFrut           ! Número de Froude, [adim].
!   Real :: NLV             ! Número de velocidad del líquido de Duns & Ross, [adim].
    Real :: HL              ! Colgamiento en flujo horizontal.
        Real:: HL_seg,HL_int
    Real :: HL_C            ! Colgamiento corregido por el ángulo de inclinación.
    Real :: HL_P            ! Colgamiento corregido por el método de Payne.
!------------------------------------------------------------------------------------------------------
!                           Variables para caracterizar el tipo de flujo

    Integer :: Direccion    ! Variable que indica si el flujo es ascendete o descendente:
                            ! 1 = ascendente, 2 = descendente.
    Real    :: angulo       ! Angulo de inclinacion de la tuberia en [grados].
!------------------------------------------------------------------------------------------------------
    Integer :: Patron       ! Indica el patron de flujo presente [entero].
    Real    :: AA           ! Variable para ponderar en caso de flujo intermitente.
    Real    :: F_Beggs      ! Factor de fricción característico [adim].
    Real    :: Delta_P      ! Gradiente de presión [psi/ft].
!------------------------------------------------------------------------------------------------------
!                           Variables auxiliares sin significado físico

    Real      :: dumy,criterio
    Character ::ask1,ask2,ask3

!------------------------------------------------------------------------------------------------------
end module Variables
