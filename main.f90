program main
!***************************************************************************************************************
!   Propósito:
!    Este programa permite estimar el gradiente de presión usando varias correlaciones de flujo multifásico
!    y además incluye varias correlaciones PVT para calcular las propiedades de los fluidos
!
!   Observaciones:
!     En esta parte del código solo se hacen llamadas a archivos para entrada y salida de datos y para llamar
!     subrutinas.
!
!   Bitácora de revisiones:
!       Fecha       Programador     Descripción del cambio
!     ==========    ===========     ======================
!     11/12/2015    Hernán Q. V.    Código Original
!     01/02/2015    ''              Prueba con valores PVT fijos (EXITOSA)
!
!***************************************************************************************************************

use Variables
use Pozo

 implicit none
        open (UNIT=8, file='datos.dat')

                    read(8,*)Pwh
                    read(8,*)Twh
                    read(8,*)Twf
                    read(8,*)Qo_ce
                    read(8,*)SGo
                    read(8,*)RGA
                    read(8,*)SGg
                    read(8,*)Diam
                    read(8,*)Rugo
                    read(8,*)Long
                    read(8,*)Direccion
                    read(8,*)angulo

        close (UNIT=8)

        open (UNIT=9, file='resultados.txt')

                    Write(9,2)Pwh,Twh,Twf,Qo_ce,SGo,RGA,SGg,Diam,Rugo,Long,Direccion,angulo

                    2 format('',T37 'Datos usados'                         ,//&
                    T14 'Presion en la cabeza del pozo              'F10.1   '    [PSI]'       ,/,&
                    T14 'Temperatura en la cabeza del pozo          'F10.1   '    [°F]'        ,/,&
                    T14 'Temperatura en el fondo del pozo           'F10.3   '    [°F]'        ,/,&
                    T14 'Gasto de aceite a condiciones estándar     'F10.3   '    [STB/día]'   ,/,&
                    T14 'Gravedad específca del aceite              'F10.3   '    [agua= 1]'   ,/,&
                    T14 'Relacion gas aceite producido              'F10.1   '    [SCF/STB]'   ,/,&
                    T14 'Gravedad especifica del gas total producido'F10.3   '    [aire=1]'    ,/,&
                    T14 'Diámetro de la tubería                     'F10.3   '    [in]'    ,/,&
                    T14 'Rugosidad de la tubería                    'F10.3   '    [in]'    ,/,&
                    T14 'Longitud de la tubería                     'F10.3   '    [ft]'    ,/,&
                    T14 'Dirección                                  'I10     '    [ascendente=1 descendente=2]',/,&
                    T14 'Ángulo de inclinación                      'F10.1   '    [°]'         ,/ )

                    write(*,2)Pwh,Twh,Twf,Qo_ce,SGo,RGA,SGg,Diam,Rugo,Long,Direccion,angulo


                    Write(*,*)'Se conocen las condiciones del separador(Temperatura Presion)??[y/n] '
                    Read(*,*)ask3

                    IF (ask3 =='y'.OR. ask3 == 'Y') Then
                        Write(*,16)
                        16 Format('',T5'Introduce la presion de separacion [PSI]')
                        Read(*,*)Psep

                        Write(*,17)
                        17 Format('',T5'Introduce la temperatura de separacion [ºF]')
                        Read(*,*)Tsep

                        Write(9,18)Tsep,Psep
                        18 Format('',T14'Presion de Separación                      'F10.2'    [PSI]',/,&
                                     T14'Temperatura de Separacion                  'F10.2'    [°F]',/)
                    END IF



            !    Calculo(Pwh,Twh,L,Tf,Qo_ce,RGA,Diam,Rugo,angulo,ask3,direccion,Pwf)
            Call Calculo(Pwh,Twh,Long,Twf,Qo_ce,RGA,Diam,Rugo,angulo,ask3,direccion,Pwf)
        close (UNIT=9)

        write(*,*)'presiona una tecla'
        read(*,*)

end program main
