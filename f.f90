PROGRAM EJ
	IMPLICIT NONE
		INTEGER :: VACIO = 0, GASNOBLE = 1, POSITIVO = 2, NEGATIVO = 3, CHOQUE = 4
		
		TYPE :: NODO
			REAL(8) :: valor
			INTEGER :: diag, tipo
		END TYPE NODO
	
		TYPE(NODO), ALLOCATABLE, DIMENSION(:,:) :: A
		REAL(8) :: volts
		INTEGER :: n, r_int, r_ext
		
		WRITE(*,'(A)',ADVANCE='NO') 'n = '
		READ*, n
		
		ALLOCATE(A(n,n))
		
		WRITE(*,'(A)',ADVANCE='NO') 'radio interno = '
		READ*, r_int
		
		WRITE(*,'(A)',ADVANCE='NO') 'radio externo = '
		READ*, r_ext
		
		
		WRITE(*,'(A)',ADVANCE='NO') 'volts = '
		READ*, volts
		
		CALL CREAR_CIRCULO(A, volts, r_int, r_ext)
		
		CALL SOBRE_RELAJACION(A, volts, 0.00001d0, 1.5d0, r_int, r_ext)
		
		CALL CARGAR_VECTORES(A, r_int, r_ext)
		
		CALL CARGAR_ARCHIVO(A,1.0d0, 1.0d0)
		
		CALL SYSTEM('gnuplot -persist script.txt')
		
		DEALLOCATE(A)
		
	CONTAINS
	
		SUBROUTINE MUESTRA_MAT(A)
			TYPE(NODO), DIMENSION(:,:) :: A
			INTEGER :: n, i, j
			
			n = SIZE(A,DIM=1)
			DO i = 1, n
				DO j = 1, n
					WRITE(*,'(F5.0,1X)',ADVANCE='NO') A(i,j)%valor
				END DO
				WRITE(*,*)
			END DO
			WRITE(*,*)
		END SUBROUTINE
	
		SUBROUTINE CREAR_CIRCULO(A, volts, r_int, r_ext)
			TYPE(NODO), DIMENSION(:,:) :: A
			REAL(8) :: volts
			INTEGER :: n, i, j, r_int, r_ext, medio
			
			n = SIZE(A, DIM = 1)
			
			medio = INT(n/2)
			
			DO i = 1, n
				DO j = 1, n
					IF ( ( i - medio )**2 + ( j - medio )**2 < r_int**2 )  THEN
							A(i,j)%valor = volts
					ELSE
							IF ( ( i - medio )**2 + ( j - medio )**2 > r_ext**2 ) THEN
									A(i,j)%valor = -volts
							ELSE
									A(i,j)%valor = 0.0d0
							END IF
					END IF
				END DO
			END DO
			
			! ( x - a )^2 + ( y - b )^2 = r^2
			
			! ( i - a )^2 + ( j - b )^2 = r^2
		END SUBROUTINE

		SUBROUTINE SOBRE_RELAJACION(A, volts, tol, omega, r_int, r_ext)
        TYPE(NODO), DIMENSION(:,:) :: A
        TYPE(NODO), ALLOCATABLE,DIMENSION(:,:) :: B_ant, B_aux
        REAL(8) :: tol, omega, volts, error
        INTEGER :: i, j, n, iter_max, iter, r_int, r_ext, medio

		n = SIZE(A,DIM = 1)

        ALLOCATE(B_aux(n,n),B_ant(n,n))
		
		iter = 1
		iter_max = 10000
		medio = INT(n/2)

        ! Inicializamos B con los valores de A
        B_aux = A
    
        ! Iteramos hasta que el error máximo sea menor que la tolerancia
        do while (.true.)
            error = 0.0
            do i = 2, n - 1
                do j = 2, n - 1
                    ! Guardamos el valor anterior
                    B_ant(i,j)%valor = B_aux(i,j)%valor
					
						! Sobrerelajación: cálculo del nuevo valor ajustado con omega
					IF ( ( ( i - medio )**2 + ( j - medio )**2 >= r_int**2 ) .and.  &
							( ( i - medio )**2 + ( j - medio )**2 <= r_ext**2 ) ) THEN
							
							B_aux(i,j)%valor = (1.0d0 - omega) * B_aux(i,j)%valor + omega * 0.25d0 * &
											(B_aux(i+1,j)%valor + B_aux(i-1,j)%valor+ B_aux(i,j+1)%valor + B_aux(i,j-1)%valor)
					END IF
						! Calculamos el error absoluto para verificar la convergencia
						error = max(error, abs(B_ant(i,j)%valor - B_aux(i,j)%valor))
				END DO
            end do
			iter = iter + 1

            ! Verificamos si el error es menor que la tolerancia
            IF (error < tol .or. iter > iter_max) THEN
				A = B_aux
				return
            END IF	
        end do
        DEALLOCATE(B_aux, B_ant)
		
    END SUBROUTINE
    
    SUBROUTINE CARGAR_VECTORES(A, r_int, r_ext)
		TYPE(NODO), DIMENSION(:,:) :: A
		REAL(8) :: dx, dy, Ex, Ey
		INTEGER :: n, i, j, k, r, r_int, r_ext, medio
		
		n = SIZE(A,DIM=1)
		medio = INT(n/2)
		
		OPEN(UNIT =1, FILE='vectores.txt')

		k = INT(n / 10)


		DO i = 1, n, k
			DO j = 1, n, k 
				r = RADIO(i,j,medio)
				IF (  r_int <= r .and. r <= r_ext ) THEN
					Ex = - ( A(i+1,j)%valor - A(i-1,j)%valor) / 2.0d0
					Ey = - ( A(i,j+1)%valor - A(i,j-1)%valor) / 2.0d0
					WRITE(1,*) i, j, Ex, Ey		
				END IF
			END DO
		END DO 
		

		CLOSE(1)
	!Ex = - ( A(i+1,j)%valor - A(i-1,j)%valor) / 2.0d0
	!Ey = - ( A(i,j+1)%valor - A(i,j-1)%valor) / 2.0d0
	!WRITE(1,*) i, j, Ex, Ey
						
    END SUBROUTINE


	FUNCTION RADIO(i, j, medio)
		INTEGER :: RADIO, i, j, medio

		RADIO = INT( sqrt( REAL((i-medio)**2 + (j-medio)**2 )) )
	END FUNCTION
   

	SUBROUTINE CARGAR_ARCHIVO(A, h, k)
        TYPE(NODO),DIMENSION(:,:) :: A
        REAL(8) :: dx, dy, h, k
		INTEGER :: n, i, j

        n = SIZE(A,DIM=1)

        OPEN(unit = 1, file='output.txt')
        dy = 0.0d0

        DO i = n, 1, -1
            dx = 0.0d0
            DO j = 1, n
                WRITE(1,*) dx, dy, A(i,j)%valor
                dx = dx + h
            END DO
            dy = dy + k
            WRITE(1,*)
        END DO

        CLOSE(1)

    END SUBROUTINE
	
END PROGRAM
























