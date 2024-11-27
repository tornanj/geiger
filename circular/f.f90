PROGRAM EJ
	
	USE SIMULACION

	IMPLICIT NONE
		TYPE(NODO), ALLOCATABLE, DIMENSION(:,:) :: A
		REAL(8) :: volts
		INTEGER :: n, r_int, r_ext
		
		!CALL SYSTEM('gnuplot -persist script.txt')
		
		WRITE(*,'(A)',ADVANCE='NO') 'n = '
		READ*, n
		
		ALLOCATE(A(n,n))
		
		WRITE(*,'(A)',ADVANCE='NO') 'radio interno = '
		READ*, r_int
		
		WRITE(*,'(A)',ADVANCE='NO') 'radio externo = '
		READ*, r_ext
		
		volts = 600.0d0
		
		CALL CREAR_CIRCULO(A, volts, r_int, r_ext)
		
		CALL SOBRE_RELAJACION(A, volts, 0.00001d0, 1.5d0, r_int, r_ext)
		
		CALL CARGAR_VECTORES(A, r_int, r_ext)
		
		CALL CARGAR_ARCHIVO(A,1.0d0, 1.0d0)
		
		CALL MOVIMIENTO(A, r_ext, r_int)
		
		DEALLOCATE(A)
		
	CONTAINS
	
		SUBROUTINE CREAR_CIRCULO(A, volts, r_int, r_ext)
			TYPE(NODO), DIMENSION(:,:) :: A
			REAL(8) :: volts, aleat
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
									CALL RANDOM_NUMBER(aleat)
									IF ( aleat <= 0.3d0 ) THEN
										A(i,j)%tipo = GASNOBLE
									ELSE
										A(i,j)%tipo = VACIO
										A(i,j)%analizado = 0
									END IF
									
									A(i,j)%valor = 0.0d0
									A(i,j)%grado = 0.0d0
							END IF
					END IF
				END DO
			END DO
			
		END SUBROUTINE

		SUBROUTINE SOBRE_RELAJACION(A, volts, tol, omega, r_int, r_ext)
        TYPE(NODO), DIMENSION(:,:) :: A
        TYPE(NODO), ALLOCATABLE,DIMENSION(:,:) :: B_ant
        REAL(8) :: tol, omega, volts, error
        INTEGER :: i, j, n, iter_max, iter, r_int, r_ext, medio

		n = SIZE(A,DIM = 1)

        ALLOCATE(B_ant(n,n))
		
		iter = 1
		iter_max = 10000
		medio = INT(n/2)

        ! Inicializamos B con los valores de A
        B_ant = A
    
        ! Iteramos hasta que el error máximo sea menor que la tolerancia
        do while (.true.)
            error = 0.0
            do i = 2, n - 1
                do j = 2, n - 1
                    ! Guardamos el valor anterior
                    B_ant(i,j)%valor = A(i,j)%valor
					
						! Sobrerelajación: cálculo del nuevo valor ajustado con omega
					IF (  ( ( i - medio )**2 + ( j - medio )**2 >= r_int**2 ) .and.  &
							( ( i - medio )**2 + ( j - medio )**2 <= r_ext**2 ) )  THEN
							
							A(i,j)%valor = (1.0d0 - omega) * A(i,j)%valor + omega * 0.25d0 * &
											(A(i+1,j)%valor + A(i-1,j)%valor+ A(i,j+1)%valor + A(i,j-1)%valor)
					END IF
						! Calculamos el error absoluto para verificar la convergencia
						error = max(error, abs(B_ant(i,j)%valor - A(i,j)%valor))
				END DO
            end do
			iter = iter + 1

            ! Verificamos si el error es menor que la tolerancia
            IF (error < tol .or. iter > iter_max) THEN
				PRINT*, 'Numero de iteraciones y tolerancia: ', iter, tol
				DEALLOCATE(B_ant)
				return
            END IF	
        end do
		
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
				IF (  r_int < r  .and. r < r_ext .and. i >= 2 .and. j >= 2 ) THEN
					Ex = - ( A(i+1,j)%valor - A(i-1,j)%valor) / 2.0d0
					Ey = - ( A(i,j+1)%valor - A(i,j-1)%valor) / 2.0d0
					Ex = Ex * 0.5d0
					Ey = Ey * 0.5d0
					WRITE(1,*) i, j, Ex, Ey
									
				END IF
			END DO
		END DO 		

		CLOSE(1)
						
    END SUBROUTINE
   
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
























