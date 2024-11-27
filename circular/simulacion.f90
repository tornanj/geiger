MODULE SIMULACION

	INTEGER, PARAMETER :: VACIO = 0, POSITIVO = 1, NEGATIVO = 2, GASNOBLE = 3, CHOQUE = 4
	REAL(8), PARAMETER :: PI = 3.14159265359
			
	TYPE :: MOV
		INTEGER :: HOR, VER
	END TYPE MOV		
			
	TYPE :: NODO
		TYPE(MOV) :: MOV
		REAL(8) :: valor, grado
		INTEGER :: analizado, tipo
	END TYPE NODO

	CONTAINS
	
		SUBROUTINE MOVIMIENTO(A, r_ext, r_int)
			TYPE(NODO), DIMENSION(:,:) :: A
			REAL(8) :: theta, r
			INTEGER :: n, medio, x, y, r_ext, r_int
			
						
			OPEN(UNIT = 1,FILE='electrones.txt')
			WRITE(1,*)
			CLOSE(1)
			
			OPEN(UNIT = 2, FILE='protones.txt')
			WRITE(2,*)
			CLOSE(2)
			
						
			r = REAL( r_ext - 20 )
			n = SIZE(A,DIM=1)
			medio = INT(n/2)
			
			theta = 0.0d0
			DO WHILE ( r > r_int )
			
				IF ( theta >= PI * 2.0d0 ) THEN
					theta = 0.0d0
				END IF
				
 				x = medio + NINT ( r * ( COS(theta) ) ) 
				y = medio + NINT(  r * ( SIN(theta) ) )
				WRITE(3,*) x, y
				
				IF ( A(x,y)%tipo == GASNOBLE ) THEN
					A(x,y+1)%tipo = POSITIVO
					A(x,y+1)%grado = theta
					
					A(x,y-1)%tipo = NEGATIVO
					A(x,y-1)%grado = theta
				END IF
				
				r = r * 0.97d0
				theta = theta + pi / 10.0d0
				
				OPEN(UNIT = 3, FILE='electron_caotico.txt')
				WRITE(3,*) x,y
				CLOSE(3)
				
				CALL BUSCAR_MOV(A, r_ext, r_int)
				CALL SYSTEM('gnuplot -persist script.txt')
				CALL DESMARCAR(A)
				print*, 'seguir ? '
				read*,
			END DO
			
			CLOSE(3)
			
			CALL SYSTEM('gnuplot -persist script.txt')
			
		
		END SUBROUTINE
	
		SUBROUTINE BUSCAR_MOV(A, r_ext, r_int)
			TYPE(NODO), DIMENSION(:,:) :: A
			REAL(8) :: r, theta
			INTEGER :: i, j, n, medio, x, y, tipo, r_ext, r_int
			
			n = SIZE(A, DIM = 1)
			medio = INT(n/2)
			
			OPEN(UNIT = NEGATIVO, FILE='electrones.txt')
			OPEN(UNIT = POSITIVO, FILE='protones.txt')
			
			DO i = 1, n 
				DO j = 1, n
					r = real(RADIO(i,j,medio) )
					IF ( (A(i,j)%tipo == POSITIVO .or. A(i,j)%tipo == NEGATIVO ) &
						.and. r <= r_ext .and. r >= r_int .and. A(i,j)%analizado == 0) THEN
						
						r = sqrt( real((i-medio)**2) + real((j-medio)**2) )
						IF ( A(i,j)%tipo == POSITIVO ) THEN
								r = r * 1.1d0  
								theta = A(i,j)%grado - Pi / 10.0d0
								tipo = POSITIVO
						ELSE
								r = r * 0.9d0  
								theta = A(i,j)%grado + Pi / 10.0d0
								tipo = NEGATIVO
						END IF
						
						x = medio + NINT( r * COS(theta) )
						y = medio + NINT( r * SIN(theta) )
						r = real(RADIO(x,y,medio) )
						
						IF ( r <= r_ext .and. r >= r_int  ) THEN
							IF ( tipo == NEGATIVO .and. A(x,y)%tipo == GASNOBLE  ) THEN
								A(x,y+3)%tipo = POSITIVO
								A(x,y+3)%grado = theta
					
								A(x,y-3)%tipo = NEGATIVO
								A(x,y-3)%grado = theta
								
								WRITE(NEGATIVO, *) x,y-3
								WRITE(POSITIVO, *) x,y+3
							END IF
							
							A(x,y)%grado = theta
							A(x,y)%tipo = tipo
							WRITE(tipo,*) x,y
						END IF
						
						A(x,y)%analizado = 1
						A(i,j)%tipo = VACIO
						
					END IF
				END DO
			END DO
			CLOSE(POSITIVO)
			CLOSE(NEGATIVO)
			
		END SUBROUTINE
		
	SUBROUTINE DESMARCAR(A)
		TYPE(NODO), DIMENSION(:,:) :: A
		INTEGER :: i, j, n
		
		n = SIZE(A, DIM = 1)
		DO i =1 , n
			DO j = 1, n
				A(i,j)%analizado = 0
			END DO
		END DO
	
	END SUBROUTINE
		
	FUNCTION RADIO(i, j, medio)
		INTEGER :: RADIO, i, j, medio

		RADIO = INT( sqrt( REAL((i-medio)**2 + (j-medio)**2 )) )
	END FUNCTION
	
	
END MODULE
