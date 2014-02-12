MODULE SURFACE_MODULE

    IMPLICIT NONE;
    
    REAL(8), PARAMETER :: MAXSIZE = 1., MINERROR = 1E-16, PI = 3.141592653589793238, MINLOCERROR = 1E-2, MINNORMERROR = 5E-2

    TYPE POINTPHI
        REAL(8) :: X
        REAL(8) :: Y
        REAL(8) :: Z
        
        REAL(8) :: PHI
        REAL(8) :: TEMPPHI
        
        REAL(8) :: REINITIAL_SIGN
        REAL(8) :: K1
        REAL(8) :: K2
        
        REAL(8) :: B_RATE
        
        REAL(8) :: U
        REAL(8) :: V
        REAL(8) :: W
        
        LOGICAL :: B_RATE_FIXED
        
        TYPE(POINTPHI), POINTER :: LEFT
        TYPE(POINTPHI), POINTER :: RIGHT
        TYPE(POINTPHI), POINTER :: TOP
        TYPE(POINTPHI), POINTER :: BOTTOM
        TYPE(POINTPHI), POINTER :: FRONT
        TYPE(POINTPHI), POINTER :: BACK
        
        TYPE(POINTPHI), POINTER :: NEXT
        TYPE(POINTPHI), POINTER :: BEFORE
    END TYPE

    TYPE OCTREE
        INTEGER :: DEPTH
        REAL(8) :: LENGTH
        
        TYPE(POINTPHI), POINTER :: PHI_LEFTTOPFRONT
        TYPE(POINTPHI), POINTER :: PHI_RIGHTTOPFRONT
        TYPE(POINTPHI), POINTER :: PHI_LEFTBOTTOMFRONT
        TYPE(POINTPHI), POINTER :: PHI_RIGHTBOTTOMFRONT
        TYPE(POINTPHI), POINTER :: PHI_LEFTTOPBACK
        TYPE(POINTPHI), POINTER :: PHI_RIGHTTOPBACK
        TYPE(POINTPHI), POINTER :: PHI_LEFTBOTTOMBACK
        TYPE(POINTPHI), POINTER :: PHI_RIGHTBOTTOMBACK
        
        TYPE(OCTREE), POINTER :: TREE_LEFTTOPFRONT
        TYPE(OCTREE), POINTER :: TREE_RIGHTTOPFRONT
        TYPE(OCTREE), POINTER :: TREE_LEFTBOTTOMFRONT
        TYPE(OCTREE), POINTER :: TREE_RIGHTBOTTOMFRONT
        TYPE(OCTREE), POINTER :: TREE_LEFTTOPBACK
        TYPE(OCTREE), POINTER :: TREE_RIGHTTOPBACK
        TYPE(OCTREE), POINTER :: TREE_LEFTBOTTOMBACK
        TYPE(OCTREE), POINTER :: TREE_RIGHTBOTTOMBACK
        
        TYPE(OCTREE), POINTER :: PARENT
        TYPE(OCTREE), POINTER :: NEXT
        TYPE(OCTREE), POINTER :: BEFORE
    END TYPE
    
    REAL(8) :: TIME_STEP
    INTEGER :: MAX_TREE_DEPTH
    
    TYPE(POINTPHI), POINTER :: FLUID_POINTPHI_HEAD
    TYPE(POINTPHI), POINTER :: FLUID_POINTPHI_TAIL
    TYPE(OCTREE), POINTER :: FLUID_TREE_HEAD
    TYPE(OCTREE), POINTER :: FLUID_TREE_TAIL
    
    TYPE(POINTPHI), POINTER :: PROPEL_POINTPHI_HEAD
    TYPE(POINTPHI), POINTER :: PROPEL_POINTPHI_TAIL
    TYPE(OCTREE), POINTER :: PROPEL_TREE_HEAD
    TYPE(OCTREE), POINTER :: PROPEL_TREE_TAIL
    
    TYPE(POINTPHI), POINTER :: CASE_POINTPHI_HEAD
    TYPE(POINTPHI), POINTER :: CASE_POINTPHI_TAIL
    TYPE(OCTREE), POINTER :: CASE_TREE_HEAD
    TYPE(OCTREE), POINTER :: CASE_TREE_TAIL
    
    REAL(8) :: FLUID_SIZE
    REAL(8), ALLOCATABLE :: FLUID_SURFACE_POINTS
    INTEGER :: FLUID_SURFACE_POINTS_NUM
    REAL(8), ALLOCATABLE :: FLUID_SURFACE_FACES
    INTEGER :: FLUID_SURFACE_FACES_NUM
    
    REAL(8) :: PROPEL_SIZE
    REAL(8), ALLOCATABLE :: PROPEL_SURFACE_POINTS
    INTEGER :: PROPEL_SURFACE_POINTS_NUM
    REAL(8), ALLOCATABLE :: PROPEL_SURFACE_FACES
    INTEGER :: PROPEL_SURFACE_FACES_NUM
    
    REAL(8) :: CASE_SIZE
    REAL(8), ALLOCATABLE :: CASE_SURFACE_POINTS
    INTEGER :: CASE_SURFACE_POINTS_NUM
    REAL(8), ALLOCATABLE :: CASE_SURFACE_FACES
    INTEGER :: CASE_SURFACE_FACES_NUM

    CONTAINS
    
    SUBROUTINE POINTPHI_CONSTRUCT(CURRENT)
        TYPE(POINTPHI) :: CURRENT
        NULLIFY(CURRENT%LEFT)
        NULLIFY(CURRENT%RIGHT)
        NULLIFY(CURRENT%TOP)
        NULLIFY(CURRENT%BOTTOM)
        NULLIFY(CURRENT%FRONT)
        NULLIFY(CURRENT%BACK)
        NULLIFY(CURRENT%NEXT)
        NULLIFY(CURRENT%BEFORE)
    END SUBROUTINE POINTPHI_CONSTRUCT
    
    SUBROUTINE OCTREE_CONSTRUCT1(CURRENT)
        TYPE(OCTREE) :: CURRENT
        
        NULLIFY(CURRENT%PHI_LEFTBOTTOMBACK)
        NULLIFY(CURRENT%PHI_RIGHTBOTTOMBACK)
        NULLIFY(CURRENT%PHI_LEFTTOPBACK)
        NULLIFY(CURRENT%PHI_RIGHTTOPBACK)
        NULLIFY(CURRENT%PHI_LEFTBOTTOMFRONT)
        NULLIFY(CURRENT%PHI_RIGHTBOTTOMFRONT)
        NULLIFY(CURRENT%PHI_LEFTTOPFRONT)
        NULLIFY(CURRENT%PHI_RIGHTTOPFRONT)
        
        NULLIFY(CURRENT%TREE_LEFTBOTTOMBACK)
        NULLIFY(CURRENT%TREE_RIGHTBOTTOMBACK)
        NULLIFY(CURRENT%TREE_LEFTTOPBACK)
        NULLIFY(CURRENT%TREE_RIGHTTOPBACK)
        NULLIFY(CURRENT%TREE_LEFTBOTTOMFRONT)
        NULLIFY(CURRENT%TREE_RIGHTBOTTOMFRONT)
        NULLIFY(CURRENT%TREE_LEFTTOPFRONT)
        NULLIFY(CURRENT%TREE_RIGHTTOPFRONT)
        
        NULLIFY(CURRENT%PARENT)
        
        NULLIFY(CURRENT%NEXT)
        NULLIFY(CURRENT%BEFORE)
        
    END SUBROUTINE OCTREE_CONSTRUCT1
    
    SUBROUTINE OCTREE_CONSTRUCT2(CURRENT, TREE_PARENT, TREE_DEPTH, TREE_LENGTH, TREE_PHI_LEFTBOTTOMBACK, TREE_PHI_LEFTTOPBACK, &
    TREE_PHI_RIGHTBOTTOMBACK, TREE_PHI_RIGHTTOPBACK, TREE_PHI_LEFTBOTTOMFRONT, TREE_PHI_LEFTTOPFRONT, &
    TREE_PHI_RIGHTBOTTOMFRONT, TREE_PHI_RIGHTTOPFRONT)
        TYPE(OCTREE) :: CURRENT
        INTEGER :: TREE_DEPTH
        REAL(8) :: TREE_LENGTH
        TYPE(POINTPHI), POINTER :: TREE_PHI_LEFTBOTTOMBACK, TREE_PHI_LEFTTOPBACK, TREE_PHI_RIGHTBOTTOMBACK, TREE_PHI_RIGHTTOPBACK, &
        TREE_PHI_LEFTBOTTOMFRONT, TREE_PHI_LEFTTOPFRONT, TREE_PHI_RIGHTBOTTOMFRONT, TREE_PHI_RIGHTTOPFRONT
        TYPE(OCTREE), POINTER :: TREE_PARENT
        
        CURRENT%DEPTH = TREE_DEPTH
        CURRENT%LENGTH = TREE_LENGTH
        
        CURRENT%PHI_LEFTBOTTOMBACK => TREE_PHI_LEFTBOTTOMBACK
        CURRENT%PHI_RIGHTBOTTOMBACK => TREE_PHI_RIGHTBOTTOMBACK
        CURRENT%PHI_LEFTTOPBACK => TREE_PHI_LEFTTOPBACK
        CURRENT%PHI_RIGHTTOPBACK => TREE_PHI_RIGHTTOPBACK
        CURRENT%PHI_LEFTBOTTOMFRONT => TREE_PHI_LEFTBOTTOMFRONT
        CURRENT%PHI_RIGHTBOTTOMFRONT => TREE_PHI_RIGHTBOTTOMFRONT
        CURRENT%PHI_LEFTTOPFRONT => TREE_PHI_LEFTTOPFRONT
        CURRENT%PHI_RIGHTTOPFRONT => TREE_PHI_RIGHTTOPFRONT
        
        NULLIFY(CURRENT%TREE_LEFTBOTTOMBACK)
        NULLIFY(CURRENT%TREE_RIGHTBOTTOMBACK)
        NULLIFY(CURRENT%TREE_LEFTTOPBACK)
        NULLIFY(CURRENT%TREE_RIGHTTOPBACK)
        NULLIFY(CURRENT%TREE_LEFTBOTTOMFRONT)
        NULLIFY(CURRENT%TREE_RIGHTBOTTOMFRONT)
        NULLIFY(CURRENT%TREE_LEFTTOPFRONT)
        NULLIFY(CURRENT%TREE_RIGHTTOPFRONT)
        
        CURRENT%PARENT => TREE_PARENT
        
        NULLIFY(CURRENT%NEXT)
        NULLIFY(CURRENT%BEFORE)
        
    END SUBROUTINE OCTREE_CONSTRUCT2
    
    SUBROUTINE DELETE_ALL_POINTPHI(TAIL)
        TYPE(POINTPHI), POINTER :: TAIL, TEMP, TEMP2
        TEMP => TAIL
        DO WHILE(ASSOCIATED(TEMP))
            TEMP2 => TEMP%BEFORE
            DEALLOCATE(TEMP)
            TEMP => TEMP2
        END DO
    END SUBROUTINE DELETE_ALL_POINTPHI
    
    RECURSIVE SUBROUTINE DELETE_ALL_OCTREE(CURRENT)
        TYPE(OCTREE), POINTER :: CURRENT
        IF(ASSOCIATED(CURRENT)) THEN
            CALL DELETE_ALL_OCTREE(CURRENT%TREE_LEFTBOTTOMBACK)
            CALL DELETE_ALL_OCTREE(CURRENT%TREE_RIGHTBOTTOMBACK)
            CALL DELETE_ALL_OCTREE(CURRENT%TREE_LEFTTOPBACK)
            CALL DELETE_ALL_OCTREE(CURRENT%TREE_RIGHTTOPBACK)
            CALL DELETE_ALL_OCTREE(CURRENT%TREE_LEFTBOTTOMFRONT)
            CALL DELETE_ALL_OCTREE(CURRENT%TREE_RIGHTBOTTOMFRONT)
            CALL DELETE_ALL_OCTREE(CURRENT%TREE_LEFTTOPFRONT)
            CALL DELETE_ALL_OCTREE(CURRENT%TREE_RIGHTTOPFRONT)
            
            DEALLOCATE(CURRENT)
        END IF
    END SUBROUTINE DELETE_ALL_OCTREE
    
    SUBROUTINE FIND_LR_BT_BF_FACE(CURRENT, X, Y, Z, LR, BT, BF,         TEMP)
        TYPE(OCTREE), POINTER :: CURRENT, TEMP
        REAL(8) :: X, Y, Z, X_TEMP, Y_TEMP, Z_TEMP, X_ERROR, Y_ERROR, Z_ERROR
        LOGICAL :: LR, BT, BF
        
        TEMP => CURRENT
        
        DO WHILE(ASSOCIATED(TEMP%TREE_LEFTBOTTOMBACK))
            X_TEMP = TEMP%PHI_LEFTBOTTOMBACK%X + TEMP%LENGTH/2.
            Y_TEMP = TEMP%PHI_LEFTBOTTOMBACK%Y + TEMP%LENGTH/2.
            Z_TEMP = TEMP%PHI_LEFTBOTTOMBACK%Z + TEMP%LENGTH/2.
            
            IF(LR) THEN
                X_ERROR = MINERROR
            ELSE
                X_ERROR = -MINERROR
            END IF
        
            IF(BT) THEN
                Y_ERROR = MINERROR
            ELSE
                Y_ERROR = -MINERROR
            END IF
        
            IF(BF) THEN
                Z_ERROR = MINERROR
            ELSE
                Z_ERROR = -MINERROR
            END IF
        
		    IF(X <= X_TEMP + X_ERROR) THEN
			    IF(Y <= Y_TEMP + Y_ERROR) THEN
				    IF(Z <= Z_TEMP + Z_ERROR) THEN
                        TEMP => TEMP%TREE_LEFTBOTTOMBACK
                    ELSE
                        TEMP => TEMP%TREE_LEFTBOTTOMFRONT
                    END IF
			    ELSE
				    IF(Z <= Z_TEMP + Z_ERROR) THEN
                        TEMP => TEMP%TREE_LEFTTOPBACK
                    ELSE
                        TEMP => TEMP%TREE_LEFTTOPFRONT
                    END IF
                END IF
            
		    ELSE
			    IF(Y <= Y_TEMP + Y_ERROR) THEN
				    IF(Z <= Z_TEMP + Z_ERROR) THEN
                        TEMP => TEMP%TREE_RIGHTBOTTOMBACK
                    ELSE
                        TEMP => TEMP%TREE_RIGHTBOTTOMFRONT
                    END IF
			    ELSE
				    IF(Z <= Z_TEMP + Z_ERROR) THEN
                        TEMP => TEMP%TREE_RIGHTTOPBACK
                    ELSE
                        TEMP => TEMP%TREE_RIGHTTOPFRONT
                    END IF
                END IF
            END IF
        END DO
    END SUBROUTINE FIND_LR_BT_BF_FACE
    
    
    SUBROUTINE FIND_VERTEX(CURRENT, DIRECTION, X, Y, Z,         TEMP)
        TYPE(POINTPHI), POINTER :: CURRENT, TEMP
        INTEGER :: DIRECTION
        REAL(8) :: X,Y,Z
        
        TEMP => CURRENT
        IF(DIRECTION == 0) THEN
            DO WHILE(TEMP%X > X + MINERROR)
                TEMP => TEMP%LEFT
            END DO
            IF(TEMP%X > X - MINERROR) THEN
                RETURN
            ELSE
                NULLIFY(TEMP)
                RETURN
            END IF
        ELSE IF(DIRECTION == 1) THEN
            DO WHILE(TEMP%Y > Y + MINERROR)
                TEMP => TEMP%BOTTOM
            END DO
            IF(TEMP%Y > Y - MINERROR) THEN
                RETURN
            ELSE
                NULLIFY(TEMP)
                RETURN
            END IF
        ELSE IF(DIRECTION == 2) THEN
            DO WHILE(TEMP%Z > Z + MINERROR)
                TEMP => TEMP%BACK
            END DO
            IF(TEMP%Z > Z - MINERROR) THEN
                RETURN
            ELSE
                NULLIFY(TEMP)
                RETURN
            END IF
        ELSE IF(DIRECTION == 3) THEN
            DO WHILE(TEMP%X > X + MINERROR)
                TEMP => TEMP%LEFT
            END DO
            IF(TEMP%X > X - MINERROR) THEN
                DO WHILE(ASSOCIATED(TEMP) .AND. TEMP%Y > Y + MINERROR)
                    TEMP => TEMP%BOTTOM
                END DO
                IF(.NOT. ASSOCIATED(TEMP) .OR. TEMP%Y > Y - MINERROR) THEN
                    RETURN
                ELSE
                    NULLIFY(TEMP)
                    RETURN
                END IF
            ELSE
                NULLIFY(TEMP)
                RETURN
            END IF
        ELSE IF(DIRECTION == 4) THEN
            DO WHILE(TEMP%X > X + MINERROR)
                TEMP => TEMP%LEFT
            END DO
            IF(TEMP%X > X - MINERROR) THEN
                DO WHILE(ASSOCIATED(TEMP) .AND. TEMP%Z > Z + MINERROR)
                    TEMP => TEMP%BACK
                END DO
                IF(.NOT. ASSOCIATED(TEMP) .OR. TEMP%Z > Z - MINERROR) THEN
                    RETURN
                ELSE
                    NULLIFY(TEMP)
                    RETURN
                END IF
            ELSE
                NULLIFY(TEMP)
                RETURN
            END IF
        ELSE IF(DIRECTION == 5) THEN
            DO WHILE(TEMP%Y > Y + MINERROR)
                TEMP => TEMP%BOTTOM
            END DO
            IF(TEMP%Y > Y - MINERROR) THEN
                DO WHILE(ASSOCIATED(TEMP) .AND. TEMP%Z > Z + MINERROR)
                    TEMP => TEMP%BACK
                END DO
                IF(.NOT. ASSOCIATED(TEMP) .OR. TEMP%Z > Z - MINERROR) THEN
                    RETURN
                ELSE
                    NULLIFY(TEMP)
                    RETURN
                END IF
            ELSE
                NULLIFY(TEMP)
                RETURN
            END IF
        END IF
        NULLIFY(TEMP)
    END SUBROUTINE FIND_VERTEX
        
END MODULE SURFACE_MODULE
    
PROGRAM MAIN
USE SURFACE_MODULE
USE PROPA_RECONST_REINITIAL
    IMPLICIT NONE
    INTEGER :: POINT_NUM
    REAL(8), ALLOCATABLE :: POINT(:,:) ! POINT(DIMENSION,POINT_NUM)
    INTEGER :: SURF_NUM
    INTEGER, ALLOCATABLE :: SURF(:,:) ! SURF(3,SURF_NUM)
    
    REAL(8) :: X,Y,Z
    REAL(8) :: RET
    
    INTEGER :: DEPTH
    REAL(8) :: TIMESTEP
    
    INTEGER :: I,J,K,N, NUM
    REAL(8), ALLOCATABLE :: SAVING(:)
    
    REAL(8) :: BRATE(0)
    
    POINT_NUM = 8
    ALLOCATE(POINT(3,POINT_NUM))
    SURF_NUM = 12
    ALLOCATE(SURF(3,SURF_NUM))
    
	POINT(1,1) = 0.2;
	POINT(2,1) = 0.2;
	POINT(3,1) = 0.2;
		
	POINT(1,2) = 0.8;
	POINT(2,2) = 0.2;
	POINT(3,2) = 0.2;

	POINT(1,3) = 0.2;
	POINT(2,3) = 0.8;
	POINT(3,3) = 0.2;
		
	POINT(1,4) = 0.8;
	POINT(2,4) = 0.8;
	POINT(3,4) = 0.2;

	POINT(1,5) = 0.2;
	POINT(2,5) = 0.2;
	POINT(3,5) = 0.8;
		
	POINT(1,6) = 0.8;
	POINT(2,6) = 0.2;
	POINT(3,6) = 0.8;
		
	POINT(1,7) = 0.2;
	POINT(2,7) = 0.8;
	POINT(3,7) = 0.8;

	POINT(1,8) = 0.8;
	POINT(2,8) = 0.8;
	POINT(3,8) = 0.8;

	SURF(1,1) = 1;
	SURF(2,1) = 3;
	SURF(3,1) = 2;
		
	SURF(1,2) = 2;
	SURF(2,2) = 3;
	SURF(3,2) = 4;
		
	SURF(1,3) = 1;
	SURF(2,3) = 2;
	SURF(3,3) = 5;
		
	SURF(1,4) = 2;
	SURF(2,4) = 6;
	SURF(3,4) = 5;
    
	SURF(1,5) = 2;
	SURF(2,5) = 4;
	SURF(3,5) = 6;
		
	SURF(1,6) = 4;
	SURF(2,6) = 8;
	SURF(3,6) = 6;
		
	SURF(1,7) = 3;
	SURF(2,7) = 7;
	SURF(3,7) = 8;
		
	SURF(1,8) = 3;
	SURF(2,8) = 8;
	SURF(3,8) = 4;
		
	SURF(1,9) = 1;
	SURF(2,9) = 5;
	SURF(3,9) = 3;
		
	SURF(1,10) = 5;
	SURF(2,10) = 7;
	SURF(3,10) = 3;
		
	SURF(1,11) = 6;
	SURF(2,11) = 8;
	SURF(3,11) = 7;
		
	SURF(1,12) = 7;
	SURF(2,12) = 5;
	SURF(3,12) = 6;
    
    X = 0.25
    Y = 0.25
    Z = 0.25
    
    !CALL DISTANCE_SURF_POINT(POINT_NUM, POINT, SURF_NUM, SURF, X, Y, Z, RET)
    !CALL DISTANCE_SURF_POINT(POINT_NUM, POINT, SURF_NUM, SURF, X, Y, Z, RET)
    !CALL DISTANCE_SURF_POINT(POINT_NUM, POINT, SURF_NUM, SURF, X, Y, Z, RET)
    !CALL DISTANCE_SURF_POINT(POINT_NUM, POINT, SURF_NUM, SURF, X, Y, Z, RET)
    !CALL DISTANCE_SURF_POINT(POINT_NUM, POINT, SURF_NUM, SURF, X, Y, Z, RET)
    !CALL DISTANCE_SURF_POINT(POINT_NUM, POINT, SURF_NUM, SURF, X, Y, Z, RET)
    !CALL DISTANCE_SURF_POINT(POINT_NUM, POINT, SURF_NUM, SURF, X, Y, Z, RET)
    !CALL DISTANCE_SURF_POINT(POINT_NUM, POINT, SURF_NUM, SURF, X, Y, Z, RET)
    !CALL DISTANCE_SURF_POINT(POINT_NUM, POINT, SURF_NUM, SURF, X, Y, Z, RET)
    
    PRINT *, RET
    
    DEPTH = 6
    TIMESTEP = 0.1
    
    NUM = 2**DEPTH
    
    
    CALL SURFACE_CONSTRUCT(POINT_NUM, POINT, SURF_NUM, SURF, DEPTH, TIMESTEP)
    
    OPEN(UNIT=21, FILE = "tree_3dlevel0.txt", ACTION = "WRITE", STATUS = "REPLACE")
    
    ALLOCATE(SAVING((NUM+1)*(NUM+1)*(NUM+1)))
    N = 1
    DO I = 0,NUM
        DO J = 0,NUM
            DO K = 0,NUM
                X = MAXSIZE * I/NUM
                Y = MAXSIZE * J/NUM
                Z = MAXSIZE * K/NUM
                IF(N==32) THEN
                    DEPTH = 5
                END IF
                CALL LEVEL_COMPUTEPHI_TYPE(X,Y,Z, 0, SAVING(N))
                N = N + 1
            END DO
        END DO
    END DO
    
    DO I = 1,(NUM+1)*(NUM+1)*(NUM+1)
        WRITE(21,*) SAVING(I)
    END DO
    
    DEALLOCATE(SAVING)
    
    CLOSE(21)
    
    
    !FLUID_SURFACE_POINTS_NUM = POINT_NUM
    !ALLOCATE(FLUID_SURFACE_POINTS(POINT_NUM))
    !FLUID_SURFACE_POINTS = POINT
    !FLUID_SURFACE_FACES_NUM = SURF_NUM
    !ALLOCATE(FLUID_SURFACE_FACES(SURF_NUM))
    !FLUID_SURFACE_FACES = SURF
    
    CALL TREE_PROPAGATE_SURFACE(BRATE, 0)
    
    OPEN(UNIT=21, FILE = "tree_3dlevel1.txt", ACTION = "WRITE", STATUS = "REPLACE")
    
    ALLOCATE(SAVING((NUM+1)*(NUM+1)*(NUM+1)))
    N = 1
    DO I = 0,NUM
        DO J = 0,NUM
            DO K = 0,NUM
                X = MAXSIZE * I/NUM
                Y = MAXSIZE * J/NUM
                Z = MAXSIZE * K/NUM
                IF(N==32) THEN
                    DEPTH = 5
                END IF
                CALL LEVEL_COMPUTEPHI_TYPE(X,Y,Z, 0, SAVING(N))
                N = N + 1
            END DO
        END DO
    END DO
    
    DO I = 1,(NUM+1)*(NUM+1)*(NUM+1)
        WRITE(21,*) SAVING(I)
    END DO
    
    DEALLOCATE(SAVING)
    
    CLOSE(21)
    
    !DEALLOCATE(FLUID_SURFACE_POINTS)
    !DEALLOCATE(FLUID_SURFACE_FACES)
    
    
END PROGRAM MAIN
    