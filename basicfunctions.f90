    SUBROUTINE DISTANCE_FACE_POINT(V,W1,W2,W3,     D)
    
        REAL(8) :: V(3),W1(3),W2(3),W3(3),S(3),VEC(3)
        REAL(8) :: U,L2,T,NORM,INNER,D
    
        CALL VEC_CURL1(W2-W1,W3-W1,S)
        CALL VEC_NORM(S,T)
        CALL VEC_INNER_PRODUCT1(S,V-W1,U)
        D = ABS(U)/T
    
    END SUBROUTINE
    
    SUBROUTINE DISTANCE_LINE_POINT(V,W1,W2,     D) 
        REAL(8) :: V(3),W1(3),W2(3),VEC(3)
        REAL(8) :: L2,INNER,D,NORM
    
        CALL VEC_NORM(W2-W1,NORM)
        L2 = NORM**2
        CALL VEC_INNER_PRODUCT1(W1-V,W2-W1,INNER)
        DO I=1,3
           VEC(I) = W1(I)- V(I) - INNER/L2 * (W2(I)-W1(I))
        END DO
        CALL VEC_NORM(VEC,D)
    
    END SUBROUTINE DISTANCE_LINE_POINT
    
    SUBROUTINE MINMOD(A,B,      T)
        REAL(8) :: A,B,T
        INTEGER :: S1,S2
    
        CALL SIGN1(A,S1)
        CALL SIGN1(B,S2)
        T = (S1+S2)/2 * MIN(ABS(A),ABS(B))
        
    END SUBROUTINE MINMOD
    
    SUBROUTINE SIGN1(A,       T)
        REAL(8) :: A 
        INTEGER :: T
    
        IF(A.LT.0.0) THEN 
            T=-1
    
        ELSE IF(A.EQ.0.0) THEN
            T=0
  
        ELSE 
            T=1
        END IF
  
    END SUBROUTINE SIGN1
    
    SUBROUTINE SIGN2(R,DX,       T)
        REAL(8) :: R,DX
        REAL(8) :: T
        
        T = R/SQRT(R**2 + DX**2)
  
    END SUBROUTINE SIGN2
    
    SUBROUTINE VEC_CURL1(V,W,       R) ! DIM=3 
        REAL(8) :: V(3),W(3),R(3)

        R(1) = V(2)*W(3) - V(3)*W(2)
        R(2) = V(3)*W(1) - V(1)*W(3)
        R(3) = V(1)*W(2) - V(2)*W(1)

    END SUBROUTINE VEC_CURL1
    
    SUBROUTINE VEC_CURL2(V1,V2,W1,W2,       R) ! DIM=3 
        REAL(8) :: V1(3),V2(3),W1(3),W2(3),R(3)

        R(1) = (V2(2)-V1(2))*(W2(3)-W1(3)) - (V2(3)-V1(3))*(W2(2)-W1(2))
        R(2) = (V2(3)-V1(3))*(W2(1)-W1(1)) - (V2(1)-V1(1))*(W2(3)-W1(3))
        R(3) = (V2(1)-V1(1))*(W2(2)-W1(2)) - (V2(2)-V1(2))*(W2(1)-W1(1))

    END SUBROUTINE VEC_CURL2
    
    SUBROUTINE VEC_INNER_PRODUCT1(V1,W1,        R) 
        REAL(8) :: R
        REAL(8) :: V1(3)
        REAL(8) :: W1(3)
    
        R=0.0
        DO I=1,3
            R = R+V1(I)*W1(I)
        END DO
        
    END SUBROUTINE VEC_INNER_PRODUCT1
    
    SUBROUTINE VEC_INNER_PRODUCT2(V1,V2,W1,W2,      R) 
        REAL(8) :: R
        REAL(8) :: V1(3)
        REAL(8) :: V2(3)
        REAL(8) :: W1(3)
        REAL(8) :: W2(3)
    
        R=0.0
        DO I=1,3
            R = R+(V2(I)-V1(I))*(W2(I)-W1(I))
        END DO
    
    END SUBROUTINE VEC_INNER_PRODUCT2
    
    SUBROUTINE VEC_NORM(W,      R) 
        REAL(8) :: R
        REAL(8) :: W(3)
   
        R=0.0
        DO I=1,3
            R=R+W(I)*W(I)
        END DO
    
        R=SQRT(R)
    
    END SUBROUTINE VEC_NORM
    
    SUBROUTINE VEC_NORMALIZATION(W)
        REAL(8) :: W(3)
        REAL(8) :: S1
    
        CALL VEC_NORM(W,S1)
        DO I=1,3
        W(I) = W(I)/S1
        END DO
    
    END SUBROUTINE VEC_NORMALIZATION
    
    SUBROUTINE INIT_RANDOM_SEED()

        IMPLICIT NONE
        
        INTEGER :: I, N, CLOCK
        INTEGER, DIMENSION(:), ALLOCATABLE :: SEED
        CALL RANDOM_SEED(SIZE = N)
        ALLOCATE(SEED(N))
        CALL SYSTEM_CLOCK(COUNT=CLOCK)
        SEED = CLOCK + 37*(/(I-1, I=1, N)/)
        CALL RANDOM_SEED(PUT = SEED)
        DEALLOCATE(SEED)

    END SUBROUTINE INIT_RANDOM_SEED