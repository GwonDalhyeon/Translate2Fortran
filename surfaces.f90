    SUBROUTINE UNSIGNED_DISTANCE_FACE_POINT(POINT_NUM, POINT, SURF_NUM, SURF,I0,V,TEMPDIST)
    
        INTEGER :: I0
        REAL(8) :: V(3) 
    
        INTEGER :: IMIN
        REAL(8) :: TEMPDIST, POINTDIST
        REAL(8) :: A1(3),A2(3),A3(3)
        REAL(8) :: N(3), ZEROVEC(3), C1(3), C2(3)
        REAL(8) :: ROT1, ROT2, INNER1, INNER2
    
        INTEGER :: POINT_NUM
        REAL(8) :: POINT(3,POINT_NUM) ! POINT(DIMENSION,POINT_NUM)
        INTEGER :: SURF_NUM
        INTEGER :: SURF(3,SURF_NUM) ! SURF(3,SURF_NUM)
    
        IMIN = 1
        
        DO I=1,3
        
        TEMPDIST = SQRT( (V(1) - POINT(1,SURF(I,I0)) ) **2 + ( V(2)- POINT(2,SURF(I,I0)) ) **2 + ( V(3)- POINT(3,SURF(I,I0)) ) **2 )
        
        IF (I==1) THEN 
           
            POINTDIST = TEMPDIST
            
        ELSEIF(TEMPDIST < POINTDIST) THEN
            
            POINTDIST = TEMPDIST
            IMIN = I
        
        END IF
    
        END DO     
            
        IF(IMIN == 1) THEN
        
            DO I=1,3
            
                A1(I) = POINT(I,SURF(1,I0))
                A2(I) = POINT(I,SURF(2,I0))
                A3(I) = POINT(I,SURF(3,I0))
        
            END DO
        
        ELSEIF(IMIN == 2) THEN
        
            DO I=1,3
            
                A1(I) = POINT(I,SURF(2,I0))
                A2(I) = POINT(I,SURF(3,I0))
                A3(I) = POINT(I,SURF(1,I0))
        
            END DO
        
        ELSE
        
            DO I=1,3
            
                A1(I) = POINT(I,SURF(3,I0))
                A2(I) = POINT(I,SURF(1,I0))
                A3(I) = POINT(I,SURF(2,I0))
        
            END DO
            
        END IF
        
        CALL VEC_CURL2(A1,A2,A1,A3,N)
        
        ZEROVEC(1) = 0
        ZEROVEC(2) = 0
        ZEROVEC(3) = 0
    
        CALL VEC_CURL2(V,A1,A1,A2,C1)
        CALL VEC_CURL2(V,A1,A1,A3,C2)
        
        CALL VEC_INNER_PRODUCT1(C1,N,ROT1)
        CALL VEC_INNER_PRODUCT1(C2,N,ROT2)
        
        CALL VEC_INNER_PRODUCT2(V,A1,A1,A2,INNER1)
        CALL VEC_INNER_PRODUCT2(V,A1,A1,A3,INNER2)
    
    
        IF(INNER1>0.AND.INNER2>0) THEN
            TEMPDIST = POINTDIST
        ELSE IF(ROT1>0.AND.ROT2<0) THEN
            CALL DISTANCE_FACE_POINT(V,A1,A2,A3,TEMPDIST)
        ELSE IF(ROT1<0) THEN
            CALL DISTANCE_LINE_POINT(V,A1,A2,TEMPDIST)
        ELSE
            CALL DISTANCE_LINE_POINT(V,A1,A3,TEMPDIST)
        END IF
    
    END SUBROUTINE UNSIGNED_DISTANCE_FACE_POINT
    
    
    
    
    SUBROUTINE LINE_FACE_INTERSECTING(POINT_NUM, POINT, SURF_NUM, SURF, I0, L, V, SGN)
    
        INTEGER :: POINT_NUM
        REAL(8) :: POINT(3,POINT_NUM) ! POINT(DIMENSION,POINT_NUM)
        INTEGER :: SURF_NUM
        INTEGER :: SURF(3,SURF_NUM) ! SURF(3,SURF_NUM)
        INTEGER :: I0
        REAL(8) :: L(3), V(3) 
        REAL(8) :: A(3),B(3),C(3)
        REAL(8) :: M(3)
        REAL(8) :: ZEROVEC(3)
        REAL(8) :: L1(3), L2(3)
        REAL(8) :: A1,A2,B1,B2,C1,C2
        REAL(8) :: D1,D2,D3
        REAL(8) :: BC
        INTEGER :: SGN
        REAL(8) :: R,S,T
        
        REAL(8) :: TEMP
    
        ZEROVEC(1) = 0.
        ZEROVEC(2) = 0.
        ZEROVEC(3) = 0.
        
        DO I=1,3
            A(I) = POINT(I,SURF(1,I0))
            B(I) = POINT(I,SURF(2,I0))
            C(I) = POINT(I,SURF(3,I0))
        END DO
    
        CALL RANDOM_NUMBER(M)
        CALL VEC_CURL1(L,M,L1)
        CALL VEC_NORM(L1, TEMP)
        IF(TEMP<MINERROR) THEN
            SGN = -1
            RETURN
        END IF
        
        CALL VEC_NORMALIZATION(L1)
        CALL VEC_CURL1(L,L1,L2)
    
        CALL VEC_INNER_PRODUCT2(A,V,ZEROVEC,L1,A1)
        CALL VEC_INNER_PRODUCT2(A,V,ZEROVEC,L2,A2)
        CALL VEC_INNER_PRODUCT2(A,B,ZEROVEC,L1,B1)
        CALL VEC_INNER_PRODUCT2(A,B,ZEROVEC,L2,B2)
        CALL VEC_INNER_PRODUCT2(A,C,ZEROVEC,L1,C1)
        CALL VEC_INNER_PRODUCT2(A,C,ZEROVEC,L2,C2)
    
        CALL VEC_INNER_PRODUCT2(A,V,ZEROVEC,L,D1)
        CALL VEC_INNER_PRODUCT2(A,B,ZEROVEC,L,D2)
        CALL VEC_INNER_PRODUCT2(A,C,ZEROVEC,L,D3)
    
        BC = B1*C2 - B2*C1 
    
        IF ( ABS(BC) < MINERROR) THEN
            SGN = -1
            RETURN
        END IF
    
        S = (A1*C2 - A2*C1)/BC
        R = (A1*B2 - A2*B1)/(-BC)
        T = S*D2 + R*D3 - D1
    
        IF(T>0.AND.S>0.AND.R>0.AND.S+R<1) THEN
            SGN = 1
        ELSEIF(ABS(S)<MINERROR .OR. ABS(R)<MINERROR .OR. ABS(S+R)<MINERROR) THEN
            SGN = -1
        ELSE 
            SGN = 0
        END IF
    
    END SUBROUTINE LINE_FACE_INTERSECTING
    
    
    
    
    
    
    SUBROUTINE LINE_SURFACE_INTERSECTING(POINT_NUM, POINT, SURF_NUM, SURF,V,RETURNS)

        INTEGER :: POINT_NUM
        REAL(8) :: POINT(3,POINT_NUM) ! POINT(DIMENSION,POINT_NUM)
        INTEGER :: SURF_NUM
        INTEGER :: SURF(3,SURF_NUM) ! SURF(3,SURF_NUM)
        INTEGER :: NUM
        INTEGER :: SGN
        REAL(8) :: L(3), V(3)
        REAL(8) :: RETURNS
        
        SGN = -1
        
        DO WHILE(SGN == -1)

            NUM = 0
            CALL RANDOM_NUMBER(L)
            CALL VEC_NORMALIZATION(L)

            DO I=1,SURF_NUM
                CALL LINE_FACE_INTERSECTING(POINT_NUM, POINT, SURF_NUM, SURF,I,L,V,SGN)

                IF (SGN == -1) THEN

                    EXIT

                END IF

                IF (SGN == 1) THEN

                    NUM = NUM+1

                END IF
            END DO
        END DO

        IF(MOD(NUM,2) == 0 ) THEN

            RETURNS = 1

        ELSE

            RETURNS = -1

        END IF

    END SUBROUTINE LINE_SURFACE_INTERSECTING
    
    

    SUBROUTINE DISTANCE_SURF_POINT(POINT_NUM, POINT, SURF_NUM, SURF, X, Y, Z, DISTFINAL)

        INTEGER :: POINT_NUM
        REAL(8) :: POINT(3,POINT_NUM) ! POINT(DIMENSION,POINT_NUM)
        INTEGER :: SURF_NUM
        INTEGER :: SURF(3,SURF_NUM) ! SURF(3,SURF_NUM)
        REAL(8) :: DIST, TEMPDIST, DISTFINAL
        REAL(8) :: VEC(3)
        REAL(8) :: X,Y,Z, TEMP
  
        VEC(1) = X
        VEC(2) = Y
        VEC(3) = Z
  
        DO I=1,SURF_NUM
            CALL UNSIGNED_DISTANCE_FACE_POINT(POINT_NUM, POINT, SURF_NUM, SURF,I,VEC,TEMPDIST)
      
            IF (I==1) THEN
                DIST = TEMPDIST
            ELSEIF(TEMPDIST <= DIST) THEN
                DIST = TEMPDIST
            END IF
        END DO 
  
        IF (DIST==0) THEN
            DISTFINAL = 0
        END IF
  
        CALL LINE_SURFACE_INTERSECTING(POINT_NUM, POINT, SURF_NUM, SURF,VEC,TEMP)
        DISTFINAL = TEMP * DIST
  
    END SUBROUTINE DISTANCE_SURF_POINT