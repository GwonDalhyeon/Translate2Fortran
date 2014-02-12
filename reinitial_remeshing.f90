MODULE REINITIAL_REMESHING
    USE SURFACE_MODULE
    USE OPERATORS
    USE MAKING
    
    CONTAINS
    
    
    SUBROUTINE LEVEL_COMPUTENORMAL_TYPE(X,Y,Z,TYP,      NORMAL)
    
        REAL(8) :: X,Y,Z
        INTEGER :: TYP
        TYPE(OCTREE), POINTER :: CURRENT
        REAL(8) :: X1,X2,Y1,Y2,Z1,Z2,L
        REAL(8) :: PHI_X, PHI_Y, PHI_Z, N
        REAL(8) :: NORMAL(3)
        
        REAL(8) :: PHI_X1, PHI_X2, PHI_X3, PHI_X4, PHI_X5, PHI_X6, PHI_X7, PHI_X8
        REAL(8) :: PHI_Y1, PHI_Y2, PHI_Y3, PHI_Y4, PHI_Y5, PHI_Y6, PHI_Y7, PHI_Y8
        REAL(8) :: PHI_Z1, PHI_Z2, PHI_Z3, PHI_Z4, PHI_Z5, PHI_Z6, PHI_Z7, PHI_Z8
    
        CALL FIND_CELL_CONTAINING_POINT_TYPE(X,Y,Z,TYP,CURRENT)
    
        X1 = CURRENT%PHI_LEFTBOTTOMBACK%X
        X2 = CURRENT%PHI_RIGHTBOTTOMBACK%X
        Y1 = CURRENT%PHI_LEFTBOTTOMBACK%Y
        Y2 = CURRENT%PHI_LEFTTOPBACK%Y
        Z1 = CURRENT%PHI_LEFTBOTTOMBACK%Z
        Z2 = CURRENT%PHI_LEFTBOTTOMFRONT%Z
        L = CURRENT%LENGTH    
    
        CALL OCTREE_DX(CURRENT%PHI_LEFTBOTTOMBACK,PHI_X1)  ! octree_dx(current->phi_leftbottomback)*(x2-x)/l*(y2-y)/l*(z2-z)/l
        CALL OCTREE_DX(CURRENT%PHI_LEFTTOPBACK,PHI_X2) ! octree_dx(current->phi_lefttopback)*(x2-x)/l*(y-y1)/l*(z2-z)/l
        CALL OCTREE_DX(CURRENT%PHI_RIGHTBOTTOMBACK,PHI_X3) ! octree_dx(current->phi_rightbottomback)*(x-x1)/l*(y2-y)/l*(z2-z)/l
        CALL OCTREE_DX(CURRENT%PHI_RIGHTTOPBACK,PHI_X4)  ! octree_dx(current->phi_righttopback)*(x-x1)/l*(y-y1)/l*(z2-z)/l
        CALL OCTREE_DX(CURRENT%PHI_LEFTBOTTOMFRONT,PHI_X5)  ! octree_dx(current->phi_leftbottomfront)*(x2-x)/l*(y2-y)/l*(z-z1)/l    
        CALL OCTREE_DX(CURRENT%PHI_LEFTTOPFRONT,PHI_X6) ! octree_dx(current->phi_lefttopfront)*(x2-x)/l*(y-y1)/l*(z-z1)/l
        CALL OCTREE_DX(CURRENT%PHI_RIGHTBOTTOMFRONT,PHI_X7)  ! octree_dx(current->phi_rightbottomfront)*(x-x1)/l*(y2-y)/l*(z-z1)/l
        CALL OCTREE_DX(CURRENT%PHI_RIGHTTOPFRONT,PHI_X8)  ! octree_dx(current->phi_righttopfront)*(x-x1)/l*(y-y1)/l*(z-z1)/l;
    
        PHI_X = PHI_X1 * (X2-X)/L * (Y2-Y)/L * (Z2-Z)/L + PHI_X2 * (X2-X)/L * (Y-Y1)/L * (Z2-Z)/L + PHI_X3 * (X-X1)/L * (Y2-Y)/L * (Z2-Z)/L &
        + PHI_X4 * (X-X1)/L * (Y-Y1)/L * (Z2-Z)/L+ PHI_X5 * (X2-X)/L * (Y2-Y)/L * (Z-Z1)/L + PHI_X6 *  (X2-X)/L * (Y-Y1)/L * (Z-Z1)/L &
        + PHI_X7 * (X-X1)/L * (Y2-Y)/L * (Z-Z1)/L + PHI_X8 * (X-X1)/L * (Y-Y1)/L * (Z-Z1)/L
    
    
        CALL OCTREE_DY(CURRENT%PHI_LEFTBOTTOMBACK,PHI_Y1) 
        CALL OCTREE_DY(CURRENT%PHI_LEFTTOPBACK,PHI_Y2) 
        CALL OCTREE_DY(CURRENT%PHI_RIGHTBOTTOMBACK,PHI_Y3) 
        CALL OCTREE_DY(CURRENT%PHI_RIGHTTOPBACK,PHI_Y4)
        CALL OCTREE_DY(CURRENT%PHI_LEFTBOTTOMFRONT,PHI_Y5)    
        CALL OCTREE_DY(CURRENT%PHI_LEFTTOPFRONT,PHI_Y6) 
        CALL OCTREE_DY(CURRENT%PHI_RIGHTBOTTOMFRONT,PHI_Y7) 
        CALL OCTREE_DY(CURRENT%PHI_RIGHTTOPFRONT,PHI_Y8) 
    
        PHI_Y = PHI_Y1 * (X2-X)/L * (Y2-Y)/L * (Z2-Z)/L + PHI_Y2 * (X2-X)/L * (Y-Y1)/L * (Z2-Z)/L + PHI_Y3 * (X-X1)/L * (Y2-Y)/L * (Z2-Z)/L &
        + PHI_Y4 * (X-X1)/L * (Y-Y1)/L * (Z2-Z)/L+ PHI_Y5 * (X2-X)/L * (Y2-Y)/L * (Z-Z1)/L + PHI_Y6 *  (X2-X)/L * (Y-Y1)/L * (Z-Z1)/L &
        + PHI_Y7 * (X-X1)/L * (Y2-Y)/L * (Z-Z1)/L + PHI_Y8 * (X-X1)/L * (Y-Y1)/L * (Z-Z1)/L
    
        CALL OCTREE_DZ(CURRENT%PHI_LEFTBOTTOMBACK,PHI_Z1) ! octree_dz(current->phi_leftbottomback)*(x2-x)/l*(y2-y)/l*(z2-z)/l
        CALL OCTREE_DZ(CURRENT%PHI_LEFTTOPBACK,PHI_Z2) ! octree_dz(current->phi_lefttopback)*(x2-x)/l*(y-y1)/l*(z2-z)/l
        CALL OCTREE_DZ(CURRENT%PHI_RIGHTBOTTOMBACK,PHI_Z3) !octree_dz(current->phi_rightbottomback)*(x-x1)/l*(y2-y)/l*(z2-z)/l
        CALL OCTREE_DZ(CURRENT%PHI_RIGHTTOPBACK,PHI_Z4) ! octree_dz(current->phi_righttopback)*(x-x1)/l*(y-y1)/l*(z2-z)/l
        CALL OCTREE_DZ(CURRENT%PHI_LEFTBOTTOMFRONT,PHI_Z5) ! octree_dz(current->phi_leftbottomfront)*(x2-x)/l*(y2-y)/l*(z-z1)/l    
        CALL OCTREE_DZ(CURRENT%PHI_LEFTTOPFRONT,PHI_Z6) ! octree_dy(current->phi_lefttopfront)*(x2-x)/l*(y-y1)/l*(z-z1)/l
        CALL OCTREE_DZ(CURRENT%PHI_RIGHTBOTTOMFRONT,PHI_Z7) ! octree_dy(current->phi_rightbottomfront)*(x-x1)/l*(y2-y)/l*(z-z1)/l
        CALL OCTREE_DZ(CURRENT%PHI_RIGHTTOPFRONT,PHI_Z8) ! octree_dy(current->phi_righttopfront)*(x-x1)/l*(y-y1)/l*(z-z1)/l
    
        PHI_Z = PHI_Z1 * (X2-X)/L * (Y2-Y)/L * (Z2-Z)/L + PHI_Z2 * (X2-X)/L * (Y-Y1)/L * (Z2-Z)/L + PHI_Z3 * (X-X1)/L * (Y2-Y)/L * (Z2-Z)/L &
        + PHI_Z4 * (X-X1)/L * (Y-Y1)/L * (Z2-Z)/L+ PHI_Z5 * (X2-X)/L * (Y2-Y)/L * (Z-Z1)/L + PHI_Z6 *  (X2-X)/L * (Y-Y1)/L * (Z-Z1)/L &
        + PHI_Z7 * (X-X1)/L * (Y2-Y)/L * (Z-Z1)/L + PHI_Z8 * (X-X1)/L * (Y-Y1)/L * (Z-Z1)/L
    
        N = SQRT(PHI_X*PHI_X + PHI_Y*PHI_Y + PHI_Z*PHI_Z)
    
        NORMAL(1) = PHI_X/N
        NORMAL(2) = PHI_Y/N
        NORMAL(3) = PHI_Z/N
    
    END SUBROUTINE LEVEL_COMPUTENORMAL_TYPE    
       
        
    SUBROUTINE LEVEL_COMPUTECURVATURE_TYPE(X,Y,Z,TYP,       CURV)

        REAL(8) :: X,Y,Z
        INTEGER :: TYP
        TYPE(OCTREE), POINTER :: CURRENT
        REAL(8) :: X1,X2,Y1,Y2,Z1,Z2,L
        REAL(8) :: CURV
        
        REAL(8) :: CURV1, CURV2, CURV3, CURV4, CURV5, CURV6, CURV7, CURV8
    
        CALL FIND_CELL_CONTAINING_POINT_TYPE(X,Y,Z,TYP,CURRENT)
   
        X1 = CURRENT%PHI_LEFTBOTTOMBACK%X
        X2 = CURRENT%PHI_RIGHTBOTTOMBACK%X
        Y1 = CURRENT%PHI_LEFTBOTTOMBACK%Y
        Y2 = CURRENT%PHI_LEFTTOPBACK%Y
        Z1 = CURRENT%PHI_LEFTBOTTOMBACK%Z
        Z2 = CURRENT%PHI_LEFTBOTTOMFRONT%Z
        L = CURRENT%LENGTH   
    
        CALL OCTREE_CURVATURE(CURRENT%PHI_LEFTBOTTOMBACK,CURV1) 
        CALL OCTREE_CURVATURE(CURRENT%PHI_LEFTTOPBACK,CURV2) 
        CALL OCTREE_CURVATURE(CURRENT%PHI_RIGHTBOTTOMBACK,CURV3) 
        CALL OCTREE_CURVATURE(CURRENT%PHI_RIGHTTOPBACK,CURV4) 
        CALL OCTREE_CURVATURE(CURRENT%PHI_LEFTBOTTOMFRONT,CURV5)     
        CALL OCTREE_CURVATURE(CURRENT%PHI_LEFTTOPFRONT,CURV6) 
        CALL OCTREE_CURVATURE(CURRENT%PHI_RIGHTBOTTOMFRONT,CURV7) 
        CALL OCTREE_CURVATURE(CURRENT%PHI_RIGHTTOPFRONT,CURV8) 
    
        CURV = CURV1 * (X2-X)/L * (Y2-Y)/L * (Z2-Z)/L + CURV2 * (X2-X)/L * (Y-Y1)/L * (Z2-Z)/L + CURV3* (X-X1)/L * (Y2-Y)/L * (Z2-Z)/L &
        + CURV4 * (X-X1)/L * (Y-Y1)/L * (Z2-Z)/L + CURV5 * (X2-X)/L * (Y2-Y)/L * (Z-Z1)/L + CURV6  * (X2-X)/L * (Y-Y1)/L * (Z-Z1)/L &
        + CURV7 * (X-X1)/L * (Y2-Y)/L * (Z-Z1)/L + CURV8 * (X-X1)/L * (Y-Y1)/L * (Z-Z1)/L
    
    END SUBROUTINE LEVEL_COMPUTECURVATURE_TYPE

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!   CHECK_REMESH
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE CHECK_REMESH(TYPE1,      B)
        INTEGER :: TYPE1
        LOGICAL :: B

        INTEGER :: POINTSNUM
        REAL(8), ALLOCATABLE :: POINTS(:,:) 
        INTEGER :: FACESNUM
        INTEGER, ALLOCATABLE:: FACES(:,:)

        REAL(8) :: DELTA
        INTEGER :: I

        REAL(8) :: X 
        REAL(8) :: Y 
        REAL(8) :: Z 

        REAL(8) :: TEMP1
        REAL(8) :: TEMP2


        REAL(8) :: VEC1(3)
        REAL(8) :: VEC2(3)
        REAL(8) :: VEC3(3)
        REAL(8) :: VEC4(3)

        REAL(8) :: R1
        REAL(8) :: R2
        REAL(8) :: R3
        REAL(8) :: R4
        REAL(8) :: R
        REAL(8) :: RR

        REAL(8) :: VEC_C(3)
        REAL(8) :: N0(3)
        REAL(8) :: N1(3)
        
        TYPE(OCTREE), POINTER :: CURRENT

        CALL FIND_CELL_CONTAINING_POINT_TYPE(X,Y,Z,TYPE1, CURRENT)

        IF (TYPE1==0) THEN 
            POINTSNUM = FLUID_SURFACE_POINTS_NUM
            POINTS = FLUID_SURFACE_POINTS
            FACESNUM = FLUID_SURFACE_FACES_NUM
            FACES = FLUID_SURFACE_FACES
        ELSE IF(TYPE1==1)THEN
            POINTSNUM = PROPEL_SURFACE_POINTS_NUM
            POINTS = PROPEL_SURFACE_POINTS
            FACESNUM = PROPEL_SURFACE_FACES_NUM
            FACES = PROPEL_SURFACE_FACES
        ELSE IF(TYPE1==2)THEN
            POINTSNUM = CASE_SURFACE_POINTS_NUM
            POINTS = CASE_SURFACE_POINTS
            FACESNUM = CASE_SURFACE_FACES_NUM
            FACES = CASE_SURFACE_FACES
        END IF

        DELTA= MAXSIZE * 2**(-MAX_TREE_DEPTH)/2.0

        B = .FALSE.

        DO I=1,POINTSNUM
            X = POINTS(1,I);
            Y = POINTS(2,I);
            Z = POINTS(3,I);
            CALL LEVEL_COMPUTEPHI_TYPE(X,Y,Z,TYPE1, TEMP1)
            IF (TEMP1 > 3.0*DELTA) THEN
                B = .TRUE.
                RETURN
            END IF
        END DO

        IF (.NOT.B) THEN
            DO I=1, FACESNUM

            VEC1(1) = POINTS(1, FACES(1,I))
            VEC1(2) = POINTS(2, FACES(1,I))
            VEC1(3) = POINTS(3, FACES(1,I))

            VEC2(1) = POINTS(1, FACES(2,I))
            VEC2(2) = POINTS(2, FACES(2,I))
            VEC2(3) = POINTS(3, FACES(2,I))

            VEC3(1) = POINTS(1, FACES(2,I))
            VEC3(2) = POINTS(2, FACES(2,I))
            VEC3(3) = POINTS(3, FACES(2,I))

            R1 = SQRT((VEC1(1)-VEC2(1))*(VEC1(1)-VEC2(1)) + (VEC1(2)-VEC2(2))*(VEC1(2)-VEC2(2)) + (VEC1(3)-VEC2(3))*(VEC1(3)-VEC2(3)))
            R2 = SQRT((VEC2(1)-VEC3(1))*(VEC2(1)-VEC3(1)) + (VEC2(2)-VEC3(2))*(VEC2(2)-VEC3(2)) + (VEC2(3)-VEC3(3))*(VEC2(3)-VEC3(3)))
            R3 = SQRT((VEC3(1)-VEC1(1))*(VEC3(1)-VEC1(1)) + (VEC3(2)-VEC1(2))*(VEC3(2)-VEC1(2)) + (VEC3(3)-VEC1(3))*(VEC3(3)-VEC1(3)))

            R = MAX(MAX(R1, R2), R3)
            RR = MIN(MIN(R1, R2), R3);

            IF(TYPE1==0 .AND. R > 3.0* FLUID_SIZE) THEN
                B = .TRUE.
                EXIT
            END IF
            IF(TYPE1==1 .AND. R > 3.0* STRUCT_SIZE) THEN
                B = .TRUE.
                EXIT
            END IF
            IF(RR < DELTA) THEN
                B = .TRUE.
                EXIT
            END IF

            VEC_C(1) = (VEC1(1)+VEC2(1)+VEC3(1))/3.0
            VEC_C(2) = (VEC1(2)+VEC2(2)+VEC3(2))/3.0
            VEC_C(3) = (VEC1(3)+VEC2(3)+VEC3(3))/3.0

            CALL LEVEL_COMPUTEPHI_TYPE(VEC_C(1), VEC_C(2), VEC_C(3),TYPE1, TEMP1)
            IF(ABS(TEMP1) > 3.0*DELTA) THEN
                B = .TRUE.
                EXIT
            END IF

            CALL LEVEL_COMPUTENORMAL_TYPE(VEC_C(1), VEC_C(2), VEC_C(3),TYPE1, N0)
            CALL VEC_CURL2(VEC2,VEC3,VEC2,VEC1,N1)
            CALL VEC_NORMALIZATION(N1)
            CALL VEC_INNER_PRODUCT1(N0,N1,TEMP2)
            IF(TEMP2<0.8) THEN
                B= .TRUE.
                EXIT
            END IF


            END DO
        END IF
    
    END SUBROUTINE

END MODULE