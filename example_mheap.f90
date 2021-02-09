PROGRAM EXAMPLE_MHEAP

   USE MHEAP

   DOUBLE PRECISION :: POINT(3)
   TYPE(THEAP) :: HEAP
   INTEGER :: K, KARR(10)

   ! Init a heap with at max 10 elements containing double precision arrays of lenght 3
   ! and the comparison function EUCLIDEAN_CMP, that compares the points' euclidean norm 
   CALL HEAP%INIT( 10, 3, GREATER1 )

   ! Insert some points
   CALL HEAP%INSERT( [ 1.0D0, 2.0D0, 1.0D0], KARR(1) )
   CALL HEAP%INSERT( [-1.0D0, 8.0D0,-3.0D0], KARR(2) )
   CALL HEAP%INSERT( [ 5.0D0, 1.0D0, 5.0D0], KARR(3) )
   CALL HEAP%INSERT( [ 2.0D0,-1.0D0,-2.0D0], KARR(4) )
   CALL HEAP%INSERT( [ 3.0D0, 7.0D0, 1.0D0], KARR(5) )
   CALL HEAP%INSERT( [-3.0D0,-8.0D0,-2.0D0], KARR(6) )
   CALL HEAP%INSERT( [ 6.0D0, 2.0D0, 5.0D0], KARR(7) )

   ! WRITE(*,*) ' Deleting nodes'

   ! CALL HEAP%DELETE( [ 2.0D0,-1.0D0,-1.0D0] )
101 FORMAT( '   Deleting point ', i2, ': ', 2(f4.1, '(', &
         f4.1, '), '), f4.1, '(', f4.1, ')' )
   CALL HEAP%PEEK( KARR(4), POINT )
   WRITE(*,101) 4, POINT(1), 2d0, POINT(2), -1d0, POINT(3), -2d0
   CALL HEAP%DELETE( K=KARR(4) )
   ! CALL HEAP%DELETE( [ 2.0D0,-1.0D0,-2.0D0] )
   CALL HEAP%PEEK( HEAP%INDEXAT(KARR(5)), POINT )
   WRITE(*,101) 5, POINT(1), 3d0, POINT(2), 7d0, POINT(3), 1d0
   CALL HEAP%DELETE( K=KARR(5) )
   ! CALL HEAP%DELETE( [ 3.0D0, 7.0D0, 1.0D0] )

   CALL HEAP%INSERT( [ 4.0D0, -1.0D0, 3.0D0], KARR(8) )
   CALL HEAP%INSERT( [ 2.0D0, -7.0D0, 2.0D0], KARR(9) )
   CALL HEAP%INSERT( [-2.0D0,  3.0D0, 4.0D0], KARR(10) )
   
   CALL HEAP%PEEK( HEAP%INDEXAT(KARR(7)), POINT )
   WRITE(*,101) 7, POINT(1), 6d0, POINT(2), 2d0, POINT(3), 5d0
   CALL HEAP%DELETE( K=KARR(7) )
   ! CALL HEAP%DELETE( [ 6.0D0, 2.0D0, 5.0D0] )

   WRITE(*,*) 'K array = ', KARR
   ! Traversal in order of distance to the origin
   DO K = 1, HEAP%SIZE()
      CALL HEAP%POP( POINT )
      WRITE(*,*) POINT
   ENDDO

CONTAINS 

   LOGICAL FUNCTION EUCLIDEAN_CMP( POINT1, POINT2 )
      DOUBLE  PRECISION, INTENT(IN) :: POINT1(:), POINT2(:)
      EUCLIDEAN_CMP = SUM(POINT1**2.0D0) < SUM(POINT2**2.0D0)
   END FUNCTION EUCLIDEAN_CMP
   
   LOGICAL FUNCTION GREATER1( NODE1, NODE2 )
     DOUBLE  PRECISION, INTENT(IN) :: NODE1(:), NODE2(:)
     GREATER1 = NODE1(1) < NODE2(1)
   END FUNCTION GREATER1

END PROGRAM EXAMPLE_MHEAP
