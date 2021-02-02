PROGRAM EXAMPLE_MHEAP

   USE MHEAP

   DOUBLE PRECISION :: POINT(3)
   TYPE(THEAP) :: HEAP
   INTEGER :: K

   ! Init a heap with at max 10 elements containing double precision arrays of lenght 3
   ! and the comparison function EUCLIDEAN_CMP, that compares the points' euclidean norm 
   CALL HEAP%INIT( 10, 3, GREATER1 )

   ! Insert some points
   CALL HEAP%INSERT( [ 1.0D0, 2.0D0, 1.0D0] )
   CALL HEAP%INSERT( [-1.0D0, 8.0D0,-3.0D0] )
   CALL HEAP%INSERT( [ 5.0D0, 1.0D0, 5.0D0] )
   CALL HEAP%INSERT( [ 2.0D0,-1.0D0,-2.0D0] )
   CALL HEAP%INSERT( [ 3.0D0, 7.0D0, 1.0D0] )
   CALL HEAP%INSERT( [-3.0D0,-8.0D0,-2.0D0] )
   CALL HEAP%INSERT( [ 6.0D0, 2.0D0, 5.0D0] )

   ! WRITE(*,*) ' Deleting nodes'

   ! CALL HEAP%DELETE( [ 6.0D0, 2.0D0, 5.0D0] )
   ! CALL HEAP%DELETE( [ 2.0D0,-1.0D0,-1.0D0] )
   ! CALL HEAP%DELETE( [ 2.0D0,-1.0D0,-2.0D0] )
   CALL HEAP%DELETE( K=4 )
   ! CALL HEAP%DELETE( [ 3.0D0, 7.0D0, 1.0D0] )
   CALL HEAP%DELETE( K=5 )
   CALL HEAP%DELETE( K=7 )

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
