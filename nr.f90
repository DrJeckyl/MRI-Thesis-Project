MODULE nr
  interface 
     subroutine omp_set_dynamic( enable_expr )
       integer :: enable_expr
     end subroutine omp_set_dynamic
  end interface
  interface
     subroutine omp_set_num_threads( number_of_threads_expr )
       integer :: number_of_threads_expr
     end subroutine omp_set_num_threads
  end interface
  interface
     SUBROUTINE RungeKutta4(x,y,h,derivs)
       use TypesAndDefs
       real, intent(inout) :: x
       complex, dimension(:), intent(inout) :: y
       real, intent(in) :: h
       INTERFACE
          SUBROUTINE derivs(x,y,dydx)
            REAL, INTENT(IN) :: x
            complex, DIMENSION(:), INTENT(IN) :: y
            complex, DIMENSION(:), INTENT(OUT) :: dydx
          END SUBROUTINE derivs
       END INTERFACE
     END SUBROUTINE RungeKutta4
  END INTERFACE
  INTERFACE
  	SUBROUTINE EquationsOfMotion(x,y,dydx)
  		USE TypesAndDefs
  		REAL, INTENT(IN) :: x
  		complex, DIMENSION(:), INTENT(IN) :: y
  		complex, DIMENSION(:), INTENT(OUT) :: dydx
  	END SUBROUTINE EquationsOfMotion
  END INTERFACE
END MODULE nr
