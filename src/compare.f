      function eqZERO( a)
      implicit none

      logical eqZERO
      double precision a

      eqZERO = (abs( a) .LE. 1.1*epsilon( 0.0))
      return
      end function

      function neZERO( a)
      implicit none

      logical neZERO
      double precision a

      neZERO = (abs( a) .GT. 1.1*epsilon( 0.0))
      return
      end function


      function eqREAL( a, b)
      implicit none

      logical eqREAL
      double precision a, b

      eqREAL = (abs( a - b) .LE. 1.1*epsilon( 0.0))

      return
      end function


      function neREAL( a, b)
      implicit none

      logical neREAL

c      double precision eps
c      parameter (eps = 1d-14)
      double precision a, b

c      neREAL = abs( a - b) .GT. eps
      neREAL = (abs( a - b) .GT. 1.1*epsilon( 0.0))

      return
      end function




