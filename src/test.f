! this file included as part of workaround to ensure fortran files in src directory are compiled

      subroutine test(a)

      implicit none

      integer,intent(out)::a

      a = 10

      return

      end subroutine test
