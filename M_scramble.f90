module M_scramble
   implicit none
   private
   public scramble

contains

   function scramble( number_of_values ) result(array)

!@(#) M_random::scramble(3f): return integer array of random values 1 to N.
      integer,intent(in)    :: number_of_values
      integer,allocatable   :: array(:)
      integer               :: i, j, k, m, n
      integer               :: temp
      real                  :: u

      array=[(i,i=1,number_of_values)]

! The intrinsic RANDOM_NUMBER(3f) returns a real number (or an array
! of such) from the uniform distribution over the interval [0,1). (ie.
! it includes 0 but not 1.).
!
! To have a discrete uniform distribution on
! the integers {n, n+1, ..., m-1, m} carve the continuous distribution
! up into m+1-n equal sized chunks, mapping each chunk to an integer.
!
! One way is:
!   call random_number(u)
!   j = n + FLOOR((m+1-n)*u)  ! choose one from m-n+1 integers

      n=1
      m=number_of_values
      do k=1,2
         do i=1,m
            call random_number(u)
            j = n + FLOOR((m+1-n)*u)
            ! switch values
            temp=array(j)
            array(j)=array(i)
            array(i)=temp
         enddo
      enddo

   end function scramble

   end module M_scramble