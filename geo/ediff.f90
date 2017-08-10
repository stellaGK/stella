program ediff

  integer :: n, ndum, i, j
  real, dimension(4) :: r1, r2
  real :: rdiff, rmax

  open(unit=12,file='eik.out')
  open(unit=13,file='eik.out.save')

  
  do i=12,13
     read(i,*) 
     read(i,*) n, ndum, ndum
  enddo

  rmax=0.

  do k=1,2
     read(12,*) 
     read(13,*) 
     do i=1,2*n+1
        read(12,*) (r1(j),j=1,4)
        read(13,*) (r2(j),j=1,4)
        do j=1,4
           rmax=max(rmax,abs(r1(j)-r2(j))/max(abs(r1(j)),1.))
        enddo
     enddo
  enddo

  do k=1,2
     read(12,*) 
     read(13,*) 
     do i=1,2*n+1
        read(12,*) (r1(j),j=1,3)
        read(13,*) (r2(j),j=1,3)
        do j=1,3
           rmax=max(rmax,abs(r1(j)-r2(j))/max(abs(r1(j)),1.))
        enddo
     enddo
  enddo

  write(*,*) rmax

  stop
  end
