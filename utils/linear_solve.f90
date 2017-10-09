module linear_solve

  implicit none

  public :: lu_decomposition
  public :: lu_back_substitution

  interface lu_decomposition
     module procedure lu_decomposition_real
     module procedure lu_decomposition_complex
  end interface

  interface lu_back_substitution
     module procedure lu_back_substitution_real
     module procedure lu_back_substitution_complex
  end interface

contains

  subroutine lu_decomposition_real (lu, idx, d)

    implicit none
    
    real, dimension (:,:), intent (in out) :: lu
    integer, dimension (:), intent (out) :: idx
    real, intent (out) :: d

    real, parameter :: zero = 1.0e-20
    real, dimension (size(lu,1)) :: vv
    real, dimension (size(lu,2)) :: dum

    integer :: j, n, imax

    n = size(lu,1)

    d = 1.0
    vv = maxval(abs(lu),dim=2)
    if (any(vv==0.0)) &
         write (*,*) 'singular matrix in lu_decomposition'
    vv = 1.0/vv
    do j = 1, n
       imax = (j-1) + imaxloc(vv(j:n)*abs(lu(j:n,j)))
       if (j /= imax) then
          dum = lu(imax,:)
          lu(imax,:) = lu(j,:)
          lu(j,:) = dum
          d = -d
          vv(imax) = vv(j)
       end if
       idx(j) = imax
       if (lu(j,j)==0.0) lu(j,j) = zero
       lu(j+1:n,j) = lu(j+1:n,j)/lu(j,j)
       lu(j+1:n,j+1:n) = lu(j+1:n,j+1:n) - spread(lu(j+1:n,j),2,n-j) &
            * spread(lu(j,j+1:n),1,n-j)
    end do

  end subroutine lu_decomposition_real

  subroutine lu_decomposition_complex (lu, idx, d)

    implicit none
    
    complex, dimension (:,:), intent (in out) :: lu
    integer, dimension (:), intent (out) :: idx
    real, intent (out) :: d

    real, parameter :: zero = 1.0e-20
    real, dimension (size(lu,1)) :: vv
    complex, dimension (size(lu,2)) :: dum

    integer :: j, n, imax

    n = size(lu,1)

    d = 1.0
    vv = maxval(cabs(lu),dim=2)
    if (any(vv==0.0)) &
         write (*,*) 'singular matrix in lu_decomposition'
    vv = 1.0/vv
    do j = 1, n
       imax = (j-1) + imaxloc(vv(j:n)*cabs(lu(j:n,j)))
       if (j /= imax) then
          dum = lu(imax,:)
          lu(imax,:) = lu(j,:)
          lu(j,:) = dum
          d = -d
          vv(imax) = vv(j)
       end if
       idx(j) = imax
       if (lu(j,j)==0.0) lu(j,j) = zero
       lu(j+1:n,j) = lu(j+1:n,j)/lu(j,j)
       lu(j+1:n,j+1:n) = lu(j+1:n,j+1:n) - spread(lu(j+1:n,j),2,n-j) &
            * spread(lu(j,j+1:n),1,n-j)
    end do

  end subroutine lu_decomposition_complex

  subroutine lu_back_substitution_real (lu, idx, b)

    implicit none

    real, dimension (:,:), intent (in) :: lu
    integer, dimension (:), intent (in) :: idx
    real, dimension (:), intent (in out) :: b

    integer :: i, n, ii, ll
    real :: summ

    n = size(lu,1)
    ii = 0
    do i = 1, n
       ll = idx(i)
       summ = b(ll)
       b(ll) = b(i)
       if (ii /= 0) then
          summ = summ - dot_product(lu(i,ii:i-1),b(ii:i-1))
       else if (summ /= 0.0) then
          ii = i
       end if
       b(i) = summ
    end do
    do i = n, 1, -1
       b(i) = (b(i) - dot_product(lu(i,i+1:n),b(i+1:n))) / lu(i,i)
    end do

  end subroutine lu_back_substitution_real

  subroutine lu_back_substitution_complex (lu, idx, b)

    implicit none

    complex, dimension (:,:), intent (in) :: lu
    integer, dimension (:), intent (in) :: idx
    complex, dimension (:), intent (in out) :: b

    integer :: i, n, ii, ll
    complex :: summ

    n = size(lu,1)
    ii = 0
    do i = 1, n
       ll = idx(i)
       summ = b(ll)
       b(ll) = b(i)
       if (ii /= 0) then
          summ = summ - dot_product(lu(i,ii:i-1),b(ii:i-1))
       else if (summ /= 0.0) then
          ii = i
       end if
       b(i) = summ
    end do
    do i = n, 1, -1
       b(i) = (b(i) - dot_product(lu(i,i+1:n),b(i+1:n))) / lu(i,i)
    end do

  end subroutine lu_back_substitution_complex

  function imaxloc (array)
    real, dimension (:), intent (in) :: array
    integer :: imaxloc
    integer, dimension (1) :: imax
    imax = maxloc(array)
    imaxloc = imax(1)
  end function imaxloc

end module linear_solve
