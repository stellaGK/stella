# Variables
vars = %w[ NPSI NCHI R0EXP B0EXP PSI CHI Rgeom ageom q dqdpsi d2qdpsi2 p dpdpsi f fdfdpsi V rho_t shear dsheardpsi kappa delta_lower delta_upper dVdpsi dpsidrhotor GDPSI_av radius_av R_av TE DTEDPSI NE DNEDPSI TI DTIDPSI NI DNIDPSI ZEFF SIGNEO JBSBAV g11 g12 g22 g33 B dBdpsi dBdchi dPsidR dPsidZ dChidR dChidZ Jacobian R Z ].map{|s| s.downcase + '_chease'}


# Exclude npsi_chease, nchi_chease
zero_dim_vars = vars.slice(2...4)

# Exclude psi_chease, chi_chease
one_dim_vars = vars.slice(6...(vars.size-14))
two_dim_vars = vars.slice((vars.size-14)...vars.size)

#p [zero_dim_vars, one_dim_vars, two_dim_vars]

module_text = <<EOF

!> CHEASE Output Reader
!!
!! A module to read in datafiles from CHEASE.
!!
!! Written by Edmund Highcock
!! edmundhighcock@sourceforge.net
!!
!! This is free software released under the GPLv3
!!
!! Available quantities are:
!! 
!! Zero D:
!! #{zero_dim_vars.join(",")}
!! One D:
!! #{one_dim_vars.join(",").gsub(/(.{40,60},)/){"#$1\n!! "}}
!! Two D:
!! #{two_dim_vars.join(",").gsub(/(.{40,60},)/){"#$1\n!! "}}
!!
!! This module is generated automatically 
!! using 
!!    $ ruby generate_read_chease.rb
!!
!! DO NOT EDIT!
!! YOUR CHANGES WILL BE LOST

! An example chease_namelist will be on the wiki

module read_chease
  implicit none
  integer :: npsi_chease, nchi_chease 
  public :: npsi_chease, nchi_chease 
	
#{zero_dim_vars.map{|v| 
"
	real :: #{v}
	public :: #{v}
"}.join("\n")}
	
  real, dimension (:), allocatable :: psi_chease,chi_chease
	public :: psi_chease, chi_chease

#{one_dim_vars.map{ |v| 
"
	real, dimension(:), allocatable :: #{v}
	public :: #{v}"}.join("\n")
}
	
#{two_dim_vars.map{|v| 
"
	real, dimension(:,:), allocatable :: #{v}
	public :: #{v}
"}.join("\n")}

  integer :: infile=1212
  integer, parameter :: ncols = 5


  contains
  subroutine read_infile(filename)
    character (len=80) :: filename
    open(infile,file= filename)
    read(infile, *)
    read(infile, *) npsi_chease
    write(*, *) npsi_chease, "<---npsi_chease"
    read(infile, *)
    read(infile, *) nchi_chease
    write(*, *) nchi_chease, "<---nchi_chease"

#{zero_dim_vars.map{ |v|
"
    read(infile, *)
    read(infile, *) #{v}
"}.join("\n")}

    allocate(psi_chease(npsi_chease))
    read(infile, *)
    read(infile, *) psi_chease
    !write(*, *) psi_chease, "<---psi_chease"

    allocate(chi_chease(nchi_chease))
    read(infile, *)
    read(infile, *) chi_chease
    !write(*, *) chi_chease, "<---chi_chease"

#{one_dim_vars.map{ |v|
"
    allocate(#{v}(npsi_chease))
    read(infile, *)
    read(infile, *) #{v}
    !write(*, *) #{v}, \"<---#{v}\"
"}.join("\n")}


#{two_dim_vars.map{ |v|
"
    allocate(#{v}(npsi_chease,nchi_chease))
    read(infile, *)
    read(infile, *) #{v}
    !write(*, *) #{v}, \"<---#{v}\"
"}.join("\n")}
   

  end subroutine read_infile
  subroutine finish
    deallocate(psi_chease)
    deallocate(chi_chease)
		#{(one_dim_vars + two_dim_vars).map{|v|
"
    deallocate(#{v})"}.join("\n")}
  end subroutine finish



end module read_chease

!program test
!  use read_chease
!  call read_infile("ogyropsi.dat")
!
!end program test
	
EOF


File.open('read_chease.f90', 'w'){|f| f.puts module_text}

