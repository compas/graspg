      module symzf_mod
      implicit none

      character*1500 :: stringorb
      logical :: flag_sym_csf
! For the zero space
      integer :: NBlockZ
      integer :: NUM_in_BLKZ(20)
      integer :: NUMZ
      character*256, dimension(:), allocatable ::  zconf1,  zconf1b
      character*256, dimension(:), allocatable ::  zconf2,  zconf2b
      character*256, dimension(:), allocatable ::  zconf3,  zconf3b
      integer, dimension(:), allocatable       :: lzconf1, lzconf1b
      integer, dimension(:), allocatable       :: lzconf2, lzconf2b
      integer, dimension(:), allocatable       :: lzconf3, lzconf3b

! For the Full space 
      integer :: NBlockF 
      integer :: NUM_in_BLKF(20)
      integer :: NUMF
      character*256, dimension(:), allocatable ::  fconf1,  fconf1b
      character*256, dimension(:), allocatable ::  fconf2,  fconf2b
      character*256, dimension(:), allocatable ::  fconf3,  fconf3b
      integer, dimension(:), allocatable       :: lfconf1, lfconf1b
      integer, dimension(:), allocatable       :: lfconf2, lfconf2b
      integer, dimension(:), allocatable       :: lfconf3, lfconf3b
      logical, dimension(:), allocatable       :: flagwrite
 
      end module symzf_mod

