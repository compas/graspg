      module symexpand_mod
      implicit none
! CYC-2023/11/26
      integer             :: nblocksym
! Number of the labelling orbitals
      integer             :: nonsym, nonsyminf
! String for all the labelling orbitals
      character(len=1070) :: strlaborb
! CYC-2023/11/26

! TotCSFs_perblock(1,i): NonSysmbolicCSFs, NCFGEN in setham_gg and
! findtype; 
! TotCSFs_perblock(2,i): the number of CSFs in the original small list
! (name.g file),
! i.e., NCFSMALL (obtained originally from the ncsfsymperblok array) in
! setham_gg
! TotCSFs_perblock(3,i): the number of CSFs in the fully expanded
! rcsf.out (name.c) file
! i.e., NCFTOT in setham_gg and findtype
      integer :: TotCSFs_perblock(3,50), NCFGTOT
      integer :: nmaxgen, ncsfDF1, scfDF1(50)
      integer :: nsym_orb(21,3)
! nsym_orb(i,1): number of this symobolic symetery; 
! nsym_orb(i,2): beginning position in orb string;
! nsym_orb(i,3): End       position in orb string;
      integer, dimension(:,:), allocatable :: map1
      integer, dimension(:,:,:), allocatable :: map2

      ! MATCH the CSFs in rcsf.inp to those in rcsfg.inp
      ! to include the CI within the diagnoal block CSFs generated by
      ! the same symmetry-ordered-CSF (CSFG).
      integer, dimension(:), allocatable :: MAPCS

! occupation number of the normal CSF 
      integer, dimension(:,:), allocatable :: iqa_csf

      end module symexpand_mod