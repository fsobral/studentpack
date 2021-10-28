module packmod

  ! COMMON PARAMETERS
  real(kind=8), parameter :: PEN = 1.0D+2
  real(kind=8), parameter :: ERR = 1.0D-2
  
  ! COMMON SCALARS
  integer      :: nbr,ndim,nite,nfix,nregh,nregw
  real(kind=8) :: cH,cW,dd,MINDIST

  ! COMMON ARRAYS
  integer,      pointer :: start(:,:),next(:),br(:,:)
  real(kind=8), pointer :: diagb(:,:,:),fcoord(:,:),frad(:)

contains

end module packmod
