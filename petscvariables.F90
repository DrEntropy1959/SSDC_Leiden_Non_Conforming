module petscvariables


#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscis.h"

! petsc 3.10.0
! use petscvec

  Vec xpetsc
  Vec xlocpetsc

  Vec nxpetsc_Shell
  Vec nxlocpetsc_Shell

  Vec Jx_r_petsc_LGL
  Vec Jx_r_locpetsc_LGL

  Vec Jx_r_petsc_Gau_shell 
  Vec Jx_r_locpetsc_Gau_shell

  Vec upetsc
  Vec ulocpetsc

  Vec mutpetsc
  Vec mutlocpetsc

  Vec upetscWENO
  Vec ulocpetscWENO

  Vec upetscWENO_Shell
  Vec ulocpetscWENO_Shell

  Vec xpetscWENO_partner
  Vec xlocpetscWENO_partner

  Vec phipetsc
  Vec philocpetsc

  Vec uelempetsc
  Vec uelemlocpetsc

  Vec r_x_petsc
  Vec r_x_loc_petsc

  Vec xpetsc_shell
  Vec xlocpetsc_shell

end module petscvariables

