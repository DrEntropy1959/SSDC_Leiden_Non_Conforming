module petscvariables

#include "include/petsc/finclude/petscsys.h"
#include "include/petsc/finclude/petscvec.h"
#include "include/petsc/finclude/petscis.h"

  use petsc

  Vec xpetsc

  Vec upetsc
  Vec ulocpetsc

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

