module petscvariables

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscis.h"

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

