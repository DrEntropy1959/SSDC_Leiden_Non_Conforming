#include <metis.h>
#include <stddef.h>
#include <stdio.h>

int calcMetisPartitions(idx_t nelems, idx_t nvertices, idx_t nparts, idx_t *xadj, 
  idx_t *adj, idx_t ncommon, idx_t *epart, idx_t *npart, int nnz)
{
  idx_t options[METIS_NOPTIONS];

  int ierr;

  // set default options
    ierr = METIS_SetDefaultOptions(options);
  // change to fortran ordering
  options[METIS_OPTION_NUMBERING] = 0;

  idx_t objval;

  objval = 0;
  int i;

  printf(" Metis is attempting to divide the domain into %d partitions\n",nparts);

//  printf("HERE %d",xadj[594]);

  ierr = METIS_PartMeshDual(&nelems, &nvertices, xadj, adj, NULL, NULL, &ncommon,
    &nparts, NULL, NULL, &objval, epart, npart);

//  printf("HERE %d",xadj[594]);

  printf("   Returned from metis, %d\n",ierr-METIS_OK);
  printf(" ===============================================================\n");

  return 0;
}
