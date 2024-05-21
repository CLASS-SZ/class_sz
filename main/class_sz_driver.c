/** @file class.c
 * Julien Lesgourgues, 17.04.2011
 * Boris Bolliet for the class_sz part, 2015-
 */

#include "class.h"
#include "time.h"

int main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct primordial pm;       /* for primordial spectra */
  struct nonlinear nl;        /* for non-linear spectra */
  struct transfers tr;        /* for transfer functions */
  struct spectra sp;          /* for output spectra */
  struct lensing le;          /* for lensed spectra */
  struct class_sz_structure tsz;     /* BB: additional structure for class_sz*/
  struct szcount csz;         /* BB: additional structure for sz cluster counts*/
  struct output op;           /* for output files */
  ErrorMsg errmsg;            /* for error messages */

  //BB: initialize w. additional class_sz structures
  clock_t start, end;
  start = clock();
  if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&tsz,&csz,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
    return _FAILURE_;
  }
  end = clock();
  double duration = ((double)end - start)/CLOCKS_PER_SEC;
  // printf("Time taken to execute input in seconds : %.3e\n", duration);


  start = clock();
  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }
  end = clock();
  duration = ((double)end - start)/CLOCKS_PER_SEC;
  // printf("Time taken to execute background in seconds : %.3e\n", duration);

  start = clock();
  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }
  end = clock();
  duration = ((double)end - start)/CLOCKS_PER_SEC;
  // printf("Time taken to execute thermo in seconds : %.3e\n", duration);

  start = clock();
  if (perturb_init(&pr,&ba,&th,&pt) == _FAILURE_) {
    printf("\n\nError in perturb_init \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }
  end = clock();
  duration = ((double)end - start)/CLOCKS_PER_SEC;
  // printf("Time taken to execute perturb in seconds : %.3e\n", duration);

  start = clock();
  if (primordial_init(&pr,&pt,&pm) == _FAILURE_) {
    printf("\n\nError in primordial_init \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }
  end = clock();
  duration = ((double)end - start)/CLOCKS_PER_SEC;
  // printf("Time taken to execute primordial in seconds : %.3e\n", duration);

  start = clock();
  if (nonlinear_init(&pr,&ba,&th,&pt,&pm,&nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_init \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }
  end = clock();
  duration = ((double)end - start)/CLOCKS_PER_SEC;
  // printf("Time taken to execute nonlionear in seconds : %.3e\n", duration);

  start = clock();
  if (transfer_init(&pr,&ba,&th,&pt,&nl,&tr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }
  end = clock();
  duration = ((double)end - start)/CLOCKS_PER_SEC;
  // printf("Time taken to execute transfer in seconds : %.3e\n", duration);

  start = clock();
  if (spectra_init(&pr,&ba,&pt,&pm,&nl,&tr,&sp,&tsz,&th,&le) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }
  end = clock();
  duration = ((double)end - start)/CLOCKS_PER_SEC;
  // printf("Time taken to execute spectra in seconds : %.3e\n", duration);

  start = clock();
  if (lensing_init(&pr,&pt,&sp,&nl,&le) == _FAILURE_) {
    printf("\n\nError in lensing_init \n=>%s\n",le.error_message);
    return _FAILURE_;
  }
  end = clock();
  duration = ((double)end - start)/CLOCKS_PER_SEC;
  // printf("Time taken to execute lensing in seconds : %.3e\n", duration);


  start = clock();
  //BB: class_sz_cosmo module
  if (class_sz_cosmo_init(&ba,&th,&pt,&nl,&pm,&sp,&le,&tsz,&pr) == _FAILURE_) {
    printf("\nError in class_sz cosmo module\n");
    return _FAILURE_;
  }
  end = clock();
  duration = ((double)end - start)/CLOCKS_PER_SEC;

  start = clock();
  //BB: class_sz tabulate module
  if (class_sz_tabulate_init(&ba,&th,&pt,&nl,&pm,&sp,&le,&tsz,&pr) == _FAILURE_) {
    printf("\nError in class_sz tabulate module\n");
    return _FAILURE_;
  }
  end = clock();
  duration = ((double)end - start)/CLOCKS_PER_SEC;

  start = clock();
  //BB: class_sz integrate module
  if (class_sz_integrate_init(&ba,&th,&pt,&nl,&pm,&sp,&le,&tsz,&pr) == _FAILURE_) {
    printf("\nError in class_sz integrate module\n");
    return _FAILURE_;
  }
  end = clock();
  duration = ((double)end - start)/CLOCKS_PER_SEC;
  // printf("Time taken to execute class_sz in seconds : %.3e\n", duration);

  start = clock();
  // //BB: sz cluster counts module
  if (szcount_init(&ba,&nl,&pm,&tsz,&csz) == _FAILURE_) {
    printf("\nError in class_sz cluster counts module\n");
    return _FAILURE_;
  }
  end = clock();
  duration = ((double)end - start)/CLOCKS_PER_SEC;
  // printf("Time taken to execute szcountsz in seconds : %.3e\n", duration);

  start = clock();
  if (output_init(&ba,&th,&pt,&pm,&tr,&sp,&nl,&le,&op) == _FAILURE_) {
    printf("\n\nError in output_init \n=>%s\n",op.error_message);
    return _FAILURE_;
  }
  end = clock();
  duration = ((double)end - start)/CLOCKS_PER_SEC;
  // printf("Time taken to execute output in seconds : %.3e\n", duration);

  /****** all calculations done, now free the structures ******/

  // //BB: free sz cluster count module
  if (szcounts_free(&csz,&tsz) == _FAILURE_) {
    printf("\n\nError in szcounts_free\n");
    return _FAILURE_;
  }

  //BB: free class_sz module
  if (class_sz_free(&tsz) == _FAILURE_) {
    printf("\n\nError in class_sz_free\n");
    return _FAILURE_;
  }


  if (lensing_free(&le) == _FAILURE_) {
    printf("\n\nError in lensing_free \n=>%s\n",le.error_message);
    return _FAILURE_;
  }

  if (spectra_free(&sp) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  if (transfer_free(&tr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (nonlinear_free(&nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_free \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }

  if (primordial_free(&pm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  if (perturb_free(&pt) == _FAILURE_) {
    printf("\n\nError in perturb_free \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_free(&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_free \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (background_free(&ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;

}
