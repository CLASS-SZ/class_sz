/** @file input.c Documented input module.
 *
 * Julien Lesgourgues, 27.08.2010
 */

#include "input.h"

/**
 * Use this routine to extract initial parameters from files 'xxx.ini'
 * and/or 'xxx.pre'. They can be the arguments of the main() routine.
 *
 * If class is embedded into another code, you will probably prefer to
 * call directly input_init() in order to pass input parameters
 * through a 'file_content' structure.
 */

int input_init_from_arguments(
                              int argc,
                              char **argv,
                              struct precision * ppr,
                              struct background *pba,
                              struct thermo *pth,
                              struct perturbs *ppt,
                              struct transfers *ptr,
                              struct primordial *ppm,
                              struct spectra *psp,
                              struct nonlinear * pnl,
                              struct lensing *ple,
                              struct class_sz_structure *pclass_sz, //BB: added for class_sz
                              struct szcount *pcsz, //BB: added for class_sz
                              struct output *pop,
                              ErrorMsg errmsg
                              ) {

  /** Summary: */

  /** - define local variables */

  struct file_content fc;             /** - --> the final structure with all parameters */
  struct file_content fc_input;       /** - --> a temporary structure with all input parameters */
  struct file_content fc_precision;   /** - --> a temporary structure with all precision parameters */
  struct file_content fc_root;        /** - --> a temporary structure with only the root name */
  struct file_content fc_inputroot;   /** - --> sum of fc_inoput and fc_root */
  struct file_content * pfc_input;    /** - --> a pointer to either fc_root or fc_inputroot */

  char input_file[_ARGUMENT_LENGTH_MAX_];
  char precision_file[_ARGUMENT_LENGTH_MAX_];
  char tmp_file[_ARGUMENT_LENGTH_MAX_+26]; // 26 is enough to extend the file name [...] with the characters "output/[...]%02d_parameters.ini" (as done below)

  int i;
  char extension[5];
  FileArg stringoutput, inifilename;
  int flag1, filenum;

  pfc_input = &fc_input;

  /** - Initialize the two file_content structures (for input
      parameters and precision parameters) to some null content. If no
      arguments are passed, they will remain null and inform
      init_params() that all parameters take default values. */

  fc.size = 0;
  fc_input.size = 0;
  fc_precision.size = 0;
  input_file[0]='\0';
  precision_file[0]='\0';

  /** - If some arguments are passed, identify eventually some 'xxx.ini'
      and 'xxx.pre' files, and store their name. */

  if (argc > 1) {
    for (i=1; i<argc; i++) {
      strncpy(extension,(argv[i]+strlen(argv[i])-4),4);
      extension[4]='\0';
      if (strcmp(extension,".ini") == 0) {
        class_test(input_file[0] != '\0',
                   errmsg,
                   "You have passed more than one input file with extension '.ini', choose one.");
        strcpy(input_file,argv[i]);
      }
      else if (strcmp(extension,".pre") == 0) {
        class_test(precision_file[0] != '\0',
                   errmsg,
                   "You have passed more than one precision with extension '.pre', choose one.");
        strcpy(precision_file,argv[i]);
      }
      else {
        fprintf(stdout,"Warning: the file '%s' has an extension different from .ini and .pre, so it has been ignored\n",argv[i]);
      }
    }
  }

  /** - if there is an 'xxx.ini' file, read it and store its content. */

  if (input_file[0] != '\0'){

    class_call(parser_read_file(input_file,&fc_input,errmsg),
               errmsg,
               errmsg);

    /** - check whether a root name has been set */

    class_call(parser_read_string(&fc_input,"root",&stringoutput,&flag1,errmsg),
               errmsg, errmsg);

    /** - if root has not been set, use root=output/inputfilennameN_ */

    if (flag1 == _FALSE_){
      //printf("strlen-4 = %zu\n",strlen(input_file)-4);
      strncpy(inifilename, input_file, strlen(input_file)-4);
      inifilename[strlen(input_file)-4] = '\0';
      for (filenum = 0; filenum < 100; filenum++){
        sprintf(tmp_file,"output/%s%02d_cl.dat", inifilename, filenum);
        if (file_exists(tmp_file) == _TRUE_)
          continue;
        sprintf(tmp_file,"output/%s%02d_pk.dat", inifilename, filenum);
        if (file_exists(tmp_file) == _TRUE_)
          continue;
        sprintf(tmp_file,"output/%s%02d_tk.dat", inifilename, filenum);
        if (file_exists(tmp_file) == _TRUE_)
          continue;
        sprintf(tmp_file,"output/%s%02d_parameters.ini", inifilename, filenum);
        if (file_exists(tmp_file) == _TRUE_)
          continue;
        break;
      }
      class_call(parser_init(&fc_root,
                             1,
                             fc_input.filename,
                             errmsg),
                 errmsg,errmsg);
      sprintf(fc_root.name[0],"root");
      sprintf(fc_root.value[0],"output/%s%02d_",inifilename,filenum);
      fc_root.read[0] = _FALSE_;
      class_call(parser_cat(&fc_input,&fc_root,&fc_inputroot,errmsg),
                 errmsg,
                 errmsg);
      class_call(parser_free(&fc_input),errmsg,errmsg);
      class_call(parser_free(&fc_root),errmsg,errmsg);
      pfc_input = &fc_inputroot;
    }
  }

  /** - if there is an 'xxx.pre' file, read it and store its content. */

  if (precision_file[0] != '\0')

    class_call(parser_read_file(precision_file,&fc_precision,errmsg),
               errmsg,
               errmsg);

  /** - if one or two files were read, merge their contents in a
      single 'file_content' structure. */

  if ((input_file[0]!='\0') || (precision_file[0]!='\0'))

    class_call(parser_cat(pfc_input,&fc_precision,&fc,errmsg),
               errmsg,
               errmsg);

  class_call(parser_free(pfc_input),errmsg,errmsg);
  class_call(parser_free(&fc_precision),errmsg,errmsg);

  /** - Finally, initialize all parameters given the input 'file_content'
      structure.  If its size is null, all parameters take their
      default values. */

  class_call(input_init(&fc,
                        ppr,
                        pba,
                        pth,
                        ppt,
                        ptr,
                        ppm,
                        psp,
                        pnl,
                        ple,
                        pclass_sz, //BB: added for class_sz
                        pcsz, //BB: added for class_sz
                        pop,
                        errmsg),
             errmsg,
             errmsg);

  class_call(parser_free(&fc),errmsg,errmsg);

  return _SUCCESS_;
}

/**
 * Initialize each parameter, first to its default values, and then
 * from what can be interpreted from the values passed in the input
 * 'file_content' structure. If its size is null, all parameters keep
 * their default values.
 *
 */

int input_init(
               struct file_content * pfc,
               struct precision * ppr,
               struct background *pba,
               struct thermo *pth,
               struct perturbs *ppt,
               struct transfers *ptr,
               struct primordial *ppm,
               struct spectra *psp,
               struct nonlinear * pnl,
               struct lensing *ple,
               struct class_sz_structure *pclass_sz, //BB: added for class_sz
               struct szcount *pcsz, //BB: added for class_sz
               struct output *pop,
               ErrorMsg errmsg
               ) {

  int flag1;
  double param1;
  int counter, index_target, i;
  double * unknown_parameter;
  int unknown_parameters_size;
  int fevals=0;
  double xzero;
  int target_indices[_NUM_TARGETS_];
  double *dxdF, *x_inout;

  char string1[_ARGUMENT_LENGTH_MAX_];
  FILE * param_output;
  FILE * param_unused;
  char param_output_name[_LINE_LENGTH_MAX_];
  char param_unused_name[_LINE_LENGTH_MAX_];

  struct fzerofun_workspace fzw;

  /**
   * Before getting into the assignment of parameters,
   * and before the shooting, we want to already fix our precision parameters.
   *
   * No precision parameter should depend on any input parameter
   *
   */

  class_call(input_read_precisions(pfc,
                                   ppr,
                                   pba,
                                   pth,
                                   ppt,
                                   ptr,
                                   ppm,
                                   psp,
                                   pnl,
                                   ple,
                                   pclass_sz,
                                   pop,
                                   errmsg),
             errmsg,
             errmsg);



  /**
   * In CLASS, we can do something we call 'shooting', where a variable,
   *  which is not directly given is calculated by another variable
   *  through successive runs of class.
   *
   * This is needed for variables which do not immediately follow from
   *  other input parameters. An example is theta_s, the angular scale
   *  of the sound horizon giving us the horizontal peak positions.
   *  This quantity can only replace the hubble parameter h, if we
   *  run all the way into class through to thermodynamics to figure out
   *  how h and theta_s relate numerically.
   *
   * A default parameter for h is chosen, and then we shoot through
   *  CLASS, finding what the corresponding theta_s is. We adjust our
   *  initial h, and shoot again, repeating this process until a
   *  suitable value for h is found which gives the correct
   *  100*theta_s value
   *
   * These two arrays must contain the strings of names to be searched
   *  for and the corresponding new parameter
   * The third array contains the module inside of which the old
   *  parameter is calculated
   *
   * See input_try_unknown_parameters for the actual shooting
   *
   */

  char * const target_namestrings[] = {"100*theta_s","Omega_dcdmdr","omega_dcdmdr",
                                       "Omega_scf","Omega_ini_dcdm","omega_ini_dcdm","sigma8","age"};
  char * const unknown_namestrings[] = {"h","Omega_ini_dcdm","Omega_ini_dcdm",
                                        "scf_shooting_parameter","Omega_dcdmdr","omega_dcdmdr","A_s","H0"};
  enum computation_stage target_cs[] = {cs_thermodynamics, cs_background, cs_background,
                                        cs_background, cs_background, cs_background, cs_nonlinear,cs_background};

  int input_verbose = 0, int1, aux_flag, shooting_failed=_FALSE_;

  class_read_int("input_verbose",input_verbose);
  if (input_verbose >0) printf("Reading input parameters\n");

  /** - Do we need to fix unknown parameters? */
  unknown_parameters_size = 0;
  fzw.required_computation_stage = 0;
  for (index_target = 0; index_target < _NUM_TARGETS_; index_target++){
    // printf(">>>%d %s\n",index_target,target_namestrings[index_target]);
    class_call(parser_read_double(pfc,
                                  target_namestrings[index_target],
                                  &param1,
                                  &flag1,
                                  errmsg),
               errmsg,
               errmsg);
    if (flag1 == _TRUE_){
      /** - --> input_auxillary_target_conditions() takes care of the case where for
          instance Omega_dcdmdr is set to 0.0.
      */
      class_call(input_auxillary_target_conditions(pfc,
                                                   index_target,
                                                   param1,
                                                   &aux_flag,
                                                   errmsg),
                 errmsg, errmsg);
      if (aux_flag == _TRUE_){
        //printf("Found target: %s\n",target_namestrings[index_target]);
        target_indices[unknown_parameters_size] = index_target;
        fzw.required_computation_stage = MAX(fzw.required_computation_stage,target_cs[index_target]);
        unknown_parameters_size++;
      }
    }
  }

  // printf(">>> unknown_param_size = %d\n",unknown_parameters_size);

  /**
   * Case with unknown parameters...
   *
   * Here we start shooting (see above for explanation of shooting)
   *
   *  */
  if (unknown_parameters_size > 0) {

    /* Create file content structure with additional entries */
    class_call(parser_init(&(fzw.fc),
                           pfc->size+unknown_parameters_size,
                           pfc->filename,
                           errmsg),
               errmsg,errmsg);
    /* Copy input file content to the new file content structure: */
    memcpy(fzw.fc.name, pfc->name, pfc->size*sizeof(FileArg));
    memcpy(fzw.fc.value, pfc->value, pfc->size*sizeof(FileArg));
    memcpy(fzw.fc.read, pfc->read, pfc->size*sizeof(short));

    class_alloc(unknown_parameter,
                unknown_parameters_size*sizeof(double),
                errmsg);
    class_alloc(fzw.unknown_parameters_index,
                unknown_parameters_size*sizeof(int),
                errmsg);
    fzw.target_size = unknown_parameters_size;
    class_alloc(fzw.target_name,
                fzw.target_size*sizeof(enum target_names),
                errmsg);
    class_alloc(fzw.target_value,
                fzw.target_size*sizeof(double),
                errmsg);

    /** - --> go through all cases with unknown parameters: */
    for (counter = 0; counter < unknown_parameters_size; counter++){
      index_target = target_indices[counter];
      class_call(parser_read_double(pfc,
                                    target_namestrings[index_target],
                                    &param1,
                                    &flag1,
                                    errmsg),
                 errmsg,
                 errmsg);

      // store name of target parameter
      fzw.target_name[counter] = index_target;
      // store target value of target parameter
      fzw.target_value[counter] = param1;
      fzw.unknown_parameters_index[counter]=pfc->size+counter;
      // substitute the name of the target parameter with the name of the corresponding unknown parameter
      strcpy(fzw.fc.name[fzw.unknown_parameters_index[counter]],unknown_namestrings[index_target]);
      // printf("%d, %d: %s\n",counter,index_target,target_namestrings[index_target]);
    }
    // exit(0);

    if (unknown_parameters_size == 1){
      if (input_verbose > 0) {
        fprintf(
                stdout,
                "Computing unknown input parameter '%s' using input parameter '%s'\n",
                fzw.fc.name[fzw.unknown_parameters_index[0]],
                target_namestrings[fzw.target_name[0]]
                );
      }
      /* We can do 1 dimensional root finding */
      /* If shooting fails, postpone error to background module to play nice with MontePython. */
      class_call_try(input_find_root(&xzero,
                                     &fevals,
                                     &fzw,
                                     ppr->tol_shooting_1d,
                                     errmsg),
                     errmsg,
                     pba->shooting_error,
                     shooting_failed=_TRUE_);

      /* Store xzero */
      sprintf(fzw.fc.value[fzw.unknown_parameters_index[0]],"%e",xzero);
      if (input_verbose > 0) {
        fprintf(stdout," -> found '%s = %s'\n",
                fzw.fc.name[fzw.unknown_parameters_index[0]],
                fzw.fc.value[fzw.unknown_parameters_index[0]]);
      }
    }
    else{
      /* We need to do multidimensional root finding */

      if (input_verbose > 0) {
        fprintf(stdout,"Computing unknown input parameters\n");
      }
      class_alloc(x_inout,
                  sizeof(double)*unknown_parameters_size,
                  errmsg);
      class_alloc(dxdF,
                  sizeof(double)*unknown_parameters_size,
                  errmsg);
      class_call(input_get_guess(x_inout,
                                 dxdF,
                                 &fzw,
                                 errmsg),
                 errmsg, errmsg);

      class_call_try(fzero_Newton(input_try_unknown_parameters,
                                  x_inout,
                                  dxdF,
                                  unknown_parameters_size,
                                  1e-4,
                                  1e-6,
                                  &fzw,
                                  &fevals,
                                  errmsg),
                     errmsg, pba->shooting_error,shooting_failed=_TRUE_);



      /* Store xzero */
      for (counter = 0; counter < unknown_parameters_size; counter++){
        sprintf(fzw.fc.value[fzw.unknown_parameters_index[counter]],
                "%e",x_inout[counter]);
        if (input_verbose > 0) {
          fprintf(stdout," -> found '%s = %s'\n",
                  fzw.fc.name[fzw.unknown_parameters_index[counter]],
                  fzw.fc.value[fzw.unknown_parameters_index[counter]]);
        }
      }

      free(x_inout);
      free(dxdF);
    }

    if (input_verbose > 1) {
      fprintf(stdout,"Shooting completed using %d function evaluations\n",fevals);
    }


    /** - --> Read all parameters from tuned pfc */
    class_call(input_read_parameters(&(fzw.fc),
                                     ppr,
                                     pba,
                                     pth,
                                     ppt,
                                     ptr,
                                     ppm,
                                     psp,
                                     pnl,
                                     ple,
                                     pclass_sz, //BB: added for class_sz
                                     pcsz, //BB: added for class_sz
                                     pop,
                                     errmsg),
               errmsg,
               errmsg);

    /** - --> Set status of shooting */
    pba->shooting_failed = shooting_failed;

    /* all parameters read in fzw must be considered as read in
       pfc. At the same time the parameters read before in pfc (like
       theta_s,...) must still be considered as read (hence we could
       not do a memcopy) */
    for (i=0; i < pfc->size; i ++) {
      if (fzw.fc.read[i] == _TRUE_)
        pfc->read[i] = _TRUE_;
    }

    // Free tuned pfc
    parser_free(&(fzw.fc));
    /** - --> Free arrays allocated*/
    free(unknown_parameter);
    free(fzw.unknown_parameters_index);
    free(fzw.target_name);
    free(fzw.target_value);
  }
  /** - case with no unknown parameters */
  else{

    /** - --> just read all parameters from input pfc: */
    class_call(input_read_parameters(pfc,
                                     ppr,
                                     pba,
                                     pth,
                                     ppt,
                                     ptr,
                                     ppm,
                                     psp,
                                     pnl,
                                     ple,
                                     pclass_sz, //BB: added for class_sz
                                     pcsz, //BB: added for class_sz
                                     pop,
                                     errmsg),
               errmsg,
               errmsg);
  }

  /** - eventually write all the read parameters in a file, unread parameters in another file, and warnings about unread parameters */

  class_call(parser_read_string(pfc,"write parameters",&string1,&flag1,errmsg),
             errmsg,
             errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {

    sprintf(param_output_name,"%s%s",pop->root,"parameters.ini");
    sprintf(param_unused_name,"%s%s",pop->root,"unused_parameters");

    class_open(param_output,param_output_name,"w",errmsg);
    class_open(param_unused,param_unused_name,"w",errmsg);

    fprintf(param_output,"# List of input/precision parameters actually read\n");
    fprintf(param_output,"# (all other parameters set to default values)\n");
    fprintf(param_output,"# Obtained with CLASS %s (for developers: svn version %s)\n",_VERSION_,_SVN_VERSION_);
    fprintf(param_output,"#\n");
    fprintf(param_output,"# This file can be used as the input file of another run\n");
    fprintf(param_output,"#\n");

    fprintf(param_unused,"# List of input/precision parameters passed\n");
    fprintf(param_unused,"# but not used (just for info)\n");
    fprintf(param_unused,"#\n");

    for (i=0; i<pfc->size; i++) {
      if (pfc->read[i] == _TRUE_)
        fprintf(param_output,"%s = %s\n",pfc->name[i],pfc->value[i]);
      else
        fprintf(param_unused,"%s = %s\n",pfc->name[i],pfc->value[i]);
    }
    fprintf(param_output,"#\n");

    fclose(param_output);
    fclose(param_unused);
  }

  if(pclass_sz->create_ref_trispectrum_for_cobaya){

    sprintf(param_output_name,"%s%s%s%s",
            pclass_sz->path_to_ref_trispectrum_for_cobaya,
            "/tSZ_params_ref_",
            pclass_sz->append_name_trispectrum_ref,
            ".txt");

    class_open(param_output,param_output_name,"w",errmsg);

    fprintf(param_output,"# List of input/precision parameters actually read\n");
    fprintf(param_output,"# (all other parameters set to default values)\n");
    fprintf(param_output,"# Obtained with CLASS %s (for developers: svn version %s)\n",_VERSION_,_SVN_VERSION_);
    fprintf(param_output,"#\n");
    fprintf(param_output,"# This file can be used as the input file of another run\n");
    fprintf(param_output,"#\n");


    for (i=0; i<pfc->size; i++)
    fprintf(param_output,"%s = %s\n",pfc->name[i],pfc->value[i]);


    fprintf(param_output,"#\n");

    fclose(param_output);

  }

  class_call(parser_read_string(pfc,"write warnings",&string1,&flag1,errmsg),
             errmsg,
             errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {

    for (i=0; i<pfc->size; i++) {
      if (pfc->read[i] == _FALSE_)
        fprintf(stdout,"[WARNING: input line not recognized and not taken into account: '%s=%s']\n",pfc->name[i],pfc->value[i]);
    }
  }

  if (pnl->has_pk_eq == _TRUE_) {

    if (input_verbose > 0) {
      printf(" -> since you want to use Halofit with a non-zero wa_fld and the Pk_equal method,\n");
      printf("    calling background module to extract the effective w(tau), Omega_m(tau) parameters");
      printf("    required by this method\n");
    }
    class_call(input_prepare_pk_eq(ppr,pba,pth,pnl,input_verbose,errmsg),
               errmsg,
               errmsg);
  }

  return _SUCCESS_;

}
int input_read_precisions(
                          struct file_content * pfc,
                          struct precision * ppr,
                          struct background *pba,
                          struct thermo *pth,
                          struct perturbs *ppt,
                          struct transfers *ptr,
                          struct primordial *ppm,
                          struct spectra *psp,
                          struct nonlinear * pnl,
                          struct lensing *ple,
                          struct class_sz_structure *pclass_sz,
                          struct output *pop,
                          ErrorMsg errmsg
                          ) {
  /** - set all precision parameters to default values */

  /**
   * Declare initial params to read into
   * */
  class_call(input_default_precision(ppr),
             errmsg,
             errmsg);


  int int1;
  int flag1;
  double param1;
  char string1[_ARGUMENT_LENGTH_MAX_];

  /**
   * Parse all precision parameters
   * */

#define __PARSE_PRECISION_PARAMETER__
#include "precisions.h"
#undef __PARSE_PRECISION_PARAMETER__

#define __ASSIGN_DEFAULT_TSZ__
#include "class_sz_precisions.h"
#undef __ASSIGN_DEFAULT_TSZ__


#define __PARSE_TSZ_PARAMETER__
#include "class_sz_precisions.h"
#undef _PARSE_TSZ_PARAMETER__


  return _SUCCESS_;
}
int input_read_parameters(
                          struct file_content * pfc,
                          struct precision * ppr,
                          struct background *pba,
                          struct thermo *pth,
                          struct perturbs *ppt,
                          struct transfers *ptr,
                          struct primordial *ppm,
                          struct spectra *psp,
                          struct nonlinear * pnl,
                          struct lensing *ple,
                          struct class_sz_structure *pclass_sz,
                          struct szcount *pcsz,
                          struct output *pop,
                          ErrorMsg errmsg
                          ) {

  /** Summary: */

  /** - define local variables */

  int flag1,flag2,flag3;
  double param1,param2,param3;
  int N_ncdm=0,n,entries_read;
  int int1,fileentries;
  double scf_lambda;
  double fnu_factor;
  double * pointer1;
  char string1[_ARGUMENT_LENGTH_MAX_];
  char string2[_ARGUMENT_LENGTH_MAX_];
  double k1=0.;
  double k2=0.;
  double prr1=0.;
  double prr2=0.;
  double pii1=0.;
  double pii2=0.;
  double pri1=0.;
  double pri2=0.;
  double n_iso=0.;
  double f_iso=0.;
  double n_cor=0.;
  double c_cor=0.;
  double stat_f_idr = 7./8.;

  double Omega_tot;

  int i;

  double sigma_B; /* Stefan-Boltzmann constant in \f$ W/m^2/K^4 = Kg/K^4/s^3 \f$*/

  double rho_ncdm;
  double R0,R1,R2,R3,R4;
  double PSR0,PSR1,PSR2,PSR3,PSR4;
  double HSR0,HSR1,HSR2,HSR3,HSR4;

  double z_max=0.;
  int bin;
  int input_verbose=0;

  sigma_B = 2. * pow(_PI_,5) * pow(_k_B_,4) / 15. / pow(_h_P_,3) / pow(_c_,2);

  /** - set all input parameters to default values */

  class_call(input_default_params(pba,
                                  pth,
                                  ppt,
                                  ptr,
                                  ppm,
                                  psp,
                                  pnl,
                                  ple,
                                  pclass_sz, //BB: added for class_sz
                                  pcsz, //BB: added for class_sz
                                  pop),
             errmsg,
             errmsg);

  /** - if entries passed in file_content structure, carefully read
      and interpret each of them, and tune the relevant input
      parameters accordingly*/

  class_read_int("input_verbose",input_verbose);

  /** Knowing the gauge from the very beginning is useful (even if
      this could be a run not requiring perturbations at all: even in
      that case, knowing the gauge is important e.g. for fixing the
      sampling in momentum space for non-cold dark matter) */

  class_call(parser_read_string(pfc,"gauge",&string1,&flag1,errmsg),
             errmsg,
             errmsg);

  if (flag1 == _TRUE_) {

    if ((strstr(string1,"newtonian") != NULL) || (strstr(string1,"Newtonian") != NULL) || (strstr(string1,"new") != NULL)) {
      ppt->gauge = newtonian;
    }

    if ((strstr(string1,"synchronous") != NULL) || (strstr(string1,"sync") != NULL) || (strstr(string1,"Synchronous") != NULL)) {
      ppt->gauge = synchronous;
    }
  }


  /** (a) background parameters */

  /** - scale factor today (arbitrary) */
  class_read_double("a_today",pba->a_today);

  /** - h (dimensionless) and [\f$ H_0/c\f$] in \f$ Mpc^{-1} = h / 2997.9... = h * 10^5 / c \f$ */
  class_call(parser_read_double(pfc,"H0",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"h",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),
             errmsg,
             "In input file, you cannot enter both h and H0, choose one");
  if (flag1 == _TRUE_) {
    pba->H0 = param1 * 1.e3 / _c_;
    pba->h = param1 / 100.;
  }
  if (flag2 == _TRUE_) {
    pba->H0 = param2 *  1.e5 / _c_;
    pba->h = param2;
  }


  // printf(">>>> we have: hubble %.3e %.3e\n",pba->H0,pba->h);

  /** - Omega_0_g (photons) and T_cmb */
  class_call(parser_read_double(pfc,"T_cmb",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"Omega_g",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"omega_g",&param3,&flag3,errmsg),
             errmsg,
             errmsg);
  class_test(class_at_least_two_of_three(flag1,flag2,flag3),
             errmsg,
             "In input file, you can only enter one of T_cmb, Omega_g or omega_g, choose one");

  if (class_none_of_three(flag1,flag2,flag3)) {
    pba->Omega0_g = (4.*sigma_B/_c_*pow(pba->T_cmb,4.)) / (3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_);
  }
  else {

    if (flag1 == _TRUE_) {
      /** - Omega0_g = rho_g / rho_c0, each of them expressed in \f$ Kg/m/s^2 \f$*/
      /** - rho_g = (4 sigma_B / c) \f$ T^4 \f$*/
      /** - rho_c0 \f$ = 3 c^2 H_0^2 / (8 \pi G) \f$*/
      pba->Omega0_g = (4.*sigma_B/_c_*pow(param1,4.)) / (3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_);
      pba->T_cmb=param1;
    }

    if (flag2 == _TRUE_) {
      pba->Omega0_g = param2;
      pba->T_cmb=pow(pba->Omega0_g * (3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_) / (4.*sigma_B/_c_),0.25);
    }

    if (flag3 == _TRUE_) {
      pba->Omega0_g = param3/pba->h/pba->h;
      pba->T_cmb = pow(pba->Omega0_g * (3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_) / (4.*sigma_B/_c_),0.25);
    }
  }

  Omega_tot = pba->Omega0_g;

  /** - Omega_0_b (baryons) */
  class_call(parser_read_double(pfc,"Omega_b",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"omega_b",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  class_test(((flag1 == _TRUE_) && (flag2 == _TRUE_)),
             errmsg,
             "In input file, you can only enter one of Omega_b or omega_b, choose one");
  if (flag1 == _TRUE_)
    pba->Omega0_b = param1;
  if (flag2 == _TRUE_)
    pba->Omega0_b = param2/pba->h/pba->h;

  Omega_tot += pba->Omega0_b;

  /** - Omega_0_ur (ultra-relativistic species / massless neutrino) */

  /* (a) try to read N_ur */
  class_call(parser_read_double(pfc,"N_ur",&param1,&flag1,errmsg),
             errmsg,
             errmsg);

  /* these lines have been added for compatibility with deprecated syntax 'N_eff' instead of 'N_ur', in the future they could be suppressed */
  class_call(parser_read_double(pfc,"N_eff",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),
             errmsg,
             "In input file, you can only enter one of N_eff (deprecated syntax) or N_ur (up-to-date syntax), since they botgh describe the same, i.e. the contribution ukltra-relativistic species to the effective neutrino number");
  if (flag2 == _TRUE_) {
    param1 = param2;
    flag1 = _TRUE_;
    flag2 = _FALSE_;
  }
  /* end of lines for deprecated syntax */

  /* (b) try to read Omega_ur */
  class_call(parser_read_double(pfc,"Omega_ur",&param2,&flag2,errmsg),
             errmsg,
             errmsg);

  /* (c) try to read omega_ur */
  class_call(parser_read_double(pfc,"omega_ur",&param3,&flag3,errmsg),
             errmsg,
             errmsg);

  /* (d) infer the unpassed ones from the passed one */
  class_test(class_at_least_two_of_three(flag1,flag2,flag3),
             errmsg,
             "In input file, you can only enter one of N_eff, Omega_ur or omega_ur, choose one");

  if (class_none_of_three(flag1,flag2,flag3)) {
    pba->Omega0_ur = 3.046*7./8.*pow(4./11.,4./3.)*pba->Omega0_g;
  }
  else {

    if (flag1 == _TRUE_) {
      pba->Omega0_ur = param1*7./8.*pow(4./11.,4./3.)*pba->Omega0_g;
    }
    if (flag2 == _TRUE_) {
      pba->Omega0_ur = param2;
    }
    if (flag3 == _TRUE_) {
      pba->Omega0_ur = param3/pba->h/pba->h;
    }
  }

  class_call(parser_read_double(pfc,"ceff2_ur",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  if (flag1 == _TRUE_) ppt->three_ceff2_ur = 3.*param1;

  class_call(parser_read_double(pfc,"cvis2_ur",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  if (flag1 == _TRUE_) ppt->three_cvis2_ur = 3.*param1;

  Omega_tot += pba->Omega0_ur;

  /** - Omega_0_idr (interacting dark radiation) */
  /* Can take both the ethos parameters, and the NADM parameters */

  class_read_double("stat_f_idr",stat_f_idr);

  class_call(parser_read_double(pfc,"N_idr",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"N_dg",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"xi_idr",&param3,&flag3,errmsg),
             errmsg,
             errmsg);
  class_test(class_at_least_two_of_three(flag1,flag2,flag3),
             errmsg,
             "In input file, you can only enter one of N_idr, N_dg or xi_idr, choose one");

  if (flag1 == _TRUE_) {
    pba->T_idr = pow(param1/stat_f_idr*(7./8.)/pow(11./4.,(4./3.)),(1./4.)) * pba->T_cmb;
    if (input_verbose > 1)
      printf("You passed N_idr = N_dg = %e, this is equivalent to xi_idr = %e in the ETHOS notation. \n", param2, pba->T_idr/pba->T_cmb);
  }
  else if (flag2 == _TRUE_) {
    pba->T_idr = pow(param2/stat_f_idr*(7./8.)/pow(11./4.,(4./3.)),(1./4.)) * pba->T_cmb;
    if (input_verbose > 2)
      printf("You passed N_dg = N_idr = %e, this is equivalent to xi_idr = %e in the ETHOS notation. \n", param2, pba->T_idr/pba->T_cmb);
  }
  else if (flag3 == _TRUE_) {
    pba->T_idr = param3 * pba->T_cmb;
    if (input_verbose > 1)
      printf("You passed xi_idr = %e, this is equivalent to N_idr = N_dg = %e in the NADM notation. \n", param3, stat_f_idr*pow(param3,4.)/(7./8.)*pow(11./4.,(4./3.)));
  }

  pba->Omega0_idr = stat_f_idr*pow(pba->T_idr/pba->T_cmb,4.)*pba->Omega0_g;

  Omega_tot += pba->Omega0_idr;

  /** - Omega_0_cdm (CDM) */
  class_call(parser_read_double(pfc,"Omega_cdm",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"omega_cdm",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  class_test(((flag1 == _TRUE_) && (flag2 == _TRUE_)),
             errmsg,
             "In input file, you can only enter one of Omega_cdm or omega_cdm, choose one");
  if (flag1 == _TRUE_)
    pba->Omega0_cdm = param1;
  if (flag2 == _TRUE_)
    pba->Omega0_cdm = param2/pba->h/pba->h;

  if ((ppt->gauge == synchronous) && (pba->Omega0_cdm==0)) pba->Omega0_cdm = ppr->Omega0_cdm_min_synchronous;

  Omega_tot += pba->Omega0_cdm;

  /** - Omega_0_icdm_dr (DM interacting with DR) */
  class_call(parser_read_double(pfc,"Omega_idm_dr",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"omega_idm_dr",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"f_idm_dr",&param3,&flag3,errmsg),
             errmsg,
             errmsg);
  class_test(class_at_least_two_of_three(flag1,flag2,flag3),
             errmsg,
             "In input file, you can only enter one of Omega_idm_dr, omega_idm_dr or f_idm_dr, choose one");

  /* ---> if user passes directly the density of idm_dr */
  if (flag1 == _TRUE_)
    pba->Omega0_idm_dr = param1;
  if (flag2 == _TRUE_)
    pba->Omega0_idm_dr = param2/pba->h/pba->h;

  /* ---> if user passes density of idm_dr as a fraction of the CDM one */
  if (flag3 == _TRUE_) {
    class_test((param3 < 0.) || (param3 > 1.),
               errmsg,
               "The fraction of interacting DM with DR must be between 0 and 1, you asked for f_idm_dr=%e",param3);
    class_test((param3 > 0.) && (pba->Omega0_cdm == 0.),
               errmsg,
               "If you want a fraction of interacting DM with DR, to be consistent, you should not set the fraction of CDM to zero");

    pba->Omega0_idm_dr = param3 * pba->Omega0_cdm;
    /* readjust Omega0_cdm */
    pba->Omega0_cdm -= pba->Omega0_idm_dr;
    /* to be consistent, remove same amount from Omega_tot */
    Omega_tot -= pba->Omega0_idm_dr;
    /* avoid Omega0_cdm =0 in synchronous gauge */
    if ((ppt->gauge == synchronous) && (pba->Omega0_cdm==0)) {
      pba->Omega0_cdm += ppr->Omega0_cdm_min_synchronous;
      Omega_tot += ppr->Omega0_cdm_min_synchronous;
      pba->Omega0_idm_dr -= ppr->Omega0_cdm_min_synchronous;
    }
  }

  Omega_tot += pba->Omega0_idm_dr;

  if (pba->Omega0_idm_dr > 0.) {

    class_test(pba->Omega0_idr == 0.0,
               errmsg,
               "You have requested interacting DM ith DR, this requires a non-zero density of interacting DR. Please set either N_idr or xi_idr");

    class_call(parser_read_double(pfc,"a_idm_dr",&param1,&flag1,errmsg),
               errmsg,
               errmsg);
    class_call(parser_read_double(pfc,"a_dark",&param2,&flag2,errmsg),
               errmsg,
               errmsg);
    class_call(parser_read_double(pfc,"Gamma_0_nadm",&param3,&flag3,errmsg),
               errmsg,
               errmsg);
    class_test(class_at_least_two_of_three(flag1,flag2,flag3),
               errmsg,
               "In input file, you can only enter one of a_idm_dr, a_dark or Gamma_0_nadm, choose one");

    if (flag1 == _TRUE_){
      pth->a_idm_dr = param1;
      if (input_verbose > 1)
        printf("You passed a_idm_dr = a_dark = %e, this is equivalent to Gamma_0_nadm = %e in the NADM notation. \n", param1, param1*(4./3.)*(pba->h*pba->h*pba->Omega0_idr));
    }
    else if (flag2 == _TRUE_){
      pth->a_idm_dr = param2;
      if (input_verbose > 1)
        printf("You passed a_dark = a_idm_dr = %e, this is equivalent to Gamma_0_nadm = %e in the NADM notation. \n", param2, param2*(4./3.)*(pba->h*pba->h*pba->Omega0_idr));
    }
    else if (flag3 == _TRUE_){
      pth->a_idm_dr = param3*(3./4.)/(pba->h*pba->h*pba->Omega0_idr);
      if (input_verbose > 1)
        printf("You passed Gamma_0_nadm = %e, this is equivalent to a_idm_dr = a_dark = %e in the ETHOS notation. \n", param3, pth->a_idm_dr);
    }

    /** - Load the rest of the parameters for idm and idr */

    if (flag3 == _TRUE_){ /* If the user passed Gamma_0_nadm, assume they want nadm parameterisation*/
      pth->nindex_idm_dr = 0;
      ppt->idr_nature = idr_fluid;
      if (input_verbose > 1)
        printf("NADM requested. Defaulting on nindex_idm_dr = %e and idr_nature = fluid \n", pth->nindex_idm_dr);
    }

    else{

      class_read_double_one_of_two("nindex_dark","nindex_idm_dr",pth->nindex_idm_dr);

      class_call(parser_read_string(pfc,"idr_nature",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);

      if (flag1 == _TRUE_) {
        if ((strstr(string1,"free_streaming") != NULL) || (strstr(string1,"Free_Streaming") != NULL) || (strstr(string1,"Free_streaming") != NULL) || (strstr(string1,"FREE_STREAMING") != NULL)) {
          ppt->idr_nature = idr_free_streaming;
        }
        if ((strstr(string1,"fluid") != NULL) || (strstr(string1,"Fluid") != NULL) || (strstr(string1,"FLUID") != NULL)) {
          ppt->idr_nature = idr_fluid;
        }
      }
    }

    class_read_double_one_of_two("m_idm","m_dm",pth->m_idm);

    class_read_double_one_of_two("b_dark","b_idr",pth->b_idr);

    /* Read alpha_idm_dr or alpha_dark */

    class_call(parser_read_list_of_doubles(pfc,"alpha_idm_dr",&entries_read,&(ppt->alpha_idm_dr),&flag1,errmsg),
               errmsg,
               errmsg);

    /* try with the other syntax */
    if (flag1 == _FALSE_) {
      class_call(parser_read_list_of_doubles(pfc,"alpha_dark",&entries_read,&(ppt->alpha_idm_dr),&flag1,errmsg),
                 errmsg,
                 errmsg);
    }

    if(flag1 == _TRUE_){
      if(entries_read != (ppr->l_max_idr-1)){
        class_realloc(ppt->alpha_idm_dr,ppt->alpha_idm_dr,(ppr->l_max_idr-1)*sizeof(double),errmsg);
        for(n=entries_read; n<(ppr->l_max_idr-1); n++) ppt->alpha_idm_dr[n] = ppt->alpha_idm_dr[entries_read-1];
      }
    }
    else{
      class_alloc(ppt->alpha_idm_dr,(ppr->l_max_idr-1)*sizeof(double),errmsg);
      for(n=0; n<(ppr->l_max_idr-1); n++) ppt->alpha_idm_dr[n] = 1.5;
    }

    /* Read alpha_idm_dr or alpha_dark */

    class_call(parser_read_list_of_doubles(pfc,"beta_idr",&entries_read,&(ppt->beta_idr),&flag1,errmsg),
               errmsg,
               errmsg);

    /* try with the other syntax */
    if (flag1 == _FALSE_) {
      class_call(parser_read_list_of_doubles(pfc,"beta_dark",&entries_read,&(ppt->beta_idr),&flag1,errmsg),
                 errmsg,
                 errmsg);
    }

    if(flag1 == _TRUE_){
      if(entries_read != (ppr->l_max_idr-1)){
        class_realloc(ppt->beta_idr,ppt->beta_idr,(ppr->l_max_idr-1)*sizeof(double),errmsg);
        for(n=entries_read; n<(ppr->l_max_idr-1); n++) ppt->beta_idr[n] = ppt->beta_idr[entries_read-1];
      }
    }
    else{
      class_alloc(ppt->beta_idr,(ppr->l_max_idr-1)*sizeof(double),errmsg);
      for(n=0; n<(ppr->l_max_idr-1); n++) ppt->beta_idr[n] = 1.5;
    }
  }
  else {
    ppt->alpha_idm_dr = NULL;
    ppt->beta_idr = NULL;
  }

  /** - Omega_0_dcdmdr (DCDM) */
  class_call(parser_read_double(pfc,"Omega_dcdmdr",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"omega_dcdmdr",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  class_test(((flag1 == _TRUE_) && (flag2 == _TRUE_)),
             errmsg,
             "In input file, you can only enter one of Omega_dcdmdr or omega_dcdmdr, choose one");
  if (flag1 == _TRUE_)
    pba->Omega0_dcdmdr = param1;
  if (flag2 == _TRUE_)
    pba->Omega0_dcdmdr = param2/pba->h/pba->h;

  if (pba->Omega0_dcdmdr > 0) {

    Omega_tot += pba->Omega0_dcdmdr;

    /** - Read Omega_ini_dcdm or omega_ini_dcdm */
    class_call(parser_read_double(pfc,"Omega_ini_dcdm",&param1,&flag1,errmsg),
               errmsg,
               errmsg);
    class_call(parser_read_double(pfc,"omega_ini_dcdm",&param2,&flag2,errmsg),
               errmsg,
               errmsg);
    class_test(((flag1 == _TRUE_) && (flag2 == _TRUE_)),
               errmsg,
               "In input file, you can only enter one of Omega_ini_dcdm or omega_ini_dcdm, choose one");
    if (flag1 == _TRUE_)
      pba->Omega_ini_dcdm = param1;
    if (flag2 == _TRUE_)
      pba->Omega_ini_dcdm = param2/pba->h/pba->h;

    /** - Read Gamma in same units as H0, i.e. km/(s Mpc)*/
    class_read_double("Gamma_dcdm",pba->Gamma_dcdm);
    /* Convert to Mpc */
    pba->Gamma_dcdm *= (1.e3 / _c_);

  }

  /** - non-cold relics (ncdm) */
  class_read_int("N_ncdm",N_ncdm);
  if ((flag1 == _TRUE_) && (N_ncdm > 0)){
    pba->N_ncdm = N_ncdm;

    if (ppt->gauge == synchronous)
      ppr->tol_ncdm = ppr->tol_ncdm_synchronous;
    if (ppt->gauge == newtonian)
      ppr->tol_ncdm = ppr->tol_ncdm_newtonian;

    /* Quadrature modes, 0 is qm_auto. */
    class_read_list_of_integers_or_default("Quadrature strategy",pba->ncdm_quadrature_strategy,0,N_ncdm);
    /* Number of momentum bins */
    class_read_list_of_integers_or_default("Number of momentum bins",pba->ncdm_input_q_size,-1,N_ncdm);

    /* qmax, if relevant */
    class_read_list_of_doubles_or_default("Maximum q",pba->ncdm_qmax,15,N_ncdm);

    /* Read temperatures: */
    class_read_list_of_doubles_or_default("T_ncdm",pba->T_ncdm,pba->T_ncdm_default,N_ncdm);

    /* Read chemical potentials: */
    class_read_list_of_doubles_or_default("ksi_ncdm",pba->ksi_ncdm,pba->ksi_ncdm_default,N_ncdm);

    /* Read degeneracy of each ncdm species: */
    class_read_list_of_doubles_or_default("deg_ncdm",pba->deg_ncdm,pba->deg_ncdm_default,N_ncdm);

    /* Read mass of each ncdm species: */
    class_read_list_of_doubles_or_default("m_ncdm",pba->m_ncdm_in_eV,0.0,N_ncdm);

    /* Read Omega of each ncdm species: */
    class_read_list_of_doubles_or_default("Omega_ncdm",pba->Omega0_ncdm,0.0,N_ncdm);

    /* Read omega of each ncdm species: (Use pba->M_ncdm temporarily)*/
    class_read_list_of_doubles_or_default("omega_ncdm",pba->M_ncdm,0.0,N_ncdm);

    /* Check for duplicate Omega/omega entries, missing mass definition and
       update pba->Omega0_ncdm:*/
    for(n=0; n<N_ncdm; n++){
      /* pba->M_ncdm holds value of omega */
      if (pba->M_ncdm[n]!=0.0){
        class_test(pba->Omega0_ncdm[n]!=0,errmsg,
                   "Nonzero values for both Omega and omega for ncdm species %d are specified!",n);
        pba->Omega0_ncdm[n] = pba->M_ncdm[n]/pba->h/pba->h;
      }
      if ((pba->Omega0_ncdm[n]==0.0) && (pba->m_ncdm_in_eV[n]==0.0)) {
        /* this is the right place for passing the default value of
           the mass (all parameters must have a default value; most of
           them are defined in input_default_params{}, but the ncdm mass
           is a bit special and there is no better place for setting its
           default value). We put an arbitrary value m << 10^-3 eV,
           i.e. the ultra-relativistic limit.*/
        pba->m_ncdm_in_eV[n]=1.e-5;
      }
    }

    /* Check if filenames for interpolation tables are given: */
    class_read_list_of_integers_or_default("use_ncdm_psd_files",pba->got_files,_FALSE_,N_ncdm);

    if (flag1==_TRUE_){
      for(n=0,fileentries=0; n<N_ncdm; n++){
        if (pba->got_files[n] == _TRUE_) fileentries++;
      }

      if (fileentries > 0) {

        /* Okay, read filenames.. */
        class_call(parser_read_list_of_strings(pfc,"ncdm_psd_filenames",
                                               &entries_read,&(pba->ncdm_psd_files),&flag2,errmsg),
                   errmsg,
                   errmsg);
        class_test(flag2 == _FALSE_,errmsg,
                   "Input use_ncdm_files is found, but no filenames found!");
        class_test(entries_read != fileentries,errmsg,
                   "Number of filenames found, %d, does not match number of _TRUE_ values in use_ncdm_files, %d",
                   entries_read,fileentries);
      }
    }
    /* Read (optional) p.s.d.-parameters:*/
    parser_read_list_of_doubles(pfc,
                                "ncdm_psd_parameters",
                                &entries_read,
                                &(pba->ncdm_psd_parameters),
                                &flag2,
                                errmsg);

    class_call(background_ncdm_init(ppr,pba),
               pba->error_message,
               errmsg);

    /* We must calculate M from omega or vice versa if one of them is missing.
       If both are present, we must update the degeneracy parameter to
       reflect the implicit normalization of the distribution function.*/
    for (n=0; n < N_ncdm; n++){
      if (pba->m_ncdm_in_eV[n] != 0.0){
        /* Case of only mass or mass and Omega/omega: */
        pba->M_ncdm[n] = pba->m_ncdm_in_eV[n]/_k_B_*_eV_/pba->T_ncdm[n]/pba->T_cmb;
        class_call(background_ncdm_momenta(pba->q_ncdm_bg[n],
                                           pba->w_ncdm_bg[n],
                                           pba->q_size_ncdm_bg[n],
                                           pba->M_ncdm[n],
                                           pba->factor_ncdm[n],
                                           0.,
                                           NULL,
                                           &rho_ncdm,
                                           NULL,
                                           NULL,
                                           NULL),
                   pba->error_message,
                   errmsg);
        if (pba->Omega0_ncdm[n] == 0.0){
          pba->Omega0_ncdm[n] = rho_ncdm/pba->H0/pba->H0;
        }
        else{
          fnu_factor = (pba->H0*pba->H0*pba->Omega0_ncdm[n]/rho_ncdm);
          pba->factor_ncdm[n] *= fnu_factor;
          /* dlnf0dlnq is already computed, but it is
             independent of any normalization of f0.
             We don't need the factor anymore, but we
             store it nevertheless:*/
          pba->deg_ncdm[n] *=fnu_factor;
        }
      }
      else{
        /* Case of only Omega/omega: */
        class_call(background_ncdm_M_from_Omega(ppr,pba,n),
                   pba->error_message,
                   errmsg);
        //printf("M_ncdm:%g\n",pba->M_ncdm[n]);
        pba->m_ncdm_in_eV[n] = _k_B_/_eV_*pba->T_ncdm[n]*pba->M_ncdm[n]*pba->T_cmb;
      }
      pba->Omega0_ncdm_tot += pba->Omega0_ncdm[n];
      //printf("Adding %g to total Omega..\n",pba->Omega0_ncdm[n]);
    }
  }
  Omega_tot += pba->Omega0_ncdm_tot;

  /** - Omega_0_k (effective fractional density of curvature) */
  class_read_double("Omega_k",pba->Omega0_k);
  /** - Set curvature parameter K */
  pba->K = -pba->Omega0_k*pow(pba->a_today*pba->H0,2);
  /** - Set curvature sign */
  if (pba->K > 0.) pba->sgnK = 1;
  else if (pba->K < 0.) pba->sgnK = -1;

  /** - Omega_0_lambda (cosmological constant), Omega0_fld (dark energy fluid), Omega0_scf (scalar field) */

  class_call(parser_read_double(pfc,"Omega_Lambda",&param1,&flag1,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"Omega_fld",&param2,&flag2,errmsg),
             errmsg,
             errmsg);
  class_call(parser_read_double(pfc,"Omega_scf",&param3,&flag3,errmsg),
             errmsg,
             errmsg);

  class_test((flag1 == _TRUE_) && (flag2 == _TRUE_) && ((flag3 == _FALSE_) || (param3 >= 0.)),
             errmsg,
             "In input file, either Omega_Lambda or Omega_fld must be left unspecified, except if Omega_scf is set and <0.0, in which case the contribution from the scalar field will be the free parameter.");

  /** - --> (flag3 == _FALSE_) || (param3 >= 0.) explained:
   *  it means that either we have not read Omega_scf so we are ignoring it
   *  (unlike lambda and fld!) OR we have read it, but it had a
   *  positive value and should not be used for filling.
   *  We now proceed in two steps:
   *  1) set each Omega0 and add to the total for each specified component.
   *  2) go through the components in order {lambda, fld, scf} and
   *     fill using first unspecified component.
   */

  /* Step 1 */
  if (flag1 == _TRUE_){
    pba->Omega0_lambda = param1;
    Omega_tot += pba->Omega0_lambda;
  }
  if (flag2 == _TRUE_){
    pba->Omega0_fld = param2;
    Omega_tot += pba->Omega0_fld;
  }
  if ((flag3 == _TRUE_) && (param3 >= 0.)){
    pba->Omega0_scf = param3;
    Omega_tot += pba->Omega0_scf;
  }
  /* Step 2 */
  if (flag1 == _FALSE_) {
    //Fill with Lambda
    pba->Omega0_lambda= 1. - pba->Omega0_k - Omega_tot;
    if (input_verbose > 0) printf(" -> matched budget equations by adjusting Omega_Lambda = %e\n",pba->Omega0_lambda);
  }
  else if (flag2 == _FALSE_) {
    // Fill up with fluid
    pba->Omega0_fld = 1. - pba->Omega0_k - Omega_tot;
    if (input_verbose > 0) printf(" -> matched budget equations by adjusting Omega_fld = %e\n",pba->Omega0_fld);
  }
  else if ((flag3 == _TRUE_) && (param3 < 0.)){
    // Fill up with scalar field
    pba->Omega0_scf = 1. - pba->Omega0_k - Omega_tot;
    if (input_verbose > 0) printf(" -> matched budget equations by adjusting Omega_scf = %e\n",pba->Omega0_scf);
  }

  /*
    fprintf(stderr,"%e %e %e %e %e\n",
    pba->Omega0_lambda,
    pba->Omega0_fld,
    pba->Omega0_scf,
    pba->Omega0_k,
    Omega_tot);
  */

  /** - Test that the user have not specified Omega_scf = -1 but left either
      Omega_lambda or Omega_fld unspecified:*/
  class_test(((flag1 == _FALSE_)||(flag2 == _FALSE_)) && ((flag3 == _TRUE_) && (param3 < 0.)),
             errmsg,
             "It looks like you want to fulfil the closure relation sum Omega = 1 using the scalar field, so you have to specify both Omega_lambda and Omega_fld in the .ini file");

  if (pba->Omega0_fld != 0.) {

    class_call(parser_read_string(pfc,
                                  "use_ppf",
                                  &string1,
                                  &flag1,
                                  errmsg),
               errmsg,
               errmsg);

    if (flag1 == _TRUE_){
      if((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)){
        pba->use_ppf = _TRUE_;
        class_read_double("c_gamma_over_c_fld",pba->c_gamma_over_c_fld);
      }
      else {
        pba->use_ppf = _FALSE_;
      }
    }

    class_call(parser_read_string(pfc,"fluid_equation_of_state",&string1,&flag1,errmsg),
               errmsg,
               errmsg);

    if (flag1 == _TRUE_) {

      if ((strstr(string1,"CLP") != NULL) || (strstr(string1,"clp") != NULL)) {
        pba->fluid_equation_of_state = CLP;
      }

      else if ((strstr(string1,"EDE") != NULL) || (strstr(string1,"ede") != NULL)) {
        pba->fluid_equation_of_state = EDE;
      }

      else {
        class_stop(errmsg,"incomprehensible input '%s' for the field 'fluid_equation_of_state'",string1);
      }
    }

    if (pba->fluid_equation_of_state == CLP) {
      class_read_double("w0_fld",pba->w0_fld);
      class_read_double("wa_fld",pba->wa_fld);
      class_read_double("cs2_fld",pba->cs2_fld);
    }

    if (pba->fluid_equation_of_state == EDE) {
      class_read_double("w0_fld",pba->w0_fld);
      class_read_double("Omega_EDE",pba->Omega_EDE);
      class_read_double("cs2_fld",pba->cs2_fld);
    }
  }

  /* Additional SCF parameters: */
  if (pba->Omega0_scf != 0.){
    /** - Read parameters describing scalar field potential */
    class_call(parser_read_list_of_doubles(pfc,
                                           "scf_parameters",
                                           &(pba->scf_parameters_size),
                                           &(pba->scf_parameters),
                                           &flag1,
                                           errmsg),
               errmsg,errmsg);
    class_read_int("scf_tuning_index",pba->scf_tuning_index);
    class_test(pba->scf_tuning_index >= pba->scf_parameters_size,
               errmsg,
               "Tuning index scf_tuning_index = %d is larger than the number of entries %d in scf_parameters. Check your .ini file.",pba->scf_tuning_index,pba->scf_parameters_size);
    /** - Assign shooting parameter */
    class_read_double("scf_shooting_parameter",pba->scf_parameters[pba->scf_tuning_index]);

    scf_lambda = pba->scf_parameters[0];
    if ((fabs(scf_lambda) <3.)&&(pba->background_verbose>1))
      printf("lambda = %e <3 won't be tracking (for exp quint) unless overwritten by tuning function\n",scf_lambda);

    class_call(parser_read_string(pfc,
                                  "attractor_ic_scf",
                                  &string1,
                                  &flag1,
                                  errmsg),
               errmsg,
               errmsg);

    if (flag1 == _TRUE_){
      if((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)){
        pba->attractor_ic_scf = _TRUE_;
      }
      else{
        pba->attractor_ic_scf = _FALSE_;
        class_test(pba->scf_parameters_size<2,
                   errmsg,
                   "Since you are not using attractor initial conditions, you must specify phi and its derivative phi' as the last two entries in scf_parameters. See explanatory.ini for more details.");
        pba->phi_ini_scf = pba->scf_parameters[pba->scf_parameters_size-2];
        pba->phi_prime_ini_scf = pba->scf_parameters[pba->scf_parameters_size-1];
      }
    }
  }

  /** (b) assign values to thermodynamics cosmological parameters */

  /** - primordial helium fraction */
  class_call(parser_read_string(pfc,"YHe",&string1,&flag1,errmsg),
             errmsg,
             errmsg);

  if (flag1 == _TRUE_) {

    if ((strstr(string1,"BBN") != NULL) || (strstr(string1,"bbn") != NULL)) {
      pth->YHe = _BBN_;
    }
    else {
      class_read_double("YHe",pth->YHe);
    }

  }

  /** - recombination parameters */
  class_call(parser_read_string(pfc,"recombination",&string1,&flag1,errmsg),
             errmsg,
             errmsg);

  if (flag1 == _TRUE_) {

    if ((strstr(string1,"HYREC") != NULL) || (strstr(string1,"hyrec") != NULL) || (strstr(string1,"HyRec") != NULL)) {
      pth->recombination = hyrec;
    }

  }

  /** - reionization parametrization */
  class_call(parser_read_string(pfc,"reio_parametrization",&string1,&flag1,errmsg),
             errmsg,
             errmsg);

  if (flag1 == _TRUE_) {
    flag2=_FALSE_;
    if (strcmp(string1,"reio_none") == 0) {
      pth->reio_parametrization=reio_none;
      flag2=_TRUE_;
    }
    if (strcmp(string1,"reio_camb") == 0) {
      pth->reio_parametrization=reio_camb;
      flag2=_TRUE_;
    }
    if (strcmp(string1,"reio_bins_tanh") == 0) {
      pth->reio_parametrization=reio_bins_tanh;
      flag2=_TRUE_;
    }
    if (strcmp(string1,"reio_half_tanh") == 0) {
      pth->reio_parametrization=reio_half_tanh;
      flag2=_TRUE_;
    }
    if (strcmp(string1,"reio_many_tanh") == 0) {
      pth->reio_parametrization=reio_many_tanh;
      flag2=_TRUE_;
    }
    if (strcmp(string1,"reio_inter") == 0) {
      pth->reio_parametrization=reio_inter;
      flag2=_TRUE_;
    }

    class_test(flag2==_FALSE_,
               errmsg,
               "could not identify reionization_parametrization value, check that it is one of 'reio_none', 'reio_camb', 'reio_bins_tanh', 'reio_half_tanh', 'reio_many_tanh', 'reio_inter'...");
  }

  /** - reionization parameters if reio_parametrization=reio_camb */
  if ((pth->reio_parametrization == reio_camb) || (pth->reio_parametrization == reio_half_tanh)){
    class_call(parser_read_double(pfc,"z_reio",&param1,&flag1,errmsg),
               errmsg,
               errmsg);
    class_call(parser_read_double(pfc,"tau_reio",&param2,&flag2,errmsg),
               errmsg,
               errmsg);
    class_test(((flag1 == _TRUE_) && (flag2 == _TRUE_)),
               errmsg,
               "In input file, you can only enter one of z_reio or tau_reio, choose one");
    if (flag1 == _TRUE_) {
      pth->z_reio=param1;
      pth->reio_z_or_tau=reio_z;
    }
    if (flag2 == _TRUE_) {
      pth->tau_reio=param2;
      pth->reio_z_or_tau=reio_tau;

    }

    class_read_double("reionization_exponent",pth->reionization_exponent);
    class_read_double("reionization_width",pth->reionization_width);
    class_read_double("helium_fullreio_redshift",pth->helium_fullreio_redshift);
    class_read_double("helium_fullreio_width",pth->helium_fullreio_width);

  }

  /** - reionization parameters if reio_parametrization=reio_bins_tanh */
  if (pth->reio_parametrization == reio_bins_tanh) {
    class_read_int("binned_reio_num",pth->binned_reio_num);
    class_read_list_of_doubles("binned_reio_z",pth->binned_reio_z,pth->binned_reio_num);
    class_read_list_of_doubles("binned_reio_xe",pth->binned_reio_xe,pth->binned_reio_num);
    class_read_double("binned_reio_step_sharpness",pth->binned_reio_step_sharpness);
  }

  /** - reionization parameters if reio_parametrization=reio_many_tanh */
  if (pth->reio_parametrization == reio_many_tanh) {
    class_read_int("many_tanh_num",pth->many_tanh_num);
    class_read_list_of_doubles("many_tanh_z",pth->many_tanh_z,pth->many_tanh_num);
    class_read_list_of_doubles("many_tanh_xe",pth->many_tanh_xe,pth->many_tanh_num);
    class_read_double("many_tanh_width",pth->many_tanh_width);
  }

  /** - reionization parameters if reio_parametrization=reio_many_tanh */
  if (pth->reio_parametrization == reio_inter) {
    class_read_int("reio_inter_num",pth->reio_inter_num);
    class_read_list_of_doubles("reio_inter_z",pth->reio_inter_z,pth->reio_inter_num);
    class_read_list_of_doubles("reio_inter_xe",pth->reio_inter_xe,pth->reio_inter_num);
  }

  /** - energy injection parameters from CDM annihilation/decay */

  class_read_double("annihilation",pth->annihilation);

  if (pth->annihilation > 0.) {

    class_read_double("annihilation_variation",pth->annihilation_variation);
    class_read_double("annihilation_z",pth->annihilation_z);
    class_read_double("annihilation_zmax",pth->annihilation_zmax);
    class_read_double("annihilation_zmin",pth->annihilation_zmin);
    class_read_double("annihilation_f_halo",pth->annihilation_f_halo);
    class_read_double("annihilation_z_halo",pth->annihilation_z_halo);

    class_call(parser_read_string(pfc,
                                  "on the spot",
                                  &(string1),
                                  &(flag1),
                                  errmsg),
               errmsg,
               errmsg);

    if (flag1 == _TRUE_) {
      if ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)) {
        pth->has_on_the_spot = _TRUE_;
      }
      else {
        if ((strstr(string1,"n") != NULL) || (strstr(string1,"N") != NULL)) {
          pth->has_on_the_spot = _FALSE_;
        }
        else {
          class_stop(errmsg,"incomprehensible input '%s' for the field 'on the spot'",string1);
        }
      }
    }
  }

  class_read_double("decay",pth->decay);

  class_call(parser_read_string(pfc,
                                "compute damping scale",
                                &(string1),
                                &(flag1),
                                errmsg),
             errmsg,
             errmsg);

  if (flag1 == _TRUE_) {
    if ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)) {
      pth->compute_damping_scale = _TRUE_;
    }
    else {
      if ((strstr(string1,"n") != NULL) || (strstr(string1,"N") != NULL)) {
        pth->compute_damping_scale = _FALSE_;
      }
      else {
        class_stop(errmsg,"incomprehensible input '%s' for the field 'compute damping scale'",string1);
      }
    }
  }

  /** (c) define which perturbations and sources should be computed, and down to which scale */

  ppt->has_perturbations = _FALSE_;
  ppt->has_cls = _FALSE_;

  class_call(parser_read_string(pfc,"output",&string1,&flag1,errmsg),
             errmsg,
             errmsg);

  if (flag1 == _TRUE_) {

    if ((strstr(string1,"tCl") != NULL) || (strstr(string1,"TCl") != NULL) || (strstr(string1,"TCL") != NULL)) {
      ppt->has_cl_cmb_temperature = _TRUE_;
      ppt->has_perturbations = _TRUE_;
      ppt->has_cls = _TRUE_;
    }

    if ((strstr(string1,"pCl") != NULL) || (strstr(string1,"PCl") != NULL) || (strstr(string1,"PCL") != NULL)) {
      ppt->has_cl_cmb_polarization = _TRUE_;
      ppt->has_perturbations = _TRUE_;
      ppt->has_cls = _TRUE_;
    }

    if ((strstr(string1,"lCl") != NULL) || (strstr(string1,"LCl") != NULL) || (strstr(string1,"LCL") != NULL)) {
      ppt->has_cl_cmb_lensing_potential = _TRUE_;
      ppt->has_perturbations = _TRUE_;
      ppt->has_cls = _TRUE_;
    }

    if ((strstr(string1,"nCl") != NULL) || (strstr(string1,"NCl") != NULL) || (strstr(string1,"NCL") != NULL) ||
        (strstr(string1,"dCl") != NULL) || (strstr(string1,"DCl") != NULL) || (strstr(string1,"DCL") != NULL)) {
      ppt->has_cl_number_count = _TRUE_;
      ppt->has_perturbations = _TRUE_;
      ppt->has_cls = _TRUE_;
    }

    if ((strstr(string1,"sCl") != NULL) || (strstr(string1,"SCl") != NULL) || (strstr(string1,"SCL") != NULL)) {
      ppt->has_cl_lensing_potential=_TRUE_;
      ppt->has_perturbations = _TRUE_;
      ppt->has_cls = _TRUE_;
    }

    if ((strstr(string1,"mPk") != NULL) || (strstr(string1,"MPk") != NULL) || (strstr(string1,"MPK") != NULL)) {
      ppt->has_pk_matter=_TRUE_;
      ppt->has_perturbations = _TRUE_;

      /*if (pba->Omega0_ncdm_tot != 0.0){
        class_call(parser_read_string(pfc,"pk_only_cdm_bar",&string1,&flag1,errmsg),
        errmsg,
        errmsg);
        if (flag1 == _TRUE_){
        if((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)){
        ppt->pk_only_cdm_bar = _TRUE_;
        }
        else {
        ppt->pk_only_cdm_bar = _FALSE_;
        }
        }
        }*/

    }

      /** (sz) SZ parameters */

      //BB: read SZ parameters from ini file
      class_read_int("class_sz_verbose",pclass_sz->sz_verbose);
      class_read_int("nlSZ",pclass_sz->nlSZ);
      class_read_int("convert_cls_to_gamma",pclass_sz->convert_cls_to_gamma);

      class_read_int("N_kSZ2_gal_multipole_grid",pclass_sz->N_kSZ2_gal_multipole_grid);
      class_read_int("N_kSZ2_gal_theta_grid",pclass_sz->N_kSZ2_gal_theta_grid);


      class_read_double("ell_max",pclass_sz->ell_max_mock);
      class_read_double("ell_min",pclass_sz->ell_min_mock);
      class_read_double("dlogell",pclass_sz->dlogell);
      class_read_double("dell",pclass_sz->dell);

      class_read_double("freq_max",pclass_sz->freq_max);
      class_read_double("freq_min",pclass_sz->freq_min);
      class_read_double("dlogfreq",pclass_sz->dlogfreq);
      class_read_double("dfreq",pclass_sz->dfreq);

      class_read_double("k_min_for_pk_hm",pclass_sz->k_min_for_pk_hm);
      class_read_double("k_max_for_pk_hm",pclass_sz->k_max_for_pk_hm);
      class_read_double("dlnk_for_pk_hm",pclass_sz->dlnk_for_pk_hm);

      class_read_double("z_for_pk_hm",pclass_sz->z_for_pk_hm);

      class_read_double("kstar_damping_1h_term (1/Mpc)",pclass_sz->kstar_damping_1h_term_Mpc);
      class_read_int("damping_1h_term",pclass_sz->damping_1h_term);
      class_read_double("photo_z_params",pclass_sz->photo_z_params);
      class_read_double("dndz_shift_source_gal",pclass_sz->dndz_shift_source_gal);
      class_read_double("dndz_shift_gal",pclass_sz->dndz_shift_gal);
      class_read_double("dndz_stretch_source_gal",pclass_sz->dndz_stretch_source_gal);
      class_read_double("dndz_stretch_gal",pclass_sz->dndz_stretch_gal);
      class_read_double("shear_calibration_m",pclass_sz->shear_calibration_m);


      //Redshift limits for the integration
      class_read_double("z_min",pclass_sz->z1SZ);
      class_read_double("z_max",pclass_sz->z2SZ);
      pclass_sz->z1SZ_dndlnM = pclass_sz->z1SZ;
      pclass_sz->z2SZ_dndlnM = pclass_sz->z2SZ;


      class_read_double("redshift_epsrel",pclass_sz->redshift_epsrel);
      class_read_double("redshift_epsabs",pclass_sz->redshift_epsabs);

      class_read_double("L_sat_epsabs",pclass_sz->epsabs_L_sat);
      class_read_double("L_sat_epsrel",pclass_sz->epsrel_L_sat);

      class_read_double("max redshift for cluster counts",pcsz->z_max);

      class_read_double("shape_noise_siggamma2",pclass_sz->shape_noise_siggamma2);
      class_read_double("ns_gal_per_arcmin2",pclass_sz->ns_gal_per_arcmin2);
      class_read_double("cl_gal_gal_A_sn",pclass_sz->cl_gal_gal_A_sn);
      class_read_double("csat_over_cdm",pclass_sz->csat_over_cdm);


      //for tabulation of sigma(m,z)
      class_read_int("ndim_redshifts",pclass_sz->ndim_redshifts);//number of z in the interpolation for sigma
      pclass_sz->n_z_dndlnM = pclass_sz->ndim_redshifts;
      // for tabulation of hmf:
      class_read_int("n_z_dndlnM",pclass_sz->n_z_dndlnM);


      class_read_int("has_b_custom1",pclass_sz->has_b_custom1);

      // class_read_int("ndim_redshifts_for_integral",pclass_sz->ndim_redshifts_for_integral);//number of z in the integration

      // class_read_int("n_m_matter_density_profile",pclass_sz->n_m_matter_density_profile);
      // printf("%d \n",pclass_sz->n_m_matter_density_profile);
      // exit(0);


      //mass limits: h^-1 Msun
      class_read_double("M_min",pclass_sz->M1SZ);
      class_read_double("M_max",pclass_sz->M2SZ);
      pclass_sz->M1SZ_dndlnM = pclass_sz->M1SZ;
      pclass_sz->M2SZ_dndlnM = pclass_sz->M2SZ;

     pclass_sz->m_min_counter_terms = pclass_sz->M1SZ;
     pclass_sz->m_max_counter_terms = pclass_sz->M2SZ;  // default M2SZ, if passed, it will be set to something different and fixed, e.g.,
     class_read_double("m_min_counter_terms",pclass_sz->m_min_counter_terms );
     class_read_double("m_max_counter_terms",pclass_sz->m_max_counter_terms );



     class_read_int("include_y_counterterms_in_yk",pclass_sz->include_y_counterterms_in_yk);
     class_read_int("include_g_counterterms_in_gk",pclass_sz->include_g_counterterms_in_gk);
     class_read_int("include_gk_counterterms_in_gk",pclass_sz->include_gk_counterterms_in_gk);
     class_read_int("include_k_counterterms_in_gk",pclass_sz->include_k_counterterms_in_gk);


      pclass_sz->M_min_ng_bar = pclass_sz->M1SZ;
      pclass_sz->M_max_ng_bar = pclass_sz->M2SZ;

      class_read_double("M_min for ng_bar",pclass_sz->M_min_ng_bar);
      class_read_double("M_max for ng_bar",pclass_sz->M_max_ng_bar);

      //mass limits: h^-1 Msun
      class_read_double("M1SZ_dndlnM",pclass_sz->M1SZ_dndlnM);
      class_read_double("M2SZ_dndlnM",pclass_sz->M2SZ_dndlnM);

      class_read_double("z1SZ_dndlnM",pclass_sz->z1SZ_dndlnM);
      class_read_double("z2SZ_dndlnM",pclass_sz->z2SZ_dndlnM);


      class_read_double("theta_ej_bcm",pclass_sz->theta_ej_bcm);
      class_read_double("delta_bcm",pclass_sz->delta_bcm);
      class_read_double("gamma_bcm",pclass_sz->gamma_bcm);
      class_read_double("mu_bcm",pclass_sz->mu_bcm);
      class_read_double("eta_star_bcm",pclass_sz->eta_star_bcm);
      class_read_double("log10Mc_bcm",pclass_sz->log10Mc_bcm);
      class_read_double("nu_log10Mc_bcm",pclass_sz->nu_log10Mc_bcm);



      //number of mass bins for cov_Y-N:
      class_read_double("number of mass bins for cluster covariances",pclass_sz->nbins_M);

      //Pressure profile is considered between x_in and x_out
      class_read_double("x_inSZ",pclass_sz->x_inSZ);
      class_read_double("x_outSZ",pclass_sz->x_outSZ);

      class_read_double("x_out_custom1",pclass_sz->x_out_custom1);

      class_read_double("delta_alpha",pclass_sz->delta_alpha);
      class_read_double("alpha_p",pclass_sz->alpha_p);

      //Hydrostatic Equilibrium Mass Bias, Piffaretti & Valdarnini [arXiv:0808.1111]

      class_read_double("B",pclass_sz->HSEbias);
      class_read_double("bias_sz",pclass_sz->HSEbias); // soliket/clusters.py notations

      class_read_double("Ap",pclass_sz->Ap);
      class_read_double("alpha_b",pclass_sz->alpha_b);
      class_read_int("mass_dependent_bias",pclass_sz->mass_dependent_bias);

      class_read_int("experiment",pclass_sz->experiment);
      class_read_int("use_planck_binned_proba",pclass_sz->use_planck_binned_proba);
      class_read_double("bin_z_min_cluster_counts",pclass_sz->bin_z_min_cluster_counts);
      class_read_double("bin_z_max_cluster_counts",pclass_sz->bin_z_max_cluster_counts);
      class_read_double("bin_dz_cluster_counts",pclass_sz->bin_dz_cluster_counts);

      class_read_double("bin_dlog10_snr",pclass_sz->bin_dlog10_snr);
      class_read_double("log10_snr_min",pclass_sz->log10_snr_min);
      class_read_double("log10_snr_max",pclass_sz->log10_snr_max);

      class_read_double("lnymin",pclass_sz->lnymin);
      class_read_double("lnymax",pclass_sz->lnymax);
      class_read_double("dlny",pclass_sz->dlny);

      class_read_double("dlnM_cluster_count_completeness_grid",pclass_sz->dlnM_cluster_count_completeness_grid);

      class_read_double("dz_cluster_count_completeness_grid_low_z",pclass_sz->dz_cluster_count_completeness_grid_low_z);
      class_read_double("dz_cluster_count_completeness_grid_mid_z",pclass_sz->dz_cluster_count_completeness_grid_mid_z);
      class_read_double("dz_cluster_count_completeness_grid_high_z",pclass_sz->dz_cluster_count_completeness_grid_high_z);

      class_read_double("cluster_count_completeness_grid_z_cutoff_low",pclass_sz->cluster_count_completeness_grid_z_cutoff_low);
      class_read_double("cluster_count_completeness_grid_z_cutoff_mid",pclass_sz->cluster_count_completeness_grid_z_cutoff_mid);

      class_read_double("mass_epsrel_cluster_counts",pclass_sz->mass_epsrel_cluster_counts);
      class_read_double("mass_epsabs_cluster_counts",pclass_sz->mass_epsabs_cluster_counts);
      class_read_double("redshift_epsrel_cluster_counts",pclass_sz->redshift_epsrel_cluster_counts);
      class_read_double("redshift_epsabs_cluster_counts",pclass_sz->redshift_epsabs_cluster_counts);



      class_read_double("sky area in deg2",pclass_sz->sky_area_deg2);

      class_read_int("apply_relativistic_correction_to_y_m",pclass_sz->apply_relativistic_correction_to_y_m);
      class_read_int("has_selection_function",pcsz->has_completeness);
      class_read_int("use_m500c_in_ym_relation",pclass_sz->use_m500c_in_ym_relation);
      class_read_int("use_m200c_in_ym_relation",pclass_sz->use_m200c_in_ym_relation);
      class_read_int("mass_range",pcsz->mass_range);


      class_read_int("use_maniyar_cib_model",pclass_sz->use_maniyar_cib_model);

      class_read_int("y_m_relation",pclass_sz->y_m_relation);
      // class_read_double("ystar",pcsz->ystar);
      // class_read_double("alpha",pcsz->alpha);
      // class_read_double("sigmaM",pcsz->sigmaM);
      //pcsz->ystar = pow(10.,pcsz->ystar)/pow(2., pcsz->alpha)*0.00472724; ////8.9138435358806980e-004;

      class_read_double("ystar_ym",pclass_sz->ystar_ym);
      class_read_double("alpha_ym",pclass_sz->alpha_ym);
      class_read_double("beta_ym",pclass_sz->beta_ym);




      // class_read_double("ystar_nika2",pclass_sz->ystar_nika2);
      // class_read_double("alpha_nika2",pclass_sz->alpha_nika2);
      // class_read_double("beta_nika2",pclass_sz->beta_nika2);
      // class_read_double("sigmaM_nika2",pclass_sz->sigmaM_nika2);

      class_read_double("alpha_thm",pclass_sz->alpha_theta);

      class_read_double("sigmaM_ym",pclass_sz->sigmaM_ym);
      class_read_double("scatter_sz",pclass_sz->sigmaM_ym); // soliket/clusters.py notations

      class_read_double("A_ym",pclass_sz->A_ym);
      class_read_double("tenToA0",pclass_sz->A_ym); // soliket/clusters.py notations

      class_read_double("B_ym",pclass_sz->B_ym);
      class_read_double("B0",pclass_sz->B_ym); // soliket/clusters.py notations

      class_read_double("C_ym",pclass_sz->C_ym);
      class_read_double("C0",pclass_sz->C_ym); // soliket/clusters.py notations

      class_read_double("m_pivot_ym_[Msun]",pclass_sz->m_pivot_ym);

      //For the tabulation of sigma2
      class_read_int("ndim_masses",pclass_sz->ndim_masses);
      pclass_sz->n_m_dndlnM = pclass_sz->ndim_masses;

      // for the tabulation of hmf itself
      class_read_int("n_m_dndlnM",pclass_sz->n_m_dndlnM);

      class_read_int("use_class_sz_fast_mode",pclass_sz->use_class_sz_fast_mode);

      class_read_double("delta_cSZ",pclass_sz->delta_cSZ);

      //Multplicity function Tinker 2010

      class_read_double("alphaSZ",pclass_sz->alphaSZ);
      class_read_double("beta0SZ",pclass_sz->beta0SZ);
      class_read_double("gamma0SZ",pclass_sz->gamma0SZ);

      class_read_double("phi0SZ",pclass_sz->phi0SZ);
      class_read_double("eta0SZ",pclass_sz->eta0SZ);

      class_read_int("no_b2",pclass_sz->no_b2);

      //Multplicity function Bocquet 2015

      class_read_double("Ap0",pclass_sz->Ap0);
      class_read_double("a0",pclass_sz->a0);
      class_read_double("b0",pclass_sz->b0);
      class_read_double("c0",pclass_sz->c0);



        class_call(parser_read_string(pfc,"component of tSZ power spectrum",&string1,&flag1,errmsg),
                   errmsg,
                   errmsg);
        if (flag1 == _TRUE_) {
          if ((strstr(string1,"total") != NULL)){
            pclass_sz->which_ps_sz = 0;
            pclass_sz->has_completeness_for_ps_SZ = 0;
          }
          else if ((strstr(string1,"unresolved") != NULL)){
            pclass_sz->which_ps_sz = 2;
            pclass_sz->has_completeness_for_ps_SZ = 1;
          }
          else if ((strstr(string1,"resolved") != NULL)){
            pclass_sz->which_ps_sz = 1;
            pclass_sz->has_completeness_for_ps_SZ = 1;
          }

        }

      class_read_int("load_Planck_noise_map",pclass_sz->has_completeness_for_ps_SZ);

      class_read_double("signal-to-noise_cut-off_for_survey_cluster_completeness",pclass_sz->sn_cutoff);
      class_read_double("signal-to-noise_cut-off_for_survey_cluster_completeness",pcsz->sn_cutoff);

      //Integration scheme for the mass integral:
      class_call(parser_read_string(pfc,"integration method (mass)",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      if (flag1 == _TRUE_) {
        if ((strstr(string1,"patterson") != NULL)){
          pclass_sz->integration_method_mass=0;

          class_read_int("patterson_show_neval",pclass_sz->patterson_show_neval);
        }
        else if ((strstr(string1,"gsl_qags") != NULL)){
          pclass_sz->integration_method_mass=1;
        }
        else if ((strstr(string1,"gsl_qag") != NULL)){
          pclass_sz->integration_method_mass=2;
        }
        else if ((strstr(string1,"gsl_romberg") != NULL)){
          pclass_sz->integration_method_mass=3;
        }
        else if ((strstr(string1,"gsl_qng") != NULL)){
          pclass_sz->integration_method_mass=4;
        }
      }

      class_read_double("mass_epsrel",pclass_sz->mass_epsrel);
      class_read_double("mass_epsabs",pclass_sz->mass_epsabs);

      class_read_double("redshift_for_dndm",pcsz->redshift_for_dndm);
      class_read_double("size_logM_for_dndm",pcsz->size_logM);

      class_read_double("f_sky",pclass_sz->f_sky);

      class_read_int("bispec_conf_id",pclass_sz->bispec_conf_id);


      //Foreground Nuisance parameters
      class_read_double("A_cib",pclass_sz->A_cib);
      class_read_double("A_rs",pclass_sz->A_rs);
      class_read_double("A_ir",pclass_sz->A_ir);
      class_read_double("A_cn",pclass_sz->A_cn);


      class_read_double("k_per_decade_class_sz",pclass_sz->k_per_decade_for_tSZ);
      class_read_double("k_min_for_pk_class_sz",pclass_sz->k_min_for_pk_in_tSZ);
      class_read_double("k_max_for_pk_class_sz",pclass_sz->k_max_for_pk_in_tSZ);


      //BB: read the quantities to be computed by class_sz
      class_call(parser_read_string(pfc,"output",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);

      if ((strstr(string1,"tSZ_tSZ_1h") != NULL) || (strstr(string1,"tSZ_1h") != NULL) || (strstr(string1,"tSZCl") != NULL) || (strstr(string1,"tszCL") != NULL)) {
        pclass_sz->has_sz_ps =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;

      }

      if ((strstr(string1,"tSZ_2h") != NULL) || (strstr(string1,"tSZ_tSZ_2h") != NULL)) {
        pclass_sz->has_sz_2halo =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"m_y_y_2h") != NULL) ) {
        pclass_sz->has_sz_m_y_y_2h =_TRUE_;
        pclass_sz->has_sz_2halo =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }
      if ((strstr(string1,"m_y_y_1h") != NULL) ) {
        pclass_sz->has_sz_m_y_y_1h =_TRUE_;
        pclass_sz->has_sz_ps =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }
      if ((strstr(string1,"te_y_y") != NULL) ) {
        pclass_sz->has_sz_te_y_y =_TRUE_;
        pclass_sz->has_sz_ps =_TRUE_; //ps is necessary in this case (Te = "Te_y_y/y_y")
        pclass_sz->has_sz_2halo =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"tSZ_cov_Y_N") != NULL) ) {
        pclass_sz->has_sz_cov_Y_N =_TRUE_;
        pclass_sz->has_sz_cov_N_N =_TRUE_;
        pclass_sz->has_sz_ps =_TRUE_;
        pclass_sz->has_sz_trispec =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"tSZ_cov_hsv") != NULL) ) {
        pclass_sz->has_sz_cov_Y_N_next_order =_TRUE_;
        pclass_sz->has_sz_cov_N_N_hsv =_TRUE_;
        pclass_sz->has_sigma2_hsv = _TRUE_;
        pclass_sz->has_sz_cov_Y_N =_TRUE_;
        pclass_sz->has_sz_cov_Y_Y_ssc =_TRUE_;
        pclass_sz->has_sz_cov_N_N =_TRUE_;
        pclass_sz->has_sz_ps =_TRUE_;
        pclass_sz->has_sz_trispec =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"cov(N,N)") != NULL) ) {
        pclass_sz->has_sz_cov_N_N_hsv =_TRUE_;
        pclass_sz->has_sigma2_hsv = _TRUE_;
        pclass_sz->has_sz_cov_N_N =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"tSZ_Trispectrum") != NULL) ) {
        pclass_sz->has_sz_trispec =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"hmf") != NULL) ) {
        pclass_sz->has_hmf =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        // pclass_sz->has_completeness_for_ps_SZ = 1;
      }
      if ((strstr(string1,"mean_y") != NULL) ) {
        pclass_sz->has_mean_y =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"pk_gg_at_z_1h") != NULL) ) {
        pclass_sz->has_pk_gg_at_z_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"pk_gg_at_z_2h") != NULL) ) {
        pclass_sz->has_pk_gg_at_z_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"pk_bb_at_z_1h") != NULL) ) {
        pclass_sz->has_pk_bb_at_z_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"pk_bb_at_z_2h") != NULL) ) {
        pclass_sz->has_pk_bb_at_z_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"pk_b_at_z_2h") != NULL) ) {
        pclass_sz->has_pk_b_at_z_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"pressure_profile_2h") != NULL) ) {
        pclass_sz->has_gas_pressure_profile_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        // printf("ok %d\n",pclass_sz->has_gas_pressure_profile_2h);
        // exit(0);
      }

      if ((strstr(string1,"gas_density_profile_2h") != NULL) ) {
        pclass_sz->has_gas_density_profile_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        // printf("ok %d\n",pclass_sz->has_gas_density_profile_2h);
        // exit(0);
      }


      if ((strstr(string1,"pk_em_at_z_1h") != NULL) ) {
        pclass_sz->has_pk_em_at_z_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"pk_em_at_z_2h") != NULL) ) {
        pclass_sz->has_pk_em_at_z_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"pk_HI_at_z_1h") != NULL) ) {
        pclass_sz->has_pk_HI_at_z_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"pk_HI_at_z_2h") != NULL) ) {
        pclass_sz->has_pk_HI_at_z_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }



      if ((strstr(string1,"pk_at_z_1h") != NULL) ) {
        pclass_sz->has_pk_at_z_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"pk_at_z_2h") != NULL) ) {
        pclass_sz->has_pk_at_z_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"bk_ttg_at_z_1h") != NULL) ) {
        pclass_sz->has_bk_ttg_at_z_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_vrms2 = 1;
      }

      if ((strstr(string1,"bk_ttg_at_z_2h") != NULL) ) {
        pclass_sz->has_bk_ttg_at_z_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_vrms2 = 1;
      }

      if ((strstr(string1,"bk_ttg_at_z_3h") != NULL) ) {
        pclass_sz->has_bk_ttg_at_z_3h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_vrms2 = 1;
      }


      if ((strstr(string1,"bk_at_z_1h") != NULL) ) {
        pclass_sz->has_bk_at_z_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"bk_at_z_2h") != NULL) ) {
        pclass_sz->has_bk_at_z_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"bk_at_z_3h") != NULL) ) {
        pclass_sz->has_bk_at_z_3h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"bk_at_z_hf") != NULL) ) {
        pclass_sz->has_bk_at_z_hf =_TRUE_;
        pclass_sz->has_knl = _TRUE_;
        pclass_sz->has_nl_index = _TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_sigma = 1;
        pclass_sz->has_pk = _TRUE_;
      }

      if ((strstr(string1,"bk_ttg_at_z_hf") != NULL) ) {
        pclass_sz->has_bk_ttg_at_z_hf =_TRUE_;
        pclass_sz->has_mean_galaxy_bias =_TRUE_;
        pclass_sz->has_knl = _TRUE_;
        pclass_sz->has_nl_index = _TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_sigma = 1;
        pclass_sz->has_vrms2 = 1;
        pclass_sz->has_pk = _TRUE_;
      }

      if ((strstr(string1,"vrms2") != NULL) ) {
        pclass_sz->has_vrms2 = _TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
      }

      if ((strstr(string1,"kSZ_kSZ_lensmag_1h") != NULL) ) {
        pclass_sz->has_kSZ_kSZ_lensmag_1halo =_TRUE_;
        pclass_sz->has_vrms2 = _TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"kSZ_kSZ_gal_1h") != NULL) ) {
        pclass_sz->has_kSZ_kSZ_gal_1h =_TRUE_;
        pclass_sz->has_vrms2 = _TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"kSZ_kSZ_gal fft (1h)") != NULL) ) {
        pclass_sz->has_kSZ_kSZ_gal_1h_fft =_TRUE_;
        pclass_sz->has_vrms2 = _TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"kSZ_kSZ_gal_2h") != NULL) ) {
        pclass_sz->has_kSZ_kSZ_gal_2h =_TRUE_;
        pclass_sz->has_vrms2 = _TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"kSZ_kSZ_gal fft (2h)") != NULL) ) {
        pclass_sz->has_kSZ_kSZ_gal_2h_fft =_TRUE_;
        pclass_sz->has_vrms2 = _TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }



      if ((strstr(string1,"kSZ_kSZ_gal_3h") != NULL) ) {
        pclass_sz->has_kSZ_kSZ_gal_3h =_TRUE_;
        pclass_sz->has_vrms2 = _TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }




      if ((strstr(string1,"kSZ_kSZ_gal_lensing_term") != NULL) ) {
        ppt->has_scalars = _TRUE_;
        ppt->has_cl_cmb_temperature = _TRUE_;
        ppt->has_cls = _TRUE_;
        // pclass_sz->has_gal_lens_1h = _TRUE_;
        // pclass_sz->has_gal_lens_2h = _TRUE_;
        pclass_sz->has_kSZ_kSZ_gal_lensing_term =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        ppt->l_scalar_max = 10000;

      }



      if ((strstr(string1,"kSZ_kSZ_gal fft (3h)") != NULL) ) {
        pclass_sz->has_kSZ_kSZ_gal_3h_fft =_TRUE_;
        pclass_sz->has_vrms2 = _TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"kSZ_kSZ_gal_hf") != NULL) ) {
        pclass_sz->has_mean_galaxy_bias = _TRUE_;
        pclass_sz->has_kSZ_kSZ_gal_hf =_TRUE_;
        pclass_sz->has_vrms2 = _TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->has_knl = _TRUE_;
        pclass_sz->has_nl_index = _TRUE_;
        pclass_sz->need_hmf = 1; // need sigma at R,z
        pclass_sz->has_pk = _TRUE_;
        pclass_sz->has_galaxy = _TRUE_;


      }

      if ((strstr(string1,"kSZ_kSZ_lens_hf") != NULL) ) {

        pclass_sz->has_kSZ_kSZ_lens_hf =_TRUE_;
        pclass_sz->has_vrms2 = _TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->has_knl = _TRUE_;
        pclass_sz->has_nl_index = _TRUE_;
        pclass_sz->need_hmf = 1; // need sigma at R,z
        pclass_sz->has_pk = _TRUE_;
        pclass_sz->has_lensing = _TRUE_;
      }

      if ((strstr(string1,"gal_gal_lens fft (1h)") != NULL) ) {
        pclass_sz->has_gal_gal_lens_1h_fft =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }
      if ((strstr(string1,"gal_gal_lens fft (2h)") != NULL) ) {
        pclass_sz->has_gal_gal_lens_2h_fft =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }
      if ((strstr(string1,"gal_gal_lens fft (3h)") != NULL) ) {
        pclass_sz->has_gal_gal_lens_3h_fft =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"kSZ_kSZ_lens fft (1h)") != NULL) ) {
        pclass_sz->has_kSZ_kSZ_lens_1h_fft =_TRUE_;
        pclass_sz->has_vrms2 = _TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"kSZ_kSZ_lens fft (2h)") != NULL) ) {
        pclass_sz->has_kSZ_kSZ_lens_2h_fft =_TRUE_;
        pclass_sz->has_vrms2 = _TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"kSZ_kSZ_lens fft (3h)") != NULL) ) {
        pclass_sz->has_kSZ_kSZ_lens_3h_fft =_TRUE_;
        pclass_sz->has_vrms2 = _TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"kSZ_kSZ_lens_covmat") != NULL) ) {
        ppt->has_scalars = _TRUE_;
        ppt->has_cl_cmb_temperature = _TRUE_;
        ppt->has_cl_cmb_lensing_potential = _TRUE_;
        ppt->has_cls = _TRUE_;
        ple->has_lensed_cls = _TRUE_;
        pclass_sz->has_lens_lens_1h = _TRUE_;
        pclass_sz->has_lens_lens_2h = _TRUE_;
        pclass_sz->has_kSZ_kSZ_lens_covmat =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        ppt->l_scalar_max = 10000;
        pclass_sz->need_ksz_template = 1;
        pclass_sz->need_tt_noise = 1;
        pclass_sz->need_lensing_noise = 1;

      }

      if ((strstr(string1,"kSZ_kSZ_lens_lensing_term") != NULL) ) {
        ppt->has_scalars = _TRUE_;
        ppt->has_cl_cmb_temperature = _TRUE_;
        ppt->has_cls = _TRUE_;
        pclass_sz->has_lens_lens_1h = _TRUE_;
        pclass_sz->has_lens_lens_2h = _TRUE_;
        pclass_sz->has_kSZ_kSZ_lens_lensing_term =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        ppt->l_scalar_max = 10000;
      }


      if ((strstr(string1,"kSZ_kSZ_gallens_hf") != NULL) ) {

        pclass_sz->has_kSZ_kSZ_gallens_hf =_TRUE_;
        pclass_sz->has_vrms2 = _TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->has_knl = _TRUE_;
        pclass_sz->has_nl_index = _TRUE_;
        pclass_sz->need_hmf = 1; // need sigma at R,z
        pclass_sz->has_lensing = _TRUE_;
        pclass_sz->has_pk = _TRUE_;
      }


      if ((strstr(string1,"kSZ_kSZ_gallens fft (1h)") != NULL) ) {
        pclass_sz->has_kSZ_kSZ_gallens_1h_fft =_TRUE_;
        pclass_sz->has_vrms2 = _TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"kSZ_kSZ_gallens fft (2h)") != NULL) ) {
        pclass_sz->has_kSZ_kSZ_gallens_2h_fft =_TRUE_;
        pclass_sz->has_vrms2 = _TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"kSZ_kSZ_gallens fft (3h)") != NULL) ) {
        pclass_sz->has_kSZ_kSZ_gallens_3h_fft =_TRUE_;
        pclass_sz->has_vrms2 = _TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"kSZ_kSZ_gallens_covmat") != NULL) ) {
        ppt->has_scalars = _TRUE_;
        ppt->has_cl_cmb_temperature = _TRUE_;
        ppt->has_cl_cmb_lensing_potential = _TRUE_;
        ppt->has_cls = _TRUE_;
        ple->has_lensed_cls = _TRUE_;
        pclass_sz->has_gallens_gallens_1h = _TRUE_;
        pclass_sz->has_gallens_gallens_2h = _TRUE_;
        pclass_sz->has_kSZ_kSZ_gallens_covmat =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        ppt->l_scalar_max = 10000;
        pclass_sz->need_ksz_template = 1;
        pclass_sz->need_tt_noise = 1;

      }

      if ((strstr(string1,"kSZ_kSZ_gallens_lensing_term") != NULL) ) {
        ppt->has_scalars = _TRUE_;
        ppt->has_cl_cmb_temperature = _TRUE_;
        ppt->has_cls = _TRUE_;
        pclass_sz->has_gallens_lens_1h = _TRUE_;
        pclass_sz->has_gallens_lens_2h = _TRUE_;
        pclass_sz->has_kSZ_kSZ_gallens_lensing_term =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        ppt->l_scalar_max = 10000;
      }



      if ((strstr(string1,"kSZ_kSZ_1h") != NULL) ) {
        pclass_sz->has_kSZ_kSZ_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_vrms2 = _TRUE_;
      }
      if ((strstr(string1,"kSZ_kSZ_2h") != NULL) ) {
        pclass_sz->has_kSZ_kSZ_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_vrms2 = _TRUE_;
      }
      if ((strstr(string1,"kSZ_kSZ_tSZ_1h") != NULL) ) {
        pclass_sz->has_kSZ_kSZ_tSZ_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_vrms2 = _TRUE_;
      }
      if ((strstr(string1,"kSZ_kSZ_tSZ_2h") != NULL) ) {
        pclass_sz->has_kSZ_kSZ_tSZ_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_vrms2 = _TRUE_;
      }

      if ((strstr(string1,"kSZ_kSZ_tSZ_3h") != NULL) ) {
        pclass_sz->has_kSZ_kSZ_tSZ_3h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_vrms2 = _TRUE_;
      }

      if ((strstr(string1,"tSZ_tSZ_tSZ_1h") != NULL) ) {
        pclass_sz->has_tSZ_tSZ_tSZ_1halo =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }
      if ((strstr(string1,"tSZ_tSZ_tSZ_2h") != NULL) ) {
        pclass_sz->has_tSZ_tSZ_tSZ_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }
      if ((strstr(string1,"tSZ_tSZ_tSZ_3h") != NULL) ) {
        pclass_sz->has_tSZ_tSZ_tSZ_3h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"galn_galn_1h") != NULL) ) {
        pclass_sz->has_ngal_ngal_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"galn_galn_2h") != NULL) ) {
        pclass_sz->has_ngal_ngal_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"cib_cib_1h") != NULL) ) {
        pclass_sz->has_cib_cib_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }
      if ((strstr(string1,"cib_monopole") != NULL) ) {
        pclass_sz->has_cib_monopole =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }
      if ((strstr(string1,"cib_shotnoise") != NULL) ) {
        pclass_sz->has_cib_shotnoise =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"dcib0dz") != NULL) ) {
        pclass_sz->has_dcib0dz =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }
      if ((strstr(string1,"dydz") != NULL) ) {
        pclass_sz->has_dydz =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"cib_cib_2h") != NULL) ) {
        pclass_sz->has_cib_cib_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"gal_gal_1h") != NULL) ) {
        pclass_sz->has_gal_gal_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"mean_galaxy_bias") != NULL) ) {
        pclass_sz->has_gal_gal_1h =_TRUE_;
        pclass_sz->has_mean_galaxy_bias =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"gal_gal_2h") != NULL) ) {
        pclass_sz->has_gal_gal_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"gal_gal_hf") != NULL) ) {
        pclass_sz->has_gal_gal_hf =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->has_pk = _TRUE_;
        // pclass_sz->need_sigma = 1;

        // printf("has_pk = %d\n",pclass_sz->has_pk);
        // exit(0);


        // pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"galn_galn_hf") != NULL) ) {
        pclass_sz->has_ngal_ngal_hf =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->has_pk = _TRUE_;

        // class_read_double("effective_galaxy_bias",pclass_sz->effective_galaxy_bias);

      }



        class_read_double("effective_galaxy_bias",pclass_sz->effective_galaxy_bias);
        class_read_double("use_bg_eff_in_ksz2g_eff",pclass_sz->use_bg_eff_in_ksz2g_eff);


      if ((strstr(string1,"kSZ_kSZ_gal_covmat") != NULL) ) {
        ppt->has_scalars = _TRUE_;
        ppt->has_cl_cmb_temperature = _TRUE_;
        ppt->has_cl_cmb_lensing_potential = _TRUE_;
        ppt->has_cls = _TRUE_;
        ple->has_lensed_cls = _TRUE_;
        // pclass_sz->has_gal_gal_1h = _TRUE_;
        // pclass_sz->has_gal_gal_2h = _TRUE_;
        pclass_sz->has_kSZ_kSZ_gal_covmat =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        ppt->l_scalar_max = 10000;
        pclass_sz->need_ksz_template = 1;
        pclass_sz->need_tt_noise = 1;

      }

      if (pclass_sz->has_kSZ_kSZ_gal_covmat ==_TRUE_){
        if ((pclass_sz->has_gal_gal_hf== _FALSE_) && ((pclass_sz->has_gal_gal_1h+pclass_sz->has_gal_gal_2h)==_FALSE_) ){
          printf("you need to request computation of gal gal to get the covmat.\n");
          exit(0);
        }
      }



      if ((strstr(string1,"galn_lens_1h") != NULL) ) {
        pclass_sz->has_ngal_lens_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"galn_lens_2h") != NULL) ) {
        pclass_sz->has_ngal_lens_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"galn_tsz_1h") != NULL) ) {
        pclass_sz->has_ngal_tsz_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"galn_tsz_2h") != NULL) ) {
        pclass_sz->has_ngal_tsz_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"galn_lens_hf") != NULL) ) {
        pclass_sz->has_ngal_lens_hf =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->has_pk = _TRUE_;
        // pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"galn_gallens_1h") != NULL) ) {
        pclass_sz->has_ngal_gallens_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"galn_gallens_2h") != NULL) ) {
        pclass_sz->has_ngal_gallens_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"galn_IA_2h") != NULL) ) {
        pclass_sz->has_ngal_IA_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"lensmagn_gallens_1h") != NULL) ) {
        pclass_sz->has_nlensmag_gallens_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"lensmagn_gallens_2h") != NULL) ) {
        pclass_sz->has_nlensmag_gallens_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"lensmagn_tsz_1h") != NULL) ) {
        pclass_sz->has_nlensmag_tsz_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"lensmagn_tsz_2h") != NULL) ) {
        pclass_sz->has_nlensmag_tsz_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"tau_gal_1h") != NULL) ) {
        pclass_sz->has_tau_gal_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"tau_gal_2h") != NULL) ) {
        pclass_sz->has_tau_gal_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"tau_tau_1h") != NULL) ) {
        pclass_sz->has_tau_tau_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"tau_tau_2h") != NULL) ) {
        pclass_sz->has_tau_tau_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"gal_lens_1h") != NULL) ) {
        pclass_sz->has_gal_lens_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"gal_lens_2h") != NULL) ) {
        pclass_sz->has_gal_lens_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"gal_lens_hf") != NULL) ) {
        pclass_sz->has_gal_lens_hf =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->has_pk = _TRUE_;
        // pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"lens_lens_hf") != NULL) ) {
        pclass_sz->has_lens_lens_hf =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->has_pk = _TRUE_;
        // pclass_sz->need_hmf = 1;
      }



      if ((strstr(string1,"gal_lensmag_hf") != NULL) ) {
        pclass_sz->has_gal_lensmag_hf =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->has_pk = _TRUE_;
        // pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"galn_lensmagn_hf") != NULL) ) {
        pclass_sz->has_ngal_nlensmag_hf =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->has_pk = _TRUE_;
        // pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"lens_lensmag_hf") != NULL) ) {
        pclass_sz->has_lens_lensmag_hf =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->has_pk = _TRUE_;
        // pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"lensmag_lensmag_hf") != NULL) ) {
        pclass_sz->has_lensmag_lensmag_hf =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->has_pk = _TRUE_;
        // pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"tSZ_lensmag_1h") != NULL) ) {
        pclass_sz->has_tSZ_lensmag_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"tSZ_lensmag_2h") != NULL) ) {
        pclass_sz->has_tSZ_lensmag_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"gallens_lensmag_1h") != NULL) ) {
        pclass_sz->has_gallens_lensmag_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"gallens_lensmag_2h") != NULL) ) {
        pclass_sz->has_gallens_lensmag_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }




      if ((strstr(string1,"gal_lensmag_1h") != NULL) ) {
        pclass_sz->has_gal_lensmag_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"gal_lensmag_2h") != NULL) ) {
        pclass_sz->has_gal_lensmag_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"gallens_gallens_1h") != NULL) ) {
        pclass_sz->has_gallens_gallens_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"gallens_gallens_2h") != NULL) ) {
        pclass_sz->has_gallens_gallens_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"gallens_lens_1h") != NULL) ) {
        pclass_sz->has_gallens_lens_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"gallens_lens_2h") != NULL) ) {
        pclass_sz->has_gallens_lens_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"gal_gallens_1h") != NULL) ) {
        pclass_sz->has_gal_gallens_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"gal_gallens_2h") != NULL) ) {
        pclass_sz->has_gal_gallens_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"tSZ_gallens_1h") != NULL) ) {
        pclass_sz->has_tSZ_gallens_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"tSZ_gallens_2h") != NULL) ) {
        pclass_sz->has_tSZ_gallens_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }



      if ((strstr(string1,"gamma_gal_gallens_1h") != NULL) ) {
        pclass_sz->has_gal_gallens_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->convert_cls_to_gamma = 1;
      }

      if ((strstr(string1,"gamma_gal_gallens_2h") != NULL) ) {
        pclass_sz->has_gal_gallens_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->convert_cls_to_gamma = 1;
      }



      if ((strstr(string1,"lensmag_lensmag_1h") != NULL) ) {
        pclass_sz->has_lensmag_lensmag_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"lensmag_lensmag_2h") != NULL) ) {
        pclass_sz->has_lensmag_lensmag_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"lens_lensmag_1h") != NULL) && (strstr(string1,"gallens_lensmag_1h") == NULL)) {
        pclass_sz->has_lens_lensmag_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"lens_lensmag_2h") != NULL) && (strstr(string1,"gallens_lensmag_2h") == NULL)) {
        pclass_sz->has_lens_lensmag_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"tSZ_gal_1h") != NULL) ) {
        pclass_sz->has_tSZ_gal_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"IA_gal_2h") != NULL) ) {
        pclass_sz->has_IA_gal_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"tSZ_gal_2h") != NULL) ) {
        pclass_sz->has_tSZ_gal_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"tSZ_cib_1h") != NULL) ) {
        pclass_sz->has_tSZ_cib_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"tSZ_cib_2h") != NULL) ) {
        pclass_sz->has_tSZ_cib_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"gallens_cib_1h") != NULL) ) {
        pclass_sz->has_gallens_cib_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"gallens_cib_2h") != NULL) ) {
        pclass_sz->has_gallens_cib_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"gal_cib_1h") != NULL) ) {
        pclass_sz->has_gal_cib_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"gal_cib_2h") != NULL) ) {
        pclass_sz->has_gal_cib_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"lens_cib_1h") != NULL) ) {
        pclass_sz->has_lens_cib_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"lens_cib_2h") != NULL) ) {
        pclass_sz->has_lens_cib_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }


      if ((strstr(string1,"lens_lens_1h") != NULL) ) {
        pclass_sz->has_lens_lens_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }
      if ((strstr(string1,"lens_lens_2h") != NULL) ) {
        pclass_sz->has_lens_lens_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"custom1_custom1_1h") != NULL) ) {
        pclass_sz->has_custom1_custom1_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_custom1 = _TRUE_;
      }

      // if ((strstr(string1,"has_b_custom1") != NULL) ) {
      //   pclass_sz->has_b_custom1 = _TRUE_;
      // }

      if ((strstr(string1,"custom1_custom1_2h") != NULL) ) {
        pclass_sz->has_custom1_custom1_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_custom1 = _TRUE_;
      }

      if ((strstr(string1,"custom1_lens_1h") != NULL) ) {
        pclass_sz->has_custom1_lens_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_custom1 = _TRUE_;
        pclass_sz->has_lensing = _TRUE_;
      }

      if ((strstr(string1,"custom1_lens_2h") != NULL) ) {
        pclass_sz->has_custom1_lens_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_custom1 = _TRUE_;
        pclass_sz->has_lensing = _TRUE_;
      }

      if ((strstr(string1,"custom1_tSZ_1h") != NULL) ) {
        pclass_sz->has_custom1_tSZ_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_custom1 = _TRUE_;
        pclass_sz->has_electron_pressure = _TRUE_;
      }

      if ((strstr(string1,"custom1_tSZ_2h") != NULL) ) {
        pclass_sz->has_custom1_tSZ_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_custom1 = _TRUE_;
        pclass_sz->has_electron_pressure = _TRUE_;
      }


      if ((strstr(string1,"custom1_cib_1h") != NULL) ) {
        pclass_sz->has_custom1_cib_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_custom1 = _TRUE_;
        pclass_sz->has_cib = _TRUE_;
      }

      if ((strstr(string1,"custom1_cib_2h") != NULL) ) {
        pclass_sz->has_custom1_cib_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_custom1 = _TRUE_;
        pclass_sz->has_cib = _TRUE_;
      }

      if ((strstr(string1,"custom1_gal_1h") != NULL) ) {
        pclass_sz->has_custom1_gal_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_custom1 = _TRUE_;
        pclass_sz->has_galaxy = _TRUE_;
      }

      if ((strstr(string1,"custom1_gal_2h") != NULL) ) {
        pclass_sz->has_custom1_gal_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_custom1 = _TRUE_;
        pclass_sz->has_galaxy = _TRUE_;
      }

      if ((strstr(string1,"custom1_gallens_1h") != NULL) ) {
        pclass_sz->has_custom1_gallens_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_custom1 = _TRUE_;
        pclass_sz->has_lensing = _TRUE_;
      }

      if ((strstr(string1,"custom1_gallens_2h") != NULL) ) {
        pclass_sz->has_custom1_gallens_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_custom1 = _TRUE_;
        pclass_sz->has_lensing = _TRUE_;
      }

      if ((strstr(string1,"tSZ_lens_1h") != NULL) ) {
        pclass_sz->has_tSZ_lens_1h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"tSZ_lens_2h") != NULL) ) {
        pclass_sz->has_tSZ_lens_2h =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

      if ((strstr(string1,"isw_lens") != NULL) ) {
        pclass_sz->has_isw_lens =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }
      if ((strstr(string1,"isw_tsz") != NULL) ) {
        pclass_sz->has_isw_tsz =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }
      if ((strstr(string1,"isw_auto") != NULL) ) {
        pclass_sz->has_isw_auto =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }
      if ((strstr(string1,"sz_cluster_counts") != NULL) ) {
        // printf("counts\n");
        pclass_sz->has_sz_counts =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_500c = 1;

      }

      if ((strstr(string1,"sz_cluster_counts_fft") != NULL) ) {
        // printf("counts\n");
        pclass_sz->has_sz_counts =_TRUE_;
        pclass_sz->has_sz_counts_fft =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_500c = 1;

      }


      if ((strstr(string1,"sz_unbinned_cluster_counts") != NULL) ) {
        // printf("counts\n");
        pclass_sz->has_sz_counts =_TRUE_;
        pclass_sz->has_sz_rates =_TRUE_;
        pclass_sz->has_hmf =_TRUE_;
        ppt->has_density_transfers=_TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
        pclass_sz->has_500c = 1;
        // pclass_sz->has_completeness_for_ps_SZ = 1;

      }

      if (strstr(string1, "dndM") || strstr(string1, "dndm") || strstr(string1, "dndlnM") || strstr(string1, "dndlnm")) {
        pclass_sz->has_dndlnM = _TRUE_;
        ppt->has_density_transfers = _TRUE_;
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->need_hmf = 1;
      }

    if (pclass_sz->need_hmf == 1)
        pclass_sz->need_sigma = 1;


      if ((strstr(string1,"m200c_to_m500c") != NULL) ) {
          pclass_sz->need_m200c_to_m500c = 1;
          ppt->has_density_transfers=_TRUE_;
          ppt->has_pk_matter = _TRUE_;
          ppt->has_perturbations = _TRUE_;
          pnl->has_pk_cb = _TRUE_;
          pnl->has_pk_m = _TRUE_;
          pclass_sz->need_sigma = 1;
          }
      if ((strstr(string1,"m500c_to_m200c") != NULL) ) {
          pclass_sz->need_m500c_to_m200c = 1;
          ppt->has_density_transfers=_TRUE_;
          ppt->has_pk_matter = _TRUE_;
          ppt->has_perturbations = _TRUE_;
          pnl->has_pk_cb = _TRUE_;
          pnl->has_pk_m = _TRUE_;
          pclass_sz->need_sigma = 1;
          }
      if ((strstr(string1,"m200m_to_m500c") != NULL) ) {
          pclass_sz->need_m200m_to_m500c = 1;
          ppt->has_density_transfers=_TRUE_;
          ppt->has_pk_matter = _TRUE_;
          ppt->has_perturbations = _TRUE_;
          pnl->has_pk_cb = _TRUE_;
          pnl->has_pk_m = _TRUE_;
          pclass_sz->need_sigma = 1;
          }
      if ((strstr(string1,"m200m_to_m200c") != NULL) ) {
          pclass_sz->need_m200m_to_m200c = 1;
          ppt->has_density_transfers=_TRUE_;
          ppt->has_pk_matter = _TRUE_;
          ppt->has_perturbations = _TRUE_;
          pnl->has_pk_cb = _TRUE_;
          pnl->has_pk_m = _TRUE_;
          pclass_sz->need_sigma = 1;
          }
      if ((strstr(string1,"m200c_to_m200m") != NULL) ) {
          pclass_sz->need_m200c_to_m200m = 1;
          ppt->has_density_transfers=_TRUE_;
          ppt->has_pk_matter = _TRUE_;
          ppt->has_perturbations = _TRUE_;
          pnl->has_pk_cb = _TRUE_;
          pnl->has_pk_m = _TRUE_;
          pclass_sz->need_sigma = 1;
          }
      if ((strstr(string1,"tabulate_rhob_xout_at_m_and_z") != NULL) ) {
          pclass_sz->tabulate_rhob_xout_at_m_and_z = 1;
          // ppt->has_density_transfers=_TRUE_;
          // ppt->has_pk_matter = _TRUE_;
          // ppt->has_perturbations = _TRUE_;
          // pnl->has_pk_cb = _TRUE_;
          // pnl->has_pk_m = _TRUE_;
          // pclass_sz->need_sigma = 1;
          }
      if ((strstr(string1,"scale_dependent_bias") != NULL) ) {
        ppt->has_density_transfers=_TRUE_;
        ppt->has_perturbations = _TRUE_;
        pclass_sz->need_ng_bias = 1;
      }

      if ((strstr(string1,"n5k") != NULL) ) {
        ppt->has_pk_matter = _TRUE_;
        ppt->has_perturbations = _TRUE_;
        pnl->has_pk_cb = _TRUE_;
        pnl->has_pk_m = _TRUE_;
        pclass_sz->has_n5k = _TRUE_;
      }

      class_call(parser_read_string(pfc,"include_ssc",&string1,&flag1,errmsg),
                   errmsg,
                   errmsg);
        if (flag1 == _TRUE_) {
            if ((strstr(string1,"yes") != NULL)){
            // pclass_sz->has_sz_cov_N_N_hsv = _TRUE_;
            // pclass_sz->has_sz_cov_Y_N_next_order = _TRUE_;
            pclass_sz->has_sz_cov_Y_Y_ssc = _TRUE_;
            pclass_sz->need_hmf = 1;
        }
      }

      /* multipoles SZ */
      class_call(parser_read_string(pfc,"multipoles_sz",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
     if (flag1 == _TRUE_) {
        if ((strstr(string1,"DKS") != NULL))
          pclass_sz->ell_sz=0;
        else  if ((strstr(string1,"P15") != NULL))
          pclass_sz->ell_sz=1;
        else  if ((strstr(string1,"KS02") != NULL))
          pclass_sz->ell_sz=2;
        else  if ((strstr(string1,"low-ell") != NULL))
          pclass_sz->ell_sz=3;
        else  if ((strstr(string1,"ell_mock") != NULL))
          pclass_sz->ell_sz=4;
        else  if ((strstr(string1,"Pl15_no_low_ell") != NULL))
          pclass_sz->ell_sz=5;
     }






      //units for tSZ spectrum
      class_call(parser_read_string(pfc,"units for tSZ spectrum",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      if (flag1 == _TRUE_) {
        if ((strstr(string1,"muK2") != NULL))
          pclass_sz->exponent_unit=0.;
        else  if ((strstr(string1,"dimensionless") != NULL))
          pclass_sz->exponent_unit=2.;
      }

      //frequency for tSZ spectrum
     //class_call(parser_read_string(pfc,"Frequency for y-distortion",&string1,&flag1,errmsg),
     class_read_double("Frequency for y-distortion in GHz",pclass_sz->nu_y_dist_GHz);

     class_read_int("Frequency_id nu for cib in GHz (to save in file)",pclass_sz->id_nu_cib_to_save);
     class_read_int("Frequency_id nu^prime for cib in GHz (to save in file)",pclass_sz->id_nu_prime_cib_to_save);


    class_read_int("cib_frequency_list_num",pclass_sz->cib_frequency_list_num);
    if (pclass_sz->has_cib_cib_1h
      + pclass_sz->has_cib_cib_2h
      + pclass_sz->has_cib_shotnoise
      + pclass_sz->has_tSZ_cib_1h
      + pclass_sz->has_tSZ_cib_2h
      + pclass_sz->has_gallens_cib_1h
      + pclass_sz->has_gallens_cib_2h
      + pclass_sz->has_gal_cib_1h
      + pclass_sz->has_gal_cib_2h
      + pclass_sz->has_lens_cib_1h
      + pclass_sz->has_lens_cib_2h
      != _FALSE_){
    class_read_list_of_doubles("cib_frequency_list_in_GHz",pclass_sz->cib_frequency_list,pclass_sz->cib_frequency_list_num);
  }


    class_read_int("has_cib_flux_cut",pclass_sz->has_cib_flux_cut);
    if (pclass_sz->has_cib_flux_cut != _FALSE_){
      class_read_list_of_doubles("cib_Snu_cutoff_list [mJy]",pclass_sz->cib_Snu_cutoff_list_in_mJy,pclass_sz->cib_frequency_list_num);
    }

    class_read_int("galaxy_samples_list_num",pclass_sz->galaxy_samples_list_num);
    if ((pclass_sz->has_ngal_ngal_1h
      +  pclass_sz->has_ngal_ngal_2h
      +  pclass_sz->has_ngal_ngal_hf
      +  pclass_sz->has_ngal_lens_1h
      +  pclass_sz->has_ngal_lens_2h
      +  pclass_sz->has_ngal_tsz_1h
      +  pclass_sz->has_ngal_tsz_2h
      +  pclass_sz->has_ngal_lens_hf
      +  pclass_sz->has_ngal_nlensmag_hf
      +  pclass_sz->has_ngal_gallens_1h
      +  pclass_sz->has_ngal_gallens_2h
      +  pclass_sz->has_ngal_IA_2h
      +  pclass_sz->has_nlensmag_gallens_1h
      +  pclass_sz->has_nlensmag_gallens_2h
      +  pclass_sz->has_nlensmag_tsz_1h
      +  pclass_sz->has_nlensmag_tsz_2h
    )
      != _FALSE_){
    class_read_list_of_integers("galaxy_samples_list",pclass_sz->galaxy_samples_list,pclass_sz->galaxy_samples_list_num);

    if (pclass_sz->has_ngal_ngal_hf
       +pclass_sz->has_ngal_lens_hf
       +pclass_sz->has_ngal_nlensmag_hf){
    class_alloc(pclass_sz->effective_galaxy_bias_ngal,sizeof(double *)*pclass_sz->galaxy_samples_list_num,pclass_sz->error_message);
    class_alloc(pclass_sz->effective_galaxy_bias_nl_ngal,sizeof(double *)*pclass_sz->galaxy_samples_list_num,pclass_sz->error_message);


    int index_g;
    // printf("reading effective bias.\n");

    for (index_g = 0;index_g<pclass_sz->galaxy_samples_list_num;index_g++){
      // printf("reading effective bias %d.\n",index_g);
      char input_param_name[_ARGUMENT_LENGTH_MAX_];
      sprintf(input_param_name,"%s%d","effective_galaxy_bias_ngal_",index_g);
      class_read_double(input_param_name,pclass_sz->effective_galaxy_bias_ngal[index_g]);
      sprintf(input_param_name,"%s%d","effective_galaxy_bias_nl_ngal_",index_g);
      class_read_double(input_param_name,pclass_sz->effective_galaxy_bias_nl_ngal[index_g]);
      // printf("reading effective bias %.3e.\n",pclass_sz->effective_galaxy_bias_ngal[index_g]);
      // if (pclass_sz->sz_verbose>1)
        // printf("reading effective bias b[%d] = %.3e\n",index_g,pclass_sz->effective_galaxy_bias_ngal[index_g]);

        }
    // exit(0);
    }

    if (pclass_sz->has_ngal_ngal_1h
      + pclass_sz->has_ngal_ngal_2h
      + pclass_sz->has_ngal_lens_1h
      + pclass_sz->has_ngal_lens_2h
      + pclass_sz->has_ngal_tsz_1h
      + pclass_sz->has_ngal_tsz_2h
      + pclass_sz->has_ngal_gallens_1h
      + pclass_sz->has_ngal_gallens_2h
      + pclass_sz->has_ngal_IA_2h

    ){

    class_alloc(pclass_sz->sigma_log10M_HOD_ngal,sizeof(double *)*pclass_sz->galaxy_samples_list_num,pclass_sz->error_message);
    class_alloc(pclass_sz->alpha_s_HOD_ngal,sizeof(double *)*pclass_sz->galaxy_samples_list_num,pclass_sz->error_message);
    class_alloc(pclass_sz->M1_prime_HOD_ngal,sizeof(double *)*pclass_sz->galaxy_samples_list_num,pclass_sz->error_message);
    class_alloc(pclass_sz->M_min_HOD_ngal,sizeof(double *)*pclass_sz->galaxy_samples_list_num,pclass_sz->error_message);
    class_alloc(pclass_sz->M0_HOD_ngal,sizeof(double *)*pclass_sz->galaxy_samples_list_num,pclass_sz->error_message);
    class_alloc(pclass_sz->x_out_truncated_nfw_profile_satellite_galaxies_ngal,sizeof(double *)*pclass_sz->galaxy_samples_list_num,pclass_sz->error_message);
    class_alloc(pclass_sz->f_cen_HOD_ngal,sizeof(double *)*pclass_sz->galaxy_samples_list_num,pclass_sz->error_message);
    class_alloc(pclass_sz->centrals_only_ngal,sizeof(double *)*pclass_sz->galaxy_samples_list_num,pclass_sz->error_message);
    class_alloc(pclass_sz->satellites_only_ngal,sizeof(double *)*pclass_sz->galaxy_samples_list_num,pclass_sz->error_message);
    class_alloc(pclass_sz->photo_z_params_ngal,sizeof(double *)*pclass_sz->galaxy_samples_list_num,pclass_sz->error_message);
    class_alloc(pclass_sz->dndz_shift_ngal,sizeof(double *)*pclass_sz->galaxy_samples_list_num,pclass_sz->error_message);
    class_alloc(pclass_sz->dndz_stretch_ngal,sizeof(double *)*pclass_sz->galaxy_samples_list_num,pclass_sz->error_message);

    int index_g;
    for (index_g = 0;index_g<pclass_sz->galaxy_samples_list_num;index_g++){
        char input_param_name[_ARGUMENT_LENGTH_MAX_];
        pclass_sz->sigma_log10M_HOD_ngal[index_g] = 1.;
        pclass_sz->alpha_s_HOD_ngal[index_g] = 1.;
        pclass_sz->M1_prime_HOD_ngal[index_g] = 1.;
        pclass_sz->M_min_HOD_ngal[index_g] = 1e11;
        pclass_sz->M0_HOD_ngal[index_g] = 1e11;
        pclass_sz->x_out_truncated_nfw_profile_satellite_galaxies_ngal[index_g] = 1.;
        pclass_sz->f_cen_HOD_ngal[index_g] = 1.;
        pclass_sz->centrals_only_ngal[index_g] = 0.;
        pclass_sz->satellites_only_ngal[index_g] = 0.;
        pclass_sz->photo_z_params_ngal[index_g] = 0.;
        pclass_sz->dndz_shift_ngal[index_g] = 0.;
        pclass_sz->dndz_stretch_ngal[index_g] = 1.;

        sprintf(input_param_name,"%s%d","sigma_log10M_HOD_ngal_",index_g);
        class_read_double(input_param_name,pclass_sz->sigma_log10M_HOD_ngal[index_g]);
        sprintf(input_param_name,"%s%d","alpha_s_HOD_ngal_",index_g);
        class_read_double(input_param_name,pclass_sz->alpha_s_HOD_ngal[index_g]);
        sprintf(input_param_name,"%s%d","M1_prime_HOD_ngal_",index_g);
        class_read_double(input_param_name,pclass_sz->M1_prime_HOD_ngal[index_g]);
        sprintf(input_param_name,"%s%d","M_min_HOD_ngal_",index_g);
        class_read_double(input_param_name,pclass_sz->M_min_HOD_ngal[index_g]);
        sprintf(input_param_name,"%s%d","M0_HOD_ngal_",index_g);
        class_read_double(input_param_name,pclass_sz->M0_HOD_ngal[index_g]);
        sprintf(input_param_name,"%s%d","x_out_truncated_nfw_profile_satellite_galaxies_ngal_",index_g);
        class_read_double(input_param_name,pclass_sz->x_out_truncated_nfw_profile_satellite_galaxies_ngal[index_g]);
        sprintf(input_param_name,"%s%d","f_cen_HOD_ngal_",index_g);
        class_read_double(input_param_name,pclass_sz->f_cen_HOD_ngal[index_g]);
        sprintf(input_param_name,"%s%d","centrals_only_ngal_",index_g);
        class_read_double(input_param_name,pclass_sz->centrals_only_ngal[index_g]);
        sprintf(input_param_name,"%s%d","satellites_only_ngal_",index_g);
        class_read_double(input_param_name,pclass_sz->satellites_only_ngal[index_g]);
        sprintf(input_param_name,"%s%d","photo_z_params_ngal_",index_g);
        class_read_double(input_param_name,pclass_sz->photo_z_params_ngal[index_g]);
        sprintf(input_param_name,"%s%d","dndz_shift_ngal_",index_g);
        class_read_double(input_param_name,pclass_sz->dndz_shift_ngal[index_g]);
        sprintf(input_param_name,"%s%d","dndz_stretch_ngal_",index_g);
        class_read_double(input_param_name,pclass_sz->dndz_stretch_ngal[index_g]);
          }
      }

  }


     // Table 1  of MM20
     class_read_double("Redshift_evolution_of_dust_temperature",pclass_sz->alpha_cib);  //redshift evolution of dust temperature
     class_read_double("Dust_temperature_today_in_Kelvins",pclass_sz->T0_cib);  // dust temperature today in Kelvins
     class_read_double("Emissivity_index_of_sed",pclass_sz->beta_cib); // emissivity index of sed
     class_read_double("Power_law_index_of_SED_at_high_frequency",pclass_sz->gamma_cib); // Power law index of SED at high frequency
     class_read_double("Redshift_evolution_of_L_-_M_normalisation",pclass_sz->delta_cib); // Redshift evolution of L − M normalisation
     class_read_double("Most_efficient_halo_mass_in_Msun",pclass_sz->m_eff_cib); // Most efficient halo mass in Msun
     class_read_double("Normalisation_of_L_-_M_relation_in_[JyMPc2/Msun]",pclass_sz->L0_cib); // Normalisation of L − M relation in [Jy MPc2/Msun]
     class_read_double("Size_of_halo_masses_sourcing_CIB_emission",pclass_sz->sigma2_LM_cib); // Size of of halo masses sourcing CIB emission
     class_read_double("z_obs_(CIB)",pclass_sz->z_obs_cib);
     class_read_double("z_plateau_cib",pclass_sz->z_plateau_cib);
     class_read_double("M_min_subhalo_in_Msun",pclass_sz->M_min_subhalo_in_Msun);

     class_read_double("maniyar_cib_etamax",pclass_sz->maniyar_cib_etamax);
     class_read_double("maniyar_cib_zc",pclass_sz->maniyar_cib_zc);
     class_read_double("maniyar_cib_tau",pclass_sz->maniyar_cib_tau);
     class_read_double("maniyar_cib_fsub",pclass_sz->maniyar_cib_fsub);

     class_read_int("use_redshift_dependent_M_min",pclass_sz->use_redshift_dependent_M_min);
     class_read_int("use_nc_1_for_all_halos_cib_HOD",pclass_sz->use_nc_1_for_all_halos_cib_HOD);

     class_read_int("cib_nu0_norm",pclass_sz->cib_nu0_norm);
     class_read_int("use scale dependent bias (from non Gaussianity)",pclass_sz->has_ng_in_bh);
     class_read_double("fNL",pclass_sz->fNL);
     class_read_int("cosmo_model",pclass_sz->cosmo_model);


     class_read_int("use_Amod",pclass_sz->use_Amod);
     class_read_double("Amod",pclass_sz->Amod);

     class_read_int("use_nlbias",pclass_sz->use_nl_bias);
     class_read_double("bnl",pclass_sz->bnl);

     // FMcC edit: read in p_fNL parameter
     // class_read_double("p_fNL",pclass_sz->p_fNL);
     // end FMcC edit

     if(pclass_sz->has_ng_in_bh){
       ppt->has_density_transfers=_TRUE_;
       ppt->has_perturbations = _TRUE_;
       pclass_sz->need_ng_bias = 1;
     }


      /* concentration parameter SZ */
      class_call(parser_read_string(pfc,"concentration parameter",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      if (flag1 == _FALSE_) {
        class_call(parser_read_string(pfc,"concentration_parameter",&string1,&flag1,errmsg),errmsg,errmsg);
      } 
     if (flag1 == _TRUE_) {
        if ((strstr(string1,"D08") != NULL))
          pclass_sz->concentration_parameter=0;
        else  if ((strstr(string1,"S00") != NULL))
          pclass_sz->concentration_parameter=1;
        else  if ((strstr(string1,"K10") != NULL))
          pclass_sz->concentration_parameter=2;
        else  if ((strstr(string1,"SC14") != NULL))
          pclass_sz->concentration_parameter=3;
        else  if ((strstr(string1,"Z09") != NULL))
          pclass_sz->concentration_parameter=4;
        else  if ((strstr(string1,"DM14") != NULL))
          pclass_sz->concentration_parameter=5; // Dutton and Maccio 2014 (https://arxiv.org/pdf/1402.7073.pdf)
        else  if ((strstr(string1,"B13") != NULL))
          pclass_sz->concentration_parameter=6; // https://arxiv.org/pdf/1112.5479.pdf
        else  if ((strstr(string1,"fixed") != NULL))
          pclass_sz->concentration_parameter=7; // https://arxiv.org/pdf/1112.5479.pdf
          // in this case the concentration is read (see class_sz_precisions.h)
          }






      /* pressure profile for thermal SZ effect */
      class_call(parser_read_string(pfc,"pressure profile",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      if (flag1 == _FALSE_) {
        class_call(parser_read_string(pfc,"pressure_profile",&string1,&flag1,errmsg),errmsg,errmsg);
      } 


      if (flag1 == _TRUE_) {

        if ((strstr(string1,"P13") != NULL)
         || (strstr(string1,"Planck2013") != NULL)
         || (strstr(string1,"Planck13") != NULL))
        
          pclass_sz->pressure_profile=0;
        
        else  if ((strstr(string1,"A10") != NULL) 
               || (strstr(string1,"Arnaud10") != NULL) 
               || (strstr(string1,"Arnaud2010") != NULL))

          pclass_sz->pressure_profile=2;
        
        // else  if ((strstr(string1,"A10XRAY") != NULL))
        
        //   pclass_sz->pressure_profile=5;

        else  if ((strstr(string1,"Custom. GNFW") != NULL)
               || (strstr(string1,"custom_gnfw") != NULL)
               || (strstr(string1,"custom_GNFW") != NULL)
               || (strstr(string1,"Custom_GNFW") != NULL))

          pclass_sz->pressure_profile=3;
        
        else  if ((strstr(string1,"B12") != NULL)
               || (strstr(string1,"Battaglia12") != NULL)
               || (strstr(string1,"Battaglia2012") != NULL))

          pclass_sz->pressure_profile=4;

        }




      /* integration method pressure profile SZ */
      class_call(parser_read_string(pfc,"integration method (pressure profile)",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
     if (flag1 == _TRUE_) {
        if ((strstr(string1,"patterson") != NULL))
          pclass_sz->integration_method_pressure_profile=0;
        else  if ((strstr(string1,"gsl") != NULL))
          pclass_sz->integration_method_pressure_profile=1;
       else  if ((strstr(string1,"spline") != NULL))
          pclass_sz->integration_method_pressure_profile=2;
          }

      class_read_double("pressure_profile_epsrel",pclass_sz->pressure_profile_epsrel);
      class_read_double("pressure_profile_epsabs",pclass_sz->pressure_profile_epsabs);

      class_read_double("nfw_profile_epsrel",pclass_sz->nfw_profile_epsrel);
      class_read_double("nfw_profile_epsabs",pclass_sz->nfw_profile_epsabs);


      class_read_double("fEDE",pclass_sz->fEDE);
      class_read_double("log10z_c",pclass_sz->log10z_c);
      class_read_double("thetai_scf",pclass_sz->thetai_scf);

      class_read_double("M_min_HOD",pclass_sz->M_min_HOD);
      class_read_double("M_min_HOD_cib",pclass_sz->M_min_HOD_cib);
      class_read_double("f_cen_HOD",pclass_sz->f_cen_HOD);
      class_read_double("Delta_z_lens",pclass_sz->Delta_z_lens);
      class_read_double("Delta_z_source",pclass_sz->Delta_z_source);

      class_read_double("bispectrum_lambda_2",pclass_sz->bispectrum_lambda_k2);
      class_read_double("bispectrum_lambda_3",pclass_sz->bispectrum_lambda_k3);


      class_call(parser_read_string(pfc,"M0 equal M_min (HOD)",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      if (flag1 == _FALSE_) {
        class_call(parser_read_string(pfc, "M0_equal_M_min_HOD", &string1, &flag1, errmsg), //new format without spaces
                  errmsg,
                  errmsg);
      }

      if (flag1 == _TRUE_) {
        if ((strstr(string1,"no") != NULL))
          pclass_sz->M0_Mmin_flag=0;
        else  if ((strstr(string1,"yes") != NULL))
          pclass_sz->M0_Mmin_flag=1;
              }
      if (pclass_sz->M0_Mmin_flag == 1){
      pclass_sz->M0_HOD = pclass_sz->M_min_HOD;
      }
      else {
      class_read_double("M0_HOD",pclass_sz->M0_HOD);
      }


      class_read_double("M1_prime_HOD",pclass_sz->M1_prime_HOD);


      class_read_double("alpha_s_HOD",pclass_sz->alpha_s_HOD);
      class_read_double("sigma_log10M_HOD",pclass_sz->sigma_log10M_HOD);
      class_read_double("rho_y_gal",pclass_sz->rho_y_gal);


      class_read_double("x_out_truncated_density_profile (electrons)",pclass_sz->x_out_truncated_nfw_profile_electrons);
      // class_read_double("x_out_truncated_density_profile",pclass_sz->x_out_truncated_density_profile);

      class_read_double("x_out_truncated_nfw_profile",pclass_sz->x_out_truncated_nfw_profile);
      class_read_double("x_out_truncated_nfw_profile_satellite_galaxies",pclass_sz->x_out_truncated_nfw_profile_satellite_galaxies);

      // printf("xout sat = %.5e\n",pclass_sz->x_out_truncated_nfw_profile_satellite_galaxies);
      // printf("xout nfw = %.5e\n",pclass_sz->x_out_truncated_nfw_profile);
      // exit(0);

      class_read_double("cvir_tau_profile_factor",pclass_sz->cvir_tau_profile_factor);
      // class_read_double("x_out_nfw_profile",pclass_sz->x_out_nfw_profile);

      class_read_int("pk_nonlinear_for_vrms2",pclass_sz->pk_nonlinear_for_vrms2);

      class_read_int("hm_consistency",pclass_sz->hm_consistency);
      pclass_sz->hm_consistency_ngbar = pclass_sz->hm_consistency; // by default set same for ng bar as for everything else
      class_read_int("hm_consistency_ngbar",pclass_sz->hm_consistency_ngbar); // read if ng bar has to be treated differently


      class_read_int("T10_alpha_fixed",pclass_sz->T10_alpha_fixed);


      class_read_int("check_consistency_conditions",pclass_sz->check_consistency_conditions);


      /* Noise for covmat (y,y)*/
      class_call(parser_read_string(pfc,"include noise curve in cov(y,y)",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      if (flag1 == _TRUE_) {
        if ((strstr(string1,"no") != NULL))
          pclass_sz->include_noise_cov_y_y=0;
        else  if ((strstr(string1,"yes") != NULL)){
          pclass_sz->include_noise_cov_y_y=1;
          class_read_string("full path to noise curve for yxy",pclass_sz->full_path_to_noise_curve_for_y_y);
          class_read_int("nl_yy_is_binned",pclass_sz->nl_yy_is_binned);
        }

      }
      // class_read_string("full_path_to_noise_curve_for_t_t",pclass_sz->full_path_to_noise_curve_for_t_t);
        // printf("-> File Name: %s\n",pclass_sz->full_path_to_noise_curve_for_t_t);

      // class_read_int("nl_yy_is_binned",pclass_sz->nl_yy_is_binned);
      // if (pclass_sz->sz_verbose >= 1)



      /* temperature mass relation SZ */
      class_call(parser_read_string(pfc,"temperature mass relation",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      if (flag1 == _TRUE_) {
        if ((strstr(string1,"standard") != NULL))
          pclass_sz->temperature_mass_relation=0;
        else  if ((strstr(string1,"lee et al 2019") != NULL))
          pclass_sz->temperature_mass_relation=1;
      }


        class_call(parser_read_string(pfc,"create reference trispectrum for likelihood code",&string1,&flag1,errmsg),
                   errmsg,
                   errmsg);
        if (flag1 == _TRUE_) {
            if ((strstr(string1,"YES") != NULL)){
                pclass_sz->create_ref_trispectrum_for_cobaya=1;
                pclass_sz->has_sz_ps = _TRUE_;
                pclass_sz->has_sz_trispec =_TRUE_;
                ppt->has_density_transfers=_TRUE_;
                ppt->has_pk_matter = _TRUE_;
                ppt->has_perturbations = _TRUE_;
                pnl->has_pk_cb = _TRUE_;
                pnl->has_pk_m = _TRUE_;
              }

            else
                pclass_sz->create_ref_trispectrum_for_cobaya=0;
        }


    /* HMF prescription for massive neutrinos SZ */
      class_call(parser_read_string(pfc,"HMF_prescription_NCDM",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
     if (flag1 == _TRUE_) {
        if ((strstr(string1,"Matter") != NULL))
          pclass_sz->HMF_prescription_NCDM=0;
        else  if ((strstr(string1,"CDM") != NULL))
          pclass_sz->HMF_prescription_NCDM=1;
        else  if ((strstr(string1,"No-pres") != NULL))
          pclass_sz->HMF_prescription_NCDM=2;

            }


//
// class_call(parser_read_string(pfc,"sigma_derivative",&string1,&flag1,errmsg),
//            errmsg,
//            errmsg);
// if (flag1 == _TRUE_) {
//   if ((strstr(string1,"Matter") != NULL))
//     pclass_sz->HMF_prescription_NCDM=0;
//   else  if ((strstr(string1,"CDM") != NULL))
//     pclass_sz->HMF_prescription_NCDM=1;
//   else  if ((strstr(string1,"No-pres") != NULL))
//     pclass_sz->HMF_prescription_NCDM=2;
//
//       }


class_read_double("A_IA",pclass_sz->A_IA);
class_read_double("eta_IA",pclass_sz->eta_IA);
class_read_double("C1_IA",pclass_sz->C1_IA);



    if (pclass_sz->pressure_profile==0){
    //Planck pressure profile 2013, Table 1 of Planck  [arXiv:1207.4061].
     pclass_sz->P0GNFW = 6.41;
     pclass_sz->c500 = 1.81;
     pclass_sz->gammaGNFW = 0.31;
     pclass_sz->alphaGNFW = 1.33;
     pclass_sz->betaGNFW = 4.13;
    }

    else if (pclass_sz->pressure_profile==2){
    //Arnaud et al pressure profile 2010, Eq B2 of [arXiv:0910.1234].
     pclass_sz->P0GNFW = 8.130;
     pclass_sz->c500 = 1.156;
     pclass_sz->gammaGNFW = 0.3292;
     pclass_sz->alphaGNFW = 1.0620;
     pclass_sz->betaGNFW = 5.4807;
    }

    // else if (pclass_sz->pressure_profile==5){
    // //Arnaud et al pressure profile 2010, Eq 12 of https://www.aanda.org/articles/aa/pdf/2010/09/aa13416-09.pdf
    //  pclass_sz->P0GNFW = 8.403;
    //  pclass_sz->c500 = 1.177;
    //  pclass_sz->gammaGNFW = 0.3081;
    //  pclass_sz->alphaGNFW = 1.0510;
    //  pclass_sz->betaGNFW = 5.4905;
    // }


     else if (pclass_sz->pressure_profile==3){
       //Custom GNFW pressure profile
       class_read_double("P0GNFW",pclass_sz->P0GNFW);
       class_read_double("c500",pclass_sz->c500);
       class_read_double("gammaGNFW",pclass_sz->gammaGNFW);
       class_read_double("alphaGNFW",pclass_sz->alphaGNFW);
       class_read_double("betaGNFW",pclass_sz->betaGNFW);
     }





     else if (pclass_sz->pressure_profile==4){
       //Battaglia et al pressure profile [arXiv:1109.3711]


       class_read_double("P0_B12",pclass_sz->P0_B12);
       class_read_double("beta_B12",pclass_sz->beta_B12);
       class_read_double("alpha_B12",pclass_sz->alpha_B12);
       class_read_double("gamma_B12",pclass_sz->gamma_B12);

       // pclass_sz->P0_B12 = 18.1;
       class_read_double("xc_B12",pclass_sz->xc_B12);// = 0.497;
       // pclass_sz->beta_B12 = 4.35;

       class_read_double("alpha_m_P0_B12",pclass_sz->alpha_m_P0_B12);// = 0.154;
       class_read_double("alpha_m_xc_B12",pclass_sz->alpha_m_xc_B12);// = -0.00865;
       class_read_double("alpha_m_beta_B12",pclass_sz->alpha_m_beta_B12);// = 0.0393;

       class_read_double("alpha_z_P0_B12",pclass_sz->alpha_z_P0_B12);// = -0.758;
       class_read_double("alpha_z_xc_B12",pclass_sz->alpha_z_xc_B12);// = 0.731;
       class_read_double("alpha_z_beta_B12",pclass_sz->alpha_z_beta_B12);// = 0.415;

       class_read_double("cp_B12",pclass_sz->c_B12);// = 1.e14;
       class_read_double("cp_B16",pclass_sz->c_B16);// = 1.e14;

       class_read_double("mcut_B12",pclass_sz->mcut_B12);// = 1.e14;
       class_read_double("alphap_m_P0_B12",pclass_sz->alphap_m_P0_B12);// = 0.154;
       class_read_double("alphap_m_xc_B12",pclass_sz->alphap_m_xc_B12);// = -0.00865;
       class_read_double("alphap_m_beta_B12",pclass_sz->alphap_m_beta_B12);// = 0.0393;

       class_read_double("alpha_c_P0_B12",pclass_sz->alpha_c_P0_B12);// = 0.;
       class_read_double("alpha_c_xc_B12",pclass_sz->alpha_c_xc_B12);// = 0.;
       class_read_double("alpha_c_beta_B12",pclass_sz->alpha_c_beta_B12);// = 0.;
       class_read_int("truncate_gas_pressure_wrt_rvir",pclass_sz->truncate_gas_pressure_wrt_rvir);
      if (pclass_sz->truncate_gas_pressure_wrt_rvir){
        pclass_sz->has_vir = 1;
      }
     }



class_read_int("use_websky_m200m_to_m200c_conversion",pclass_sz->use_websky_m200m_to_m200c_conversion);
class_read_int("no_tt_noise_in_kSZ2X_cov",pclass_sz->no_tt_noise_in_kSZ2X_cov);
      /* mass function */
      class_call(parser_read_string(pfc,"sub_halo_mass_function",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
     if (flag1 == _TRUE_) {
        if ((strstr(string1,"TW10") != NULL)){
        pclass_sz->SHMF=1;
        }
        else  if ((strstr(string1,"JvdB14") != NULL)){
        pclass_sz->SHMF=2;
        }
      }

      /* mass function SZ */
      class_call(parser_read_string(pfc,"mass function",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      if (flag1 == _FALSE_) {
        class_call(parser_read_string(pfc,"mass_function",&string1,&flag1,errmsg),errmsg,errmsg);
      }
     if (flag1 == _TRUE_) {
        if ((strstr(string1,"T10M200m") != NULL)){
            pclass_sz->MF=1;
            pclass_sz->integrate_wrt_m200m = 1;
        }

        else  if ((strstr(string1,"B15") != NULL)){
          pclass_sz->MF=2;
          pclass_sz->integrate_wrt_m200m = 1;
        }
        else  if ((strstr(string1,"J01") != NULL)){
          pclass_sz->MF=3;
          // pclass_sz->integrate_wrt_m180m = 1;
        }
        else  if ((strstr(string1,"T08M200c") != NULL)) { // needs to be before T08!
          pclass_sz->MF=8;
          pclass_sz->integrate_wrt_m200c = 1;
        }
        else  if ((strstr(string1,"T08M200m") != NULL)){
          pclass_sz->MF=4;
          pclass_sz->integrate_wrt_m200m = 1;

          // avoid problems here
          if  (strstr(string1,"T08M500c") != NULL){
          pclass_sz->MF=5;
          pclass_sz->integrate_wrt_m500c = 1;
          pclass_sz->integrate_wrt_m200c = 0;
          }
        }
        else  if ((strstr(string1,"T08M500c") != NULL)){
          pclass_sz->MF=5;
          pclass_sz->integrate_wrt_m500c = 1;
        }
        else  if ((strstr(string1,"M1600") != NULL)){
          pclass_sz->MF=6;
          // pclass_sz->integrate_wrt_m1600m = 1;
        }
        else  if ((strstr(string1,"B16M500c") != NULL)){ //not working yet
          pclass_sz->MF=7;
          pclass_sz->integrate_wrt_m500c = 1;
        }

     }

      // class_read_int("integrate_wrt_mvir",pclass_sz->integrate_wrt_mvir);
      // class_read_int("integrate_wrt_m500c",pclass_sz->integrate_wrt_m500c);
      // class_read_int("integrate_wrt_m200m",pclass_sz->integrate_wrt_m200m);
      // class_read_int("integrate_wrt_m200c",pclass_sz->integrate_wrt_m200c);
     if (pclass_sz->integrate_wrt_mvir == 1) {
        pclass_sz->integrate_wrt_m500c = 0;
        pclass_sz->integrate_wrt_m200m = 0;
        pclass_sz->integrate_wrt_m200c = 0;
      }
      if (pclass_sz->integrate_wrt_m500c == 1){
        pclass_sz->integrate_wrt_mvir = 0;
        pclass_sz->integrate_wrt_m200m = 0;
        pclass_sz->integrate_wrt_m200c = 0;

        pclass_sz->delta_def_cib = 2;
        pclass_sz->delta_def_galaxies = 2;
        pclass_sz->delta_def_matter_density = 2;
        pclass_sz->delta_def_electron_density = 2;
        pclass_sz->delta_def_electron_pressure = 2;
        pclass_sz->delta_def_custom1 = 2;
      }
      // m200c is the default
      // if (pclass_sz->integrate_wrt_m200c == 1){
      //   pclass_sz->integrate_wrt_mvir = 0;
      //   pclass_sz->integrate_wrt_m200m = 0;
      //   pclass_sz->integrate_wrt_m500c = 0;
      // }
      if (pclass_sz->integrate_wrt_m200m == 1){
        pclass_sz->integrate_wrt_mvir = 0;
        pclass_sz->integrate_wrt_m200c = 0;
        pclass_sz->integrate_wrt_m500c = 0;


        pclass_sz->delta_def_cib = 0;
        pclass_sz->delta_def_galaxies = 0;
        pclass_sz->delta_def_matter_density = 0;
        pclass_sz->delta_def_electron_density = 0;
        pclass_sz->delta_def_electron_pressure = 0;
        pclass_sz->delta_def_custom1 = 0;
      }

      class_call(parser_read_string(pfc,"profile_for_matter_density",&string1,&flag1,errmsg),errmsg,errmsg);

        if ((strstr(string1,"nfw_with_power_law") != NULL)){
          pclass_sz->profile_matter_density=1;
          class_read_double("matter_nfw_power_law_index",pclass_sz->matter_nfw_power_law_index);
          if (pclass_sz->sz_verbose>0){
            printf("using matter_nfw_power_law with n=%.3e\n",pclass_sz->matter_nfw_power_law_index);
          }
        }

        else if (flag1 == _TRUE_) {
          if ((strstr(string1,"nfw") != NULL)){
            pclass_sz->profile_matter_density=0;
            if (pclass_sz->sz_verbose>0){
              printf("using standard nfw profile for matter profile.\n",pclass_sz->matter_nfw_power_law_index);
            }
          }
          }


      class_call(parser_read_string(pfc,"delta for galaxies",&string1,&flag1,errmsg),errmsg,errmsg);
      if (flag1 == _FALSE_) {
        class_call(parser_read_string(pfc,"delta_for_galaxies",&string1,&flag1,errmsg),errmsg,errmsg);
      }
      if (flag1 == _TRUE_) {
        if ((strstr(string1,"200m") != NULL))
          pclass_sz->delta_def_galaxies=0;
        else  if ((strstr(string1,"200c") != NULL))
          pclass_sz->delta_def_galaxies=1;
        else  if ((strstr(string1,"500c") != NULL))
          pclass_sz->delta_def_galaxies=2;
          }


      class_call(parser_read_string(pfc,"delta for cib",&string1,&flag1,errmsg),errmsg,errmsg);
      if (flag1 == _FALSE_) {
        class_call(parser_read_string(pfc,"delta_for_cib",&string1,&flag1,errmsg),errmsg,errmsg);
      }
      if (flag1 == _TRUE_) {
        if ((strstr(string1,"200m") != NULL))
          pclass_sz->delta_def_cib=0;
        else  if ((strstr(string1,"200c") != NULL))
          pclass_sz->delta_def_cib=1;
        else  if ((strstr(string1,"500c") != NULL))
          pclass_sz->delta_def_cib=2;
          }

      class_call(parser_read_string(pfc,"delta for matter density",&string1,&flag1,errmsg),errmsg,errmsg);
      if (flag1 == _FALSE_) {
        class_call(parser_read_string(pfc,"delta_for_matter_density",&string1,&flag1,errmsg),errmsg,errmsg);
      }
      if (flag1 == _FALSE_) {
        class_call(parser_read_string(pfc,"delta_matter_density",&string1,&flag1,errmsg),errmsg,errmsg);
      }
      if (flag1 == _FALSE_) {
        class_call(parser_read_string(pfc,"delta_matter",&string1,&flag1,errmsg),errmsg,errmsg);
      }  
      if (flag1 == _FALSE_) {
        class_call(parser_read_string(pfc,"delta_def_matter_density",&string1,&flag1,errmsg),errmsg,errmsg);
      }
      if (flag1 == _FALSE_) {
        class_call(parser_read_string(pfc,"delta_def_matter",&string1,&flag1,errmsg),errmsg,errmsg);
      }     
      if (flag1 == _TRUE_) {
        if ((strstr(string1,"200m") != NULL))
          pclass_sz->delta_def_matter_density=0;
        else  if ((strstr(string1,"200c") != NULL))
          pclass_sz->delta_def_matter_density=1;
        else  if ((strstr(string1,"500c") != NULL))
          pclass_sz->delta_def_matter_density=2;
        else  if ((strstr(string1,"vir") != NULL))
          pclass_sz->delta_def_matter_density=3;
          }

      class_call(parser_read_string(pfc,"delta_for_custom1",&string1,&flag1,errmsg),errmsg,errmsg);
      if (flag1 == _TRUE_) {
        if ((strstr(string1,"200m") != NULL))
          pclass_sz->delta_def_custom1=0;
        else  if ((strstr(string1,"200c") != NULL))
          pclass_sz->delta_def_custom1=1;
        else  if ((strstr(string1,"500c") != NULL))
          pclass_sz->delta_def_custom1=2;
        else  if ((strstr(string1,"vir") != NULL))
          pclass_sz->delta_def_custom1=3;
          }


      class_call(parser_read_string(pfc,"delta for electron density",&string1,&flag1,errmsg),errmsg,errmsg);
      if (flag1 == _FALSE_) {
        class_call(parser_read_string(pfc,"delta_for_electron_density",&string1,&flag1,errmsg),errmsg,errmsg);
      } 
      if (flag1 == _TRUE_) {
        if ((strstr(string1,"200m") != NULL))
          pclass_sz->delta_def_electron_density=0;
        else  if ((strstr(string1,"200c") != NULL))
          pclass_sz->delta_def_electron_density=1;
        else  if ((strstr(string1,"500c") != NULL))
          pclass_sz->delta_def_electron_density=2;
          }

      class_call(parser_read_string(pfc,"delta for electron pressure",&string1,&flag1,errmsg),errmsg,errmsg);
      if (flag1 == _FALSE_) {
        class_call(parser_read_string(pfc,"delta_for_electron_pressure",&string1,&flag1,errmsg),errmsg,errmsg);
      } 
      if (flag1 == _TRUE_) {
        if ((strstr(string1,"200m") != NULL))
          pclass_sz->delta_def_electron_pressure=0;
        else  if ((strstr(string1,"200c") != NULL))
          pclass_sz->delta_def_electron_pressure=1;
        else  if ((strstr(string1,"500c") != NULL))
          pclass_sz->delta_def_electron_pressure=2;
          }



      class_read_string("UNWISE_dndz_file",pclass_sz->UNWISE_dndz_file);
      class_read_string("UNWISE_fdndz_file",pclass_sz->UNWISE_dndz_file);
      // class_read_string("path_to_class",pclass_sz->path_to_class);
      class_read_string("sz_selection_function_thetas_file",pclass_sz->SO_thetas_file);
      class_read_string("sz_selection_function_skyfracs_file",pclass_sz->SO_skyfracs_file);
      class_read_string("sz_selection_function_ylims_file",pclass_sz->SO_ylims_file);

      class_read_string("append_name_trispectrum_ref",pclass_sz->append_name_trispectrum_ref);
      class_read_string("path to reference trispectrum for likelihood code",pclass_sz->path_to_ref_trispectrum_for_cobaya);
      class_read_string("root",pclass_sz->root);
      class_read_string("root",pcsz->root);


      class_read_double("f_free",pclass_sz->f_free);
      class_read_double("f_b_gas",pclass_sz->f_b_gas);

      class_read_int("compute_ksz2ksz2",pclass_sz->compute_ksz2ksz2);

      class_read_int("overwrite_clpp_with_limber",psp->overwrite_clpp_with_limber);

      class_call(parser_read_string(pfc,"write sz results to files",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);

       if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {

         pclass_sz->write_sz = _TRUE_;

       }



    /* // Eq. 15 or 16 of KA20 */
    class_call(parser_read_string(pfc,"use_hod",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
   if (flag1 == _TRUE_) {
      if ((strstr(string1,"yes") != NULL))
        pclass_sz->use_hod=1;
      else  if ((strstr(string1,"no") != NULL))
        pclass_sz->use_hod=0;
      }

    /* // Use analytical truncated nfw rather than numerical fourier transform */
    class_call(parser_read_string(pfc,"use_analytical_truncated_nfw",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
   if (flag1 == _TRUE_) {
      if ((strstr(string1,"yes") != NULL))
        pclass_sz->use_analytical_truncated_nfw=1;
      else  if ((strstr(string1,"no") != NULL))
        pclass_sz->use_analytical_truncated_nfw=0;
      }

      /* tau profile kSZ */
      class_call(parser_read_string(pfc,"gas profile",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      if (flag1 == _TRUE_) {
        if ((strstr(string1,"nfw") != NULL))
          pclass_sz->tau_profile=0;
        else  if ((strstr(string1,"B16") != NULL)){
          pclass_sz->tau_profile=1;
          pclass_sz->use_analytical_truncated_nfw=0; // overwrite previous
        }
        else  if ((strstr(string1,"BCM") != NULL)){
          pclass_sz->tau_profile=2;
          // pclass_sz->use_analytical_truncated_nfw=0; // overwrite previous
        }

        }

      /* tau profile mode kSZ */
      class_call(parser_read_string(pfc,"gas profile mode",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      if (flag1 == _TRUE_) {
        if ((strstr(string1,"agn") != NULL)){
          pclass_sz->tau_profile_mode=0;
            pclass_sz->A_rho0 = 4.e3;
            pclass_sz->A_alpha = 0.88;
            pclass_sz->A_beta = 3.83;

            pclass_sz->alpha_m_rho0 = 0.29;
            pclass_sz->alpha_m_alpha = -0.03;
            pclass_sz->alpha_m_beta = 0.04;

            pclass_sz->alpha_z_rho0 = -0.66;
            pclass_sz->alpha_z_alpha = 0.19;
            pclass_sz->alpha_z_beta = -0.025;
            pclass_sz->xc_B16 = 0.5;

            pclass_sz->c_B16 = 0.;
	          pclass_sz->mcut = 1.e14;
            pclass_sz->alphap_m_rho0 = 0.29;
            pclass_sz->alphap_m_alpha = -0.03;
            pclass_sz->alphap_m_beta = 0.04;

            pclass_sz->alpha_c_rho0 = 0.;
            pclass_sz->alpha_c_alpha = 0.;
            pclass_sz->alpha_c_beta = 0.;
        }
        else  if ((strstr(string1,"shock") != NULL)){
          pclass_sz->tau_profile_mode=1;
          //   // shock heating
            pclass_sz->A_rho0 = 1.9e4;
            pclass_sz->A_alpha = 0.70;
            pclass_sz->A_beta = 4.43;

            pclass_sz->alpha_m_rho0 = 0.09;
            pclass_sz->alpha_m_alpha = -0.017;
            pclass_sz->alpha_m_beta = 0.005;

            pclass_sz->alpha_z_rho0 = -0.95;
            pclass_sz->alpha_z_alpha = 0.27;
            pclass_sz->alpha_z_beta = 0.037;
            pclass_sz->xc_B16 = 0.5;

            pclass_sz->c_B16 = 0.;
	          pclass_sz->mcut = 1.e14;
            pclass_sz->alphap_m_rho0 = 0.09;
            pclass_sz->alphap_m_alpha = -0.017;
            pclass_sz->alphap_m_beta = 0.005;

            pclass_sz->alpha_c_rho0 = 0.;
            pclass_sz->alpha_c_alpha = 0.;
            pclass_sz->alpha_c_beta = 0.;
        }
        else if ((strstr(string1,"custom") != NULL)){
          pclass_sz->tau_profile_mode=2;
          class_read_double("A_rho0",pclass_sz->A_rho0);
          class_read_double("A_alpha",pclass_sz->A_alpha);
          class_read_double("A_beta",pclass_sz->A_beta);

          class_read_double("alpha_m_rho0",pclass_sz->alpha_m_rho0);
          class_read_double("alpha_m_alpha",pclass_sz->alpha_m_alpha);
          class_read_double("alpha_m_beta",pclass_sz->alpha_m_beta);

          class_read_double("alpha_z_rho0",pclass_sz->alpha_z_rho0);
          class_read_double("alpha_z_alpha",pclass_sz->alpha_z_alpha);
          class_read_double("alpha_z_beta",pclass_sz->alpha_z_beta);


	        // B.H.
          class_read_double("cp_B16",pclass_sz->c_B16);
          class_read_double("mcut",pclass_sz->mcut);
          class_read_double("alphap_m_rho0",pclass_sz->alphap_m_rho0);
          class_read_double("alphap_m_alpha",pclass_sz->alphap_m_alpha);
          class_read_double("alphap_m_beta",pclass_sz->alphap_m_beta);

          class_read_double("alpha_c_rho0",pclass_sz->alpha_c_rho0);
          class_read_double("alpha_c_alpha",pclass_sz->alpha_c_alpha);
          class_read_double("alpha_c_beta",pclass_sz->alpha_c_beta);

          class_read_double("gamma_B16",pclass_sz->gamma_B16);
          class_read_double("xc_B16",pclass_sz->xc_B16);
        }
}



    /* galaxy sample */
    class_call(parser_read_string(pfc,"galaxy_sample",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
   if (flag1 == _TRUE_) {
      if ((strstr(string1,"WIxSC") != NULL))
        pclass_sz->galaxy_sample=0;
      else  if ((strstr(string1,"unwise") != NULL))
        pclass_sz->galaxy_sample=1;
      else  if ((strstr(string1,"custom") != NULL))
        pclass_sz->galaxy_sample=2;
      }

      if (pclass_sz->sz_verbose > 3) {
        switch(pclass_sz->galaxy_sample) {
          case 0:
            printf("Galaxy sample: WIxSC\n");
            break;
          case 1:
            printf("Galaxy sample: unwise\n");
            break;
          case 2:
            printf("Galaxy sample: custom\n");
            break;
          default:
            printf("Galaxy sample: unknown\n");
        }
      }


      /* unwise galaxy sample id */
      class_call(parser_read_string(pfc,"galaxy_sample_id",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);
     if (flag1 == _TRUE_) {
        if ((strstr(string1,"red") != NULL))
          pclass_sz->unwise_galaxy_sample_id=0;
        else  if ((strstr(string1,"green") != NULL))
          pclass_sz->unwise_galaxy_sample_id=1;
        else  if ((strstr(string1,"gr_shallow") != NULL))
          pclass_sz->unwise_galaxy_sample_id=2;
        else  if ((strstr(string1,"blue") != NULL))
          pclass_sz->unwise_galaxy_sample_id=3;
        }


    if (pclass_sz->has_pk)
        pclass_sz->need_sigma = _TRUE_;
    // end of class_sz parameters


    if ((strstr(string1,"mTk") != NULL) || (strstr(string1,"MTk") != NULL) || (strstr(string1,"MTK") != NULL) ||
        (strstr(string1,"dTk") != NULL) || (strstr(string1,"DTk") != NULL) || (strstr(string1,"DTK") != NULL)) {
      ppt->has_density_transfers=_TRUE_;
      ppt->has_perturbations = _TRUE_;
    }

    if ((strstr(string1,"vTk") != NULL) || (strstr(string1,"VTk") != NULL) || (strstr(string1,"VTK") != NULL)) {
      ppt->has_velocity_transfers=_TRUE_;
      ppt->has_perturbations = _TRUE_;
    }

  }

  /* The following lines make sure that if perturbations are not computed, idm_dr and idr parameters are still freed */

  if(ppt->has_perturbations == _FALSE_) {

    if (ppt->alpha_idm_dr != NULL)
      free(ppt->alpha_idm_dr);

    if (ppt->beta_idr != NULL)
      free(ppt->beta_idr);
  }

  if (ppt->has_density_transfers == _TRUE_) {
    class_call(parser_read_string(pfc,"extra metric transfer functions",&string1,&flag1,errmsg),
               errmsg,
               errmsg);

    if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"y") != NULL))) {
      ppt->has_metricpotential_transfers = _TRUE_;
    }
  }

  if (ppt->has_cl_cmb_temperature == _TRUE_) {

    class_call(parser_read_string(pfc,"temperature contributions",&string1,&flag1,errmsg),
               errmsg,
               errmsg);

    if (flag1 == _TRUE_) {

      ppt->switch_sw = 0;
      ppt->switch_eisw = 0;
      ppt->switch_lisw = 0;
      ppt->switch_dop = 0;
      ppt->switch_pol = 0;

      if ((strstr(string1,"tsw") != NULL) || (strstr(string1,"TSW") != NULL))
        ppt->switch_sw = 1;
      if ((strstr(string1,"eisw") != NULL) || (strstr(string1,"EISW") != NULL))
        ppt->switch_eisw = 1;
      if ((strstr(string1,"lisw") != NULL) || (strstr(string1,"LISW") != NULL))
        ppt->switch_lisw = 1;
      if ((strstr(string1,"dop") != NULL) || (strstr(string1,"Dop") != NULL))
        ppt->switch_dop = 1;
      if ((strstr(string1,"pol") != NULL) || (strstr(string1,"Pol") != NULL))
        ppt->switch_pol = 1;

      class_test((ppt->switch_sw == 0) && (ppt->switch_eisw == 0) && (ppt->switch_lisw == 0) && (ppt->switch_dop == 0) && (ppt->switch_pol == 0),
                 errmsg,
                 "In the field 'output', you selected CMB temperature, but in the field 'temperature contributions', you removed all contributions");

      class_read_double("early/late isw redshift",ppt->eisw_lisw_split_z);

    }

  }

  if (ppt->has_cl_number_count == _TRUE_) {

    class_call(parser_read_string(pfc,"number count contributions",&string1,&flag1,errmsg),
               errmsg,
               errmsg);

    if (flag1 == _TRUE_) {

      if (strstr(string1,"density") != NULL)
        ppt->has_nc_density = _TRUE_;
      if (strstr(string1,"rsd") != NULL)
        ppt->has_nc_rsd = _TRUE_;
      if (strstr(string1,"lensing") != NULL)
        ppt->has_nc_lens = _TRUE_;
      if (strstr(string1,"gr") != NULL)
        ppt->has_nc_gr = _TRUE_;

      class_test((ppt->has_nc_density == _FALSE_) && (ppt->has_nc_rsd == _FALSE_) && (ppt->has_nc_lens == _FALSE_) && (ppt->has_nc_gr == _FALSE_),
                 errmsg,
                 "In the field 'output', you selected number count Cl's, but in the field 'number count contributions', you removed all contributions");

    }

    else {
      /* default: only the density contribution */
      ppt->has_nc_density = _TRUE_;
    }
  }

  if (ppt->has_perturbations == _TRUE_) {

    /* perturbed recombination */
    class_call(parser_read_string(pfc,
                                  "perturbed recombination",
                                  &(string1),
                                  &(flag1),
                                  errmsg),
               errmsg,
               errmsg);

    if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {
      ppt->has_perturbed_recombination = _TRUE_;
    }

    /* modes */
    class_call(parser_read_string(pfc,"modes",&string1,&flag1,errmsg),
               errmsg,
               errmsg);

    if (flag1 == _TRUE_) {

      /* if no modes are specified, the default is has_scalars=_TRUE_;
         but if they are specified we should reset has_scalars to _FALSE_ before reading */
      ppt->has_scalars=_FALSE_;

      if ((strstr(string1,"s") != NULL) || (strstr(string1,"S") != NULL))
        ppt->has_scalars=_TRUE_;

      if ((strstr(string1,"v") != NULL) || (strstr(string1,"V") != NULL))
        ppt->has_vectors=_TRUE_;

      if ((strstr(string1,"t") != NULL) || (strstr(string1,"T") != NULL))
        ppt->has_tensors=_TRUE_;

      class_test(class_none_of_three(ppt->has_scalars,ppt->has_vectors,ppt->has_tensors),
                 errmsg,
                 "You wrote: modes='%s'. Could not identify any of the modes ('s', 'v', 't') in such input",string1);
    }

    if (ppt->has_scalars == _TRUE_) {

      class_call(parser_read_string(pfc,"ic",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);

      if (flag1 == _TRUE_) {

        /* if no initial conditions are specified, the default is has_ad=_TRUE_;
           but if they are specified we should reset has_ad to _FALSE_ before reading */
        ppt->has_ad=_FALSE_;

        if ((strstr(string1,"ad") != NULL) || (strstr(string1,"AD") != NULL))
          ppt->has_ad=_TRUE_;

        if ((strstr(string1,"bi") != NULL) || (strstr(string1,"BI") != NULL))
          ppt->has_bi=_TRUE_;

        if ((strstr(string1,"cdi") != NULL) || (strstr(string1,"CDI") != NULL))
          ppt->has_cdi=_TRUE_;

        if ((strstr(string1,"nid") != NULL) || (strstr(string1,"NID") != NULL))
          ppt->has_nid=_TRUE_;

        if ((strstr(string1,"niv") != NULL) || (strstr(string1,"NIV") != NULL))
          ppt->has_niv=_TRUE_;

        class_test(ppt->has_ad==_FALSE_ && ppt->has_bi ==_FALSE_ && ppt->has_cdi ==_FALSE_ && ppt->has_nid ==_FALSE_ && ppt->has_niv ==_FALSE_,
                   errmsg,
                   "You wrote: ic='%s'. Could not identify any of the initial conditions ('ad', 'bi', 'cdi', 'nid', 'niv') in such input",string1);

      }
    }

    else {

      class_test(ppt->has_cl_cmb_lensing_potential == _TRUE_,
                 errmsg,
                 "Inconsistency: you want C_l's for cmb lensing potential, but no scalar modes\n");

      class_test(ppt->has_pk_matter == _TRUE_,
                 errmsg,
                 "Inconsistency: you want P(k) of matter, but no scalar modes\n");

    }

    if (ppt->has_vectors == _TRUE_){

      class_test((ppt->has_cl_cmb_temperature == _FALSE_) && (ppt->has_cl_cmb_polarization == _FALSE_),
                 errmsg,
                 "inconsistent input: you asked for vectors, so you should have at least one non-zero tensor source type (temperature or polarization). Please adjust your input.");

    }

    if (ppt->has_tensors == _TRUE_){

      class_test((ppt->has_cl_cmb_temperature == _FALSE_) && (ppt->has_cl_cmb_polarization == _FALSE_),
                 errmsg,
                 "inconsistent input: you asked for tensors, so you should have at least one non-zero tensor source type (temperature or polarization). Please adjust your input.");

    }
  }

  /** (d) define the primordial spectrum */

  class_call(parser_read_string(pfc,"P_k_ini type",&string1,&flag1,errmsg),
             errmsg,
             errmsg);

  if (flag1 == _TRUE_) {
    flag2=_FALSE_;
    if (strcmp(string1,"analytic_Pk") == 0) {
      ppm->primordial_spec_type = analytic_Pk;
      flag2=_TRUE_;
    }
    if (strcmp(string1,"two_scales") == 0) {
      ppm->primordial_spec_type = two_scales;
      flag2=_TRUE_;
    }
    if (strcmp(string1,"inflation_V") == 0) {
      ppm->primordial_spec_type = inflation_V;
      flag2=_TRUE_;
    }
    if (strcmp(string1,"inflation_H") == 0) {
      ppm->primordial_spec_type = inflation_H;
      flag2=_TRUE_;
    }
    if (strcmp(string1,"inflation_V_end") == 0) {
      ppm->primordial_spec_type = inflation_V_end;
      flag2=_TRUE_;
    }
    if (strcmp(string1,"external_Pk") == 0) {
      ppm->primordial_spec_type = external_Pk;
      flag2=_TRUE_;
    }
    class_test(flag2==_FALSE_,
               errmsg,
               "could not identify primordial spectrum type, check that it is one of 'analytic_pk', 'two_scales', 'inflation_V', 'inflation_H', 'external_Pk'...");
  }

  class_read_double("k_pivot",ppm->k_pivot);

  if (ppm->primordial_spec_type == two_scales) {

    class_read_double("k1",k1);
    class_read_double("k2",k2);
    class_test(k1<=0.,errmsg,"enter strictly positive scale k1");
    class_test(k2<=0.,errmsg,"enter strictly positive scale k2");

    if (ppt->has_scalars == _TRUE_) {

      class_read_double("P_{RR}^1",prr1);
      class_read_double("P_{RR}^2",prr2);
      class_test(prr1<=0.,errmsg,"enter strictly positive scale P_{RR}^1");
      class_test(prr2<=0.,errmsg,"enter strictly positive scale P_{RR}^2");

      ppm->n_s = log(prr2/prr1)/log(k2/k1)+1.;
      ppm->A_s = prr1*exp((ppm->n_s-1.)*log(ppm->k_pivot/k1));

      if ((ppt->has_bi == _TRUE_) ||
          (ppt->has_cdi == _TRUE_) ||
          (ppt->has_nid == _TRUE_) ||
          (ppt->has_niv == _TRUE_)) {

        class_read_double("P_{II}^1",pii1);
        class_read_double("P_{II}^2",pii2);
        class_read_double("P_{RI}^1",pri1);
        class_read_double("|P_{RI}^2|",pri2);

        class_test(pii1 <= 0.,
                   errmsg,
                   "since you request iso modes, you should have P_{ii}^1 strictly positive");
        class_test(pii2 < 0.,
                   errmsg,
                   "since you request iso modes, you should have P_{ii}^2 positive or eventually null");
        class_test(pri2 < 0.,
                   errmsg,
                   "by definition, you should have |P_{ri}^2| positive or eventually null");

        flag1 = _FALSE_;

        class_call(parser_read_string(pfc,"special iso",&string1,&flag1,errmsg),
                   errmsg,
                   errmsg);

        /* axion case, only one iso parameter: piir1  */
        if ((flag1 == _TRUE_) && (strstr(string1,"axion") != NULL)) {
          n_iso = 1.;
          n_cor = 0.;
          c_cor = 0.;
        }
        /* curvaton case, only one iso parameter: piir1  */
        else if ((flag1 == _TRUE_) && (strstr(string1,"anticurvaton") != NULL)) {
          n_iso = ppm->n_s;
          n_cor = 0.;
          c_cor = 1.;
        }
        /* inverted-correlation-curvaton case, only one iso parameter: piir1  */
        else if ((flag1 == _TRUE_) && (strstr(string1,"curvaton") != NULL)) {
          n_iso = ppm->n_s;
          n_cor = 0.;
          c_cor = -1.;
        }
        /* general case, but if pii2 or pri2=0 the code interprets it
           as a request for n_iso=n_ad or n_cor=0 respectively */
        else {
          if (pii2 == 0.) {
            n_iso = ppm->n_s;
          }
          else {
            class_test((pii1==0.) || (pii2 == 0.) || (pii1*pii2<0.),errmsg,"should NEVER happen");
            n_iso = log(pii2/pii1)/log(k2/k1)+1.;
          }
          class_test(pri1==0,errmsg,"the general isocurvature case requires a non-zero P_{RI}^1");
          if (pri2 == 0.) {
            n_cor = 0.;
          }
          else {
            class_test((pri1==0.) || (pri2 <= 0.) || (pii1*pii2<0),errmsg,"should NEVER happen");
            n_cor = log(pri2/fabs(pri1))/log(k2/k1)-0.5*(ppm->n_s+n_iso-2.);
          }
          class_test((pii1*prr1<=0.),errmsg,"should NEVER happen");
          class_test(fabs(pri1)/sqrt(pii1*prr1)>1,errmsg,"too large ad-iso cross-correlation in k1");
          class_test(fabs(pri1)/sqrt(pii1*prr1)*exp(n_cor*log(k2/k1))>1,errmsg,"too large ad-iso cross-correlation in k2");
          c_cor = -pri1/sqrt(pii1*prr1)*exp(n_cor*log(ppm->k_pivot/k1));
        }
        /* formula for f_iso valid in all cases */
        class_test((pii1==0.) || (prr1 == 0.) || (pii1*prr1<0.),errmsg,"should NEVER happen");
        f_iso = sqrt(pii1/prr1)*exp(0.5*(n_iso-ppm->n_s)*log(ppm->k_pivot/k1));

      }

      if (ppt->has_bi == _TRUE_) {
        ppm->f_bi = f_iso;
        ppm->n_bi = n_iso;
        ppm->c_ad_bi = c_cor;
        ppm->n_ad_bi = n_cor;
      }

      if (ppt->has_cdi == _TRUE_) {
        ppm->f_cdi = f_iso;
        ppm->n_cdi = n_iso;
        ppm->c_ad_cdi = c_cor;
        ppm->n_ad_cdi = n_cor;
      }

      if (ppt->has_nid == _TRUE_) {
        ppm->f_nid = f_iso;
        ppm->n_nid = n_iso;
        ppm->c_ad_nid = c_cor;
        ppm->n_ad_nid = n_cor;
      }

      if (ppt->has_niv == _TRUE_) {
        ppm->f_niv = f_iso;
        ppm->n_niv = n_iso;
        ppm->c_ad_niv = c_cor;
        ppm->n_ad_niv = n_cor;
      }
    }

    ppm->primordial_spec_type = analytic_Pk;

  }

  else if (ppm->primordial_spec_type == analytic_Pk) {

    if (ppt->has_scalars == _TRUE_) {

      class_call(parser_read_double(pfc,"A_s",&param1,&flag1,errmsg),
                 errmsg,
                 errmsg);
      class_call(parser_read_double(pfc,"ln10^{10}A_s",&param2,&flag2,errmsg),
                 errmsg,
                 errmsg);
      class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),
                 errmsg,
                 "In input file, you cannot enter both A_s and ln10^{10}A_s, choose one");
      if (flag1 == _TRUE_)
        ppm->A_s = param1;
      else if (flag2 == _TRUE_)
        ppm->A_s = exp(param2)*1.e-10;

      if (ppt->has_ad == _TRUE_) {

        class_read_double("n_s",ppm->n_s);
        class_read_double("alpha_s",ppm->alpha_s);

      }

      if (ppt->has_bi == _TRUE_) {

        class_read_double("f_bi",ppm->f_bi);
        class_read_double("n_bi",ppm->n_bi);
        class_read_double("alpha_bi",ppm->alpha_bi);

      }

      if (ppt->has_cdi == _TRUE_) {

        class_read_double("f_cdi",ppm->f_cdi);
        class_read_double("n_cdi",ppm->n_cdi);
        class_read_double("alpha_cdi",ppm->alpha_cdi);

      }

      if (ppt->has_nid == _TRUE_) {

        class_read_double("f_nid",ppm->f_nid);
        class_read_double("n_nid",ppm->n_nid);
        class_read_double("alpha_nid",ppm->alpha_nid);

      }

      if (ppt->has_niv == _TRUE_) {

        class_read_double("f_niv",ppm->f_niv);
        class_read_double("n_niv",ppm->n_niv);
        class_read_double("alpha_niv",ppm->alpha_niv);

      }

      if ((ppt->has_ad == _TRUE_) && (ppt->has_bi == _TRUE_)) {
        class_read_double_one_of_two("c_ad_bi","c_bi_ad",ppm->c_ad_bi);
        class_read_double_one_of_two("n_ad_bi","n_bi_ad",ppm->n_ad_bi);
        class_read_double_one_of_two("alpha_ad_bi","alpha_bi_ad",ppm->alpha_ad_bi);
      }

      if ((ppt->has_ad == _TRUE_) && (ppt->has_cdi == _TRUE_)) {
        class_read_double_one_of_two("c_ad_cdi","c_cdi_ad",ppm->c_ad_cdi);
        class_read_double_one_of_two("n_ad_cdi","n_cdi_ad",ppm->n_ad_cdi);
        class_read_double_one_of_two("alpha_ad_cdi","alpha_cdi_ad",ppm->alpha_ad_cdi);
      }

      if ((ppt->has_ad == _TRUE_) && (ppt->has_nid == _TRUE_)) {
        class_read_double_one_of_two("c_ad_nid","c_nid_ad",ppm->c_ad_nid);
        class_read_double_one_of_two("n_ad_nid","n_nid_ad",ppm->n_ad_nid);
        class_read_double_one_of_two("alpha_ad_nid","alpha_nid_ad",ppm->alpha_ad_nid);
      }

      if ((ppt->has_ad == _TRUE_) && (ppt->has_niv == _TRUE_)) {
        class_read_double_one_of_two("c_ad_niv","c_niv_ad",ppm->c_ad_niv);
        class_read_double_one_of_two("n_ad_niv","n_niv_ad",ppm->n_ad_niv);
        class_read_double_one_of_two("alpha_ad_niv","alpha_niv_ad",ppm->alpha_ad_niv);
      }

      if ((ppt->has_bi == _TRUE_) && (ppt->has_cdi == _TRUE_)) {
        class_read_double_one_of_two("c_bi_cdi","c_cdi_bi",ppm->c_bi_cdi);
        class_read_double_one_of_two("n_bi_cdi","n_cdi_bi",ppm->n_bi_cdi);
        class_read_double_one_of_two("alpha_bi_cdi","alpha_cdi_bi",ppm->alpha_bi_cdi);
      }

      if ((ppt->has_bi == _TRUE_) && (ppt->has_nid == _TRUE_)) {
        class_read_double_one_of_two("c_bi_nid","c_nid_bi",ppm->c_bi_nid);
        class_read_double_one_of_two("n_bi_nid","n_nid_bi",ppm->n_bi_nid);
        class_read_double_one_of_two("alpha_bi_nid","alpha_nid_bi",ppm->alpha_bi_nid);
      }

      if ((ppt->has_bi == _TRUE_) && (ppt->has_niv == _TRUE_)) {
        class_read_double_one_of_two("c_bi_niv","c_niv_bi",ppm->c_bi_niv);
        class_read_double_one_of_two("n_bi_niv","n_niv_bi",ppm->n_bi_niv);
        class_read_double_one_of_two("alpha_bi_niv","alpha_niv_bi",ppm->alpha_bi_niv);
      }

      if ((ppt->has_cdi == _TRUE_) && (ppt->has_nid == _TRUE_)) {
        class_read_double_one_of_two("c_cdi_nid","c_nid_cdi",ppm->c_cdi_nid);
        class_read_double_one_of_two("n_cdi_nid","n_nid_cdi",ppm->n_cdi_nid);
        class_read_double_one_of_two("alpha_cdi_nid","alpha_nid_cdi",ppm->alpha_cdi_nid);
      }

      if ((ppt->has_cdi == _TRUE_) && (ppt->has_niv == _TRUE_)) {
        class_read_double_one_of_two("c_cdi_niv","c_niv_cdi",ppm->c_cdi_niv);
        class_read_double_one_of_two("n_cdi_niv","n_niv_cdi",ppm->n_cdi_niv);
        class_read_double_one_of_two("alpha_cdi_niv","alpha_niv_cdi",ppm->alpha_cdi_niv);
      }

      if ((ppt->has_nid == _TRUE_) && (ppt->has_niv == _TRUE_)) {
        class_read_double_one_of_two("c_nid_niv","c_niv_nid",ppm->c_nid_niv);
        class_read_double_one_of_two("n_nid_niv","n_niv_nid",ppm->n_nid_niv);
        class_read_double_one_of_two("alpha_nid_niv","alpha_niv_nid",ppm->alpha_nid_niv);
      }

    }

    if (ppt->has_tensors == _TRUE_) {

      class_read_double("r",ppm->r);

      if (ppt->has_scalars == _FALSE_) {
        class_read_double("A_s",ppm->A_s);
      }

      if (ppm->r <= 0) {
        ppt->has_tensors = _FALSE_;
      }
      else {

        class_call(parser_read_string(pfc,"n_t",&string1,&flag1,errmsg),
                   errmsg,
                   errmsg);

        if ((flag1 == _TRUE_) && !((strstr(string1,"SCC") != NULL) || (strstr(string1,"scc") != NULL))) {
          class_read_double("n_t",ppm->n_t);
        }
        else {
          /* enforce single slow-roll self-consistency condition (order 2 in slow-roll) */
          ppm->n_t = -ppm->r/8.*(2.-ppm->r/8.-ppm->n_s);
        }

        class_call(parser_read_string(pfc,"alpha_t",&string1,&flag1,errmsg),
                   errmsg,
                   errmsg);

        if ((flag1 == _TRUE_) && !((strstr(string1,"SCC") != NULL) || (strstr(string1,"scc") != NULL))) {
          class_read_double("alpha_t",ppm->alpha_t);
        }
        else {
          /* enforce single slow-roll self-consistency condition (order 2 in slow-roll) */
          ppm->alpha_t = ppm->r/8.*(ppm->r/8.+ppm->n_s-1.);
        }
      }
    }
  }

  else if ((ppm->primordial_spec_type == inflation_V) || (ppm->primordial_spec_type == inflation_H)) {

    if (ppm->primordial_spec_type == inflation_V) {

      class_call(parser_read_string(pfc,"potential",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);

      /* only polynomial coded so far: no need to interpret string1 **/

      class_call(parser_read_string(pfc,"PSR_0",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);

      if (flag1 == _TRUE_) {

        PSR0=0.;
        PSR1=0.;
        PSR2=0.;
        PSR3=0.;
        PSR4=0.;

        class_read_double("PSR_0",PSR0);
        class_read_double("PSR_1",PSR1);
        class_read_double("PSR_2",PSR2);
        class_read_double("PSR_3",PSR3);
        class_read_double("PSR_4",PSR4);

        class_test(PSR0 <= 0.,
                   errmsg,
                   "inconsistent parametrization of polynomial inflation potential");
        class_test(PSR1 <= 0.,
                   errmsg,
                   "inconsistent parametrization of polynomial inflation potential");

        R0 = PSR0;
        R1 = PSR1*16.*_PI_;
        R2 = PSR2*8.*_PI_;
        R3 = PSR3*pow(8.*_PI_,2);
        R4 = PSR4*pow(8.*_PI_,3);

        ppm->V0 = R0*R1*3./128./_PI_;
        ppm->V1 = -sqrt(R1)*ppm->V0;
        ppm->V2 = R2*ppm->V0;
        ppm->V3 = R3*ppm->V0*ppm->V0/ppm->V1;
        ppm->V4 = R4*ppm->V0/R1;
      }

      else {

        class_call(parser_read_string(pfc,"R_0",&string1,&flag1,errmsg),
                   errmsg,
                   errmsg);

        if (flag1 == _TRUE_) {

          R0=0.;
          R1=0.;
          R2=0.;
          R3=0.;
          R4=0.;

          class_read_double("R_0",R0);
          class_read_double("R_1",R1);
          class_read_double("R_2",R2);
          class_read_double("R_3",R3);
          class_read_double("R_4",R4);

          class_test(R0 <= 0.,
                     errmsg,
                     "inconsistent parametrization of polynomial inflation potential");
          class_test(R1 <= 0.,
                     errmsg,
                     "inconsistent parametrization of polynomial inflation potential");

          ppm->V0 = R0*R1*3./128./_PI_;
          ppm->V1 = -sqrt(R1)*ppm->V0;
          ppm->V2 = R2*ppm->V0;
          ppm->V3 = R3*ppm->V0*ppm->V0/ppm->V1;
          ppm->V4 = R4*ppm->V0/R1;
        }

        else {

          class_read_double("V_0",ppm->V0);
          class_read_double("V_1",ppm->V1);
          class_read_double("V_2",ppm->V2);
          class_read_double("V_3",ppm->V3);
          class_read_double("V_4",ppm->V4);

        }
      }
    }

    else {

      class_call(parser_read_string(pfc,"HSR_0",&string1,&flag1,errmsg),
                 errmsg,
                 errmsg);

      if (flag1 == _TRUE_) {

        HSR0=0.;
        HSR1=0.;
        HSR2=0.;
        HSR3=0.;
        HSR4=0.;

        class_read_double("HSR_0",HSR0);
        class_read_double("HSR_1",HSR1);
        class_read_double("HSR_2",HSR2);
        class_read_double("HSR_3",HSR3);
        class_read_double("HSR_4",HSR4);

        ppm->H0 = sqrt(HSR0*HSR1*_PI_);
        ppm->H1 = -sqrt(4.*_PI_*HSR1)*ppm->H0;
        ppm->H2 = 4.*_PI_*HSR2*ppm->H0;
        ppm->H3 = 4.*_PI_*HSR3*ppm->H0*ppm->H0/ppm->H1;
        ppm->H4 = 4.*_PI_*HSR4*ppm->H0*ppm->H0*ppm->H0/ppm->H1/ppm->H1;

      }
      else {

        class_read_double("H_0",ppm->H0);
        class_read_double("H_1",ppm->H1);
        class_read_double("H_2",ppm->H2);
        class_read_double("H_3",ppm->H3);
        class_read_double("H_4",ppm->H4);
      }

      class_test(ppm->H0 <= 0.,
                 errmsg,
                 "inconsistent parametrization of polynomial inflation potential");

    }
  }

  else if (ppm->primordial_spec_type == inflation_V_end) {

    class_call(parser_read_string(pfc,"full_potential",&string1,&flag1,errmsg),
               errmsg,
               errmsg);

    if (flag1 == _TRUE_) {
      if (strcmp(string1,"polynomial") == 0) {
        ppm->potential = polynomial;
      }
      else if (strcmp(string1,"higgs_inflation") == 0) {
        ppm->potential = higgs_inflation;
      }
      else {
        class_stop(errmsg,"did not recognize input parameter 'potential': should be one of 'polynomial' or 'higgs_inflation'");
      }
    }

    class_read_double("phi_end",ppm->phi_end);
    class_read_double("Vparam0",ppm->V0);
    class_read_double("Vparam1",ppm->V1);
    class_read_double("Vparam2",ppm->V2);
    class_read_double("Vparam3",ppm->V3);
    class_read_double("Vparam4",ppm->V4);

    class_call(parser_read_string(pfc,"ln_aH_ratio",&string1,&flag1,errmsg),
               errmsg,
               errmsg);

    class_call(parser_read_string(pfc,"N_star",&string2,&flag2,errmsg),
               errmsg,
               errmsg);

    class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),
               errmsg,
               "In input file, you can only enter one of ln_aH_ratio or N_star, the two are not compatible");

    if (flag1 == _TRUE_) {
      if ((strstr(string1,"auto") != NULL) || (strstr(string1,"AUTO") != NULL)) {
        ppm->phi_pivot_method = ln_aH_ratio_auto;
      }
      else {
        ppm->phi_pivot_method = ln_aH_ratio;
        class_read_double("ln_aH_ratio",ppm->phi_pivot_target);
      }
    }

    if (flag2 == _TRUE_) {
      ppm->phi_pivot_method = N_star;
      class_read_double("N_star",ppm->phi_pivot_target);
    }

    class_call(parser_read_string(pfc,"inflation_behavior",&string1,&flag1,errmsg),
               errmsg,
               errmsg);

    if (flag1 == _TRUE_) {
      if (strstr(string1,"numerical") != NULL) {
        ppm->behavior = numerical;
      }
      else if (strstr(string1,"analytical") != NULL) {
        ppm->behavior = analytical;
      }
      else {
        class_stop(errmsg,"Your entry for 'inflation behavior' could not be understood");
      }
    }
  }

  else if (ppm->primordial_spec_type == external_Pk) {
    class_call(parser_read_string(pfc, "command", &(string1), &(flag1), errmsg),
               errmsg, errmsg);
    class_test(strlen(string1) == 0,
               errmsg,
               "You omitted to write a command for the external Pk");

    ppm->command = (char *) malloc (strlen(string1) + 1);
    strcpy(ppm->command, string1);
    class_read_double("custom1",ppm->custom1);
    class_read_double("custom2",ppm->custom2);
    class_read_double("custom3",ppm->custom3);
    class_read_double("custom4",ppm->custom4);
    class_read_double("custom5",ppm->custom5);
    class_read_double("custom6",ppm->custom6);
    class_read_double("custom7",ppm->custom7);
    class_read_double("custom8",ppm->custom8);
    class_read_double("custom9",ppm->custom9);
    class_read_double("custom10",ppm->custom10);
  }

  /* Tests moved from primordial module: */
  if ((ppm->primordial_spec_type == inflation_V) || (ppm->primordial_spec_type == inflation_H) || (ppm->primordial_spec_type == inflation_V_end)) {

    class_test(ppt->has_scalars == _FALSE_,
               errmsg,
               "inflationary module cannot work if you do not ask for scalar modes");

    class_test(ppt->has_vectors == _TRUE_,
               errmsg,
               "inflationary module cannot work if you ask for vector modes");

    class_test(ppt->has_tensors == _FALSE_,
               errmsg,
               "inflationary module cannot work if you do not ask for tensor modes");

    class_test(ppt->has_bi == _TRUE_ || ppt->has_cdi == _TRUE_ || ppt->has_nid == _TRUE_ || ppt->has_niv == _TRUE_,
               errmsg,
               "inflationary module cannot work if you ask for isocurvature modes");
  }

  /** (e) parameters for final spectra */

  if (ppt->has_cls == _TRUE_) {

    if (ppt->has_scalars == _TRUE_) {
      if ((ppt->has_cl_cmb_temperature == _TRUE_) ||
          (ppt->has_cl_cmb_polarization == _TRUE_) ||
          (ppt->has_cl_cmb_lensing_potential == _TRUE_))
        class_read_double("l_max_scalars",ppt->l_scalar_max);

      if ((ppt->has_cl_lensing_potential == _TRUE_) || (ppt->has_cl_number_count == _TRUE_))
        class_read_double("l_max_lss",ppt->l_lss_max);
    }

    if (ppt->has_vectors == _TRUE_) {
      class_read_double("l_max_vectors",ppt->l_vector_max);
    }

    if (ppt->has_tensors == _TRUE_) {
      class_read_double("l_max_tensors",ppt->l_tensor_max);
    }
  }

  class_call(parser_read_string(pfc,
                                "lensing",
                                &(string1),
                                &(flag1),
                                errmsg),
             errmsg,
             errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {

    if ((ppt->has_scalars == _TRUE_) &&
        ((ppt->has_cl_cmb_temperature == _TRUE_) || (ppt->has_cl_cmb_polarization == _TRUE_)) &&
        (ppt->has_cl_cmb_lensing_potential == _TRUE_)) {
      ple->has_lensed_cls = _TRUE_;
    }
    else {
      class_stop(errmsg,"you asked for lensed CMB Cls, but this requires a minimal number of options: 'modes' should include 's', 'output' should include 'tCl' and/or 'pCL', and also, importantly, 'lCl', the CMB lensing potential spectrum. You forgot one of those in your input.");
    }
  }

  if ((ppt->has_scalars == _TRUE_) &&
      (ppt->has_cl_cmb_lensing_potential == _TRUE_)) {

    class_read_double("lcmb_rescale",ptr->lcmb_rescale);
    class_read_double("lcmb_tilt",ptr->lcmb_tilt);
    class_read_double("lcmb_pivot",ptr->lcmb_pivot);

  }

  if ((ppt->has_pk_matter == _TRUE_) || (ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_)) {

    class_call(parser_read_double(pfc,"P_k_max_h/Mpc",&param1,&flag1,errmsg),
               errmsg,
               errmsg);
    class_call(parser_read_double(pfc,"P_k_max_1/Mpc",&param2,&flag2,errmsg),
               errmsg,
               errmsg);
    class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),
               errmsg,
               "In input file, you cannot enter both P_k_max_h/Mpc and P_k_max_1/Mpc, choose one");
    if (flag1 == _TRUE_) {
      ppt->k_max_for_pk=param1*pba->h;
    }
    if (flag2 == _TRUE_) {
      ppt->k_max_for_pk=param2;
    }

    class_call(parser_read_list_of_doubles(pfc,
                                           "z_pk",
                                           &(int1),
                                           &(pointer1),
                                           &flag1,
                                           errmsg),
               errmsg,
               errmsg);

    if (flag1 == _TRUE_) {
      class_test(int1 > _Z_PK_NUM_MAX_,
                 errmsg,
                 "you want to write some output for %d different values of z, hence you should increase _Z_PK_NUM_MAX_ in include/output.h to at least this number",
                 int1);
      pop->z_pk_num = int1;
      for (i=0; i<int1; i++) {
        pop->z_pk[i] = pointer1[i];
      }
      free(pointer1);
    }
  }

  /** Do we want density and velocity transfer functions in Nbody gauge? */
  if ((ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_)){
    class_call(parser_read_string(pfc,"Nbody gauge transfer functions",&string1,&flag1,errmsg),
	       errmsg,
	       errmsg);

    if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"y") != NULL))) {
      ppt->has_Nbody_gauge_transfers = _TRUE_;
    }
  }

  /* deal with selection functions */
  if ((ppt->has_cl_number_count == _TRUE_) || (ppt->has_cl_lensing_potential == _TRUE_)) {

    class_call(parser_read_string(pfc,
                                  "selection",
                                  &(string1),
                                  &(flag1),
                                  errmsg),
               errmsg,
               errmsg);

    if (flag1 == _TRUE_) {
      if (strstr(string1,"gaussian") != NULL) {
        ppt->selection=gaussian;
      }
      else if (strstr(string1,"tophat") != NULL) {
        ppt->selection=tophat;
      }
      else if (strstr(string1,"dirac") != NULL) {
        ppt->selection=dirac;
      }
      else {
        class_stop(errmsg,"In selection function input: type '%s' is unclear",string1);
      }
    }

    class_call(parser_read_list_of_doubles(pfc,
                                           "selection_mean",
                                           &(int1),
                                           &(pointer1),
                                           &flag1,
                                           errmsg),
               errmsg,
               errmsg);

    if ((flag1 == _TRUE_) && (int1>0)) {

      class_test(int1 > _SELECTION_NUM_MAX_,
                 errmsg,
                 "you want to compute density Cl's for %d different bins, hence you should increase _SELECTION_NUM_MAX_ in include/perturbations.h to at least this number",
                 int1);

      ppt->selection_num = int1;
      for (i=0; i<int1; i++) {
        class_test((pointer1[i] < 0.) || (pointer1[i] > 1000.),
                   errmsg,
                   "input of selection functions: you asked for a mean redshift equal to %e, sounds odd",
                   pointer1[i]);
        ppt->selection_mean[i] = pointer1[i];
      }
      free(pointer1);
      /* first set all widths to default; correct eventually later */
      for (i=1; i<int1; i++) {
        class_test(ppt->selection_mean[i]<=ppt->selection_mean[i-1],
                   errmsg,
                   "input of selection functions: the list of mean redshifts must be passed in growing order; you entered %e before %e",ppt->selection_mean[i-1],ppt->selection_mean[i]);
        ppt->selection_width[i] = ppt->selection_width[0];
        ptr->selection_bias[i] = ptr->selection_bias[0];
        ptr->selection_magnification_bias[i] = ptr->selection_magnification_bias[0];
      }

      class_call(parser_read_list_of_doubles(pfc,
                                             "selection_width",
                                             &(int1),
                                             &(pointer1),
                                             &flag1,
                                             errmsg),
                 errmsg,
                 errmsg);

      if ((flag1 == _TRUE_) && (int1>0)) {

        if (int1==1) {
          for (i=0; i<ppt->selection_num; i++) {
            ppt->selection_width[i] = pointer1[0];
          }
        }
        else if (int1==ppt->selection_num) {
          for (i=0; i<int1; i++) {
            ppt->selection_width[i] = pointer1[i];
          }
        }
        else {
          class_stop(errmsg,
                     "In input for selection function, you asked for %d bin centers and %d bin widths; number of bins unclear; you should pass either one bin width (common to all bins) or %d bin widths",
                     ppt->selection_num,int1,ppt->selection_num);
        }
        free(pointer1);
      }

      class_call(parser_read_list_of_doubles(pfc,
                                             "selection_bias",
                                             &(int1),
                                             &(pointer1),
                                             &flag1,
                                             errmsg),
                 errmsg,
                 errmsg);

      if ((flag1 == _TRUE_) && (int1>0)) {

        if (int1==1) {
          for (i=0; i<ppt->selection_num; i++) {
            ptr->selection_bias[i] = pointer1[0];
          }
        }
        else if (int1==ppt->selection_num) {
          for (i=0; i<int1; i++) {
            ptr->selection_bias[i] = pointer1[i];
          }
        }
        else {
          class_stop(errmsg,
                     "In input for selection function, you asked for %d bin centers and %d bin biases; number of bins unclear; you should pass either one bin bias (common to all bins) or %d bin biases",
                     ppt->selection_num,int1,ppt->selection_num);
        }
        free(pointer1);
      }

      class_call(parser_read_list_of_doubles(pfc,
                                             "selection_magnification_bias",
                                             &(int1),
                                             &(pointer1),
                                             &flag1,
                                             errmsg),
                 errmsg,
                 errmsg);

      if ((flag1 == _TRUE_) && (int1>0)) {

        if (int1==1) {
          for (i=0; i<ppt->selection_num; i++) {
            ptr->selection_magnification_bias[i] = pointer1[0];
          }
        }
        else if (int1==ppt->selection_num) {
          for (i=0; i<int1; i++) {
            ptr->selection_magnification_bias[i] = pointer1[i];
          }
        }
        else {
          class_stop(errmsg,
                     "In input for selection function, you asked for %d bin centers and %d bin biases; number of bins unclear; you should pass either one bin bias (common to all bins) or %d bin biases",
                     ppt->selection_num,int1,ppt->selection_num);
        }
        free(pointer1);
      }

    }

    if (ppt->selection_num>1) {
      class_read_int("non_diagonal",psp->non_diag);
      if ((psp->non_diag<0) || (psp->non_diag>=ppt->selection_num))
        class_stop(errmsg,
                   "Input for non_diagonal is %d, while it is expected to be between 0 and %d\n",
                   psp->non_diag,ppt->selection_num-1);
    }

    class_call(parser_read_string(pfc,
                                  "dNdz_selection",
                                  &(string1),
                                  &(flag1),
                                  errmsg),
               errmsg,
               errmsg);

    if ((flag1 == _TRUE_)) {
      if ((strstr(string1,"analytic") != NULL)){
        ptr->has_nz_analytic = _TRUE_;
      }
      else{
        ptr->has_nz_file = _TRUE_;
        class_read_string("dNdz_selection",ptr->nz_file_name);
      }
    }

    class_call(parser_read_string(pfc,
                                  "dNdz_evolution",
                                  &(string1),
                                  &(flag1),
                                  errmsg),
               errmsg,
               errmsg);

    if ((flag1 == _TRUE_)) {
      if ((strstr(string1,"analytic") != NULL)){
        ptr->has_nz_evo_analytic = _TRUE_;
      }
      else{
        ptr->has_nz_evo_file = _TRUE_;
        class_read_string("dNdz_evolution",ptr->nz_evo_file_name);
      }
    }

    flag1 = _FALSE_;
    class_call(parser_read_double(pfc,"bias",&param1,&flag1,errmsg),
               errmsg,
               errmsg);
    class_test(flag1 == _TRUE_,
               errmsg,
               "the input parameter 'bias' is obsolete, because you can now pass an independent light-to-mass bias for each bin/selection function. The new input name is 'selection_bias'. It can be set to a single number (common bias for all bins) or as many numbers as the number of bins");

    flag1 = _FALSE_;
    class_call(parser_read_double(pfc,"s_bias",&param1,&flag1,errmsg),
               errmsg,
               errmsg);
    class_test(flag1 == _TRUE_,
               errmsg,
               "the input parameter 's_bias' is obsolete, because you can now pass an independent magnitude bias for each bin/selection function. The new input name is 'selection_magnitude_bias'. It can be set to a single number (common magnitude bias for all bins) or as many numbers as the number of bins");

  }
  /* end of selection function section */

  /* deal with z_max issues */
  if ((ppt->has_pk_matter == _TRUE_) || (ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_) || (ppt->has_cl_number_count == _TRUE_) || (ppt->has_cl_lensing_potential == _TRUE_)) {

    class_call(parser_read_double(pfc,"z_max_pk",&param1,&flag1,errmsg),
               errmsg,
               errmsg);

    if (flag1==_TRUE_) {
      ppt->z_max_pk = param1;
    }
    else {
      // ppt->z_max_pk = 0.; //class_sz modif

      if ((ppt->has_pk_matter == _TRUE_) || (ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_)) {
        for (i=0; i<pop->z_pk_num; i++) {
          ppt->z_max_pk = MAX(ppt->z_max_pk,pop->z_pk[i]);
        }
      }

      if ((ppt->has_cl_number_count == _TRUE_) || (ppt->has_cl_lensing_potential == _TRUE_)) {

        for (bin=0; bin<ppt->selection_num; bin++) {

          /* the few lines below should be consistent with their counterpart in transfer.c, in transfer_selection_times() */
          if (ppt->selection==gaussian) {
            z_max = ppt->selection_mean[bin]+ppt->selection_width[bin]*ppr->selection_cut_at_sigma;
          }
          if (ppt->selection==tophat) {
            z_max = ppt->selection_mean[bin]+(1.+ppr->selection_cut_at_sigma*ppr->selection_tophat_edge)*ppt->selection_width[bin];
          }
          if (ppt->selection==dirac) {
            z_max = ppt->selection_mean[bin];
          }
          ppt->z_max_pk = MAX(ppt->z_max_pk,z_max);
        }
      }
    }
    psp->z_max_pk = ppt->z_max_pk;
  }

    if (pclass_sz->has_sz_ps
      + pclass_sz->has_sz_counts
      + pclass_sz->has_sz_counts_fft
      + pclass_sz->has_sz_rates
      + pclass_sz->has_hmf
      + pclass_sz->has_pk_at_z_1h
      + pclass_sz->has_pk_at_z_2h
      + pclass_sz->has_pk_gg_at_z_1h
      + pclass_sz->has_pk_gg_at_z_2h
      + pclass_sz->has_pk_bb_at_z_1h
      + pclass_sz->has_pk_bb_at_z_2h
      + pclass_sz->has_pk_b_at_z_2h
      + pclass_sz->has_gas_pressure_profile_2h
      + pclass_sz->has_gas_density_profile_2h
      + pclass_sz->has_pk_em_at_z_1h
      + pclass_sz->has_pk_em_at_z_2h
      + pclass_sz->has_pk_HI_at_z_1h
      + pclass_sz->has_pk_HI_at_z_2h
      + pclass_sz->has_bk_at_z_1h
      + pclass_sz->has_bk_at_z_2h
      + pclass_sz->has_bk_at_z_3h
      + pclass_sz->has_bk_ttg_at_z_1h
      + pclass_sz->has_bk_ttg_at_z_2h
      + pclass_sz->has_bk_ttg_at_z_3h
      + pclass_sz->has_bk_at_z_hf
      + pclass_sz->has_mean_y
      + pclass_sz->has_cib_monopole
      + pclass_sz->has_cib_shotnoise
      + pclass_sz->has_dcib0dz
      + pclass_sz->has_dydz
      + pclass_sz->has_sz_2halo
      + pclass_sz->has_sz_trispec
      + pclass_sz->has_sz_m_y_y_1h
      + pclass_sz->has_sz_m_y_y_2h
      + pclass_sz->has_sz_te_y_y
      + pclass_sz->has_sz_cov_N_N
      + pclass_sz->has_sz_cov_N_N_hsv
      + pclass_sz->has_tSZ_tSZ_tSZ_1halo
      + pclass_sz->has_tSZ_tSZ_tSZ_2h
      + pclass_sz->has_tSZ_tSZ_tSZ_3h
      + pclass_sz->has_kSZ_kSZ_1h
      + pclass_sz->has_kSZ_kSZ_2h
      + pclass_sz->has_kSZ_kSZ_tSZ_1h
      + pclass_sz->has_kSZ_kSZ_tSZ_2h
      + pclass_sz->has_kSZ_kSZ_tSZ_3h
      + pclass_sz->has_kSZ_kSZ_gal_1h
      + pclass_sz->has_kSZ_kSZ_gal_1h_fft
      + pclass_sz->has_kSZ_kSZ_gal_2h_fft
      + pclass_sz->has_kSZ_kSZ_gal_3h_fft
      + pclass_sz->has_kSZ_kSZ_gallens_1h_fft
      + pclass_sz->has_kSZ_kSZ_gallens_2h_fft
      + pclass_sz->has_kSZ_kSZ_gallens_3h_fft
      + pclass_sz->has_kSZ_kSZ_gallens_hf
      + pclass_sz->has_kSZ_kSZ_lens_1h_fft
      + pclass_sz->has_kSZ_kSZ_lens_2h_fft
      + pclass_sz->has_kSZ_kSZ_lens_3h_fft
      + pclass_sz->has_gal_gal_lens_1h_fft
      + pclass_sz->has_gal_gal_lens_2h_fft
      + pclass_sz->has_gal_gal_lens_3h_fft
      + pclass_sz->has_kSZ_kSZ_lens_hf
      + pclass_sz->has_gallens_gallens_1h
      + pclass_sz->has_gallens_gallens_2h
      + pclass_sz->has_gallens_lens_1h
      + pclass_sz->has_gallens_lens_2h
      + pclass_sz->has_kSZ_kSZ_gal_2h
      + pclass_sz->has_kSZ_kSZ_gal_3h
      + pclass_sz->has_kSZ_kSZ_gal_hf
      + pclass_sz->has_kSZ_kSZ_lensmag_1halo
      + pclass_sz->has_tSZ_gal_1h
      + pclass_sz->has_tSZ_gal_2h
      + pclass_sz->has_IA_gal_2h
      + pclass_sz->has_tSZ_cib_1h
      + pclass_sz->has_tSZ_cib_2h
      + pclass_sz->has_gallens_cib_1h
      + pclass_sz->has_gallens_cib_2h
      + pclass_sz->has_gal_cib_1h
      + pclass_sz->has_gal_cib_2h
      + pclass_sz->has_lens_cib_1h
      + pclass_sz->has_lens_cib_2h
      + pclass_sz->has_cib_cib_1h
      + pclass_sz->has_ngal_ngal_1h
      + pclass_sz->has_ngal_ngal_2h
      + pclass_sz->has_ngal_ngal_hf
      + pclass_sz->has_ngal_lens_1h
      + pclass_sz->has_ngal_lens_2h
      + pclass_sz->has_ngal_lens_hf
      + pclass_sz->has_ngal_nlensmag_hf
      + pclass_sz->has_ngal_gallens_1h
      + pclass_sz->has_ngal_gallens_2h
      + pclass_sz->has_ngal_IA_2h
      + pclass_sz->has_nlensmag_gallens_1h
      + pclass_sz->has_nlensmag_gallens_2h
      + pclass_sz->has_ngal_tsz_1h
      + pclass_sz->has_ngal_tsz_2h
      + pclass_sz->has_nlensmag_tsz_1h
      + pclass_sz->has_nlensmag_tsz_2h
      + pclass_sz->has_cib_cib_2h
      + pclass_sz->has_gal_gal_1h
      + pclass_sz->has_gal_gal_2h
      + pclass_sz->has_gal_gal_hf
      + pclass_sz->has_n5k
      + pclass_sz->has_tau_gal_1h
      + pclass_sz->has_tau_gal_2h
      + pclass_sz->has_tau_tau_1h
      + pclass_sz->has_tau_tau_2h
      + pclass_sz->has_gal_lens_1h
      + pclass_sz->has_gal_lens_2h
      + pclass_sz->has_gal_lens_hf
      + pclass_sz->has_lens_lens_hf
      + pclass_sz->has_gal_lensmag_1h
      + pclass_sz->has_gal_lensmag_2h
      + pclass_sz->has_gal_gallens_1h
      + pclass_sz->has_gal_gallens_2h
      + pclass_sz->has_tSZ_gallens_1h
      + pclass_sz->has_tSZ_gallens_2h
      + pclass_sz->has_gal_lensmag_hf
      + pclass_sz->has_tSZ_lensmag_1h
      + pclass_sz->has_tSZ_lensmag_2h
      + pclass_sz->has_gallens_lensmag_1h
      + pclass_sz->has_gallens_lensmag_2h
      + pclass_sz->has_lensmag_lensmag_1h
      + pclass_sz->has_lensmag_lensmag_2h
      + pclass_sz->has_lensmag_lensmag_hf
      + pclass_sz->has_lens_lensmag_1h
      + pclass_sz->has_lens_lensmag_2h
      + pclass_sz->has_lens_lensmag_hf
      + pclass_sz->has_lens_lens_1h
      + pclass_sz->has_lens_lens_2h
      + pclass_sz->has_tSZ_lens_1h
      + pclass_sz->has_tSZ_lens_2h
      + pclass_sz->has_isw_lens
      + pclass_sz->has_isw_tsz
      + pclass_sz->has_isw_auto
      + pclass_sz->has_vrms2
      + pclass_sz->has_dndlnM != _FALSE_){
        ppt->z_max_pk = pclass_sz->z2SZ;
        psp->z_max_pk = ppt->z_max_pk;
      }
  /* end of z_max section */

  class_call(parser_read_string(pfc,"root",&string1,&flag1,errmsg),
             errmsg,
             errmsg);
  if (flag1 == _TRUE_){
    class_test(strlen(string1)>_FILENAMESIZE_-32,errmsg,"Root directory name is too long. Please install in other directory, or increase _FILENAMESIZE_ in common.h");
    strcpy(pop->root,string1);
  }

  class_call(parser_read_string(pfc,
                                "headers",
                                &(string1),
                                &(flag1),
                                errmsg),
             errmsg,
             errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL))) {
    pop->write_header = _FALSE_;
  }

  class_call(parser_read_string(pfc,"format",&string1,&flag1,errmsg),
             errmsg,
             errmsg);

  if (flag1 == _TRUE_) {

    if ((strstr(string1,"class") != NULL) || (strstr(string1,"CLASS") != NULL))
      pop->output_format = class_format;
    else {
      if ((strstr(string1,"camb") != NULL) || (strstr(string1,"CAMB") != NULL))
        pop->output_format = camb_format;
      else
        class_stop(errmsg,
                   "You wrote: format='%s'. Could not identify any of the possible formats ('class', 'CLASS', 'camb', 'CAMB')",string1);
    }
  }

  /** (f) parameter related to the non-linear spectra computation */

  class_call(parser_read_string(pfc,
                                "non_linear",
                                &(string1),
                                &(flag1),
                                errmsg),
             errmsg,
             errmsg);

  if (flag1 == _TRUE_) {

    class_test(ppt->has_perturbations == _FALSE_, errmsg, "You requested non linear computation but no linear computation. You must set output to tCl or similar.");

    if ((strstr(string1,"halofit") != NULL) || (strstr(string1,"Halofit") != NULL) || (strstr(string1,"HALOFIT") != NULL)) {
      pnl->method=nl_halofit;
      ppt->k_max_for_pk = MAX(ppt->k_max_for_pk,MAX(ppr->halofit_min_k_max,ppr->nonlinear_min_k_max));
      ppt->has_nl_corrections_based_on_delta_m = _TRUE_;
    }
    if ((strstr(string1,"hmcode") != NULL) || (strstr(string1,"HMCODE") != NULL) || (strstr(string1,"HMcode") != NULL) || (strstr(string1,"Hmcode") != NULL)) {
      pnl->method=nl_HMcode;
      ppt->k_max_for_pk = MAX(ppt->k_max_for_pk,MAX(ppr->hmcode_min_k_max,ppr->nonlinear_min_k_max));
      ppt->has_nl_corrections_based_on_delta_m = _TRUE_;
      class_read_int("extrapolation_method",pnl->extrapolation_method);

      class_call(parser_read_string(pfc,
                                    "feedback model",
                                    &(string1),
                                    &(flag1),
                                    errmsg),
                 errmsg,
                 errmsg);

      if (flag1 == _TRUE_) {

		if (strstr(string1,"emu_dmonly") != NULL) {
          pnl->feedback = nl_emu_dmonly;
		}
		if (strstr(string1,"owls_dmonly") != NULL) {
          pnl->feedback = nl_owls_dmonly;
		}
		if (strstr(string1,"owls_ref") != NULL) {
          pnl->feedback = nl_owls_ref;
		}
		if (strstr(string1,"owls_agn") != NULL) {
          pnl->feedback = nl_owls_agn;
		}
		if (strstr(string1,"owls_dblim") != NULL) {
          pnl->feedback = nl_owls_dblim;
		}
      }

      class_call(parser_read_double(pfc,"eta_0",&param2,&flag2,errmsg),
                 errmsg,
                 errmsg);
      class_call(parser_read_double(pfc,"c_min",&param3,&flag3,errmsg),
                 errmsg,
                 errmsg);

      class_test(((flag1 == _TRUE_) && ((flag2 == _TRUE_) || (flag3 == _TRUE_))),
                 errmsg,
                 "In input file, you cannot enter both a baryonic feedback model and a choice of baryonic feedback parameters, choose one of both methods");

      if ((flag2 == _TRUE_) && (flag3 == _TRUE_)) {
		pnl->feedback = nl_user_defined;
		class_read_double("eta_0", pnl->eta_0);
		class_read_double("c_min", pnl->c_min);
      }
      else if ((flag2 == _TRUE_) && (flag3 == _FALSE_)) {
		pnl->feedback = nl_user_defined;
		class_read_double("eta_0", pnl->eta_0);
		pnl->c_min = (0.98 - pnl->eta_0)/0.12;
      }
      else if ((flag2 == _FALSE_) && (flag3 == _TRUE_)) {
		pnl->feedback = nl_user_defined;
		class_read_double("c_min", pnl->c_min);
		pnl->eta_0 = 0.98 - 0.12*pnl->c_min;
      }

      class_call(parser_read_double(pfc,"z_infinity",&param1,&flag1,errmsg),
                 errmsg,
                 errmsg);

      if (flag1 == _TRUE_) {
        class_read_double("z_infinity", pnl->z_infinity);
      }
    }
  }

  /** (g) amount of information sent to standard output (none if all set to zero) */

  class_read_int("background_verbose",
                 pba->background_verbose);

  class_read_int("thermodynamics_verbose",
                 pth->thermodynamics_verbose);

  class_read_int("perturbations_verbose",
                 ppt->perturbations_verbose);

  class_read_int("transfer_verbose",
                 ptr->transfer_verbose);

  class_read_int("primordial_verbose",
                 ppm->primordial_verbose);

  class_read_int("spectra_verbose",
                 psp->spectra_verbose);

  class_read_int("nonlinear_verbose",
                 pnl->nonlinear_verbose);

  class_read_int("lensing_verbose",
                 ple->lensing_verbose);

  class_read_int("output_verbose",
                 pop->output_verbose);


  if (ppt->has_tensors == _TRUE_) {
    /** - ---> Include ur and ncdm shear in tensor computation? */
    class_call(parser_read_string(pfc,"tensor method",&string1,&flag1,errmsg),
               errmsg,
               errmsg);
    if (flag1 == _TRUE_) {
      if (strstr(string1,"photons") != NULL)
        ppt->tensor_method = tm_photons_only;
      if (strstr(string1,"massless") != NULL)
        ppt->tensor_method = tm_massless_approximation;
      if (strstr(string1,"exact") != NULL)
        ppt->tensor_method = tm_exact;
    }
  }

  /** - ---> derivatives of baryon sound speed only computed if some non-minimal tight-coupling schemes is requested */
  if ((ppr->tight_coupling_approximation == (int)first_order_CLASS) || (ppr->tight_coupling_approximation == (int)second_order_CLASS)) {
    pth->compute_cb2_derivatives = _TRUE_;
  }

  class_test(ppr->ur_fluid_trigger_tau_over_tau_k==ppr->radiation_streaming_trigger_tau_over_tau_k,
             errmsg,
             "please choose different values for precision parameters ur_fluid_trigger_tau_over_tau_k and radiation_streaming_trigger_tau_over_tau_k, in order to avoid switching two approximation schemes at the same time");

  if (pba->N_ncdm>0) {

    class_test(ppr->ncdm_fluid_trigger_tau_over_tau_k==ppr->radiation_streaming_trigger_tau_over_tau_k,
               errmsg,
               "please choose different values for precision parameters ncdm_fluid_trigger_tau_over_tau_k and radiation_streaming_trigger_tau_over_tau_k, in order to avoid switching two approximation schemes at the same time");

    class_test(ppr->ncdm_fluid_trigger_tau_over_tau_k==ppr->ur_fluid_trigger_tau_over_tau_k,
               errmsg,
               "please choose different values for precision parameters ncdm_fluid_trigger_tau_over_tau_k and ur_fluid_trigger_tau_over_tau_k, in order to avoid switching two approximation schemes at the same time");

  }
  if (pba->Omega0_idr != 0.){
    class_test(ppr->idr_streaming_trigger_tau_over_tau_k==ppr->radiation_streaming_trigger_tau_over_tau_k,
               errmsg,
               "please choose different values for precision parameters dark_radiation_trigger_tau_over_tau_k and radiation_streaming_trigger_tau_over_tau_k, in order to avoid switching two approximation schemes at the same time");

    class_test(ppr->idr_streaming_trigger_tau_over_tau_k==ppr->ur_fluid_trigger_tau_over_tau_k,
               errmsg,
               "please choose different values for precision parameters dark_radiation_trigger_tau_over_tau_k and ur_fluid_trigger_tau_over_tau_k, in order to avoid switching two approximation schemes at the same time");

    class_test(ppr->idr_streaming_trigger_tau_over_tau_k==ppr->ncdm_fluid_trigger_tau_over_tau_k,
               errmsg,
               "please choose different values for precision parameters dark_radiation_trigger_tau_over_tau_k and ncdm_fluid_trigger_tau_over_tau_k, in order to avoid switching two approximation schemes at the same time");
  }


  /**
   * Here we can place all obsolete (deprecated) names for the precision parameters,
   * so they will still get read.
   * The new parameter names should be used preferrably
   * */
  class_read_double("k_scalar_min_tau0",ppr->k_min_tau0); // obsolete precision parameter: read for compatibility with old precision files
  class_read_double("k_scalar_max_tau0_over_l_max",ppr->k_max_tau0_over_l_max); // obsolete precision parameter: read for compatibility with old precision files
  class_read_double("k_scalar_step_sub",ppr->k_step_sub); // obsolete precision parameter: read for compatibility with old precision files
  class_read_double("k_scalar_step_super",ppr->k_step_super); // obsolete precision parameter: read for compatibility with old precision files
  class_read_double("k_scalar_step_transition",ppr->k_step_transition); // obsolete precision parameter: read for compatibility with old precision files
  class_read_double("k_scalar_k_per_decade_for_pk",ppr->k_per_decade_for_pk); // obsolete precision parameter: read for compatibility with old precision files
  class_read_double("k_scalar_k_per_decade_for_bao",ppr->k_per_decade_for_bao); // obsolete precision parameter: read for compatibility with old precision files
  class_read_double("k_scalar_bao_center",ppr->k_bao_center); // obsolete precision parameter: read for compatibility with old precision files
  class_read_double("k_scalar_bao_width",ppr->k_bao_width); // obsolete precision parameter: read for compatibility with old precision files

  class_read_double("k_step_trans_scalars",ppr->q_linstep); // obsolete precision parameter: read for compatibility with old precision files
  class_read_double("k_step_trans_tensors",ppr->q_linstep); // obsolete precision parameter: read for compatibility with old precision files
  class_read_double("k_step_trans",ppr->q_linstep); // obsolete precision parameter: read for compatibility with old precision files
  class_read_double("q_linstep_trans",ppr->q_linstep); // obsolete precision parameter: read for compatibility with old precision files
  class_read_double("q_logstep_trans",ppr->q_logstep_spline); // obsolete precision parameter: read for compatibility with old precision files

  class_call(parser_read_string(pfc,
                                "l_switch_limber_for_cl_density_over_z",
                                &string1,
                                &flag1,
                                errmsg),
             errmsg,
             errmsg);

  class_test(flag1 == _TRUE_,
             errmsg,
             "You passed in input a precision parameter called l_switch_limber_for_cl_density_over_z. This syntax is deprecated since v2.5.0. Please use instead the two precision parameters l_switch_limber_for_nc_local_over_z, l_switch_limber_for_nc_los_over_z, defined in include/common.h, and allowing for better performance.");

  /** (i) Write values in file */
  if (ple->has_lensed_cls == _TRUE_)
    ppt->l_scalar_max+=ppr->delta_l_max;


  // class_call(parser_read_string(pfc,"sBBN_file",&string1,&flag1,errmsg),
  //            errmsg,
  //            errmsg);
  // class_read_string("sBBN_file",ppr->sBBN_file);
  class_read_string("ksz_filter_file",pclass_sz->ksz_filter_file);
  class_read_string("projected_field_filter_file",pclass_sz->ksz_filter_file);
  class_read_string("full_path_to_dndz_gal",pclass_sz->full_path_to_dndz_gal);
  class_read_string("full_path_and_prefix_to_dndz_ngal",pclass_sz->full_path_and_prefix_to_dndz_ngal);
  class_read_string("full_path_to_redshift_dependent_M_min",pclass_sz->full_path_to_redshift_dependent_M_min);
  class_read_string("full_path_to_source_dndz_gal",pclass_sz->full_path_to_source_dndz_gal);
  // printf("-> File Name: %s\n",pclass_sz->ksz_filter_file);
    // printf("-> File Name: %s\n",pclass_sz->full_path_to_dndz_gal);
    // exit(0);
  class_read_string("cmb_lensing_noise_file",pclass_sz->cmb_lensing_noise_file);
  // printf("-> File Name: %s\n",pclass_sz->cmb_lensing_noise_file);
  // exit(0);
  class_read_string("A10_file",pclass_sz->A10_file);
  class_read_string("P13_file",pclass_sz->P13_file);
  // ppr->sBBN_file = string1;

  /** - (i.1.) shall we write background quantities in a file? */

  class_call(parser_read_string(pfc,"write background",&string1,&flag1,errmsg),
             errmsg,
             errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {

    pop->write_background = _TRUE_;

  }

  /** - (i.2.) shall we write thermodynamics quantities in a file? */

  class_call(parser_read_string(pfc,"write thermodynamics",&string1,&flag1,errmsg),
             errmsg,
             errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {

    pop->write_thermodynamics = _TRUE_;

  }

  /** - (i.3.) shall we write perturbation quantities in files? */

  class_call(parser_read_list_of_doubles(pfc,
                                         "k_output_values",
                                         &(int1),
                                         &(pointer1),
                                         &flag1,
                                         errmsg),
             errmsg,
             errmsg);

  if (flag1 == _TRUE_) {
    class_test(int1 > _MAX_NUMBER_OF_K_FILES_,
               errmsg,
               "you want to write some output for %d different values of k, hence you should increase _MAX_NUMBER_OF_K_FILES_ in include/perturbations.h to at least this number",
               int1);
    ppt->k_output_values_num = int1;

    for (i=0; i<int1; i++) {
      ppt->k_output_values[i] = pointer1[i];
    }
    free(pointer1);

    /* Sort the k_array using qsort */
    qsort (ppt->k_output_values, ppt->k_output_values_num, sizeof(double), compare_doubles);

    ppt->store_perturbations = _TRUE_;
    pop->write_perturbations = _TRUE_;
  }

  /** - (i.4.) shall we write primordial spectra in a file? */

  class_call(parser_read_string(pfc,"write primordial",&string1,&flag1,errmsg),
             errmsg,
             errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {

    pop->write_primordial = _TRUE_;

  }

  /** - (i.5) special steps if we want Halofit with wa_fld non-zero:
      so-called "Pk_equal method" of 0810.0190 and 1601.07230 */

  if (pnl->method == nl_halofit) {

    class_call(parser_read_string(pfc,"pk_eq",&string1,&flag1,errmsg),
               errmsg,
               errmsg);

    if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {

      if ((pba->Omega0_fld != 0.) && (pba->wa_fld != 0.)){

        pnl->has_pk_eq = _TRUE_;
      }
    }
  }

  return _SUCCESS_;

}

/**
 * All default parameter values (for input parameters)
 *
 * @param pba Input: pointer to background structure
 * @param pth Input: pointer to thermodynamics structure
 * @param ppt Input: pointer to perturbation structure
 * @param ptr Input: pointer to transfer structure
 * @param ppm Input: pointer to primordial structure
 * @param psp Input: pointer to spectra structure
 * @param pnl Input: pointer to nonlinear structure
 * @param ple Input: pointer to lensing structure
 * @param pop Input: pointer to output structure
 * @return the error status
 */

int input_default_params(
                         struct background *pba,
                         struct thermo *pth,
                         struct perturbs *ppt,
                         struct transfers *ptr,
                         struct primordial *ppm,
                         struct spectra *psp,
                         struct nonlinear * pnl,
                         struct lensing *ple,
                         struct class_sz_structure *pclass_sz, //BB: added for class_sz
                         struct szcount *pcsz, //BB: added for class_sz
                         struct output *pop
                         ) {

  double sigma_B; /* Stefan-Boltzmann constant in \f$ W/m^2/K^4 = Kg/K^4/s^3 \f$*/

  sigma_B = 2. * pow(_PI_,5) * pow(_k_B_,4) / 15. / pow(_h_P_,3) / pow(_c_,2);

  /** Define all default parameter values (for input parameters) for each structure:*/
  /** - background structure */

  /* 5.10.2014: default parameters matched to Planck 2013 + WP
     best-fitting model, with ones small difference: the published
     Planck 2013 + WP bestfit is with h=0.6704 and one massive
     neutrino species with m_ncdm=0.06eV; here we assume only massless
     neutrinos in the default model; for the CMB, taking m_ncdm = 0 or
     0.06 eV makes practically no difference, provided that we adapt
     the value of h in order ot get the same peak scale, i.e. the same
     100*theta_s. The Planck 2013 + WP best-fitting model with
     h=0.6704 gives 100*theta_s = 1.042143 (or equivalently
     100*theta_MC=1.04119). By taking only massless neutrinos, one
     gets the same 100*theta_s provided that h is increased to
     0.67556. Hence, we take h=0.67556, N_ur=3.046, N_ncdm=0, and all
     other parameters from the Planck2013 Cosmological Parameter
     paper. */

  pba->h = 0.67556;
  pba->H0 = pba->h * 1.e5 / _c_;
  pba->T_cmb = 2.7255;
  pba->Omega0_g = (4.*sigma_B/_c_*pow(pba->T_cmb,4.)) / (3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_);
  pba->Omega0_ur = 3.046*7./8.*pow(4./11.,4./3.)*pba->Omega0_g;
  pba->Omega0_idr = 0.0;
  pba->Omega0_idm_dr = 0.0;
  pba->T_idr = 0.0;
  pba->Omega0_b = 0.022032/pow(pba->h,2);
  pba->Omega0_cdm = 0.12038/pow(pba->h,2);
  pba->Omega0_dcdmdr = 0.0;
  pba->Omega0_dcdm = 0.0;
  pba->Gamma_dcdm = 0.0;
  pba->N_ncdm = 0;
  pba->Omega0_ncdm_tot = 0.;
  pba->ksi_ncdm_default = 0.;
  pba->ksi_ncdm = NULL;
  pba->T_ncdm_default = 0.71611; /* this value gives m/omega = 93.14 eV b*/
  pba->T_ncdm = NULL;
  pba->deg_ncdm_default = 1.;
  pba->deg_ncdm = NULL;
  pba->ncdm_psd_parameters = NULL;
  pba->ncdm_psd_files = NULL;

  pba->Omega0_scf = 0.; /* Scalar field defaults */
  pba->attractor_ic_scf = _TRUE_;
  pba->scf_parameters = NULL;
  pba->scf_parameters_size = 0;
  pba->scf_tuning_index = 0;
  //MZ: initial conditions are as multiplicative factors of the radiation attractor values
  pba->phi_ini_scf = 1;
  pba->phi_prime_ini_scf = 1;

  pba->Omega0_k = 0.;
  pba->K = 0.;
  pba->sgnK = 0;
  pba->Omega0_lambda = 1.-pba->Omega0_k-pba->Omega0_g-pba->Omega0_ur-pba->Omega0_b-pba->Omega0_cdm-pba->Omega0_ncdm_tot-pba->Omega0_dcdmdr-pba->Omega0_idm_dr-pba->Omega0_idr;
  pba->Omega0_fld = 0.;
  pba->a_today = 1.;
  pba->use_ppf = _TRUE_;
  pba->c_gamma_over_c_fld = 0.4;
  pba->fluid_equation_of_state = CLP;
  pba->w0_fld = -1.;
  pba->wa_fld = 0.;
  pba->Omega_EDE = 0.;
  pba->cs2_fld = 1.;

  pba->shooting_failed = _FALSE_;

  /** - thermodynamics structure */

  pth->YHe=_BBN_;
  pth->recombination=recfast;
  pth->reio_parametrization=reio_camb;
  pth->reio_z_or_tau=reio_z;
  pth->z_reio=11.357;
  pth->tau_reio=0.0925;
  pth->reionization_exponent=1.5;
  pth->reionization_width=0.5;
  pth->helium_fullreio_redshift=3.5;
  pth->helium_fullreio_width=0.5;

  pth->binned_reio_num=0;
  pth->binned_reio_z=NULL;
  pth->binned_reio_xe=NULL;
  pth->binned_reio_step_sharpness = 0.3;

  pth->annihilation = 0.;
  pth->decay = 0.;

  pth->annihilation_variation = 0.;
  pth->annihilation_z = 1000.;
  pth->annihilation_zmax = 2500.;
  pth->annihilation_zmin = 30.;
  pth->annihilation_f_halo = 0.;
  pth->annihilation_z_halo = 30.;
  pth->has_on_the_spot = _TRUE_;

  pth->compute_cb2_derivatives=_FALSE_;

  pth->compute_damping_scale = _FALSE_;

  pth->a_idm_dr = 0.;
  pth->b_idr = 0.;
  pth->nindex_idm_dr = 4.;
  pth->m_idm = 1.e11;

  /** - perturbation structure */

  ppt->has_cl_cmb_temperature = _FALSE_;
  ppt->has_cl_cmb_polarization = _FALSE_;
  ppt->has_cl_cmb_lensing_potential = _FALSE_;
  ppt->has_cl_number_count = _FALSE_;
  ppt->has_cl_lensing_potential = _FALSE_;
  ppt->has_pk_matter = _FALSE_;
  ppt->has_density_transfers = _FALSE_;
  ppt->has_velocity_transfers = _FALSE_;
  ppt->has_metricpotential_transfers = _FALSE_;

  ppt->has_nl_corrections_based_on_delta_m = _FALSE_;

  ppt->has_nc_density = _FALSE_;
  ppt->has_nc_rsd = _FALSE_;
  ppt->has_nc_lens = _FALSE_;
  ppt->has_nc_gr = _FALSE_;

  //ppt->pk_only_cdm_bar=_FALSE_;

  ppt->switch_sw = 1;
  ppt->switch_eisw = 1;
  ppt->switch_lisw = 1;
  ppt->switch_dop = 1;
  ppt->switch_pol = 1;
  ppt->eisw_lisw_split_z = 120;

  ppt->has_ad=_TRUE_;
  ppt->has_bi=_FALSE_;
  ppt->has_cdi=_FALSE_;
  ppt->has_nid=_FALSE_;
  ppt->has_niv=_FALSE_;

  ppt->has_perturbed_recombination=_FALSE_;
  ppt->tensor_method = tm_massless_approximation;
  ppt->evolve_tensor_ur = _FALSE_;
  ppt->evolve_tensor_ncdm = _FALSE_;

  ppt->has_scalars=_TRUE_;
  ppt->has_vectors=_FALSE_;
  ppt->has_tensors=_FALSE_;

  ppt->l_scalar_max=2500;
  ppt->l_vector_max=500;
  ppt->l_tensor_max=500;
  ppt->l_lss_max=300;
  ppt->k_max_for_pk=60.; //BB: changed in class_sz; was set to 1 in original class version, now 10. [TBC]

  ppt->gauge=synchronous;

  ppt->idr_nature=idr_free_streaming;

  ppt->has_Nbody_gauge_transfers = _FALSE_;

  ppt->k_output_values_num=0;
  ppt->store_perturbations = _FALSE_;

  ppt->three_ceff2_ur=1.;
  ppt->three_cvis2_ur=1.;

  ppt->z_max_pk=0.;

  ppt->selection_num=1;
  ppt->selection=gaussian;
  ppt->selection_mean[0]=1.;
  ppt->selection_width[0]=0.1;

  /** - primordial structure */

  ppm->primordial_spec_type = analytic_Pk;
  ppm->k_pivot = 0.05;
  ppm->A_s = 2.215e-9;
  ppm->n_s = 0.9619;
  ppm->alpha_s = 0.;
  ppm->f_bi = 1.;
  ppm->n_bi = 1.;
  ppm->alpha_bi = 0.;
  ppm->f_cdi = 1.;
  ppm->n_cdi = 1.;
  ppm->alpha_cdi = 0.;
  ppm->f_nid = 1.;
  ppm->n_nid = 1.;
  ppm->alpha_nid = 0.;
  ppm->f_niv = 1.;
  ppm->n_niv = 1.;
  ppm->alpha_niv = 0.;
  ppm->c_ad_bi = 0.;
  ppm->n_ad_bi = 0.;
  ppm->alpha_ad_bi = 0.;
  ppm->c_ad_cdi = 0.;
  ppm->n_ad_cdi = 0.;
  ppm->alpha_ad_cdi = 0.;
  ppm->c_ad_nid = 0.;
  ppm->n_ad_nid = 0.;
  ppm->alpha_ad_nid = 0.;
  ppm->c_ad_niv = 0.;
  ppm->n_ad_niv = 0.;
  ppm->alpha_ad_niv = 0.;
  ppm->c_bi_cdi = 0.;
  ppm->n_bi_cdi = 0.;
  ppm->alpha_bi_cdi = 0.;
  ppm->c_bi_nid = 0.;
  ppm->n_bi_nid = 0.;
  ppm->alpha_bi_nid = 0.;
  ppm->c_bi_niv = 0.;
  ppm->n_bi_niv = 0.;
  ppm->alpha_bi_niv = 0.;
  ppm->c_cdi_nid = 0.;
  ppm->n_cdi_nid = 0.;
  ppm->alpha_cdi_nid = 0.;
  ppm->c_cdi_niv = 0.;
  ppm->n_cdi_niv = 0.;
  ppm->alpha_cdi_niv = 0.;
  ppm->c_nid_niv = 0.;
  ppm->n_nid_niv = 0.;
  ppm->alpha_nid_niv = 0.;
  ppm->r = 1.;
  ppm->n_t = -ppm->r/8.*(2.-ppm->r/8.-ppm->n_s);
  ppm->alpha_t = ppm->r/8.*(ppm->r/8.+ppm->n_s-1.);
  ppm->potential=polynomial;
  ppm->phi_end=0.;
  ppm->phi_pivot_method = N_star;
  ppm->phi_pivot_target = 60;
  ppm->V0=1.25e-13;
  ppm->V1=-1.12e-14;
  ppm->V2=-6.95e-14;
  ppm->V3=0.;
  ppm->V4=0.;
  ppm->H0=3.69e-6;
  ppm->H1=-5.84e-7;
  ppm->H2=0.;
  ppm->H3=0.;
  ppm->H4=0.;
  ppm->behavior=numerical;
  ppm->command="write here your command for the external Pk";
  ppm->custom1=0.;
  ppm->custom2=0.;
  ppm->custom3=0.;
  ppm->custom4=0.;
  ppm->custom5=0.;
  ppm->custom6=0.;
  ppm->custom7=0.;
  ppm->custom8=0.;
  ppm->custom9=0.;
  ppm->custom10=0.;

  /** - nonlinear structure */

  pnl->method = nl_none;
  pnl->extrapolation_method = extrap_max_scaled;
  pnl->has_pk_eq = _FALSE_;

  pnl->feedback = nl_emu_dmonly;
  pnl->z_infinity = 10.;

  /** - transfer structure */

  ptr->selection_bias[0]=1.;
  ptr->selection_magnification_bias[0]=0.;
  ptr->lcmb_rescale=1.;
  ptr->lcmb_pivot=0.1;
  ptr->lcmb_tilt=0.;
  ptr->initialise_HIS_cache=_FALSE_;
  ptr->has_nz_analytic = _FALSE_;
  ptr->has_nz_file = _FALSE_;
  ptr->has_nz_evo_analytic = _FALSE_;
  ptr->has_nz_evo_file = _FALSE_;

  /** - spectra structure */

  psp->z_max_pk = pop->z_pk[0];
  psp->non_diag=0;

  /** - lensing structure */

  ple->has_lensed_cls = _FALSE_;

  /** - output structure */

  pop->z_pk_num = 1;
  pop->z_pk[0] = 0.;
  sprintf(pop->root,"output/");
  pop->write_header = _TRUE_;
  pop->output_format = class_format;
  pop->write_background = _FALSE_;
  pop->write_thermodynamics = _FALSE_;
  pop->write_perturbations = _FALSE_;
  pop->write_primordial = _FALSE_;

  /** - all verbose parameters */

  pba->background_verbose = 0;
  pth->thermodynamics_verbose = 0;
  ppt->perturbations_verbose = 0;
  ptr->transfer_verbose = 0;
  ppm->primordial_verbose = 0;
  psp->spectra_verbose = 0;
  pnl->nonlinear_verbose = 0;
  ple->lensing_verbose = 0;
  pop->output_verbose = 0;



  psp->overwrite_clpp_with_limber = 0;
  //BB: SZ parameters default values



  pclass_sz->write_sz = _FALSE_;
  pclass_sz->ell_sz = 4;
  pclass_sz->dlogell = .3;
  pclass_sz->dell = 0.;
  pclass_sz->ell_max_mock = 4000.;
  pclass_sz->ell_min_mock = 100.;

  pclass_sz->dlogfreq = .3;
  pclass_sz->dfreq = 0.;
  pclass_sz->freq_max = 1000.;
  pclass_sz->freq_min = 100.;
  //pclass_sz->nlSZ = 18;

  pclass_sz->bispectrum_lambda_k2 = 1.;
  pclass_sz->bispectrum_lambda_k3 = 1.;

  pclass_sz->M0_Mmin_flag = 0;
  pclass_sz->f_cen_HOD = 1.;
  pclass_sz->Delta_z_lens = 0.; // DES photo-z errors
  pclass_sz->Delta_z_source = 0.;
  pclass_sz->photo_z_params = 0.; //
  pclass_sz->dndz_shift_gal = 0; // shift and stretch params
  pclass_sz->dndz_shift_source_gal = 0; //as in https://arxiv.org/pdf/2210.08633.pdf
  pclass_sz->dndz_stretch_gal = 1.;
  pclass_sz->dndz_stretch_source_gal = 1.;
  pclass_sz->shear_calibration_m = 0.;
  pclass_sz->cosmo_model = 0; // 0 lcdm, 1 mnu, 2 neff, 3 wcdm, 4 ede, 5 mnu-3states, 6 ede-3states
  pclass_sz->use_Amod = 0;
  pclass_sz->Amod = 0.;
  pclass_sz->M_min_HOD_mass_factor_unwise = 1.;
  pclass_sz->M0_HOD = 0.; //DES-like HOD see https://arxiv.org/pdf/2106.08438.pdf

  pclass_sz->x_out_truncated_density_profile = 1.; // the numerical one.
  pclass_sz->x_out_truncated_nfw_profile = 1.; //the analytical one
  pclass_sz->x_out_truncated_nfw_profile_electrons = 1.; //the analytical one
  pclass_sz->x_out_truncated_nfw_profile_satellite_galaxies =1.;

  pclass_sz->x_out_custom1 = 1.;
  // pclass_sz->x_out_nfw_profile = 2.5;
  pclass_sz->cvir_tau_profile_factor =  1.;
  pclass_sz->M1_prime_HOD_factor = 15.;
  pclass_sz->M_min_HOD_satellite_mass_factor_unwise = 0.1;

  pclass_sz->hm_consistency = 0; //0: nothing 1: counter terms 2: alpha(z)
  pclass_sz->check_consistency_conditions = 0;
  pclass_sz->damping_1h_term = 1;

  pclass_sz->use_class_sz_fast_mode = 0;

  pclass_sz->N_kSZ2_gal_multipole_grid = 20;

  pclass_sz->concentration_parameter=6;
  //pclass_sz->hod_model=-1;
  pclass_sz->effective_temperature=0;
  pclass_sz->create_ref_trispectrum_for_cobaya=0;


  pclass_sz->effective_galaxy_bias = 1.;
  pclass_sz->use_bg_eff_in_ksz2g_eff = 0;
  pclass_sz->compute_ksz2ksz2  = 0;
  // sprintf(pclass_sz->path_to_ref_trispectrum_for_cobaya,"output/");
  // sprintf(pclass_sz->append_name_cobaya_ref,"for_cobaya");
  // sprintf(pclass_sz->full_path_to_noise_curve_for_y_y,"sz_auxiliary_files/my_noise_curve_yxy.txt");

  pclass_sz->nl_yy_is_binned = 0;

  pcsz->redshift_for_dndm = 1.e-5;





  pclass_sz->Tcmb_gNU = pba->T_cmb*((_h_P_*150.0e9/(_k_B_*pba->T_cmb))*(1./tanh((_h_P_*150.0e9/(_k_B_*pba->T_cmb))/2.))-4.);




  pclass_sz->has_completeness_for_ps_SZ = 0;

  //default: total power spectrum (no completeness cut)
  pclass_sz->which_ps_sz = 0; //0: total, 1: resolved, 2: unresolved

  pclass_sz->use_analytical_truncated_nfw = 1;

  pclass_sz->use_m500c_in_ym_relation = 1;
  pclass_sz->use_m200c_in_ym_relation = 0;
  pclass_sz->use_hod = 1;
  pclass_sz->galaxy_sample = 2; // WIxSC
  pclass_sz->unwise_galaxy_sample_id = -1; // red
  //pclass_sz->unwise_m_min_cut = 1e10; // Msun/h

  pclass_sz->sn_cutoff = 5.;
  pcsz->sn_cutoff = 5.;
  //Redshift limits for the integration
  pclass_sz->z1SZ = 1.e-4;
  pclass_sz->z2SZ = 6.;
  ppt->z_max_pk = pclass_sz->z2SZ;
  psp->z_max_pk = ppt->z_max_pk;

  pclass_sz->z1SZ_dndlnM = 5.e-3;
  pclass_sz->z2SZ_dndlnM = 6.;
  pclass_sz->N_redshift_dndlnM = 50;

  pclass_sz->M1SZ_dndlnM = 1.e8;
  pclass_sz->M2SZ_dndlnM = 1.e17;
  pclass_sz->N_mass_dndlnM = 200;

  pclass_sz->kstar_damping_1h_term_Mpc = 0.01; //same as hmvec default value in inverse Mpc

  pclass_sz->ndim_redshifts_for_integral = 30; //used in redshift integral



  // pclass_sz->n_m_matter_density_profile =  100;

  pclass_sz->include_y_counterterms_in_yk = 1;
  pclass_sz->include_g_counterterms_in_gk = 1;
  pclass_sz->include_k_counterterms_in_gk = 1;
  pclass_sz->include_gk_counterterms_in_gk = 1;
  //mass limits: h^-1 Msun
  pclass_sz->M1SZ = 1.e10;
  pclass_sz->M2SZ = 5.e15;

  pclass_sz->n_z_dndlnM = 500;
  pclass_sz->n_m_dndlnM = 500;


  // pclass_sz->M1SZ_L_sat = 1.e9;
  // pclass_sz->M2SZ_L_sat = 1.e17;
  // pclass_sz->z1SZ_L_sat = 1.e-3;
  // pclass_sz->z2SZ_L_sat = 6.;
  // pclass_sz->n_z_L_sat = 101;
  // pclass_sz->n_m_L_sat = 102;
  // pclass_sz->n_nu_L_sat = 103;
  // pclass_sz->epsabs_L_sat = 1e-15;
  // pclass_sz->epsrel_L_sat = 1e-6;

  pclass_sz->convert_cls_to_gamma = 0;




  pclass_sz->n_z_W_lensmag = 500;
  pclass_sz->n_z_W_gallens_sources = 500;

  //Set pressure profile to A10
  pclass_sz->pressure_profile=2;


  pclass_sz->tau_profile = 0; // scaled nfw
  pclass_sz->tau_profile_mode = 0; // agn feedback

  //Pressure profile is considered between x_in and x_out
  //See Komatsu
  pclass_sz->x_inSZ = 1.e-5; //KS02
  pclass_sz->x_outSZ = 6.; //KS02
  pclass_sz->ln_x_size_for_pp = 1000;
  pclass_sz->x_size_for_pp = 2000;

  pclass_sz->f_sky = 1.; // full sky
  pclass_sz->fsky_from_skyfracs = 1.;
  pclass_sz->szunbinned_loglike = -1000.;
  pclass_sz->Omega_survey = 4.*_PI_*pclass_sz->f_sky;

  // integration_method_pressure_profile
  // 1: 0 -> Patterson rule
  // 2: 1 -> GSL
  // 3: 2 -> spline integral
  pclass_sz->integration_method_pressure_profile = 1;

  //P13 UPP parameters
  // pclass_sz->P0GNFW = 6.41;
  // pclass_sz->c500 = 1.81;
  // pclass_sz->gammaGNFW = 0.31;
  // pclass_sz->alphaGNFW = 1.33;
  // pclass_sz->betaGNFW = 4.13;

  pclass_sz->A_IA = 0.5; // see https://arxiv.org/pdf/2106.08438.pdf
  pclass_sz->eta_IA = -1.0; // see https://arxiv.org/pdf/2106.08438.pdf
  pclass_sz->C1_IA =  5e-14;
  //A10 UPP parameters
  pclass_sz->P0GNFW = 8.130;
  pclass_sz->c500 = 1.156;
  pclass_sz->gammaGNFW = 0.3292;
  pclass_sz->alphaGNFW = 1.0620;
  pclass_sz->betaGNFW = 5.4807;

  pclass_sz->delta_alpha = 0.;
  pclass_sz->alpha_p = 0.12;
  //Hydrostatic Equilibrium Mass Bias, Piffaretti & Valdarnini [arXiv:0808.1111]

   pclass_sz->truncate_gas_pressure_wrt_rvir = 0;
   pclass_sz->no_tt_noise_in_kSZ2X_cov = 0;
  // battaglia pressure profile:
  pclass_sz->gamma_B12 = -0.3;
  pclass_sz->alpha_B12 = 1.;
   pclass_sz->P0_B12 = 18.1;
   pclass_sz->xc_B12 = 0.497;
   pclass_sz->beta_B12 = 4.35;

   pclass_sz->alpha_m_P0_B12 = 0.154;
   pclass_sz->alpha_m_xc_B12 = -0.00865;
   pclass_sz->alpha_m_beta_B12 = 0.0393;

   pclass_sz->alpha_z_P0_B12 = -0.758;
   pclass_sz->alpha_z_xc_B12 = 0.731;
   pclass_sz->alpha_z_beta_B12 = 0.415;

   pclass_sz->c_B12 = 0.;
   pclass_sz->mcut_B12 = 1.e14;
   // repeated??
   pclass_sz->alphap_m_P0_B12 = 0.154;
   pclass_sz->alphap_m_xc_B12 = -0.00865;
   pclass_sz->alphap_m_beta_B12 = 0.0393;

   pclass_sz->alpha_c_P0_B12 = 0.;
   pclass_sz->alpha_c_xc_B12 = 0.;
   pclass_sz->alpha_c_beta_B12 = 0.;



   pclass_sz->use_websky_m200m_to_m200c_conversion = 0;


   //battaglia density profile
   // default set to agn paremeters
   pclass_sz->A_rho0 = 4.e3;
   pclass_sz->A_alpha = 0.88;
   pclass_sz->A_beta = 3.83;

   pclass_sz->alpha_m_rho0 = 0.29;
   pclass_sz->alpha_m_alpha = -0.03;
   pclass_sz->alpha_m_beta = 0.04;

   pclass_sz->alpha_z_rho0 = -0.66;
   pclass_sz->alpha_z_alpha = 0.19;
   pclass_sz->alpha_z_beta = -0.025;
   pclass_sz->gamma_B16 = -0.2;
   pclass_sz->xc_B16 = 0.5;


   pclass_sz->c_B16 = 0.;
	 pclass_sz->mcut = 1.e14;
   pclass_sz->alphap_m_rho0 = 0.29;
   pclass_sz->alphap_m_alpha = -0.03;
   pclass_sz->alphap_m_beta = 0.04;

   pclass_sz->alpha_c_rho0 = 0.;
   pclass_sz->alpha_c_alpha = 0.;
   pclass_sz->alpha_c_beta = 0.;
  //units
  pclass_sz->nu_y_dist_GHz = 150.;
  pclass_sz->exponent_unit = 2; // 2: dimensionless, 0:  'muK' (micro Kelvin)

  pclass_sz->id_nu_cib_to_save = 0;
  pclass_sz->id_nu_prime_cib_to_save = 0;

  // Table 1  of MM20
  pclass_sz->alpha_cib = 0.36; //redshift evolution of dust temperature
  pclass_sz->T0_cib = 24.4; // dust temperature today in Kelvins
  pclass_sz->beta_cib = 1.75; // emissivity index of sed
  pclass_sz->gamma_cib = 1.7; // Power law index of SED at high frequency
  pclass_sz->delta_cib = 3.6; // Redshift evolution of L − M normalisation
  pclass_sz->m_eff_cib = pow(10.,12.6); // Most efficient halo mass in Msun
  pclass_sz->L0_cib = 6.4e-8; // Normalisation of L − M relation in [Jy MPc2/Msun]
  pclass_sz->sigma2_LM_cib = 0.5; // Size of of halo masses sourcing CIB emission
  pclass_sz->z_obs_cib = 0.;
  pclass_sz->z_plateau_cib = 1e100; // see 5.2.1 of https://arxiv.org/pdf/1208.5049.pdf
  pclass_sz->M_min_subhalo_in_Msun = 0;
  pclass_sz->use_redshift_dependent_M_min = 0;
  pclass_sz->cib_nu0_norm = 1;
  pclass_sz->use_nc_1_for_all_halos_cib_HOD = 0;
  pclass_sz->maniyar_cib_etamax = 0.4028353504978569; // see https://github.com/abhimaniyar/halomodel_cib_tsz_cibxtsz/blob/master/input_var.py
  pclass_sz->maniyar_cib_zc = 1.5; // see https://github.com/abhimaniyar/halomodel_cib_tsz_cibxtsz/blob/master/input_var.py
  pclass_sz->maniyar_cib_tau = 1.2040244128818796; // see https://github.com/abhimaniyar/halomodel_cib_tsz_cibxtsz/blob/master/input_var.py
  pclass_sz->maniyar_cib_fsub = 0.134; // see https://github.com/abhimaniyar/halomodel_cib_tsz_cibxtsz/blob/master/Cell_cib.py
  pclass_sz->fNL = 0.;


  // ## fiducial BCM parameters (matches BAHAMAS)
  // logMc = 13.25
  // theta_ej = 4.711
  // eta_star = 0.2
  // eta_cga = 0.297 /// not needed.
  // mu = 1.0
  // gamma = 2.5
  // delta = 7.0
  // nu_Mc = 0.038
  pclass_sz->log10Mc_bcm = 13.25;
  pclass_sz->theta_ej_bcm = 4.711;
  pclass_sz->eta_star_bcm = 0.2;
  pclass_sz->delta_bcm = 7.0;
  pclass_sz->gamma_bcm = 2.5;
  pclass_sz->mu_bcm = 1.0;
  pclass_sz->nu_log10Mc_bcm = -0.038;


  //# Table 1 of https://arxiv.org/pdf/1309.0382.pdf
  pclass_sz->has_cib_flux_cut  = 0;
  pclass_sz->cib_Snu_cutoff_list_in_mJy = NULL;


  pclass_sz->cib_frequency_list_num=1;
  pclass_sz->cib_frequency_list=NULL;


  pclass_sz->HSEbias = 1.;
  pclass_sz->Ap = 0.1;
  pclass_sz->alpha_b = 0.7;
  pclass_sz->mass_dependent_bias = 0;

  pclass_sz->bin_z_min_cluster_counts = 0.;
  pclass_sz->bin_z_max_cluster_counts = 1.;
  pclass_sz->bin_dz_cluster_counts = 0.1;

  pclass_sz->use_planck_binned_proba = 0;


  pclass_sz->bin_dlog10_snr = 0.25;
  pclass_sz->log10_snr_min = 0.6;
  pclass_sz->log10_snr_max = 2.0;


  pclass_sz->lnymin = -11.5;// -11.5 in planck
  pclass_sz->lnymax = 10.; // 10 in planck
  pclass_sz->dlny = 0.05; // 0.05 in planck
  pcsz->has_completeness = 0;
  pcsz->mass_range = 1;//szcount masses
  pclass_sz->experiment = 0; //planck
  pclass_sz->sky_area_deg2 = 599.;
  pclass_sz->apply_relativistic_correction_to_y_m = 0;

  pclass_sz->y_m_relation = 1; //0: planck, 1: act/so


  pclass_sz->use_maniyar_cib_model = 0;


  // pcsz->ystar = -0.186; //-0.186 ref. value in SZ_plus_priors.ini (cosmomc)
  // pcsz->alpha = 1.789; //1.789 ref. value in SZ_plus_priors.ini (cosmomc)
  // pcsz->ystar = pow(10.,pcsz->ystar)/pow(2., pcsz->alpha)*0.00472724; ////8.9138435358806980e-004;
  // pcsz->sigmaM = 0.075; //in log10 see tab 1 of planck cc 2015 paper
  // pcsz->beta = 0.66;
  // pcsz->thetastar = 6.997;
  // pcsz->alpha_theta = 1./3.;

  pclass_sz->ystar_ym = -0.186; //-0.186 ref. value in SZ_plus_priors.ini (cosmomc)
  // pclass_sz->alpha_ym = 1.78; //1.789 ref. value in SZ_plus_priors.ini (cosmomc)
  // pclass_sz->ystar_ym = pow(10.,-0.19)/pow(2., pclass_sz->alpha_ym)*0.00472724; ////8.9138435358806980e-004;
  pclass_sz->alpha_ym = 5./3.+0.12;//
  pclass_sz->sigmaM_ym = 0.173; //in natural log, see tab 1 of planck cc 2015 paper, intrinsic scatter
  //0.173/np.log(10.) = 0.075, in szcounts.f90. But it's wrong. It should be in natural log.

  pclass_sz->beta_ym = 2./3.;
  pclass_sz->thetastar = 6.997;
  pclass_sz->alpha_theta = 1./3.;

  // values in Hasselfield et al 2013
  pclass_sz->B_ym = 0.08;
  pclass_sz->A_ym = 4.95e-5;
  pclass_sz->C_ym = -0.025;

  pclass_sz->m_pivot_ym = 3e14;

  pclass_sz->temperature_mass_relation=0;

  //For the computation of sigma2

  //Array size
  pclass_sz->ndim_redshifts = 100;//number of z in the sigma Interpolation

  pclass_sz->ndim_masses = 100;
  pclass_sz->logR1SZ = -10; // 0.0034Mpc/h, 1.8e4  solar mass
  pclass_sz->logR2SZ = 10.; //default =4 , i.e., 54.9Mpc/h, 7.5e16 solar mass

  pclass_sz->delta_cSZ = (3./20.)*pow(12.*_PI_,2./3.); // this is = 1.686470199841145


  pclass_sz->k_per_decade_for_tSZ = 128.; //#default 40
  pclass_sz->k_min_for_pk_in_tSZ = 1.e-4;
  pclass_sz->k_max_for_pk_in_tSZ = 5.e1;

  pclass_sz->z_for_pk_hm = 1.;
  pclass_sz->k_min_for_pk_hm = 1e-4;
  pclass_sz->k_max_for_pk_hm = 1e2;
  pclass_sz->dlnk_for_pk_hm = 0.1;


  //Multplicity function Tinker 2010
  // https://arxiv.org/pdf/1001.3162.pdf
  // this is Table 4 of the T10 paper
  // for Delta = 200
  pclass_sz->T10_alpha_fixed = 0; // allpha is computed at each z.
  pclass_sz->alphaSZ = 0.368;
  pclass_sz->beta0SZ = 0.589;
  pclass_sz->gamma0SZ = 0.864;

  pclass_sz->phi0SZ = -0.729;
  pclass_sz->eta0SZ = -0.243;

  pclass_sz->no_b2 = 0;

  //Multplicity function Bocquet 2015 "Hydro"

  pclass_sz->Ap0 = 0.228;
  pclass_sz->a0 = 2.15;
  pclass_sz->b0 = 1.69;
  pclass_sz->c0 = 1.30;

  pclass_sz->pk_nonlinear_for_vrms2 = 0;

  pclass_sz->MF = 8; //Tinker et al 2008 @ M200c
  pclass_sz->SHMF = 1;

  //////////////////////////////////
  //Integration method and parameters (mass)
  //patterson
  pclass_sz->integration_method_mass = 0;
  pclass_sz->patterson_show_neval = 0;


  pclass_sz->redshift_epsrel = 1e-6;
  pclass_sz->redshift_epsabs = 1e-40;

  pclass_sz->mass_epsrel_cluster_counts = 1e-5;
  pclass_sz->mass_epsabs_cluster_counts = 1e-40;

  pclass_sz->redshift_epsrel_cluster_counts = 1e-5;
  pclass_sz->redshift_epsabs_cluster_counts = 1e-40;


  pclass_sz->dlnM_cluster_count_completeness_grid = 0.05;
  pclass_sz->dz_cluster_count_completeness_grid_low_z = 1.e-3;
  pclass_sz->dz_cluster_count_completeness_grid_mid_z = 1.e-2;
  pclass_sz->dz_cluster_count_completeness_grid_high_z = 1.e-1;

  pclass_sz->cluster_count_completeness_grid_z_cutoff_low = 0.2;
  pclass_sz->cluster_count_completeness_grid_z_cutoff_mid = 1.;



  pclass_sz->mass_epsrel = 1e-6;
  pclass_sz->mass_epsabs = 1e-40;

  pclass_sz->pressure_profile_epsrel = 1e-2;
  pclass_sz->pressure_profile_epsabs = 1e-4;

  pclass_sz->nfw_profile_epsrel = 1e-9;
  pclass_sz->nfw_profile_epsabs = 1e-10;
  //trapezoidal
  pclass_sz->number_of_mass_bins = 60;
  /////////////////////////////////

  //number of mass bins for cov(Y,N)
  pclass_sz->nbins_M = 20;

  pclass_sz->JMAX = 30; //for the pressure profile
  pclass_sz->EPS = 1.e-5; //for the pressure profile (default 1e-5)
  pclass_sz->K = 5; //for the pressure profile (default 5)


  pclass_sz->JMAX_sigma = 20; //not used
  pclass_sz->EPS_sigma = 1.e-6; //not used
  pclass_sz->K_sigma = 5; //not used


  //Nuisance

  pclass_sz->A_cib = 0.29;
  pclass_sz->A_rs = 0.01;
  pclass_sz->A_ir = 1.97;
  pclass_sz->A_cn = 1.0;

  pclass_sz->bispec_conf_id = 0;

  //pclass_sz->has_class_sz_structure = _FALSE_;
  pclass_sz->has_sz_counts = _FALSE_;
  pclass_sz->has_sz_counts_fft = _FALSE_;
  pclass_sz->has_sz_rates = _FALSE_;
  pclass_sz->has_isw_lens = _FALSE_;
  pclass_sz->has_isw_tsz = _FALSE_;
  pclass_sz->has_isw_auto = _FALSE_;
  pclass_sz->has_dndlnM = _FALSE_;
  pclass_sz->has_gal_gal_1h = _FALSE_;
  pclass_sz->has_gal_gal_2h = _FALSE_;
  pclass_sz->has_gal_gal_hf = _FALSE_;
  pclass_sz->has_n5k = _FALSE_;
  pclass_sz->has_tau_gal_1h = _FALSE_;
  pclass_sz->has_tau_gal_2h = _FALSE_;
  pclass_sz->has_tau_tau_1h = _FALSE_;
  pclass_sz->has_tau_tau_2h = _FALSE_;
  pclass_sz->has_gal_lens_1h = _FALSE_;
  pclass_sz->has_gal_lens_2h = _FALSE_;
  pclass_sz->has_gal_lens_hf = _FALSE_;
  pclass_sz->has_lens_lens_hf = _FALSE_;
  pclass_sz->has_gal_lensmag_1h = _FALSE_;
  pclass_sz->has_gal_lensmag_2h = _FALSE_;
  pclass_sz->has_gal_gallens_1h = _FALSE_;
  pclass_sz->has_gal_gallens_2h = _FALSE_;
  pclass_sz->has_tSZ_gallens_1h = _FALSE_;
  pclass_sz->has_tSZ_gallens_2h = _FALSE_;
  pclass_sz->has_gal_lensmag_hf = _FALSE_;
  pclass_sz->has_tSZ_lensmag_1h = _FALSE_;
  pclass_sz->has_tSZ_lensmag_2h = _FALSE_;
  pclass_sz->has_gallens_lensmag_1h = _FALSE_;
  pclass_sz->has_gallens_lensmag_2h = _FALSE_;
  pclass_sz->has_lensmag_lensmag_1h = _FALSE_;
  pclass_sz->has_lensmag_lensmag_2h = _FALSE_;
  pclass_sz->has_lensmag_lensmag_hf = _FALSE_;
  pclass_sz->has_lens_lensmag_1h = _FALSE_;
  pclass_sz->has_lens_lensmag_2h = _FALSE_;
  pclass_sz->has_lens_lensmag_hf = _FALSE_;
  pclass_sz->has_tSZ_gal_1h = _FALSE_;
  pclass_sz->has_IA_gal_2h = _FALSE_;
  pclass_sz->has_tSZ_gal_2h = _FALSE_;
  pclass_sz->has_tSZ_cib_1h = _FALSE_;
  pclass_sz->has_tSZ_cib_2h = _FALSE_;
  pclass_sz->has_gallens_cib_1h = _FALSE_;
  pclass_sz->has_gallens_cib_2h = _FALSE_;
  pclass_sz->has_gal_cib_1h = _FALSE_;
  pclass_sz->has_gal_cib_2h = _FALSE_;
  pclass_sz->has_lens_cib_1h = _FALSE_;
  pclass_sz->has_lens_cib_2h = _FALSE_;
  pclass_sz->has_ngal_ngal_1h = _FALSE_;
  pclass_sz->has_ngal_ngal_2h = _FALSE_;
  pclass_sz->has_ngal_ngal_hf = _FALSE_;
  pclass_sz->has_ngal_lens_1h = _FALSE_;
  pclass_sz->has_ngal_lens_2h = _FALSE_;
  pclass_sz->has_ngal_lens_hf = _FALSE_;
  pclass_sz->has_ngal_nlensmag_hf = _FALSE_;
  pclass_sz->has_ngal_tsz_1h = _FALSE_;
  pclass_sz->has_ngal_tsz_2h = _FALSE_;
  pclass_sz->has_nlensmag_tsz_1h = _FALSE_;
  pclass_sz->has_nlensmag_tsz_2h = _FALSE_;
  pclass_sz->has_ngal_gallens_1h = _FALSE_;
  pclass_sz->has_ngal_gallens_2h = _FALSE_;
  pclass_sz->has_ngal_IA_2h = _FALSE_;
  pclass_sz->has_nlensmag_gallens_1h = _FALSE_;
  pclass_sz->has_nlensmag_gallens_2h = _FALSE_;
  pclass_sz->has_cib_cib_1h = _FALSE_;
  pclass_sz->has_cib_cib_2h = _FALSE_;
  pclass_sz->has_pk_at_z_1h = _FALSE_;
  pclass_sz->has_pk_at_z_2h = _FALSE_;
  pclass_sz->has_pk_gg_at_z_1h = _FALSE_;
  pclass_sz->has_pk_gg_at_z_2h = _FALSE_;
  pclass_sz->has_pk_bb_at_z_1h = _FALSE_;
  pclass_sz->has_pk_bb_at_z_2h = _FALSE_;
  pclass_sz->has_pk_b_at_z_2h = _FALSE_;
  pclass_sz->has_gas_pressure_profile_2h = _FALSE_;
  pclass_sz->has_gas_density_profile_2h = _FALSE_;
  pclass_sz->has_pk_em_at_z_1h = _FALSE_;
  pclass_sz->has_pk_em_at_z_2h = _FALSE_;
  pclass_sz->has_pk_HI_at_z_1h = _FALSE_;
  pclass_sz->has_pk_HI_at_z_2h = _FALSE_;
  pclass_sz->has_bk_at_z_1h = _FALSE_;
  pclass_sz->has_bk_at_z_2h = _FALSE_;
  pclass_sz->has_bk_at_z_3h = _FALSE_;
  pclass_sz->has_bk_ttg_at_z_1h = _FALSE_;
  pclass_sz->has_bk_ttg_at_z_2h = _FALSE_;
  pclass_sz->has_bk_ttg_at_z_3h = _FALSE_;
  pclass_sz->has_bk_at_z_hf = _FALSE_;
  pclass_sz->has_lens_lens_1h = _FALSE_;
  pclass_sz->has_lens_lens_2h = _FALSE_;
  pclass_sz->has_custom1 = _FALSE_;
  pclass_sz->has_b_custom1 = 0;
  pclass_sz->has_custom1_custom1_1h = _FALSE_;
  pclass_sz->has_custom1_custom1_2h = _FALSE_;
  pclass_sz->has_custom1_lens_1h = _FALSE_;
  pclass_sz->has_custom1_lens_2h = _FALSE_;
  pclass_sz->has_custom1_tSZ_1h = _FALSE_;
  pclass_sz->has_custom1_tSZ_2h = _FALSE_;
  pclass_sz->has_custom1_cib_1h = _FALSE_;
  pclass_sz->has_custom1_cib_2h = _FALSE_;
  pclass_sz->has_custom1_gal_1h = _FALSE_;
  pclass_sz->has_custom1_gal_2h = _FALSE_;
  pclass_sz->has_custom1_gallens_1h = _FALSE_;
  pclass_sz->has_custom1_gallens_2h = _FALSE_;
  pclass_sz->has_tSZ_lens_1h = _FALSE_;
  pclass_sz->has_tSZ_lens_2h = _FALSE_;
  pclass_sz->has_kSZ_kSZ_gal_1h = _FALSE_;
  pclass_sz->has_kSZ_kSZ_gal_1h_fft = _FALSE_;
  pclass_sz->has_kSZ_kSZ_gal_2h_fft = _FALSE_;
  pclass_sz->has_kSZ_kSZ_gal_3h_fft = _FALSE_;
  pclass_sz->has_kSZ_kSZ_gal_covmat = _FALSE_;
  pclass_sz->has_kSZ_kSZ_gal_lensing_term = _FALSE_;
  pclass_sz->has_kSZ_kSZ_gal_2h = _FALSE_;
  pclass_sz->has_kSZ_kSZ_gal_3h = _FALSE_;
  pclass_sz->has_kSZ_kSZ_gal_hf = _FALSE_;
  pclass_sz->has_kSZ_kSZ_gallens_1h_fft = _FALSE_;
  pclass_sz->has_kSZ_kSZ_gallens_2h_fft = _FALSE_;
  pclass_sz->has_kSZ_kSZ_gallens_3h_fft = _FALSE_;
  pclass_sz->has_kSZ_kSZ_gallens_covmat = _FALSE_;
  pclass_sz->has_kSZ_kSZ_gallens_lensing_term = _FALSE_;
  pclass_sz->has_kSZ_kSZ_gallens_hf = _FALSE_;
  pclass_sz->has_kSZ_kSZ_lens_1h_fft = _FALSE_;
  pclass_sz->has_kSZ_kSZ_lens_2h_fft = _FALSE_;
  pclass_sz->has_kSZ_kSZ_lens_3h_fft = _FALSE_;
  pclass_sz->has_gal_gal_lens_1h_fft = _FALSE_;
  pclass_sz->has_gal_gal_lens_2h_fft = _FALSE_;
  pclass_sz->has_gal_gal_lens_3h_fft = _FALSE_;
  pclass_sz->has_kSZ_kSZ_lens_covmat = _FALSE_;
  pclass_sz->has_kSZ_kSZ_lens_lensing_term = _FALSE_;
  pclass_sz->has_kSZ_kSZ_lens_hf = _FALSE_;
  pclass_sz->has_gallens_gallens_1h = _FALSE_;
  pclass_sz->has_gallens_gallens_2h = _FALSE_;
  pclass_sz->has_gallens_lens_1h = _FALSE_;
  pclass_sz->has_gallens_lens_2h = _FALSE_;
  pclass_sz->has_mean_galaxy_bias = _FALSE_;
  pclass_sz->has_kSZ_kSZ_lensmag_1halo = _FALSE_;
  pclass_sz->has_tSZ_tSZ_tSZ_1halo = _FALSE_;
  pclass_sz->has_tSZ_tSZ_tSZ_2h = _FALSE_;
  pclass_sz->has_tSZ_tSZ_tSZ_3h = _FALSE_;
  pclass_sz->has_kSZ_kSZ_1h = _FALSE_;
  pclass_sz->has_kSZ_kSZ_2h = _FALSE_;
  pclass_sz->has_kSZ_kSZ_tSZ_1h = _FALSE_;
  pclass_sz->has_kSZ_kSZ_tSZ_2h = _FALSE_;
  pclass_sz->has_kSZ_kSZ_tSZ_3h = _FALSE_;
  pclass_sz->has_sz_te_y_y = _FALSE_;
  pclass_sz->has_sz_m_y_y_1h = _FALSE_;
  pclass_sz->has_sz_m_y_y_2h = _FALSE_;
  pclass_sz->has_sz_ps = _FALSE_;
  pclass_sz->has_sz_2halo = _FALSE_;
  pclass_sz->has_sz_trispec = _FALSE_;
  pclass_sz->has_hmf = _FALSE_;
  pclass_sz->has_cib_monopole = _FALSE_;
  pclass_sz->has_cib_shotnoise = _FALSE_;
  pclass_sz->has_dcib0dz = _FALSE_;
  pclass_sz->has_dydz = _FALSE_;
  pclass_sz->has_mean_y = _FALSE_;
  pclass_sz->has_sz_cov_Y_N = _FALSE_;
  pclass_sz->has_sz_cov_Y_Y_ssc = _FALSE_;
  pclass_sz->has_sz_cov_Y_N_next_order = _FALSE_;
  pclass_sz->has_sz_cov_N_N = _FALSE_;
  pclass_sz->has_sz_cov_N_N_hsv = _FALSE_;

  pclass_sz->has_vrms2 = _FALSE_;
  pclass_sz->has_knl = _FALSE_;
  pclass_sz->has_nl_index = _FALSE_;
  pclass_sz->has_sigma2_hsv = _FALSE_;
  pclass_sz->has_ng_in_bh = _FALSE_;
  pclass_sz->need_ng_bias = 0;

  pclass_sz->index_md_hmf = 0;
  pclass_sz->index_md_mean_y = 1;
  pclass_sz->index_md_sz_ps = 2;
  pclass_sz->index_md_trispectrum = 3;
  pclass_sz->index_md_2halo = 4;
  pclass_sz->index_md_te_y_y = 5;
  pclass_sz->index_md_cov_Y_N = 6;
  pclass_sz->index_md_cov_N_N = 7;
  pclass_sz->index_md_cov_Y_N_next_order = 8;
  pclass_sz->index_md_cov_N_N_hsv = 9;
  pclass_sz->index_md_kSZ_kSZ_gal_1h = 10;
  pclass_sz->index_md_tSZ_lens_1h = 11;
  pclass_sz->index_md_tSZ_lens_2h = 12;
  pclass_sz->index_md_isw_lens = 13;
  pclass_sz->index_md_isw_tsz = 14;
  pclass_sz->index_md_isw_auto = 15;
  pclass_sz->index_md_dndlnM = 16;
  pclass_sz->index_md_cov_Y_Y_ssc = 17;
  pclass_sz->index_md_m_y_y_1h = 18;
  pclass_sz->index_md_m_y_y_2h = 19;
  pclass_sz->index_md_tSZ_gal_1h = 20;
  pclass_sz->index_md_tSZ_tSZ_tSZ_1halo = 21;
  pclass_sz->index_md_gal_gal_1h = 22;
  pclass_sz->index_md_gal_gal_2h = 23;
  pclass_sz->index_md_gal_lens_2h = 24;
  pclass_sz->index_md_gal_lens_1h = 25;
  pclass_sz->index_md_tSZ_gal_2h = 26;
  pclass_sz->index_md_lens_lens_2h = 27;
  pclass_sz->index_md_lens_lens_1h = 28;
  pclass_sz->index_md_tSZ_cib_1h = 29;
  pclass_sz->index_md_tSZ_cib_2h = 30;
  pclass_sz->index_md_cib_cib_1h = 31;
  pclass_sz->index_md_cib_cib_2h = 32;
  pclass_sz->index_md_lens_cib_1h = 33;
  pclass_sz->index_md_lens_cib_2h = 34;
  pclass_sz->index_md_pk_at_z_1h = 35;
  pclass_sz->index_md_pk_at_z_2h = 36;
  pclass_sz->index_md_gal_lensmag_2h = 37;
  pclass_sz->index_md_gal_lensmag_1h = 38;
  pclass_sz->index_md_lensmag_lensmag_2h = 39;
  pclass_sz->index_md_lensmag_lensmag_1h = 40;
  pclass_sz->index_md_lens_lensmag_2h = 41;
  pclass_sz->index_md_lens_lensmag_1h = 42;
  pclass_sz->index_md_kSZ_kSZ_gal_2h = 43;
  pclass_sz->index_md_kSZ_kSZ_gal_3h = 44;
  pclass_sz->index_md_tSZ_lensmag_2h = 45;
  pclass_sz->index_md_tSZ_lensmag_1h = 46;
  pclass_sz->index_md_bk_at_z_1h = 47;
  pclass_sz->index_md_bk_at_z_2h = 48;
  pclass_sz->index_md_bk_at_z_3h = 49;

  pclass_sz->index_md_kSZ_kSZ_gal_hf = 50;
  pclass_sz->index_md_gal_gal_hf = 51;
  pclass_sz->index_md_gal_lens_hf = 52;
  pclass_sz->index_md_gal_lensmag_hf = 53;
  pclass_sz->index_md_lens_lensmag_hf = 54;
  pclass_sz->index_md_lensmag_lensmag_hf = 55;
  pclass_sz->index_md_pk_gg_at_z_1h = 56;
  pclass_sz->index_md_pk_gg_at_z_2h = 57;
  pclass_sz->index_md_gal_cib_1h = 58;
  pclass_sz->index_md_gal_cib_2h = 59;
  pclass_sz->index_md_kSZ_kSZ_gal_1h_fft = 60;
  pclass_sz->index_md_kSZ_kSZ_gal_2h_fft = 61;
  pclass_sz->index_md_kSZ_kSZ_gal_3h_fft = 62;
  pclass_sz->index_md_cib_monopole = 63;
  pclass_sz->index_md_dcib0dz = 64;
  pclass_sz->index_md_dydz = 65;

  pclass_sz->index_md_bk_ttg_at_z_1h = 66;
  pclass_sz->index_md_bk_ttg_at_z_2h = 67;
  pclass_sz->index_md_bk_ttg_at_z_3h = 68;

  pclass_sz->index_md_gal_gallens_2h = 69;
  pclass_sz->index_md_gal_gallens_1h = 70;

  pclass_sz->index_md_kSZ_kSZ_tSZ_1h = 71;
  pclass_sz->index_md_kSZ_kSZ_tSZ_2h = 72;
  pclass_sz->index_md_kSZ_kSZ_tSZ_3h = 73;
  // pclass_sz->index_md_bk_at_z_hf = 51;
  pclass_sz->index_md_kSZ_kSZ_1h = 74;
  pclass_sz->index_md_kSZ_kSZ_2h = 75;

  pclass_sz->index_md_pk_bb_at_z_1h = 76;
  pclass_sz->index_md_pk_bb_at_z_2h = 77;

  pclass_sz->index_md_kSZ_kSZ_gallens_1h_fft = 78;
  pclass_sz->index_md_kSZ_kSZ_gallens_2h_fft = 79;
  pclass_sz->index_md_kSZ_kSZ_gallens_3h_fft = 80;
  pclass_sz->index_md_kSZ_kSZ_gallens_hf = 81;
  pclass_sz->index_md_gallens_gallens_1h = 82;
  pclass_sz->index_md_gallens_gallens_2h = 83;
  pclass_sz->index_md_gallens_lens_1h = 84;
  pclass_sz->index_md_gallens_lens_2h = 85;

  pclass_sz->index_md_kSZ_kSZ_lens_1h_fft = 86;
  pclass_sz->index_md_kSZ_kSZ_lens_2h_fft = 87;
  pclass_sz->index_md_kSZ_kSZ_lens_3h_fft = 88;
  pclass_sz->index_md_kSZ_kSZ_lens_hf = 89;

  pclass_sz->index_md_lens_lens_hf = 90;

  pclass_sz->index_md_pk_em_at_z_1h = 91;
  pclass_sz->index_md_pk_em_at_z_2h = 92;
  pclass_sz->index_md_szrates = 93;
  pclass_sz->index_md_pk_HI_at_z_2h = 94;

  pclass_sz->index_md_tSZ_tSZ_tSZ_2h = 95;
  pclass_sz->index_md_tSZ_tSZ_tSZ_3h = 96;

  pclass_sz->index_md_cib_shotnoise= 97;

  pclass_sz->index_md_ngal_ngal_1h = 98;
  pclass_sz->index_md_ngal_ngal_2h = 99;
  pclass_sz->index_md_ngal_ngal_hf = 100;

  pclass_sz->index_md_ngal_lens_1h = 101;
  pclass_sz->index_md_ngal_lens_2h = 102;
  pclass_sz->index_md_ngal_lens_hf = 103;

  pclass_sz->index_md_pk_b_at_z_2h = 104;

  pclass_sz->index_md_tSZ_gallens_2h = 105;
  pclass_sz->index_md_tSZ_gallens_1h = 106;

  pclass_sz->index_md_gallens_cib_1h = 107;
  pclass_sz->index_md_gallens_cib_2h = 108;

  pclass_sz->index_md_tau_gal_2h = 109;
  pclass_sz->index_md_tau_gal_1h = 110;

  pclass_sz->index_md_IA_gal_2h = 111;

  pclass_sz->index_md_gal_gal_lens_1h_fft = 112;
  pclass_sz->index_md_gal_gal_lens_2h_fft = 113;
  pclass_sz->index_md_gal_gal_lens_3h_fft = 114;

  pclass_sz->index_md_custom1_custom1_1h = 115;
  pclass_sz->index_md_custom1_custom1_2h = 116;
  pclass_sz->index_md_custom1_lens_1h = 117;
  pclass_sz->index_md_custom1_lens_2h = 118;
  pclass_sz->index_md_custom1_tSZ_1h = 119;
  pclass_sz->index_md_custom1_tSZ_2h = 120;
  pclass_sz->index_md_custom1_cib_1h = 121;
  pclass_sz->index_md_custom1_cib_2h = 122;
  pclass_sz->index_md_custom1_gal_1h = 123;
  pclass_sz->index_md_custom1_gal_2h = 124;
  pclass_sz->index_md_custom1_gallens_1h = 125;
  pclass_sz->index_md_custom1_gallens_2h = 126;


  pclass_sz->index_md_gallens_lensmag_1h = 127;
  pclass_sz->index_md_gallens_lensmag_2h = 128;

  pclass_sz->index_md_ngal_nlensmag_hf = 129;

  pclass_sz->index_md_tau_tau_2h = 130;
  pclass_sz->index_md_tau_tau_1h = 131;

  pclass_sz->index_md_ngal_tsz_1h = 132;
  pclass_sz->index_md_ngal_tsz_2h = 133;

  pclass_sz->index_md_ngal_gallens_1h = 134;
  pclass_sz->index_md_ngal_gallens_2h = 135;
  pclass_sz->index_md_nlensmag_gallens_1h = 136;
  pclass_sz->index_md_nlensmag_gallens_2h = 137;
  pclass_sz->index_md_ngal_IA_2h = 138;
  pclass_sz->index_md_nlensmag_tsz_1h = 139;
  pclass_sz->index_md_nlensmag_tsz_2h = 140;

  pclass_sz->integrate_wrt_mvir = 0;
  pclass_sz->integrate_wrt_m500c = 0;
  pclass_sz->integrate_wrt_m200c = 1;

  pclass_sz->delta_def_cib = 1;
  pclass_sz->delta_def_galaxies = 1;
  pclass_sz->delta_def_matter_density = 1;
  pclass_sz->delta_def_electron_density = 1;
  pclass_sz->delta_def_electron_pressure = 1;
  pclass_sz->delta_def_custom1 = 1;


  pclass_sz->profile_matter_density = 0;
  pclass_sz->matter_nfw_power_law_index = 0.;


  pclass_sz->integrate_wrt_m200m = 0;

  pclass_sz->need_m200c_to_m200m = 0;
  pclass_sz->need_m200m_to_mvir = 0;
  pclass_sz->need_m200m_to_m200c = 0;
  pclass_sz->need_m200m_to_m500c = 0;
  pclass_sz->need_m200c_to_m500c = 0;
  pclass_sz->need_m500c_to_m200c = 0;
  pclass_sz->need_hmf = 0;
  pclass_sz->need_sigma = 0;
  pclass_sz->need_ksz_template = 0;
  pclass_sz->need_tt_noise = 0;
  pclass_sz->need_lensing_noise=0;
  pclass_sz->has_electron_pressure = 0;
  pclass_sz->has_electron_density = 0;
  pclass_sz->has_HI_density = 0;
  pclass_sz->has_galaxy = 0;
  pclass_sz->has_matter_density = 0;
  pclass_sz->has_lensing = 0;
  pclass_sz->has_cib = 0;
  pclass_sz->has_isw = 0;

  pclass_sz->has_vir = 0;
  pclass_sz->has_500c = 0;
  pclass_sz->has_200m = 0;
  pclass_sz->has_200c = 0;


  pclass_sz->HMF_prescription_NCDM=1; //CDM


  pcsz->size_logM = 100;

  pcsz->z_max = 1.;

  pclass_sz->sz_verbose = 0;

  pclass_sz->f_free  = 1.; //  Ionization state of Helium (0.86 = neutral, 0.93 = singly ionized, 1 = completely ionized for Y_p = 0.24)
  pclass_sz->f_b_gas  = -1.;
  pclass_sz->mu_e = 1.14;

  pclass_sz->csat_over_cdm = 1.;
  pclass_sz->cl_gal_gal_A_sn = 0.;
  pclass_sz->shape_noise_siggamma2 = 0.3;

  //HOD
  pclass_sz->M_min_HOD = pow(10,11.5); //Msun/h
  pclass_sz->M_min_HOD_cib = pow(10,11.5); //Msun/h
  pclass_sz->M1_prime_HOD = pow(10,12.6); //Msun/h
  pclass_sz->alpha_s_HOD = 1.;
  pclass_sz->sigma_log10M_HOD = 0.15;
  pclass_sz->rho_y_gal = -0.6;


  pclass_sz->include_noise_cov_y_y = 0;


  return _SUCCESS_;

}

/**
 * Initialize the precision parameter structure.
 *
 * All precision parameters used in the other modules are listed here
 * and assigned here a default value.
 *
 * @param ppr Input/Output: a precision_params structure pointer
 * @return the error status
 *
 */

int input_default_precision ( struct precision * ppr ) {

  /**
   * - automatic estimate of machine precision
   */
  ppr->smallest_allowed_variation=DBL_EPSILON;

  //get_machine_precision(&(ppr->smallest_allowed_variation));

  class_test(ppr->smallest_allowed_variation < 0,
             ppr->error_message,
             "smallest_allowed_variation = %e < 0",ppr->smallest_allowed_variation);


#define __ASSIGN_DEFAULT_PRECISION__
#include "precisions.h"
#undef __ASSIGN_DEFAULT_PRECISION__


return _SUCCESS_;

}

int class_version(
                  char * version
                  ) {

  sprintf(version,"%s",_VERSION_);
  return _SUCCESS_;
}

/**
 * Automatically computes the machine precision.
 *
 * @param smallest_allowed_variation a pointer to the smallest allowed variation
 *
 * Returns the smallest
 * allowed variation (minimum epsilon * _TOLVAR_)
 */

int get_machine_precision(double * smallest_allowed_variation) {
  double one, meps, sum;

  one = 1.0;
  meps = 1.0;
  do {
    meps /= 2.0;
    sum = one + meps;
  } while (sum != one);
  meps *= 2.0;

  *smallest_allowed_variation = meps * _TOLVAR_;

  return _SUCCESS_;

}

int input_fzerofun_1d(double input,
                      void* pfzw,
                      double *output,
                      ErrorMsg error_message){

  class_call(input_try_unknown_parameters(&input,
                                          1,
                                          pfzw,
                                          output,
                                          error_message),
             error_message,
             error_message);

  return _SUCCESS_;
}

int class_fzero_ridder(int (*func)(double x, void *param, double *y, ErrorMsg error_message),
                       double x1,
                       double x2,
                       double xtol,
                       void *param,
                       double *Fx1,
                       double *Fx2,
                       double *xzero,
                       int *fevals,
                       ErrorMsg error_message){
  /**Using Ridders' method, return the root of a function func known to
     lie between x1 and x2. The root, returned as zriddr, will be found to
     an approximate accuracy xtol.
  */
  int j,MAXIT=1000;
  double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;
  if ((Fx1!=NULL)&&(Fx2!=NULL)){
    fl = *Fx1;
    fh = *Fx2;
  }
  else{
    class_call((*func)(x1, param, &fl, error_message),
               error_message, error_message);
    class_call((*func)(x2, param, &fh, error_message),
               error_message, error_message);

    *fevals = (*fevals)+2;
  }
  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    xl=x1;
    xh=x2;
    ans=-1.11e11;
    for (j=1;j<=MAXIT;j++) {
      xm=0.5*(xl+xh);
      class_call((*func)(xm, param, &fm, error_message),
                 error_message, error_message);
      *fevals = (*fevals)+1;
      s=sqrt(fm*fm-fl*fh);
      if (s == 0.0){
        *xzero = ans;
        //printf("Success 1\n");
        return _SUCCESS_;
      }
      xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
      if (fabs(xnew-ans) <= xtol) {
        *xzero = ans;
        return _SUCCESS_;
      }
      ans=xnew;
      class_call((*func)(ans, param, &fnew, error_message),
                 error_message, error_message);
      *fevals = (*fevals)+1;
      if (fnew == 0.0){
        *xzero = ans;
        //printf("Success 2, ans=%g\n",ans);
        return _SUCCESS_;
      }
      if (NRSIGN(fm,fnew) != fm) {
        xl=xm;
        fl=fm;
        xh=ans;
        fh=fnew;
      } else if (NRSIGN(fl,fnew) != fl) {
        xh=ans;
        fh=fnew;
      } else if (NRSIGN(fh,fnew) != fh) {
        xl=ans;
        fl=fnew;
      } else return _FAILURE_;
      if (fabs(xh-xl) <= xtol) {
        *xzero = ans;
        //        printf("Success 3\n");
        return _SUCCESS_;
      }
    }
    class_stop(error_message,"zriddr exceed maximum iterations");
  }
  else {
    if (fl == 0.0) return x1;
    if (fh == 0.0) return x2;
    class_stop(error_message,"root must be bracketed in zriddr.");
  }
  class_stop(error_message,"Failure in int.");
}

int input_try_unknown_parameters(double * unknown_parameter,
                                 int unknown_parameters_size,
                                 void * voidpfzw,
                                 double * output,
                                 ErrorMsg errmsg){
  /** Summary:
   * - Call the structures*/

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct nonlinear nl;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct class_sz_structure tsz;     /* for tsz spectrum */ //BB: added for class_sz
  struct szcount csz;         /* for sz cluster count */ //BB: added for class_sz
  struct output op;           /* for output files */

  int i;
  double rho_dcdm_today, rho_dr_today;
  struct fzerofun_workspace * pfzw;
  int input_verbose;
  int flag;
  int param;
  short compute_sigma8 = _FALSE_;

  pfzw = (struct fzerofun_workspace *) voidpfzw;
  /** - Read input parameters */
  for (i=0; i < unknown_parameters_size; i++) {
    sprintf(pfzw->fc.value[pfzw->unknown_parameters_index[i]],
            "%e",unknown_parameter[i]);
  }

  class_call(input_read_precisions(&(pfzw->fc),
                                   &pr,
                                   &ba,
                                   &th,
                                   &pt,
                                   &tr,
                                   &pm,
                                   &sp,
                                   &nl,
                                   &le,
                                   &tsz,
                                   &op,
                                   errmsg),
             errmsg,
             errmsg);

  class_call(input_read_parameters(&(pfzw->fc),
                                   &pr,
                                   &ba,
                                   &th,
                                   &pt,
                                   &tr,
                                   &pm,
                                   &sp,
                                   &nl,
                                   &le,
                                   &tsz,
                                   &csz,
                                   &op,
                                   errmsg),
             errmsg,
             errmsg);

  class_call(parser_read_int(&(pfzw->fc),
                             "input_verbose",
                             &param,
                             &flag,
                             errmsg),
             errmsg,
             errmsg);

  if (flag == _TRUE_)
    input_verbose = param;
  else
    input_verbose = 0;

  /** - Optimise flags for sigma8 calculation.*/
  for (i=0; i < unknown_parameters_size; i++) {
    if (pfzw->target_name[i] == sigma8) {
      compute_sigma8 = _TRUE_;
    }
  }
  if (compute_sigma8 == _TRUE_) {
    pt.k_max_for_pk=10.0; // increased in June 2020 for higher accuracy
    pt.has_pk_matter=_TRUE_;
    pt.has_perturbations = _TRUE_;
    pt.has_cl_cmb_temperature = _FALSE_;
    pt.has_cls = _FALSE_;
    pt.has_cl_cmb_polarization = _FALSE_;
    pt.has_cl_cmb_lensing_potential = _FALSE_;
    pt.has_cl_number_count = _FALSE_;
    pt.has_cl_lensing_potential=_FALSE_;
    pt.has_density_transfers=_FALSE_;
    pt.has_velocity_transfers=_FALSE_;
    nl.has_pk_eq=_FALSE_; // not needed since sigma8 is derived from linear P(k)
    nl.method=nl_none;    // not needed since sigma8 is derived from linear P(k)
  }

  /** - Shoot forward into class up to required stage */
  if (pfzw->required_computation_stage >= cs_background){
    if (input_verbose>2)
      printf("Stage 1: background\n");
    ba.background_verbose = 0;
    class_call(background_init(&pr,&ba),
               ba.error_message,
               errmsg
               );
  }

  if (pfzw->required_computation_stage >= cs_thermodynamics){
    if (input_verbose>2)
      printf("Stage 2: thermodynamics\n");
    pr.recfast_Nz0 = 10000;
    th.thermodynamics_verbose = 0;
    class_call_except(thermodynamics_init(&pr,&ba,&th),
                      th.error_message,
                      errmsg,
                      background_free(&ba)
                      );
  }

  if (pfzw->required_computation_stage >= cs_perturbations){
    if (input_verbose>2)
      printf("Stage 3: perturbations\n");
    pt.perturbations_verbose = 0;
    class_call_except(perturb_init(&pr,&ba,&th,&pt),
                      pt.error_message,
                      errmsg,
                      thermodynamics_free(&th);background_free(&ba)
                      );
  }

  if (pfzw->required_computation_stage >= cs_primordial){
    if (input_verbose>2)
      printf("Stage 4: primordial\n");
    pm.primordial_verbose = 0;
    class_call_except(primordial_init(&pr,&pt,&pm),
                      pm.error_message,
                      errmsg,
                      perturb_free(&pt);thermodynamics_free(&th);background_free(&ba)
                      );
  }

  if (pfzw->required_computation_stage >= cs_nonlinear){
    if (input_verbose>2)
      printf("Stage 5: nonlinear\n");
    nl.nonlinear_verbose = 0;
    class_call_except(nonlinear_init(&pr,&ba,&th,&pt,&pm,&nl),
                      nl.error_message,
                      errmsg,
                      primordial_free(&pm);perturb_free(&pt);thermodynamics_free(&th);background_free(&ba)
                      );
  }

  if (pfzw->required_computation_stage >= cs_transfer){
    if (input_verbose>2)
      printf("Stage 6: transfer\n");
    tr.transfer_verbose = 0;
    class_call_except(transfer_init(&pr,&ba,&th,&pt,&nl,&tr),
                      tr.error_message,
                      errmsg,
                      nonlinear_free(&nl);primordial_free(&pm);perturb_free(&pt);thermodynamics_free(&th);background_free(&ba)
                      );
  }

  if (pfzw->required_computation_stage >= cs_spectra){
    if (input_verbose>2)
      printf("Stage 7: spectra\n");
    sp.spectra_verbose = 0;
    class_call_except(spectra_init(&pr,&ba,&pt,&pm,&nl,&tr,&sp,&tsz,&th,&le),
                      sp.error_message,
                      errmsg,
                      transfer_free(&tr);nonlinear_free(&nl);primordial_free(&pm);perturb_free(&pt);thermodynamics_free(&th);background_free(&ba)
                      );
  }

  // printf(">>> shooting:\n");
  // exit(0);
  /** - Get the corresponding shoot variable and put into output */
  for (i=0; i < pfzw->target_size; i++) {
    switch (pfzw->target_name[i]) {
    case theta_s:
      output[i] = 100.*th.rs_rec/th.ra_rec-pfzw->target_value[i];
      break;
    case Omega_dcdmdr:
      rho_dcdm_today = ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_rho_dcdm];
      if (ba.has_dr == _TRUE_)
        rho_dr_today = ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_rho_dr];
      else
        rho_dr_today = 0.;
      output[i] = (rho_dcdm_today+rho_dr_today)/(ba.H0*ba.H0)-pfzw->target_value[i];
      break;
    case omega_dcdmdr:
      rho_dcdm_today = ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_rho_dcdm];
      if (ba.has_dr == _TRUE_)
        rho_dr_today = ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_rho_dr];
      else
        rho_dr_today = 0.;
      output[i] = (rho_dcdm_today+rho_dr_today)/(ba.H0*ba.H0)-pfzw->target_value[i]/ba.h/ba.h;
      break;
    case Omega_scf:
      /** - In case scalar field is used to fill, pba->Omega0_scf is not equal to pfzw->target_value[i].*/
      output[i] = ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_rho_scf]/(ba.H0*ba.H0)
        -ba.Omega0_scf;
      break;
    case Omega_ini_dcdm:
    case omega_ini_dcdm:
      rho_dcdm_today = ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_rho_dcdm];
      if (ba.has_dr == _TRUE_)
        rho_dr_today = ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_rho_dr];
      else
        rho_dr_today = 0.;
      output[i] = -(rho_dcdm_today+rho_dr_today)/(ba.H0*ba.H0)+ba.Omega0_dcdmdr;
      break;
    case sigma8:
      output[i] = nl.sigma8[nl.index_pk_m]-pfzw->target_value[i];
      break;
    case age:
      output[i] = ba.age-pfzw->target_value[i];
      break;

    }
  }


  /** - Free structures */
  if (pfzw->required_computation_stage >= cs_spectra){
    class_call(spectra_free(&sp), sp.error_message, errmsg);
  }
  if (pfzw->required_computation_stage >= cs_transfer){
    class_call(transfer_free(&tr), tr.error_message, errmsg);
  }
  if (pfzw->required_computation_stage >= cs_nonlinear){
    class_call(nonlinear_free(&nl), nl.error_message, errmsg);
  }
  if (pfzw->required_computation_stage >= cs_primordial){
    class_call(primordial_free(&pm), pm.error_message, errmsg);
  }
  if (pfzw->required_computation_stage >= cs_perturbations){
    class_call(perturb_free(&pt), pt.error_message, errmsg);
  }
  if (pfzw->required_computation_stage >= cs_thermodynamics){
    class_call(thermodynamics_free(&th), th.error_message, errmsg);
  }
  if (pfzw->required_computation_stage >= cs_background){
    class_call(background_free(&ba), ba.error_message, errmsg);
  }

  /** - Set filecontent to unread */
  for (i=0; i<pfzw->fc.size; i++) {
    pfzw->fc.read[i] = _FALSE_;
  }

  return _SUCCESS_;
}

int input_get_guess(double *xguess,
                    double *dxdy,
                    struct fzerofun_workspace * pfzw,
                    ErrorMsg errmsg){

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct nonlinear nl;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct class_sz_structure tsz;     //BB: added for class_sz
  struct szcount csz;         //BB: added for class_sz
  struct output op;           /* for output files */
  int i;

  double Omega_M, a_decay, gamma, Omega0_dcdmdr=1.0;
  int index_guess;

  /* Cheat to read only known parameters: */
  pfzw->fc.size -= pfzw->target_size;

  class_call(input_read_precisions(&(pfzw->fc),
                                   &pr,
                                   &ba,
                                   &th,
                                   &pt,
                                   &tr,
                                   &pm,
                                   &sp,
                                   &nl,
                                   &le,
                                   &tsz, //BB: added for class_sz
                                   //&csz, //BB: added for class_sz
                                   &op,
                                   errmsg),
             errmsg,
             errmsg);

  class_call(input_read_parameters(&(pfzw->fc),
                                   &pr,
                                   &ba,
                                   &th,
                                   &pt,
                                   &tr,
                                   &pm,
                                   &sp,
                                   &nl,
                                   &le,
                                   &tsz, //BB: added for class_sz
                                   &csz, //BB: added for class_sz
                                   &op,
                                   errmsg),
             errmsg,
             errmsg);

  pfzw->fc.size += pfzw->target_size;
  /** Summary: */
  /** - Here we should write reasonable guesses for the unknown parameters.
      Also estimate dxdy, i.e. how the unknown parameter responds to the known.
      This can simply be estimated as the derivative of the guess formula.*/

  for (index_guess=0; index_guess < pfzw->target_size; index_guess++) {
    switch (pfzw->target_name[index_guess]) {
    case theta_s:
      xguess[index_guess] = 3.54*pow(pfzw->target_value[index_guess],2)-5.455*pfzw->target_value[index_guess]+2.548;
      dxdy[index_guess] = (7.08*pfzw->target_value[index_guess]-5.455);
      /** - Update pb to reflect guess */
      ba.h = xguess[index_guess];
      ba.H0 = ba.h *  1.e5 / _c_;
      break;
    case Omega_dcdmdr:
      Omega_M = ba.Omega0_cdm+ba.Omega0_idm_dr+ba.Omega0_dcdmdr+ba.Omega0_b;
      /* This formula is exact in a Matter + Lambda Universe, but only
         for Omega_dcdm, not the combined.
         sqrt_one_minus_M = sqrt(1.0 - Omega_M);
         xguess[index_guess] = pfzw->target_value[index_guess]*
         exp(2./3.*ba.Gamma_dcdm/ba.H0*
         atanh(sqrt_one_minus_M)/sqrt_one_minus_M);
         dxdy[index_guess] = 1.0;//exp(2./3.*ba.Gamma_dcdm/ba.H0*atanh(sqrt_one_minus_M)/sqrt_one_minus_M);
      */
      gamma = ba.Gamma_dcdm/ba.H0;
      if (gamma < 1)
        a_decay = 1.0;
      else
        a_decay = pow(1+(gamma*gamma-1.)/Omega_M,-1./3.);
      xguess[index_guess] = pfzw->target_value[index_guess]/a_decay;
      dxdy[index_guess] = 1./a_decay;
      //printf("x = Omega_ini_guess = %g, dxdy = %g\n",*xguess,*dxdy);
      break;
    case omega_dcdmdr:
      Omega_M = ba.Omega0_cdm+ba.Omega0_idm_dr+ba.Omega0_dcdmdr+ba.Omega0_b;
      /* This formula is exact in a Matter + Lambda Universe, but only
         for Omega_dcdm, not the combined.
         sqrt_one_minus_M = sqrt(1.0 - Omega_M);
         xguess[index_guess] = pfzw->target_value[index_guess]*
         exp(2./3.*ba.Gamma_dcdm/ba.H0*
         atanh(sqrt_one_minus_M)/sqrt_one_minus_M);
         dxdy[index_guess] = 1.0;//exp(2./3.*ba.Gamma_dcdm/ba.H0*atanh(sqrt_one_minus_M)/sqrt_one_minus_M);
      */
      gamma = ba.Gamma_dcdm/ba.H0;
      if (gamma < 1)
        a_decay = 1.0;
      else
        a_decay = pow(1+(gamma*gamma-1.)/Omega_M,-1./3.);
      xguess[index_guess] = pfzw->target_value[index_guess]/ba.h/ba.h/a_decay;
      dxdy[index_guess] = 1./a_decay/ba.h/ba.h;
      //printf("x = Omega_ini_guess = %g, dxdy = %g\n",*xguess,*dxdy);
      break;
    case Omega_scf:

      /** - This guess is arbitrary, something nice using WKB should be implemented.
       *
       * - Version 2: use a fit: `xguess[index_guess] = 1.77835*pow(ba.Omega0_scf,-2./7.);
       * dxdy[index_guess] = -0.5081*pow(ba.Omega0_scf,-9./7.)`;
       *
       * - Version 3: use attractor solution */

      if (ba.scf_tuning_index == 0){
        xguess[index_guess] = sqrt(3.0/ba.Omega0_scf);
        dxdy[index_guess] = -0.5*sqrt(3.0)*pow(ba.Omega0_scf,-1.5);
      }
      else{
        /* Default: take the passed value as xguess and set dxdy to 1. */
        xguess[index_guess] = ba.scf_parameters[ba.scf_tuning_index];
        dxdy[index_guess] = 1.;
      }
      break;
    case omega_ini_dcdm:
      Omega0_dcdmdr = 1./(ba.h*ba.h);
    case Omega_ini_dcdm:
      /** - This works since correspondence is
          Omega_ini_dcdm -> Omega_dcdmdr and
          omega_ini_dcdm -> omega_dcdmdr */
      Omega0_dcdmdr *=pfzw->target_value[index_guess];
      Omega_M = ba.Omega0_cdm+ba.Omega0_idm_dr+Omega0_dcdmdr+ba.Omega0_b;
      gamma = ba.Gamma_dcdm/ba.H0;
      if (gamma < 1)
        a_decay = 1.0;
      else
        a_decay = pow(1+(gamma*gamma-1.)/Omega_M,-1./3.);
      xguess[index_guess] = pfzw->target_value[index_guess]*a_decay;
      dxdy[index_guess] = a_decay;
      if (gamma > 100)
        dxdy[index_guess] *= gamma/100;

      //printf("x = Omega_ini_guess = %g, dxdy = %g\n",*xguess,*dxdy);
      break;

    case sigma8:
      /* Assume linear relationship between A_s and sigma8 and fix coefficient
         according to vanilla LambdaCDM. Should be good enough... */
      xguess[index_guess] = 2.43e-9/0.87659*pfzw->target_value[index_guess];
      dxdy[index_guess] = 2.43e-9/0.87659;
      break;
    case age:
      xguess[index_guess] = 0.96/pfzw->target_value[index_guess]/_Gyr_over_Mpc_*_c_/1.e3;
      dxdy[index_guess] = -0.96/pfzw->target_value[index_guess]/pfzw->target_value[index_guess]/_Gyr_over_Mpc_*_c_/1.e3;
      /** - Update pb to reflect guess */
      // ba.H0 = xguess[index_guess]; //in Mpc^-1
      // ba.h = ba.H0*_c_/1.e5;
      // printf(">>> getting guess for H0, got: %.3e inv Mpc or %.3e\n",
      // xguess[index_guess],
      // xguess[index_guess]);
      break;
    }
    //printf("xguess = %g\n",xguess[index_guess]);
  }

  for (i=0; i<pfzw->fc.size; i++) {
    pfzw->fc.read[i] = _FALSE_;
  }

  /** - Deallocate everything allocated by input_read_parameters */
  background_free_input(&ba);

  return _SUCCESS_;
}

int input_find_root(double *xzero,
                    int *fevals,
                    struct fzerofun_workspace *pfzw,
                    double tol,
                    ErrorMsg errmsg){
  double x1, x2, f1, f2, dxdy, dx;
  int iter, iter2;
  int return_function;
  /** Summary: */

  /** - Fisrt we do our guess */
  class_call(input_get_guess(&x1, &dxdy, pfzw, errmsg),
             errmsg, errmsg);
  //      printf("x1= %g\n",x1);

  class_call(input_fzerofun_1d(x1,
                               pfzw,
                               &f1,
                               errmsg),
             errmsg, errmsg);

  (*fevals)++;
  //printf("x1= %g, f1= %g\n",x1,f1);

  dx = 1.5*f1*dxdy;

  /** - Then we do a linear hunt for the boundaries */
  for (iter=1; iter<=15; iter++){
    //x2 = x1 + search_dir*dx;
    x2 = x1 - dx;

    for (iter2=1; iter2 <= 3; iter2++) {
      return_function = input_fzerofun_1d(x2,pfzw,&f2,errmsg);
      (*fevals)++;
      //printf("x2= %g, f2= %g\n",x2,f2);
      //fprintf(stderr,"iter2=%d\n",iter2);

      if (return_function ==_SUCCESS_) {
        break;
      }
      else if (iter2 < 3) {
        dx*=0.5;
        x2 = x1-dx;
      }
      else {
        //fprintf(stderr,"get here\n");
        class_stop(errmsg,errmsg);
      }
    }

    if (f1*f2<0.0){
      /** - root has been bracketed */
      if (0==1){
        printf("Root has been bracketed after %d iterations: [%g, %g].\n",iter,x1,x2);
      }
      break;
    }

    x1 = x2;
    f1 = f2;
  }

  /** - Find root using Ridders method. (Exchange for bisection if you are old-school.)*/
  class_call(class_fzero_ridder(input_fzerofun_1d,
                                x1,
                                x2,
                                tol*MAX(fabs(x1),fabs(x2)),
                                pfzw,
                                &f1,
                                &f2,
                                xzero,
                                fevals,
                                errmsg),
             errmsg,errmsg);

  return _SUCCESS_;
}

int file_exists(const char *fname){
  FILE *file = fopen(fname, "r");
  if (file != NULL){
    fclose(file);
    return _TRUE_;
  }
  return _FALSE_;
}

int input_auxillary_target_conditions(struct file_content * pfc,
                                      enum target_names target_name,
                                      double target_value,
                                      int * aux_flag,
                                      ErrorMsg errmsg){
  *aux_flag = _TRUE_;
  /*
    double param1;
    int int1, flag1;
    int input_verbose = 0;
    class_read_int("input_verbose",input_verbose);
  */
  switch (target_name){
  case Omega_dcdmdr:
  case omega_dcdmdr:
  case Omega_scf:
  case Omega_ini_dcdm:
  case omega_ini_dcdm:
    /* Check that Omega's or omega's are nonzero: */
    if (target_value == 0.)
      *aux_flag = _FALSE_;
    break;
  default:
    /* Default is no additional checks */
    *aux_flag = _TRUE_;
    break;
  }
  return _SUCCESS_;
}

int compare_integers (const void * elem1, const void * elem2) {
  int f = *((int*)elem1);
  int s = *((int*)elem2);
  if (f > s) return  1;
  if (f < s) return -1;
  return 0;
}

int compare_doubles(const void *a,const void *b) {
  double *x = (double *) a;
  double *y = (double *) b;
  if (*x < *y)
    return -1;
  else if
    (*x > *y) return 1;
  return 0;
}


/**
 * Perform preliminary steps fur using the method called Pk_equal,
 * described in 0810.0190 and 1601.07230, extending the range of
 * validity of HALOFIT from constant w to (w0,wa) models. In that
 * case, one must compute here some effective values of w0_eff(z_i)
 * and Omega_m_eff(z_i), that will be interpolated later at arbitrary
 * redshift in the non-linear module.
 *
 * Returns table of values [z_i, tau_i, w0_eff_i, Omega_m_eff_i]
 * stored in nonlinear structure.
 *
 * @param ppr           Input: pointer to precision structure
 * @param pba           Input: pointer to background structure
 * @param pth           Input: pointer to thermodynamics structure
 * @param pnl    Input/Output: pointer to nonlinear structure
 * @param input_verbose Input: verbosity of this input module
 * @param errmsg  Input/Ouput: error message
 */

int input_prepare_pk_eq(
                        struct precision * ppr,
                        struct background *pba,
                        struct thermo *pth,
                        struct nonlinear *pnl,
                        int input_verbose,
                        ErrorMsg errmsg
                        ) {

  /** Summary: */

  /** - define local variables */

  double tau_of_z;
  double delta_tau;
  double error;
  double delta_tau_eq;
  double * pvecback;
  int last_index=0;
  int index_pk_eq_z;
  int index_eq;
  int true_background_verbose;
  int true_thermodynamics_verbose;
  double true_w0_fld;
  double true_wa_fld;
  double * z;

  /** - store the true cosmological parameters (w0, wa) somwhere before using temporarily some fake ones in this function */

  true_background_verbose = pba->background_verbose;
  true_thermodynamics_verbose = pth->thermodynamics_verbose;
  true_w0_fld = pba->w0_fld;
  true_wa_fld = pba->wa_fld;

  /** - the fake calls of the background and thermodynamics module will be done in non-verbose mode */

  pba->background_verbose = 0;
  pth->thermodynamics_verbose = 0;

  /** - allocate indices and arrays for storing the results */

  pnl->pk_eq_tau_size = 10;
  class_alloc(pnl->pk_eq_tau,pnl->pk_eq_tau_size*sizeof(double),errmsg);
  class_alloc(z,pnl->pk_eq_tau_size*sizeof(double),errmsg);

  index_eq = 0;
  class_define_index(pnl->index_pk_eq_w,_TRUE_,index_eq,1);
  class_define_index(pnl->index_pk_eq_Omega_m,_TRUE_,index_eq,1);
  pnl->pk_eq_size = index_eq;
  class_alloc(pnl->pk_eq_w_and_Omega,pnl->pk_eq_tau_size*pnl->pk_eq_size*sizeof(double),errmsg);
  class_alloc(pnl->pk_eq_ddw_and_ddOmega,pnl->pk_eq_tau_size*pnl->pk_eq_size*sizeof(double),errmsg);

  /** - call the background module in order to fill a table of tau_i[z_i] */

  class_call(background_init(ppr,pba), pba->error_message, errmsg);
  for (index_pk_eq_z=0; index_pk_eq_z<pnl->pk_eq_tau_size; index_pk_eq_z++) {
    z[index_pk_eq_z] = exp(log(1.+ppr->pk_eq_z_max)/(pnl->pk_eq_tau_size-1)*index_pk_eq_z)-1.;
    class_call(background_tau_of_z(pba,z[index_pk_eq_z],&tau_of_z),
               pba->error_message, errmsg);
    pnl->pk_eq_tau[index_pk_eq_z] = tau_of_z;
  }
  class_call(background_free_noinput(pba), pba->error_message, errmsg);

  /** - loop over z_i values. For each of them, we will call the
      background and thermodynamics module for fake models. The goal is
      to find, for each z_i, and effective w0_eff[z_i] and
      Omega_m_eff{z_i], such that: the true model with (w0,wa) and the
      equivalent model with (w0_eff[z_i],0) have the same conformal
      distance between z_i and z_recombination, namely chi = tau[z_i] -
      tau_rec. It is thus necessary to call both the background and
      thermodynamics module for each fake model and to re-compute
      tau_rec for each of them. Once the eqauivalent model is found we
      compute and store Omega_m_effa(z_i) of the equivalent model */

  for (index_pk_eq_z=0; index_pk_eq_z<pnl->pk_eq_tau_size; index_pk_eq_z++) {

    if (input_verbose > 2)
      printf("    * computing Pk_equal parameters at z=%e\n",z[index_pk_eq_z]);

    /* get chi = (tau[z_i] - tau_rec) in true model */

    pba->w0_fld = true_w0_fld;
    pba->wa_fld = true_wa_fld;

    class_call(background_init(ppr,pba), pba->error_message, errmsg);
    class_call(thermodynamics_init(ppr,pba,pth), pth->error_message, errmsg);

    delta_tau = pnl->pk_eq_tau[index_pk_eq_z] - pth->tau_rec;

    /* launch iterations in order to coverge to effective model with wa=0 but the same chi = (tau[z_i] - tau_rec) */

    pba->wa_fld=0.;

    do {
      class_call(background_free_noinput(pba), pba->error_message, errmsg);
      class_call(thermodynamics_free(pth), pth->error_message, errmsg);

      class_call(background_init(ppr,pba), pba->error_message, errmsg);
      class_call(background_tau_of_z(pba,z[index_pk_eq_z],&tau_of_z), pba->error_message, errmsg);
      class_call(thermodynamics_init(ppr,pba,pth), pth->error_message, errmsg);

      delta_tau_eq = tau_of_z - pth->tau_rec;

      error = 1.-delta_tau_eq/delta_tau;
      pba->w0_fld = pba->w0_fld*pow(1.+error,10.);

    }
    while(fabs(error) > ppr->pk_eq_tol);

    /* Equivalent model found. Store w0(z) in that model. Find Omega_m(z) in that model and store it. */

    pnl->pk_eq_w_and_Omega[pnl->pk_eq_size*index_pk_eq_z+pnl->index_pk_eq_w] = pba->w0_fld;

    class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
    class_call(background_at_tau(pba,
                                 tau_of_z,
                                 pba->long_info,
                                 pba->inter_normal,
                                 &last_index,
                                 pvecback),
               pba->error_message, errmsg);
    pnl->pk_eq_w_and_Omega[pnl->pk_eq_size*index_pk_eq_z+pnl->index_pk_eq_Omega_m] = pvecback[pba->index_bg_Omega_m];
    free(pvecback);

    class_call(background_free_noinput(pba), pba->error_message, errmsg);
    class_call(thermodynamics_free(pth), pth->error_message, errmsg);

  }

  /** - restore cosmological parameters (w0, wa) to their true values before main call to CLASS modules */

  pba->background_verbose = true_background_verbose;
  pth->thermodynamics_verbose = true_thermodynamics_verbose;
  pba->w0_fld = true_w0_fld;
  pba->wa_fld = true_wa_fld;

  /* in verbose mode, report the results */

  if (input_verbose > 1) {

    fprintf(stdout,"    Effective parameters for Pk_equal:\n");

    for (index_pk_eq_z=0; index_pk_eq_z<pnl->pk_eq_tau_size; index_pk_eq_z++) {

      fprintf(stdout,"    * at z=%e, tau=%e w=%e Omega_m=%e\n",
              z[index_pk_eq_z],
              pnl->pk_eq_tau[index_pk_eq_z],
              pnl->pk_eq_w_and_Omega[pnl->pk_eq_size*index_pk_eq_z+pnl->index_pk_eq_w],
              pnl->pk_eq_w_and_Omega[pnl->pk_eq_size*index_pk_eq_z+pnl->index_pk_eq_Omega_m]
              );
    }
  }

  free(z);

  /** - spline the table for later interpolation */

  class_call(array_spline_table_lines(
                                      pnl->pk_eq_tau,
                                      pnl->pk_eq_tau_size,
                                      pnl->pk_eq_w_and_Omega,
                                      pnl->pk_eq_size,
                                      pnl->pk_eq_ddw_and_ddOmega,
                                      _SPLINE_NATURAL_,
                                      errmsg),
             errmsg,errmsg);

  return _SUCCESS_;

}
