template<typename T>
int GMRESLinSolver<T>::
operator()(JacobiImpApply MV, JacobiImpApplyPrecond PC, T *computed_solution, T *rhs, 
        T reltol, T abstol, size_t N, int maxIters, int restartIters, int vid) const
{
  /*---------------------------------------------------------------------------
  * Allocate storage for the ?par parameters and the solution/rhs vectors
  *---------------------------------------------------------------------------*/
  MKL_INT ipar[GMRESParSize];
  double dpar[GMRESParSize];
  int tmp_N = N * (2 * restartIters + 1) + (restartIters * (restartIters + 9)) / 2 + 1;
  /*
  double tmp[tmp_N];
  double b[N];
  double residual[N];
  */
  double *tmp      = new double[tmp_N];
  double *b        = new double[N];
  double *residual = new double[N];
  /*---------------------------------------------------------------------------
  * Some additional variables to use with the RCI (P)FGMRES solver
  *---------------------------------------------------------------------------*/
  MKL_INT itercount;
  MKL_INT RCI_request, i, ivar;
  double dvar;
  
  COUT("Solving the linear system using RCI FGMRES solver\n");
  
  /*---------------------------------------------------------------------------
  * Initialize variables and the right hand side through matrix-vector product
  *---------------------------------------------------------------------------*/
  ivar = N;
  /*---------------------------------------------------------------------------
  * Save the right-hand side in vector b for future use
  *---------------------------------------------------------------------------*/
  i = 1;
  dcopy (&ivar, rhs, &i, b, &i);
  /*---------------------------------------------------------------------------
  * Initialize the solver
  *---------------------------------------------------------------------------*/
  dfgmres_init (&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp);
  if (RCI_request != 0)
    goto FAILED;
  /*---------------------------------------------------------------------------
  * Set the desired parameters:
  * LOGICAL parameters:
  * set the maximum number of iterations
  * do residual stopping test
  * do not request for the user defined stopping test
  * do the check of the norm of the next generated vector automatically
  * do the restart after restartIters iterations
  * DOUBLE PRECISION parameters
  * set the relative tolerance to tol instead of default value 1.0D-6
  * set the absulote tolerance to 1.0E-07
  *---------------------------------------------------------------------------*/
  ipar[14] = restartIters;
  ipar[4] = maxIters;
  ipar[7] = 0;
  ipar[10] = 1;
  //ipar[8] = 1;
  //ipar[9] = 0;
  //ipar[11] = 1;
  dpar[0] = reltol;
  dpar[1] = abstol;
  /*---------------------------------------------------------------------------
  * Check the correctness and consistency of the newly set parameters
  *---------------------------------------------------------------------------*/
  dfgmres_check (&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp);
  if (RCI_request != 0)
    goto FAILED;
  /*---------------------------------------------------------------------------
  * Print the info about the RCI FGMRES method
  *---------------------------------------------------------------------------*/
  /*
  printf ("Some info about the current run of RCI FGMRES method:\n\n");
  if (ipar[7])
    {
      printf ("As ipar[7]=%d, the automatic test for the maximal number of ", ipar[7]);
      printf ("iterations will be\nperformed\n");
    }
  else
    {
      printf ("As ipar[7]=%d, the automatic test for the maximal number of ", ipar[7]);
      printf ("iterations will be\nskipped\n");
    }
  printf ("+++\n");
  if (ipar[8])
    {
      printf ("As ipar[8]=%d, the automatic residual test will be performed\n", ipar[8]);
    }
  else
    {
      printf ("As ipar[8]=%d, the automatic residual test will be skipped\n", ipar[8]);
    }
  printf ("+++\n");
  if (ipar[9])
    {
      printf ("As ipar[9]=%d, the user-defined stopping test will be ", ipar[9]);
      printf ("requested via\nRCI_request=2\n");
    }
  else
    {
      printf ("As ipar[9]=%d, the user-defined stopping test will not be ", ipar[9]);
      printf ("requested, thus,\nRCI_request will not take the value 2\n");
    }
  printf ("+++\n");
  if (ipar[10])
    {
      printf ("As ipar[10]=%d, the Preconditioned FGMRES iterations will be ", ipar[10]);
      printf ("performed, thus,\nthe preconditioner action will be requested via");
      printf ("RCI_request=3\n");
    }
  else
    {
      printf ("As ipar[10]=%d, the Preconditioned FGMRES iterations will not ", ipar[10]);
      printf ("be performed,\nthus, RCI_request will not take the value 3\n");
    }
  printf ("+++\n");
  if (ipar[11])
    {
      printf ("As ipar[11]=%d, the automatic test for the norm of the next ", ipar[11]);
      printf ("generated vector is\nnot equal to zero up to rounding and ");
      printf ("computational errors will be performed,\nthus, RCI_request will not take the value 4\n");
    }
  else
    {
      printf ("As ipar[11]=%d, the automatic test for the norm of the next ", ipar[11]);
      printf ("generated vector is\nnot equal to zero up to rounding and ");
      printf ("computational errors will be skipped,\nthus, the user-defined test ");
      printf ("will be requested via RCI_request=4\n");
    }
  printf ("+++\n\n");
  */
  /*---------------------------------------------------------------------------
  * Compute the solution by RCI (P)FGMRES solver without preconditioning
  * Reverse Communication starts here
  *---------------------------------------------------------------------------*/
ONE:dfgmres (&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp);
  /*---------------------------------------------------------------------------
  * If RCI_request=0, then the solution was found with the required precision
  *---------------------------------------------------------------------------*/
  if (RCI_request == 0)
    goto COMPLETE;
  /*---------------------------------------------------------------------------
  * If RCI_request=1, then compute the vector A*tmp[ipar[21]-1]
  * and put the result in vector tmp[ipar[22]-1]
  *---------------------------------------------------------------------------
  * NOTE that ipar[21] and ipar[22] contain FORTRAN style addresses,
  * therefore, in C code it is required to subtract 1 from them to get C style
  * addresses
  *---------------------------------------------------------------------------*/
  if (RCI_request == 1)
    {
      MV(this, &tmp[ipar[21] - 1], &tmp[ipar[22] - 1], vid);
      goto ONE;
    }
  if (RCI_request == 2)
    {
      ipar[12] = 1;
      dfgmres_get (&ivar, computed_solution, b, &RCI_request, ipar, dpar, tmp, &itercount);
      MV(this, b, residual, vid);
      dvar = -1.0E0;
      i = 1;
      daxpy (&ivar, &dvar, rhs, &i, residual, &i);
      dvar = dnrm2 (&ivar, residual, &i);
      
      //COUT("rank: "<<myrank<<". Current iteration: "<<itercount<<". residual: "<<dvar<<"\n");
      
      if (dvar <= dpar[3])
        goto COMPLETE;
      if (itercount > ipar[4])
        goto FAILED;
      else
        goto ONE;
    }
  if (RCI_request == 3)
    {
      PC(this, &tmp[ipar[21] - 1], &tmp[ipar[22] - 1]);
      goto ONE;
    }
  if (RCI_request == 4)
    {
      if (dpar[6] < 1.0E-12)
        goto COMPLETE;
      else
        goto ONE;
    }
  /*---------------------------------------------------------------------------
  * If RCI_request=anything else, then dfgmres subroutine failed
  * to compute the solution vector: computed_solution[N]
  *---------------------------------------------------------------------------*/
  else
    {
      goto FAILED;
    }
  /*---------------------------------------------------------------------------
  * Reverse Communication ends here
  * Get the current iteration number and the FGMRES solution (DO NOT FORGET to
  * call dfgmres_get routine as computed_solution is still containing
  * the initial guess!)
  *---------------------------------------------------------------------------*/
COMPLETE:ipar[12] = 0;
  dfgmres_get (&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp, &itercount);
  /*---------------------------------------------------------------------------
  * Print solution vector: computed_solution[N] and the number of iterations: itercount
  *--------------------------------------------------------------------------- */
  COUT("The system has been solved. Number of iterations: "<<itercount<<".\n");
  /*-------------------------------------------------------------------------*/
  /* Release internal Intel(R) MKL memory that might be used for computations         */
  /* NOTE: It is important to call the routine below to avoid memory leaks   */
  /* unless you disable Intel(R) MKL Memory Manager                                   */
  /*-------------------------------------------------------------------------*/
  delete[] tmp;
  delete[] b;
  delete[] residual;
  MKL_Free_Buffers ();
  return 0;
  /*-------------------------------------------------------------------------*/
  /* Release internal Intel(R) MKL memory that might be used for computations         */
  /* NOTE: It is important to call the routine below to avoid memory leaks   */
  /* unless you disable Intel(R) MKL Memory Manager                                   */
  /*-------------------------------------------------------------------------*/

FAILED:COUT("Solving the system FAILED as the solver has returned the ERROR code "<<RCI_request<<"\n");
  delete[] tmp;
  delete[] b;
  delete[] residual;
  MKL_Free_Buffers ();
  return 1;
}

template<typename T>
Error_t GMRESLinSolver<T>::
SetContext(const void *ctx)
{
    ctx_ = ctx;
    return ErrorEvent::Success;
}

template<typename T>
Error_t GMRESLinSolver<T>::
Context(const void **ctx) const
{
    *ctx = ctx_;
    return ErrorEvent::Success;
}
