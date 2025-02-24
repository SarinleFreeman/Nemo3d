25c25
< $Source: /repo/eigen/src/patch.c,v $
---
> $Source: /repo/eigen/src/patch.c,v $
407a408
>     int		*n_nosuccess,	/*   out number no success		*/
419,420c420,421
< /*const nml_dscalar	bias = 1e-80;*/	/* PLEASE EXPLAIN ME		*/
< /*const nml_dscalar	gain = 1e-8;*/	/* PLEASE EXPLAIN ME		*/
---
>   const nml_dscalar	bias = 1e-80;	/* PLEASE EXPLAIN ME			*/
>   const nml_dscalar	gain = 1e-8;	/* PLEASE EXPLAIN ME			*/
432,433c433
< /*  nml_dscalar	abstol = 1.0e-13;*/	/* eigenvalue absolute tolerance*/
<     nml_dscalar	abstol = 0.01*tolerance;/* eigenvalue absolute tolerance*/
---
>     nml_dscalar	abstol = 1.0e-13;	/* eigenvalue absolute tolerance*/
459a460
>   *n_nosuccess = 0;
464,465c465,466
< /*    if ((0 < i) && (fabs(eval[i-1] - eval[i]) < (bias + eval[i])*gain))    */
<       if ((0 < i) && (fabs(eval[i-1] - eval[i]) < resolution))
---
>       if ((0 < i) &&
>         (fabs(eval[i-1] - eval[i]) < (bias + eval[i])*gain))
502a504
> 	      *n_nosuccess = 0;
520,521c522
<     for (iteration = 0; iteration < iterations; ++iteration) {
<       /*
---
>     for (iteration = 0; iteration < iterations; ++iteration)
526,530d526
<       */
<       if ((emin/MTRX_SCALE < eval[iteration])
<        && (eval[iteration] < emax/MTRX_SCALE)
<        && ((0 == iteration)
< 	|| (resolution < fabs((eval[iteration-1] - eval[iteration]))))) {
534d529
<       }
536a532,533
>   ++(*n_nosuccess);
> 
551a549
>     int		 nsuc,		/* in maximum number no sucess		*/
555a554
>     nml_dscalar	 shift,		/* in eigenvalue shift			*/
568a568
>   int	n_nosuccess = 0;	/* number of no success			*/
639c639
< 	    &n_old, eigv_old, fp_trace, verbose, fp_instant,
---
> 	    &n_nosuccess, &n_old, eigv_old, fp_trace, verbose, fp_instant,
643c643,644
<       no_convergence = ((iterations + 1 < imax) && (eigenvalues < requested));
---
>       no_convergence = ((iterations + 1 < imax)
> 		     && (eigenvalues < requested) && (n_nosuccess < nsuc));
680a682
>     nml_extent	 nsuc,		/* in maximum number no sucess		*/
684a687
>     nml_dscalar	 shift,		/* in eigenvalue shift			*/
697,698c700,701
<       ConvCheckStartIter, ConvCheckSkipRate, imax,
<       emin, emax, tolerance, resolution, matmul, argument,
---
>       ConvCheckStartIter, ConvCheckSkipRate, imax, nsuc,
>       emin, emax, tolerance, resolution, shift, matmul, argument,
959a963
>     int		 nsuc,		/* in maximum number no sucess		*/
970a975
>     nml_dscalar	 shift,		/* in eigenvalue shift			*/
1039,1040c1044,1045
< 	ConvCheckStartIter, ConvCheckSkipRate, imax,
< 	emin, emax, tolerance, 100.0*tolerance, matmul, (const int**)argument,
---
> 	ConvCheckStartIter, ConvCheckSkipRate, imax, nsuc,
> 	emin, emax, tolerance, 100.0*tolerance, shift, matmul, (const int**)argument,
