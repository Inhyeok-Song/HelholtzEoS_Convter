#include <stdlib.h>
#include <stdio.h>
#include "helmholtz.h"



extern double *logrho;
extern double *logtemp;
extern double *yes;



void set_bins( int rho_bs, int temp_bs, int ye_bs, double eos_rhomin, double eos_rhomax,
               double eos_tempmin, double eos_tempmax, double eos_yemin, double eos_yemax )
{

// Allocate memory for logrho, logtemp, and ye
   if (  ! ( logrho  = (double*)malloc( rho_bs  * sizeof(double) ) )  )
      fprintf( stderr, "Cannot allocate memery for logtemp\n" );
   if (  ! ( logtemp = (double*)malloc( temp_bs * sizeof(double) ) )  )
      fprintf(stderr, "Cannot allocate memery for logtemp\n");
   if (  ! ( yes     = (double*)malloc( ye_bs   * sizeof(double) ) )  )
      fprintf(stderr, "Cannot allocate memery for logtemp\n");


// Set logrho
   fprintf( stdout, "*******************************\n"  );
   fprintf( stdout, "Setting rho\n"                      );
   fprintf( stdout, "rho_binsize:  %d\n",    rho_bs      );
   fprintf( stdout, "rho_min:      %9.3E\n", eos_rhomin  );
   fprintf( stdout, "rho_max:      %9.3E\n", eos_rhomax  );

   double dlrho = LOG10( eos_rhomax/eos_rhomin ) / ( rho_bs-1 );
   for (int i=0; i<rho_bs-1; i++)
      logrho[i] = LOG10( eos_rhomin ) + dlrho*i;
   logrho[rho_bs-1] = LOG10( eos_rhomax );


// Set logtemp
   fprintf( stdout, "*******************************\n"  );
   fprintf( stdout, "Setting temp\n"                     );
   fprintf( stdout, "temp_binsize: %d\n",    temp_bs     );
   fprintf( stdout, "temp_min:     %9.3E\n", eos_tempmin );
   fprintf( stdout, "temp_max:     %9.3E\n", eos_tempmax );
   double dltemp = LOG10( eos_tempmax/eos_tempmin ) / ( temp_bs-1 );
   for (int i=0; i<temp_bs-1; i++)
      logtemp[i] = LOG10(eos_tempmin) + dltemp * i;
   logtemp[temp_bs - 1] = LOG10(eos_tempmax);


// Set Ye
   fprintf( stdout, "*******************************\n"  );
   fprintf( stdout, "Setting ye\n"                       );
   fprintf( stdout, "ye_binsize:   %d\n",    ye_bs       );
   fprintf( stdout, "ye_min:       %9.3E\n", eos_yemin   );
   fprintf( stdout, "ye_max:       %9.3E\n", eos_yemax   );
   double dye = ( eos_yemax-eos_yemin ) / ( ye_bs-1 );
   for (int i=0; i<ye_bs-1; i++)
      yes[i] = eos_yemin + dye * i;
   yes[ye_bs - 1] = eos_yemax;


   return;
}
