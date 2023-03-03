#include <stdio.h>
#include <math.h>
#include "helmholtz.h"

extern int     imax;
extern int     jmax;

extern double *helmholtz_table;
extern double *helmholtz_dd;
extern double *helmholtz_dt;
extern double *helm_dens;
extern double *helm_temp;

void Helmholtz_eos_C_ReadTable( char *helmeos_table_name )
{

   FILE *FILE = fopen( helmeos_table_name, "r" );


// allocate memory for tables
   if (  ! ( helmholtz_table = (double*)malloc(imax*jmax*21       *sizeof(double)) )  )
      fprintf( stderr, "cannot allocate memory for EOS table !!\n" );

   if (  ! ( helmholtz_dd    = (double*)malloc((imax-1)*5         *sizeof(double)) )  )
      fprintf( stderr, "cannot allocate memory for EOS table !!\n" );

   if (  ! ( helmholtz_dt    = (double*)malloc((jmax-1)*5         *sizeof(double)) )  )
      fprintf( stderr, "cannot allocate memory for EOS table !!\n" );

   if (  ! ( helm_dens       = (double*)malloc(imax               *sizeof(double)) )  )
      fprintf( stderr, "cannot allocate memory for EOS table !!\n" );

   if (  ! ( helm_temp       = (double*)malloc(jmax               *sizeof(double)) )  )
      fprintf( stderr, "cannot allocate memory for EOS table !!\n" );



   const double tlo   = 3.0;
   const double thi   = 13.0;
   const double tstp  = (thi - tlo) / (jmax - 1);
   const double tstpi = 1.0 / tstp;
   const double dlo   = -12.0;
   const double dhi   = 15.0;
   const double dstp  = (dhi - dlo) / (imax - 1);
   const double dstpi = 1.0 / dstp;


//  Set up arrays to store the Helmholtz free energy and its derivatives
   double *f      = (double *)malloc( imax * jmax * sizeof(double) );
   double *fd     = (double *)malloc( imax * jmax * sizeof(double) );
   double *ft     = (double *)malloc( imax * jmax * sizeof(double) );
   double *fdd    = (double *)malloc( imax * jmax * sizeof(double) );
   double *ftt    = (double *)malloc( imax * jmax * sizeof(double) );
   double *fdt    = (double *)malloc( imax * jmax * sizeof(double) );
   double *fddt   = (double *)malloc( imax * jmax * sizeof(double) );
   double *fdtt   = (double *)malloc( imax * jmax * sizeof(double) );
   double *fddtt  = (double *)malloc( imax * jmax * sizeof(double) );
   double *dpdf   = (double *)malloc( imax * jmax * sizeof(double) );
   double *dpdfd  = (double *)malloc( imax * jmax * sizeof(double) );
   double *dpdft  = (double *)malloc( imax * jmax * sizeof(double) );
   double *dpdfdt = (double *)malloc( imax * jmax * sizeof(double) );
   double *ef     = (double *)malloc( imax * jmax * sizeof(double) );
   double *efd    = (double *)malloc( imax * jmax * sizeof(double) );
   double *eft    = (double *)malloc( imax * jmax * sizeof(double) );
   double *efdt   = (double *)malloc( imax * jmax * sizeof(double) );
   double *xf     = (double *)malloc( imax * jmax * sizeof(double) );
   double *xfd    = (double *)malloc( imax * jmax * sizeof(double) );
   double *xft    = (double *)malloc( imax * jmax * sizeof(double) );
   double *xfdt   = (double *)malloc( imax * jmax * sizeof(double) );

   double *dd_sav   = (double *)malloc( (imax-1) * sizeof(double) );
   double *dd2_sav  = (double *)malloc( (imax-1) * sizeof(double) );
   double *ddi_sav  = (double *)malloc( (imax-1) * sizeof(double) );
   double *dd2i_sav = (double *)malloc( (imax-1) * sizeof(double) );
   double *dd3i_sav = (double *)malloc( (imax-1) * sizeof(double) );

   double *dt_sav   = (double *)malloc( (jmax-1) * sizeof(double) );
   double *dt2_sav  = (double *)malloc( (jmax-1) * sizeof(double) );
   double *dti_sav  = (double *)malloc( (jmax-1) * sizeof(double) );
   double *dt2i_sav = (double *)malloc( (jmax-1) * sizeof(double) );
   double *dt3i_sav = (double *)malloc( (jmax-1) * sizeof(double) );


   for (int j=0; j<jmax; j++)
   {
      double tsav = tlo + j*tstp;
      helm_temp[j] = pow( 10.0, tsav );
      for (int i=0; i<imax; i++)
      {
         double dsav = dlo + i*dstp;
         helm_dens[i] = pow( 10.0, dsav );
         fscanf( FILE, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                 &f[i*jmax+j],    &fd[i*jmax+j],   &ft[i*jmax+j],
                 &fdd[i*jmax+j],  &ftt[i*jmax+j],  &fdt[i*jmax+j],
                 &fddt[i*jmax+j], &fdtt[i*jmax+j], &fddtt[i*jmax+j] );
      }
   }


   for (int j=0; j<jmax; j++)
   for (int i=0; i<imax; i++)
      fscanf( FILE, "%lf %lf %lf %lf",
              &dpdf[i*jmax+j], &dpdfd[i*jmax+j], &dpdft[i*jmax+j], &dpdfdt[i*jmax+j] );

   for (int j=0; j<jmax; j++)
   for (int i=0; i<imax; i++)
      fscanf( FILE, "%lf %lf %lf %lf",
              &ef[i*jmax+j], &efd[i*jmax+j], &eft[i*jmax+j], &efdt[i*jmax+j] );

   for (int j=0; j<jmax; j++)
   for (int i=0; i<imax; i++)
      fscanf( FILE, "%lf %lf %lf %lf",
              &xf[i*jmax+j], &xfd[i*jmax+j], &xft[i*jmax+j], &xfdt[i*jmax+j] );


// calculate the derivatives
   for (int i=0; i<imax-1; i++)
   {
      double dd   = helm_dens[i+1] - helm_dens[i];
      double dd2  = dd * dd;
      double ddi  = 1.0/dd;
      double dd2i = 1.0/dd2;
      double dd3i = dd2i * ddi;

      dd_sav  [i] = dd;
      dd2_sav [i] = dd2;
      ddi_sav [i] = ddi;
      dd2i_sav[i] = dd2i;
      dd3i_sav[i] = dd3i;
   }

   for (int j=0; j<jmax-1; j++)
   {
      double dth  = helm_temp[j+1] - helm_temp[j];
      double dt2  = dth * dth;
      double dti  = 1.0/dth;
      double dt2i = 1.0/dt2;
      double dt3i = dt2i * dti;

      dt_sav  [j] = dth;
      dt2_sav [j] = dt2;
      dti_sav [j] = dti;
      dt2i_sav[j] = dt2i;
      dt3i_sav[j] = dt3i;
   }


// store the variables
   for (int i=0; i<imax; i++)
   {
      for (int j=0; j<jmax; j++)
      {
         helmholtz_table[ ( i*jmax*21) + (j*21) + 0  ] = f     [i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 1  ] = fd    [i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 2  ] = ft    [i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 3  ] = fdd   [i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 4  ] = ftt   [i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 5  ] = fdt   [i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 6  ] = fddt  [i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 7  ] = fdtt  [i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 8  ] = fddtt [i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 9  ] = dpdf  [i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 10 ] = dpdfd [i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 11 ] = dpdft [i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 12 ] = dpdfdt[i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 13 ] = ef    [i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 14 ] = efd   [i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 15 ] = eft   [i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 16 ] = efdt  [i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 17 ] = xf    [i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 18 ] = xfd   [i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 19 ] = xft   [i*jmax + j];
         helmholtz_table[ ( i*jmax*21) + (j*21) + 20 ] = xfdt  [i*jmax + j];
      }
   }

   for (int i=0; i<imax-1; i++)
   {
      helmholtz_dd[ i*5 + 0 ] = dd_sav[i];
      helmholtz_dd[ i*5 + 1 ] = dd2_sav[i];
      helmholtz_dd[ i*5 + 2 ] = ddi_sav[i];
      helmholtz_dd[ i*5 + 3 ] = dd2i_sav[i];
      helmholtz_dd[ i*5 + 4 ] = dd3i_sav[i];
   }

   for (int j=0; j<jmax-1; j++)
   {
      helmholtz_dt[ j*5 + 0 ] = dt_sav[j];
      helmholtz_dt[ j*5 + 1 ] = dt2_sav[j];
      helmholtz_dt[ j*5 + 2 ] = dti_sav[j];
      helmholtz_dt[ j*5 + 3 ] = dt2i_sav[j];
      helmholtz_dt[ j*5 + 4 ] = dt3i_sav[j];
   }



// free memory
   free(f     ); free(fd  ); free(ft   ); free(fdd  );
   free(ftt   ); free(fdt ); free(fddt ); free(fdtt );
   free(fddtt ); free(dpdf); free(dpdfd); free(dpdft);
   free(dpdfdt); free(ef  ); free(efd  ); free(eft  );
   free(efdt  ); free(xf  ); free(xfd  ); free(xft  );
   free(xfdt  ); 
   
   free(dd_sav); free(dd2_sav ); free(ddi_sav); free(dd2i_sav); free(dd3i_sav);
   free(dt_sav); free(dt2_sav ); free(dti_sav); free(dt2i_sav); free(dt3i_sav);

// close the file
   fclose( FILE );

}