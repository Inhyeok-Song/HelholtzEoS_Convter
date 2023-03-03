#include <stdio.h>
#include "helmholtz.h"



int     imax  = 541;
int     jmax  = 201;
int     nrho  = 100;
int     ntemp = 100;
int     nye   = 100;

double *alltables;
double *logrho;
double *logtemp;
double *yes;
double *helm_dens;
double *helm_temp;
double *helmholtz_table;
double *helmholtz_dd;
double *helmholtz_dt;



int main( int argc, char* argv[] ) {

// read helmholtz table
   fprintf( stdout, "Read Helmholtz EOS table\n" );
   Helmholtz_eos_C_ReadTable( "../helm_table.dat" );


// set eos limits
   const double eos_rhomin  = 1e-3;
   const double eos_rhomax  = 1e+6;

   const double eos_tempmin = 1e-4;
   const double eos_tempmax = 1e-1;

   const double eos_yemin   = 0.4;
   const double eos_yemax   = 0.85;

   set_bins( nrho, ntemp, nye, eos_rhomin, eos_rhomax,
             eos_tempmin, eos_tempmax, eos_yemin, eos_yemax );


// write HDF5 table
   fprintf( stdout, "Read Helmholtz EOS table\n" );
   write_eos_table( "../Helmholtz.h5" );


   return 0;
}
