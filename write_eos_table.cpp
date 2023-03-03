#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "helmholtz.h"



// Catch HDF5 errors
#define HDF5_ERROR( fn_call )                                        \
   do                                                                \
   {                                                                 \
      int _error_code = fn_call;                                     \
      if ( _error_code<0 )                                           \
      {                                                              \
         fprintf( stderr,                                            \
                  "HDF5 call '%s' returned error code %d",           \
                   #fn_call, _error_code );                          \
         abort();                                                    \
      }                                                              \
   } while (0)



void write_eos_table( char *helmeos_table_name )
{


   fprintf( stdout,"*********************************\n" );
   fprintf( stdout,"Writing Helmtholz EoS table file:\n" );
   fprintf( stdout,"%s\n", helmeos_table_name            );
   fprintf( stdout,"*********************************\n" );


   hid_t     file;
   hid_t     dataset;
   hid_t     space;
   hsize_t   dims[3] = {nye, ntemp, nrho};


   // Allocating eos table
   double ****eos_table = NULL;  // Write buffer
   if (  ! ( eos_table = (double****)malloc(nye * sizeof(double***)) )  )
   {
      fprintf( stderr, "Cannot allocate memory for EOS table\n" );
      abort();
   }
   for (int i=0; i<nye; i++)
   {
      if (  ! ( eos_table[i] = (double***)malloc(ntemp * sizeof(double**)) )  )
      {
         fprintf( stderr, "Cannot allocate memory for EOS table\n" );
         abort();
      }
      for (int j=0; j<ntemp; j++)
      {
         if (  ! ( eos_table[i][j] = (double**)malloc(nrho * sizeof(double*)) )  )
         {
            fprintf( stderr, "Cannot allocate memory for EOS table\n" );
            abort();
         }
         for (int k=0; k<nrho; k++)
         {
            if (  ! ( eos_table[i][j][k] = (double*)malloc((NTABLES) * sizeof(double)) )  )
            {
               fprintf( stderr, "Cannot allocate memory for EOS table\n" );
               abort();
            }
         }
      }
   }


   double *wdata;
   if (  ! ( wdata = (double*)malloc(nye * ntemp * nrho * sizeof(double)) )  )
   {
      fprintf(stderr, "Cannot allocate memory for EOS table\n");
      abort();
   }


   // Find variables. The last row is the elements in the dataspace, i, j
   // and k are the elements within the array datatype
   for (int i=0; i<nye;   i++)
   for (int j=0; j<ntemp; j++)
   for (int k=0; k<nrho;  k++) {

      const int  NTarget = 4;
      double In[3], Out[NTarget+1];
      int  TargetIdx[NTarget] = { NUC_VAR_IDX_PRES, NUC_VAR_IDX_ENGY, NUC_VAR_IDX_ENTR, NUC_VAR_IDX_CSQR };

      In[0] = POW( 10.0, logrho[k]  );
      In[1] = POW( 10.0, logtemp[j] );
      In[2] = yes[i];

      Helmholtz_eos( Out, In, NTarget, TargetIdx, imax, jmax, helm_dens, helm_temp,
                     helmholtz_table, helmholtz_dd, helmholtz_dt );


      eos_table[i][j][k][NUC_VAR_IDX_PRES] = Out[NUC_VAR_IDX_PRES];
      eos_table[i][j][k][NUC_VAR_IDX_ENGY] = Out[NUC_VAR_IDX_ENGY];
      eos_table[i][j][k][NUC_VAR_IDX_ENTR] = Out[NUC_VAR_IDX_ENTR];
      eos_table[i][j][k][NUC_VAR_IDX_CSQR] = Out[NUC_VAR_IDX_CSQR];
   }



// Create a new file using the default properties
   HDF5_ERROR(  file = H5Fcreate( helmeos_table_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT )  );


// Use these defines to easily write a lot of variables in the same way
#define WRITE_EOSTABLE_HDF5( NAME, OFF  )                                                             \
   do                                                                                                 \
   {                                                                                                  \
      for (int i=0; i<nye; i++)                                                                       \
         for (int j=0; j<ntemp; j++)                                                                  \
            for (int k=0; k<nrho; k++)                                                                \
               wdata[i*ntemp*nrho + j*nrho + k] = eos_table[i][j][k][OFF];                            \
               HDF5_ERROR(  H5LTmake_dataset( file, NAME, 3, dims, H5T_NATIVE_DOUBLE, wdata )  );     \
   } while (0)

#define WRITE_EOS_HDF5( NAME, var )                                                                   \
   do                                                                                                 \
   {                                                                                                  \
      hsize_t dims[1] = {1};                                                                          \
      HDF5_ERROR(  space = H5Screate_simple( 1, dims, NULL )  );                                      \
      HDF5_ERROR(  dataset = H5Dcreate( file, NAME, H5T_INTEL_I32, space, H5P_DEFAULT,                \
                   H5P_DEFAULT, H5P_DEFAULT )  );                                                     \
      HDF5_ERROR(  H5Dwrite( dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, var )  );        \
      HDF5_ERROR(  H5Dclose( dataset )  );                                                            \
      HDF5_ERROR(  H5Sclose( space )  );                                                              \
   } while (0)

#define WRITE_EOS_INIT_HDF5( NAME, nvar, var )                                                        \
   do                                                                                                 \
   {                                                                                                  \
      hsize_t dims[1] = {nvar};                                                                       \
      HDF5_ERROR(  space = H5Screate_simple( 1, dims, NULL )  );                                      \
      HDF5_ERROR(  dataset = H5Dcreate( file, NAME, H5T_IEEE_F64LE, space, H5P_DEFAULT,               \
                   H5P_DEFAULT, H5P_DEFAULT )  );                                                     \
      HDF5_ERROR(  H5Dwrite( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, var )  );     \
      HDF5_ERROR(  H5Dclose( dataset )  );                                                            \
      HDF5_ERROR(  H5Sclose( space ) );                                                               \
   } while (0)


   // Write size of tables
   WRITE_EOS_HDF5( "pointsrho",  &nrho  );
   WRITE_EOS_HDF5( "pointstemp", &ntemp );
   WRITE_EOS_HDF5( "pointsye",   &nye   );

   // writes alltables
   WRITE_EOSTABLE_HDF5( "logpress",  0 );
   WRITE_EOSTABLE_HDF5( "logenergy", 1 );
   WRITE_EOSTABLE_HDF5( "entropy",   2 );
   WRITE_EOSTABLE_HDF5( "cs2",       3 );

   // Write additional tables and variable
   WRITE_EOS_INIT_HDF5( "logrho",  nrho,  logrho  );
   WRITE_EOS_INIT_HDF5( "logtemp", ntemp, logtemp );
   WRITE_EOS_INIT_HDF5( "ye",      nye,   yes     );

   // Close and release resources
   HDF5_ERROR(  H5Fclose( file )  );


   // Free memory
   free( wdata );
   wdata = NULL;

   for (int i=0; i<nye; i++)
   {
      for (int j=0; j<ntemp; j++)
      {
         for (int k=0; k<nrho; k++)
         {
            free( eos_table[i][j][k] );
            eos_table[i][j][k] = NULL;
         }
         free( eos_table[i][j] );
         eos_table[i][j] = NULL;
      }
      free( eos_table[i] );
      eos_table[i] = NULL;
   }
   free( eos_table );
   eos_table = NULL;


   return;
}
