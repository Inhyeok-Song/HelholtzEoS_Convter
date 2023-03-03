#ifndef _HEADER_H_
#define _HEADER_H_

#include <math.h>


#define MAX( a, b )          (  ( (a) > (b) ) ? (a) : (b)  )
#define MIN( a, b )          (  ( (a) < (b) ) ? (a) : (b)  )

#define FABS( a )        fabs( a )
#define SQRT( a )        sqrt( a )
#define LOG10( a )      log10( a )
#define LOG( a )          log( a )
#define POW( a, b )       pow( a, b )


#define NTABLES               4


// variable indices in the Nuclear EoS table and auxiliary table
#define NUC_VAR_IDX_PRES      0     // pressure
#define NUC_VAR_IDX_ENGY      1     // internal energy
#define NUC_VAR_IDX_ENTR      2     // entropy
#define NUC_VAR_IDX_CSQR      3     // sound speed squared


// variable indices for the Helmholtz EoS
#define PSI0                  0     // HelmholtzEoS
#define DPSI0                 1     // HelmholtzEoS
#define DDPSI0                2     // HelmholtzEoS
#define PSI1                  3     // HelmholtzEoS
#define DPSI1                 4     // HelmholtzEoS
#define DDPSI1                5     // HelmholtzEoS
#define PSI2                  6     // HelmholtzEoS
#define DPSI2                 7     // HelmholtzEoS
#define DDPSI2                8     // HelmholtzEoS
#define XPSI0                 9     // HelmholtzEoS
#define XDPSI0               10     // HelmholtzEoS
#define XPSI1                11     // HelmholtzEoS
#define XDPSI1               12     // HelmholtzEoS


const double Const_c             = 2.99792458e10;           // speed of light
const double Const_amu           = 1.660539040e-24;         // atomic mass unit
const double Const_kB            = 1.38064852e-16;          // Boltzmann constant in erg/K
const double Const_Planck        = 1.054571800e-27;         // reduced Planck constant in erg*s
const double Const_NA            = 6.022140857e23;          // Avogadro constant (in 1/mole)
const double Kelvin2MeV          = 8.6173303e-11;


extern int imax, jmax;
extern int nrho;
extern int ntemp;
extern int nye;

extern double *alltables;
extern double *logrho;
extern double *logtemp;
extern double *yes;
extern double *helm_dens;
extern double *helm_temp;
extern double *helmholtz_table;
extern double *helmholtz_dd;
extern double *helmholtz_dt;


// prototype
void set_bins( int rho_bs, int temp_bs, int ye_bs, double eos_rhomin, double eos_rhomax,
               double eos_tempmin, double eos_tempmax, double eos_yemin, double eos_yemax );

void Helmholtz_eos_C_ReadTable( char *helmeos_table_name );

void Helmholtz_eos( double *Out, const double *In, const int NTarget, const int *TargetIdx,
                    const int imax, const int jmax, const double *helm_dens, const double *helm_temp,
                    const double *helmholtz_table, const double *helmholtz_dd, const double *helmholtz_dt );

void write_eos_table( char *helmeos_table_name );


// Tolerance
const double Tolerance = 1e-10;

#endif
