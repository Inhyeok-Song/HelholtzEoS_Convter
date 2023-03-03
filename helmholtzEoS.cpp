#include <stdio.h>
#include "helmholtz.h"



extern int     imax;
extern int     jmax;

extern double *helmholtz_table;
extern double *helmholtz_dd;
extern double *helmholtz_dt;
extern double *helm_dens;
extern double *helm_temp;



static
double hermite_poly( const double z, const int TargetIdx );

static
double cubic_h5( const double *fi, const double w0t, const double w1t,
                 const double w2t, const double w0mt, const double w1mt,
                 const double w2mt, const double w0d, const double w1d,
                 const double w2d, const double w0md, const double w1md,
                 const double w2md );

static
double biquintic_h3( const double *fi, const double w0t, const double w1t,
                     const double w0mt, const double w1mt, const double w0d,
                     const double w1d, const double w0md, const double w1md );



void Helmholtz_eos( double *Out, const double *In, const int NTarget, const int *TargetIdx,
                    const int imax, const int jmax, const double *helm_dens, const double *helm_temp,
                    const double *helmholtz_table, const double *helmholtz_dd, const double *helmholtz_dt )
{

   double       dvardlt    = 0.0;
   int          max_iter   = 1;

   double       pres;
   double       ener;
   double       entr;
   double       entr_kb;
   double       sound;

   const double kB         = 1.38064852e-16;
   const double kB_eV      = 8.6173303e-5;
   const double c          = 2.99792458e10;
   const double MeV2Kelvin = 1e6 / kB_eV;
   const double amu        = 1.660539040e-24;
   const double Planck     = 1.054571800e-27;
   const double Na         = 6.022140857e23;

   const double sigma      = 5.670367e-5;
   const double rad_a      = 4*sigma / c;
   const double rad_ai3    = rad_a/3.0;
   const double qe         = 4.803204673e-10;
   const double esqu       = qe*qe;
   const double sioncon    = ( amu * kB ) / ( 2.0*M_PI*Planck*Planck );

   const double   a1       = -0.898004;
   const double   b1       =  0.96786;
   const double   c1       =  0.220703;
   const double   d1       = -0.86097;
   const double   e1       =  2.5269;
   const double   a2       =  0.29561;
   const double   b2       =  1.9885;
   const double   c2       =  0.288675;

   const double   tlo      = 3.0;
   const double   thi      = 13.0;
   const double   tstp     = (thi - tlo)/(jmax-1);
   const double   tstpi    = 1.0/tstp;
   const double   dlo      = -12.0;
   const double   dhi      = 15.0;
   const double   dstp     = (dhi - dlo)/(imax-1);
   const double   dstpi    = 1.0/dstp;


   const double Dens  = In[0];
   const double Temp  = In[1] * MeV2Kelvin;
   const double Ye    = In[2];



// enter the table with Ye*Dens
   double YeDens = Ye*Dens;
   const double Densi = 1.0/Dens;

   const double xmass[3] = { 0.75, 0.23, 0.02 };
   const double aion[3]  = { 1.0,  4.0,  12.0 };
   const double zion[3]  = { 1.0,  2.0,  6.0  };

   double sum_x_a  = 0.0;
   double sum_xz_a = 0.0;
   for (int i=0; i<3; i++) {
      sum_x_a  += xmass[i]/aion[i];
      sum_xz_a += xmass[i]*zion[i]/aion[i];
   }
   // const double abar  = 1.0/sum_x_a;
   // assume abar = 4.0
   const double abar  = 4.0;
   const double zbar  = abar * sum_xz_a;
   const double ytot1 = 1.0/abar;
   // const double Ye    = MAX( Tolerance, ytot1 * zbar );


   const double Tempi = 1.0/Temp;

   const double kt    = Temp*kB;
   const double ktinv = 1.0/kt;

   const double prad  = rad_ai3 * Temp * Temp * Temp * Temp;
   const double dpraddd = 0.0;
   const double dpraddt = 4.0 * prad * Tempi;

   const double erad    = 3.0 * prad * Densi;
   const double deraddt = 3.0 * dpraddt * Densi;

   const double srad    = ( prad*Densi + erad ) * Tempi;

// ion section
   const double xni     = Na * ytot1 * Dens;
   const double dxnidd  = Na * ytot1;
   const double dxnida  = -xni * ytot1;

   const double pion    = xni * kt;
   const double dpiondd = dxnidd * kt;
   const double dpiondt = xni * kB;
   const double dpionda = dxnida * kt;
   const double dpiondz = 0.0;

   const double eion    = 1.5 * pion * Densi;
   const double deiondt = 1.5 * dpiondt * Densi;

// sackur-tetrode equation for the ion entropy of
// a single ideal gas characterized by abar
   double x = abar*abar*SQRT(abar) * Densi/Na;
   double s = sioncon * Temp;
   double z = x * s * SQRT(s);
   double y = LOG(z);

   const double sion = ( pion*Densi + eion ) * Tempi + kB * Na * ytot1 * y;


   int jat = int( (LOG10(Temp) - tlo)*tstpi );
      jat = MAX( 0, MIN(jat, jmax-2) );
   int iat = int( (LOG10(YeDens) - dlo)*dstpi );
      iat = MAX( 0, MIN(iat, imax-2) );

   double fi[36];
   int arr_idx_9[9] = { 0, 2, 4, 1, 3, 5, 6, 7, 8 };

   for (int k=0; k<9; k++)
   for (int j=0; j<2; j++)
   for (int i=0; i<2; i++)
      fi[4*k + 2*j + i] = helmholtz_table[ ((iat+i)*jmax*21) + ((jat+j)*21) + arr_idx_9[k] ];


// various differences
   double xt  = MAX( (Temp - helm_temp[jat])  *helmholtz_dt[jat*5 + 2], 0.0 );
   double xd  = MAX( (YeDens - helm_dens[iat])*helmholtz_dd[iat*5 + 2],  0.0 );
   double mxt = 1.0 - xt;
   double mxd = 1.0 - xd;


// the six density and six temperature basis functions
   double si0t    =  hermite_poly( xt,  PSI0   );
   double si1t    =  hermite_poly( xt,  PSI1   ) * helmholtz_dt[jat*5 + 0];
   double si2t    =  hermite_poly( xt,  PSI2   ) * helmholtz_dt[jat*5 + 1];

   double si0mt   =  hermite_poly( mxt, PSI0   );
   double si1mt   = -hermite_poly( mxt, PSI1   ) * helmholtz_dt[jat*5 + 0];
   double si2mt   =  hermite_poly( mxt, PSI2   ) * helmholtz_dt[jat*5 + 1];

   double si0d    =  hermite_poly( xd,  PSI0   );
   double si1d    =  hermite_poly( xd,  PSI1   ) * helmholtz_dd[iat*5 + 0];
   double si2d    =  hermite_poly( xd,  PSI2   ) * helmholtz_dd[iat*5 + 1];

   double si0md   =  hermite_poly( mxd, PSI0   );
   double si1md   = -hermite_poly( mxd, PSI1   ) * helmholtz_dd[iat*5 + 0];
   double si2md   =  hermite_poly( mxd, PSI2   ) * helmholtz_dd[iat*5 + 1];

// derivatives of the weight functions
   double dsi0t   =  hermite_poly( xt,  DPSI0  ) * helmholtz_dt[jat*5 + 2];
   double dsi1t   =  hermite_poly( xt,  DPSI1  );
   double dsi2t   =  hermite_poly( xt,  DPSI2  ) * helmholtz_dt[jat*5 + 0];

   double dsi0mt  = -hermite_poly( mxt, DPSI0  ) * helmholtz_dt[jat*5 + 2];
   double dsi1mt  =  hermite_poly( mxt, DPSI1  );
   double dsi2mt  = -hermite_poly( mxt, DPSI2  ) * helmholtz_dt[jat*5 + 0];

   double dsi0d   =  hermite_poly( xd,  DPSI0  ) * helmholtz_dd[iat*5 + 2];
   double dsi1d   =  hermite_poly( xd,  DPSI1  );
   double dsi2d   =  hermite_poly( xd,  DPSI2  ) * helmholtz_dd[iat*5 + 0];

   double dsi0md  = -hermite_poly( mxd, DPSI0  ) * helmholtz_dd[iat*5 + 2];
   double dsi1md  =  hermite_poly( mxd, DPSI1  );
   double dsi2md  = -hermite_poly( mxd, DPSI2  ) * helmholtz_dd[iat*5 + 0];

// second derivatives of the weight functions
   double ddsi0t  =  hermite_poly( xt,  DDPSI0 ) * helmholtz_dt[jat*5 + 3];
   double ddsi1t  =  hermite_poly( xt,  DDPSI1 ) * helmholtz_dt[jat*5 + 2];
   double ddsi2t  =  hermite_poly( xt,  DDPSI2 );

   double ddsi0mt =  hermite_poly( mxt, DDPSI0 ) * helmholtz_dt[jat*5 + 3];
   double ddsi1mt = -hermite_poly( mxt, DDPSI1 ) * helmholtz_dt[jat*5 + 2];
   double ddsi2mt =  hermite_poly( mxt, DDPSI2 );


// the free energy
   double free  = cubic_h5( fi, si0t, si1t, si2t, si0mt, si1mt, si2mt,
                              si0d, si1d, si2d, si0md, si1md, si2md );

// derivative with respect to density
   double df_d  = cubic_h5( fi, si0t,  si1t,  si2t,  si0mt,  si1mt,  si2mt,
                              dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md );

// derivative with respect to temperature
   double df_t  = cubic_h5( fi, dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt,
                              si0d, si1d, si2d, si0md, si1md, si2md );


// derivative with respect to temperature**2
   double df_tt = cubic_h5( fi, ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt,
                              si0d, si1d, si2d, si0md, si1md, si2md );

// derivative with respect to temperature and density
   double df_dt = cubic_h5( fi, dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt,
                              dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md );

// now get the pressure derivative with density, chemical potential, and
// electron positron number densities
// get the interpolation weight functions
   si0t   =  hermite_poly( xt,  XPSI0  );
   si1t   =  hermite_poly( xt,  XPSI1  ) * helmholtz_dt[jat*5 + 0];;

   si0mt  =  hermite_poly( mxt, XPSI0  );
   si1mt  = -hermite_poly( mxt, XPSI1  ) * helmholtz_dt[jat*5 + 0];;

   si0d   =  hermite_poly( xd,  XPSI0  );
   si1d   =  hermite_poly( xd,  XPSI1  ) * helmholtz_dd[iat*5 + 0];

   si0md  =  hermite_poly( mxd, XPSI0  );
   si1md  = -hermite_poly( mxd, XPSI1  ) * helmholtz_dd[iat*5 + 0];


// derivatives of weight functions
   dsi0t  =  hermite_poly( xt,  XDPSI0 ) * helmholtz_dt[jat*5 + 2];
   dsi1t  =  hermite_poly( xt,  XDPSI1 );

   dsi0mt = -hermite_poly( mxt, XDPSI0 ) * helmholtz_dt[jat*5 + 2];
   dsi1mt =  hermite_poly( mxt, XDPSI1 );

   dsi0d  =  hermite_poly( xd,  XDPSI0 ) * helmholtz_dd[iat*5 + 2];
   dsi1d  =  hermite_poly( xd,  XDPSI1 );

   dsi0md = -hermite_poly( mxd, XDPSI0 ) * helmholtz_dd[iat*5 + 2];
   dsi1md =  hermite_poly( mxd, XDPSI1 );


// # look in the pressure derivative only once
   int arr_idx_4[4] = { 9, 11, 10, 12 };

   for (int k=0; k<4; k++)
   for (int j=0; j<2; j++)
   for (int i=0; i<2; i++)
      fi[4*k + 2*j + i] = helmholtz_table[ ((iat+i)*jmax*21) + ((jat+j)*21) + arr_idx_4[k] ];

// pressure derivative with density
   double dpepdd  = biquintic_h3( fi, si0t, si1t, si0mt, si1mt,
                                 si0d, si1d, si0md, si1md );
         dpepdd  = MAX( Ye * dpepdd, Tolerance );


// look in the electron chemical potential table only once
   arr_idx_4[0] = 13;  arr_idx_4[1] = 15;  arr_idx_4[2] = 14;  arr_idx_4[3] = 16;

   for (int k=0; k<4; k++)
   for (int j=0; j<2; j++)
   for (int i=0; i<2; i++)
      fi[4*k + 2*j + i] = helmholtz_table[ ((iat+i)*jmax*21) + ((jat+j)*21) + arr_idx_4[k] ];


// derivative with respect to density
   x = biquintic_h3( fi, si0t, si1t, si0mt, si1mt,
                     dsi0d, dsi1d, dsi0md, dsi1md );


// look in the number density table only once
   arr_idx_4[0] = 17;  arr_idx_4[1] = 19;  arr_idx_4[2] = 18;  arr_idx_4[3] = 20;

   for (int k=0; k<4; k++)
   for (int j=0; j<2; j++)
   for (int i=0; i<2; i++)
      fi[4*k + 2*j + i] = helmholtz_table[ ((iat+i)*jmax*21) + ((jat+j)*21) + arr_idx_4[k] ];

// the desired electron-positron thermodynamic quantities

// dpepdd at high temperatures and low densities is below the
// floating point limit of the subtraction of two large terms.
// since dpresdd doesn't enter the maxwell relations at all, use the
// bicubic interpolation done above instead of the formally correct expression
         x      = YeDens * YeDens;
   double pele   = x * df_d;
   double dpepdt = x * df_dt;
         s      = dpepdd/Ye - 2.0 * YeDens * df_d;

         x       = Ye * Ye;
   double sele    = -df_t * Ye;
   double dsepdt  = -df_tt * Ye;

   double eele    = Ye*free + Temp * sele;
   double deepdt  = Temp * dsepdt;

//  coulomb section:
//
// uniform background corrections only
// from yakovlev & shalybkov 1989
// lami is the average ion seperation
// plasg is the plasma coupling parameter
   const double third = 1.0/3.0;

                  z    = 4.0/3.0 * M_PI;
                  s    = z * xni;
   const double dsdd = z * dxnidd;
   const double dsda = z * dxnida;

   const double lami     = 1.0 / POW( s, third );
   const double inv_lami = 1.0/lami;
                  z        = -third * lami;
   const double lamidd   = z * dsdd/s;
   const double lamida   = z * dsda/s;

   const double plasg   = zbar*zbar*esqu*ktinv*inv_lami;
                  z       = -plasg * inv_lami;
   const double plasgdd = z * lamidd;
   const double plasgda = z * lamida;
   const double plasgdt = -plasg*ktinv * kB;
   const double plasgdz = 2.0 * plasg/zbar;

   double pcoul    = 0.0;  double dpcouldd = 0.0;  double dpcouldt = 0.0;  double dpcoulda = 0.0;
   double dpcouldz = 0.0;  double ecoul    = 0.0;  double decouldd = 0.0;  double decouldt = 0.0;
   double decoulda = 0.0;  double decouldz = 0.0;  double scoul    = 0.0;  double dscouldz = 0.0;

   if ( plasg >= 1.0 ) {
      x        = POW( plasg, 0.25 );
      y        = Na * ytot1 * kB;
      ecoul    = y * Temp * (a1*plasg + b1*x + c1/x + d1);
      pcoul    = third * Dens * ecoul;
      scoul    = -y * (3.0*b1*x - 5.0*c1/x + d1*(LOG(plasg) - 1.0) - e1);

      y        = Na*ytot1*kt*(a1 + 0.25/plasg*(b1*x - c1/x));
      decouldd = y * plasgdd;
      decouldt = y * plasgdt + ecoul/Temp;
      decoulda = y * plasgda - ecoul/abar;
      decouldz = y * plasgdz;

      y        = third * Dens;
      dpcouldd = third * ecoul + y*decouldd;
      dpcouldt = y * decouldt;
      dpcoulda = y * decoulda;
      dpcouldz = y * decouldz;

      y        = -Na*kB/(abar*plasg)*(0.75*b1*x+1.25*c1/x+d1);
      dscouldz = y * plasgdz;
   }

   else if ( plasg < 1.0 ) {
      x        = plasg*SQRT( plasg );
      y        = POW( plasg, b2 );
      z        = c2 * x - third * a2 * y;
      pcoul    = -pion * z;
      ecoul    = 3.0 * pcoul/Dens;
      scoul    = -Na/abar*kB*( c2*x -a2 * (b2-1.0) / b2*y );

      s        = 1.5*c2*x/plasg - third*a2*b2*y/plasg;
      dpcouldd = -dpiondd*z - pion*s*plasgdd;
      dpcouldt = -dpiondt*z - pion*s*plasgdt;
      dpcoulda = -dpionda*z - pion*s*plasgda;
      dpcouldz = -dpiondz*z - pion*s*plasgdz;

      s        = 3.0/Dens;
      decouldd = s * dpcouldd - ecoul/Dens;
      decouldt = s * dpcouldt;
      decoulda = s * dpcoulda;
      decouldz = s * dpcouldz;

      s        = -Na*kB/(abar*plasg)*(1.5*c2*x-a2*(b2-1.0)*y);
      dscouldz = s * plasgdz;
   }

// bomb proof
   x = prad + pion + pele + pcoul;
   y = erad + eion + eele + ecoul;
   z = srad + sion + sele + scoul;

   if ( x <= 0.0  ||  y <= 0.0 ) {
      pcoul    = 0.0;  dpcouldd = 0.0;  dpcouldt = 0.0;  dpcoulda = 0.0;
      dpcouldz = 0.0;  ecoul    = 0.0;  decouldd = 0.0;  decouldt = 0.0;
      decoulda = 0.0;  decouldz = 0.0;  scoul    = 0.0;  dscouldz = 0.0;
   }

// sum all the gas components
   const double pgas = pion + pele + pcoul;
   const double egas = eion + eele + ecoul;
   const double sgas = sion + sele + scoul;

   const double dpgasdd = dpiondd + dpepdd + dpcouldd;
   const double dpgasdt = dpiondt + dpepdt + dpcouldt;
   const double degasdt = deiondt + deepdt + decouldt;

// add in radiation to get the total
   pres    = prad + pgas;
   ener    = erad + egas;
   entr    = srad + sgas;
   entr_kb = entr / (kB*Na);

   const double dpresdd = dpraddd + dpgasdd;
   const double dpresdt = dpraddt + dpgasdt;
   const double denerdt = deraddt + degasdt;


// for the gas
// the temperature and density exponents (c&g 9.81 9.82)
// the specific heat at constant volume (c&g 9.92)
// the third adiabatic exponent (c&g 9.93)
// the first adiabatic exponent (c&g 9.97)
// the second adiabatic exponent (c&g 9.105)
// the specific heat at constant pressure (c&g 9.98)
// and relativistic formula for the sound speed (c&g 14.29)
   // double zz        = pgas*Densi;
   // double zzi       = Dens/pgas;
   // double chit_gas  = Temp/pgas * dpgasdt;
   // double chid_gas  = dpgasdd*zzi;
   // double cv_gas    = degasdt;
   //              x         = zz * chit_gas/(Temp * cv_gas);
   // double gam1_gas  = chit_gas*x + chid_gas;
   //              z         = 1.0 + ( egas + c*c ) * zzi;
   // double sound_gas = c * SQRT(gam1_gas/z);

// // for the totals
   const double zz    = pres*Densi;
   const double zzi   = Dens/pres;
   const double chit  = Temp/pres * dpresdt;
   const double chid  = dpresdd*zzi;
         double cv    = denerdt;
                x     = zz * chit/(Temp * cv);
   const double gam1  = chit*x + chid;
                z     = 1.0 + (ener + c*c)*zzi;
                sound = c * SQRT(gam1/z);



   for (int i=0; i<NTarget; i++)
   {
      switch ( TargetIdx[i] )
      {
         case NUC_VAR_IDX_PRES :
            Out[i] = LOG10( pres );
            break;

         case NUC_VAR_IDX_ENGY :
            Out[i] = LOG10( ener );
            break;

         case NUC_VAR_IDX_ENTR :
            Out[i] = entr_kb;
            break;

         case NUC_VAR_IDX_CSQR :
            Out[i] = sound * sound;
            break;

//       unsupported output variables
         default :
         {
            fprintf( stderr, "unspported output variables\n" );
            abort();
         }

      }
   }


} // FUNCTION : Helmholtz_eos



double hermite_poly( const double z, const int TargetIdx )
{

   double res = 0.0;
   switch ( TargetIdx )
   {
   case PSI0   :
      res = z*z*z * ( z * (-6.0 * z + 15.0) - 10.0 ) + 1.0;
      break;

   case DPSI0  :
      res = z * z * ( z * (-30.0 * z + 60.0 ) - 30.0 );
      break;

   case DDPSI0 :
      res = z* ( z * ( -120.0*z + 180.0 ) -60.0 );
      break;

   case PSI1:
      res = z * ( z*z * ( z * ( -3.0*z + 8.0 ) - 6.0 ) + 1.0 );
      break;

   case DPSI1  :
      res = z*z * ( z * ( -15.0*z + 32.0 ) - 18.0 ) + 1.0;
      break;

   case DDPSI1 :
      res = z * ( z * (-60.0*z + 96.0 ) - 36.0 );
      break;

   case PSI2   :
      res = 0.5*z*z * ( z * ( z * ( -z + 3.0 ) - 3.0 ) + 1.0 );
      break;

   case DPSI2  :
      res = 0.5*z * ( z * ( z * ( -5.0*z + 12.0 ) - 9.0 ) + 2.0 );
      break;

   case DDPSI2 :
      res = 0.5 * ( z * ( z * ( -20.0*z + 36.0 ) - 18.0 ) + 2.0 );
      break;

   case XPSI0  :
      res = z*z * ( 2.0*z - 3.0 ) + 1.0;
      break;

   case XDPSI0 :
      res = z * ( 6.0*z - 6.0 );
      break;

   case XPSI1  :
      res = z * ( z * ( z - 2.0 ) + 1.0 );
      break;

   case XDPSI1 :
      res = z * ( 3.0*z - 4.0 ) + 1.0;
      break;
   }

   return res;

}

double cubic_h5( const double *fi, const double w0t, const double w1t,
                 const double w2t, const double w0mt, const double w1mt,
                 const double w2mt, const double w0d, const double w1d,
                 const double w2d, const double w0md, const double w1md,
                 const double w2md )
{

   double res;

   res = fi[0 ] * w0d * w0t   + fi[1 ] * w0md * w0t
       + fi[2 ] * w0d * w0mt  + fi[3 ] * w0md * w0mt
       + fi[4 ] * w0d * w1t   + fi[5 ] * w0md * w1t
       + fi[6 ] * w0d * w1mt  + fi[7 ] * w0md * w1mt
       + fi[8 ] * w0d * w2t   + fi[9 ] * w0md * w2t
       + fi[10] * w0d * w2mt  + fi[11] * w0md * w2mt
       + fi[12] * w1d * w0t   + fi[13] * w1md * w0t
       + fi[14] * w1d * w0mt  + fi[15] * w1md * w0mt
       + fi[16] * w2d * w0t   + fi[17] * w2md * w0t
       + fi[18] * w2d * w0mt  + fi[19] * w2md * w0mt
       + fi[20] * w1d * w1t   + fi[21] * w1md * w1t
       + fi[22] * w1d * w1mt  + fi[23] * w1md * w1mt
       + fi[24] * w2d * w1t   + fi[25] * w2md * w1t
       + fi[26] * w2d * w1mt  + fi[27] * w2md * w1mt
       + fi[28] * w1d * w2t   + fi[29] * w1md * w2t
       + fi[30] * w1d * w2mt  + fi[31] * w1md * w2mt
       + fi[32] * w2d * w2t   + fi[33] * w2md * w2t
       + fi[34] * w2d * w2mt  + fi[35] * w2md * w2mt;

   return res;
}

double biquintic_h3( const double *fi, const double w0t, const double w1t,
                     const double w0mt, const double w1mt, const double w0d,
                     const double w1d, const double w0md, const double w1md )
{
   double res;

   res = fi[0 ]*w0d*w0t  + fi[1 ]*w0md*w0t
       + fi[2 ]*w0d*w0mt + fi[3 ]*w0md*w0mt
       + fi[4 ]*w0d*w1t  + fi[5 ]*w0md*w1t
       + fi[6 ]*w0d*w1mt + fi[7 ]*w0md*w1mt
       + fi[8 ]*w1d*w0t  + fi[9 ]*w1md*w0t
       + fi[10]*w1d*w0mt + fi[11]*w1md*w0mt
       + fi[12]*w1d*w1t  + fi[13]*w1md*w1t
       + fi[14]*w1d*w1mt + fi[15]*w1md*w1mt;

    return res;
}
