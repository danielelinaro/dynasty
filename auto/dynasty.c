#include "auto_f2c.h"

#define B 0.17
#define D 0.42
#define G 0.09
#define H 0.1
#define Q 0.4

int func (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *f, doublereal *dfdu, doublereal *dfdp)
{
        /* System generated locals */
        integer dfdu_dim1 = ndim;
        integer dfdp_dim1 = ndim;

        double x, y, z;
        double r, e, b, d, g, h, q;
        double denx, deny, denxsquare, denysquare;
        
        x = u[0];
        y = u[1];
        z = u[2];

        r = par[0];
        e = par[1];
        b = par[2];
        d = par[3];
        g = par[4];
        h = par[5];
        q = par[6];
        
        f[0] = x*(1-x-y/(b+x)-h*z);
        f[1] = q*y*(e*x/(b+x)-1-z/(d+y));
        f[2] = r*(x*y/(b+x)-g*z);

        if (ijac == 0) {
                return 0;
        }

        denx = 1/(b+x);
        deny = 1/(d+y);
        denxsquare = denx*denx;
        denysquare = deny*deny;
        
        ARRAY2D(dfdu, 0, 0) = 1 - 2*x - b*y*denxsquare - h*z;
        ARRAY2D(dfdu, 0, 1) = - x*denx;
        ARRAY2D(dfdu, 0, 2) = - h*x;
        ARRAY2D(dfdu, 1, 0) = q*e*b*y*denxsquare;
        ARRAY2D(dfdu, 1, 1) = q*e*x*denx - q*d*z*denysquare - q;
        ARRAY2D(dfdu, 1, 2) = -q*y*deny;
        ARRAY2D(dfdu, 2, 0) = r*b*y*denxsquare;
        ARRAY2D(dfdu, 2, 1) = r*x*denx;
        ARRAY2D(dfdu, 2, 2) = -r*g;

        if (ijac == 1) {
                return 0;
        }

        return 0;
}

int stpnt (integer ndim, doublereal t,
           doublereal *u, doublereal *par)
{
        par[0] = 0.1;
        par[1] = 1.5;
        par[2] = B;
        par[3] = D;
        par[4] = G;
        par[5] = H;
        par[6] = Q;

        u[0] = 0.97691713;
        u[1] = 0.01269461;
        u[2] = 0.12014408;

        return 0;
}

int pvls (integer ndim, const doublereal *u,
          doublereal *par)
{

        return 0;
}

int bcnd (integer ndim, const doublereal *par, const integer *icp,
          integer nbc, const doublereal *u0, const doublereal *u1, integer ijac,
          doublereal *fb, doublereal *dbc)
{
        return 0;
}

int icnd (integer ndim, const doublereal *par, const integer *icp,
          integer nint, const doublereal *u, const doublereal *uold,
          const doublereal *udot, const doublereal *upold, integer ijac,
          doublereal *fi, doublereal *dint)
{
        return 0;
}

int fopt (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *fs, doublereal *dfdu, doublereal *dfdp)
{
        return 0;
}

