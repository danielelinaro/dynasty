#include <stdio.h>
#include <math.h>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */

#define Ith(v,i)    NV_Ith_S(v,i)         /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) SM_ELEMENT_D(A,i,j) /* IJth numbers rows,cols 1..NEQ */

/* Problem Constants */

#define NEQ   3                /* number of equations  */
#define NPARS 7                /* number of equations  */
#define Y0    RCONST(0.8)      /* initial y components */
#define Y1    RCONST(0.1)
#define Y2    RCONST(0.1)
#define RTOL  RCONST(1.0e-10)   /* scalar relative tolerance            */
#define ATOL1 RCONST(1.0e-6)   /* vector absolute tolerance components */
#define ATOL2 RCONST(1.0e-6)
#define ATOL3 RCONST(1.0e-6)
#define T0    RCONST(0.0)      /* initial time           */
#define T1    RCONST(0.4)      /* first output time      */
#define TMULT RCONST(10.0)     /* output time factor     */
#define NOUT  12               /* number of output times */

#define ZERO  RCONST(0.0)

#define EPS   RCONST(1.0E-3)

#define R     RCONST(0.1)
#define E     RCONST(2.7)
#define B     RCONST(0.17)
#define D     RCONST(0.42)
#define G     RCONST(0.09)
#define H     RCONST(0.1)
#define Q     RCONST(0.4)

/* Functions Called by the Solver */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int g(realtype t, N_Vector y, realtype *gout, void *user_data);

static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, 
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */

static int check_retval(void *returnvalue, const char *funcname, int opt);


/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main()
{
        realtype reltol, t, tout, tstep;
        N_Vector y, pars, abstol;
        SUNMatrix A;
        SUNLinearSolver LS;
        void *cvode_mem;
        int retval, retvalr;
        int rootsfound[2];
        int nev;
        
        y = pars = abstol = NULL;
        A = NULL;
        LS = NULL;
        cvode_mem = NULL;
        
        /* Create serial vector of length NEQ for I.C. and abstol */
        y = N_VNew_Serial(NEQ);
        if (check_retval((void *)y, "N_VNew_Serial", 0)) return(1);
        pars = N_VNew_Serial(NPARS);
        if (check_retval((void *)y, "N_VNew_Serial", 0)) return(1);
        abstol = N_VNew_Serial(NEQ); 
        if (check_retval((void *)abstol, "N_VNew_Serial", 0)) return(1);
        
        /* Initialize y */
        Ith(y,0) = Y0;
        Ith(y,1) = Y1;
        Ith(y,2) = Y2;

        /* Initialize pars */
        Ith(pars,0) = R;
        Ith(pars,1) = E;
        Ith(pars,2) = B;
        Ith(pars,3) = D;
        Ith(pars,4) = G;
        Ith(pars,5) = H;
        Ith(pars,6) = Q;
        
        /* Set the scalar relative tolerance */
        reltol = RTOL;
        /* Set the vector absolute tolerance */
        Ith(abstol,0) = ATOL1;
        Ith(abstol,1) = ATOL2;
        Ith(abstol,2) = ATOL3;
        
        /* Call CVodeCreate to create the solver memory and specify the 
         * Backward Differentiation Formula */
        cvode_mem = CVodeCreate(CV_BDF);
        if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);
        
        /* Call CVodeInit to initialize the integrator memory and specify the
         * user's right hand side function in y'=f(t,y), the inital time T0, and
         * the initial dependent variable vector y. */
        retval = CVodeInit(cvode_mem, f, T0, y);
        if (check_retval(&retval, "CVodeInit", 1)) return(1);
        
        /* Call CVodeSVtolerances to specify the scalar relative tolerance
         * and vector absolute tolerances */
        retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
        if (check_retval(&retval, "CVodeSVtolerances", 1)) return(1);
        
        /* Call CVodeRootInit to specify the root function g with 1 component */
        retval = CVodeRootInit(cvode_mem, 1, g);
        if (check_retval(&retval, "CVodeRootInit", 1)) return(1);
        
        /* Create dense SUNMatrix for use in linear solves */
        A = SUNDenseMatrix(NEQ, NEQ);
        if(check_retval((void *)A, "SUNDenseMatrix", 0)) return(1);
        
        /* Create dense SUNLinearSolver object for use by CVode */
        LS = SUNLinSol_Dense(y, A);
        if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) return(1);
        
        /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
        retval = CVodeSetLinearSolver(cvode_mem, LS, A);
        if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);
        
        /* Set the user-supplied Jacobian routine Jac */
        retval = CVodeSetJacFn(cvode_mem, Jac);
        if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);
        
        retval = CVodeSetUserData(cvode_mem, (void *) pars);
        if(check_retval(&retval, "CVodeSetUserData", 1)) return(1);

        /* In loop, call CVode, print results, and test for error.
           Break out of loop when NOUT preset output times have been reached.  */
        
        t = 0.0;
        tstep = 0.05;
        //tout = 100000;
        tout = tstep;
        nev = 0;
        while(nev < 500 && t < tout) {
                retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
                if (retval == CV_ROOT_RETURN) {
                        retvalr = CVodeGetRootInfo(cvode_mem, rootsfound);
                        if (check_retval(&retvalr, "CVodeGetRootInfo", 1)) return(1);
                        printf("%g %g %g %g 1\n", t, Ith(y,0), Ith(y,1), Ith(y,2));
                        nev++;
                }
                else {
                        printf("%g %g %g %g 0\n", t, Ith(y,0), Ith(y,1), Ith(y,2));
                }
                if (check_retval(&retval, "CVode", 1)) {
                        break;
                }
                tout += tstep;
        }
        
        /* Free y and abstol vectors */
        N_VDestroy(y);
        N_VDestroy(pars);
        N_VDestroy(abstol);
        
        /* Free integrator memory */
        CVodeFree(&cvode_mem);
        
        /* Free the linear solver memory */
        SUNLinSolFree(LS);
        
        /* Free the matrix memory */
        SUNMatDestroy(A);
        
        return(retval);
}


/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/*
 * f routine. Compute function f(t,y). 
 */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
        realtype y1, y2, y3;
        realtype r, e, b, d, g, h, q;
        N_Vector pars = (N_Vector) user_data;
        
        r = Ith(pars,0);
        e = Ith(pars,1);
        b = Ith(pars,2);
        d = Ith(pars,3);
        g = Ith(pars,4);
        h = Ith(pars,5);
        q = Ith(pars,6);
        
        y1 = Ith(y,0);
        y2 = Ith(y,1);
        y3 = Ith(y,2);
        
        if (y1 < 0.) {
                if (y1 > -EPS)
                        y1 = 0.;
                else {
                        fprintf(stderr, "y1 = %g\n", y1);
                        return CV_ERR_FAILURE;
                }
        }
        if (y2 < 0.) {
                if (y2 > -EPS)
                        y2 = 0.;
                else {
                        fprintf(stderr, "y2 = %g\n", y2);
                        return CV_ERR_FAILURE;
                }
        }
        if (y3 < 0.) {
                if (y3 > -EPS)
                        y3 = 0.;
                else {
                        fprintf(stderr, "y3 = %g\n", y3);
                        return CV_ERR_FAILURE;
                }
        }
                
        //if(y1 < 0.0 || y2 < 0.0 || y3 < 0.0) {
        //        // small negative components are trimmed to zero
        //        if(y1 > -EPS)
        //                y1 = 0.0;
        //        else {
        //                fprintf(stderr, "y1 = %g\n", y1);
        //                return CV_ERR_FAILURE;
        //        }
        //        if(y2 > -EPS)
        //                y2 = 0.0;
        //        else {
        //                fprintf(stderr, "y2 = %g\n", y2);
        //                return CV_ERR_FAILURE;
        //        }
        //        if(y3 > -EPS)
        //                y3 = 0.0;
        //        else {
        //                fprintf(stderr, "y3 = %g\n", y3);
        //                return CV_ERR_FAILURE;
        //        }
        //}
        
        Ith(ydot,0) = y1*(1-y1-y2/(b+y1)-h*y3);
        Ith(ydot,1) = q*y2*(e*y1/(b+y1)-1-y3/(d+y2));
        Ith(ydot,2) = r*(y1*y2/(b+y1)-g*y3);

        return(0);
}

/*
 * g routine. Compute functions g_i(t,y) for i = 0,1. 
 */

static int g(realtype t, N_Vector y, realtype *gout, void *user_data)
{
        realtype y1, y2, y3;
        realtype b, h;
        N_Vector pars = (N_Vector) user_data;
        
        b = Ith(pars,2);
        h = Ith(pars,5);
        
        y1 = Ith(y,0);
        y2 = Ith(y,1);
        y3 = Ith(y,2);
        
        gout[0] = y1*(1-y1-y2/(b+y1)-h*y3);
        return(0);
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */

static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, 
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
        realtype y1, y2, y3;
        realtype r, e, b, d, g, h, q;
        realtype deny1, deny2, deny1square, deny2square;
        N_Vector pars = (N_Vector) user_data;
        
        r = Ith(pars,0);
        e = Ith(pars,1);
        b = Ith(pars,2);
        d = Ith(pars,3);
        g = Ith(pars,4);
        h = Ith(pars,5);
        q = Ith(pars,6);
        
        y1 = Ith(y,0);
        y2 = Ith(y,1);
        y3 = Ith(y,2);
        
        deny1 = 1/(b+y1);
        deny2 = 1/(d+y2);
        deny1square = deny1*deny1;
        deny2square = deny2*deny2;
        
        IJth (J, 0, 0) = 1 - 2*y1 - b*y2*deny1square - h*y3;
        IJth (J, 0, 1) = - y1*deny1;
        IJth (J, 0, 2) = - h*y1;
        IJth (J, 1, 0) = q*e*b*y2*deny1square;
        IJth (J, 1, 1) = q*e*y1*deny1 - q*d*y3*deny2square - q;
        IJth (J, 1, 2) = -q*y2*deny2;
        IJth (J, 2, 0) = r*b*y2*deny1square;
        IJth (J, 2, 1) = r*y1*deny1;
        IJth (J, 2, 2) = -r*g;

        return(0);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_retval(void *returnvalue, const char *funcname, int opt)
{
        int *retval;
        
        /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
        if (opt == 0 && returnvalue == NULL) {
                fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
                return(1);
        }
        
        /* Check if retval < 0 */
        else if (opt == 1) {
                retval = (int *) returnvalue;
                if (*retval < 0) {
                fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n", funcname, *retval);
                return(1);
                }
        }
        
        /* Check if function returned NULL pointer - no memory allocated */
        else if (opt == 2 && returnvalue == NULL) {
                fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
                return(1);
        }
        
        return(0);
}

