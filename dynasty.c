#include <stdlib.h>
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
#define NPARS 7                /* number of parameters */
#define Y0    RCONST(0.8)      /* initial y components */
#define Y1    RCONST(0.1)
#define Y2    RCONST(0.1)
#define RTOL  RCONST(1.0e-12)  /* scalar relative tolerance            */
#define ATOL1 RCONST(1.0e-10)  /* vector absolute tolerance components */
#define ATOL2 RCONST(1.0e-10)
#define ATOL3 RCONST(1.0e-10)
#ifndef PYTHON
#define TTRAN RCONST(1.0e4)
#define TEND  RCONST(1.0e5)
#endif

#define EPS   RCONST(1.0E-3)

#ifndef PYTHON
#define R     RCONST(0.07)
#define E     RCONST(2.0)
#define B     RCONST(0.17)
#define D     RCONST(0.42)
#define G     RCONST(0.09)
#define H     RCONST(0.1)
#define Q     RCONST(0.4)
#endif

#define MAX_STEPS 50000
#ifndef PYTHON
#define MAX_NEV   64
#endif

//#define DEBUG

/* Functions Called by the Solver */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int g(realtype t, N_Vector y, realtype *gout, void *user_data);

static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, 
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */

static int check_retval(void *returnvalue, const char *funcname, int opt);

#ifdef TRIM_ZERO
static int trim_zero(realtype *y1, realtype *y2, realtype *y3);
#endif


/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

#ifdef PYTHON

#include <Python.h>

int integrate(double R,
	      double E,
	      double B,
	      double D,
	      double G,
              double H,
              double Q,
	      double TTRAN,
	      double TEND,
	      size_t MAX_NEV)

#else

int main()

#endif
{
        realtype reltol, t;
        N_Vector y, pars, abstol;
        SUNMatrix A;
        SUNLinearSolver LS;
        void *cvode_mem;
        int retval, retvalr;
        int rootsfound[2], rootdir;
        size_t nev = 0;

#ifndef PYTHON
	size_t i;
	realtype buf[NEQ+1], pars_buf[NPARS];
        FILE *fid;
#endif
        
        y = pars = abstol = NULL;
        A = NULL;
        LS = NULL;
        cvode_mem = NULL;
        
         /* Create serial vector of length NEQ for I.C. and abstol */
        y = N_VNew_Serial(NEQ);
        if (check_retval((void *) y, "N_VNew_Serial", 0)) {
                return(1);
        }

        pars = N_VNew_Serial(NPARS);
        if (check_retval((void *) pars, "N_VNew_Serial", 0)) {
                N_VDestroy(y);
                return(1);
        }

        abstol = N_VNew_Serial(NEQ); 
        if (check_retval((void *) abstol, "N_VNew_Serial", 0)) {
                N_VDestroy(pars);
                N_VDestroy(y);
                return(1);
        }
        
#ifndef PYTHON
        fid = fopen("dynasty.dat", "wb");
        if (fid == NULL) {
                N_VDestroy(abstol);
                N_VDestroy(pars);
                N_VDestroy(y);
                return(1);
        }
#endif

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
        if ( check_retval((void *)cvode_mem, "CVodeCreate", 0) ) {
                retval = 1;
                goto free_vec;
        }
        
        /* Call CVodeInit to initialize the integrator memory and specify the
         * user's right hand side function in y'=f(t,y), the inital time, and
         * the initial dependent variable vector y. */
        retval = CVodeInit(cvode_mem, f, 0.0, y);
        if ( check_retval(&retval, "CVodeInit", 1) )
                goto free_mem;
        
        /* Call CVodeSVtolerances to specify the scalar relative tolerance
         * and vector absolute tolerances */
        retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
        if ( check_retval(&retval, "CVodeSVtolerances", 1) )
                goto free_mem;
        
        /* Create dense SUNMatrix for use in linear solvers */
        A = SUNDenseMatrix(NEQ, NEQ);
        if ( check_retval((void *)A, "SUNDenseMatrix", 0) ) {
                retval = 1;
                goto free_matrix;
        }
        
        /* Create dense SUNLinearSolver object for use by CVode */
        LS = SUNLinSol_Dense(y, A);
        if ( check_retval((void *)LS, "SUNLinSol_Dense", 0) ) {
                retval = 1;
                goto free_solver;
        }
        
        /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
        retval = CVodeSetLinearSolver(cvode_mem, LS, A);
        if( check_retval( &retval, "CVodeSetLinearSolver", 1) )
                goto free_solver;

        /* Set the user-supplied Jacobian routine Jac */
        retval = CVodeSetJacFn(cvode_mem, Jac);
        if( check_retval(&retval, "CVodeSetJacFn", 1) )
                goto free_solver;
        
        /* Set the parameters of the system */
        retval = CVodeSetUserData(cvode_mem, (void *) pars);
        if( check_retval(&retval, "CVodeSetUserData", 1) )
                goto free_solver;

        retval = CVodeSetMaxNumSteps(cvode_mem, MAX_STEPS);
        if( check_retval(&retval, "CVodeSetMaxNumSteps", 1) )
                goto free_solver;

        t = 0.0;
        while(t < TTRAN) {
                retval = CVode(cvode_mem, TTRAN, y, &t, CV_NORMAL);
                if (check_retval(&retval, "CVode", 1))
                        goto free_solver;
        }

        /* Call CVodeRootInit to specify the root function g with 1 component */
        retval = CVodeRootInit(cvode_mem, 1, g);
        if ( check_retval(&retval, "CVodeRootInit", 1) )
                goto free_solver;

        /* Call CVodeSetRootDirection to specify that only negative crossings of *
         * the Poincare' section should be reported                              */
        rootdir = -1;
        retval = CVodeSetRootDirection(cvode_mem, &rootdir);
        if ( check_retval(&retval, "CVodeRootInit", 1) )
                goto free_solver;

#ifndef PYTHON
        fprintf(fid, "NEQ: %d\n", NEQ);
        fprintf(fid, "NPARS: %d\n", NPARS);
        fprintf(fid, "NEV: %d\n", MAX_NEV);
        fprintf(fid, "NINT: %d\n", 1);
        fprintf(fid, "sizeof(size_t): %lu\n", sizeof(size_t));
        fprintf(fid, "sizeof(realtype): %lu\n", sizeof(realtype));
        fprintf(fid, "Binary:\n");

        fwrite(&nev, 1, sizeof(size_t), fid);

        for (i=0; i<NPARS; i++)
                pars_buf[i] = Ith(pars, i);
        fwrite(pars_buf, NPARS, sizeof(realtype), fid);
#endif

        while (nev < MAX_NEV && t < TEND) {
                retval = CVode(cvode_mem, TEND, y, &t, CV_NORMAL);
                if (retval == CV_ROOT_RETURN) {
                        retvalr = CVodeGetRootInfo(cvode_mem, rootsfound);
                        if (check_retval(&retvalr, "CVodeGetRootInfo", 1))
                                break;
                        ++nev;
#ifndef PYTHON
                        buf[0] = t;
                        buf[1] = Ith(y,0);
                        buf[2] = Ith(y,1);
                        buf[3] = Ith(y,2);
                        fwrite(&nev, 1, sizeof(size_t), fid);
                        fwrite(buf, NEQ+1, sizeof(realtype), fid);
#endif
#ifdef DEBUG
                        printf("[%03zu/%03zu]  %8.2f  %7.5f  %7.5f  %7.5f\n", nev, MAX_NEV, t, Ith(y,0), Ith(y,1), Ith(y,2));
#endif
                }
                if (check_retval(&retval, "CVode", 1))
                        break;
        }
        
free_solver:
        /* Free the linear solver memory */
        SUNLinSolFree(LS);
        
free_matrix:
        /* Free the matrix memory */
        SUNMatDestroy(A);
        
free_mem:
        /* Free integrator memory */
        CVodeFree(&cvode_mem);

free_vec:
#ifndef PYTHON
        /* Close the data file */
        fclose(fid);
#endif

        /* Free y, pars and abstol vectors */
        N_VDestroy(abstol);
        N_VDestroy(pars);
        N_VDestroy(y);

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

        r = exp(Ith(pars,0));
        e = Ith(pars,1);
        b = Ith(pars,2);
        d = Ith(pars,3);
        g = Ith(pars,4);
        h = Ith(pars,5);
        q = Ith(pars,6);
        
        y1 = Ith(y,0);
        y2 = Ith(y,1);
        y3 = Ith(y,2);
        
#ifdef TRIM_ZERO
        if ( trim_zero(&y1, &y2, &y3) )
                return CV_ERR_FAILURE;
#endif

        Ith(ydot,0) = y1 * ( 1 - y1 - y2 / ( b + y1 ) - h * y3 );
        Ith(ydot,1) = q * y2 * ( e * y1 / ( b + y1 ) - 1 - y3 / ( d + y2 ) );
        Ith(ydot,2) = r * ( y1 * y2 / ( b + y1 ) - g * y3 );

        return(0);
}

/*
 * g routine. Compute functions g_i(t,y) for i = 0,1. 
 */

static int g(realtype t, N_Vector y, realtype *gout, void *user_data)
{
        realtype y1, y2, y3;
        realtype q, e, b, d;
        N_Vector pars = (N_Vector) user_data;
        
        e = Ith(pars,1);
        b = Ith(pars,2);
        d = Ith(pars,3);
        q = Ith(pars,6);

        y1 = Ith(y,0);
        y2 = Ith(y,1);
        y3 = Ith(y,2);
        
#ifdef TRIM_ZERO
        if ( trim_zero(&y1, &y2, &y3) )
                return CV_ERR_FAILURE;
#endif

        gout[0] = q * y2 * ( e * y1 / ( b + y1 ) - 1 - y3 / ( d + y2 ) );

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
        
        r = exp(Ith(pars,0));
        e = Ith(pars,1);
        b = Ith(pars,2);
        d = Ith(pars,3);
        g = Ith(pars,4);
        h = Ith(pars,5);
        q = Ith(pars,6);
        
        y1 = Ith(y,0);
        y2 = Ith(y,1);
        y3 = Ith(y,2);
        
#ifdef TRIM_ZERO
        if ( trim_zero(&y1, &y2, &y3) )
                return CV_ERR_FAILURE;
#endif

        deny1 = 1 / ( b + y1 );
        deny2 = 1 / ( d + y2 );
        deny1square = deny1 * deny1;
        deny2square = deny2 * deny2;

        IJth (J, 0, 0) = 1 - 2 * y1 - b * y2 * deny1square - h * y3;
        IJth (J, 0, 1) = - y1 * deny1;
        IJth (J, 0, 2) = - h * y1;
        IJth (J, 1, 0) = q * e * b * y2 * deny1square;
        IJth (J, 1, 1) = q * e * y1 * deny1 - q * d * y3 * deny2square - q;
        IJth (J, 1, 2) = - q * y2 * deny2;
        IJth (J, 2, 0) = r * b * y2 * deny1square;
        IJth (J, 2, 1) = r * y1 * deny1;
        IJth (J, 2, 2) = -r * g;

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

int trim_zero(realtype *y1, realtype *y2, realtype *y3) {
        /* small negative components are trimmed to zero */
        if(*y1 < 0.0 || *y2 < 0.0 || *y3 < 0.0) {
                if(*y1 > -EPS) {
                        *y1 = 0.0;
                }
                else {
                        fprintf(stderr, "y1 = %g\n", *y1);
                        return 1;
                }
                if(*y2 > -EPS) {
                        *y2 = 0.0;
                }
                else {
                        fprintf(stderr, "y2 = %g\n", *y2);
                        return 1;
                }
                if(*y3 > -EPS) {
                        *y3 = 0.0;
                }
                else {
                        fprintf(stderr, "y3 = %g\n", *y3);
                        return 1;
                }
        }
        return 0;
}

