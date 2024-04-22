#pragma once

typedef double (*NewuoaFunction)(long n, const double *x);
typedef double (*NewuoaClosureFunction)(void *data, long n, const double *x);
typedef double (*NewuoaClosureFunctionConst)(const void *data, long n, const double *x);

typedef struct {
    void *data;
    NewuoaClosureFunction function;
} NewuoaClosure;

typedef struct {
    const void *data;
    NewuoaClosureFunctionConst function;
} NewuoaClosureConst;

/* This subroutine seeks the least value of a function of many variables,
 * by a trust region method that forms quadratic models by interpolation.
 * There can be some freedom in the interpolation conditions, which is
 * taken up by minimizing the Frobenius norm of the change to the second
 * derivative of the quadratic model, beginning with a zero matrix. The
 * arguments of the subroutine are as follows. */

/* N must be set to the number of variables and must be at least two.
 * NPT is the number of interpolation conditions. Its value must be in the
 *   interval [N+2,(N+1)(N+2)/2].
 * Initial values of the variables must be set in X(1),X(2),...,X(N). They
 *   will be changed to the values that give the least calculated F.
 * RHOBEG and RHOEND must be set to the initial and final values of a trust
 *   region radius, so both must be positive with RHOEND<=RHOBEG. Typically
 *   RHOBEG should be about one tenth of the greatest expected change to a
 *   variable, and RHOEND should indicate the accuracy that is required in
 *   the final values of the variables.
 * MAXFUN must be set to an upper bound on the number of calls of CALFUN.
 * The array W will be used for working space. Its length must be at least
 * (NPT+13)*(NPT+N)+3*N*(N+3)/2. */

#ifdef __cplusplus
extern "C" {
#endif

double newuoa(
    NewuoaFunction function,
    long n,
    long npt,
    double *x,
    double rhobeg,
    double rhoend,
    long maxfun,
    double *w,
    int *counter
);

double newuoa_closure(
    NewuoaClosure *closure,
    long n,
    long npt,
    double *x,
    double rhobeg,
    double rhoend,
    long maxfun,
    double *w,
    int *counter
);

double newuoa_closure_const(
    const NewuoaClosureConst* closure,
    long n,
    long npt,
    double *x,
    double rhobeg,
    double rhoend,
    long maxfun,
    double *w,
    int *counter
);

#ifdef __cplusplus
}
#endif

#define NEWUOA_WORKING_SPACE_SIZE(n, npt) \
    (npt + 13)*(npt + n) + 3*n*(n + 3)/2

