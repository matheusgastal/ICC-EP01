#ifndef __GS_H__
#define __GS_H__

#include "edo.h"

int gauss_seidel_tridiag(Tridiag *sl, real_t *x, int maxit, real_t tol, int *iters, real_t *normL2);

#endif // __GS_H__
