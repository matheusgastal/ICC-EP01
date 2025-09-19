#include <math.h>
#include <stdlib.h>
#include <fenv.h>

#include "gs.h"

//Matheus Gastal Magalhães
//GRR20244620

int gauss_seidel_tridiag(Tridiag *sl, real_t *x, int maxit, real_t tol, int *iters, real_t *normL2){
    int n = sl->n;
    real_t *Di = sl->Di; // sub-diagonal
    real_t *D  = sl->D;  // diagonal
    real_t *Ds = sl->Ds; // super-diagonal
    real_t *B  = sl->B;  // RHS

    // inicializa x com zeros
    for (int i = 0; i < n; ++i) x[i] = 0.0;

    int k = 0;
    real_t erro = tol + 1.0;

    while (erro > tol && k < maxit) {
        k++;

        // primeira equação
        x[0] = (B[0] - Ds[0] * x[1]) / D[0];

        // equações do meio
        for (int i = 1; i < n-1; ++i)
            x[i] = (B[i] - Di[i] * x[i-1] - Ds[i] * x[i+1]) / D[i];

        // última equação
        x[n-1] = (B[n-1] - Di[n-1] * x[n-2]) / D[n-1];

        erro = 0.0;
        for (int i = 0; i < n; ++i) {
            real_t r = D[i]*x[i];
            if (i > 0) r += Di[i]*x[i-1];
            if (i < n-1) r += Ds[i]*x[i+1];
            r -= B[i];
            erro += r*r;
        }
        erro = sqrt(erro);
    }

    *iters = k;
    *normL2 = erro;

    return 0;
}

