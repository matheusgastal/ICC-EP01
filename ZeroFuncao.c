#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <fenv.h>
#include "utils.h"
#include "ZeroFuncao.h"
#include "DoubleType.h"

// ---------------- DIFERENÇA EM ULPs ----------------
long ulps_diff(double a, double b) {
    Double_t da, db;
    da.f = a;
    db.f = b;

    if (da.i >> 63) da.i = 0x8000000000000000ULL - da.i;
    if (db.i >> 63) db.i = 0x8000000000000000ULL - db.i;

    return llabs((long long)da.i - (long long)db.i);
}

// ---------------- NEWTON-RAPHSON ----------------
real_t newtonRaphson(Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz, int tipo_calc) {
    real_t px, dpx, dpx_antigo = 1.0;
    real_t x_novo = x0;
    real_t erro = 0.0;

    do {
        (*it)++;
        real_t x_ant = x_novo;

        if (tipo_calc == 1)
            calcPolinomio_lento(p, x_ant, &px, &dpx);
        else
            calcPolinomio_rapido(p, x_ant, &px, &dpx);

        if (fabs(dpx) > ZERO) {
            dpx_antigo = dpx;
            x_novo = x_ant - (px / dpx);
        } else {
            x_novo = x_ant - (px / dpx_antigo);
        }

        // erro de acordo com o critério
        if (criterioParada == 0)
            erro = (fabs(x_novo) > ZERO) ? fabs(x_novo - x_ant)/fabs(x_novo) : fabs(x_novo - x_ant);
        else if (criterioParada == 1)
            erro = fabs(px);
        else // ULPS
            erro = (real_t) ulps_diff(x_novo, x_ant);

        // condição de parada
        if ((criterioParada == 0 && erro < EPS) ||
            (criterioParada == 1 && fabs(px) <= DBL_EPSILON) ||
            (criterioParada == 2 && ulps_diff(x_novo, x_ant) <= ULPS))
            break;

    } while (*it < MAXIT);

    *raiz = x_novo;
    return erro;
}

// ---------------- BISSECCAO ----------------
real_t bisseccao(Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz, int tipo_calc) {
    real_t px_a, px_b, x_ant, x_novo = a;
    real_t dpx_a = 0.0, dpx_b = 0.0; 
    real_t erro = 0.0;

    if (tipo_calc == 1) {
        calcPolinomio_lento(p, a, &px_a, &dpx_a);
        calcPolinomio_lento(p, b, &px_b, &dpx_b);
    } else {
        calcPolinomio_rapido(p, a, &px_a, &dpx_a);
        calcPolinomio_rapido(p, b, &px_b, &dpx_b);
    }

    if(fabs(px_a) <= ZERO) {
        *raiz = a;
        (*it)++;
        return erro;
    }

    do {
        (*it)++;
        x_ant = x_novo;
        x_novo = (a + b) / 2.0;

        real_t px_m, dpx_m = 0.0;
        if (tipo_calc == 1) calcPolinomio_lento(p, x_novo, &px_m, &dpx_m);
        else calcPolinomio_rapido(p, x_novo, &px_m, &dpx_m);

        if (px_a * px_m < 0) {
            b = x_novo;
            px_b = px_m;
            dpx_b = dpx_m;
        } else {
            a = x_novo;
            px_a = px_m;
            dpx_a = dpx_m;
        }

        // erro de acordo com o critério
        if (criterioParada == 0)
            erro = (fabs(x_novo) > ZERO) ? fabs(x_novo - x_ant)/fabs(x_novo) : fabs(x_novo - x_ant);
        else if (criterioParada == 1)
            erro = fabs(px_m);
        else // ULPS
            erro = (real_t) ulps_diff(x_novo, x_ant);

        // condição de parada
        if ((criterioParada == 0 && erro < EPS) ||
            (criterioParada == 1 && fabs(px_m) <= DBL_EPSILON) ||
            (criterioParada == 2 && ulps_diff(x_novo, x_ant) <= ULPS))
            break;

    } while (*it < MAXIT);

    *raiz = x_novo;
    return erro;
}

// ---------------- CÁLCULO POLINÔMIO ----------------
void calcPolinomio_rapido(Polinomio p, double x, double *px, double *dpx) {
    *px = p.p[p.grau];
    *dpx = 0.0;
    for (int i = p.grau - 1; i >= 0; i--) {
        *dpx = (*dpx) * x + (*px);
        *px = (*px) * x + p.p[i];
    }
}

void calcPolinomio_lento(Polinomio p, double x, double *px, double *dpx) {
    *px = 0.0;
    *dpx = 0.0;
    for (int i = p.grau; i >= 0; i--) {
        *px += p.p[i] * pow(x, i);
        if (i > 0) *dpx += i * p.p[i] * pow(x, i - 1);
    }
}
