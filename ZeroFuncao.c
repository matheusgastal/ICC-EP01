#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>

#include "utils.h"
#include "ZeroFuncao.h"

// Função para calcular a diferença em ULPs entre dois doubles
long ulps_diff(double a, double b) {
    union { double d; long long i; } ua, ub;
    ua.d = a;
    ub.d = b;
    return llabs(ua.i - ub.i);
}

// ---------------- NEWTON-RAPHSON ----------------
real_t newtonRaphson(Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz, int tipo_calc)
{
    double px, dpx;
    double erro = 0.0;
    *it = 0;
    *raiz = x0;

    do {
        (*it)++;
        double x_ant = x0;

        if (tipo_calc == 1)
            calcPolinomio_lento(p, x0, &px, &dpx);
        else
            calcPolinomio_rapido(p, x0, &px, &dpx);

        if (fabs(dpx) < 1e-12) break;

        x0 = x0 - px / dpx;
        *raiz = x0;

        erro = fabs(x0 - x_ant); // erro acumulado

        // Critérios de parada
        if (criterioParada == EPS && fabs(x0 - x_ant) < EPS) return erro;
        if (criterioParada == DBL_EPSILON && fabs(x0 - x_ant) < DBL_EPSILON) return erro;
        if (criterioParada == ULPS && ulps_diff(x0, x_ant) <= ULPS) return erro;

    } while (*it < MAXIT);

    return erro; // retorna erro acumulado final
}

// ---------------- BISSECCAO ----------------
real_t bisseccao(Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz, int tipo_calc)
{
    double px_a, dpx_a, px_b, dpx_b, px_m, dpx_m;
    double m;
    *it = 0;

    if (tipo_calc == 1) {
        calcPolinomio_lento(p, a, &px_a, &dpx_a);
        calcPolinomio_lento(p, b, &px_b, &dpx_b);
    } else {
        calcPolinomio_rapido(p, a, &px_a, &dpx_a);
        calcPolinomio_rapido(p, b, &px_b, &dpx_b);
    }

    if (px_a * px_b > 0) {
        printf("Intervalo invalido\n");
        *raiz = a;
        return 0;
    }

    double erro = b - a; // erro inicial é tamanho do intervalo

    do {
        (*it)++;
        double m_ant = m;
        m = (a + b) / 2.0;

        if (tipo_calc == 1)
            calcPolinomio_lento(p, m, &px_m, &dpx_m);
        else
            calcPolinomio_rapido(p, m, &px_m, &dpx_m);

        *raiz = m;

        if (px_a * px_m < 0) {
            b = m;
            px_b = px_m;
        } else if (px_m * px_b < 0) {
            a = m;
            px_a = px_m;
        } else {
            break;
        }

        erro = fabs(b - a); // erro acumulado

        if (criterioParada == EPS && fabs(m_ant - m) < EPS) return erro;
        if (criterioParada == DBL_EPSILON && fabs(m_ant - m) < DBL_EPSILON) return erro;
        if (criterioParada == ULPS && ulps_diff(m_ant, m)  <= ULPS) return erro;
    } while (*it < MAXIT);

    return erro;
}


// ---------------- CÁLCULO POLINÔMIO ----------------
void calcPolinomio_rapido(Polinomio p, double x, double *px, double *dpx) {
    *px = p.p[p.grau];  
    *dpx = 0.0;
    for (int i = p.grau - 1; i >= 0; i--) {
        *dpx = (*dpx) * x + (*px);   
        *px  = (*px)  * x + p.p[i];  
    }
}

void calcPolinomio_lento(Polinomio p, double x, double *px, double *dpx) {
    *px = 0.0;   
    *dpx = 0.0;  
    for (int i = p.grau; i >= 0; i--) {
        *px  += p.p[i] * pow(x, i);               
        if (i > 0) {
            *dpx += i * p.p[i] * pow(x, i - 1);  
        }
    }
}
