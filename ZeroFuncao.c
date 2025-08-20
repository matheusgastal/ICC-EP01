#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"



// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t newtonRaphson (Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz)
{

}


// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t bisseccao(Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz, int tipo_calc)
{
    double px_a, dpx_a;
    double px_b, dpx_b;
    double px_m, dpx_m;
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
        return 0;
    }

    do {
        (*it)++;
        m = (a + b) / 2.0;

        if (tipo_calc == 1) {
            calcPolinomio_lento(p, m, &px_m, &dpx_m);
        } else {
            calcPolinomio_rapido(p, m, &px_m, &dpx_m);
        }

        if (px_a * px_m < 0) {
            b = m;
            px_b = px_m;
        } else if (px_m * px_b < 0) {
            a = m;
            px_a = px_m;
        } else {
            break; 
        }

    } while (*it < 600 && fabs(px_m) > criterioParada);

    *raiz = m;
    return fabs(px_m); 
}




void calcPolinomio_rapido(Polinomio p, double x, double *px, double *dpx)
{
    *px = p.p[p.grau];  
    *dpx = 0.0;

    // Itera do grau-1 até o termo constante
    for (int i = p.grau - 1; i >= 0; i--) {
        *dpx = (*dpx) * x + (*px);   // regra de Horner para derivada
        *px  = (*px)  * x + p.p[i];  // regra de Horner para polinômio
    }
}


void calcPolinomio_lento(Polinomio p, double x, double *px, double *dpx)
{
    *px = 0.0;   
    *dpx = 0.0;  
    for (int i = p.grau; i >= 0; i--) {
        *px  += p.p[i] * pow(x, i);               
        if (i > 0) {
            *dpx += i * p.p[i] * pow(x, i - 1);  
        }
    }
}
