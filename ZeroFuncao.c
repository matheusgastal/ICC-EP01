#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"



// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t newtonRaphson (Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz)
{

}


// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t bisseccao (Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz)
{
    double anterior, novo;
    novo = (a + b) / 2;
    real_t resultado_novo = 0;
    real_t resultado_a = 0;
    for(int i = p.grau; i > 0; i--){
        resultado_novo += p.p[i] * pow(novo,i);
        resultado_a += p.p[i] * pow(a,i);
    }

    if(resultado_a * resultado_novo > 0)
        a = novo;
    else if(resultado_a * resultado_novo < 0)
        b = novo;
    else{
        return novo;
    }
}


void calcPolinomio_rapido(Polinomio p, real_t x, real_t *px, real_t *dpx)
{
}   

void calcPolinomio_lento(Polinomio p, real_t x, real_t *px, real_t *dpx)
{

}
