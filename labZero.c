#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <fenv.h>
#include "utils.h"
#include "ZeroFuncao.h"

int main() {
    fesetround(FE_DOWNWARD); // arredondamento para baixo

    real_t a, b;
    Polinomio pol;

    // Lê grau do polinômio
    scanf("%d", &pol.grau);

    // Aloca vetor de coeficientes
    pol.p = malloc((pol.grau + 1) * sizeof(double));

    // Lê coeficientes
    for (int i = pol.grau; i >= 0; --i)
        scanf("%lf", &pol.p[i]);

    // Lê intervalo
    scanf("%lf %lf", &a, &b);

    int it = 1;
    real_t raiz, erro;
    double t0, t1, tempo;
    int criterio = 0;
    // ---------- RAPIDO ----------
    printf("\nRAPIDO\n\n");
        for(int i = 0; i< 3; i++){
          it = 0; 
          
          t0 = timestamp();
          erro = bisseccao(pol, a, b, i, &it, &raiz, 0); // rápido
          t1 = timestamp();
          tempo = t1 - t0;
          printf("bissec  %+1.15e % .15e %4d % .8f\n", raiz, erro, it, tempo);
        }


      for(int i = 0; i<3; i++){
          it = 0; 
          t0 = timestamp();
          erro = newtonRaphson(pol, (a+b)/2, i, &it, &raiz, 0);
          t1 = timestamp();
          tempo = t1 - t0;
          printf("newton  %+1.14e % .15e %4d % .8f\n", raiz, erro, it, tempo);
      }
          


 
    // ---------- LENTO ----------
    printf("\nLENTO\n\n");

    for(int i = 0; i< 3; i++){
          it = 0;  
          t0 = timestamp();
          erro = bisseccao(pol, a, b, i, &it, &raiz, 1); // rápido
          t1 = timestamp();
          tempo = t1 - t0;
          printf("bissec  %+1.14e % .15e %4d % .8f\n", raiz, erro, it, tempo);
          
         
        } 
  for(int i = 0; i<3; i++){
       it = 0; 
       t0 = timestamp();
       erro = newtonRaphson(pol, (a+b)/2, i, &it, &raiz, 1);
       t1 = timestamp();
       tempo = t1 - t0;
       printf("newton  %+1.14e % .15e %4d % .8f\n", raiz, erro, it, tempo);
  }
    free(pol.p);
    return 0;
}
