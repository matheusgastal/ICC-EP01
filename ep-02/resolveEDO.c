#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fenv.h>
#include "utils.h"
#include <likwid.h>
#include "edo.h"
#include "gs.h"

//Matheus Gastal Magalhães
//GRR20244620

int main(void)
{
    LIKWID_MARKER_INIT;
    fesetround(FE_DOWNWARD);

    int n;
    scanf("%d", &n); 

    double a, b;
    scanf("%lf %lf", &a, &b);

    double ya, yb;
    scanf("%lf %lf", &ya, &yb);

    double p, q;
    scanf("%lf %lf", &p, &q);

    // ler linhas restantes com 4 coeficientes até EOF
    double r1, r2, r3, r4;
    while ( (scanf("%lf %lf %lf %lf", &r1, &r2, &r3, &r4)) == 4 ) {

        EDo eq;
        eq.n = n;
        eq.a = a; 
        eq.b = b;
        eq.ya = ya; 
        eq.yb = yb;
        eq.p = p; 
        eq.q = q;
        eq.r1 = r1; 
        eq.r2 = r2; 
        eq.r3 = r3; 
        eq.r4 = r4;

        // imprime SL (matriz aumentada)
        prnEDOsl(&eq);

        // gerar estrutura tridiagonal
        Tridiag *sl = genTridiag(&eq);

        // preparar vetor solução
        real_t *x = (real_t *) calloc(n, sizeof(real_t));

        // medir tempo
        rtime_t t0 = timestamp();

        int iters = 0;
        real_t normL2 = 0.0;
        
       
        string_t marker = markerName("GS", n);
        LIKWID_MARKER_START(marker);

        gauss_seidel_tridiag(sl, x, 100, 1e-5, &iters, &normL2);

        LIKWID_MARKER_STOP(marker);
        free(marker);


        rtime_t t_elapsed = timestamp() - t0;

        // imprimir solução (linha única, n elementos)
        for (int i = 0; i < n; ++i) {
            printf(FORMAT, x[i]);
        }
        printf("\n");

        // número de iterações (apenas o inteiro em linha separada)
        printf("%d\n", iters);

        // norma L2 (15 casas)
        printf(FORMAT "\n", normL2);

        // tempo em milisegundos com 8 casas decimais
        printf("%.8e\n\n", t_elapsed);

        // liberar
        free(sl->D);
        free(sl->Di);
        free(sl->Ds);
        free(sl->B);
        free(sl);
        free(x);
    }
    LIKWID_MARKER_CLOSE;
    return 0;
}
