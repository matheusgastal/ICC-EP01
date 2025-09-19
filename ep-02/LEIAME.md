Matheus Gastal Magalhães
GRR: 20244620


Estruturas de dados utilizadas:

EDo: Estrutura que armazena todas as informações da EDO, como o número de pontos internos, o intervalo do problema, as condições de contorno e os coeficientes da EDO genérica. É essencial para ler o input e para posteriormente gerar a matriz tridiagonal que representa o sistema linear.

Tridiag: Estrutura utilizada para representar sistemas lineares tridiagonais. Armazena os vetores da diagonal inferior, diagonal principal, diagonal superior e os termos independentes, além do tamanho da matriz. É gerada a partir da estrutura EDo e serve como entrada para o método de Gauss-Seidel.

Descrição das funções feitas por mim:

int gauss_seidel_tridiag(Tridiag *sl, real_t *x, int maxit, real_t tol, int *iters, real_t *normL2): 
Parâmetros:
*sl que armazena o sistema linear, com diagonal inferior, principal, superior e de termos independentes, além do tamanho do sistema. 
*x é o vetor de soluções que inicialmente é preenchido com zeros.
int maxit que é o número máximo de iterações para o método parar
real_t tol que é a tolerância de convergência do método, que para quando a norma L2 do resíduo é menor ou igual a essa tolerância.
int *iters que é a variável que vai retornar pra main as iterações feitas
real_t *normL2 que é a variavel que vai retornar a norma L2 do resíduo final, ou seja, a precisão da solução.

Explicando a função gauss_seidel_tridiag:
EA função gauss_seidel_tridiag resolve sistemas lineares tridiagonais usando o método de Gauss-Seidel. Ela recebe o sistema sl, que contém a diagonal principal, a diagonal inferior, a diagonal superior e o vetor de termos independentes, além do tamanho do sistema. O vetor x é inicializado com zeros e, a cada iteração, seus elementos são atualizados usando os valores das diagonais e do vetor de termos independentes, aproveitando a estrutura tridiagonal para reduzir operações. Após cada iteração, é calculada a norma L2 do resíduo para verificar se a solução já atingiu a precisão desejada (tol) ou se deve continuar até o máximo de iterações (maxit). Ao final, iters recebe o número de iterações realizadas, normL2 a norma do resíduo, e x contém a solução aproximada do sistema.

int main():
A função main lê os dados do problema: o tamanho do sistema, os limites da EDO, as condições de contorno e os coeficientes da equação. Para cada conjunto de coeficientes lidos, ela cria uma estrutura EDo e imprime a matriz aumentada. Em seguida, gera a matriz tridiagonale aloca dinamicamente o vetor solução x. O método de Gauss-Seidel é aplicado com  a função que descrevi acima, calculando a solução aproximada, o número de iterações realizadas e a norma L2 do resíduo.Depois disso imprime a solução, a quantidade de iterações, a norma L2 e o tempo gasto.


Funções entregues pelo professor:

Tridiag *genTridiag (EDo *edoeq):
Parâmetros:
edoeq: Equação diferencial com todos os dados necessários.

Explicação da função *genTridiag:
A função genTridiag recebe os dados da EDO e monta a estrutura tridiagonal que representa o sistema linear correspondente. Ela cria dinamicamente os vetores da diagonal principal, das diagonais inferior e superior, e o vetor de termos independentes. Cada elemento é calculado usando os coeficientes da EDO e o passo h. As condições de contorno também são aplicadas ajustando os primeiros e últimos elementos de B. No final, ela retorna o ponteiro para a matriz tridiagonal.

void prnEDOsl (EDo *edoeq):

Parâmetros:
*edoeq: Equação diferencial com todos os dados necessários.

Explicação da função prnEDOsl:
Já a função prnEDOsl é responsável por imprimir o sistema linear tridiagonal na saída padrão. Ela percorre cada linha, calcula os valores das diagonais e do vetor de termos independentes, e organiza tudo no formato de matriz aumentada. 




