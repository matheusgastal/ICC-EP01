#!/bin/bash

# Arquivo de entrada
INPUT="teste.dat"  # substitua pelo seu arquivo, se necessário

# Ajuste os núcleos que deseja usar (pinning de threads)
CORES="0-3"
export LIKWID_THREADS=$CORES

# Executa o programa normalmente (opcional, para conferir saída)
./resolveEDO < "$INPUT"

# Executa LIKWID para medir apenas o contador FP_ARITH_INST_RETIRED_SCALAR_DOUBLE
likwid-perfctr -g FLOPS_DP -m -C $CORES --marker GS_ ./resolveEDO < "$INPUT" | \
awk '/FP_ARITH_INST_RETIRED_SCALAR_DOUBLE/ {print $1","$2}'
