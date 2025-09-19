# PROGRAMA
PROG = resolveEDO
OBJS = $(PROG).o edo.o utils.o gs.o
VERIF = verificaEP02

# Compilador
CC     = gcc

# Acrescentar onde apropriado as opções para incluir uso da biblioteca LIKWID
CFLAGS = -O0 -std=c11 -Wall -DLIKWID_PERFMON -I${LIKWID_INCLUDE}
LFLAGS = -lm -lrt -L${LIKWID_LIB} -llikwid

# Lista de arquivos para distribuição. Acrescentar mais arquivos se necessário.
DISTFILES = *.c *.h LEIAME* Makefile *.dat
DISTDIR = ${USER}

.PHONY: clean purge dist all

# Alvo principal
all: $(PROG) $(VERIF)

# Regras de compilação de objetos
%.o: %.c %.h utils.h
	$(CC) -c $(CFLAGS) -o $@ $<

# Linkagem do programa principal
$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

# Linkagem do verificador
$(VERIF): $(VERIF).c
	$(CC) -Wno-format -o $@ $<

# Teste de formato
testeFormato: $(PROG) $(VERIF)
	@cat teste.dat | ./$(PROG) | ./$(VERIF)

# Limpeza básica
clean:
	@echo "Limpando sujeira ..."
	@rm -f *~ *.bak

# Limpeza completa
purge:  clean
	@echo "Limpando tudo ..."
	@rm -f core a.out $(OBJS)
	@rm -f $(PROG) $(VERIF) $(DISTDIR) $(DISTDIR).tgz

# Gerar arquivo de distribuição
dist: purge
	@echo "Gerando arquivo de distribuição ($(DISTDIR).tgz) ..."
	@ln -s . $(DISTDIR)
	@tar -chzvf $(DISTDIR).tgz $(addprefix ./$(DISTDIR)/, $(DISTFILES))
	@rm -f $(DISTDIR)
