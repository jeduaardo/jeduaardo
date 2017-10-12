import copy #Para utilizar deepcopy
from random import *

## Parâmetros Ajustáveis
nc = 4 # Número de Colunas
nr = 2 # Número de Linhas
lb = 2 # Levels-Back
ni = 2 # Número de Entradas
no = 2 # Número de Saídas
inputA = [0,0,1,1]
inputB = [0,1,0,1]
max_entradas = 2 # Número de máximos de conexões por porta lógica
ur = 1 # Taxa de mutação (em %)
ee_lambda = 4 #Lambda para estratégias de evolução

## Parâmetros Globais
#FT = [100, 110, 120, 130, 900] # Tabela de Funções
FT = [100, 110]
Ln = nc * nr # Número máximo de Nós
ug = int((Ln * ur) / 100) # Transforma porcentagem em quantidade de genes
if ug == 0:
    ug = 1

## Arrays
G = [] # Array de Genótipos
NU = [] # Array de Nós Ativos/Inativos
aux_AI = [] # Array auxiliar de nós ativos/inativos
GENS_EVOL = [] # Array que armazena todos os genótipos de todos os indivíduos de cada geração (Pai + ee_lambda descendentes)
NOS_EVOL = [] # Array que armazena todos os nós ativos de todos os indivíduos da geração anterior

## Arrays auxiliares
fila = []
myarray = []
MEUARRAY = []
TABELA_V = []
FITNESS_EE = []
ATIVO_MOM = []


#INPUTS = [[0,0,0,0,1,1,1,1],[0,0,1,1,0,0,1,1],[0,1,0,1,0,1,0,1]]
#OUTPUTS = [[1,1,0,0,0,0,0,1],[1,0,1,1,1,0,0,0]]

INPUTS = [[0,0,1,1],[0,1,0,1]]
OUTPUTS = [[0,0,0,1],[1,0,0,0]]

def formato_ee():
    for i in range(ee_lambda + 1):
        GENS_EVOL.append([])
        NOS_EVOL.append([])

def limpa_tabela_verdade():
    TABELA_V.clear()

def formato_tabela_verdade():
    for i in range(Ln + ni + no):
        TABELA_V.append([])

def formato_array_fitness():
    for j in range(ee_lambda + 1):
        FITNESS_EE.append([])
        for k in range(no):
            FITNESS_EE[j].append([])
                        

def gera_formato_nos():
    for i in range(Ln):
        NU.append([])
    for k in range(ee_lambda + 1):
        MEUARRAY.append([])
        for j in range(no):
            MEUARRAY[k].append([])

def gera_formato_aux():
    for i in range(nc):
        aux_AI.append([])
        aux_AI[i] = "-"

def gera_formato_genotipo(): #Define o formato da matriz genótipos G
    for i in range(Ln + ni):
        if i <= (ni -1):
            G.append([]) #Aloca espaços para as entradas
        else:
            G.append([]) #Aloca espaços para os nós
            for k in range(3):
                G[i].append([])
    for i in range(no):
        G.append([]) #Aloca espaço para as saidas

def populacao_inicial():
    for j in range(Ln):
        i = j + ni # Pula as posições das entradas
        n_porta = randint(0,len(FT)-1) # Sorteia uma porta dentre as disponíveis em FT
        porta = FT[n_porta]
        G[i][max_entradas] = porta # A porta ocupa o último locus do gene
        colunaatual = int(j/nr) + 1 # Determina a coluna atual (sem contar entrada)
        if colunaatual == 1: #Neste caso as portas lógicas só recebem inputs
            for k in range(max_entradas):
                entrada = randint(0, (ni-1))
                G[i][k] = entrada
        if colunaatual - lb < 1 and colunaatual != 1: #Neste caso as portas lógicas recebem inputs e nós
            #valorpossivel = (ni-1) + (nr * (colunaatual-2)) + (nr-1) + ni -1 # Determina todos os valores de entradas possiveis
            valorpossivel = (nr * (colunaatual-2)) + (nr-1) + ni
            print(valorpossivel)
            for k in range(max_entradas):
                sorteado = randint(0, valorpossivel)
                G[i][k] = sorteado
        if colunaatual - lb >= 1: #Neste caso as portas logicas recebem somente nós
            primeiro_elemento = nr * (colunaatual-lb-1) + ni # O fator + ni serve para deslocar os elementos para frente, contando a quantidade de inputs.
            ultimo_elemento = (nr * (colunaatual-2)) + (nr-1) + ni
            for k in range(max_entradas):
                #entrada_ou_porta = randint(0, 1) #Determina se a entrada da porta lógica atual será uma entrada do circuito ou a saída de alguma porta lógica já existente (0: entrada; 1: porta lógica)
                elementos_portas = ultimo_elemento - primeiro_elemento + 1 #Determina quantas portas possíveis existem
                qtd_elementos = elementos_portas + ni - 1 #Determina a quantidade de valores distintos que a porta pode assumir como entrada. O -1 tem como função apenas possibilitar o uso de randint(0, qtd_elementos)
                entrada_ou_porta = randint(0, qtd_elementos) #Dá a mesma chance para todos os possiveis inputs
                if entrada_ou_porta >= ni:
                    G[i][k] = randint(primeiro_elemento, ultimo_elemento)
                else:
                    G[i][k] = randint(0, ni-1)

def define_saida():
    ultimo_elemento = ni + Ln - 1 #Último elemento possivel que a saída pode assumir
    primeiro_elemento = ultimo_elemento - (lb * nr) + 1 #Primeiro valor possivel que a saída pode assumir
    elementos_portas = ultimo_elemento - primeiro_elemento + 1 #Quantidade de valores existentes entre o primeiro e último possiveis valores
    qtd_elementos = elementos_portas + ni - 1 #Determina a quantidade de valores distintos que a porta pode assumir como entrada. O -1 tem como função apenas possibilitar o uso de randint(0, qtd_elementos)
    for i in range(no):
        entrada_ou_porta = randint(0, qtd_elementos) #Dá a mesma chance para todos os possiveis inputs
        if entrada_ou_porta >= ni: #G[Ln+ni+i] representa cada saida, sequencialmente
            G[Ln+ni+i] = randint(primeiro_elemento, ultimo_elemento)
        else:
            G[Ln+ni+i] = randint(0, ni-1)

def nos_ativos():
    fila = [] #Array que funcionará como uma fila
    myarray = [] #Array que armazenará os genes dos nós da fila
    vetor = []
    for k in range(Ln+ni):
        vetor.append("-")
    #vetor = ["-", "-", "-", "-", "-", "-", "-", "-"] #Este vetor temporário (existe somente nesta função), recebe "X" quando o nó é ativo e "-" quando o nó é inativo.
    for w in range(ee_lambda + 1): #Percorre o pai e seus lambda filhos
        for i in range(no): #Realiza o procedimento para o número de saídas
            saida = int(Ln+ni+i) #Determina a saída a ser processada
            fila.append(GENS_EVOL[w][saida])
            noatual = GENS_EVOL[w][GENS_EVOL[w][saida]] #Obtém as informações do gene responsável pela saída
            if type(noatual) is not int: #len(noatual) == 3: #Analisa se a saída é o resultado de um nó ou é uma entrada do programa. (3 significa: entrada, entrada e porta)
                fila.append(noatual[0]) #Adiciona primeira entrada do nó saída
                fila.append(noatual[1]) #Adiciona segunda entrada do nó saída
                myarray.append(GENS_EVOL[w][noatual[0]]) #Adiciona o gene responsável pela primeira entrada
                myarray.append(GENS_EVOL[w][noatual[1]]) #Adiciona o gene responsável pela segunda entrada
            else:
                fila.append(noatual) #Se tamanho é diferente de 3, então é uma entrada. Logo a fila recebe o valor da posição da entrada.
            while len(myarray) != 0: #Procedimento para esvaziar a fila
                noatual = myarray[0]
                if type(noatual) is int: #(len(noatual) != 3: #Considerando tamanho 3 do gene, isto é, duas entradas mais a função
                    fila.append("entrada")
                    myarray.pop(0)
                else:
                    fila.append(noatual[0])
                    fila.append(noatual[1])
                    myarray.append(GENS_EVOL[w][noatual[0]])
                    myarray.append(GENS_EVOL[w][noatual[1]])
                    myarray.pop(0)
            #print("TAMANHO FILA:", len(fila))
            for j in range(len(fila)): # Coloca X nas posições ativas do genótipo
                if type(fila[j]) is int:
                    if vetor[fila[j]] == "-":
                        vetor[fila[j]] = "X"
            print(fila)
            print(vetor)
            MEUARRAY[w][i] = copy.deepcopy(vetor)
            vetor.clear()
            for k in range(Ln+ni):
                vetor.append("-")
            #vetor = ["-", "-", "-", "-", "-", "-", "-", "-"] #Retorna o formato padrão
            fila = [] #Limpa o vetor fila
            myarray = [] #Limpa o vetor myarray
               
def mutacao_pontual():
    for k in range(ee_lambda): #Realiza mutações e gera "ee_lambda" variações do pai
        print("ESTAMOS EM K IGUAL A", k)
        #print("GENOTIPO:", G)
        qtd_mutacao = randint(1, ug) #Determina quantos nós serão mutados, dentro do intervalo recomendado pela literatura
        for i in range(qtd_mutacao): #Realiza "qtd_mutacao" de mutações NESTE descendente
            G = copy.deepcopy(GENS_EVOL[0]) #Busca o pai original de volta (sempre alocado em GENS_EVOL[0])
            gene_mutado = randint(ni, (ni+no+Ln-1)) #O gene que será mutado pode receber qualquer valor maior que o último espaço ocupado pela entrada, isto é da posição ni (já que começa por 0) até a posição da última saida
            print(gene_mutado)
            if gene_mutado >= (Ln+ni): #Garante que o gene selecionado é uma saída
                ultimo_elemento = ni + Ln - 1 #Último elemento possivel que a saída pode assumir
                primeiro_elemento = ultimo_elemento - (lb * nr) + 1 #Primeiro valor possivel que a saída pode assumir
                elementos_portas = ultimo_elemento - primeiro_elemento + 1 #Quantidade de valores existentes entre o primeiro e último possiveis valores
                qtd_elementos = elementos_portas + ni - 1 #Determina a quantidade de valores distintos que a porta pode assumir como entrada. O -1 tem como função apenas possibilitar o uso de randint(0, qtd_elementos)
                entrada_ou_porta = randint(0, qtd_elementos) #Dá a mesma chance para todos os possiveis inputs
                if entrada_ou_porta >= ni: #G[Ln+ni+i] representa cada saida, sequencialmente
                    nova_saida = randint(primeiro_elemento, ultimo_elemento)
                    while G[gene_mutado] == nova_saida: #Garante que a mutação mudará o gene_mutado para um valor diferente do inicial
                        nova_saida = randint(primeiro_elemento, ultimo_elemento)
                    G[gene_mutado] = nova_saida
                else:
                    nova_saida = randint(0, ni-1)
                    while G[gene_mutado] == nova_saida: #Garante que a mutação mudará o gene_mutado para um valor diferente do inicial
                        nova_saida = randint(0, ni-1)
                    G[gene_mutado] = nova_saida    
            else:
                colunaatual = int((gene_mutado - ni)/nr) + 1 #Determina a coluna na qual o gene selecionado está
                alelo_mutado = randint(0, max_entradas) #Determina qual será o alelo a ser mutado
                print("Mutação em nó", gene_mutado, " alelo ", alelo_mutado, " coluna ", colunaatual)
                print("Nó correspondente: ", G[gene_mutado])
    
                if alelo_mutado == max_entradas: #Significa que a mutação ocorrerá no tipo da porta
                    porta_logica = randint(0, (len(FT)-1)) #Determina qual será a nova porta
                    while G[gene_mutado][alelo_mutado] == FT[porta_logica]: #Garante que a mutação mudará o gene_mutado para um valor diferente do inicial
                        porta_logica = randint(0, (len(FT)-1))
                    G[gene_mutado][alelo_mutado] = FT[porta_logica] #Atribui a nova porta ao gene
                else:
                    if colunaatual == 1: #Neste caso as portas lógicas só recebem inputs
                        entrada = randint(0, (ni-1))
                        while G[gene_mutado][alelo_mutado] == entrada:#Garante que a mutação mudará o gene_mutado para um valor diferente do inicial
                            entrada = randint(0, (ni-1))
                        G[gene_mutado][alelo_mutado] = entrada

                    
                    if colunaatual - lb < 1 and colunaatual != 1: #Neste caso as portas lógicas recebem inputs e nós       
                        valorpossivel = (nr * (colunaatual-2)) + (nr-1) + ni
                        sorteado = randint(0, valorpossivel)
                        while G[gene_mutado][alelo_mutado] == sorteado:#Garante que a mutação mudará o gene_mutado para um valor diferente do inicial
                            sorteado = randint(0, valorpossivel)
                        G[gene_mutado][alelo_mutado] = sorteado

            
                    if colunaatual - lb >= 1: #Neste caso as portas logicas recebem somente nós
                        primeiro_elemento = nr * (colunaatual-lb-1) + ni # O fator + ni serve para deslocar os elementos para frente, contando a quantidade de inputs.
                        ultimo_elemento = (nr * (colunaatual-2)) + (nr-1) + ni 
                        elementos_portas = ultimo_elemento - primeiro_elemento + 1 #Determina quantas portas possíveis existem
                        qtd_elementos = elementos_portas + ni - 1 #Determina a quantidade de valores distintos que a porta pode assumir como entrada. O -1 tem como função apenas possibilitar o uso de randint(0, qtd_elementos)
                        entrada_ou_porta = randint(0, qtd_elementos) #Dá a mesma chance para todos os possiveis inputs
                        if entrada_ou_porta >= ni:
                            novo_valor_mutacao = randint(primeiro_elemento, ultimo_elemento)
                            while G[gene_mutado][alelo_mutado] == novo_valor_mutacao:#Garante que a mutação mudará o gene_mutado para um valor diferente do inicial
                                novo_valor_mutacao = randint(primeiro_elemento, ultimo_elemento)
                            G[gene_mutado][alelo_mutado] = novo_valor_mutacao
                        else:
                            novo_valor_mutacao = randint(0, ni-1)
                            while G[gene_mutado][alelo_mutado] == novo_valor_mutacao:#Garante que a mutação mudará o gene_mutado para um valor diferente do inicial
                                novo_valor_mutacao = randint(0, ni-1)
                            G[gene_mutado][alelo_mutado] = novo_valor_mutacao
    
                print("Nó correspondente: ", G[gene_mutado])

        GENS_EVOL[k+1] = copy.deepcopy(G) #Passa o genótipo mutado para a matriz de genótipos
        print("EVOL K+1", GENS_EVOL[k+1])
        #atualiza_G() #Retorna com o pai para sofrer nova mutação


def SAM_ativo(posicao):
    fila = [] #Array que funcionará como uma fila
    myarray = [] #Array que armazenará os genes dos nós da fila
    vetor = []
    for k in range(Ln+ni):
        vetor.append("-")
    #vetor = ["-", "-", "-", "-", "-", "-", "-", "-"] #Este vetor temporário (existe somente nesta função), recebe "X" quando o nó é ativo e "-" quando o nó é inativo.
    for w in range(1): #Percorre o pai e seus lambda filhos
        for i in range(no): #Realiza o procedimento para o número de saídas
            saida = int(Ln+ni+i) #Determina a saída a ser processada
            fila.append(GENS_EVOL[posicao][saida])
            noatual = GENS_EVOL[posicao][GENS_EVOL[posicao][saida]] #Obtém as informações do gene responsável pela saída
            if type(noatual) is not int: #len(noatual) == 3: #Analisa se a saída é o resultado de um nó ou é uma entrada do programa. (3 significa: entrada, entrada e porta)
                fila.append(noatual[0]) #Adiciona primeira entrada do nó saída
                fila.append(noatual[1]) #Adiciona segunda entrada do nó saída
                myarray.append(GENS_EVOL[posicao][noatual[0]]) #Adiciona o gene responsável pela primeira entrada
                myarray.append(GENS_EVOL[posicao][noatual[1]]) #Adiciona o gene responsável pela segunda entrada
            else:
                fila.append(noatual) #Se tamanho é diferente de 3, então é uma entrada. Logo a fila recebe o valor da posição da entrada.
            while len(myarray) != 0: #Procedimento para esvaziar a fila
                noatual = myarray[0]
                if type(noatual) is int: #(len(noatual) != 3: #Considerando tamanho 3 do gene, isto é, duas entradas mais a função
                    fila.append("entrada")
                    myarray.pop(0)
                else:
                    fila.append(noatual[0])
                    fila.append(noatual[1])
                    myarray.append(GENS_EVOL[posicao][noatual[0]])
                    myarray.append(GENS_EVOL[posicao][noatual[1]])
                    myarray.pop(0)
            #print("TAMANHO FILA:", len(fila))
            for j in range(len(fila)): # Coloca X nas posições ativas do genótipo
                if type(fila[j]) is int:
                    if vetor[fila[j]] == "-":
                        vetor[fila[j]] = "X"
            print(fila)
            print(vetor)
            MEUARRAY[posicao][i] = copy.deepcopy(vetor)
            vetor.clear()
            for k in range(Ln+ni):
                vetor.append("-")
            #vetor = ["-", "-", "-", "-", "-", "-", "-", "-"] #Retorna o formato padrão
            fila = [] #Limpa o vetor fila
            myarray = [] #Limpa o vetor myarray
    

def mutacao_SAM(posicao):
    ATIVO_MOM = copy.deepcopy(MEUARRAY[posicao])
    while ATIVO_MOM == MEUARRAY[posicao]:
        #FICA FAZENDO MUTAÇÃO
        
        qtd_mutacao = randint(1, ug) #Determina quantos nós serão mutados, dentro do intervalo recomendado pela literatura
        for i in range(qtd_mutacao): #Realiza "qtd_mutacao" de mutações NESTE descendente
            G = copy.deepcopy(GENS_EVOL[0]) #Busca o pai original de volta (sempre alocado em GENS_EVOL[0])
            gene_mutado = randint(ni, (ni+no+Ln-1)) #O gene que será mutado pode receber qualquer valor maior que o último espaço ocupado pela entrada, isto é da posição ni (já que começa por 0) até a posição da última saida
            print(gene_mutado)
            if gene_mutado >= (Ln+ni): #Garante que o gene selecionado é uma saída
                ultimo_elemento = ni + Ln - 1 #Último elemento possivel que a saída pode assumir
                primeiro_elemento = ultimo_elemento - (lb * nr) + 1 #Primeiro valor possivel que a saída pode assumir
                elementos_portas = ultimo_elemento - primeiro_elemento + 1 #Quantidade de valores existentes entre o primeiro e último possiveis valores
                qtd_elementos = elementos_portas + ni - 1 #Determina a quantidade de valores distintos que a porta pode assumir como entrada. O -1 tem como função apenas possibilitar o uso de randint(0, qtd_elementos)
                entrada_ou_porta = randint(0, qtd_elementos) #Dá a mesma chance para todos os possiveis inputs
                if entrada_ou_porta >= ni: #G[Ln+ni+i] representa cada saida, sequencialmente
                    nova_saida = randint(primeiro_elemento, ultimo_elemento)
                    while G[gene_mutado] == nova_saida: #Garante que a mutação mudará o gene_mutado para um valor diferente do inicial
                        nova_saida = randint(primeiro_elemento, ultimo_elemento)
                    G[gene_mutado] = nova_saida
                else:
                    nova_saida = randint(0, ni-1)
                    while G[gene_mutado] == nova_saida: #Garante que a mutação mudará o gene_mutado para um valor diferente do inicial
                        nova_saida = randint(0, ni-1)
                    G[gene_mutado] = nova_saida    
            else:
                colunaatual = int((gene_mutado - ni)/nr) + 1 #Determina a coluna na qual o gene selecionado está
                alelo_mutado = randint(0, max_entradas) #Determina qual será o alelo a ser mutado
                print("Mutação em nó", gene_mutado, " alelo ", alelo_mutado, " coluna ", colunaatual)
                print("Nó correspondente: ", G[gene_mutado])
    
                if alelo_mutado == max_entradas: #Significa que a mutação ocorrerá no tipo da porta
                    porta_logica = randint(0, (len(FT)-1)) #Determina qual será a nova porta
                    while G[gene_mutado][alelo_mutado] == FT[porta_logica]: #Garante que a mutação mudará o gene_mutado para um valor diferente do inicial
                        porta_logica = randint(0, (len(FT)-1))
                    G[gene_mutado][alelo_mutado] = FT[porta_logica] #Atribui a nova porta ao gene
                else:
                    if colunaatual == 1: #Neste caso as portas lógicas só recebem inputs
                        entrada = randint(0, (ni-1))
                        while G[gene_mutado][alelo_mutado] == entrada:#Garante que a mutação mudará o gene_mutado para um valor diferente do inicial
                            entrada = randint(0, (ni-1))
                        G[gene_mutado][alelo_mutado] = entrada

                    
                    if colunaatual - lb < 1 and colunaatual != 1: #Neste caso as portas lógicas recebem inputs e nós       
                        valorpossivel = (nr * (colunaatual-2)) + (nr-1) + ni
                        sorteado = randint(0, valorpossivel)
                        while G[gene_mutado][alelo_mutado] == sorteado:#Garante que a mutação mudará o gene_mutado para um valor diferente do inicial
                            sorteado = randint(0, valorpossivel)
                        G[gene_mutado][alelo_mutado] = sorteado

            
                    if colunaatual - lb >= 1: #Neste caso as portas logicas recebem somente nós
                        primeiro_elemento = nr * (colunaatual-lb-1) + ni # O fator + ni serve para deslocar os elementos para frente, contando a quantidade de inputs.
                        ultimo_elemento = (nr * (colunaatual-2)) + (nr-1) + ni 
                        elementos_portas = ultimo_elemento - primeiro_elemento + 1 #Determina quantas portas possíveis existem
                        qtd_elementos = elementos_portas + ni - 1 #Determina a quantidade de valores distintos que a porta pode assumir como entrada. O -1 tem como função apenas possibilitar o uso de randint(0, qtd_elementos)
                        entrada_ou_porta = randint(0, qtd_elementos) #Dá a mesma chance para todos os possiveis inputs
                        if entrada_ou_porta >= ni:
                            novo_valor_mutacao = randint(primeiro_elemento, ultimo_elemento)
                            while G[gene_mutado][alelo_mutado] == novo_valor_mutacao:#Garante que a mutação mudará o gene_mutado para um valor diferente do inicial
                                novo_valor_mutacao = randint(primeiro_elemento, ultimo_elemento)
                            G[gene_mutado][alelo_mutado] = novo_valor_mutacao
                        else:
                            novo_valor_mutacao = randint(0, ni-1)
                            while G[gene_mutado][alelo_mutado] == novo_valor_mutacao:#Garante que a mutação mudará o gene_mutado para um valor diferente do inicial
                                novo_valor_mutacao = randint(0, ni-1)
                            G[gene_mutado][alelo_mutado] = novo_valor_mutacao
    
                print("Nó correspondente: ", G[gene_mutado])
        #print("EXECUTANDO")
        SAM_ativo(posicao)
        GENS_EVOL[posicao] = copy.deepcopy(G) #Passa o genótipo mutado para a matriz de genótipos

    print("EVOL K+1", GENS_EVOL[posicao])
    G = copy.deepcopy(GENS_EVOL[0])
    #atualiza_G() #Retorna com o pai para sofrer nova mutação    



def atualiza_G(): #Apenas atualiza o genótipo para o "Pai" da geração corrente (buscando ele no array GENS_EVOL (sempre na posição 0))
    #G = []
    G = copy.deepcopy(GENS_EVOL[0])
    print("COPIANDO DE VOLTA:", G)

def monta_tabela(posicao):
    auxiliar = []
    for i in range(no):
        print("ARRAY: ", MEUARRAY[posicao][i])
        for j in range(len(MEUARRAY[posicao][i])):
            if MEUARRAY[posicao][i][j] == "X":
                print("ESTA POSIÇÃO EQUIVALE A: ", GENS_EVOL[posicao][j])
                print("JOTA É: ", j)
                if j < ni:
                    print("É UMA ENTRADA")
                    TABELA_V[j] = copy.deepcopy(INPUTS[j])
                else:
                    print("NÃO É UMA ENTRADA")
                    entradaA = GENS_EVOL[posicao][j][0]
                    entradaB = GENS_EVOL[posicao][j][1]
                    funcao = GENS_EVOL[posicao][j][2]
                    print("A, B, FUNCAO: ", entradaA, entradaB, funcao)
                    if funcao == 100:
                        for w in range(len(TABELA_V[entradaA])):
                            auxiliar.append((TABELA_V[entradaA][w]) and (TABELA_V[entradaB][w]))
                        TABELA_V[j] = copy.deepcopy(auxiliar)
                        auxiliar = []
                    if funcao == 110:
                        for w in range(len(TABELA_V[entradaA])):
                            auxiliar.append((TABELA_V[entradaA][w]) or (TABELA_V[entradaB][w]))
                        TABELA_V[j] = copy.deepcopy(auxiliar)
                        auxiliar = []

        saida_ckt = GENS_EVOL[posicao][Ln+ni+i] #Determina a saida atual do circuito, que está sendo processada
        tabela_saida = TABELA_V[saida_ckt]
        if len(TABELA_V[saida_ckt]) == 0:
            print("TABELA VAZIA")
        fitness_atual = 0
        for s in range(len(OUTPUTS[i])):
            if OUTPUTS[i][s] == TABELA_V[saida_ckt][s]:
                print("OUTPUT IGUAL")
                fitness_atual = fitness_atual + 1
        FITNESS_EE[posicao][i] = copy.deepcopy(fitness_atual)
        TABELA_V.clear()
        formato_tabela_verdade()

def atualiza_pai():
    soma = 0
    fitness_soma = []
    aux = []
    for k in range(ee_lambda + 1):
        fitness_soma.append([])
        aux.append("-")
    for i in range(len(FITNESS_EE)):
        for j in range(len(FITNESS_EE[i])):
            soma = soma + FITNESS_EE[i][j]
        fitness_soma[i] = copy.deepcopy(soma)
        soma = 0
    print(fitness_soma)



def main():
    formato_ee()
    formato_array_fitness()
    gera_formato_nos()
    formato_tabela_verdade()
    gera_formato_genotipo()
    populacao_inicial()
    define_saida()
    for j in range(ni):
        G[j] = copy.deepcopy(INPUTS[j])
    GENS_EVOL[0] = copy.deepcopy(G)        
    for i in range(5):
        mutacao_pontual()
        nos_ativos()
        for w in range(ee_lambda + 1):
            monta_tabela(w)
        print(FITNESS_EE)


def main_teste():
    formato_ee()
    formato_array_fitness()
    gera_formato_nos()
    formato_tabela_verdade()
    gera_formato_genotipo()
    populacao_inicial()
    define_saida()
    for j in range(ni):
        G[j] = copy.deepcopy(j)
        #G[j] = copy.deepcopy(INPUTS[j])
    GENS_EVOL[0] = copy.deepcopy(G)
    for i in range(1):
        mutacao_pontual()
        nos_ativos()
        for w in range(ee_lambda + 1):
            monta_tabela(w)
        print(FITNESS_EE)


def main_mutacao_SAM():
    formato_ee()
    formato_array_fitness()
    gera_formato_nos()
    formato_tabela_verdade()
    gera_formato_genotipo()
    populacao_inicial()
    define_saida()
    for j in range(ni):
        G[j] = copy.deepcopy(j)
    GENS_EVOL[0] = copy.deepcopy(G)
    mutacao_pontual()
    nos_ativos()
    for i in range(5):
        for k in range(ee_lambda):
            mutacao_SAM(k+1)


        

