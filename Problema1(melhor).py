# -*- coding: utf-8 -*-
"""
PROBLEMA 1, Python (spyder)

@author: José Rodrigues, nº 2019246536
"""

from math import exp, sqrt, pi, factorial as fct
from scipy.special import assoc_laguerre
import matplotlib.pyplot as plt
import numpy as np

n=3 #nº quântico principal
l=[0,1,2]
z=1 #nº atómico
a=0.52917721092 #raio de bohr: a0 [a]
e=-1 #carga do eletrão [eV]

def R(r,n,l):
    y=sqrt((2/(n*a))**3*fct(n-l-1)/(2*n*fct(n+l)))\
        *exp(-r/(n*a))*(2*r/(n*a))**l\
        *assoc_laguerre((2*r/(n*a)),n-l-1,2*l+1)
    return y


"""
1.
"""


r_a=[i for i in np.arange(0.1,30,0.1)] #r/a0
R30=[]
R31=[]
R32=[]
Rnl=[R30,R31,R32]

for i in l:
    for j in r_a:
        Rnl[i]+=[(4*pi*(j*a)**2*(R(j*a,n,i))**2)]

plt.suptitle('DENSIDADE DE PROBABILIDADE RADIAL DO HIDROGÉNIO PARA n=3')
plt.plot(r_a,R30, label='l=0')
plt.legend(loc='upper right')
plt.plot(r_a,R31, label='l=1',color='red')
plt.legend(loc='upper right')
plt.plot(r_a,R32, label='l=2',color='green')
plt.legend(loc='upper right')
plt.xlabel('r/a0')
plt.ylabel('4\u03C0r^2*|R_3l|^2')
plt.show()


'COMENTÁRIO'
print("\n\n")

print('"""1."""\n')
print('     A partir do gráfico da densidade de probabilidade radial do hidrog\
énio, verifica-se que existem nodos (pontos no qual a densidade de probabilida\
de radial é nula) entre pontos de densidade máxima, sendo que, para l=0 existem 2 nod\
os entre 3 máximos, para l=1 existe 1 nodo entre 2 máximos, e l=3 não tem nodo\
s mas tem um máximo. Verificam-se também pontos nos quais a densidade de proba\
bilidade radial é igual para 2 estados degenerados diferentes.')

print("\n\n")



"""
2.
"""

v1=[]
v2=[]
v=[]
for i in r_a:
    v2+=[-11*e**2/(i*a)]
    v1+=[-e**2/(i*a)]
for i in r_a:
    if i>=4:
        z=1
    else:
        z=11
    v+=[-z*e**2/(i*a)]
    
fig, axs = plt.subplots(2, sharex= True)
fig.suptitle('Potencial do Sódio (Z=11)')
axs[0].plot(r_a,v1,label='Z=1')
axs[0].legend(loc='lower right')
axs[0].plot(r_a,v2,label='Z=11')
axs[0].legend(loc='lower right')
axs[0].set(ylabel='V(r) [V]')
axs[1].plot(r_a,v,label='Z=11, r<4a0; Z=1, r>=4a0', color='g')
axs[1].legend(loc='lower right')
axs[1].set(xlabel='r/a0')
axs[1].set(ylabel='V(r) [V]')
plt.show()

import matplotlib.pyplot as plt
"""
3.
"""
print('"""3."""\n')
print('     O estado do eletrão exterior no átomo de sódio menos afetado pelo \
facto de o potencial não ser um potencial puro de Coulomb é o estado 3d, uma v\
ez que a alteração do potencial ocorre em r=4a0 e este estado é definido apenas p\
or um máximo perto de r=9a0 (zona na qual o potencial não sofre qualquer altera\
ção). Por outro lado, o estado mais afetado por esta condição do potencial é o estado 3s, pois\
  neste estado a densidade de probabilidade radial atinge um máximo local perto de r=a\
0 (no qual está sujeito a um potencial de carga +11), perto de r=4.5a0 atinge \
o seu 2º máximo local e na vizinhança de r=13a0 atinge o seu máximo absoluto (\
sendo y=[i for i in np.arage(0,)]que nestes dois pontos a energia potencial do eletrão é consequente apen\
as de uma carga +1).')

print('\n\n')



"""
4.
"""
print('"""4."""\n')
print('     Os valores indicados para a energia de ligação do eletrão nos es\
tados 3s, 3p e 3d, vão de encontro às conclusões acerca da penetração relativa\
  dos estados, já que o estado onde o eletrão tem maior energia de ligação é o e\
stado no qual, na alínea anterior, se verificou uma zona de densidade de probabilid\
ade radial máxima para um potencial devido à carga total do núcleo de Na (Z=+11), o estado 3s (dest\
es 3 estados, é o estado mais penetrante).\n    Estes valores de energia diminuem com l u\
ma vez que os estados são sucessivamente menos penetrantes, pois o estado 3p ainda \
tem uma zona de densidade de probabilidade radial máxima perto de r=3a0 (cf. gráfico potencial do sódio), pelo \
que a enegia de ligação do eletrão nesse estado vai ser maior que no estado 3d.')

print('\n\n')

"""
5.
"""
print('"""5."""\n')

Eh=-13.60569301 #[eV]
En=[5.12, 2.10, 1.50] #[eV]

def zef(n,En):
    y=sqrt(-n**2*En/Eh)
    return y

j=0
for i in En:
    n=3
    orbital=['3s','3p','3d']
    print(f"Z_eff para {orbital[j]} = {zef(n,i)}")
    j+=1

print('\n')
print("     Estes resultados vêm reforçar as conclusões das alíneas anteriores,\
  a energia de ligação eletrão-núcleo diminui com o nº quântico orbital, \
o que só pode ser explicado pelo facto do potencial efetivo diminuir também com l, \
que é exatamente o que se verifica nesta alínea.")

print('\n\n')



"""
6.1
"""
z=37
n=5

r_a=[i for i in np.arange(0.1,80,0.1)] 
R50=[]
R51=[]
R52=[]
Rnl=[R50,R51,R52]

for i in l:
    for j in r_a:
        Rnl[i]+=[(4*pi*(j*a)**2*(R(j*a,n,i))**2)]

plt.suptitle('DENSIDADE DE PROBABILIDADE RADIAL DO HIDROGÉNIO PARA n=5')
plt.plot(r_a,R50, label='l=0')
plt.legend(loc='upper right')
plt.plot(r_a,R51, label='l=1',color='red')
plt.legend(loc='upper right')
plt.plot(r_a,R52, label='l=2',color='green')
plt.legend(loc='upper right')
plt.xlabel('r/a0')
plt.ylabel('4\u03C0r^2*|R_3l|^2')
plt.show()

'COMENTÁRIO'
print('"""6."""\n')
print('     Tendo em conta ambos os gráficos da densidade de probabilidade rad\
ial do hidrogénio, verfica-se que para r<4a0 o comportamento das curvas (à exc\
eção das amplitudes dos máximos) é muito semelhante para ambos os casos, n=3 e n=5, pe\
lo que a análise feita na alínea 1. é também válida para esta situação (fazend\
o um paralelismo entre os estados, 3s e 5s, 3p e 5p, e 3d e 5d).')

print('\n\n')

"""
6.2
"""
v1=[]
v2=[]
v=[]
for i in r_a:
    v2+=[-z*e**2/(i*a)]
    v1+=[-e**2/(i*a)]
for i in r_a:
    if i>=4:
        z=1
    else:
        z=37
    v+=[-z*e**2/(i*a)]
    
fig, axs = plt.subplots(2, sharex= True)
fig.suptitle('Potencial do Rubídio (Z=37)')
axs[0].plot(r_a,v1,label='Z=1')
axs[0].legend(loc='lower right')
axs[0].plot(r_a,v2,label='Z=11')
axs[0].legend(loc='lower right')
axs[0].set(ylabel='V(r) [V]')
axs[1].plot(r_a,v,label='Z=1, r>=4a0; Z=11, r<4a0',color='g')
axs[1].legend(loc='lower right')
axs[1].set(xlabel='r/a0')
axs[1].set(ylabel='V(r) [V]')
plt.show()



"""
6.3, 6.4 e 6.5
"""
print('"""6.3, 6.4 e 6.5"""\n')
Eh=-13.60569301 #[eV]
En=[4.18,2.62,0.99] #[eV], https://nvlpubs.nist.gov/nistpubs/jres/048/jresv48n1p61_A1b.pdf
    
j=0
for i in En:
    n=5
    orbital=['5s','5p','5d']
    print(f"Z_eff para {orbital[j]} = {zef(n,i)}")
    j+=1

print('\n')
print('     A transposição de informação da alínea 3 para este comentário tamb\
ém é válida, pelas razões referidas anteriormente, pelo que o estado mais afet\
tado é o 5s e o menos afetado o 5d.\n   Utilizando os valores para as energias\
de ligação do eletrão nos estados 5s, 5p e 5d do Rb apresentados em Moore e R\
ussel (1952)- 4.18 eV para estado 5s, 2.62 eV para o estado 5p e 0.99 eV para \
o estado 5d- as conclusões são em tudo idênticas às apresentadas na alínea 4: \
a energia de ligação diminui com o nº quântico orbital.\n   Do mesmo modo os p\
otênciais efetivos calculados estão de acordo com as observações feitas acima,\
 i.e, quanto menor é a energia de ligação do eletrão, menor será o potencial e\
fetivo que interage com esse eletrão.')