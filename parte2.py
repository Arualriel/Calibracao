#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 12:31:29 2020

@author: laura
"""



# bibliotecas

from random import *
import numpy as np
import matplotlib.pyplot as plt
import time



inicio=time.time()
global g,l,x0,y0



g,l,m=9.81,0.552,0.060314




# funcao f

def F(X,B):
    x=X[0,0]
    y=X[1,0]
    f1 = y
    f2 =(-g/l)*x-B*y
    saida=np.matrix([[f1],[f2]])
    return saida


# solucao
    
def G(t,x0,y0,B):
    w=((-B**2.0+4.0*(g/l))**(0.5))/2.0
    c=x0*np.cos(w*t)
    s=((y0/w)+(x0*B/(2.0*w))*np.sin(w*t))
    solucao=(np.exp(-(B/(2.0))*t))*(c+s)
    return solucao
    


# arquivo de dados





with open('data3.dat','r') as arquivo:
    linhas = arquivo.read().splitlines()
    
tamanho=len(linhas)
tempo=np.zeros(tamanho, int)
dados=np.zeros(tamanho)
erro_dados=np.zeros(tamanho)
cont=0


j=0
for i in linhas:
    tempo[j]=int(float(i[0:7]))
    dados[j]=float(i[10:14])
    j=j+1
         
err=0.0

(_, caps, _) = plt.errorbar(
    tempo, dados, xerr=erro_dados,fmt='o-', markersize=5, capsize=10)

for cap in caps:
    cap.set_markeredgewidth(err)

plt.show()


temponovo=np.zeros(tamanho,int)
dadosnovo=np.zeros(tamanho)

err=0.0
for j in range(tamanho-10,tamanho):
    err=err+np.abs(dados[j]-dados[j-1])
    
err=(err/(20.0))*(((1.0/9.0)*np.pi)/1.14)



i=0
indice=0

for j in range(len(tempo)):
    tempo[j]=tempo[j]-5027515
    if(tempo[j]>=0):
        temponovo[i]=tempo[j]
        dadosnovo[i]=(dados[j]-1.31-err)*(((1.0/9.0)*np.pi)/1.14)
        if(dadosnovo[i]==dadosnovo[i-1])and(dadosnovo[i]==dadosnovo[i-2])and(indice==0):
            
            indice=i
        i=i+1

for k in range(i):
    if(k>=indice):
        dadosnovo[k]=0.0

td=np.zeros(indice,int)
dd=np.zeros(indice)
ed=np.zeros(indice)
for i in range(len(td)):
    td[i]=temponovo[i]
    dd[i]=dadosnovo[i]

print(td,dd)




for j in range(len(td)):
    ed[j]=err
altura=1.0

(_, caps, _) = plt.errorbar(
    td, dd,xerr=ed, label="Dados experimentais",fmt='o', markersize=3, capsize=10)

for cap in caps:
    cap.set_markeredgewidth(altura)

plt.show()

###################

num=len(td) #tempo total do experimento em segundos
N=100
tam=int(td[num-1])+1
X0=np.zeros((2,1))
y0=0.0
x0=dd[0]
X0[0,0],X0[1,0]=x0,y0

# variaveis auxiliares




X1=np.zeros((2,1)) #vetor do passo seguinte

#taxa de variacao do tempo
dt=10.0**(-3.0)
R=np.zeros((2,tam)) #matriz de resultados
L=np.zeros((1,tam))
erro=np.zeros(N)
B=np.zeros(N)

# intervalo em que B esta [i1,i2]

i1=0.4
i2=0.7
R[0,0],R[1,0]=X0[0,0],X0[1,0]


epsilon=10.0**(-5.0)
E=1.0
a=0

seed()
while(E>=epsilon)and(a<N):
    B[a]=uniform(i1,i2) #normalvariate((i1+i2)/2.0,0.02) 
    X0[0,0],X0[1,0]=x0,y0
    X1=X0
    soma=0.0
    for k in range(tam):
        K1=F(X0,B[a])
        K2=F(X0+dt*K1*0.5,B[a])
        K3=F(X0+dt*K2*0.5,B[a])
        K4=F(X0+dt*K3,B[a])
        X1=X0+(dt/6.0)*(K1+2*K2+2*K3+K4)
        
        R[0,k]=X1[0,0]
        R[1,k]=X1[1,0]
        X0=X1
    for k in range(len(td)):
        soma=soma+(dd[k]-R[0,td[k]])**2.0
        
        
    erro[a]=soma#((soma)**(0.5))/num
    E=erro[a]
    a=a+1
menor=erro[0]
ind=0
for k in range(a):
    if(erro[k]<=menor):
        menor=erro[k]
        ind=k

#### B best=B[ind]!
Bbest=B[ind]
print('B=',B,'erro=',erro)
#print('menor erro e indice',menor,ind)
#print('melhor b',Bbest)
tr=np.zeros(num)
X0[0,0],X0[1,0]=x0,y0
for k in range(tam):
        
    K1=F(X0,Bbest)
    K2=F(X0+dt*K1*0.5,Bbest)
    K3=F(X0+dt*K2*0.5,Bbest)
    K4=F(X0+dt*K3,Bbest)
    X1=X0+(dt/6.0)*(K1+2*K2+2*K3+K4)
        
    R[0,k]=X1[0,0]
    R[1,k]=X1[1,0]
    X0=X1
for k in range(num):
    tr[k]=np.abs(R[0,k]-dd[k])


xr=R[0,:]
yr=R[1,:]
t=np.zeros(tam)
for i in range(tam):
    t[i]=i
    L[0,i]=G(t[i]*10.0**(-3.0),x0,y0,Bbest)
    

tl=L[0,:]


   
linha1,=plt.plot(t,xr,label="Solução numérica",color='blue',ls='-.')

linha2,=plt.plot(t,tl,label="Solução exata",color='purple',ls='-')

data=plt.scatter(td,dd,marker='o',color='r',label="Dados")
plt.xlabel("t")
plt.ylabel("teta(t)")
plt.legend(handles=[linha1,linha2,data],loc='best')

plt.show() 


ERRO=plt.scatter(B,erro,marker='x',color='purple',label="Erro")
plt.xlabel("Bi")
plt.ylabel("Erro relacionado ao Bi")
plt.legend(handles=[ERRO],loc='best')
plt.show()


fim=time.time()
print('TEMPO DE EXECUCAO=',fim-inicio,'melhor B=',Bbest,'menor erro=',menor)




