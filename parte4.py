#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 30 01:48:50 2020

@author: laura
"""


# bibliotecas


import numpy as np
import matplotlib.pyplot as plt


global g,l,x0,y0,b0,b1,b



g,l,m,b=9.81,0.552,0.060314,0.5460896761464202#m1=0.060314  m2=0.0402140 m3=0.0203469

b1,b0=(b/10.0)+b,b-(b/10.0)


alfacuts=5

def F(X,B):
    x=X[0,0]
    y=X[1,0]
    f1 = y
    f2 =(-g/l)*x-B*y
    saida=np.matrix([[f1],[f2]])
    return saida


def pertinencia(x):
    if((x>=b0)and(x<=b)):
        return (x-b0)/(b-b0)
    elif ((x>b)and(x<=b1)):
        return (b1-x)/(b1-b)
    



    

with open('m2.txt','r') as arquivo:
    linhas = arquivo.read().splitlines()

tamanho=len(linhas)
tempo=np.zeros(tamanho, int)
dados=np.zeros(tamanho)
erro_dados=np.zeros(tamanho)



j=0
for i in linhas:
    tempo[j]=int(float(i[0:5]))
    dados[j]=float(i[8:12])
    j=j+1
         
err=0.0

(_, caps, _) = plt.errorbar(
    tempo, dados, xerr=erro_dados,fmt='o-',color='g', markersize=5, capsize=10)

for cap in caps:
    cap.set_markeredgewidth(err)
plt.xlabel("Tempo")
plt.ylabel("Volts")
dv="Dados da massa 2"
plt.legend([dv],loc='best')
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
    tempo[j]=tempo[j]-60491 #60491
    if(tempo[j]>=0):
        temponovo[i]=tempo[j]
        dadosnovo[i]=(dados[j]-1.31-err)*(((1.0/9.0)*np.pi)/1.14)
        if((dadosnovo[i]==dadosnovo[i-1])and(dadosnovo[i]==dadosnovo[i-2]))and((indice==0)):#and(dadosnovo[i]==dadosnovo[i+3])):
            indice=i
        i=i+1
indice=indice+25
for k in range(i):
    if(k>=indice):
        dadosnovo[k]=0.0

td=np.zeros(indice,int)
dd=np.zeros(indice)
ed=np.zeros(indice)
for i in range(len(td)):
    td[i]=temponovo[i]
    dd[i]=dadosnovo[i]






for j in range(len(td)):
    ed[j]=err

altura=1.0
(_, caps, _) = plt.errorbar(
    td, dd,xerr=ed,color='g', label="Dados experimentais",fmt='o', markersize=3, capsize=10)

for cap in caps:
    cap.set_markeredgewidth(altura)
plt.xlabel("Tempo")
plt.ylabel("Theta")
dl="Dados da massa 2"
plt.legend([dl],loc='best')
plt.show()
   
    
    


x0,y0=dd[0],0.0# m1= -0.05224119401595079,0.0 # mudar para m2

db=(b-b0)/(alfacuts-1.0)
alfa=np.zeros(2*alfacuts-1)
alfa[0],alfa[2*alfacuts-2]=b0,b1
pert=np.zeros(2*alfacuts-1)

for i in range(1,alfacuts):
    if(i==alfacuts-1):
        alfa[i]=b
        #pert[i]=pertinencia(alfa[i])
    else:
        alfa[i]=alfa[i-1]+db #uniform(b0,b)
        alfa[2*alfacuts-2-i]=alfa[2*alfacuts-1-i]-db
    pert[i]=pertinencia(alfa[i])
    
    pert[2*alfacuts-2-i]=pertinencia(alfa[2*alfacuts-2-i])
    

 
Fuzzy=plt.scatter(alfa,pert,marker='.',color='blue',label="Pertinencia(bi)")
plt.xlabel("bi")
plt.ylabel("Pertinencia de bi")
plt.legend(handles=[Fuzzy],loc='best')
plt.show()

tam=18000#  7392+1#7392 #te2 #m3=7848? m2=7393 m1 = 8304
tf=6638 #t2 6638 t1 6944
tf1=8303+1
t1=6643
X0=np.zeros((2,1))
X1=np.zeros((2,1))

dt=10.0**(-3.0)
R=np.zeros((len(alfa)+1,tam+2,2))

    
t=np.zeros(tam)
tb=np.zeros(len(alfa))
erro=np.zeros(len(alfa))
memb=np.zeros(len(alfa))
for i in range(tam):
    t[i]=i
    

listacores=["purple","magenta","red","orange","yellow","orange","red","magenta",
            "purple"]
    
for i in range(len(alfa)):
    X0[0,0],X0[1,0]=x0,y0
    X1=X0
    soma=0.0
    for k in range(tam):
        K1=F(X0,alfa[i])
        K2=F(X0+dt*K1*0.5,alfa[i])
        K3=F(X0+dt*K2*0.5,alfa[i])
        K4=F(X0+dt*K3,alfa[i])
        X1=X0+(dt/6.0)*(K1+2*K2+2*K3+K4)
        
        R[i,k,0]=X1[0,0]
        R[i,k,1]=X1[1,0]
        X0=X1
    for k in range(len(td)):
        soma=soma+(dd[k]-R[i,td[k],0])**2.0
        
    erro[i]=soma
    E=erro[i]
    if(alfa[i]==b):
        nf=R[i,tf,0]
        tb[i]=tf
        memb[i]=pert[i]
    
    xr=R[i,0:tam,0]
    linha1,=plt.plot(t,xr,label="Solução numérica para cada bi",color=listacores[i],ls='-')


# dadosm2=plt.scatter(td,dd, marker='o',color='g',label="Dados para a massa m2")

plt.xlabel("t")
plt.ylabel("teta(t)")
# plt.legend(handles=[dadosm2],loc='best')
    
plt.show() 

# ERRO=plt.scatter(alfa,erro,marker='x',color='orange',label="Erro")
# plt.xlabel("Bi")
# plt.ylabel("Erro relacionado ao Bi")
# plt.legend(handles=[ERRO],loc='best')
# plt.show()
menor=1.0
for i in range(len(erro)):
    if(menor>=erro[i]):
        menor=erro[i]
        
        bmin=alfa[i]


print('erro minimo m2',menor,bmin)
epsilon=10.0**(-4.0)#18
for i in range(len(alfa)):
    cont=0
    for j in range(tam):
        if((((np.abs(R[i,j,0]-0.0)+np.abs(R[i,j+1,0])+np.abs(R[i,j+2,0])+
              np.abs(R[i,j-1,0])+np.abs(R[i,j-2,0]))<epsilon)and(cont==0))and(j>tf+5000)):
            cont=1
            tb[i]=j
            memb[i]=pert[i]
tb.sort()
Fuzzyt2=plt.scatter(tb,memb,marker='*',color='green',label="Pertinência(t2 (bi))")
plt.xlabel("t(bi)")
plt.ylabel("Pertinência de t(bi)")
plt.legend(handles=[Fuzzyt2],loc='best')
plt.show()

print('epsilon2=',np.abs(tb[alfacuts-1]-tf))

print(tb)

