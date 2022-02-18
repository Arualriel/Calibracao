#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 19:29:30 2020

@author: laura
"""


# bibliotecas


import numpy as np
import matplotlib.pyplot as plt


global g,l,x0,y0,b0,b1,b



g,l,m,b=9.81,0.552,0.060314,0.5460896761464202

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



with open('data3.dat','r') as arquivo:
    linhas1 = arquivo.read().splitlines()
    
tamanho1=len(linhas1)
tempo1=np.zeros(tamanho1, int)
dados1=np.zeros(tamanho1)
erro_dados1=np.zeros(tamanho1)



j=0
for i in linhas1:
    tempo1[j]=int(float(i[0:7]))
    dados1[j]=float(i[10:14])
    j=j+1
         
err1=0.0

temponovo1=np.zeros(tamanho1,int)
dadosnovo1=np.zeros(tamanho1)

err1=0.0
for j in range(tamanho1-10,tamanho1):
    err1=err1+np.abs(dados1[j]-dados1[j-1])
    
err1=(err1/(20.0))*(((1.0/9.0)*np.pi)/1.14)



i=0
indice1=0

for j in range(len(tempo1)):
    tempo1[j]=tempo1[j]-5027665
    if(tempo1[j]>=0):
        temponovo1[i]=tempo1[j]
        dadosnovo1[i]=(dados1[j]-1.31-err1)*(((1.0/9.0)*np.pi)/1.14)
        if(dadosnovo1[i]==dadosnovo1[i-1])and(dadosnovo1[i]==dadosnovo1[i-2])and(indice1==0):
            
            indice1=i
        i=i+1

for k in range(i):
    if(k>=indice1):
        dadosnovo1[k]=0.0

td1=np.zeros(indice1,int)
dd1=np.zeros(indice1)
ed1=np.zeros(indice1)
for i in range(len(td1)):
    td1[i]=temponovo1[i]
    dd1[i]=dadosnovo1[i]




x0,y0=dd1[0],0.0 # mudar para m2

db=(b-b0)/(alfacuts-1.0)
alfa=np.zeros(2*alfacuts-1)
alfa[0],alfa[2*alfacuts-2]=b0,b1
pert=np.zeros(2*alfacuts-1)

for i in range(1,alfacuts):
    if(i==alfacuts-1):
        alfa[i]=b
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

tam=18000#8304 #te1
tf=6643 #7096 #t1 
X0=np.zeros((2,1))
X1=np.zeros((2,1))

dt=10.0**(-3.0)
R=np.zeros((len(alfa)+1,tam+2,2))

    
t=np.zeros(tam)
tb=np.zeros(len(alfa))
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

    if(alfa[i]==b):
        nf=R[i,tf,0]
        tb[i]=tf
        memb[i]=pert[i]
    
    xr=R[i,0:tam,0]
    linha1,=plt.plot(t,xr,label="Solução numérica para cada bi",color=listacores[i],ls='-')

dadosm1=plt.scatter(td1,dd1, marker='o',color='r',label="Dados para a massa m1")

plt.xlabel("t")
plt.ylabel("teta(t)")
plt.legend(handles=[dadosm1],loc='best')
plt.show() 


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



Fuzzyt=plt.scatter(tb,memb,marker='*',color='red',label="Pertinência(t1 (bi))")
plt.xlabel("t(bi)")
plt.ylabel("Pertinência de t(bi)")
plt.legend(handles=[Fuzzyt],loc='best')
plt.show()

print(tb)
print('epsilon1=',np.abs(tb[alfacuts-1]-tf))



