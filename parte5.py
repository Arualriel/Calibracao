#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 15:12:58 2020

@author: laura
"""



# bibliotecas


import numpy as np
import matplotlib.pyplot as plt


global g,l,x0,y0,b0,b1,b



g,l,m,b=9.81,0.552,0.0203469,0.5460896761464202

b1,b0=(b/10.0)+b,b-(b/10.0)

x0,y0=(-0.06147401609012629-0.05224119401595079)/2.0,0.0 # mudar para m2

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
tf=6600 #7096 #t1 
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


plt.xlabel("t")
plt.ylabel("teta(t)")
    
plt.show() 

Epsilon=5700.0
epsilon=10.0**(-4.0)#18

for i in range(len(alfa)):
    cont=0
    for j in range(tam):
        if((((np.abs(R[i,j,0]-0.0)+np.abs(R[i,j+1,0])+np.abs(R[i,j+2,0])+
              np.abs(R[i,j-1,0])+np.abs(R[i,j-2,0]))<epsilon)and(cont==0))and(j>tf+5000)):
            cont=1
            tb[i]=j
            memb[i]=pert[i]
            print(b,alfa[i],tb[i],i,j,np.abs(R[i,j,0]))
            

tb.sort()
print(tb,memb,nf)


Fuzzyt=plt.scatter(tb,memb,marker='*',color='navy',label="Pertinência(t3 (bi))")
plt.xlabel("t(bi)")
plt.ylabel("Pertinência de t(bi)")
plt.legend(handles=[Fuzzyt],loc='best')
plt.show()


print('estimativa para t3=',np.abs(tb[alfacuts-1]-Epsilon))

