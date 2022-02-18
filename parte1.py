#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 12:31:29 2020

@author: laura
"""



#bibliotecas

import numpy as np

global g,l,p,x0,y0,eta,r,mi,n, b

g,l,m,eta,r=9.81,0.8,0.125,0.01820,0.025#########corrigir!!!!########
p=6.0*np.pi*eta*r
y0=0.0
x0=np.pi/9.0
mi,n=0.2,0.005
b=-(n*mi)/(m*l)

#funcao f

def F(X):
    x=X[0,0]
    y=X[1,0]
    # B=0.5#0.25436054126242197
    f1 = y
    f2 =(-g/l)*x-(p/m+b)*y
    saida=np.matrix([[f1],[f2]])
    return saida


#solucao
    
def G(t,x0,y0):
    w=((-(p/m+b)**2.0+4.0*(g/l))**(0.5))/2.0
    c=x0*np.cos(w*t)
    s=((y0/w)+(x0*(p/m+b)/(2.0*w))*np.sin(w*t))
    solucao=(np.exp(-((p/m+b)/(2.0))*t))*(c+s)
    return solucao
    

#condicao inicial


X0=np.zeros((2,1))
X0[0,0],X0[1,0]=x0,y0


#variaveis auxiliares


num=50000 #numero de iteracoes

X1=np.zeros((2,1)) #vetor do passo seguinte

#taxa de variacao do tempo
dt=0.001
R=np.zeros((2,num)) #matriz de resultados
L=np.zeros((1,num))
erro=np.zeros(num)


R[:,0]=X0[:,0]



#aplicando o metodo

for k in range(num):
    K1=F(X0)
    K2=F(X0+dt*K1*0.5)
    K3=F(X0+dt*K2*0.5)
    K4=F(X0+dt*K3)
    X1=X0+(dt/6.0)*(K1+2*K2+2*K3+K4)
    
    R[0,k]=X1[0,0]
    R[1,k]=X1[1,0]
    X0=X1

import matplotlib.pyplot as plt

xr=R[0,:]
yr=R[1,:]
t=np.zeros(num)
for i in range(num-1):
    t[i]=dt*i
    L[:,i]=G(t[i],x0,y0)
   
tl=L[0,:]

for i in range(num):
    erro[i]=(R[0,i]-L[0,i])**2.0
    


   
linha1,=plt.plot(t[0:num-1],xr[0:num-1],label="solução numérica",color='blue',ls='-.')
linha2,=plt.plot(t[0:num-1],tl[0:num-1],label="solução exata",color='red',ls='-')
# linha3,=plt.plot(t[0:num-1],erro[0:num-1],label="erro",color='green')
plt.xlabel("t")
plt.ylabel("teta(t)")
plt.legend(handles=[linha1,linha2],loc='best')
plt.show() 

linha1,=plt.plot(t[0:num-1],xr[0:num-1],label="solução numérica",color='blue',ls='-.')

plt.xlabel("t")
plt.ylabel("teta(t)")
plt.legend(handles=[linha1],loc='best')
plt.show() 

linha2,=plt.plot(t[0:num-1],tl[0:num-1],label="solução exata",color='red',ls='-')

plt.xlabel("t")
plt.ylabel("teta(t)")
plt.legend(handles=[linha2],loc='best')
plt.show() 

linha3,=plt.plot(t[0:num-1],erro[0:num-1],label="erro",color='green')

plt.xlabel("t")
plt.ylabel("erro(t)")
plt.legend(handles=[linha3],loc='best')
plt.show() 
   