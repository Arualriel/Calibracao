#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 12:31:29 2020

@author: laura
"""



#bibliotecas

import numpy as np

# global g,l,x0,y0, b


global g,l,p,x0,y0,eta,r,mi,n, b

g,l,m,eta,r=9.81,0.8,0.125,0.01820,0.025#########corrigir!!!!########
p=6.0*np.pi*eta*r
y0=0.0
x0=np.pi/6.0
mi,n=0.2,0.005
b=-(n*mi)/(m*l)

# g,l,m=9.81,0.552,0.060314#########corrigir!!!!########

# y0=0.0
# x0=np.pi/9.0

##b= calibrar

#funcao f

def F(X):
    x=X[0,0]
    y=X[1,0]
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
    


#arquivo de dados
import matplotlib.pyplot as plt


num=20000 #milissegundos 20*10^3

with open('data3.dat','r') as arquivo:
    linhas = arquivo.read().splitlines()
    
tamanho=len(linhas)
tempo=np.zeros(tamanho)
dados=np.zeros(tamanho)
erro_dados=np.zeros(tamanho)
cont=0


j=0
for i in linhas:
    tempo[j]=float(i[0:7])
    dados[j]=float(i[10:16])
    j=j+1
         




err=0.0

(_, caps, _) = plt.errorbar(
    tempo, dados, xerr=erro_dados,fmt='o-', markersize=5, capsize=10)

for cap in caps:
    cap.set_markeredgewidth(err)

plt.show()


temponovo=np.zeros(tamanho)
dadosnovo=np.zeros(tamanho)

err=0.0
for j in range(tamanho-10,tamanho):
    err=err+np.abs(dados[j]-dados[j-1])
    
err=(err/(20.0))*(((1.0/9.0)*np.pi)/1.14)



i=0
indice=0
for j in range(len(tempo)):
    tempo[j]=tempo[j]-5027665.0
    if(tempo[j]>=0):
        temponovo[i]=tempo[j]
        dadosnovo[i]=(dados[j]-1.31-err)*(((1.0/9.0)*np.pi)/1.14)
        if(dadosnovo[i]==dadosnovo[i-1])and(dadosnovo[i]==dadosnovo[i-2])and(indice==0):
            # print('entrou!')
            indice=i
        i=i+1

for k in range(i):
    if(k>=indice):
        dadosnovo[k]=0.0

td=np.zeros(indice)
dd=np.zeros(indice)
ed=np.zeros(indice)
for i in range(len(td)):
    td[i]=temponovo[i]
    dd[i]=dadosnovo[i]


# print('dados',dd)



for j in range(len(td)):
    ed[j]=err
altura=1.0

(_, caps, _) = plt.errorbar(
    td, dd,xerr=ed,fmt='o', markersize=3, capsize=10)

for cap in caps:
    cap.set_markeredgewidth(altura)

plt.show()




X0=np.zeros((2,1))
X0[0,0],X0[1,0]=x0,y0


#variaveis auxiliares




X1=np.zeros((2,1)) #vetor do passo seguinte

#taxa de variacao do tempo
dt=0.01
R=np.zeros((6,num)) #matriz de resultados
L=np.zeros((1,num))
erro=np.zeros(3)
graficoerro=np.zeros(num)
graficot=np.zeros(num)
R[0,0],R[1,0],R[2,0],R[3,0],R[4,0],R[5,0]=X0[0,0],X0[1,0],X0[0,0],X0[1,0],X0[0,0],X0[1,0]


t=np.zeros(num)
#aplicando o metodo
for j in range(2,4):
    dt=0.01 #10**(-j)
    X0[0,0],X1[1,0]=x0,y0
    X1=X0
    for k in range(num-1):
        
        K1=F(X0)
        K2=F(X0+dt*K1*0.5)
        K3=F(X0+dt*K2*0.5)
        K4=F(X0+dt*K3)
        X1=X0+(dt/6.0)*(K1+2*K2+2*K3+K4)
        
        R[j,k]=X1[0,0]
        R[j+1,k]=X1[1,0]
        X0=X1
        t[k]=k*dt
    xr=R[0,:]
    linha1,=plt.plot(t[0:num-1],xr[0:num-1],label="solução numérica",color='blue',ls='-.')

    plt.xlabel("t")
    plt.ylabel("teta(t)")
    plt.legend(handles=[linha1],loc='best')
    plt.show() 

yr=R[1,:]

for i in range(num-1):
    t[i]=i*10**(-3)
    L[0,i]=G(t[i],x0,y0)
    
   
tl=L[0,:]
emax=0.0
eind=0

for i in range(3):
    soma=0.0
    for j in range(num-1):
        soma=soma+(R[i,j]-L[0,j])**2.0
    erro[i]=((soma)**(0.5))/num
print(erro)
   
# linha1,=plt.plot(t[0:num-1],xr[0:num-1],label="solução numérica",color='blue',ls='-.')
# linha2,=plt.plot(t[0:num-1],tl[0:num-1],label="solução exata",color='red',ls='-')
# linha3,=plt.plot(t[0:num-1],erro[0:num-1],label="erro",color='green')
# plt.xlabel("t")
# plt.ylabel("teta(t)")
# plt.legend(handles=[linha1,linha2,linha3],loc='best')
# plt.show() 

# linha1,=plt.plot(t[0:num-1],xr[0:num-1],label="solução numérica",color='blue',ls='-.')

# plt.xlabel("t")
# plt.ylabel("teta(t)")
# plt.legend(handles=[linha1],loc='best')
# plt.show() 

# linha2,=plt.plot(t[0:num-1],tl[0:num-1],label="solução exata",color='red',ls='-')

# plt.xlabel("t")
# plt.ylabel("teta(t)")
# plt.legend(handles=[linha2],loc='best')
# # plt.show() 

# linha3,=plt.plot(t[0:num-1],erro[0:num-1],label="erro",color='green')

# plt.xlabel("t")
# plt.ylabel("erro(t)")
# plt.legend(handles=[linha3],loc='best')
# plt.show() 



