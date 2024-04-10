#STAB Mean-square and asymptotic stability test for E-M
#
# SDE is  dX = gamma*X dt + mu*X dW,  X(0) = Xzero,
#      where gamma and mu are constants and Xzero = 1
#
# Adapted from 
# Desmond J. Higham "An Algorithmic Introduction to Numerical Simulation of 
#                    Stochastic Differential Equations"
#
# http://www.caam.rice.edu/~cox/stoch/dhigham.pdf


import numpy as np
import matplotlib.pyplot as plt

T=20; M=50000; Xzero = 3
# b) Simulate (over the interval [0,20]) this stochastic process using an implicit method of the form
# 𝑋𝑛+1 = 𝑋𝑛 + (1 − 𝜃)𝛥𝑡𝑓(𝑋𝑛) + 𝜃𝛥𝑡𝑓(𝑋𝑛+1) + √𝛥𝑡𝛼𝑛𝑔(𝑋𝑛)
# For what values of 𝜇 𝑎𝑛𝑑 𝜎 is the SDE asymptotically stable.

ax = plt.subplot(211)
mu = 2; sigma = 0.1; theta = 0.5
for k in range(1,4):
    Dt=2**(1-k)
    N = int(T/Dt)
    Xms=np.zeros((N,1)); Xtemp=Xzero*np.ones((M,1))
    Xms[0]=Xzero
    for j in range(1,int(N)):
        alphan = np.random.randn(M,1)
        Xtemp = Xtemp + (1-theta)*mu*Xtemp*Dt + theta*mu*Xtemp*Dt + sigma*Xtemp*alphan*np.sqrt(Dt)
        Xms[j]=np.mean(Xtemp**2)
    ax.plot(np.arange(0,T,Dt), Xms)
    plt.yscale('log')

ax.legend(("$\Delta t = 1$", "$\Delta t = 1/2$", "$\Delta t = 1/4$"),loc=2)
plt.title("Mean square: $\mu = 2, \sigma = 0.1$",fontsize=12)
plt.ylabel("$E[X^2]$",fontsize=10), plt.axis([0,T,1e-20,1e+20])


bx = plt.subplot(212)
T=20
mu = 0.5; sigma = np.sqrt(6)
for k in range(1,4):
    Dt=2**(1-k)
    N = int(T/Dt)
    Xemabs=np.zeros((N,1)); Xtemp=Xzero
    Xemabs[0]=Xzero
    for j in range(1,int(N)):
        alphan = np.random.randn(1)
        Xtemp = Xtemp + (1-theta)*mu*Xtemp*Dt + theta*mu*Xtemp*Dt + sigma*Xtemp*alphan*np.sqrt(Dt)
        Xemabs[j]=np.abs(Xtemp)
    bx.plot(np.arange(0,T,Dt), Xemabs)
    plt.yscale('log')
bx.legend(("$\Delta t = 1$", "$\Delta t = 1/2$", "$\Delta t = 1/4$"),loc=2)
plt.title("Single path: $\mu = 0.5, \sigma = \sqrt{6}$",fontsize=12)
plt.ylabel("$|X|$",fontsize=10), plt.axis([0,T,1e-50,1e+100])

plt.suptitle("Mean-square and asymptotic stability test for E-M",fontsize=16)
# increase spacing between subplots
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

plt.savefig("imgs/4python.png")
plt.close()
plt.clf()
# For what values of 𝜃 is the implicit method mean-square stable.

# 𝑋𝑛+1 = 𝑋𝑛 + (1 − 𝜃)𝛥𝑡𝑓(𝑋𝑛) + 𝜃𝛥𝑡𝑓(𝑋𝑛+1) + √𝛥𝑡𝛼𝑛𝑔(𝑋𝑛)
def implicit_mean_square_stability(theta):
    T=20; M=50000; Xzero = 3
    mu = 2; sigma = 0.1
    Dt=1/4
    N = int(T/Dt)
    Xms=np.zeros((N,1)); Xtemp=Xzero*np.ones((M,1))
    Xms[0]=Xzero
    for j in range(1,int(N)):
        alphan = np.random.randn(M,1)
        Xtemp = Xtemp + (1-theta)*mu*Xtemp*Dt + theta*mu*Xtemp*Dt + sigma*Xtemp*alphan*np.sqrt(Dt)
        Xms[j]=np.mean(Xtemp**2)
    return max(Xms)

theta = np.linspace(0,1,100)
y = [implicit_mean_square_stability(t) for t in theta]
plt.plot(theta,y)
plt.title("Mean-square stability test for E-M")
plt.xlabel("$\\theta$")
plt.ylabel("$E[X^2]$")
plt.savefig("imgs/4python2.png")
plt.close()
plt.clf()