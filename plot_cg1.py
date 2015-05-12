import matplotlib.pyplot as plt
import math
import scipy as sp
import numpy as np


data_a = np.loadtxt('efgh_C40_aniso.txt');

N = 30;
n1 = 550;
Nav = 100; 


ka = data_a[0:N,0];
EaY = data_a[n1*N:(n1+1)*N,1]; 
Ea0 = data_a[n1*N:(n1+1)*N,2];
EaX = data_a[n1*N:(n1+1)*N,3];  

FaY = data_a[n1*N:(n1+1)*N,4]; 
Fa0 = data_a[n1*N:(n1+1)*N,5];
FaX = data_a[n1*N:(n1+1)*N,6];  


GaY = data_a[n1*N:(n1+1)*N,7]; 
Ga0 = data_a[n1*N:(n1+1)*N,8];
GaX = data_a[n1*N:(n1+1)*N,9];  

HaY = data_a[n1*N:(n1+1)*N,10]; 
Ha0 = data_a[n1*N:(n1+1)*N,11];
HaX = data_a[n1*N:(n1+1)*N,12];  



for n in range(Nav): 
    EaY = EaY + data_a[(n1+n)*N:(n1+n+1)*N,1];
    Ea0 = Ea0 + data_a[(n1+n)*N:(n1+n+1)*N,2];
    EaX = EaX + data_a[(n1+n)*N:(n1+n+1)*N,3];
    FaY = FaY + data_a[(n1+n)*N:(n1+n+1)*N,4];
    Fa0 = Fa0 + data_a[(n1+n)*N:(n1+n+1)*N,5];
    FaX = FaX + data_a[(n1+n)*N:(n1+n+1)*N,6];
    GaY = GaY + abs(data_a[(n1+n)*N:(n1+n+1)*N,7]);
    Ga0 = Ga0 + abs(data_a[(n1+n)*N:(n1+n+1)*N,8]);
    GaX = GaX + abs(data_a[(n1+n)*N:(n1+n+1)*N,9]);
    HaY = HaY + abs(data_a[(n1+n)*N:(n1+n+1)*N,10]);
    Ha0 = Ha0 + abs(data_a[(n1+n)*N:(n1+n+1)*N,11]);
    HaX = HaX + abs(data_a[(n1+n)*N:(n1+n+1)*N,12]);


Ea = EaY + 2*Ea0 +EaX;
Fa = FaY + 2*Fa0 +FaX;
Ga = GaY + 2*Ga0 +GaX;
Ha = HaY + 2*Ha0 +HaX;

#--------------- PLOTTING


plt.subplot(2,1,1);
plt.loglog(ka,Fa,'b')
plt.loglog(ka[1:19],ka[1:19]**(-11./3.)*math.exp(6.5),'g')
plt.loglog(ka[19:27],ka[19:27]**(-4.8)*math.exp(9),'r')
plt.loglog(ka[27:60],ka[27:60]**(-3.2)*math.exp(1.5),'c')

 
plt.subplot(2,1,2);
plt.loglog(ka,Ea,'b')
plt.loglog(ka[1:20],ka[1:20]**(-5./3.)*math.exp(5.5),'g')
plt.loglog(ka[20:60],ka[20:60]**(-3.2)*math.exp(9),'r')
plt.show()


