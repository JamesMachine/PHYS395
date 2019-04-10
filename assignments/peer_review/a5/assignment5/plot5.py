#!/usr/bin/env python

from numpy import loadtxt
import matplotlib.pyplot as plt

x_pt, psi_x =loadtxt('data_q5_1.dat',unpack=True)
# x_pt_neg, x1_neg,psi_x_neg,psi_prime_x_neg =loadtxt('data_q2_2.dat',unpack=True)

#plotting negative data file
plt.plot(x_pt, psi_x,color='k', linewidth=1.0) 
#plotting positive data file
# plt.plot(x_pt_neg, psi_x_neg,color='k', linewidth=1.0) 
plt.title('Problem 5 plot of different_energy level of Schrodinger equation')
plt.xlabel('x range')
plt.ylabel('\psi(x)')
plt.xlim([-10,10])
# plt.show()
plt.savefig('Problem5.pdf')