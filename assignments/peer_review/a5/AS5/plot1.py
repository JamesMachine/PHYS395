#!/usr/bin/env python

from numpy import loadtxt
import matplotlib.pyplot as plt

x_pt, x1,psi_x,psi_prime_x =loadtxt('data_q1_1.dat',unpack=True)
x_pt_neg, x1_neg,psi_x_neg,psi_prime_x_neg =loadtxt('data_q1_2.dat',unpack=True)
x_pt_1, x1,psi_x_1,psi_prime_x_1 =loadtxt('data_q1_3.dat',unpack=True)
x_pt_neg_2, x1_neg_2,psi_x_neg_2,psi_prime_x_neg_2 =loadtxt('data_q1_4.dat',unpack=True)

#plotting negative data file
plt.figure(0)
plt.plot(x_pt, psi_x,color='k', linewidth=1.0) 
#plotting positive data file
plt.plot(x_pt_neg, psi_x_neg,color='k', linewidth=1.0) 
plt.title('Problem 1 plot of first exicited state of Schrodinger equation')
plt.xlabel('x range')
plt.ylabel('\psi(x)')
# plt.show()
plt.savefig('Problem1_odd.pdf')

plt.figure(1)
plt.plot(x_pt_1, psi_x_1,color='k', linewidth=1.0) 
#plotting positive data file
plt.plot(x_pt_neg_2, psi_x_neg_2,color='k', linewidth=1.0) 
plt.title('Problem 1 plot of ground state of Schrodinger equation')
plt.xlabel('x range')
plt.ylabel('\psi(x)')
# plt.show()
plt.savefig('Problem1_even.pdf')