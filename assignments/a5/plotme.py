import numpy as np
import matplotlib.pyplot as plt


plt.figure(figsize=(10,10))
numpsi= 10

E = np.loadtxt('./output/output_q4/output_q4_E')
data = np.loadtxt('./output/output_q4/output_q4_psi') # transpose => first column is x, next each column is psi

x=data[:,0]
psi=data[:,1]
psisq=data[:,2]

for i in range(numpsi):
  start_i = 150*i
  end_i = 150*i+149

  plt.title("psi(x), psi(x)^2 for Energy level E"+str(i)+": " + str(E[i]))
  plt.xlabel("x")
  plt.ylabel("psi(x) or psi(x)^2")
  plt.xlim(-5,5)
  plt.ylim(-0.2,0.2)
  plt.grid()

  plt.plot(x[start_i:end_i],psi[start_i:end_i])
  plt.plot(x[start_i:end_i],psisq[start_i:end_i])

  plt.gca().legend(("psi(x)","psi(x)^2"))

  plt.savefig('./output/output_q4/output_q4_psi_E{}.png'.format(str(i)))

  plt.clf()
