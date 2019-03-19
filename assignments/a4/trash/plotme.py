import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, ArtistAnimation
import matplotlib.animation as animation

#initialization
fig = plt.figure(figsize=(10,10))
ims = []

# read data
t,x0,y0,x1,y1,x2,y2,E = np.loadtxt('./output/results').T
rows = t.shape[0]

for i in range(rows):
  plt.title('motion animation')
  plt.xlabel("x")
  plt.ylabel("y")
  plt.xlim(-2,2)
  plt.ylim(-2,2)

  #px1, py1 = [0, 0], [x1[i], y1[i]]
  #px2, py2 = [x1[i], y1[i]], [x2[i], y2[i]]
  #plt.plot(0,0,'o')
  #im = plt.plot(px1, py1, px2, py2, marker = 'o')

  im = plt.plot(x1[i], y1[i], 'ro',x2[i],y2[i],'bo')
  ims.append(im)


ani = ArtistAnimation(fig, ims, interval=10)

try:
  plt.show()
except:
  print("Aborted!!!")

ani.save("./output/motion.gif", writer="imagemagick") #takes too much time



