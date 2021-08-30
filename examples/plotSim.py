import numpy as np
import matplotlib.pyplot as plt
import sys
import time
from dataclasses import dataclass


testdir  = sys.argv[1]
filebase = sys.argv[2]
wallfile = sys.argv[3]
tsave    = int(sys.argv[4])
tfinal   = int(sys.argv[5])

R = 1.0

NT    = int(tfinal/tsave);
theta = np.linspace(0,2*np.pi,100);
xcirc = R*np.cos(theta);
ycirc = R*np.sin(theta);

fig = plt.figure(1);
plt.ion();
plt.show()
for i in range(NT):
    xy = np.genfromtxt(testdir + filebase + str(tsave*(i+1)) + '.csv', delimiter=',')
    for j in range(xy.shape[0]):
        plt.plot(xcirc + xy[j,0] , ycirc + xy[j,1] , 'b');
    plt.gca().set_aspect('equal');
    plt.xlim([-10,10])
    plt.ylim([-4,16])
    fig.canvas.draw()
    fig.canvas.flush_events()
    time.sleep(0.01)
    plt.clf()

plt.ioff(); plt.show();
