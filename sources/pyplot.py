import matplotlib.pyplot as plt
import numpy as np
with open('out1.dat', 'r') as f:
    lines = f.readlines()
    y = [float(line.split()[0]) for line in lines]
    x1 = [float(line.split()[1]) for line in lines]
    x2 = [float(line.split()[2]) for line in lines]
#    x3 = [float(line.split()[3]) for line in lines]

fig, panel = plt.subplots(1, 1) 
ax = plt.gca()
ax.invert_yaxis()   
panel.plot(x1, y, color='red', label = 'linear') # this is the  "linear" problem
panel.plot(x2, y, color='green', label = 'Jaupart') # this is the "Jaupart" non-linear problem etc. etc.
#panel.plot(x3, y, color='blue', label = 'Chapman') # this is the "Chapman" non-linear problem etc. etc.
panel.set_ylabel('Depth(m)')
panel.set_xlabel('Temp. (°C)')
panel.legend(loc='upper right',fontsize='medium')
plt.show()


with open('out2.dat', 'r') as f:
    lines = f.readlines()
    y = [float(line.split()[0]) for line in lines]
    x1 = [float(line.split()[1]) for line in lines]
    x2 = [float(line.split()[2]) for line in lines]
#    x3 = [float(line.split()[3]) for line in lines]

fig, panel = plt.subplots(1, 1) 
ax = plt.gca()
ax.invert_yaxis()   
panel.plot(x1, y, color='red', label = 'linear') # this is the  "linear" problem
panel.plot(x2, y, color='green', label = 'Jaupart') # this is the "Jaupart" non-linear problem etc. etc.
#panel.plot(x3, y, color='blue', label = 'Chapman') # this is the "Chapman" non-linear problem etc. etc.
panel.set_ylabel('Depth(m)')
panel.set_xlabel('Heat flow W/m2')
panel.legend(loc='upper right',fontsize='medium')
plt.show()

with open('out3.dat', 'r') as f:
    lines = f.readlines()
    y = [float(line.split()[0]) for line in lines]
    x1 = [float(line.split()[1]) for line in lines]
    x2 = [float(line.split()[2]) for line in lines]
#    x3 = [float(line.split()[3]) for line in lines]

fig, panel = plt.subplots(1, 1) 
ax = plt.gca()
ax.invert_yaxis()   
panel.plot(x1, y, color='red', label = 'linear') # this is the  "linear" problem
panel.plot(x2, y, color='green', label = 'Jaupart') # this is the "Jaupart" non-linear problem etc. etc.
#panel.plot(x3, y, color='blue', label = 'Chapman') # this is the "Chapman" non-linear problem etc. etc.
panel.set_ylabel('Depth(m)')
panel.set_xlabel('Thermal Conductivity W/(mK)')
panel.legend(loc='upper right',fontsize='medium')
plt.show()
