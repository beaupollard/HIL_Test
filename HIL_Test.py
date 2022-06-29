from mujoco_py import load_model_from_path, MjSim, MjViewer
import math
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mujoco_py.generated import const
from scipy.signal import butter, lfilter, freqz


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

## Get the ground truth knowledge ##
_xml_path = 'HIL.xml'
model  = load_model_from_path(_xml_path)
sim = MjSim(model)
# viewer = MjViewer(sim)
sim.step()
k2=sim.model.jnt_stiffness[-1]
b2=sim.model.dof_damping[-1]
m2=sim.model.body_mass[-1]
x0=sim.data.body_xpos[-1][0]
force = []
qpeg = []
vpeg=[]
qpeg_f=[]
while sim.data.time<15.:
    sim.data.ctrl[0]=3.
    sim.step()
    force.append(sim.data.sensordata[0])
    qpeg.append(sim.data.body_xpos[-1][0])
    vpeg.append(sim.data.body_xvelp[-1][0])


## Run the impedance matching code ##
_xml_path = 'HIL_control.xml'
model  = load_model_from_path(_xml_path)
sim = MjSim(model)
sim.step()
forcec = []
qpegc = []
vpegc=[]
dFrec=[]
k1=sim.model.jnt_stiffness[-1]
b1=sim.model.dof_damping[-1]
m1=sim.model.body_mass[-1]
x0=sim.data.body_xpos[-1][0]
dF=0
count=0
kp=26.
kd=50.
dt=sim.model.opt.timestep
ctrl=[]

x2_t1=0
x2_t2=0
while sim.data.time<15.:
    sim.data.ctrl[0]=3.
    x2=sim.data.body_xpos[-1][0]
    dx2dt=sim.data.body_xvelp[-1][0]
    
    ## Match impedances if there is a force ##
    if sim.data.sensordata[0]!=0:
        F2=-(k2*(x2-x0)+b2*dx2dt+m2*sim.data.cacc[-1][3])                       # Expected force given the space impedance model
        dF=sim.data.sensordata[0]-F2                                            # Difference in measured force and expected
        xdes=1/(m1/dt**2+b1/dt+k1)*(dF+m1/dt**2*(2*x2_t1-x2_t2)+b1/dt*x2_t1)    # Determine the delta x
        x2_t2=x2_t1                                                             # delta x[t-2]
        x2_t1=xdes                                                              # delta x[t-1]       
        dxdtdes=(xdes)/dt                                                       # Determine the velocity required to move delta x in delta t
        sim.data.ctrl[1]=-(kp*(x2-xdes)+dxdtdes*kd)                             # Position controller 
    
    ## Step the Simulation Forward ##
    sim.step()

    ## Record the data from the simulation ##
    forcec.append(sim.data.sensordata[0])
    qpegc.append(x2)
    vpegc.append(dx2dt)
    ctrl.append(sim.data.ctrl[1])


y = butter_lowpass_filter(np.array(forcec), 10, 1/dt, 6)
plt.plot(qpeg,'b')
plt.plot(qpegc,'r')
plt.show()

# while sim.data.time<15.:
#     # sim.data.ctrl[0]=3.
#     x2=sim.data.body_xpos[-1][0]
#     dx2=sim.data.body_xvelp[-1][0]
#     # if sim.data.sensordata[0]!=0:
#     #     F2=k2*(x2-x0)/m2+b2*dx2/m2
#     #     dF=sim.data.sensordata[0]-F2
#     #     # sim.data.ctrl[1]=dF
#     #     dx=1/(b1*dx2+k1*x2)*dF
#     dx=sim.data.body_xpos[-1][0]-qpeg[count]
#     dxdt=sim.data.body_xvelp[-1][0]-vpeg[count]
#     sim.data.ctrl[1]=-(kp*dx+dxdt*kd)
#     count=count+1
#     # print(sim.data.sensordata)
#     sim.step()
#     forcec.append(sim.data.sensordata[0])
#     qpegc.append(x2)
#     vpegc.append(dx2)
#     # dFrec.append(dF)
#     # viewer.render()
# plt.plot(qpeg,'b')
# plt.plot(qpegc,'r')
# plt.show()
# print('hey')