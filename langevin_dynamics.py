import math
import sys
from math import exp
import numpy as np

## To switch from PES-1 to PES-3 comment the lines under 'PES-1' and uncomment
## the lines under 'PES-3'. Only the potential and dV_dx/dV_dy functions are altered.

def potential(pes_type,x,y,well=4):
    # For PES-1
    if pes_type == 1:
        return 0.02*(x**4+y**4) - well*exp(-((x+2)**2 + (y+2)**2)) - well*exp(-((x-2)**2 + (y-2)**2)) + 0.3*(x-y)**2 + 2*well*exp(-8)
    # For PES-2
    elif pes_type == 2:
        return 0.03*(x**4+y**4) - well*exp(-((x+2)**2 + (y+2)**2)) - well*exp(-((x-2)**2 + (y-2)**2)) + 0.4*(x-y)**2 + 4*exp(-(x**2+y**2)) - 2.1245
    # For invalid PES
    else:
        sys.exit("{} is an invalid choice of PES".format(pes_type))
    # For PES-3
    #return 0.02*(x**4+y**4) - 3.73*exp(-((x+2)**2/8 + (y+2)**2/8)) - 3.73*exp(-((x-2)**2/8 + (y-2)**2/8)) + 3*exp(-((x)**2/2 + (y)**2/15)) + 2*exp(-((x)**2/2 + (y)**2/2)) - 0.5085

def kinetic_energy(phasepoint):
    px = phasepoint[2]
    py = phasepoint[3]
    return 0.5*(px*px+py*py) 
 
def force(x,y,px,py,dt,beta,gamma,pes_type,well=4):
    std_dev = math.sqrt(2.0*gamma/(beta*dt))
    # For PES-1
    if pes_type == 1:
        dV_dx = 0.08*x**3 + 2*well*(x-2)*exp(-(x-2)**2-(y-2)**2) + 2*well*(x+2)*exp(-(x+2)**2-(y+2)**2) + 0.6*(x-y)
        dV_dy = -0.6*(x-y) + 2*well*(y-2)*exp(-(x-2)**2-(y-2)**2) + 2*well*(y+2)*exp(-(x+2)**2-(y+2)**2) + 0.08*y**3
    # For PES-2
    elif pes_type == 2:
        dV_dx = 0.12*x**3 + 2*well*(x-2)*exp(-(x-2)**2-(y-2)**2) + 2*well*(x+2)*exp(-(x+2)**2-(y+2)**2) + 0.8*(x-y) - 8*x*exp(-(x**2+y**2))
        dV_dy = 0.12*y**3 + 2*well*(y-2)*exp(-(x-2)**2-(y-2)**2) + 2*well*(y+2)*exp(-(x+2)**2-(y+2)**2) - 0.8*(x-y) - 8*y*exp(-(x**2+y**2))
    else:
        sys.exit("{} is an invalid choice of PES".format(pes_type))
    # For PES-3
    #dV_dx = 0.08*x**3 - 2*x*exp(-x**2/2-y**2/2) - 3*x*exp(-x**2/2-y**2/15) - 0.9325*(-x-2)*exp(-(1/8)*(x+2)**2-(1/8)*(y+2)**2) - 0.9325*(2-x)*exp(-(1/8)*(x-2)**2 - (1/8)*(y-2)**2)
    #dV_dy = -2*y*exp(-x**2/2 - y**2/2) - (2/5)*y*exp(-x**2/2 - y**2/15) - 0.9325*(-y-2)*exp(-(1/8)*(x+2)**2 - (1/8)*(y+2)**2) - 0.9325*(2-y)*exp(-(1/8)*(x-2)**2 - (1/8)*(y-2)**2) + 0.08*y**3
    fx = -dV_dx - gamma*px + np.random.normal(0,std_dev)
    fy = -dV_dy - gamma*py + np.random.normal(0,std_dev)
    return fx,fy
 
def vv_step(phasepoint,dt,beta,gamma,pes_type,well):
    x = phasepoint[0]
    y = phasepoint[1]
    px = phasepoint[2]
    py = phasepoint[3]
    fx = phasepoint[4]
    fy = phasepoint[5]
    px = px + (1/2)*dt*fx
    py = py + (1/2)*dt*fy
    x = x + dt*px
    y = y + dt*py
    fx,fy = force(x,y,px,py,dt,beta,gamma,pes_type,well=well)
    px = px + (1/2)*dt*fx
    py = py + (1/2)*dt*fy
    return np.asarray([x,y,px,py,fx,fy])
 
def calc_op(op_xtype,op_ytype,op_xcoef,op_ycoef,x,y):
    if len(op_xcoef) != op_xtype:
        sys.exit("Invalid op_xtype and op_xcoef pairing")
    if len(op_ycoef) != op_ytype:
        sys.exit("Invalid op_ytype and op_ycoef pairing")
    x_vec = [x**i for i in range(1,op_xtype+1)]
    y_vec = [y**i for i in range(1,op_ytype+1)]
#     if op_type == 1:
#         return x
#     elif op_type == 2:
#         return y
#     elif op_type == 3:
#         return x + y
#     else:
#         sys.exit("Invalid choice of OP")
    return np.dot(op_xcoef,x_vec) + np.dot(op_ycoef,y_vec)
    

