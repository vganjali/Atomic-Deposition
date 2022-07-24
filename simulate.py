###############################################
# Atomic Deposition Project, EE216, Fall 2017 #
# Developed by Yucheng Li, Md Nafiz Amin and  #
# Vahid Ganjalizadeh                          #
###############################################

import sys
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np
import threading as trd
import time


name = '>> EE216-Project Copper Atom Deposition'
resolution = 1
screen_size = [1200, 900]
grid_size = [40 * resolution, 40 * resolution, 50 * 1]

# Copper
radius = 2.*resolution
mass = 3.0758e5
#mass = 1.055e-25

# L-J Potential
sigma = 3.348
rm = 1.122 * sigma * resolution
eps = 2.1407e31
#eps = 3.589e-20

#gravity = 9.89
gravity = 1.413e11
mu_k = 0.0

coarse_step = 1.0e-7
fine_step = 2.0e-15
#coarse_step = 1.0e-2
#fine_step = 1.0e-4

time_step = .001e1
speed = 2000.
timer = 0.0
run_time = 0.0
prev_time = 0.0

particles_count = 1000
particle_active = False
current_particle = np.zeros(9,dtype=float)
current_position = np.zeros(3,dtype=float)
#critical_distance = 2*np.round(rm).astype(int)
critical_distance = 7
minimum_distance = 0.
tolerance_z = 1.5*mass*gravity
tolerance = 1.0e40
tolerance_v = 1.0e12
position_change = 1.0e-5
tolerance_v_d = 2e1
damp_factor = np.ones(3,dtype=float)
force_damp_factor = 0.02
velocity_damp_factor = 0.
osc_v = np.zeros(3,dtype=float)
timeout = 0
timeout_init = 20
stuck_timeout = 0
stuck_timeout_init = 20
eq_point = np.zeros(3,dtype=float)
initial_velocity = -1.0e5

global wd
global ht
global mouseX
global mouseY
global aff
global nrange

global s2
global p2
global axrng
axrng = 5.0
p2 = 0
s2 = 0

wd = 600
ht = 600
mouseX = wd/2
mouseY = ht/2
brake = 128.0
aff = [1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0]

colors = [[1.,1.,1.,1.],[0.,.1,.5,1.0],[0.0,0.5,0.0,1.0],[0.5,0.5,0.0,1.0],[0.5,0.0,0.5,1.0],[0.0,0.5,0.5,1.0]]
#colors = [[1.,1.,1.,1.],[0.,.1,.5,1.0]]
add_del = True  # add/delete new point if True/False
particles = np.zeros((particles_count,4),dtype=int)
p_lattice = np.zeros((grid_size[0]+1,grid_size[1]+1,grid_size[2]+1,5),dtype=int)
p_extended = np.zeros((grid_size[0]+1+2*critical_distance,grid_size[1]+1+2*critical_distance,grid_size[2]+1+2*critical_distance),dtype=int)
p_sub = np.zeros((2*critical_distance+1,2*critical_distance+1,2*critical_distance+1),dtype=int)
fx = np.zeros((2*critical_distance+1,2*critical_distance+1,2*critical_distance+1),dtype=float)
fy = np.zeros((2*critical_distance+1,2*critical_distance+1,2*critical_distance+1),dtype=float)
fz = np.zeros((2*critical_distance+1,2*critical_distance+1,2*critical_distance+1),dtype=float)
f_lut = np.zeros((2*critical_distance+1,2*critical_distance+1,2*critical_distance+1),dtype=float)
r_lut = np.zeros((2*critical_distance+1,2*critical_distance+1,2*critical_distance+1),dtype=float)
vlj_lut = np.zeros((2*critical_distance+1,2*critical_distance+1,2*critical_distance+1),dtype=float)
theta_lut = np.zeros((2*critical_distance+1,2*critical_distance+1,2*critical_distance+1),dtype=float)
phi_lut = np.zeros((2*critical_distance+1,2*critical_distance+1,2*critical_distance+1),dtype=float)
#print(particles.shape)
count = 0
initial_count = 0
shooted_count = 0



grid_box_verticies = (
    (0, 0, 0),
    (grid_size[0], 0, 0),
    (grid_size[0], grid_size[1], 0),
    (0, grid_size[1], 0),
    (0, 0, grid_size[2]),
    (grid_size[0], 0, grid_size[2]),
    (grid_size[0], grid_size[1], grid_size[2]),
    (0, grid_size[1], grid_size[2])
    )

grid_box_lines = (
    (0,1),
    (0,3),
    (0,4),
    (2,1),
    (2,3),
    (2,6),
    (5,1),
    (5,4),
    (5,6),
    (7,3),
    (7,4),
    (7,6)
    )
def force():
    global fx,fy,fz,r_lut,vlj_lut,f_lut
    for x in range(-critical_distance, critical_distance + 1, 1):
        for y in range(-critical_distance, critical_distance + 1, 1):
            for z in range(-critical_distance, critical_distance + 1, 1):
                r = np.sqrt(float(x**2 + y**2 + z**2))
                if ( (r > 0) & (r <= critical_distance)):
                    r_lut[x+critical_distance,y+critical_distance,z+critical_distance] = r
                    vlj_lut[x+critical_distance,y+critical_distance,z+critical_distance] = eps*((rm/r)**(12) - 2*(rm/r)**(6))
                    f_lut[x+critical_distance,y+critical_distance,z+critical_distance] = eps*(12*(r/rm)**(-7)-12*(r/rm)**(-13))/rm
                    if (z > 0):
                        if (x > 0):
                            phi_lut[x+critical_distance,y+critical_distance,z+critical_distance] = np.arctan(float(y)/float(x))
                            theta_lut[x+critical_distance,y+critical_distance,z+critical_distance] = np.arctan(np.sqrt(float(x**2+y**2))/float(z))
                        elif (x < 0):
                            phi_lut[x+critical_distance,y+critical_distance,z+critical_distance] = np.arctan(float(y)/float(x))+np.pi
                            theta_lut[x+critical_distance,y+critical_distance,z+critical_distance] = np.arctan(np.sqrt(float(x**2+y**2))/float(z))
                        else:
                            phi_lut[x+critical_distance,y+critical_distance,z+critical_distance] = np.sign(y)*np.pi/2
                            theta_lut[x+critical_distance,y+critical_distance,z+critical_distance] = np.arctan(float(np.abs(y))/float(z))
                    elif (z < 0):
                        if (x > 0):
                            phi_lut[x+critical_distance,y+critical_distance,z+critical_distance] = np.arctan(float(y)/float(x))
                            theta_lut[x+critical_distance,y+critical_distance,z+critical_distance] = np.arctan(np.sqrt(float(x**2+y**2))/float(z))+np.pi
                        elif (x < 0):
                            phi_lut[x+critical_distance,y+critical_distance,z+critical_distance] = np.arctan(float(y)/float(x))+np.pi
                            theta_lut[x+critical_distance,y+critical_distance,z+critical_distance] = np.arctan(np.sqrt(float(x**2+y**2))/float(z))+np.pi
                        else:
                            phi_lut[x+critical_distance,y+critical_distance,z+critical_distance] = np.sign(y)*np.pi/2
                            theta_lut[x+critical_distance,y+critical_distance,z+critical_distance] = np.arctan(float(np.abs(y))/float(z))+np.pi
                    else:
                        theta_lut[x+critical_distance,y+critical_distance,z+critical_distance] = np.pi/2
                        if (x > 0):
                            phi_lut[x+critical_distance,y+critical_distance,z+critical_distance] = np.arctan(float(y)/float(x))
                        elif (x < 0):
                            phi_lut[x+critical_distance,y+critical_distance,z+critical_distance] = np.arctan(float(y)/float(x))+np.pi
                        else:
                            phi_lut[x+critical_distance,y+critical_distance,z+critical_distance] = np.sign(y)*np.pi/2
                    #print(theta_lut[x+critical_distance,y+critical_distance])
                    fx[x+critical_distance,y+critical_distance,z+critical_distance] = \
                        f_lut[x+critical_distance,y+critical_distance,z+critical_distance]* \
                        np.sin(theta_lut[x+critical_distance,y+critical_distance,z+critical_distance])* \
                        np.cos(phi_lut[x+critical_distance,y+critical_distance,z+critical_distance])
                    fy[x+critical_distance,y+critical_distance,z+critical_distance] = \
                        f_lut[x+critical_distance,y+critical_distance,z+critical_distance]* \
                        np.sin(theta_lut[x+critical_distance,y+critical_distance,z+critical_distance])* \
                        np.sin(phi_lut[x+critical_distance,y+critical_distance,z+critical_distance])
                    fz[x+critical_distance,y+critical_distance,z+critical_distance] = \
                        f_lut[x+critical_distance,y+critical_distance,z+critical_distance]* \
                        np.cos(theta_lut[x+critical_distance,y+critical_distance,z+critical_distance])
                    if (x == 0):
                        fx[x+critical_distance,y+critical_distance,z+critical_distance] = 0.0
                    if (y == 0):
                        fy[x+critical_distance,y+critical_distance,z+critical_distance] = 0.0
                    if (z == 0):
                        fz[x+critical_distance,y+critical_distance,z+critical_distance] = 0.0
##    print(np.amax(fx))
##    print(np.amin(np.abs(fx[np.nonzero(fx)])))
##    print("theta: ")
##    print(theta_lut)
##    print("phi: ")
##    print(phi_lut)
##    print("fx: ")
##    print(fx)
##    print("fy: ")
##    print(fy)
##    print("fz: ")
##    print(fz[critical_distance,:,:])
##    #plt.plot(r_lut[critical_distance+4:,critical_distance,critical_distance],f_lut[critical_distance+4:,critical_distance,critical_distance])
##    plt.plot(vlj_lut[critical_distance,critical_distance,critical_distance+1:])
##    print(np.amin(vlj_lut[critical_distance,critical_distance,critical_distance+1:]))
##    plt.show()
def velocity(f):
    global current_particle
    return current_particle[3:6]+time_step/(2*mass)*(current_particle[6:]+f)

def init():
    glClearColor(0.,0.,0.,1.)
    glShadeModel(GL_SMOOTH)
    #glutFullScreen()
    
def main():
    global particles, count, particle_active, p_lattice, initial_count, vlj_lut, r_lut
    glutInit(sys.argv)
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
    glutInitWindowSize(screen_size[0], screen_size[1])
    glutCreateWindow(name)
    init()
    force()
    p_lattice[:,:,0,3] = np.ones((grid_size[0]+1,grid_size[1]+1))
    for x in range(grid_size[0]+1):
        for y in range(grid_size[1]+1):
            if (((float(x) % round(rm)) == 0) and ((float(y) % round(rm)) == 0)):
                particles[count] = np.array([x,y,0,0])
                count += 1
    initial_count = count
    glEnable(GL_CULL_FACE)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_LIGHTING)
    lightZeroPosition = [0.,0.,10*np.max(grid_size[:2]),10000.]
    lightZeroColor = [1.,1.,1.,1.]
    glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition)
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor)
    glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 0.1/grid_size[2])
    glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.1/grid_size[2])
    glEnable(GL_LIGHT0)
    glutDisplayFunc(display)
    glMatrixMode(GL_PROJECTION)
    #gluPerspective(6.,float(screen_size[0])/float(screen_size[1]),0.1,10000.)
    glOrtho(-1.1*grid_size[0], grid_size[0], -grid_size[1], 1*grid_size[1], -1.0, 100000.0)
    glMatrixMode(GL_MODELVIEW)
    gluLookAt(2.5*np.amax(grid_size),2.2*np.amax(grid_size),3*np.amax(grid_size),
              0.5*np.amax(grid_size[:2]),0.2*np.amax(grid_size[:2]),0,
              0,0,1)
    glPushMatrix()
    glutIdleFunc(display)
    #glutMotionFunc(motion)
    #glutPassiveMotionFunc(mousemotion)
    #glutReshapeFunc(reshape)
    glutMainLoop()

def collide(p2):
    global current_particle
    speed = np.sqrt(np.sum(current_particle[3:6]**2))
    diff = p2[:3] - current_particle[:3]
    angle1 = np.arctan(diff[1]/(np.sqrt((diff[0]**2)+(diff[2]**2))))
    if diff[0] != 0 and diff[2] != 0:
        if diff[0] > 0:
            if diff[2] > 0:   angle2 = np.arctan(diff[2]/diff[0])
            elif diff[2] < 0: angle2 = np.arctan(diff[2]/diff[0])
        elif diff[0] < 0:
            if diff[2] > 0:   angle2 =  np.pi + np.arctan(diff[2]/diff[0])
            elif diff[2] < 0: angle2 = -np.pi + np.arctan(diff[2]/diff[0])
        xspeed = -speed*np.cos(angle2)*np.cos(angle1)
        yspeed = -speed*np.sin(angle1)
        zspeed = -speed*np.sin(angle2)*np.cos(angle1)
    else:
        if diff[0] == 0 and diff[2] == 0:
            angle2 = 0
        if diff[0] == 0 and diff[2] != 0:
            if diff[2] > 0:   angle2 = -np.pi/2
            else:           angle2 =  np.pi/2
        if diff[0] != 0 and diff[2] == 0:
            if diff[0] < 0:   angle2 =  0.0
            else:           angle2 =  np.pi
        xspeed = speed*np.cos(angle2)*np.cos(angle1)
        yspeed = speed*np.sin(angle1)
        zspeed = speed*np.sin(angle2)*np.cos(angle1)
    current_particle[3:6] = np.array([xspeed,yspeed,zspeed])*damp_factor

def add_particle(x,y,z):
    global particles, count, particle_active, p_lattice, p_extended, p_sub, current_particle, timeout, stuck_timeout, damp_factor
    if (p_lattice[int(round(x)),int(round(y)),int(round(z)),3] == 0):
        p_extended[critical_distance:-critical_distance,critical_distance:-critical_distance,critical_distance:-critical_distance] = p_lattice[:,:,:,3]
        current_particle = np.array([x,y,z,0.,0.,initial_velocity,0.,0.,-mass*gravity])
        particle_active = True
        count += 1
        timeout = timeout_init
        stuck_timeout = stuck_timeout_init
        damp_factor = np.ones(3,dtype=float)
        #print("Particle added")
    #print count

def move_particle():
    global particles, count, particle_active, p_sub, p_lattice, current_particle, r_lut, timeout, shooted_count, osc_v, time_step, eq_point, stuck_timeout, current_position, run_time, damp_factor
    #print(current_particle[:3])
##    p_sub = p_extended[int(round(current_particle[0])):int(round(current_particle[0]))+2*critical_distance+1, \
##                       int(round(current_particle[1])):int(round(current_particle[1]))+2*critical_distance+1, \
##                       int(round(current_particle[2])):int(round(current_particle[2]))+2*critical_distance+1]
    #print(p_extended.shape)
    #print(p_sub.shape)
    near_particles = np.multiply(r_lut,p_sub.astype(float))
    #print(near_particles)
    nearest_particle = 0
    if (np.count_nonzero(near_particles) == 0):
        time_step = coarse_step
    #print(p_sub)
    #print(nearest_particle)
        current_particle[:3] = current_particle[:3]+time_step*current_particle[3:6]+0.5*((time_step**2)/mass)*current_particle[6:]
        if (current_particle[2] <= radius):
            current_particle[2] = radius
            current_particle[5] = 0.
            current_particle[8] = 0.
        if ((current_particle[2] > grid_size[2]) or \
            (current_particle[0] <= 0) or (current_particle[0] >= grid_size[0]) or \
            (current_particle[1] <= 0) or (current_particle[1] >= grid_size[1])):
                particle_active = False
                count -= 1
                #print("shooted!")
                #time.sleep(5)
                shooted_count += 1
                return
        p_sub = p_extended[int(round(current_particle[0])):int(round(current_particle[0]))+2*critical_distance+1, \
                           int(round(current_particle[1])):int(round(current_particle[1]))+2*critical_distance+1, \
                           int(round(current_particle[2])):int(round(current_particle[2]))+2*critical_distance+1]
        f = np.array([0.,0.,-mass*gravity])
        v = velocity(f)
        #v = velocity_damp_factor*v
        #print("V: " + str(v))
        #current_particle[6:] = f
        #print("No neighbor")
        #print(current_particle[:3])
    #if (len(idx[0]) and (tmp[idx[0][0]] > minimum_distance)):
     #   current_particle[:3] += 0.5*minimum_distance*(p_sub[idx[0][0],:3] - current_particle[:3])
        
        
    else:
##        delta = time_step*current_particle[3:6]+0.5*(time_step**2)/mass*current_particle[6:]
##        if (np.sum(delta) > radius):
##            delta -= damp_factor*delta
##        current_particle[:3] += delta
        time_step = fine_step
        current_position = current_particle[:3]+time_step*current_particle[3:6]+(time_step**2)/mass*current_particle[6:]
        if ((current_position[2] <= 0) or (current_position[2] > grid_size[2]) or \
            (current_position[0] <= 0) or (current_position[0] >= grid_size[0]) or \
            (current_position[1] <= 0) or (current_position[1] >= grid_size[1])):
                particle_active = False
                #print("shooted!")
                #time.sleep(5)
                count -= 1
                shooted_count += 1 
                return
        p_sub = p_extended[int(round(current_position[0])):int(round(current_position[0]))+2*critical_distance+1, \
                           int(round(current_position[1])):int(round(current_position[1]))+2*critical_distance+1, \
                           int(round(current_position[2])):int(round(current_position[2]))+2*critical_distance+1]
        f = np.array([np.sum(np.multiply(fx,p_sub.astype(float))),np.sum(np.multiply(fy,p_sub.astype(float))),np.sum(np.multiply(fz,p_sub.astype(float)))])
        #print(fz[critical_distance,:,:])
        #print(np.multiply(fz[critical_distance,:,:],p_sub[critical_distance,:,:].astype(float)))
        f += np.array([0.,0.,-mass*gravity])
        if ((np.max(current_particle[3:5]) > 0) and (f[2] == 0)):
            f -= np.array([current_particle[3]/abs(np.max(current_particle[3:5])),current_particle[4]/abs(np.max(current_particle[3:5])),0])*mu_k*abs(f[2])
        v = velocity(f)
        #v = velocity_damp_factor*v
        near_particles = np.multiply(r_lut,p_sub.astype(float))
        if (np.count_nonzero(near_particles) > 0):
            #nearest_particle = np.argmin(near_particles[np.nonzero(near_particles)])
            nearest_particle = np.where(near_particles == np.amin(near_particles[np.nonzero(near_particles)]))[:3]
            #collide(np.array(nearest_particle).T[0])
            #print("nearest")
            #print(near_particles[nearest_particle])
            if ((near_particles[nearest_particle].all() > 0) and (near_particles[nearest_particle].any() < 2*rm)):
                c_factor = np.log((coarse_step/fine_step)-1)
                time_step = coarse_step / (1 + np.exp(c_factor - 2.0*c_factor/(2*rm-1.0)*(near_particles[nearest_particle][0] - 1.0)))
                #print(time_step)
            
        #print(np.amin(near_particles[np.nonzero(near_particles)]))
##        if (np.count_nonzero(near_particles) > 0):
##            #nearest_particle = np.argmin(near_particles[np.nonzero(near_particles)])
##            nearest_particle = np.where(near_particles == np.amin(near_particles[np.nonzero(near_particles)]))
##            #print("nearest")
##            if (np.alen(np.array(nearest_particle)[0]) > 0):
##                if ((near_particles[nearest_particle][0] > 0) and (near_particles[nearest_particle][0] <= minimum_distance)):
##                    print("Collision detected")
##                    #np.dot(v,(nearest_particle - current_particle[:3]))
##                    collide(np.array(nearest_particle).T[0])
        #print("neighbor")
        #if (abs(current_particle[6]) >= tolerance_z):
            #current_particle[2] = np.ceil(current_particle[2]/rm)*rm
            #current_particle[3] = 0.
            #f[0] = 0.
        #if (abs(current_particle[7]) >= tolerance_z):
            #current_particle[2] = np.ceil(current_particle[2]/rm)*rm
            #current_particle[4] = 0.
            #f[1] = 0.
        #if (abs(current_particle[8]) >= tolerance_z):
            #current_particle[2] = np.ceil(current_particle[2]/rm)*rm
            #current_particle[5] = 0.
            #current_particle[8] = 0.
            #f[2] = 0.
##        if (current_particle[2] <= radius):
##            current_particle[2] = radius
##            current_particle[5] = 0.
##            current_particle[8] = 0.
        
        #print(current_particle[:3])
        #print(current_particle[8])
        #print(osc_v)
##        for n in range(3):
##            if ((np.sign(v[n])*np.sign(current_particle[3+n])) < 0):
##                damp_factor[n] = 0.0
##                #print("Bounced")
##                #print(np.sum(v**2))
##        #print(v)
##        for n in range(3):
##            if ((np.sign(f[n])*np.sign(current_particle[6+n])) < 0):
##                damp_factor[n] = force_damp_factor
        #print(f)
        current_particle[:3] = current_position
        #if (((np.sign(f[2])*np.sign(current_particle[6+2])) <= 0) and (np.sum(current_particle[3:5]**2) <= tolerance_v**2)):
        if (((np.sign(f[2])*np.sign(current_particle[6+2])) <= 0) and (np.sum(current_particle[3:5]**2) <= tolerance_v**2)):
            if (np.sum((eq_point-current_particle[:3])**2) < position_change/10):
                stuck_timeout -= 1
            else:
                eq_point = current_particle[:3]
                stuck_timeout = stuck_timeout_init
                damp_factor = np.ones(3,dtype=float)
            print (stuck_timeout)
            if (stuck_timeout == 0):
                print("Stuck")
                if (p_lattice[int(round(current_position[0])),int(round(current_position[1])),int(round(current_position[2])),3] == 0):
                    p_lattice[int(round(current_position[0])),int(round(current_position[1])),int(round(current_position[2])),3] = 1
                    particles[count-1] = [int(round(current_particle[0])),int(round(current_particle[1])),int(round(current_particle[2])), \
                                          (int(round((current_particle[2]/rm)) % len(colors)))]
                    particle_active = False
                    #print("Atom not there")
                else:
                    stuck_timeout += 1
        if (((np.sign(f[0])*np.sign(current_particle[6+0])) <= 0) and \
            ((np.sign(f[1])*np.sign(current_particle[6+1])) <= 0) and \
            ((np.sign(f[2])*np.sign(current_particle[6+2])) <= 0)):
            if (np.sum((eq_point-current_particle[:3])**2) < position_change):
                timeout -= 1
            else:
                eq_point = current_particle[:3]
                timeout = timeout_init
                damp_factor = np.ones(3,dtype=float)
            #v = v*velocity_damp_factor
            print (timeout)
            if (timeout == 0):
                print("Force Balanced")
                if (p_lattice[int(round(current_particle[0])),int(round(current_particle[1])),int(round(current_particle[2])),3] == 0):
                    p_lattice[int(round(current_particle[0])),int(round(current_particle[1])),int(round(current_particle[2])),3] = 1
                    particles[count-1] = [int(round(current_particle[0])),int(round(current_particle[1])),int(round(current_particle[2])), \
                                          (int(round((current_particle[2]/rm)) % len(colors)))]
                    particle_active = False
                    #print("Atom not there")
                else:
                    timeout += 1
        
        #if (p_lattice[int(round(current_particle[0])),int(round(current_particle[1])),int(round(current_particle[2])),3] == 1):
        #    current_particle[3:6] = -1.0*current_particle[3:6]
        
    #if (len(idx[0]) and (tmp[idx[0][0]] > minimum_distance)):
     #   current_particle[:3] += 0.5*minimum_distance*(p_sub[idx[0][0],:3] - current_particle[:3])
        #print(current_particle[:3])
        #if ((np.sum(current_particle[3:6]**2) < tolerance) and (np.sum(current_particle[6:]**2) < tolerance) and \
        #    (p_lattice[int(round(current_particle[0])),int(round(current_particle[1])),int(round(current_particle[2])),3] == 0)):
        #if (p_lattice[int(round(current_particle[0])),int(round(current_particle[1])),int(round(current_particle[2])),3] == 0):
         #   p_lattice[int(round(current_particle[0])),int(round(current_particle[1])),int(round(current_particle[2])),3] = 1
          #  particles[count-1] = [int(round(current_particle[0])),int(round(current_particle[1])),int(round(current_particle[2]))]
           # particle_active = False
           
        #print(run_time)
           
        #print("Velocity: " + str(current_particle[3:6]))
        #print("Force: " + str(current_particle[6:]))
        #print("Pos: " + str(current_particle[:3]))
        #print("Velo: " + str(current_particle[3:6]))
    
    for n in range(3):
        if ((np.sign(v[n])*np.sign(current_particle[3+n])) < 0):
            damp_factor[n] = (1.0 - force_damp_factor)
            #print("Bounced")
            #print(np.sum(v**2))
        if ((np.sign(f[n])*np.sign(current_particle[6+n])) <= 0):
            damp_factor[n] = (1.0 - force_damp_factor)
    v = np.multiply(damp_factor, v)
    current_particle[6:] = f
    current_particle[3:6] = v
##    print("time step: " + str(time_step))
##    print("Z: " + str(current_particle[2]))
##    print("Velocity: " + str(current_particle[3:6]))
##    print("Force: " + str(current_particle[6:]))
    #print(count)
    
    run_time += time_step
def update_particles():
    global particles, count, particle_active
    if ((count < particles_count) and (not particle_active)):
        x,y,z = float(grid_size[0])*np.random.random(),float(grid_size[1])*np.random.random(),float(grid_size[2])
        add_particle(x,y,z)
    elif (particle_active):
        move_particle()
    #print count
def save_image(name):
    frame_buffer = ( GLubyte * (3*screen_size[0]*screen_size[1]))(0)
    glReadPixels(0, 0, screen_size[0], screen_size[1], GL_RGB, GL_UNSIGNED_BYTE, frame_buffer)
    image = Image.frombytes(mode="RGB", size=(screen_size[0], screen_size[1]), data=frame_buffer)     
    image = image.transpose(Image.FLIP_TOP_BOTTOM)
    addr = r'D:\Documents\Classes\EE216-Nanomaterials\Images\frame-'+name+'.png'
    image.save(addr, 'PNG')
def mousemotion(x,y):
    global mouseX
    global mouseY
    mouseX = x
    mouseY = y
def chaptrack():
    global mouseX
    global mouseY
    global wd
    global ht
    global aff
    dx = (mouseX-wd/2)/brake
    dy = (mouseY-ht/2)/brake
    glMatrixMode(GL_MODELVIEW)
    glPushMatrix()
    glLoadIdentity()
    glRotatef(dx,0,1.0,0.0)
    glRotatef(dy,1.0,0.0,0.0)
    glMultMatrixf(aff)
    aff = glGetFloatv(GL_MODELVIEW_MATRIX)
    glPopMatrix()
def motion():
    return 0
def reshape(width, height):
    global wd
    global ht
    glClearColor(0.0, 0.0, 0.0, 0.0)
    if height == 0:
        height = 1
    wd = width
    ht = height
    glViewport(0,0,wd,ht)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    if wd<=ht:
        glOrtho(-axrng,axrng,-axrng*ht/wd,axrng*ht/wd,-
        axrng,axrng)
    else:
        glOrtho(-axrng*wd/ht,axrng*wd/ht,-axrng,axrng,-
        axrng,axrng)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
def display():
    global particles, count, particle_active, timer, shooted_count, run_time, prev_time
    chaptrack()
    t0 = time.time()
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
    glPushMatrix()
    color = [1.0,0.,0.,100.0]
    glMaterialfv(GL_FRONT,GL_DIFFUSE,color)
    glLineWidth(1)
    glDisable(GL_LIGHTING)
    glColor3f(.5,.0,.0)
    glBegin(GL_LINES)
    for line in grid_box_lines:
        for vertex in line:
            glVertex3fv(grid_box_verticies[vertex])
    glEnd()
    glPopMatrix()
    glEnable(GL_LIGHTING)
    color = [1.,1.,1.,1.0]
    glMaterialfv(GL_FRONT,GL_DIFFUSE,color)
    update_particles()
    for n in range(initial_count):
        glPushMatrix()
        glTranslatef(particles[n,0],particles[n,1],particles[n,2])
        glutSolidSphere(radius,20,20)
        #print(particles[n,:3])
        glPopMatrix() 
    #color = [0.,.1,.5,1.0]
    #glMaterialfv(GL_FRONT,GL_DIFFUSE,color)
    for m in range(initial_count,(count-1)):
        glMaterialfv(GL_FRONT,GL_DIFFUSE,colors[particles[m,3]])
        glPushMatrix()
        glTranslatef(particles[m,0],particles[m,1],particles[m,2])
        glutSolidSphere(radius,20,20)
        #print(particles[n,:3])
        glPopMatrix()
    if (particle_active):
        color = [1.,.1,0.,1.0]
        glMaterialfv(GL_FRONT,GL_DIFFUSE,color)
        glPushMatrix()
        glTranslatef(current_particle[0],current_particle[1],current_particle[2])
        glutSolidSphere(radius,20,20)
        glDisable(GL_LIGHTING)
        glColor3f(.2,.2,.0)
        glutWireSphere(radius+critical_distance,20,20)
        glColor3f(1.0,.0,.0)
        glBegin(GL_LINES)
        glVertex3fv((0,0,0))
        glVertex3fv((1e-22*current_particle[6],0,0))
        glEnd()
        glColor3f(.0,1.0,.0)
        glBegin(GL_LINES)
        glVertex3fv((0,0,0))
        glVertex3fv((0,1e-22*current_particle[7],0))
        glEnd()
        glColor3f(.0,.0,1.0)
        glBegin(GL_LINES)
        glVertex3fv((0,0,0))
        glVertex3fv((0,0,1e-22*current_particle[8]))
        glEnd()
        glColor3f(1.0,1.0,1.0)
        glBegin(GL_LINES)
        glVertex3fv((0,0,0))
        glVertex3fv((1e-22*current_particle[6],1e-22*current_particle[7],1e-22*current_particle[8]))
        glEnd()
        glEnable(GL_LIGHTING)
        glPopMatrix()
    #glTranslatef(10*np.random.random(),10*np.random.random(),1*np.random.random())
    #glutSolidSphere(radius,10,10)
    glMatrixMode(GL_PROJECTION)
    matrix = glGetDouble( GL_PROJECTION_MATRIX )
    glLoadIdentity()
    glOrtho(0, screen_size[0], 0, screen_size[1], -1.0, 1.0)
    glMatrixMode(GL_MODELVIEW)
    glPushMatrix()
    glLoadIdentity()
    glDisable(GL_LIGHTING)
    glColor3f(.0,1.,.0)
    y = screen_size[1]-20
    glRasterPos2f(10, y)
    for c in name:
        glutBitmapCharacter(GLUT_BITMAP_9_BY_15, ord(c))
    y -= 20
    glColor3f(.5,1.0,.0)
    glRasterPos2f(10, y)
    for c in ">> Particles Count: " + str(count-initial_count):
        glutBitmapCharacter(GLUT_BITMAP_9_BY_15, ord(c))
    y -= 20
    glColor3f(1.,.0,.0)
    glRasterPos2f(10, y)
    for c in ">> Shooted Particles Count: " + str(shooted_count):
        glutBitmapCharacter(GLUT_BITMAP_9_BY_15, ord(c))
    y -= 20 
    glColor3f(1.,1.0,.0)
    glRasterPos2f(10, y)
    if (count < particles_count):
        timer = time.clock()
        #run_time += time_step
    for c in ">> Time Elapsed: " + "%.3f" % (run_time * 1e6) + " [usec]":
        glutBitmapCharacter(GLUT_BITMAP_9_BY_15, ord(c))
##    for c in ">> Time Elapsed: " + "%.3f" % timer + " [sec]":
##        glutBitmapCharacter(GLUT_BITMAP_9_BY_15, ord(c))
    glPopMatrix()
    glMatrixMode(GL_PROJECTION)
    glLoadMatrixd( matrix )
    glMatrixMode(GL_MODELVIEW)
    glEnable(GL_LIGHTING)
    t1 = time.time()
    glutSwapBuffers()
    #if ((run_time * 1e6 - prev_time) >= 1):
     #   print("frame-" + str(int(run_time * 1e6)) + " saved")
      #  save_image(str(int(run_time * 1e6)))
       # prev_time += 1
##    if (time_step/speed-(t1-t0) > 0):
##        time.sleep(time_step/speed-(t1-t0))
##    else:
##        time.sleep(0)

if __name__ == '__main__': main()
