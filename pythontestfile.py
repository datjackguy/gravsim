import numpy as np
import random

G=6.67430e-11
epsilon = 0.1

class particle():
    def __init__(self,mass,pos,vel,frames):
        self.mass = mass
        self.pos = np.array(pos)
        self.vel = np.array(vel)
        # self.allPos = np.empty((frames,3),dtype=float)
        self.allPos = np.full((frames,3),np.nan)
        # self.allMass = np.empty((frames+1,1),dtype=float)
        self.allMass = np.full((frames,1),np.nan)
        self.merged = False
        self.prevVel = np.array(vel)
    def getMass(self):
        return self.mass
    def getAllPos(self):
        return self.allPos
    def getPos(self):
        return self.pos
    def getVel(self):
        return self.vel
    def calcPos(self,dt,frame):
        self.allPos[frame]=self.pos
        self.pos = self.pos + self.vel * dt
    def calcVel(self,dt,netForce):
        self.prevVel = self.vel
        self.vel = self.vel + netForce * dt / self.mass
    def merge(self):
        self.merged = True
    def getMerged(self):
        return self.merged
    def setPos(self,pos):
        self.pos = pos
    def setVel(self,vel):
        self.prevVel = vel
        self.vel = vel
    def setMass(self,mass):
        self.mass = mass
    def setCurrentMass(self,mass,frame):
        self.allMass[frame]=mass
    def getAllMass(self):
        return self.allMass
    def getPrevVel(self):
        return self.prevVel

def initial_mass_function(m):
    if m < 0.08:
        return 0.0
    elif m < 0.5:
        return m ** -1.3
    elif 0.5 < m < 1:
        return m ** -2.3
    else:
        return m ** -2.7
    
def generate_masses(N):
    masses = []
    while len(masses) < N:
        #Random Mass
        x_candidate = np.random.uniform(0, 20)
        #Probability
        y_candidate = np.random.uniform(0, 1)
        #Check if mass accepted
        if y_candidate <= initial_mass_function(x_candidate):
            masses.append(x_candidate)
    return np.array(masses)

def randomSph(maxRad):
    rho=random.random()*maxRad
    # rho=maxRad
    theta=random.random()*np.pi
    phi=random.random()*2*np.pi
    x=rho*np.sin(theta)*np.cos(phi)
    y=rho*np.sin(theta)*np.sin(phi)
    z=rho*np.cos(theta)
    return x,y,z

def calcDistance(p1,p2):
    distance=((p1.getPos()[0]-p2.getPos()[0])**2+(p1.getPos()[1]-p2.getPos()[1])**2+(p1.getPos()[2]-p2.getPos()[2])**2)**(1/2)
    return distance

def calcForce(p1,p2):
    distance=((p1.getPos()[0]-p2.getPos()[0])**2+(p1.getPos()[1]-p2.getPos()[1])**2+(p1.getPos()[2]-p2.getPos()[2])**2)**(1/2)
    """if distance < 5e20:
        force=-G*p1.getMass()*p2.getMass()*(p1.getPos()-p2.getPos())/(distance**2+(epsilon)**2)**(3/2)
    else:
        force=-G*p1.getMass()*p2.getMass()*(p1.getPos()-p2.getPos())/(distance)**3"""
    # if distance < epsilon:
    #     force=-G*p1.getMass()*p2.getMass()*(p1.getPos()-p2.getPos())/(distance**2+(epsilon)**2)**(3/2)
    # else:
    #     force=-G*p1.getMass()*p2.getMass()*(p1.getPos()-p2.getPos())/(distance**2)**(3/2)
    force=-G*p1.getMass()*p2.getMass()*(p1.getPos()-p2.getPos())/(distance**2+(epsilon)**2)**(3/2)
    return force
    
print("IMF Check")
print(initial_mass_function(0.03))
print(initial_mass_function(0.3))
print(initial_mass_function(0.7))
print(initial_mass_function(1.2))

print("Masses Check")
print(generate_masses(10))

print("Random Sphere Check")
print(randomSph(10))

print("Distance Calc Check")
testParticle1 = particle(10, [1, 1, 1], [0, 0.2, 0], 10)
testParticle2 = particle(20, [3, 0, -2], [0.4, 0, 0.1], 10)
print(calcDistance(testParticle1,testParticle2))

print("Force Calc Check")
print(calcForce(testParticle1,testParticle2))

