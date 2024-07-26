import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math
import random
from fractions import Fraction

class particle():
    def __init__(self,mass,velocity,position,number,color):
        self.m = mass
        self.velocity = velocity
        self.position = position
        self.number = number
        self.color = color

class vector():
    def __init__(self,x,y,z):
        self.x = Fraction(x)
        self.y = Fraction(y)
        self.z = Fraction(z)
    def __repr__(self):
        return "(" + str(self.x) + "," + str(self.y) + "," + str(self.z) + ")"
    def length(self):
        return Math.norm(self)
    def __add__(self, other):
        return vector(self.x + other.x,self.y + other.y, self.z + other.z)   
    def __sub__(self,other):
        return vector(self.x - other.x, self.y - other.y, self.z -other.z)
    def __mul__(self,other):
        if type(other) != type(self):
            return vector(self.x*other, self.y*other, self.z*other)
        else:
            return Math.scalar(self,other)
    def __getitem__(self,num):
        if num == 0:
            return self.x
        elif num == 1:
            return self.y
        elif num == 2:
            return self.z
        
        
    
class Math():
    def scalar(a,b):
        res = a.x*b.x + a.y*b.y + a.z*b.z
        return res
    def norm(a):
        res = Math.scalar(a,a)
        return Fraction(np.sqrt(float(res)))
    def crossproduct(a,b):
        x = a.y*b.z - a.z*b.y
        y = a.z*b.x - a.x*b.z
        z = a.x*b.y - a.y*b.x
        return vector(x,y,z)
    def unitvector(v):
        n = Math.norm(v)
        x = v.x/n
        y = v.y/n
        z = v.z/n
        return vector(x,y,z)
    
class Universe():
    def __init__(self):
        self.particles = create_universe()
        self.step = Fraction(350000)
        self.tend = Fraction(4*1000*17300000)
        self.e = Fraction(1000000)
        self.framepositions = []
        self.framepositions2 = []
        self.colors = []
        for particle in self.particles:
            self.colors.append(particle.color)
        
    
    def Force(self,p1, p2):
        G = Fraction(6.67*10**(-11))
        m1 = p1.m
        m2 = p2.m
        distance = p1.position - p2.position
        Force = distance*1
        factor = -(G*m1*m2)/(((distance.length())**2) + self.e**2)**(3/2)
        Force = Force * factor
        return Force

    def timeevolution(self):
        t = 0
        positions = []
        while t < self.tend:
            position = []
            t += self.step
            self.leapfrog(Fraction(self.step))
            for particle in self.particles:
                coordinates = (particle.position[0],particle.position[1],particle.position[2])
                position.append(coordinates)
            positions.append(position)
        self.framepositions = positions
        print("done: " + str(t))

    def remove_particle(self):
        amount = len(self.particles)
        particle_number = random.randint(0,amount)
        del self.particles[particle_number]
        print(particle_number)

    
    def leapfrog(self,t):
        for particle in self.particles:
            R = particle.position + particle.velocity*(t/2)
            particle.position = R
        for particle in self.particles:
            A = vector(0,0,0)
            for number, p in enumerate(self.particles):
                if number == particle.number:
                    continue
                else:
                    A += (self.Force(particle, p)*(1/particle.m))
            V = particle.velocity + A*t
            particle.velocity = V
        for particle in self.particles:
            R = particle.position + particle.velocity*(t/2)
            particle.position = R
        
    def update_frame(frame):
        scatter.set_offsets(self.framepositions[frame])
        return scatter
            
        
    def get_coordinates(particle):
        x = (particle.x,particle.y,particle.z)
        return x
    def timeevolute_back(self):
        t = math.ceil(self.tend/self.step)*self.step
        print(t)
        positions = []
        while t > 0:
            position = []
            t -= self.step
            self.leapfrog(Fraction(-self.step))
            for particle in self.particles:
                coordinates = (particle.position[0], particle.position[1], particle.position[2])
                position.append(coordinates)
            positions.append(position)
        self.framepositions2 = positions
        print("done: " + str(t))
def create_cluster(n,colour,mass,vx,seed):
    particles = []
    np.random.seed(seed)
    for n in range(n):
        
        # Random position within a cubical volume (length 6*10**(12))
        radius = 6*10**(12)
        x = random.uniform(-radius, radius)
        y = random.uniform(-radius, radius)
        z = random.uniform(-radius, radius)
        position = vector(x, y, z)

        """if vx < 0:
            vx = vx #random.uniform(vx, 0)
            vy = random.uniform(-radius2, radius2)
            vz = random.uniform(-radius2, radius2)
            velocity = vector(vx, vy, vz)
        else:
            vx = vx #random.uniform(0, vx)
            vy = random.uniform(-radius2, radius2)
            vz = random.uniform(-radius2, radius2)
            velocity = vector(vx, vy, vz) """
        velocity = vector(vx,0,0)
        Particle = particle(mass, velocity, position, n,colour)
        particles.append(Particle)
    if is_bound(particles):
        return particles
    else:
        return create_cluster(n, colour)
    return particles     
        
def is_bound(particles):
    T = 0
    U = 0
    for particle in particles:
        T += 0.5*particle.m*particle.velocity.length()**2
    for n,particle in enumerate(particles):
        for i,particle2 in enumerate(particles):
            if n == particle2.number:
                continue
            else:
                val = particle.m*particle2.m
                radius = particle.position - particle2.position
                radius = radius.length()
                energy = val/radius
                U += energy
    if T - U < 0:
        print(T)
        print(U)
        return True
    else:
        return False

def shift_cluster(cluster,shift):
    for particle in cluster:
        particle.position.x = particle.position.x + shift
    return cluster

def create_universe():
    cluster1 = create_cluster(13,"red",2*(10**27),-5,80)
    cluster2 = create_cluster(13,"blue",2*(10**27),5,90)
    cluster1 = shift_cluster(cluster1,6.5*10**(12))
    cluster2 = shift_cluster(cluster2,-6.5*10**(12))
    #plot_particle_positions(cluster1 + cluster2)
    return cluster1 + cluster2

def plot_particle_positions(particles):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    colour = particles[0].color

    # Extract x, y, z coordinates from particles
    x_positions = [float(particle.position[0]) for particle in particles]
    y_positions = [float(particle.position[1]) for particle in particles]
    z_positions = [float(particle.position[2]) for particle in particles]
    colors = [particle.color for particle in particles]

    ax.scatter(x_positions, y_positions, z_positions, marker='o', label='Particles',color = colors)
    ax.scatter(1.6*10**13,1.6*10**13,1.6*10**13,marker="o",alpha=0)
    ax.scatter(-1.6*10**13,-1.6*10**13,-1.6*10**13,marker="o",alpha=0)
    
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Particle Positions')
    plt.legend()
    plt.show()

def plot_particle_positions2(frames, colors):
    first_frame = frames[0]
    last_frame = frames[-1]

    fig = plt.figure(figsize=(14, 6))

    # Plotting the first frame
    ax1 = fig.add_subplot(121, projection='3d')
    for (x, y, z), color in zip(first_frame, colors):
        ax1.scatter(float(x), float(y), float(z), color=color, label=color)
    ax1.set_title('First Frame')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')

    # Plotting the last frame
    ax2 = fig.add_subplot(122, projection='3d')
    for (x, y, z), color in zip(last_frame, colors):
        ax2.scatter(float(x), float(y), float(z), color=color, label=color)
    ax2.set_title('Last Frame')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')

    plt.show()
    
def Animate(frames,color,title):
    n = len(frames[0])
    colors = color
    print(len(frames))
    def update(frame):
        frame_coords = frames[frame]
        for i in range(n):
            x = float(frame_coords[i][0])
            y = float(frame_coords[i][1])
            z = float(frame_coords[i][2])
            line = lines[i]
            line.set_data([x], [y])
            line.set_3d_properties([z])
            line.set_color(colors[i])
        return lines
    
    def init():
        for line in lines:
            line.set_data([], [])
            line.set_3d_properties([])
        return lines
    
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.set(xlim3d=(-10**14, 10**14), xlabel='X')
    ax.set(ylim3d=(-10**14, 10**14), ylabel='Y')
    ax.set(zlim3d=(-10**14, 10**14), zlabel='Z')
    ax.set_title(title)
    lines = []
    for i in range(n):
        line, = ax.plot([], [], [], 'o')
        lines.append(line)
    ani = FuncAnimation(fig, update, frames=range(len(frames)), init_func=init, blit=True,interval=1)
    plt.show()

Uni = Universe()
Uni.timeevolution()
Uni.timeevolute_back()
frames = 1*Uni.framepositions
frames2 = 1*Uni.framepositions2
colours1 = Uni.colors
Uni = Universe()
Uni.timeevolution()
Uni.remove_particle()
Uni.timeevolute_back()
frames3 = 1*Uni.framepositions
frames4 = 1*Uni.framepositions2
colours2 = Uni.colors
#Successful
Animate(frames,colours1,"Forward")
Animate(frames2,colours1,"Backwards")
plot_particle_positions2(frames,colours1)
plot_particle_positions2(frames2,colours1)
#Unsuccesful
Animate(frames3,colours2,"Forward")
Animate(frames4,colours2,"Backwards")
plot_particle_positions2(frames3,colours2)
plot_particle_positions2(frames4,colours2)