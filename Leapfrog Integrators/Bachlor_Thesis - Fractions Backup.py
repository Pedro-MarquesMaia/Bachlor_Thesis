import numpy as np
import matplotlib.pyplot as plt
import math
import random
import time
import os
import imageio
import copy as c
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
        self.scale = Fraction(1000000000000)
        self.step = 2**17  #262144
        self.tend = 5.5*20000000#5*1000*17300000
        self.e = 1000000
        self.framepositions = []
        self.framepositions2 = []
        self.colors = []
        self.energy = []
        self.momentum = []
        self.angularmomentum = []
        self.G = Fraction(6.67*(10**(-11)))
        for particle in self.particles:
            self.colors.append(particle.color)
    def calc_energy(self):
        kenergy = 0
        penergy = 0
        for p in self.particles:
            kenergy += 0.5*p.m*(p.velocity.length())**2
        for n,part in enumerate(self.particles):
            pot = 0
            for p in self.particles:
                if p.number != part.number:
                    lower = ((p.position - part.position).length())
                    pot += (self.G*part.m*p.m)/lower
                else:
                    continue
            penergy += pot
        return 0.5*penergy + kenergy

    def create_copy(self):
        x = Universe()
        newlist = []
        newlist2 = []
        x.framepositions = c.deepcopy(self.framepositions)
        x.energy = c.deepcopy(self.energy)
        x.momentum = c.deepcopy(self.momentum)
        x.angularmomentum = c.deepcopy(self.angularmomentum)
        for p in self.particles:
            part = particle(mass=c.deepcopy(p.m),velocity=p.velocity*1,position=p.position*1,number=p.number*1,color=c.deepcopy(p.color))
            newlist.append(part)
            newlist2.append(part.color)
        x.particles = newlist
        x.colors = newlist2
        return x

    def calc_momentum(self):
        moment = vector(0,0,0)
        for part in self.particles:
            moment += part.velocity*part.m
        return moment.length()

    def calc_angular(self):
        angular = vector(0,0,0)
        for part in self.particles:
            angular += Math.crossproduct(part.position,part.velocity*part.m)
        return angular.length()

    def Force(self,p1, p2):
        m2 = p2.m
        distance = p1.position - p2.position
        Force = distance*1
        factor = Fraction(-(self.G*m2)/(((distance.length())**2) + self.e**2)**(3/2))
        Force = Force * factor
        return Force

    def timeevolution(self):
        t = 0
        positions = []
        n = 1
        while t < self.tend:
            n += 1
            position = []
            t += self.step
            self.leapfrog(Fraction(self.step))
            self.energy.append(self.calc_energy())
            self.momentum.append(self.calc_momentum())
            self.angularmomentum.append(self.calc_angular())
            for particle in self.particles:
                coordinates = (particle.position[0],particle.position[1],particle.position[2])
                position.append(coordinates)
            positions.append(position)
        self.framepositions = positions
        print("done: " + str(t))

    def remove_particle(self):
        amount = len(self.particles)
        particle_number = random.randint(0,amount-1)
        del self.particles[particle_number]
        del self.colors[particle_number]
        return particle_number

    
    def leapfrog(self,t):
        for particle in self.particles:
            R = particle.position + particle.velocity*int(t/2)
            particle.position = R
        for particle in self.particles:
            A = vector(0,0,0)
            for p in self.particles:
                if p.number == particle.number:
                    continue
                else:
                    A += (self.Force(particle, p))
            V = particle.velocity + A*t
            particle.velocity = V
        for particle in self.particles:
            R = particle.position + particle.velocity*int(t/2)
            particle.position = R

    def get_coordinates(particle):
        x = (particle.x,particle.y,particle.z)
        return x
    def timeevolute_back(self):
        t = math.ceil(self.tend/self.step)*self.step
        n = 0
        positions = []
        while t > 0:
            n += 1
            position = []
            t -= self.step
            self.leapfrog(Fraction(-self.step))
            self.energy.append(self.calc_energy())
            self.momentum.append(self.calc_momentum())
            self.angularmomentum.append(self.calc_angular())
            for particle in self.particles:
                coordinates = (particle.position[0], particle.position[1], particle.position[2])
                position.append(coordinates)
            positions.append(position)
        self.framepositions2 = positions
        print("done: " + str(t))
def create_cluster(n,colour,mass,vx,seed,start):
    particles = []
    np.random.seed(seed)
    for i in range(start, start + n):
        radius = 7.5*10**(10)
        x = np.random.uniform(-radius, radius)
        y = np.random.uniform(-radius, radius)
        z = np.random.uniform(-radius, radius)
        position = vector(Fraction(x), Fraction(y), Fraction(z))
        velocity = vector(vx,0,0)
        Particle = particle(int(mass), velocity, position, i,colour)
        particles.append(Particle)
    return particles


def shift_cluster(cluster,shift):
    for particle in cluster:
        particle.position.x = particle.position.x + Fraction(shift)
    return cluster

def create_universe():
    cluster1 = create_cluster(8,"red",2*(10**27),-0,200,0)
    cluster2 = create_cluster(8,"blue",2*(10**27),0,300,8)
    cluster1 = shift_cluster(cluster1,7.6*10**(10))
    cluster2 = shift_cluster(cluster2,-7.6*10**(10))
    return cluster1 + cluster2


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

def plot_quantities_vs_steps(values,quantity):
    steps = list(range(1, len(values) + 1))
    plt.figure(figsize=(10, 6))
    plt.plot(steps, values, marker='o', linestyle='-', color='b', label=quantity)
    plt.xlabel('Iteration Step')
    plt.ylabel(quantity)
    plt.title(quantity + ' vs. Iteration Steps')
    plt.grid(True)
    plt.legend()
    plt.show()
def create_particle_video(framepositions,output_filename, skip_frames=1, fps=10):
    os.makedirs('temp_images', exist_ok=True)
    num_frames = len(framepositions)
    images = []
    for i in range(0, num_frames, skip_frames):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        frame = framepositions[i]
        frame_array = np.array(frame)
        ax.scatter(frame_array[:, 0], frame_array[:, 1], frame_array[:, 2], marker='o',color=colours1)
        ax.set_xlim([-5*(10**11), 5*(10**11)])
        ax.set_ylim([-5*(10**11), 5*(10**11)])
        ax.set_zlim([-5*(10**11), 5*(10**11)])
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        image_path = f'temp_images/frame_{i}.png'
        plt.savefig(image_path)
        plt.close(fig)
        images.append(image_path)
    with imageio.get_writer(output_filename, fps=fps) as writer:
        for image_path in images:
            image = imageio.imread(image_path)
            writer.append_data(image)
    for image_path in images:
        os.remove(image_path)
    os.rmdir('temp_images')
    print(f'Video saved as {output_filename}')
start_time = time.time()
print(f"Start time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}")
Uni = Universe()
Uni.timeevolution()
Uni2 = Uni.create_copy()
Uni.timeevolute_back()
frames = 1*Uni.framepositions
frames2 = 1*Uni.framepositions2
print("Length:")
print(len(frames))
colours1 = Uni.colors
end_time = time.time()
print(f"End time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}")
duration = end_time - start_time
print(f"Duration: {duration:.2f} seconds")
plot_particle_positions2(frames,colours1)
plot_particle_positions2(frames2,colours1)
plot_quantities_vs_steps(Uni.energy,"Total Energy")
plot_quantities_vs_steps(Uni.momentum, "Total Momentum")
plot_quantities_vs_steps(Uni.angularmomentum, "Total Angular Momentum")
create_particle_video(frames,"vid1.mp4")
create_particle_video(frames2,"vid2.mp4")
ind = Uni2.remove_particle()
Uni2.timeevolute_back()
frames3 = 1*Uni2.framepositions
frames4 = 1*Uni2.framepositions2
colours2 = Uni2.colors
plot_particle_positions2(frames3,colours1)
plot_particle_positions2(frames4,colours2)
plot_quantities_vs_steps(Uni2.energy,"Total Energy")
plot_quantities_vs_steps(Uni2.momentum, "Total Momentum")
plot_quantities_vs_steps(Uni2.angularmomentum, "Total Angular Momentum")
create_particle_video(frames3,"vid3.mp4")
del colours1[ind]
create_particle_video(frames4,"vid4.mp4")