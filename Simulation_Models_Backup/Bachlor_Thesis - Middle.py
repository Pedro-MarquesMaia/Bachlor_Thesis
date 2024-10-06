import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math
import random
import time
import os
import imageio
import copy as c
from fractions import Fraction
from decimal import Decimal, getcontext

getcontext().prec = 50

class particle():
    def __init__(self,mass,velocity,position,number,color):
        self.m = mass
        self.velocity = velocity
        self.position = position
        self.number = number
        self.color = color

class vector():
    scale = 1
    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z
    def __repr__(self):
        return "(" + str(self.x) + "," + str(self.y) + "," + str(self.z) + ")"
    def dec(self):
        self.x = Decimal(self.x)
        self.y = Decimal(self.y)
        self.z = Decimal(self.z)
    def length(self):
        return Math.norm(self)
    def length_new(self):
        return Math.norm_new(self)
    def length(self):
        return Math.norm(self)
    def __add__(self, other):
        return vector(self.x + other.x,self.y + other.y, self.z + other.z)   
    def __sub__(self,other):
        return vector(self.x - other.x, self.y - other.y, self.z -other.z)
    def __rmul__(self,other):
        if type(other) != type(self):
            return vector(self.x*other, self.y*other, self.z*other)
        else:
            return Math.scalar(self,other)
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

    def intify(self):
        self.x = int(round(self.x))
        self.y = int(round(self.y))
        self.z = int(round(self.z))
        return vector(self.x, self.y, self.z)

    def intify_div(self, var):
        return vector(self.x // var, self.y // var, self.z // var)

    def discretize(self):
        new = self.intify_div(self.scale)
        new = new * self.scale
        return new

        
    
class Math():
    def scalar(a,b):
        res = a.x*b.x + a.y*b.y + a.z*b.z
        return res
    def norm(a):
        res = Math.scalar(a,a)
        return np.sqrt(float(res))
    def crossproduct(a,b):
        x = a.y*b.z - a.z*b.y
        y = a.z*b.x - a.x*b.z
        z = a.x*b.y - a.y*b.x
        return vector(x,y,z)
    def scalar_new(a,b):
        res = int(a.x*b.x) + int(a.y*b.y) + int(a.z*b.z)
        return res
    def unitvector(v):
        n = Math.norm(v)
        x = v.x/n
        y = v.y/n
        z = v.z/n
        return vector(x,y,z)
    def norm_new(a):
        res = Math.scalar_new(a,a)
        res = math.isqrt(res)
        return res
    def crossproduct_new(a,b):
        x = int(a.y)*int(b.z) - int(a.z)*int(b.y)
        y = int(a.z)*int(b.x) - int(a.x)*int(b.z)
        z = int(a.x)*int(b.y) - int(a.y)*int(b.x)
        return vector(x,y,z)
    
class Universe():
    def __init__(self):
        self.scale = 1000000000000
        vector.scale = self.scale
        self.particles = create_universe(self.scale)
        self.scale = 1000000000000
        self.step = 200   #262144
        self.tend = 20*(10**7)
        self.e = 1000
        self.framepositions = []
        self.framepositions2 = []
        self.colors = []
        self.energy = []
        self.momentum = []
        self.angularmomentum = []
        self.G = 6.67*10**(-11)
        for particle in self.particles:
            self.colors.append(particle.color)
        
    def calc_energy(self):
        kenergy = 0
        penergy = 0
        for p in self.particles:
            kenergy += 0.5*p.m*(self.to_float(p.velocity,self.scale).length())**2
        for n,part in enumerate(self.particles):
            pot = 0
            for p in self.particles:
                if p.number != part.number:
                    lower = np.sqrt((self.to_float(p.position,self.scale) - self.to_float(part.position,self.scale)).length()**2 + self.e**2) #lower = np.sqrt(((self.to_float(p.position,self.scale) - self.to_float(part.position,self.scale)).length())**2 + self.e**2)
                    pot += (self.G*part.m*p.m)/lower
                else:
                    continue
            penergy += pot
        return penergy + kenergy

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
        ang = vector(0,0,0)
        for part in self.particles:
            ang += Math.crossproduct(part.position.intify_div(self.scale),(part.velocity*part.m).intify_div(self.scale))
        ang = self.to_float(ang, self.scale)
        return ang.length()



    def Force(self,m2,p1, p2,t):
        G = (-1)*(6.67*10**(-11))*m2*self.scale
        distance = p1 - p2
        factor = G/((distance.length()**2 + self.e**2)**(3/2))
        Force = factor*distance*t
        Force.intify()
        return Force

    def timeevolution(self):
        t = 0
        positions = []
        n = 0
        print(((self.G * 2 * (10 ** 27)) ** (1 / 2)) * ((self.e) ** (-3 / 2)) * self.step)
        while t < self.tend:
            position = []
            t += self.step
            self.leapfrog(self.step)
            self.energy.append(self.calc_energy())
            self.momentum.append(self.calc_momentum())
            self.angularmomentum.append(self.calc_angular())
            for particle in self.particles:
                p = self.to_float(particle.position,self.scale)
                coordinates = (p[0],p[1],p[2])
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
            Vchange = vector(0,0,0)
            for p in self.particles:
                if p.number == particle.number:
                    continue
                else:
                    Vchange += (self.Force(p.m,self.to_float(particle.position,self.scale),self.to_float(p.position,self.scale),t))
            V = particle.velocity + Vchange
            particle.velocity = V
        for particle in self.particles:
            R = particle.position + particle.velocity*int(t/2)
            particle.position = R
            
    def to_int(self,v,scale):
        v.x = int(v.x * scale)
        v.y = int(v.y * scale)
        v.z = int(v.z * scale)
        return v

    def to_float(self,v,scale):
        try:
            vn = v * (1 / scale)
        except:
            vn = v * (1 / Decimal(scale))
        return vn

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
            self.leapfrog(int(-self.step))
            self.energy.append(self.calc_energy())
            self.momentum.append(self.calc_momentum())
            self.angularmomentum.append(self.calc_angular())
            for particle in self.particles:
                p = self.to_float(particle.position, self.scale)
                coordinates = (p[0], p[1], p[2])
                position.append(coordinates)
            positions.append(position)
        self.framepositions2 = positions
        print("done: " + str(t))

def create_cluster(n,colour,mass,vx,seed,start,scale):
    particles = []
    np.random.seed(seed)
    for i in range(start, start + n):
        # Random position within a cubical volume (length 6*10**(12))
        radius = 6*10**(12)
        x = np.random.uniform(-radius, radius)
        y = np.random.uniform(-radius, radius)
        z = np.random.uniform(-radius, radius)
        position = vector(int(x)*scale*0, int(y)*scale*0, int(z)*scale*0)

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
        velocity = vector(0,vx*scale,0)
        Particle = particle(mass, velocity, position, i,colour)
        particles.append(Particle)
    """if is_bound(particles):
        return particles
    else:
        return create_cluster(n, colour)"""
    return particles     
        
def is_bound(particles):
    T = 0
    U = 0
    for particle in particles:
        T += 0.5*particle.m*(particle.velocity*(1/100000000)).length()**2
    for n,particle in enumerate(particles):
        for i,particle2 in enumerate(particles):
            if n == i:
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

def shift_cluster(cluster,shift,scale):
    for particle in cluster:
        particle.position.x = particle.position.x + int(shift)*scale
    return cluster

def create_universe(scale):
    cluster1 = create_cluster(1, "red", 2 * (10 ** 30), 0, 200, 0, scale)
    cluster2 = create_cluster(1, "blue", 6 * (10 ** 24), 35000, 300, 1, scale)
    cluster1 = shift_cluster(cluster1,0,scale)
    cluster2 = shift_cluster(cluster2,1.5*10**11,scale)
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


def create_particle_video(framepositions,output_filename, skip_frames=300, fps=30):
    os.makedirs('temp_images', exist_ok=True)
    num_frames = len(framepositions)
    images = []
    for i in range(0, num_frames, skip_frames):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        frame = framepositions[i]
        frame_array = np.array(frame)
        ax.scatter(frame_array[:, 0], frame_array[:, 1], frame_array[:, 2], marker='o',color=colours1)
        ax.set_xlim([-2*(10**13), 2*(10**13)])
        ax.set_ylim([-2*(10**13), 2*(10**13)])
        ax.set_zlim([-2*(10**13), 2*(10**13)])
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

def plot_fractional_change(data, title):
    ref = data[1]
    plt.figure(figsize=(8, 7))
    output = []
    x = []
    for i in range(1, len(data)):
        try:
            fractional_change = (data[i] - ref )/ ref
            output.append(fractional_change)
            x.append(i)
        except:
            print("One")
            fractional_change = 0
            output.append(fractional_change)
            x.append(i)
    plt.plot(x, output, marker='o', color='b', linestyle="None")
    plt.title(title)
    plt.xlabel('Index (Starting from 1)')
    plt.ylabel('Fractional Change')
    plt.grid(True)
    plt.show()
    
def plot_particle_paths(frames):
    num_particles = len(frames[0])
    colors = ['r', 'g', 'b', 'y']
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # Extract and plot the paths for each particle
    for i in range(num_particles):
        x = []
        y = []
        z = []
        for frame in frames:
            pos = frame[i]
            x.append(pos[0])
            y.append(pos[1])
            z.append(pos[2])
        ax.plot(x, y, z, color=colors[i], label=f'Particle {i+1}')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()
    plt.show()
start_time = time.time()
print(f"Start time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}")
Uni = Universe()
Uni.timeevolution()
Uni2 = Uni.create_copy()
Uni.timeevolute_back()
frames = 1*Uni.framepositions
frames2 = 1*Uni.framepositions2
plot_particle_paths(frames)
print("Length:")
print(len(frames))
colours1 = Uni.colors
end_time = time.time()
print(f"End time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}")
duration = end_time - start_time
print(f"Duration: {duration:.2f} seconds")
ind = Uni2.remove_particle()
Uni2.timeevolute_back()
frames3 = 1*Uni2.framepositions
frames4 = 1*Uni2.framepositions2
colours2 = Uni2.colors
plot_particle_positions2(frames,colours1)
plot_particle_positions2(frames2,colours1)
plot_quantities_vs_steps(Uni.energy,"Total Energy")
plot_quantities_vs_steps(Uni.momentum, "Total Momentum")
plot_quantities_vs_steps(Uni.angularmomentum, "Total Angular Momentum")
plot_particle_positions2(frames3,colours1)
plot_particle_positions2(frames4,colours2)
plot_quantities_vs_steps(Uni2.energy,"Total Energy")
plot_quantities_vs_steps(Uni2.momentum, "Total Momentum")
plot_quantities_vs_steps(Uni2.angularmomentum, "Total Angular Momentum")
plot_fractional_change(Uni.energy,"Fractional change Energy")
plot_fractional_change(Uni.momentum,"Fractional change Momentum")
plot_fractional_change(Uni.angularmomentum,"Fractional change Angular Momentum")
plot_fractional_change(Uni2.energy,"Fractional change Energy")
plot_fractional_change(Uni2.momentum,"Fractional change Momentum")
plot_fractional_change(Uni2.angularmomentum,"Fractional change Angular Momentum")
create_particle_video(frames,"vid1.mp4")
create_particle_video(frames2,"vid2.mp4")
create_particle_video(frames3,"vid3.mp4")
del colours1[ind]
create_particle_video(frames4,"vid4.mp4")

