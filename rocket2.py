import numpy as np
import matplotlib.pyplot as  plt
from matplotlib.animation import FuncAnimation 
import copy
import sys

# G force
G = 6.67430e-11
M = 5.972e24    # mass Earth
R = 6.371e6     # radius

rho0 = 1.225    # Density at surface
Cd = 0.6
g = 9.81
A = 75   # diameter of 2 m

class Rocket():
    def __init__(self, dt=0.01, rocket_mass=0, r=R,
                mTot=5e5,
                m1=0, t1=150, ve1=3000,
                m2=0, t2=300, ve2=4500
                ):
        # ship mass 
        self.mTot = mTot
        self.rocket_mass = mTot*rocket_mass

        # Ship variables
        self.r = r #6378e3
        self.v_r = 0
        self.theta = np.pi/2
        self.dtheta = 0
        self.t = 0
        self.alfa = 0    # Thrust angle against radial direction
        self.status = "Active"  # "Has exited", "Crashed", "Active"
        
        self.stage1 = Stage(mTot*m1, t1, ve1)
        self.stage2 = Stage(mTot*m2, t2, ve2)
        self.current_stage = self.stage1
        self.next_stage = self.stage2

    def update_rocket(self, obs, dt):
        if self.current_stage is not None:
            # Remove fuel, check active stage
            self.current_stage.update_stage(dt)
            self.mTot = self.rocket_mass + self.stage1.m + self.stage2.m
        return

    def Fg(self):
        return G*M*self.mTot/(self.r)**2
    
    def Fthrust(self):
        if self.current_stage == None:
            return 0
        else:
            return self.current_stage.get_thrust()
        
    def Fdrag(self):
        return self.getRho()*(self.v_r**2 + (self.r*self.dtheta)**2) * Cd * A/2

    def getRho(self):
        if self.current_stage == None:
            return 0
        if self.r < 11000:
            H = 8500
        elif self.r < 50e3 and self.r > 11e3:
            H = 7600
        elif self.r < 85e3 and self.r > 50e3:
            H = 6500
        elif self.r < 500e3 and self.r > 85e3:
            H = 15000
        else:
            return 0
        return rho0*np.exp(-(self.r-R)/H)
    
    def get_energy(self):
        K = (self.mTot*(self.v_r**2 + (self.r*self.dtheta)**2))/2
        #Gravitational pot energy
        U = -G*M*self.mTot/self.r
        return K + U
    
    def update_alfa(self):
        if self.current_stage == self.stage2:
            self.alfa = min((self.t-self.stage1.burnTime)/(self.stage2.burnTime)*np.pi/2, np.pi/3)
        return
            
    def accel_r(self):
        #print("Forces", self.Fg(), self.Fthrust(), self.Fdrag())
        return (-self.Fg() + self.Fthrust()*np.cos(self.alfa) - self.Fdrag()*np.cos(self.alfa))/self.mTot
    
    def accel_theta(self):
        return (self.Fthrust()*np.sin(self.alfa) - self.Fdrag()*np.sin(self.alfa))/self.mTot + 2*self.v_r*(self.r*self.dtheta)/self.r
    
class Stage:
    def __init__(self, m, t, v_e):
        self.m = m
        self.burnTime = t
        self.dmdt = self.m/self.burnTime
        self.v_e = v_e     # Exhaust velocity
        self.active = True
    
    def get_thrust(self):
        return self.dmdt*self.v_e
    
    # Update fuel, check runtime
    def update_stage(self, dt):
        self.m -= self.dmdt * dt
        if self.m < 0:
            self.m = 0
            self.active = False
        return

class Observables:
    def __init__(self):
        self.v_r = []
        self.r = [] 
        self.dtheta = []
        self.theta = []
        self.energy = []
        self.t = []
        self.stageswap_coord = []
        self.stageswap_time = []
        self.stageswap_energy = []

class Integrator:

    def __init__(self, _dt=0.01):
        self.dt = _dt

    def integrate(self, rocket, obs):

        self.step(rocket)

        obs.v_r.append(rocket.v_r)
        obs.r.append(rocket.r)
        obs.dtheta.append(rocket.dtheta)
        obs.theta.append(rocket.theta)
        obs.t.append(rocket.t)
        # removes fuelmass and sets status of stage to inactive if empty tank
        rocket.update_rocket(obs, self.dt)
        energy = rocket.get_energy()       
        obs.energy.append(energy)

        if rocket.current_stage is not None and rocket.current_stage.active == False:
            print("Current stage ended, r:", rocket.r)
            obs.stageswap_coord.append([rocket.r, rocket.theta])
            obs.stageswap_time.append(rocket.t)
            obs.stageswap_energy.append(energy)
            rocket.current_stage = rocket.next_stage
            rocket.next_stage = None

        # check crash into earth
        if rocket.r < R:
            print("Rocket crashed")
            rocket.status = "Crashed"

        if energy > 0:
            rocket.status = "Has exited"
       

    def simulate(self, rocket, obs, nn):
        #nn = self.nsteps//self.numStepsPerFrame
        #print("Integrating for"+str(nn*self.numStepsPerFrame)+"steps")
        for i in range(nn):
            self.integrate(rocket, obs)
            print(i)
            
            if rocket.status == "Crashed":
                print("Crashed")
                break
            #elif rocket.status == "Has exited":
                #for j in range(100):
                #    self.integrate(rocket, obs)
                #break
    def step(self, rocket):
        pass
# Doesnt fit the problem
class Verlet(Integrator):
    def step(self, rocket, obs):

        ddtheta = rocket.accel_theta/rocket.r
        new_theta = rocket.theta + rocket.dtheta * self.dt + ddtheta/2 * (self.dt)**2

        new_r = rocket.r + rocket.v_r * self.dt + rocket.accel_r/2 * (self.dt)**2

        new_dtheta = rocket.dtheta + rocket.accel_theta

class EulerCromer(Integrator):
    def step(self, rocket):
        rocket.update_alfa()
        accel_r = rocket.accel_r()
        accel_theta = rocket.accel_theta()

        rocket.v_r += accel_r * self.dt
        temp = rocket.r
        rocket.r += rocket.v_r * self.dt
        rocket.dtheta += accel_theta/temp * self.dt
        rocket.theta += rocket.dtheta * self.dt
# This is used
class RK4(Integrator):
    def step(self, rocket):
        # Radial
        rocket.update_alfa()
        accel_r = rocket.accel_r()

        accel_theta = rocket.accel_theta()
        a1 = accel_r * self.dt
        b1 = rocket.v_r * self.dt
        c1 = accel_theta/rocket.r * self.dt
        d1 = rocket.dtheta * self.dt

        temp_rocket = copy.deepcopy(rocket)
        temp_rocket.v_r += a1/2
        temp_rocket.r += b1/2
        temp_rocket.dtheta += c1/2
        temp_rocket.theta += d1/2
        temp_rocket.t += self.dt/2
        accel_r = temp_rocket.accel_r()
        accel_theta = temp_rocket.accel_theta()
        a2 = accel_r * self.dt
        b2 = (rocket.v_r + a1/2) * self.dt
        c2 = accel_theta/rocket.r * self.dt
        d2 = (rocket.dtheta + d1/2) * self.dt

        temp_rocket = copy.deepcopy(rocket)
        temp_rocket.v_r += a2/2
        temp_rocket.r += b2/2
        temp_rocket.dtheta += c2/2
        temp_rocket.theta += d2/2
        temp_rocket.t += self.dt/2
        accel_r = temp_rocket.accel_r()
        accel_theta = temp_rocket.accel_theta()
        a3 = accel_r * self.dt
        b3 = (rocket.v_r + a2/2) * self.dt
        c3 = accel_theta/rocket.r * self.dt
        d3 = (rocket.dtheta + d2/2) * self.dt

        temp_rocket = copy.deepcopy(rocket)
        temp_rocket.v_r += a3
        temp_rocket.r += b3
        temp_rocket.dtheta += c3
        temp_rocket.theta += d3
        temp_rocket.t += self.dt
        accel_r = temp_rocket.accel_r()
        accel_theta = temp_rocket.accel_theta()
        a4 = accel_r * self.dt
        b4 = (rocket.v_r + a3) * self.dt
        c4 = accel_theta/rocket.r * self.dt
        d4 = (rocket.dtheta + d3) * self.dt

        rocket.v_r += (a1 + 2*a2 + 2*a3 + a4)/6
        rocket.r += (b1 + 2*b2 + 2*b3 + b4)/6
        rocket.dtheta += (c1 + 2*c2 + 2*c3 + c4)/6
        rocket.theta += (d1 + 2*d2 + 2*d3 + d4)/6
        rocket.t += self.dt
         
class Simulation: 
    def reset(self, rocket, obs):
        self.rocket = rocket
        self.obs = obs
    
    def __init__(self, rocket, obs):
        self.reset(rocket, obs)
    
    def run_plot(self):
        r_list = self.obs.r
        theta_list = self.obs.theta
        print("Maximum r:", max(r_list))
        print(self.rocket.status)

        swap_list = self.obs.stageswap_coord

        x_list = r_list * np.cos(theta_list)
        y_list = r_list * np.sin(theta_list)

        x_min, x_max = min(x_list) - 1e6, max(x_list) + 1e6
        y_min, y_max = min(y_list) - 1e6, max(y_list) + 1e6

        # Stageswap
        x_stage = [coord[0]*np.cos(coord[1]) for coord in swap_list]
        y_stage = [coord[0]*np.sin(coord[1]) for coord in swap_list]

        # Create earth and outer atmosphere
        theta_space = np.linspace(0, 2*np.pi, 100)
        r_earth = R
        x_earth = r_earth*np.cos(theta_space)
        y_earth = r_earth*np.sin(theta_space)

        r_atmo = R + 1e5
        x_atmo = r_atmo*np.cos(theta_space)
        y_atmo = r_atmo*np.sin(theta_space)

        mask = (x_earth >= x_min) & (x_earth <= x_max) & (y_earth >= y_min) & (y_earth <= y_max)

        plt.figure(figsize=(8, 8))

        plt.plot(x_list, y_list, 'r-', label="Rocket Trajectory")
        plt.scatter(x_stage, y_stage, label="Change of Stage")
        plt.plot(x_atmo, y_atmo, color='lightblue', label='Atmosphere')
        plt.plot(x_earth[mask], y_earth[mask], color='black', label='Earth Surface')

        plt.axis('equal')
        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)
        plt.legend()
        plt.show()

    def run_plot_obs(self):
        energy_list = self.obs.energy
        t_list = self.obs.t

        stageswap_time = self.obs.stageswap_time
        stageswap_energy = self.obs.stageswap_energy

        print(self.obs.stageswap_time)

        plt.figure()
        plt.plot(t_list, energy_list, label="Total Energy")
        #plt.scatter(stageswap_time, stageswap_energy, color="black", label="Change of Stage")
        plt.legend()
        plt.show()

    # Takes too long time
    def run_animate(self,
                    tmax=30,
                    stepsPerFrame=1,
                    ):


        r_list = self.obs.r
       
        theta_list = self.obs.theta

        x_list = r_list * np.cos(theta_list)
        y_list = r_list * np.sin(theta_list)

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.set_xlim(-max(r_list) - 1, max(r_list) + 1)
        ax.set_ylim(-max(r_list) - 1, max(r_list) + 1)
        ax.set_aspect('equal')
        ax.legend()

        rocket_line, = ax.plot([], [], 'r-', lw=2, label='Trajectory')
        rocket_dot, = ax.plot([], [], 'bo', markersize=8, label='Rocket')

        def init():
            rocket_line.set_data([], [])
            rocket_dot.set_data([], [])
            return rocket_line, rocket_dot

        def update(frame):
            idx = frame * stepsPerFrame
            if idx >=len(x_list):
                idx = len(x_list) - 1
            
            rocket_line.set_data(x_list[:idx+1], y_list[:idx+1])
            rocket_dot.set_data(x_list[:idx], y_list[:idx])
            
            return rocket_line, rocket_dot
        
        time_steps = len(x_list) // stepsPerFrame
        ani = FuncAnimation(fig, update, frames=time_steps, init_func=init, blit=True, interval=50)
        plt.show()
def main():
    _dt = 1
    #integrator = EulerCromer(_dt)
    integrator = RK4(_dt)
    rocket = Rocket(rocket_mass=0.05, mTot=5e5, m1=0.7, t1=200, ve1=6000, m2=0.25, t2=300, ve2=6000)
    obs = Observables()
    integrator.simulate(rocket, obs, nn=100)
    print(rocket.status)

    sim = Simulation(rocket, obs)
    sim.run_plot()
    #sim.run_plot_obs()
    #sim.run_animate()

def main4():
    dt_list =[1, 0.1, 0.01, 0.01] 
    nn_list = [1000, 10000, 100000, 1000000]
    x_list = [] 
    y_list = []
    for i in range(len(dt_list)):
        #integrator = RK4(dt_list[i])    
        integrator = RK4(dt_list[i])    
        rocket = Rocket(rocket_mass=0.05, mTot=5e5, m1=0.7, t1=200, ve1=6000, m2=0.25, t2=300, ve2=6000)
        obs = Observables()
        integrator.simulate(rocket, obs, nn=nn_list[i])
        y_list.append(obs.r)
        x_list.append(obs.t)
    
    plt.figure()
    for j in range(len(y_list)):
        plt.plot(x_list[j], y_list[j], label="dt="+str(dt_list[j]))
    plt.legend()
    plt.show()

def main5():
    true_list = []
    time_list = []
    error_list = []

    integrator = RK4(0.0001)    
    rocket = Rocket(rocket_mass=0.05, mTot=5e5, m1=0.7, t1=200, ve1=6000, m2=0.25, t2=300, ve2=6000)
    obs = Observables()
    integrator.simulate(rocket, obs, nn=1000000)
    true_list = obs.energy
    time_list = obs.t
    
    dt_list = [0.1, 0.01, 0.001]
    nn_list = [1000, 10000, 100000]
    #dt_list = [1, 0.1, 0.01]
    #nn_list = [100000, 100000, 100000]
    x_list = [] 
    y_list = []
    for i in range(len(dt_list)):
        #integrator = RK4(dt_list[i])    
        integrator = EulerCromer(dt_list[i])    
        rocket = Rocket(rocket_mass=0.05, mTot=5e5, m1=0.7, t1=200, ve1=6000, m2=0.25, t2=300, ve2=6000)
        obs = Observables()
        integrator.simulate(rocket, obs, nn=nn_list[i])
        y_list.append(obs.energy)
        x_list.append(obs.t)
    print("1")
    plt.figure()
    y_inter = []
    for j in range(len(y_list)):
        y_inter.append(np.interp(time_list, x_list[j], y_list[j]))
        error_list.append(y_inter[-1]-true_list)
    print("2")
    for k in range(len(error_list)):
        plt.plot(time_list, error_list[k], label="dt="+str(dt_list[k]))
    plt.legend()
    plt.show()

def main6():
    true_list = []
    time_list = []
    error_list = []

    integrator = RK4(0.0001)    
    rocket = Rocket(rocket_mass=0.05, mTot=5e5, m1=0.7, t1=200, ve1=6000, m2=0.25, t2=300, ve2=6000)
    obs = Observables()
    integrator.simulate(rocket, obs, nn=100000)
    true_list = obs.r
    time_list = obs.t
    
    dt_list = [0.01, 0.001, 0.0001]
    nn_list = [1000, 10000, 100000]
    #dt_list = [1, 0.1, 0.01]
    #nn_list = [100000, 100000, 100000]
    x_list = [] 
    y_list = []
    for i in range(len(dt_list)):
        #integrator = RK4(dt_list[i])    
        integrator = RK4(dt_list[i])    
        rocket = Rocket(rocket_mass=0.05, mTot=5e5, m1=0.7, t1=200, ve1=6000, m2=0.25, t2=300, ve2=6000)
        obs = Observables()
        integrator.simulate(rocket, obs, nn=nn_list[i])
        y_list.append(obs.r)
        x_list.append(obs.t)
        
    print("1")
    plt.figure()
    y_inter = []
    for j in range(len(y_list)):
        y_inter.append(np.interp(time_list, x_list[j], y_list[j]))
        error_list.append(y_inter[-1]-true_list)
    print("2")
    for k in range(len(error_list)):
        plt.plot(time_list, error_list[k], label="dt="+str(dt_list[k]))
    plt.legend()
    plt.show()

main4()
    

