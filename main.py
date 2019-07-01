import math


class S:  # Needed for storing space and time objects
    def __init__(self, v, x):
        self.vel = v
        # Space
        self.x = x

    def __add__(self, external):
        return S(self.vel + external.vel, self.x + external.x)

    def __sub__(self, external):
        return S(self.vel - external.vel, self.x - external.x)

    def __mul__(self, value):
        return S(self.vel * value, self.x * value)

    def __truediv__(self, value):
        return S(self.vel / value, self.x / value)

    def __radd__(self, external):
        return S(self.vel + external.vel, self.x + external.x)

    def __rsub__(self, external):
        return S(self.vel - external.vel, self.x - external.x)

    def __rmul__(self, value):
        return S(self.vel * value, self.x * value)

    def __rtruediv__(self, value):
        return S(self.vel / value, self.x / value)

    def __str__(self):
        return f"[v={self.vel}, x={self.x}]"

    
def runge_kutta(s_0, dt, F, final_time, initial_time):
    s = s_0
    while initial_time < final_time:
        k1 = F(s, initial_time)
        k2 = F(s + (k1/2), initial_time + (dt/2))
        k3 = F(s + (k2/2), initial_time + (dt/2))
        k4 = F(s + k3, initial_time + dt)
        s = s + dt*(k1/6 + k2/3 + k3/3 + k4/6)
        initial_time += dt
    return s


def main():
    tolerance = 0.00001
    mass = 2
    k = 4
    w = math.sqrt(k/mass)
    zeta = 0.05

    # Function of the force applied to the cart
    force = lambda t: (4*t) if t <= 0.5 else (4*(1-t)) if (t>0.5 and t<1.0) else 0

    initial_velocity = 0.0
    initial_position = 1.0
    s_0 = S(initial_velocity, initial_position)

    final_time = 1.2

    # dx(t)/d(t) = velocity
    # d²x(t)/dt² + 2*zeta*w*dx(t)/dt + w²x(t) = force(t)/mass, isolating the d²x(t)/dt² we have
    F = lambda s_i, t: S(force(t) / mass - 2 * zeta * w * s_i.vel - w * w * s_i.x, s_i.vel)

    # Solving
    count = 0
    dt = 0.5
    s_i_plus_one = runge_kutta(s_0, dt, F, final_time, initial_time=0.0)
    x_old = s_i_plus_one.x
    print('First solution: ', s_i_plus_one)
    while True and count < 10000:
        dt /= 2
        s_i_plus_one = runge_kutta(s_0, dt, F, final_time, initial_time=0.0)
        print('Next solution: ', s_i_plus_one)
        x_new = s_i_plus_one.x
        count += 1
        if abs((x_new - x_old)/x_new) < tolerance:
            break
        x_old = x_new
    print('\n\nDelta t used: ', dt)
    print('Delta t adjustments: ', count)
    print('Solution for 4th degree Runge-Kutta at time 1.2s and tolerance 0.00001: ', s_i_plus_one)


if __name__ == "__main__":
    main()
