import numpy as np
from scipy.interpolate import interp1d
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
import scipy

class DelayVar:
    """
    The instances of this class are special function-like
    variables which store their past values in an interpolator and
    can be called for any past time: Y(t), Y(t-d).
    Very convenient for the integration of DDEs.
    """

    def __init__(self, g, tc=0):
        """ g(t) = expression of Y(t) for t<tc """

        self.g = g
        self.tc = tc
        # We must fill the interpolator with 2 points minimum

        self.interpolator = scipy.interpolate.interp1d(
            np.array([tc - 1, tc]),  # X
            np.array([self.g(tc), self.g(tc)]).T,  # Y
            kind="linear",
            bounds_error=False,
            fill_value=self.g(tc)
        )

    def update(self, t, Y):
        """ Add one new (ti,yi) to the interpolator """
        Y2 = Y if (Y.size == 1) else np.array([Y]).T
        self.interpolator = scipy.interpolate.interp1d(
            np.hstack([self.interpolator.x, [t]]),  # X
            np.hstack([self.interpolator.y, Y2]),  # Y
            kind="linear",
            bounds_error=False,
            fill_value=Y
        )

    def __call__(self, t=0):
        """ Y(t) will return the instance's value at time t """
        return self.g(t) if (t <= self.tc) else self.interpolator(t)


class DelayedDynamicSolver:
    def __init__(self, func, x0, t_span, dt, delays, difEqNum=None,method="newton", **kwargs):
        """
        Initialize the solver.

        Parameters:
            func: Callable
                The function f(x(t), x_delayed, t), representing the dynamic system.
            x0: numpy array
                Initial state vector.
            t_span: tuple
                Time span as (t_start, t_end).
            dt: float
                Time step size.
            delays: list of floats
                List of time delays for each variable.
            method: str
                Numerical method to use: "newton" or "runge_kutta".
            kwargs: dict
                Additional arguments for specific methods.
        """
        self._func = func
        self.x0 = x0
        self.num_d=len(self.x0) if difEqNum is None else difEqNum
        self.t_span = t_span
        self.dt = dt
        self.delays = delays
        self.method = method.lower()
        self.tol = kwargs.get("tol", 1e-6)
        self.max_iter = kwargs.get("max_iter", 100)
        self.perturbation = kwargs.get("perturbation", 1e-8)
        self.jacobian_mode = kwargs.get("jacobian_mode", "full_numerical")

        # Initialize time and history buffers
        self.t_values = np.arange(t_span[0], t_span[1] + dt, dt)
        self.max_delay = max(delays)
        self.history = [(t_span[0] - i * dt, x0) for i in range(int(np.ceil(self.max_delay / dt)) + 1)][::-1]
        # self.delay_vars = [DelayVar(g[i], tc=0) for i in range(len(delays))]

    def get_delayed_state(self, current_time):
        """
        Interpolates to get the delayed states for all variables.
        """
        times = np.array([h[0] for h in self.history])
        states = np.array([h[1] for h in self.history])
        interpolator = interp1d(times, states.T, kind="linear", fill_value="extrapolate", axis=1)
        return np.array([states[0] if delay<=current_time else interpolator(current_time - delay) for index,delay in enumerate(self.delays)]).T 

    def get_residual(self,x_next,x_prev,x_delayed,t):
        return np.hstack([
                    x_next[:self.num_d] - x_prev[:self.num_d] - self.dt * self.func(x_next, x_delayed, t)[:self.num_d],  # Differential part
                    self.func(x_next, x_delayed, t)[self.num_d:]  # Algebraic part
                ])
    
    def diffunc(self,x,x_delayed,t):
        f_d =np.array([
        (-43.2*x[0]**2 +4.2- x[2])/0.02535,
        (-195.37*x[1]**2 +4.2- x[2])/0.02535,
        -360*x[2]+(x_delayed[0,0]+x_delayed[1,1])*3280
        ])
        return f_d

    def algebraic_eqns(self,x,x_delayed,t):
        tau = 1 #这个比较重要
        f_a = np.array([
            (x[3]**2 - 3 * x[0]) / tau
        ])
        return f_a

    def func(self,x,x_delayed,t):
        difeqns=self.diffunc(x,x_delayed,t)
        alg_eqns=self.algebraic_eqns(x,x_delayed,t)
        return np.hstack([difeqns,alg_eqns])

    def compute_jacobian(self, x_next, x_prev,residual, x_delayed, t):
        """
        Compute the Jacobian matrix based on the selected mode.
        """
        n = len(x_next)
        if self.jacobian_mode == "full_numerical":
            jacobian = np.zeros((n, n))
            for i in range(n):
                x_perturbed = np.copy(x_next)
                x_perturbed[i] += self.perturbation
                perturbed_residual = self.get_residual(x_perturbed,x_prev,x_delayed,t)
                jacobian[:, i] = (perturbed_residual.flatten() - residual.flatten()) / self.perturbation
        elif self.jacobian_mode == "sparse_numerical":
            jacobian = lil_matrix((n, n))
            for i in range(n):
                x_perturbed = np.copy(x_next)
                x_perturbed[i] += self.perturbation
                perturbed_residual = self.get_residual(x_perturbed,x_prev,x_delayed,t)
                jacobian[:, i] = (perturbed_residual.flatten() - residual.flatten()) / self.perturbation
            jacobian = jacobian.tocsr()
        else:
            raise ValueError(f"Invalid Jacobian mode: {self.jacobian_mode}")
        return jacobian

    def solve_newton(self):
        """
        Solve using the Newton-Raphson method.
        """
        x_values = [self.x0]
        for t in self.t_values[:-1]:
            x_prev = x_values[-1]
            x_next = np.copy(x_prev)  # Initial guess for x_next
            x_delayed = self.get_delayed_state(t + self.dt)

            for _ in range(self.max_iter):
                residual =  self.get_residual(x_next,x_prev,x_delayed,t+self.dt)
                if np.linalg.norm(residual, ord=2) < self.tol:
                    break
                jacobian = self.compute_jacobian(x_next,x_prev, residual, x_delayed, t + self.dt)
                delta_x = spsolve(jacobian, -residual) if self.jacobian_mode == "sparse_numerical" else np.linalg.solve(jacobian, -residual)
                x_next += delta_x
            else:
                raise RuntimeError(f"Newton-Raphson failed to converge at time {t + self.dt}.")
            self.history.append((t + self.dt, x_next))
            self.history = self.history[-int(np.ceil(self.max_delay / self.dt)) - 1:]
            x_values.append(x_next)
            # print(x_next[3]/x_next[1])
        return np.array(self.t_values), np.array(x_values)

    def solve_runge_kutta(self):
        """
        Solve using the Runge-Kutta method (RK4).
        """
        x_values = [self.x0]
        for t in self.t_values[:-1]:
            x_prev = x_values[-1]
            x_delayed = self.get_delayed_state(t)
            k1 = self.func(x_prev, x_delayed, t)
            k2 = self.func(x_prev + 0.5 * self.dt * k1, x_delayed, t + 0.5 * self.dt)
            k3 = self.func(x_prev + 0.5 * self.dt * k2, x_delayed, t + 0.5 * self.dt)
            k4 = self.func(x_prev + self.dt * k3, x_delayed, t + self.dt)
            x_next = x_prev + (self.dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
            x_next[self.num_d:]=k1[self.num_d:]

            self.history.append((t + self.dt, x_next))
            self.history = self.history[-int(np.ceil(self.max_delay / self.dt)) - 1:]
            x_values.append(x_next)
        return np.array(self.t_values), np.array(x_values)

    def solve(self):
        """
        Main interface to solve the system using the specified method.
        """
        if self.method == "newton":
            return self.solve_newton()
        elif self.method == "runge_kutta":
            return self.solve_runge_kutta()
        else:
            raise ValueError(f"Invalid method: {self.method}")

# Example usage
def delayed_dynamic_system(x, x_delayed, t):
    # a, b, c, d = 1.0, 0.5, 0.3, -0.7
    # Differential equations
    f_d =np.array([
        (-43.2*x[0]**2 +4.2- x[2])/0.02535,
        (-195.37*x[1]**2 +4.2- x[2])/0.02535,
        -360*x[2]+(x_delayed[0,0]+x_delayed[1,1])*3280
    ])
     # Algebraic equations
    tau = 1 #这个比较重要
    f_a = np.array([
        (x[3]**2 - 3 * x[0]) / tau
    ])
    #替换为
    # f_a=np.array([
    #     3 * x[1]
    # ])

    return np.hstack([f_d, f_a])

x0 = np.array([0.1, 0.1,0.1,0.3])
t_span = (0, 0.1)
dt = 0.001
delays = [0.005, 0.005]  # Different delays for x1 and x2

solver = DelayedDynamicSolver(delayed_dynamic_system, x0, t_span, dt, delays, difEqNum=3,method="newton")
t_values, x_values = solver.solve()

plt.plot(t_values, x_values[:, 0], label="mass1_t")
plt.plot(t_values, x_values[:, 1], label="mass2_t")
plt.plot(t_values, x_values[:, 2], label="Pressure_t")
plt.plot(t_values, x_values[:, 3], label="3_mass2")
plt.xlabel("Time")
plt.ylabel("State")
plt.grid(True)
plt.legend()
plt.title("Delayed Dynamic System with Multiple Delays")
plt.show()
