from BondGraphTools import new, draw, simulate
model = new(name='RC')
C = new("C", value=1.0)
R = new("R", value=1.0)
se = new("Se", value=1.0)
one = new("1")
from BondGraphTools import add, connect, expose
add(model, R,C,one,se)
connect(se,one)
connect(R,one)
connect(C,one)
model.state_vars
timespan = [0, 5]
x0 = {'x_0':1}
t, x = simulate(model, timespan=timespan, x0=x0)
import matplotlib.pyplot as plt

plt.plot(t,x)
plt.show()
plt.savefig("RC_2.svg", pad_inches=0, bbox_inches="tight")