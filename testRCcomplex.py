from BondGraphTools import new, draw, simulate
from BondGraphTools import add, connect, expose
from BondGraphTools.base import Port
submodel = new(name='subR')
R1=new("R",value=1.0)
one=new("1")
C1=new("C",value=1.0)
ss=new("SS")
add(submodel, R1,one,C1,ss)
connect(one,R1)
connect( one,C1)
connect( one,ss)
expose(ss)

Se=new("Se",value=1.0)
mainmodel = new(name='RC')
add(mainmodel, submodel, Se)
# P=submodel.get_port("P")
connect(Se,submodel)

mainmodel.state_vars
timespan = [0, 5]
x0 = {'x_0':1}
t, x = simulate(mainmodel, timespan=timespan, x0=x0)
import matplotlib.pyplot as plt

plt.plot(t,x)
plt.show()
plt.savefig("RC_2.svg", pad_inches=0, bbox_inches="tight")
