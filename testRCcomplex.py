from BondGraphTools import new, draw, simulate
from BondGraphTools import add, connect, expose
submodel = new(name='subR')
R1=new("R",value=1.0)
one=new("1")
C1=new("C",value=1.0)
zero = new("0")  # 新增 0-结
add(submodel, R1,one,C1)
connect(R1, zero)
connect(zero, one) # 0-结 -> 1-结
connect( C1,one)
expose(zero,label='P')

Se=new("Se",value=1.0)
mainmodel = new(name='RC')
add(mainmodel, submodel, Se)
connect(Se,submodel.get_port("P"))

mainmodel.state_vars
timespan = [0, 5]
x0 = {'x_0':1}
t, x = simulate(mainmodel, timespan=timespan, x0=x0)
import matplotlib.pyplot as plt

plt.plot(t,x)
plt.show()
plt.savefig("RC_2.svg", pad_inches=0, bbox_inches="tight")
