import BondGraphTools as bgt
model=bgt.new(name='circuit_1')
# Parameters' values
C1_value=100*1e-6     #(100 uF)
C2_value=150*1e-6     #(150 uF)
R_value=100*1e3       #(100 k)
C1=bgt.new("C", value=C1_value)
C2=bgt.new("C", value=C2_value)
R=bgt.new("R", value=R_value)
one_junc=bgt.new("1")
bgt.add(model,C1,C2,R,one_junc)
bgt.connect(C1,one_junc)
bgt.connect(one_junc,R)
bgt.connect(one_junc,C2)
# bgt.draw(model)
timespan=[0,50]
model.state_vars
x0={"x_0":1, "x_1":0}
model.constitutive_relations
t, x = bgt.simulate(model, timespan=timespan, x0=x0)
import matplotlib.pyplot as plt
plt.plot(t,x[:,0], '-b', label='q_C1')
plt.plot(t,x[:,1], '-r', label='q_C2')
plt.xlabel("time (s)")
plt.ylabel("electric charge (Coulomb)")
plt.legend(loc='upper right')
plt.grid()
plt.show()
import numpy as np
f = np.array(x[:,0], dtype=float)
slope=np.gradient(f,0.1)
v_C1=slope

# dq_C2/dt = v_C2 (flow in C2)
import numpy as np
f = np.array(x[:,1], dtype=float)
slope=np.gradient(f,0.1)
v_C2=slope

plt.plot(t,v_C1, '-b', label='v_C1')
plt.plot(t,v_C2, '-r', label='v_C2')
plt.xlabel("time (s)")
plt.ylabel("Flow (Coulomb/s)")
plt.legend(loc='upper right')
plt.grid()
plt.show()