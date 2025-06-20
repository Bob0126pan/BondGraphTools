def testRCBasic():
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

def testRCcomplex():
    from BondGraphTools import new, draw, simulate
    from BondGraphTools import add, connect, expose
    from BondGraphTools.base import Port
    submodel = new(name='subR')
    R1=new("R",value=1.0)
    one=new("1")
    C1=new("C",value=1.0)
    add(submodel, R1,one,C1)
    connect(one,R1)
    connect( one,C1)
    expose(one,'P')

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


def testRCSubmodel():
    from BondGraphTools.actions import new,new_extended
    from BondGraphTools import new, draw, simulate
    from BondGraphTools import add, connect, expose
    # 创建基本电容
    C0 = new_extended("base.C", value=1e-5)
    C1 = new_extended("base.C", value={'C':1e-6})
    C2 = new_extended("base.C", value={'C':1e-3})
    C3 = new_extended("BioChem.Ce",value={'k':1, 'R':1, 'T':1})
    # 创建复合模块 RCBlock（通过lib）
    double_rc_model = new_extended("elecSubmodel.RCBlock", value={"R.r": 100000, "C1.C": 0.00015}, name="my_rc")
    print(double_rc_model.components[1].params)

    C2=new("C", value=100*1e-6)
    mainmodel = new(name='RC')
    add(mainmodel, double_rc_model, C2)
    connect(double_rc_model.get_port('P'),C2)

    print(mainmodel.state_vars)
    timespan = [0, 50]
    x0 = {'x_0':1, 'x_1':0}  # 初始状态变量
    print(mainmodel.constitutive_relations)
    t, x = simulate(mainmodel, timespan=timespan, x0=x0)
    
    import matplotlib.pyplot as plt

    plt.plot(t,x[:,0], '-b', label='q_C1')
    plt.plot(t,x[:,1], '-r', label='q_C2')
    # plt.plot(t,x[:,2], '-g', label='q_C3')
    plt.xlabel("time (s)")
    plt.ylabel("electric charge (Coulomb)")
    plt.legend(loc='upper right')
    plt.grid()
    plt.show()

def testRCSubSubmodel():
    from BondGraphTools.actions import new,new_extended
    from BondGraphTools import new, draw, simulate
    from BondGraphTools import add, connect, expose
    # 创建基本电容
    # 创建复合模块 RCBlock（通过lib）
    double_rc_model = new_extended("elecSubmodel.DoubleRC1", value={"rc1.R.r": 100000, "rc1.C1.C": 0.00015,"rc2.R.r": 10000, "rc2.C2.C": 0.00022}, name="my_rc")
    print(double_rc_model.components[1].params)

    C2=new("C", value=100*1e-6)
    mainmodel = new(name='RC')
    add(mainmodel, double_rc_model, C2)
    connect(double_rc_model.get_port('P'),C2)

    print(mainmodel.state_vars)
    timespan = [0, 50]
    x0 = {'x_0':1, 'x_1':0,'x_2':0}  # 初始状态变量
    print(mainmodel.constitutive_relations)
    t, x = simulate(mainmodel, timespan=timespan, x0=x0)
    
    import matplotlib.pyplot as plt

    plt.plot(t,x[:,0], '-b', label='q_C1')
    plt.plot(t,x[:,1], '-r', label='q_C2')
    plt.plot(t,x[:,2], '-g', label='q_C3')
    plt.xlabel("time (s)")
    plt.ylabel("electric charge (Coulomb)")
    plt.legend(loc='upper right')
    plt.grid()
    plt.show()


testRCSubSubmodel()