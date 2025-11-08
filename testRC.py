from BondGraphTools.modelpostPrcessor import BondGraphPost

def testRCsimple():
    from BondGraphTools import new, draw, simulate
    from BondGraphTools import add, connect
    from sympy import init_printing, SparseMatrix, Eq, Symbol, lambdify
    model = new(name='R_Se')
    R = new("R", value=1.0)
    one = new("1")
    se = new("Se", value=1.0)
    I = new("I", value=1.0)
    add(model, R, one, I, se)
    connect(se, one)
    connect(R, one)
    connect(I, one)
    model.constitutive_relations
    timespan = [0, 5]
    x0 = {'x_0':1}
    t, x = simulate(model, timespan=timespan, x0=x0)

    post = BondGraphPost(model, (t, x))
    
    # 4. 查看组件层次
    post.list_components()

    # 5. 获取子模型数据
    submodel_data = post.get_component('I3')

    post.plot_component('I3')

def testRCBasic():
    from BondGraphTools import new, draw, simulate
    model = new(name='RC')
    C = new("C", value=1.0)
    R = new("R", value=1.0)
    I = new("I", value=1.0)
    se = new("Se", value=1.0)
    one = new("1")
    from BondGraphTools import add, connect, expose
    add(model, R,C,one,se,I)
    connect(se,one)
    model.constitutive_relations
    connect(R,one)
    connect(C,one)
    connect(I,one)
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
    # import matplotlib.pyplot as plt

    # plt.plot(t,x)
    # plt.show()
    # plt.savefig("RC_2.svg", pad_inches=0, bbox_inches="tight")
    post = BondGraphPost(mainmodel, (t, x))

    
    # 4. 查看组件层次
    post.list_components()

    # 5. 获取子模型数据
    submodel_data = post.get_component('subR')

    post.plot_component('subR.C2')


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
    double_rc_model = new_extended("elecSubmodel.RCBlock", value={"R.r": 100000, "C1.C": 0.0001}, name="my_rc")
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

    post = BondGraphPost(mainmodel, (t, x))

    
    # 4. 查看组件层次
    post.list_components()

    # 5. 获取子模型数据
    submodel_data = post.get_component('my_rc')

    post.plot_component('my_rc.C1')
    post.plot_component('my_rc.C8')
    
    # import matplotlib.pyplot as plt

    # plt.plot(t,x[:,0], '-b', label='q_C1')
    # plt.plot(t,x[:,1], '-r', label='q_C2')
    # # plt.plot(t,x[:,2], '-g', label='q_C3')
    # plt.xlabel("time (s)")
    # plt.ylabel("electric charge (Coulomb)")
    # plt.legend(loc='upper right')
    # plt.grid()
    # plt.show()

def testRCSubSubmodel():
    from BondGraphTools.actions import new,new_extended
    from BondGraphTools import new, draw, simulate
    from BondGraphTools import add, connect, expose
    # 创建基本电容
    # 创建复合模块 RCBlock（通过lib）
    double_rc_model = new_extended("elecSubmodel.DoubleRC1", value={"rc1.R.r": 10000, "rc1.C1.C": 0.015,"rc2.R.r": 1000, "rc2.C1.C": 0.0022}, name="my_rc")
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

    post = BondGraphPost(mainmodel, (t, x))

    
    # 4. 查看组件层次
    post.list_components()

    # 5. 获取子模型数据
    submodel_data = post.get_component('my_rc')

    post.plot_component('my_rc.rc2.R')
    post.plot_component('my_rc.rc1.C1')
    
    # import matplotlib.pyplot as plt

    # plt.plot(t,x[:,0], '-b', label='q_C1')
    # plt.plot(t,x[:,1], '-r', label='q_C2')
    # plt.plot(t,x[:,2], '-g', label='q_C3')
    # plt.xlabel("time (s)")
    # plt.ylabel("electric charge (Coulomb)")
    # plt.legend(loc='upper right')
    # plt.grid()
    # plt.show()


def test_Tube():
    from BondGraphTools.actions import new,new_extended
    from BondGraphTools import new, draw, simulate
    from BondGraphTools import add, connect, expose

    tube_model = new_extended("fluidSubmodel.Tube", name="my_tube")
    se1 = new("Se", value=11.997e6)
    se2 = new("Se", value=10.664e6)

    mainmodel = new(name='TubeSystem')
    add(mainmodel, tube_model,se1,se2)
    connect(tube_model.get_port('one1'),se1)
    connect(tube_model.get_port('one2'),se2)

    print(mainmodel.state_vars)
    timespan = [0, 10]
    x0={"x_0":5*1e-6, "x_1":0, "x_2":0} # 初始状态变量
    print(mainmodel.constitutive_relations)
    t, x = simulate(mainmodel, timespan=timespan, x0=x0)

    post = BondGraphPost(mainmodel, (t, x))

    
    # 4. 查看组件层次
    post.list_components()

    # 5. 获取子模型数据
    submodel_data = post.get_component('my_rc')

    post.plot_component('my_tube.C')
    post.plot_component('my_tube.I1')
    
    # import matplotlib.pyplot as plt

    # plt.plot(t,x[:,0], '-b', label='q_C')
    # plt.xlabel("time (s)")
    # plt.ylabel("fluid volume (m^3)")
    # plt.legend(loc='upper right')
    # plt.grid()
    # plt.show()

def Tube():
    import BondGraphTools as bgt
    model=bgt.new(name='straight tube')
    Se1=bgt.new("Se",value=11.997e6)     #(J/m6)
    Se2=bgt.new("Se",value=10.664e6)     #(J/m6)

    C=bgt.new("C", value=0.60015e-6)     #(m6/J)

    # The amounts R-elements are assumed to be equal in a straight tube
    R1=bgt.new("R", value=10.664e-6)     #(J.s/m6)
    R2=bgt.new("R", value=10.664e-6)     #(J.s/m6)

    # The amounts of the I-elements are assumed to be equal in a straight tube
    L1=bgt.new("I", value=0.06665e6)     #(J.s2/m6)
    L2=bgt.new("I", value=0.06665e6)     #(J.s2/m6)

    zero_junc=bgt.new("0")
    one_junc_1=bgt.new("1")
    one_junc_2=bgt.new("1")

    bgt.add(model,Se1,Se2,C,R1,R2,L1,L2,zero_junc,one_junc_1,one_junc_2)
    bgt.connect(one_junc_1,Se1)
    bgt.connect(one_junc_1,R1)
    bgt.connect(one_junc_1,L1)
    bgt.connect(one_junc_1,zero_junc)
    bgt.connect(zero_junc,one_junc_2)
    bgt.connect(zero_junc,C)
    bgt.connect(one_junc_2,R2)
    bgt.connect(one_junc_2,L2)
    bgt.connect(one_junc_2,Se2)
    timespan=[0,5]
    print(model.state_vars)
    x0={"x_0":5*1e-6, "x_1":0, "x_2":0}
    print(model.constitutive_relations)
    t, x = bgt.simulate(model, timespan=timespan, x0=x0)
    import matplotlib.pyplot as plt
    plt.plot(t,x[:,0], '-b', label='q_C')
    plt.xlabel("time (s)")
    plt.ylabel("volume (m3)")  #metre3
    plt.legend(loc='upper right')
    plt.grid()
    plt.show()

test_Tube()
