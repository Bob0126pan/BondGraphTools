def test_tube_buid_eqs():
    import BondGraphTools as bgt
    from BondGraphTools import expose, connect, add,draw
    model=bgt.new(name='straight tube')
    model.para_symbols=True

    C=bgt.new("C", name="C",library="base")     #(m6/J)

    # The amounts R-elements are assumed to be equal in a straight tube
    R1=bgt.new("R", name="R1",library="base")     #(J.s/m6)
    R2=bgt.new("R", name="R2",library="base")     #(J.s/m6)

    # The amounts of the I-elements are assumed to be equal in a straight tube
    L1=bgt.new("I", name="L1",library="base")     #(J.s2/m6)
    L2=bgt.new("I", name="L2",library="base")     #(J.s2/m6)

    zero_junc=bgt.new("0")
    one_junc_1=bgt.new("1")
    one_junc_2=bgt.new("1")

    bgt.add(model,C,R1,R2,L1,L2,zero_junc,one_junc_1,one_junc_2)

    bgt.connect(one_junc_1,R1)
    bgt.connect(one_junc_1,L1)
    bgt.connect(one_junc_1,zero_junc)
    bgt.connect(zero_junc,one_junc_2)
    bgt.connect(zero_junc,C)
    bgt.connect(one_junc_2,R2)
    bgt.connect(one_junc_2,L2)
    # model.map_port('P1', ((one_junc_1, 'e'),(one_junc_1, 'f')))
    # model.map_port('P2', ((one_junc_2, 'e'),(one_junc_2, 'f')))
    expose(one_junc_1)
    expose(one_junc_2)
    print(model.constitutive_relations)
    #  前半段采用para_symbols  True 求得带参数方程
    # 后半段采用False 求不带参数方程，两个综合起来确定对应的方程。
    model.para_symbols=False
    # model.map_port('P1', ((one_junc_1, 'e'),(one_junc_1, 'f')))
    # model.map_port('P2', ((one_junc_2, 'e'),(one_junc_2, 'f')))
    expose(one_junc_1)
    expose(one_junc_2)
    print(model.constitutive_relations)
    pass