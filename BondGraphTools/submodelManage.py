import json
import BondGraphTools as bgt

from BondGraphTools import new, add, connect, expose
from BondGraphTools.exceptions import InvalidComponentException




### 支持嵌套的子模块读取
class CompositeBuilder:
    def __init__(self, comp_def: dict,default_value=None):
        self.comp_def = comp_def
        self.name = comp_def.get("name")
        self.submodel = new(name=self.name)
        self.components = {}  # 名称到组件的映射
        self.default_value=default_value
        self._build()

    def _create_component(self, name: str, comp_spec: dict):
        """创建组件（支持基本类型和复合类型）"""
        comp_type = comp_spec["type"]
        lib, comp = comp_type.rsplit('.', 1)
        
        # 1. 基本组件
        if lib in ['base','BioChem','elec','base1']:
            if self.default_value:
                value = comp_spec.get("value", None)
            else:
                value = None
            return new(comp, value=value, name=name,library=lib)
        
        # 2. 复合组件（递归构建）

        try:
            from BondGraphTools.actions import new_extended
            user_params = comp_spec.get("value", {})
            return new_extended(comp_type,value=user_params, name=name)
        
        # 3. 未知组件类型
        except Exception as e:
            raise InvalidComponentException(f"无法创建组件 {name} ({comp_type}): {str(e)}")
        #     raise ValueError(f"未知组件类型: {comp_type}")

    def _build(self):
        desc = self.comp_def
        subcomponents = desc.get("subcomponents", {})
        connections = desc.get("connections", [])
        exposed = desc.get("expose_ports", {})

        # Step 1: 创建所有子组件
        for name, comp_spec in subcomponents.items():
            component = self._create_component(name, comp_spec)
            self.components[name] = component
            add(self.submodel, component)

        # Step 2: 建立连接（只需处理当前层组件）
        for conn in connections:
            if len(conn) == 2:
                src, dest = conn
                connect(self.components[src], self.components[dest])
            else:
                raise ValueError(f"无效连接格式: {conn}")

        # Step 3: 暴露端口（只需处理当前层组件）
        for comp_name, port_info in exposed.items():
            component = self.components[comp_name]
            port_label = port_info["port"]        # 从字典中获取 "port" 的值
            direction = port_info["direction"] 
            expose(component, label=port_label,direction=direction)

    def get_model(self):
        return self.submodel 
    
# 用法示例
if __name__ == "__main__":
    # 从文件读配置（可根据需要替换为你的json路径）
    # from joblib import Parallel, delayed, cpu_count

    import os
    with open(os.path.join(os.path.dirname(__file__),"components", "elecComp.json"), "r") as f:
        data = json.load(f)

    from matplotlib import pyplot as plt

    comp_lib=data
    # 构建嵌套复合模块

    ### 测试复合模块1
    # builder = CompositeBuilder(
    #     comp_def=comp_lib["components"]["RCBlock"],
    #     comp_lib=comp_lib,
    #     # name="MyDoubleRC"
    # )
    # #test RC 模型
    # double_rc_model = builder.get_model()
 
    # C2=new("C", value=150*1e-6)
    # mainmodel = new(name='RC')
    # add(mainmodel, double_rc_model, C2)
    # connect(C2,double_rc_model.get_port('P'))

    # mainmodel.state_vars
    # timespan = [0, 50]
    # x0 = {'x_0':1, 'x_1':0}  # 初始状态变量
    # mainmodel.constitutive_relations
    # t, x = bgt.simulate(mainmodel, timespan=timespan, x0=x0)
    # import matplotlib.pyplot as plt

    # plt.plot(t,x[:,0], '-b', label='q_C1')
    # plt.plot(t,x[:,1], '-r', label='q_C2')
    # plt.xlabel("time (s)")
    # plt.ylabel("electric charge (Coulomb)")
    # plt.legend(loc='upper right')
    # plt.grid()
    # plt.show()

    ### 测试嵌套复合模块1   ## 这个测试DoubleRC有问题
    builder = CompositeBuilder(
        comp_def=comp_lib["components"]["DoubleRC1"]
    )
    #test RC 模型
    double_rc_model = builder.get_model()
    # print(double_rc_model.components[0].components[1].params)
 
    C2=new("C", value=100*1e-6)
    mainmodel = new(name='RC')
    add(mainmodel, double_rc_model, C2)
    connect(double_rc_model.get_port('P'),C2)

    print(mainmodel.state_vars)
    timespan = [0, 50]
    x0 = {'x_0':1, 'x_1':0, 'x_2':0}  # 初始状态变量
    print(mainmodel.constitutive_relations)
    t, x = bgt.simulate(mainmodel, timespan=timespan, x0=x0)
    
    import matplotlib.pyplot as plt

    plt.plot(t,x[:,0], '-b', label='q_C1')
    plt.plot(t,x[:,1], '-r', label='q_C2')
    plt.plot(t,x[:,2], '-g', label='q_C3')
    plt.xlabel("time (s)")
    plt.ylabel("electric charge (Coulomb)")
    plt.legend(loc='upper right')
    plt.grid()
    plt.show()
