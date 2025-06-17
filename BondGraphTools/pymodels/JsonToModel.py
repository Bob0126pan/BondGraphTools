import json
import BondGraphTools as bgt

from BondGraphTools import new, add, connect, expose
from BondGraphTools.exceptions import InvalidComponentException

# class CompositeBuilder:
#     def __init__(self, comp_def: dict, name: str = None):
#         self.comp_def = comp_def
#         self.name = name or comp_def.get("id", "composite_model")
#         self.submodel = new(name=self.name)
#         self.components = {}  # 映射名称到 BondGraphTools 元件对象
#         self._build()

#     def _build(self):
#         desc = self.comp_def
#         subcomponents = desc.get("subcomponents", {})
#         connections = desc.get("connections", [])
#         exposed = desc.get("expose_ports", {})

#         # Step 1: 创建所有子组件
#         for name, comp in subcomponents.items():
#             typ = comp["type"]
#             value = comp.get("value", None)
#             if value is not None:
#                 element = new(typ, value=value)
#             else:
#                 element = new(typ)
#             self.components[name] = element
#             add(self.submodel, element)

#         # Step 2: 建立连接
#         for conn in connections:
#             if len(conn) == 2:
#                 connect(self.components[conn[0]], self.components[conn[1]])
#             else:
#                 raise ValueError(f"Invalid connection format: {conn}")

#         # Step 3: 暴露端口
#         for internal_name, port_label in exposed.items():
#             expose(self.components[internal_name], label=port_label)

#     def get_model(self):
#         return self.submodel



### 支持嵌套的子模块读取
class CompositeBuilder:
    def __init__(self, comp_def: dict, comp_lib: dict, name: str = None):
        self.comp_def = comp_def
        self.comp_lib = comp_lib
        self.name = name or comp_def.get("id", "composite_model")
        self.submodel = new(name=self.name)
        self.components = {}  # 名称到组件的映射
        self._build()

    def _create_component(self, name: str, comp_spec: dict):
        """创建组件（支持基本类型和复合类型）"""
        comp_type = comp_spec["type"]
        
        # 1. 基本组件
        if comp_type in ["R", "C", "I", "Se", "Sf", "TF", "GY", "0", "1", "SS"]:
            value = comp_spec.get("value", None)
            return new(comp_type, value=value, name=name)
        
        # 2. 复合组件（递归构建）
        elif comp_type in self.comp_lib.get("components", {}):
            composite_def = self.comp_lib["components"][comp_type]
            sub_builder = CompositeBuilder(composite_def, self.comp_lib, name)
            return sub_builder.get_model()
        
        # 3. 未知组件类型
        else:
            raise ValueError(f"未知组件类型: {comp_type}")

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
        for comp_name, port_label in exposed.items():
            component = self.components[comp_name]
            expose(component, label=port_label)

    def get_model(self):
        return self.submodel 
    
# 用法示例
if __name__ == "__main__":
    # 从文件读配置（可根据需要替换为你的json路径）
    import os
    with open(os.path.join(os.path.dirname(__file__), "RC.json"), "r") as f:
        data = json.load(f)

    from matplotlib import pyplot as plt

    comp_lib=data
    # 构建嵌套复合模块
    builder = CompositeBuilder(
        comp_def=comp_lib["components"]["RCBlock"],
        comp_lib=comp_lib,
        name="MyDoubleRC"
    )
    double_rc_model = builder.get_model()

    # testsimple rc
    Se=new("Se",value=1.0)
    mainmodel = new(name='RC')
    add(mainmodel, double_rc_model, Se)
    connect(Se,double_rc_model.get_port('P'))

    mainmodel.state_vars
    timespan = [0, 5]
    x0 = {'x_0':1}
    t, x = bgt.simulate(mainmodel, timespan=timespan, x0=x0)
    import matplotlib.pyplot as plt

    plt.plot(t,x)
    plt.show()
    plt.savefig("RC_2.svg", pad_inches=0, bbox_inches="tight")

    # # 构建主模型，并连接一个 Se 元件
    # mainmodel = new(name="MainModel")
    # Se = new("Se", value=1.0)  # 恒压源
    # load = new("R", value=1.0, name="load_resistor")
    # one = new("0")  # 连接点
    # add(mainmodel, Se, double_rc_model, load, one)

    # connect(one, double_rc_model.get_port('Input'))  # 自动连接暴露端口 Input
    # connect(one,double_rc_model.get_port('Output'))  # 自动连接暴露端口 Output
    # connect(one, Se)  # 连接恒压源到输入端口
    # connect(one, load)  # 连接输出到负载电阻

    # # 初始化、仿真
    # x0 = {'x_0': 1.0, 'x_1': 0.5}  # 2 个电容器初值
    # tspan = [0, 5]

    # t, x = bgt.simulate(mainmodel, timespan=tspan, x0=x0)

    # # 绘图
    # plt.plot(t, x[:, 0], label='C1 voltage')
    # plt.plot(t, x[:, 1], label='C2 voltage')
    # plt.xlabel("Time")
    # plt.ylabel("State variables")
    # plt.title("Double RC Block Simulation")
    # plt.legend()
    # plt.grid()
    # plt.tight_layout()
    # plt.savefig("double_rc_sim.svg")
    # plt.show()
