import BondGraphTools as bgt
import numpy as np


import BondGraphTools as bgt
from abc import ABC, abstractmethod

class SubModelsBase(ABC):
    def __init__(self, name="Unnamed", params=None):
        self.name = name
        self.system = bgt.new(name=name)
        self.params = self.default_params()
        if params:
            self.params.update(params)
        self.build()

    @abstractmethod
    def default_params(self):
        """返回默认参数的字典"""
        return {}

    @abstractmethod
    def build(self):
        """子类实现建模逻辑"""
        pass

    def get_system(self):
        """获取构建完成的bond graph子系统"""
        return self.system

    def connect_to(self, other, port_self=0, port_other=0):
        """连接当前系统的端口到其他模型"""
        self.system.connect(self.system, port_self, other.system, port_other)

class ThrottleValveComponent(SubModelsBase):
    def default_params(self):
        return {
            "Cd": 0.6,
            "rho": 1000,
            "Q_nom": 2.0e-5,
            "Q_leak": 5e-8,
            "p_nom": 5e5,
            "p_vapour": 0
        }

    def build(self):
        p = self.params
        A = p["Q_nom"] / (p["Cd"] * np.sqrt((2 / p["rho"]) * p["p_nom"]))
        G_leak = p["Q_leak"] / p["p_nom"]
        R_throttle = bgt.new("R", value=1/(p["Cd"] * A), name="R_throttle")
        R_leak = bgt.new("R", value=1/G_leak, name="R_leak")
        zero = bgt.new("0", name="Node")
        bgt.add(self.system, R_throttle, R_leak, zero)
        bgt.connect(R_throttle, zero)
        bgt.connect(R_leak, zero)

if __name__ == "__main__":
    # 示例：创建一个节流阀组件并获取其bond graph系统
    throttle_valve1 = ThrottleValveComponent(name="ThrottleValve1").get_system()
    throttle_valve2 = ThrottleValveComponent(name="ThrottleValve1").get_system()
    system = bgt.new(name="MainSystem")
    zero = bgt.new("0", name="Node")
   
    from BondGraphTools import new, draw, simulate
    bgt.add(system, throttle_valve1, throttle_valve2, zero)
    bgt.connect(throttle_valve1, zero)
    bgt.connect(throttle_valve2, zero)
    bgt.draw(system)
    print(system)
    # 可以在这里进一步操作或分析system