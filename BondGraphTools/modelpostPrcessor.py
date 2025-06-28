import numpy as np
import sympy as sp
from collections import defaultdict

import numpy as np
import sympy as sp
from collections import deque, defaultdict

class BondGraphPostProcessor:
    def __init__(self, model, states_values=None, t=None):
        """
        键合图模型后处理器
        
        参数:
            model: BondGraphTools模型实例
        """
        self.model = model
        self.components = model.components
        self.bonds = model.bonds
        self.states_vars = model.state_vars
        self.processed = set()  # 已处理的元件
        self.results = {'t': t.T[0]}  # 存储所有变量的结果
        self.dt= (t[1] - t[0]).item()
        self.createBondGraph()
        self.process(states_values,t)
        
    def process(self, states_values, t=None):
        
        ### 处理的顺序 应该是储能元件,与储能相邻的元件，0-1结，tf,gy结，源点
        '''根据states_vars 找到对应是状态变量的元件，并计算其结果'''
        #  求解储能元件
        for state_str,comp_var in self.states_vars.items():
            comp,varname = comp_var
            idx = int(state_str.split('_')[1])
            self.results[f"{comp.name}_{varname}"] = states_values[:, idx]
            self.solve_component(comp, t)

        # 求解源元件
        # self.solve_source_components(t)    ## 这里可能需要调整到后面，最后来求解他，因为源的一个变量已知，另外一个变量通过结来求

        # 根据bonds的连接关系，求解其他元件变量
        for comp in self.processed:
            self.process_connected_components(comp, t)

    def solve_component(self, component, t, known_vars=None):
        """
        求解单个元件的变量
        
        参数:
            component: 元件实例
            states: 状态变量数组 (n_time, n_states)
            t: 时间数组
            known_vars: 已知变量字典 {var_name: values}
            
            results: 元件所有变量的时间序列
        """
        if component in self.processed:
            return
        # 获取元件所有变量名
        all_symbols = set()
        for eq in component.constitutive_relations:
            all_symbols |= eq.free_symbols
        
        # 检查是否已有足够信息求解
        known_vars=set() if known_vars is None else known_vars
        component.results={}  
        if component.state_vars:
            # 如果是状态变量元件，先求出其微分值表现如 dq_0
            for state_str, comp_var in component.state_vars.items():
                state_value = self.results.get(f"{component.name}_{state_str}", None)
                derivative_value=np.gradient(state_value, self.dt) if state_value is not None else None
                for symb in all_symbols:
                    if symb.name == state_str:
                        component.results[symb] = state_value
                        known_vars |= {symb}
                    elif symb.name == f"d{state_str}":
                        component.results[symb] = derivative_value
                        known_vars |= {symb}

                
        # 解析本构关系求解其他变量
        unknown_vars = all_symbols - known_vars
        # t = sp.symbols('t')
        # unknown_vars=[ sp.Function(var.name)(t) for var in unknown_vars]
        solutions=sp.solve(component.constitutive_relations, list(unknown_vars))
        f=sp.lambdify(tuple(known_vars), list(solutions.values()), modules='numpy')
        arrayvalue=tuple(component.results[var] for var in known_vars)
        values=f(*arrayvalue)
        component.results.update({var: value for var, value in zip(unknown_vars, values)})
        self.processed.add(component)
    
    def createBondGraph(self):
        """
        创建键合图的连接关系
        """
        self.bond_graph = defaultdict(list)
        for bond in self.bonds:
            head = bond.head.component
            tail = bond.tail.component
            self.bond_graph[head].append(tail)
            self.bond_graph[tail].append(head)

    def solve_source_components(self, t):
        """
        求解源元件的变量
        """
        for comp in self.components:
            if comp.template.split("/")[-1] in ['Se', 'Sei', 'Sf', 'Sf0', 'Se0', 'Sei0']:
                # 获取源元件的所有变量
                self.solve_component(comp, t)

    def process_connected_components(self, component, t):
        """
        处理与元件连接的其他元件
        """
        connected_components = self.bond_graph[component]
        for conn_comp in connected_components:
            if conn_comp in self.processed:
                continue
            # 获取连接元件的所有变量
            if conn_comp.template.split("/")[-1] in ['1']:  #  # 如果是1结，则所有的f都相等
                # 如果是源元件，直接求解
                f_0=sp.symbols('f_0')
                self.solve_component(conn_comp, t, known_vars=[f_0])
            elif conn_comp.template.split("/")[-1] in ['0']:
                self.solve_component(conn_comp, t, known_vars=['e_0'])
            elif conn_comp.template.split("/")[-1] in ['tf', 'gy']:
                # 如果是tf或gy结，求解其已知变量
                self.solve_component(conn_comp, t, known_vars=['f_0', 'e_0'])
                # 否则求解连接元件
            else:
                # 对于其他元件，求解其所有变量
                self.solve_component(conn_comp, t)

if __name__ == "__main__":
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
    
    post_process_results = BondGraphPostProcessor(model,x,t)
    
