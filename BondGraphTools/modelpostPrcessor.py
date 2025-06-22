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
        self.process(states_values,t)
        
    def process(self, states_values, t=None):
        
        '''根据states_vars 找到对应是状态变量的元件，并计算其结果'''
        for state_str,comp_var in self.states_vars.items():
            comp,varname = comp_var
            idx = int(state_str.split('_')[1])
            self.results[f"{comp.name}_{varname}"] = states_values[:, idx]
            self.solve_component(comp, t)

        # 处理非状态变量元件
        for comp in self.components:
            if comp not in self.processed:
                self.solve_component(comp, t)

    def solve_component(self, component, t, known_vars=None):
        """
        求解单个元件的变量
        
        参数:
            component: 元件实例
            states: 状态变量数组 (n_time, n_states)
            t: 时间数组
            known_vars: 已知变量字典 {var_name: values}
            
        返回:
            dict: 元件所有变量的时间序列
        """
        if component in self.processed:
            return
        # 获取元件所有变量名
        all_vars = set()
        for rel in component.constitutive_relations:
            try:
                expr = sp.sympify(rel)
                all_vars.update({str(s) for s in expr.free_symbols})
            except:
                pass
        
        # 检查是否已有足够信息求解
        known_vars = [var for var in all_vars if f"{component.name}_{var}" in self.results]
        if len(known_vars) < 1:  # 需要至少n-1个已知变量
            return
        
        
        # 解析本构关系求解其他变量
        for rel in component.constitutive_relations:
            try:
                expr = sp.sympify(rel)
                unknowns = [str(s) for s in expr.free_symbols 
                          if f"{component.name}_{str(s)}" not in self.results]
                
                # 只有一个未知量时求解
                if len(unknowns) == 1:
                    # 直接解析方程而不是逐点计算
                    param_dict = {k: v for k, v in component.params.items()}
                    known_dict = {
                        k: self.results[f"{component.name}_{k}"] 
                        for k in all_vars if k != unknowns[0]
                    }
                    
                    # 向量化求解
                    sol = self._solve_equation_vectorized(
                        expr, 
                        unknowns[0], 
                        known_dict, 
                        param_dict
                    )
                    
                    if sol is not None:
                        self.results[f"{component.name}_{unknowns[0]}"] = sol
            except Exception as e:
                print(f"Error solving {component.name}: {str(e)}")
        
        self.processed.add(component)
        
        # 计算功率和能量
        # self._calculate_power(component, results)
        
    
    def _solve_equation_vectorized(self, equation, unknown, known_dict, param_dict):
        """向量化求解方程"""
        try:
            # 求解符号表达式
            sol_expr = sp.solve(equation, unknown, dict=True)[0][sp.symbols(unknown)]
            
            # 创建替换函数
            symbols = {sp.symbols(k): v for k, v in {**known_dict, **param_dict}.items()}
            
            # 向量化计算
            return sp.lambdify(list(symbols.keys()), sol_expr, 'numpy')(**{
                str(k): v for k, v in symbols.items()
            })
        except:
            return None
    
    def _calculate_power(self, component, results):
        """
        计算元件的功率
        """
        # 识别端口功率
        power_vars = []
        for i in range(len(component.ports)):
            e_var = f"e_{i}"
            f_var = f"f_{i}"
            if e_var in results and f_var in results:
                power_var = f"P_{i}"
                results[power_var] = results[e_var] * results[f_var]
                power_vars.append(power_var)
        
        # 计算总功率
        if power_vars:
            results['P_total'] = np.zeros_like(results['t'])
            for p_var in power_vars:
                results['P_total'] += results[p_var]
        
        # 计算能量（积分功率）
        if 'P_total' in results and len(results['t']) > 1:
            results['E'] = np.cumsum(results['P_total']) * (results['t'][1] - results['t'][0])
    


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
    
