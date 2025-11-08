import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union
import matplotlib.pyplot as plt
from sympy import Symbol, lambdify

class BondGraphPost:
    """
    Bond Graph模型的通用后处理类
    
    基于BondGraphTools的符号系统，直接从模型方程提取和计算：
    - 状态变量(广义动量、广义位移)
    - 流(flow)和势(effort)
    - 能量和功率
    - 所有键合图变量
    """
    
    def __init__(self, model, simulation_results: Union[Tuple, np.ndarray] = None):
        """
        初始化后处理类
        
        Parameters:
        -----------
        model : BondGraph模型对象
            已构建的Bond Graph模型
        simulation_results : tuple or array-like, optional
            仿真结果 (t, x) 或仅 x
        """
        self.model = model
        
        # 处理仿真结果
        if simulation_results is not None:
            if isinstance(simulation_results, tuple) and len(simulation_results) == 2:
                self.t, self.x = simulation_results
                self.t = np.array(self.t).flatten() if self.t is not None else None
            else:
                self.t = None
                self.x = simulation_results
            self.x = np.array(self.x)
        else:
            self.t = None
            self.x = None
        
        # 获取模型的符号系统和映射
        self.X, self.mapping, self.A, self.F, self.G = model.system_model()
        
        # 获取完整方程
        self.full_equations = self._build_full_equations()
        self.variable_names = list(self.full_equations.keys())
        
        # 创建计算函数
        self.eval_func = self._create_eval_function()
        
        # 计算结果
        self.results_df = self._compute_results()
        
    @staticmethod
    def find_port(component, direction):
        """查找组件的端口"""
        if direction in ['f', 'forward']:
            index = 0
        elif direction in ['r', 'reverse']:
            index = 1
        return list(component.ports.keys())[index]
    
    def _rename_variables(self) -> Tuple[List[str], Dict[str, str]]:
        """
        将方程中的变量替换为有物理意义的名称
        
        Returns:
        --------
        renamed_vars : list
            重命名后的变量列表
        replacements : dict
            变量替换映射
        """
        replacements = {}
        renamed_vars = []
        
        # 处理状态变量: x_0 -> C1_q, dx_0 -> C1_dq
        for (comp, state), idx in self.mapping[0].items():
            # 获取组件名（去掉命名空间前缀）
            comp_name = comp.name.split(':')[-1] if hasattr(comp, 'name') else f'comp{idx}'
            
            # 状态变量名
            state_var = f'{comp_name}_{state}'
            replacements[f'x_{idx}'] = state_var
            
            # 状态导数名
            state_deriv = f'{comp_name}_d{state}'
            replacements[f'dx_{idx}'] = state_deriv
        
        # 处理端口变量 (bond variables): e_0, f_0 -> comp_p0_e, comp_p0_f
        for port, bond_idx in self.mapping[1].items():
            # 获取组件信息
            comp = port.component
            comp_name = comp.name.split(':')[-1] if hasattr(comp, 'name') else f'comp{bond_idx}'
            
            # 对于junction(0或1节点)，添加前缀'j'
            if hasattr(comp, 'metamodel') and comp.metamodel in ['0', '1']:
                comp_name = f'j{comp.metamodel}_{comp_name}'
            
            port_name = str(port.index) if hasattr(port, 'index') else '0'
            
            # effort 和 flow
            effort_var = f'{comp_name}_p{port_name}_e'
            flow_var = f'{comp_name}_p{port_name}_f'
            
            replacements[f'e_{bond_idx}'] = effort_var
            replacements[f'f_{bond_idx}'] = flow_var
        
        return replacements
    
    def _build_full_equations(self) -> Dict[str, any]:
        """
        构建完整的模型方程
        
        Returns:
        --------
        equations : dict
            变量名 -> 符号表达式的映射
        """
        from sympy import SparseMatrix
        
        # 获取变量重命名映射
        replacements = self._rename_variables()
        
        # 构建方程: AX + F(X) = 0 => X_i = (AX + F)_i
        # 创建符号变量列表
        symbol_vars = []
        for old_name, new_name in sorted(replacements.items(), key=lambda x: x[0]):
            if old_name.startswith('x_') or old_name.startswith('dx_'):
                symbol_vars.append(Symbol(new_name))
        
        # 如果没有符号变量，使用原始状态变量
        if not symbol_vars:
            symbol_vars = [Symbol(x) for x in self.model.state_vars.keys()]
            var_matrix = SparseMatrix(self.model.state_vars.keys())
        else:
            var_matrix = SparseMatrix([str(v) for v in symbol_vars])
        
        # 计算 AX + F
        try:
            AX_F = self.A * var_matrix + self.F
        except:
            # 如果矩阵运算失败，使用原始X
            AX_F = self.A * SparseMatrix(self.X) + self.F
        
        # 构建方程字典
        equations = {}
        for i in range(len(self.X)):
            xi = str(self.X[i])
            
            # 重命名变量
            if xi in replacements:
                var_name = replacements[xi]
            else:
                var_name = xi
            
            # 方程: var = expression
            eqn = AX_F[i, 0]
            equations[var_name] = eqn
        
        return equations
    
    def _create_eval_function(self):
        """创建可以数值计算的函数"""
        from sympy import symbols
        
        # 获取所有符号变量（包括状态变量和导数）
        all_symbols = set()
        for expr in self.full_equations.values():
            all_symbols.update(expr.free_symbols)
        
        # 按名称排序，确保顺序一致
        sorted_symbols = sorted(all_symbols, key=lambda s: str(s))
        
        # 获取表达式
        expressions = tuple(self.full_equations.values())
        
        # 如果没有找到符号，使用状态变量
        if not sorted_symbols:
            sorted_symbols = [Symbol(x) for x in self.model.state_vars.keys()]
        
        # 创建 lambda 函数
        return lambdify(sorted_symbols, expressions, modules='numpy')
    
    def _compute_results(self) -> pd.DataFrame:
        """计算所有变量的数值结果"""
        if self.x is None:
            return pd.DataFrame()
        
        # 确保 x 是二维数组
        if self.x.ndim == 1:
            x_2d = self.x.reshape(1, -1)
        else:
            x_2d = self.x
        
        # 计算导数 dx/dt
        if self.t is not None and len(self.t) > 1:
            dx_dt = np.gradient(x_2d, self.t, axis=0)
        else:
            dx_dt = np.zeros_like(x_2d)
        
        # 获取需要的符号变量
        all_symbols = set()
        for expr in self.full_equations.values():
            all_symbols.update(expr.free_symbols)
        sorted_symbols = sorted(all_symbols, key=lambda s: str(s))
        
        # 对每一行状态值计算表达式
        all_results = []
        for i, row in enumerate(x_2d):
            try:
                # 准备输入参数：状态变量 + 导数
                input_values = []
                for sym in sorted_symbols:
                    sym_str = str(sym)
                    if sym_str.startswith('dx_'):
                        # 导数变量
                        idx = int(sym_str.split('_')[1])
                        input_values.append(dx_dt[i, idx])
                    elif sym_str.startswith('x_'):
                        # 状态变量
                        idx = int(sym_str.split('_')[1])
                        input_values.append(row[idx])
                    else:
                        # 其他变量（如重命名后的变量）
                        # 尝试从状态变量中匹配
                        found = False
                        for j, state_name in enumerate(self.model.state_vars.keys()):
                            if state_name in sym_str:
                                input_values.append(row[j])
                                found = True
                                break
                        if not found:
                            input_values.append(0.0)  # 默认值
                
                result = self.eval_func(*input_values)
                all_results.append(result)
            except Exception as e:
                print(f"Warning: 计算失败 {e}")
                all_results.append([np.nan] * len(self.variable_names))
        
        # 转换为数组 (n_samples, n_vars)
        result_array = np.array(all_results)
        
        # 创建 DataFrame
        df_data = {}
        
        # 添加时间
        if self.t is not None:
            df_data['time'] = self.t
        
        # 添加状态变量
        for i, state_name in enumerate(self.model.state_vars.keys()):
            df_data[f'state_{state_name}'] = x_2d[:, i]
        
        # 添加状态导数
        for i, state_name in enumerate(self.model.state_vars.keys()):
            df_data[f'deriv_d{state_name}'] = dx_dt[:, i]
        
        # 添加计算的变量
        for i, var_name in enumerate(self.variable_names):
            df_data[var_name] = result_array[:, i]
        
        return pd.DataFrame(df_data)
    
    def get_state_variables(self) -> Dict[str, np.ndarray]:
        """获取所有状态变量"""
        states = {}
        for col in self.results_df.columns:
            if col.startswith('state_'):
                states[col.replace('state_', '')] = self.results_df[col].values
        return states
    
    def get_flows(self) -> Dict[str, np.ndarray]:
        """获取所有流变量"""
        flows = {}
        for col in self.results_df.columns:
            if '_f' in col or 'flow' in col.lower():
                flows[col] = self.results_df[col].values
        return flows
    
    def get_efforts(self) -> Dict[str, np.ndarray]:
        """获取所有势变量"""
        efforts = {}
        for col in self.results_df.columns:
            if '_e' in col or 'effort' in col.lower():
                efforts[col] = self.results_df[col].values
        return efforts
    
    def get_generalized_displacement(self) -> Dict[str, np.ndarray]:
        """获取广义位移 (对应C元件的q)"""
        displacements = {}
        for col in self.results_df.columns:
            if '_q' in col or col.startswith('state_q'):
                displacements[col] = self.results_df[col].values
        return displacements
    
    def get_generalized_momentum(self) -> Dict[str, np.ndarray]:
        """获取广义动量 (对应I元件的p)"""
        momenta = {}
        for col in self.results_df.columns:
            if '_p' in col or col.startswith('state_p'):
                momenta[col] = self.results_df[col].values
        return momenta
    
    def get_component_variables(self, component_name: str) -> pd.DataFrame:
        """获取特定组件的所有变量"""
        matching_cols = [col for col in self.results_df.columns 
                        if component_name in col]
        return self.results_df[['time'] + matching_cols] if 'time' in self.results_df.columns else self.results_df[matching_cols]
    
    def calculate_energy(self) -> Dict[str, np.ndarray]:
        """计算储能元件的能量"""
        energies = {}
        
        # 对于C元件: E = 0.5 * q^2 / C
        for col in self.results_df.columns:
            if '_q' in col and (col.startswith('state_') or 'C' in col):
                comp_name = col.split('_')[0].replace('state_', '')
                q = self.results_df[col].values
                
                # 尝试获取C值
                try:
                    C_value = self._get_component_parameter(comp_name, 'C')
                    if C_value is not None and C_value != 0:
                        energies[f'energy_{comp_name}'] = 0.5 * q**2 / C_value
                except:
                    pass
        
        # 对于I元件: E = 0.5 * p^2 / I
        for col in self.results_df.columns:
            if '_p' in col and (col.startswith('state_') or 'I' in col):
                comp_name = col.split('_')[0].replace('state_', '')
                p = self.results_df[col].values
                
                # 尝试获取I值
                try:
                    I_value = self._get_component_parameter(comp_name, 'I')
                    if I_value is not None and I_value != 0:
                        energies[f'energy_{comp_name}'] = 0.5 * p**2 / I_value
                except:
                    pass
        
        return energies
    
    def calculate_power(self) -> Dict[str, np.ndarray]:
        """计算功率 P = e * f"""
        powers = {}
        
        # 找到所有的 effort-flow 配对
        for col in self.results_df.columns:
            if '_e' in col:
                # 尝试找到对应的flow
                flow_col = col.replace('_e', '_f')
                if flow_col in self.results_df.columns:
                    effort = self.results_df[col].values
                    flow = self.results_df[flow_col].values
                    power_name = col.replace('_e', '_power')
                    powers[power_name] = effort * flow
        
        return powers
    
    def _get_component_parameter(self, comp_name: str, param_type: str) -> Optional[float]:
        """获取组件参数值"""
        try:
            # 遍历模型组件
            for attr_name in dir(self.model):
                if not attr_name.startswith('_'):
                    comp = getattr(self.model, attr_name)
                    if hasattr(comp, 'name') and comp_name in comp.name:
                        if hasattr(comp, 'value'):
                            return float(comp.value)
        except:
            pass
        return None
    
    def eval_at_state(self, *state_values, state_derivs=None):
        """
        在特定状态值处计算所有变量
        
        Parameters:
        -----------
        *state_values : float
            状态变量的值，顺序与 model.state_vars 一致
        state_derivs : array-like, optional
            状态导数的值，如果为None则设为0
        
        Returns:
        --------
        result : dict
            变量名 -> 值的映射
        """
        # 获取需要的符号变量
        all_symbols = set()
        for expr in self.full_equations.values():
            all_symbols.update(expr.free_symbols)
        sorted_symbols = sorted(all_symbols, key=lambda s: str(s))
        
        # 准备输入值
        input_values = []
        n_states = len(state_values)
        
        if state_derivs is None:
            state_derivs = [0.0] * n_states
        
        for sym in sorted_symbols:
            sym_str = str(sym)
            if sym_str.startswith('dx_'):
                idx = int(sym_str.split('_')[1])
                input_values.append(state_derivs[idx] if idx < len(state_derivs) else 0.0)
            elif sym_str.startswith('x_'):
                idx = int(sym_str.split('_')[1])
                input_values.append(state_values[idx] if idx < len(state_values) else 0.0)
            else:
                # 尝试匹配状态变量名
                found = False
                for j, state_name in enumerate(self.model.state_vars.keys()):
                    if state_name in sym_str:
                        input_values.append(state_values[j] if j < len(state_values) else 0.0)
                        found = True
                        break
                if not found:
                    input_values.append(0.0)
        
        result_values = self.eval_func(*input_values)
        return dict(zip(self.variable_names, result_values))
    
    def to_dataframe(self) -> pd.DataFrame:
        """返回完整的结果 DataFrame"""
        return self.results_df
    
    def plot(self, variables: List[str] = None, figsize=(12, 8), 
             plot_type: str = 'all'):
        """
        绘制结果
        
        Parameters:
        -----------
        variables : list of str, optional
            要绘制的变量列表，如果为None则根据plot_type选择
        figsize : tuple
            图形大小
        plot_type : str
            绘图类型: 'all', 'states', 'flows', 'efforts', 'energy', 'power'
        """
        if 'time' not in self.results_df.columns:
            print("无时间数据，无法绘图")
            return
        
        # 根据 plot_type 选择变量
        if variables is None:
            if plot_type == 'states':
                variables = [col for col in self.results_df.columns if col.startswith('state_')]
            elif plot_type == 'flows':
                variables = [col for col in self.results_df.columns if '_f' in col]
            elif plot_type == 'efforts':
                variables = [col for col in self.results_df.columns if '_e' in col]
            elif plot_type == 'energy':
                energies = self.calculate_energy()
                if energies:
                    # 临时添加能量到DataFrame用于绘图
                    temp_df = self.results_df.copy()
                    for name, values in energies.items():
                        temp_df[name] = values
                    variables = list(energies.keys())
                    self.results_df = temp_df
                else:
                    print("没有能量数据")
                    return
            elif plot_type == 'power':
                powers = self.calculate_power()
                if powers:
                    temp_df = self.results_df.copy()
                    for name, values in powers.items():
                        temp_df[name] = values
                    variables = list(powers.keys())
                    self.results_df = temp_df
                else:
                    print("没有功率数据")
                    return
            else:  # 'all'
                variables = [col for col in self.results_df.columns if col != 'time']
        
        # 过滤存在的变量
        variables = [v for v in variables if v in self.results_df.columns]
        
        if len(variables) == 0:
            print("没有可绘制的变量")
            return
        
        # 创建子图
        n_vars = len(variables)
        n_cols = min(2, n_vars)
        n_rows = (n_vars + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
        axes = axes.flatten()
        
        # 绘制每个变量
        time = self.results_df['time'].values
        for i, var in enumerate(variables):
            axes[i].plot(time, self.results_df[var].values, linewidth=2)
            axes[i].set_xlabel('Time')
            axes[i].set_ylabel(var)
            axes[i].set_title(var)
            axes[i].grid(True, alpha=0.3)
        
        # 隐藏多余的子图
        for i in range(n_vars, len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        plt.show()
    
    def summary(self):
        """打印结果摘要"""
        print("=" * 70)
        print("Bond Graph 后处理结果摘要")
        print("=" * 70)
        
        print(f"\n模型名称: {self.model.name if hasattr(self.model, 'name') else 'Unknown'}")
        
        states = self.get_state_variables()
        print(f"\n状态变量数量: {len(states)}")
        for name in states.keys():
            print(f"  - {name}")
        
        flows = self.get_flows()
        print(f"\n流变量数量: {len(flows)}")
        for name in list(flows.keys())[:5]:  # 只显示前5个
            print(f"  - {name}")
        if len(flows) > 5:
            print(f"  ... 还有 {len(flows)-5} 个")
        
        efforts = self.get_efforts()
        print(f"\n势变量数量: {len(efforts)}")
        for name in list(efforts.keys())[:5]:
            print(f"  - {name}")
        if len(efforts) > 5:
            print(f"  ... 还有 {len(efforts)-5} 个")
        
        print(f"\n总变量数量: {len(self.variable_names)}")
        
        if 'time' in self.results_df.columns:
            t = self.results_df['time'].values
            print(f"\n仿真时间范围: [{t[0]:.4f}, {t[-1]:.4f}]")
            print(f"时间步数: {len(t)}")
        
        print("=" * 70)
    
    def export_csv(self, filename: str):
        """导出结果到CSV文件"""
        self.results_df.to_csv(filename, index=False)
        print(f"结果已导出到: {filename}")


# 使用示例
if __name__ == "__main__":

    # 1. 创建和仿真Bond Graph模型
    from BondGraphTools import new, add, connect, simulate
    
    model = new(name='RC')
    C = new("C", value=1.0)
    R = new("R", value=1.0)
    se = new("Se", value=1.0)
    one = new("1")
    
    add(model, R, C, one, se)
    connect(se, one)
    connect(R, one)
    connect(C, one)
    
    timespan = [0, 5]
    x0 = {'x_0': 1}
    t, x = simulate(model, timespan=timespan, x0=x0)
    
    # 2. 创建后处理对象
    post = BondGraphPost(model, (t, x))
    
    # 3. 查看摘要
    post.summary()
    
    # 4. 获取各类变量
    states = post.get_state_variables()
    flows = post.get_flows()
    efforts = post.get_efforts()
    energies = post.calculate_energy()
    powers = post.calculate_power()
    
    # 5. 获取特定组件的变量
    c_vars = post.get_component_variables('C')
    print(c_vars)
    
    # 6. 在特定状态计算
    result = post.eval_at_state(1.0)  # 在状态x=1.0处计算
    print(result)
    
    # 7. 获取完整DataFrame
    df = post.to_dataframe()
    print(df.head())
    
    # 8. 绘图
    post.plot(plot_type='states')     # 只绘制状态变量
    post.plot(plot_type='flows')      # 只绘制流
    post.plot(plot_type='efforts')    # 只绘制势
    post.plot(plot_type='energy')     # 绘制能量
    post.plot(plot_type='power')      # 绘制功率
    
    # 9. 导出CSV
    post.export_csv('results.csv')
