import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union, Any
import matplotlib.pyplot as plt
from sympy import Symbol, lambdify
from collections import defaultdict

class ComponentData:
    """单个组件的数据容器"""
    
    def __init__(self, name: str, component_type: str, component=None):
        self.name = name
        self.type = component_type
        self.component = component
        
        # 组件参数
        self.parameters = {}
        
        # 状态变量
        self.states = {}  # {state_name: array}
        self.state_derivatives = {}  # {state_name: array}
        
        # 端口变量
        self.ports = defaultdict(dict)  # {port_id: {effort: array, flow: array}}
        
        # 能量和功率
        self.energy = None
        self.power = {}  # {port_id: array}
        
        # 子组件（如果是子模型）
        self.subcomponents = {}  # {name: ComponentData}
        
    def add_state(self, state_name: str, values: np.ndarray, derivative: np.ndarray = None):
        """添加状态变量"""
        self.states[state_name] = values
        if derivative is not None:
            self.state_derivatives[state_name] = derivative
    
    def add_port_variable(self, port_id: str, var_type: str, values: np.ndarray):
        """添加端口变量 (effort 或 flow)"""
        self.ports[port_id][var_type] = values
    
    def add_parameter(self, param_name: str, value: float):
        """添加组件参数"""
        self.parameters[param_name] = value
    
    def add_subcomponent(self, subcomp: 'ComponentData'):
        """添加子组件"""
        self.subcomponents[subcomp.name] = subcomp
    
    def calculate_energy(self):
        """计算储能元件的能量"""
        if self.type == 'C' and self.parameters.get('value'):
            # E = 0.5 * q^2 / C
            for state_name, q in self.states.items():
                if 'q' in state_name or state_name == 'x':
                    C = self.parameters['value']
                    self.energy = 0.5 * q**2 / C if C != 0 else None
                    break
        
        elif self.type == 'I' and self.parameters.get('value'):
            # E = 0.5 * p^2 / I
            for state_name, p in self.states.items():
                if 'p' in state_name:
                    I = self.parameters['value']
                    self.energy = 0.5 * p**2 / I if I != 0 else None
                    break
    
    def calculate_power(self):
        """计算各端口的功率 P = e * f"""
        for port_id, port_vars in self.ports.items():
            if 'effort' in port_vars and 'flow' in port_vars:
                self.power[port_id] = port_vars['effort'] * port_vars['flow']
    
    def to_dict(self, include_subcomponents: bool = True) -> Dict:
        """转换为字典格式，方便序列化"""
        data = {
            'name': self.name,
            'type': self.type,
            'parameters': self.parameters,
            'states': {k: v.tolist() if isinstance(v, np.ndarray) else v 
                      for k, v in self.states.items()},
            'state_derivatives': {k: v.tolist() if isinstance(v, np.ndarray) else v 
                                 for k, v in self.state_derivatives.items()},
            'ports': {port_id: {k: v.tolist() if isinstance(v, np.ndarray) else v 
                               for k, v in port_vars.items()}
                     for port_id, port_vars in self.ports.items()},
            'energy': self.energy.tolist() if isinstance(self.energy, np.ndarray) else self.energy,
            'power': {k: v.tolist() if isinstance(v, np.ndarray) else v 
                     for k, v in self.power.items()},
        }
        
        if include_subcomponents and self.subcomponents:
            data['subcomponents'] = {name: comp.to_dict(True) 
                                    for name, comp in self.subcomponents.items()}
        
        return data
    
    def get_dataframe(self, time: np.ndarray = None) -> pd.DataFrame:
        """获取该组件的DataFrame"""
        data = {}
        
        if time is not None:
            data['time'] = time
        
        # 状态变量
        for state_name, values in self.states.items():
            data[f'{self.name}_state_{state_name}'] = values
        
        # 状态导数
        for state_name, values in self.state_derivatives.items():
            data[f'{self.name}_d{state_name}_dt'] = values
        
        # 端口变量
        for port_id, port_vars in self.ports.items():
            for var_type, values in port_vars.items():
                data[f'{self.name}_port{port_id}_{var_type}'] = values
        
        # 能量
        if self.energy is not None:
            data[f'{self.name}_energy'] = self.energy
        
        # 功率
        for port_id, values in self.power.items():
            data[f'{self.name}_port{port_id}_power'] = values
        
        return pd.DataFrame(data)
    
    def __repr__(self):
        return f"ComponentData(name='{self.name}', type='{self.type}', states={len(self.states)}, ports={len(self.ports)})"


class BondGraphPost:
    """
    Bond Graph模型的层次化后处理类
    
    支持：
    - 组件级数据组织
    - 子模型嵌套
    - 层次化数据访问
    - UI友好的数据结构
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
        
        # 组件数据字典 {component_name: ComponentData}
        self.components = {}
        
        # 获取完整方程
        self.full_equations = self._build_full_equations()
        self.variable_names = list(self.full_equations.keys())
        
        # 创建计算函数
        self.eval_func = self._create_eval_function()
        
        # 处理结果并组织到组件中
        self._process_components()
        
    def _get_component_info(self, comp) -> Tuple[str, str]:
        """获取组件名称和类型"""
        name = comp.name.split(':')[-1] if hasattr(comp, 'name') else 'unknown'
        comp_type = comp.metamodel if hasattr(comp, 'metamodel') else \
                   (comp.component_type if hasattr(comp, 'component_type') else 'unknown')
        return name, comp_type
    
    def _extract_model_components(self, model=None, parent_path: str = "") -> Dict[str, Any]:
        """
        递归提取模型中的所有组件（包括子模型）
        
        Returns:
        --------
        components : dict
            {full_path: component_object}
        """
        if model is None:
            model = self.model
        
        components = {}
        
        # 方法1: 使用 components 属性
        if hasattr(model, 'components'):
            for comp in model.components:
                name, comp_type = self._get_component_info(comp)
                full_path = f"{parent_path}.{name}" if parent_path else name
                components[full_path] = comp
                
                # 如果是子模型，递归提取
                if hasattr(comp, 'components'):
                    sub_components = self._extract_model_components(comp, full_path)
                    components.update(sub_components)
        
        # 方法2: 遍历属性
        else:
            for attr_name in dir(model):
                if not attr_name.startswith('_'):
                    try:
                        comp = getattr(model, attr_name)
                        if hasattr(comp, 'metamodel') or hasattr(comp, 'component_type'):
                            name, comp_type = self._get_component_info(comp)
                            full_path = f"{parent_path}.{name}" if parent_path else name
                            components[full_path] = comp
                            
                            # 递归子模型
                            if hasattr(comp, 'components'):
                                sub_components = self._extract_model_components(comp, full_path)
                                components.update(sub_components)
                    except:
                        pass
        
        return components
    
    def _rename_variables(self) -> Dict[str, str]:
        """将方程中的变量替换为有物理意义的名称"""
        replacements = {}
        
        # 处理状态变量
        for (comp, state), idx in self.mapping[0].items():
            comp_name, comp_type = self._get_component_info(comp)
            replacements[f'x_{idx}'] = f'{comp_name}_{state}'
            replacements[f'dx_{idx}'] = f'd{comp_name}_{state}'
        
        # 处理端口变量
        for port, bond_idx in self.mapping[1].items():
            comp = port.component
            comp_name, comp_type = self._get_component_info(comp)
            
            if comp_type in ['0', '1']:
                comp_name = f'j{comp_type}_{comp_name}'
            
            port_id = str(port.index) if hasattr(port, 'index') else '0'
            
            replacements[f'e_{bond_idx}'] = f'{comp_name}_p{port_id}_e'
            replacements[f'f_{bond_idx}'] = f'{comp_name}_p{port_id}_f'
        
        return replacements
    
    def _build_full_equations(self) -> Dict[str, any]:
        """构建完整的模型方程"""
        from sympy import SparseMatrix
        
        replacements = self._rename_variables()
        
        # 创建符号变量
        symbol_vars = []
        for old_name in sorted([k for k in replacements.keys() if k.startswith('x_')]):
            symbol_vars.append(Symbol(replacements[old_name]))
        
        if not symbol_vars:
            symbol_vars = [Symbol(x) for x in self.model.state_vars.keys()]
        
        var_matrix = SparseMatrix([str(v) for v in symbol_vars])
        
        try:
            AX_F = self.A * var_matrix + self.F
        except:
            AX_F = self.A * SparseMatrix(self.X) + self.F
        
        equations = {}
        for i in range(len(self.X)):
            xi = str(self.X[i])
            var_name = replacements.get(xi, xi)
            equations[var_name] = AX_F[i, 0]
        
        return equations
    
    def _create_eval_function(self):
        """创建可以数值计算的函数"""
        all_symbols = set()
        for expr in self.full_equations.values():
            all_symbols.update(expr.free_symbols)
        
        sorted_symbols = sorted(all_symbols, key=lambda s: str(s))
        
        if not sorted_symbols:
            sorted_symbols = [Symbol(x) for x in self.model.state_vars.keys()]
        
        expressions = tuple(self.full_equations.values())
        return lambdify(sorted_symbols, expressions, modules='numpy')
    
    def _compute_all_variables(self) -> Dict[str, np.ndarray]:
        """计算所有变量的数值"""
        if self.x is None:
            return {}
        
        # 确保 x 是二维数组
        if self.x.ndim == 1:
            x_2d = self.x.reshape(1, -1)
        else:
            x_2d = self.x
        
        # 计算导数
        if self.t is not None and len(self.t) > 1:
            dx_dt = np.gradient(x_2d, self.t, axis=0)
        else:
            dx_dt = np.zeros_like(x_2d)
        
        # 获取符号变量
        all_symbols = set()
        for expr in self.full_equations.values():
            all_symbols.update(expr.free_symbols)
        sorted_symbols = sorted(all_symbols, key=lambda s: str(s))
        
        # 计算每个时间点的结果
        all_results = []
        for i, row in enumerate(x_2d):
            input_values = []
            for sym in sorted_symbols:
                sym_str = str(sym)
                
                # 匹配导数变量
                if sym_str.startswith('d') and '_' in sym_str:
                    # 尝试找到对应的状态索引
                    for j, state_name in enumerate(self.model.state_vars.keys()):
                        if state_name in sym_str:
                            input_values.append(dx_dt[i, j])
                            break
                    else:
                        input_values.append(0.0)
                
                # 匹配状态变量
                else:
                    found = False
                    for j, state_name in enumerate(self.model.state_vars.keys()):
                        if state_name in sym_str:
                            input_values.append(row[j])
                            found = True
                            break
                    if not found:
                        input_values.append(0.0)
            
            try:
                result = self.eval_func(*input_values)
                all_results.append(result)
            except Exception as e:
                print(f"Warning: 计算失败 at time index {i}: {e}")
                all_results.append([np.nan] * len(self.variable_names))
        
        # 转换为字典
        result_dict = {}
        result_array = np.array(all_results)
        
        for i, var_name in enumerate(self.variable_names):
            result_dict[var_name] = result_array[:, i]
        
        # 添加状态变量和导数
        for i, state_name in enumerate(self.model.state_vars.keys()):
            result_dict[f'state_{state_name}'] = x_2d[:, i]
            result_dict[f'deriv_{state_name}'] = dx_dt[:, i]
        
        return result_dict
    
    def _process_components(self):
        """处理结果并组织到各个组件中"""
        # 计算所有变量
        all_vars = self._compute_all_variables()
        
        if not all_vars:
            return
        
        # 提取所有组件（包括子模型）
        all_components = self._extract_model_components()
        
        # 为每个组件创建 ComponentData
        for comp_path, comp in all_components.items():
            comp_name, comp_type = self._get_component_info(comp)
            
            comp_data = ComponentData(comp_path, comp_type, comp)
            
            # 添加参数
            if hasattr(comp, 'value'):
                comp_data.add_parameter('value', float(comp.value))
            
            # 匹配该组件的变量
            for var_name, values in all_vars.items():
                # 检查变量名是否包含组件名
                if comp_name in var_name or comp_path in var_name:
                    
                    # 状态变量
                    if 'state_' in var_name:
                        state_name = var_name.split('state_')[-1]
                        comp_data.add_state(state_name, values)
                    
                    # 状态导数
                    elif 'deriv_' in var_name:
                        deriv_name = var_name.split('deriv_')[-1]
                        if deriv_name not in comp_data.state_derivatives:
                            comp_data.state_derivatives[deriv_name] = values
                    
                    # 端口变量
                    elif '_p' in var_name and ('_e' in var_name or '_f' in var_name):
                        # 解析端口ID和变量类型
                        parts = var_name.split('_p')
                        if len(parts) >= 2:
                            port_part = parts[1]
                            if '_e' in port_part:
                                port_id = port_part.split('_e')[0]
                                comp_data.add_port_variable(port_id, 'effort', values)
                            elif '_f' in port_part:
                                port_id = port_part.split('_f')[0]
                                comp_data.add_port_variable(port_id, 'flow', values)
            
            # 计算能量和功率
            comp_data.calculate_energy()
            comp_data.calculate_power()
            
            self.components[comp_path] = comp_data
        
        # 处理层次结构（将子组件关联到父组件）
        self._organize_hierarchy()
    
    def _organize_hierarchy(self):
        """组织组件的层次结构"""
        # 找出所有根组件（没有父路径的）
        root_components = {}
        nested_components = {}
        
        for path, comp_data in self.components.items():
            if '.' in path:
                nested_components[path] = comp_data
            else:
                root_components[path] = comp_data
        
        # 将嵌套组件添加到父组件
        for path, comp_data in nested_components.items():
            parts = path.split('.')
            if len(parts) >= 2:
                parent_path = '.'.join(parts[:-1])
                if parent_path in self.components:
                    self.components[parent_path].add_subcomponent(comp_data)
    
    def get_component(self, component_name: str) -> Optional[ComponentData]:
        """获取特定组件的数据"""
        # 精确匹配
        if component_name in self.components:
            return self.components[component_name]
        
        # 模糊匹配（匹配最后一段名称）
        for path, comp_data in self.components.items():
            if path.endswith(component_name) or component_name in path:
                return comp_data
        
        return None
    
    def get_all_components(self, include_junctions: bool = False) -> Dict[str, ComponentData]:
        """
        获取所有组件
        
        Parameters:
        -----------
        include_junctions : bool
            是否包含junction节点(0, 1)
        """
        if include_junctions:
            return self.components
        
        return {name: comp for name, comp in self.components.items() 
                if not comp.type in ['0', '1']}
    
    def get_components_by_type(self, comp_type: str) -> Dict[str, ComponentData]:
        """按类型获取组件"""
        return {name: comp for name, comp in self.components.items() 
                if comp.type == comp_type}
    
    def list_components(self, indent: int = 0):
        """列出所有组件的层次结构"""
        def print_component(comp: ComponentData, level: int):
            prefix = "  " * level
            print(f"{prefix}├─ {comp.name} ({comp.type})")
            
            if comp.states:
                print(f"{prefix}│  States: {list(comp.states.keys())}")
            if comp.ports:
                print(f"{prefix}│  Ports: {list(comp.ports.keys())}")
            if comp.subcomponents:
                print(f"{prefix}│  Subcomponents:")
                for sub in comp.subcomponents.values():
                    print_component(sub, level + 1)
        
        print("\n" + "="*60)
        print("Component Hierarchy")
        print("="*60)
        
        # 只打印根组件
        for comp_data in self.components.values():
            if '.' not in comp_data.name:
                print_component(comp_data, 0)
    
    def to_dataframe(self, component_name: str = None) -> pd.DataFrame:
        """
        转换为DataFrame
        
        Parameters:
        -----------
        component_name : str, optional
            特定组件名，如果为None则返回所有数据
        """
        if component_name:
            comp = self.get_component(component_name)
            if comp:
                return comp.get_dataframe(self.t)
            return pd.DataFrame()
        
        # 合并所有组件的数据
        all_dfs = []
        
        if self.t is not None:
            df = pd.DataFrame({'time': self.t})
            all_dfs.append(df)
        
        for comp_data in self.components.values():
            if '.' not in comp_data.name:  # 只包含根组件
                comp_df = comp_data.get_dataframe(None)
                all_dfs.append(comp_df)
        
        if all_dfs:
            return pd.concat(all_dfs, axis=1)
        return pd.DataFrame()
    
    def to_dict(self, include_hierarchy: bool = True) -> Dict:
        """
        转换为字典格式（JSON友好）
        
        Parameters:
        -----------
        include_hierarchy : bool
            是否包含层次结构
        """
        result = {
            'model_name': self.model.name if hasattr(self.model, 'name') else 'unknown',
            'time': self.t.tolist() if self.t is not None else None,
            'components': {}
        }
        
        for name, comp_data in self.components.items():
            if include_hierarchy:
                # 只包含根组件（子组件会递归包含）
                if '.' not in name:
                    result['components'][name] = comp_data.to_dict(True)
            else:
                result['components'][name] = comp_data.to_dict(False)
        
        return result
    
    def plot_component(self, component_name: str, figsize=(12, 6)):
        """绘制特定组件的所有变量"""
        comp = self.get_component(component_name)
        if not comp:
            print(f"Component '{component_name}' not found")
            return
        
        if self.t is None:
            print("No time data available")
            return
        
        df = comp.get_dataframe(self.t)
        
        # 移除time列
        plot_cols = [col for col in df.columns if col != 'time']
        
        if not plot_cols:
            print(f"No data to plot for component '{component_name}'")
            return
        
        n_vars = len(plot_cols)
        n_cols = min(2, n_vars)
        n_rows = (n_vars + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
        axes = axes.flatten()
        
        for i, col in enumerate(plot_cols):
            axes[i].plot(self.t, df[col].values, linewidth=2)
            axes[i].set_xlabel('Time')
            axes[i].set_ylabel(col)
            axes[i].set_title(col)
            axes[i].grid(True, alpha=0.3)
        
        for i in range(n_vars, len(axes)):
            axes[i].set_visible(False)
        
        plt.suptitle(f'Component: {component_name} ({comp.type})', fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.show()
    
    def plot_comparison(self, component_names: List[str], variable: str = 'energy', 
                       figsize=(10, 6)):
        """比较多个组件的同一变量"""
        if self.t is None:
            print("No time data available")
            return
        
        plt.figure(figsize=figsize)
        
        for comp_name in component_names:
            comp = self.get_component(comp_name)
            if comp:
                df = comp.get_dataframe(self.t)
                matching_cols = [col for col in df.columns if variable in col.lower()]
                
                for col in matching_cols:
                    plt.plot(self.t, df[col].values, label=f'{comp_name}', linewidth=2)
        
        plt.xlabel('Time')
        plt.ylabel(variable.capitalize())
        plt.title(f'{variable.capitalize()} Comparison')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
    
    def summary(self):
        """打印结果摘要"""
        print("\n" + "="*70)
        print("Bond Graph Post-Processing Summary")
        print("="*70)
        
        print(f"\nModel: {self.model.name if hasattr(self.model, 'name') else 'Unknown'}")
        
        if self.t is not None:
            print(f"\nTime range: [{self.t[0]:.4f}, {self.t[-1]:.4f}]")
            print(f"Time steps: {len(self.t)}")
        
        print(f"\nTotal components: {len(self.components)}")
        
        # 按类型统计
        type_counts = defaultdict(int)
        for comp in self.components.values():
            type_counts[comp.type] += 1
        
        print("\nComponents by type:")
        for comp_type, count in sorted(type_counts.items()):
            print(f"  {comp_type}: {count}")
        
        # 统计有状态的组件
        stateful_comps = [c for c in self.components.values() if c.states]
        print(f"\nComponents with states: {len(stateful_comps)}")
        
        print("\n" + "="*70)


# 使用示例
if __name__ == "__main__":
    """
    使用示例:
    
    from BondGraphTools import new, add, connect, simulate
    
    # 1. 创建模型
    model = new(name='RC')
    C = new("C", value=1.0, name="C1")
    R = new("R", value=1.0, name="R1")
    se = new("Se", value=1.0, name="Source")
    one = new("1")
    
    add(model, R, C, one, se)
    connect(se, one)
    connect(R, one)
    connect(C, one)
    
    # 2. 仿真
    timespan = [0, 5]
    x0 = {'x_0': 1}
    t, x = simulate(model, timespan=timespan, x0=x0)
    
    # 3. 后处理
    post = BondGraphPost(model, (t, x))
    
    # 4. 查看组件层次
    post.list_components()
    
    # 5. 查看摘要
    post.summary()
    
    # 6. 获取特定组件
    c1 = post.get_component('C1')
    print(c1)
    print("States:", c1.states)
    print("Ports:", c1.ports)
    print("Energy:", c1.energy)
    
    # 7. 按类型获取组件
    capacitors = post.get_components_by_type('C')
    resistors = post.get_components_by_type('R')
    
    # 8. 获取DataFrame
    df_all = post.to_dataframe()
    df_c1 = post.to_dataframe('C1')
    
    # 9. 转换为字典（用于JSON/UI）
    data_dict = post.to_dict(include_hierarchy=True)
    
    # 10. 绘图
    post.plot_component('C1')
    post.plot_comparison(['C1', 'R1'], variable='power')
    
    # 11. 导出
    import json
    with open('bondgraph_results.json', 'w') as f:
        json.dump(post.to_dict(), f, indent=2)
    """
    print("BondGraphPost类已定义，支持组件层次结构和子模型")