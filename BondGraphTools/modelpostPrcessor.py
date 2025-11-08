import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union, Any
import matplotlib.pyplot as plt
from sympy import Symbol, lambdify, Matrix, SparseMatrix
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
        
        # 子组件
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
            for state_name, q in self.states.items():
                if 'q' in state_name or state_name == 'x':
                    C = self.parameters['value']
                    self.energy = 0.5 * q**2 / C if C != 0 else None
                    break
        
        elif self.type == 'I' and self.parameters.get('value'):
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
        """转换为字典格式"""
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
        
        for state_name, values in self.states.items():
            data[f'{self.name}_state_{state_name}'] = values
        
        for state_name, values in self.state_derivatives.items():
            data[f'{self.name}_d{state_name}_dt'] = values
        
        for port_id, port_vars in self.ports.items():
            for var_type, values in port_vars.items():
                data[f'{self.name}_port{port_id}_{var_type}'] = values
        
        if self.energy is not None:
            data[f'{self.name}_energy'] = self.energy
        
        for port_id, values in self.power.items():
            data[f'{self.name}_port{port_id}_power'] = values
        
        return pd.DataFrame(data)
    
    def __repr__(self):
        return f"ComponentData(name='{self.name}', type='{self.type}', states={len(self.states)}, ports={len(self.ports)})"


class BondGraphPost:
    """
    Bond Graph模型的层次化后处理类 - 递归扁平化版本
    
    核心思路：
    1. 递归调用所有子模型的 system_model()
    2. 将所有状态变量、端口映射展开到扁平化结构
    3. 构建完整的代数方程求解所有变量
    """
    
    def __init__(self, model, simulation_results: Union[Tuple, np.ndarray] = None):
        """
        初始化后处理类
        
        Parameters:
        -----------
        model : BondGraph模型对象
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
        
        # 组件数据字典
        self.components = {}
        
        # 扁平化的映射关系
        self.flat_state_map = {}  # {(comp_obj, state_name): (full_path, idx)}
        self.flat_port_map = {}   # {(comp_obj, port_obj): (full_path, bond_idx)}
        
        # 递归扁平化所有子模型
        self._flatten_model()
        
        # 处理组件数据
        self._process_components()
        
    def _get_component_info(self, comp) -> Tuple[str, str]:
        """获取组件名称和类型"""
        name = comp.name.split(':')[-1] if hasattr(comp, 'name') else 'unknown'
        comp_type = comp.metamodel if hasattr(comp, 'metamodel') else \
                   (comp.component_type if hasattr(comp, 'component_type') else 'unknown')
        return name, comp_type
    
    def _flatten_model(self, model=None, parent_path: str = ""):
        """
        递归扁平化模型结构
        
        对每个子模型调用 system_model()，提取其内部的状态和端口映射
        """
        if model is None:
            model = self.model
        
        print(f"\n扁平化模型: {parent_path if parent_path else 'root'}")
        
        try:
            # 获取模型的 system_model
            X, mapping, A, F, G = model.system_model()
            
            print(f"  状态变量数: {len(X)}")
            print(f"  映射关系: {len(mapping) if mapping else 0}")
            
            # 处理状态变量映射 mapping[0]
            if mapping and len(mapping) > 0:
                state_mapping = mapping[0]
                print(f"  状态映射项数: {len(state_mapping)}")
                
                for (comp, state_name), idx in state_mapping.items():
                    comp_name, comp_type = self._get_component_info(comp)
                    full_path = f"{parent_path}.{comp_name}" if parent_path else comp_name
                    
                    # 保存到扁平化映射
                    self.flat_state_map[(comp, state_name)] = (full_path, idx)
                    print(f"    状态: {full_path}.{state_name} -> idx {idx}")
            
            # 处理端口映射 mapping[1]
            if mapping and len(mapping) > 1:
                port_mapping = mapping[1]
                print(f"  端口映射项数: {len(port_mapping)}")
                
                for port, bond_idx in port_mapping.items():
                    comp = port.component
                    comp_name, comp_type = self._get_component_info(comp)
                    full_path = f"{parent_path}.{comp_name}" if parent_path else comp_name
                    
                    # 保存到扁平化映射
                    self.flat_port_map[(comp, port)] = (full_path, bond_idx)
                    print(f"    端口: {full_path}.port{port.index if hasattr(port, 'index') else 0} -> bond {bond_idx}")
            
        except Exception as e:
            print(f"  警告: 无法获取 system_model: {e}")
        
        # 递归处理子模型
        if hasattr(model, 'components'):
            for comp in model.components:
                comp_name, comp_type = self._get_component_info(comp)
                
                # 如果组件本身是一个模型（有 components 和 system_model）
                if hasattr(comp, 'components') and hasattr(comp, 'system_model'):
                    full_path = f"{parent_path}.{comp_name}" if parent_path else comp_name
                    print(f"  发现子模型: {comp_name}")
                    self._flatten_model(comp, full_path)
    
    def _extract_all_components(self, model=None, parent_path: str = "") -> Dict[str, Any]:
        """递归提取所有组件"""
        if model is None:
            model = self.model
        
        components = {}
        
        if hasattr(model, 'components'):
            for comp in model.components:
                comp_name, comp_type = self._get_component_info(comp)
                full_path = f"{parent_path}.{comp_name}" if parent_path else comp_name
                components[full_path] = comp
                
                # 递归子模型
                if hasattr(comp, 'components'):
                    sub_components = self._extract_all_components(comp, full_path)
                    components.update(sub_components)
        
        return components
    
    def _process_components(self):
        """
        使用扁平化的映射关系处理组件数据
        """
        if self.x is None:
            return
        
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
        
        print(f"\n处理组件数据:")
        print(f"  仿真时间步数: {len(x_2d)}")
        print(f"  状态变量数: {x_2d.shape[1]}")
        
        # 提取所有组件
        all_components = self._extract_all_components()
        
        # 创建 ComponentData
        for comp_path, comp in all_components.items():
            comp_name, comp_type = self._get_component_info(comp)
            comp_data = ComponentData(comp_path, comp_type, comp)
            
            # 添加参数
            if hasattr(comp, 'value'):
                try:
                    comp_data.add_parameter('value', float(comp.value))
                except:
                    pass
            
            self.components[comp_path] = comp_data
        
        print(f"  提取了 {len(self.components)} 个组件")
        
        # 使用扁平化映射填充状态变量
        print(f"\n填充状态变量:")
        for (comp, state_name), (full_path, idx) in self.flat_state_map.items():
            if full_path in self.components:
                comp_data = self.components[full_path]
                
                if idx < x_2d.shape[1]:
                    state_values = x_2d[:, idx]
                    deriv_values = dx_dt[:, idx]
                    
                    comp_data.add_state(state_name, state_values, deriv_values)
                    comp_data.calculate_energy()
                    
                    print(f"  {full_path}.{state_name} <- x[{idx}]")
        
        # 计算端口变量
        self._compute_port_variables_flat(x_2d)
        
        # 组织层次结构
        self._organize_hierarchy()
    
    def _compute_port_variables_flat(self, x_2d: np.ndarray):
        """
        使用扁平化映射计算端口变量
        
        基本思路：
        对每个子模型，求解其代数方程 G * [e, f] = A * x + F
        """
        print(f"\n计算端口变量:")
        
        # 为每个有端口的组件初始化端口数组
        for (comp, port), (full_path, bond_idx) in self.flat_port_map.items():
            if full_path in self.components:
                comp_data = self.components[full_path]
                port_id = str(port.index) if hasattr(port, 'index') else '0'
                
                if port_id not in comp_data.ports:
                    comp_data.ports[port_id]['effort'] = np.zeros(len(x_2d))
                    comp_data.ports[port_id]['flow'] = np.zeros(len(x_2d))
        
        # 对主模型求解
        self._solve_algebraic_equations(self.model, x_2d, "")
        
        # 对每个子模型递归求解
        self._solve_submodels_algebraic(self.model, x_2d, "")
        
        # 计算功率
        for comp_data in self.components.values():
            comp_data.calculate_power()
    
    def _solve_submodels_algebraic(self, model, x_2d: np.ndarray, parent_path: str):
        """递归求解子模型的代数方程"""
        if not hasattr(model, 'components'):
            return
        
        for comp in model.components:
            comp_name, comp_type = self._get_component_info(comp)
            full_path = f"{parent_path}.{comp_name}" if parent_path else comp_name
            
            # 如果是子模型
            if hasattr(comp, 'system_model') and hasattr(comp, 'components'):
                print(f"  求解子模型: {full_path}")
                self._solve_algebraic_equations(comp, x_2d, full_path)
                
                # 递归
                self._solve_submodels_algebraic(comp, x_2d, full_path)
    
    def _solve_algebraic_equations(self, model, x_2d: np.ndarray, model_path: str):
        """
        求解特定模型的代数方程
        
        G * [e, f] = A * x + F
        """
        try:
            X, mapping, A, F, G = model.system_model()
            
            if not mapping or len(mapping) < 2:
                return
            
            port_mapping = mapping[1]
            
            if G.rows == 0 or G.cols == 0:
                print(f"    跳过 {model_path}: G 矩阵为空")
                return
            
            print(f"    G 矩阵: {G.rows}x{G.cols}, 状态数: {len(X)}")
            
            # 对每个时间步求解
            for i in range(len(x_2d)):
                x_val = x_2d[i, :len(X)]
                
                try:
                    # 计算 A * x + F
                    AX_F = A * Matrix(x_val) + F
                    
                    # 求解 G * y = AX_F
                    if G.is_square:
                        try:
                            y_solution = G.inv() * AX_F
                        except:
                            # 奇异矩阵，使用伪逆
                            y_solution = G.pinv() * AX_F
                    else:
                        # 非方阵，使用最小二乘
                        y_solution = G.pinv() * AX_F
                    
                    # 分配到各个端口
                    for port, bond_idx in port_mapping.items():
                        comp = port.component
                        comp_name = self._get_component_info(comp)[0]
                        
                        # 构建完整路径
                        if model_path:
                            full_path = f"{model_path}.{comp_name}"
                        else:
                            full_path = comp_name
                        
                        if full_path in self.components:
                            comp_data = self.components[full_path]
                            port_id = str(port.index) if hasattr(port, 'index') else '0'
                            
                            # effort 和 flow 在 y_solution 中
                            # 通常格式: [e_0, f_0, e_1, f_1, ...]
                            if 2 * bond_idx + 1 < len(y_solution):
                                try:
                                    e_val = float(y_solution[2 * bond_idx])
                                    f_val = float(y_solution[2 * bond_idx + 1])
                                    
                                    comp_data.ports[port_id]['effort'][i] = e_val
                                    comp_data.ports[port_id]['flow'][i] = f_val
                                    
                                    if i == 0:
                                        print(f"      {full_path}.port{port_id}: e={e_val:.4f}, f={f_val:.4f}")
                                except Exception as e:
                                    pass
                
                except Exception as e:
                    if i == 0:
                        print(f"    警告: 时间步 {i} 求解失败: {e}")
                    continue
        
        except Exception as e:
            print(f"    警告: 无法求解 {model_path}: {e}")
    
    def _organize_hierarchy(self):
        """组织组件的层次结构"""
        for path, comp_data in list(self.components.items()):
            if '.' in path:
                parts = path.split('.')
                parent_path = '.'.join(parts[:-1])
                if parent_path in self.components:
                    self.components[parent_path].add_subcomponent(comp_data)
    
    def get_component(self, component_name: str) -> Optional[ComponentData]:
        """获取特定组件的数据"""
        # 精确匹配
        if component_name in self.components:
            return self.components[component_name]
        
        # 模糊匹配
        for path, comp_data in self.components.items():
            if path.endswith(component_name):
                return comp_data
        
        # 更模糊的匹配
        for path, comp_data in self.components.items():
            if component_name in path:
                return comp_data
        
        return None
    
    def get_all_components(self, include_junctions: bool = False) -> Dict[str, ComponentData]:
        """获取所有组件"""
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
                port_info = []
                for pid, pdata in comp.ports.items():
                    has_e = 'effort' in pdata and np.any(pdata['effort'] != 0)
                    has_f = 'flow' in pdata and np.any(pdata['flow'] != 0)
                    port_info.append(f"p{pid}({'e' if has_e else ''}{'f' if has_f else ''})")
                print(f"{prefix}│  Ports: {port_info}")
            if comp.subcomponents:
                print(f"{prefix}│  Subcomponents:")
                for sub in comp.subcomponents.values():
                    print_component(sub, level + 1)
        
        print("\n" + "="*60)
        print("Component Hierarchy")
        print("="*60)
        
        for comp_data in self.components.values():
            if '.' not in comp_data.name:
                print_component(comp_data, 0)
    
    def to_dataframe(self, component_name: str = None) -> pd.DataFrame:
        """转换为DataFrame"""
        if component_name:
            comp = self.get_component(component_name)
            if comp:
                return comp.get_dataframe(self.t)
            return pd.DataFrame()
        
        all_dfs = []
        
        if self.t is not None:
            df = pd.DataFrame({'time': self.t})
            all_dfs.append(df)
        
        for comp_data in self.components.values():
            if '.' not in comp_data.name:
                comp_df = comp_data.get_dataframe(None)
                all_dfs.append(comp_df)
        
        if all_dfs:
            return pd.concat(all_dfs, axis=1)
        return pd.DataFrame()
    
    def to_dict(self, include_hierarchy: bool = True) -> Dict:
        """转换为字典格式（JSON友好）"""
        result = {
            'model_name': self.model.name if hasattr(self.model, 'name') else 'unknown',
            'time': self.t.tolist() if self.t is not None else None,
            'components': {}
        }
        
        for name, comp_data in self.components.items():
            if include_hierarchy:
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
            print(f"Available: {list(self.components.keys())}")
            return
        
        if self.t is None:
            print("No time data available")
            return
        
        df = comp.get_dataframe(self.t)
        plot_cols = [col for col in df.columns if col != 'time']
        
        if not plot_cols:
            print(f"No data to plot for '{component_name}'")
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
        
        type_counts = defaultdict(int)
        for comp in self.components.values():
            type_counts[comp.type] += 1
        
        print("\nComponents by type:")
        for comp_type, count in sorted(type_counts.items()):
            print(f"  {comp_type}: {count}")
        
        stateful_comps = [c for c in self.components.values() if c.states]
        print(f"\nComponents with states: {len(stateful_comps)}")
        
        port_comps = [c for c in self.components.values() if c.ports]
        print(f"Components with port data: {len(port_comps)}")
        
        print("\n" + "="*70)


# 使用示例
if __name__ == "__main__":
    print("BondGraphPost - 递归扁平化版本")
    print("\n核心改进:")
    print("1. _flatten_model() - 递归调用所有子模型的 system_model()")
    print("2. _solve_submodels_algebraic() - 递归求解每个子模型的代数方程")
    print("3. 完整的端口变量计算，支持嵌套子模型")
    print("\n使用:")
    print("post = BondGraphPost(mainmodel, (t, x))")
    print("post.list_components()")
    print("post.plot_component('subR.R1')  # 现在应该有端口数据了")