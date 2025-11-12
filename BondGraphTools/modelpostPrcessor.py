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
    Bond Graph模型的层次化后处理类 - 整合版本
    
    核心思路：
    1. 递归处理所有层次的组件
    2. 对每个层次单独调用 system_model()
    3. 使用正确的方法从 AX+F 中提取端口变量
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
        
        # 提取所有组件（包括子模型内的）
        self._extract_all_components()
        
        # 处理组件数据
        if self.x is not None:
            self._process_components()
        
    def _get_component_info(self, comp) -> Tuple[str, str]:
        """获取组件名称和类型"""
        name = comp.name.split(':')[-1] if hasattr(comp, 'name') else 'unknown'
        comp_type = comp.metamodel if hasattr(comp, 'metamodel') else \
                   (comp.component_type if hasattr(comp, 'component_type') else 'unknown')
        return name, comp_type
    
    def _extract_all_components(self, model=None, parent_path: str = ""):
        """递归提取所有组件"""
        if model is None:
            model = self.model
        
        print(f"\n提取组件: {parent_path if parent_path else 'root'}")
        
        if hasattr(model, 'components'):
            for comp in model.components:
                comp_name, comp_type = self._get_component_info(comp)
                full_path = f"{parent_path}.{comp_name}" if parent_path else comp_name
                
                print(f"  发现组件: {full_path} ({comp_type})")
                
                # 创建 ComponentData
                comp_data = ComponentData(full_path, comp_type, comp)
                
                # 添加参数
                if hasattr(comp, 'value'):
                    try:
                        comp_data.add_parameter('value', float(comp.value))
                        print(f"    参数: value = {comp.value}")
                    except:
                        pass
                
                self.components[full_path] = comp_data
                
                # 递归处理子模型
                if hasattr(comp, 'components') and hasattr(comp, 'system_model'):
                    print(f"  {full_path} 是子模型，递归处理...")
                    self._extract_all_components(comp, full_path)
    
    def _process_components(self):
        """处理组件数据：状态变量和端口变量"""
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
        
        # 处理主模型
        self._process_model_level(self.model, x_2d, dx_dt, "")
        
        # 组织层次结构
        self._organize_hierarchy()
    
    def _process_model_level(self, model, x_2d: np.ndarray, dx_dt: np.ndarray, parent_path: str):
        """
        处理特定层次的模型
        
        关键：对每个模型单独调用 system_model()，使用其返回的映射关系
        """
        print(f"\n处理模型层次: {parent_path if parent_path else 'root'}")
        
        try:
            X, mapping, A, F, G = model.system_model()
            
            if not mapping or len(mapping) < 1:
                print(f"  警告: 没有映射关系")
                return
            
            state_mapping = mapping[0]
            port_mapping = mapping[1] if len(mapping) > 1 else {}
            
            print(f"  状态变量数: {len(X)}")
            print(f"  状态映射数: {len(state_mapping)}")
            print(f"  端口映射数: {len(port_mapping)}")
            
            # 处理状态变量
            print(f"  处理状态映射:")
            for (comp, state_name), idx in state_mapping.items():
                comp_name, comp_type = self._get_component_info(comp)
                full_path = f"{parent_path}.{comp_name}" if parent_path else comp_name
                
                # 尝试匹配组件
                matched_path = None
                if full_path in self.components:
                    matched_path = full_path
                else:
                    # 通过组件对象匹配
                    for path, comp_data in self.components.items():
                        if comp_data.component is comp:
                            matched_path = path
                            break
                
                if matched_path and idx < x_2d.shape[1]:
                    comp_data = self.components[matched_path]
                    state_values = x_2d[:, idx]
                    deriv_values = dx_dt[:, idx]
                    
                    comp_data.add_state(state_name, state_values, deriv_values)
                    comp_data.calculate_energy()
                    
                    print(f"    ✓ {matched_path}.{state_name} <- x[{idx}]")
            
            # 求解端口变量 - 使用正确的方法
            self._solve_port_variables_from_equations(model, x_2d, dx_dt, parent_path, X, mapping, A, F)
            
            # 递归处理子模型
            if hasattr(model, 'components'):
                for comp in model.components:
                    comp_name, comp_type = self._get_component_info(comp)
                    full_path = f"{parent_path}.{comp_name}" if parent_path else comp_name
                    
                    # 如果是子模型
                    if hasattr(comp, 'components') and hasattr(comp, 'system_model'):
                        print(f"  发现子模型: {full_path}")
                        self._process_model_level(comp, x_2d, dx_dt, full_path)
        
        except Exception as e:
            print(f"  错误: {e}")
            import traceback
            traceback.print_exc()
    
    def _solve_port_variables_from_equations(self, model, x_2d: np.ndarray, dx_dt: np.ndarray,
                                            model_path: str, X, mapping, A, F):
        """
        从完整方程 AX + F 中提取端口变量
        
        核心思路（基于第一个正确版本）：
        1. X 包含所有变量: [dx_0, e_0, f_0, e_1, f_1, ..., x_0]
        2. 通过 AX + F 计算出所有变量的值
        3. 根据符号名称匹配提取 effort 和 flow
        """
        if not mapping or len(mapping) < 2:
            print(f"    没有端口映射")
            return
        
        state_mapping = mapping[0]
        port_mapping = mapping[1]
        
        if not port_mapping:
            print(f"    没有端口需要处理")
            return
        
        print(f"    从方程提取端口变量:")
        print(f"    X 向量长度: {len(X)}")
        
        # 初始化端口数组
        for port, bond_idx in port_mapping.items():
            comp = port.component
            comp_name = self._get_component_info(comp)[0]
            full_path = f"{model_path}.{comp_name}" if model_path else comp_name
            
            # 匹配组件
            matched_path = None
            if full_path in self.components:
                matched_path = full_path
            else:
                for path, comp_data in self.components.items():
                    if comp_data.component is comp:
                        matched_path = path
                        break
            
            if matched_path:
                comp_data = self.components[matched_path]
                port_id = str(port.index) if hasattr(port, 'index') else '0'
                
                if port_id not in comp_data.ports:
                    comp_data.ports[port_id]['effort'] = np.zeros(len(x_2d))
                    comp_data.ports[port_id]['flow'] = np.zeros(len(x_2d))
        
        # 对每个时间步计算
        for i in range(len(x_2d)):
            # 准备输入：状态变量 + 状态导数
            # 需要按照 X 中符号的顺序准备输入
            input_values = []
            
            for sym in X:
                sym_str = str(sym)
                
                if sym_str.startswith('dx_'):
                    # 状态导数
                    idx = int(sym_str.split('_')[1])
                    if idx < dx_dt.shape[1]:
                        input_values.append(dx_dt[i, idx])
                    else:
                        input_values.append(0.0)
                
                elif sym_str.startswith('x_'):
                    # 状态变量
                    idx = int(sym_str.split('_')[1])
                    if idx < x_2d.shape[1]:
                        input_values.append(x_2d[i, idx])
                    else:
                        input_values.append(0.0)
                
                else:
                    # 其他变量（e, f 等）暂时设为 0，后面会计算
                    input_values.append(0.0)
            
            try:
                # 计算 AX + F，得到所有变量的值
                X_computed = A * Matrix(input_values) + F
                
                # 从计算结果中提取端口变量
                for port, bond_idx in port_mapping.items():
                    comp = port.component
                    
                    # 匹配组件
                    matched_path = None
                    comp_name = self._get_component_info(comp)[0]
                    full_path = f"{model_path}.{comp_name}" if model_path else comp_name
                    
                    if full_path in self.components:
                        matched_path = full_path
                    else:
                        for path, comp_data in self.components.items():
                            if comp_data.component is comp:
                                matched_path = path
                                break
                    
                    if matched_path:
                        comp_data = self.components[matched_path]
                        port_id = str(port.index) if hasattr(port, 'index') else '0'
                        
                        # 在 X 向量中查找对应的 e 和 f 符号
                        e_symbol = f"e_{bond_idx}"
                        f_symbol = f"f_{bond_idx}"
                        
                        for idx, var in enumerate(X):
                            var_str = str(var)
                            
                            if var_str == e_symbol and idx < len(X_computed):
                                try:
                                    e_val = float(X_computed[idx])
                                    comp_data.ports[port_id]['effort'][i] = e_val
                                    
                                    if i == 0:
                                        print(f"      {matched_path}.port{port_id}.effort <- {e_symbol} = {e_val:.4f}")
                                except:
                                    pass
                            
                            elif var_str == f_symbol and idx < len(X_computed):
                                try:
                                    f_val = float(X_computed[idx])
                                    comp_data.ports[port_id]['flow'][i] = f_val
                                    
                                    if i == 0:
                                        print(f"      {matched_path}.port{port_id}.flow <- {f_symbol} = {f_val:.4f}")
                                except:
                                    pass
            
            except Exception as e:
                if i == 0:
                    print(f"    警告: 时间步 {i} 计算失败: {e}")
                continue
        
        # 计算功率
        for comp_data in self.components.values():
            if model_path in comp_data.name or (not model_path and '.' not in comp_data.name):
                comp_data.calculate_power()
    
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
        
        # 模糊匹配（末尾匹配）
        for path, comp_data in self.components.items():
            if path.endswith(component_name):
                return comp_data
        
        # 更模糊的匹配（包含）
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
            
            if comp.parameters:
                params_str = ", ".join([f"{k}={v}" for k, v in comp.parameters.items()])
                print(f"{prefix}│  参数: {params_str}")
            
            if comp.states:
                print(f"{prefix}│  状态: {list(comp.states.keys())}")
            
            if comp.ports:
                port_info = []
                for pid, pdata in comp.ports.items():
                    has_e = 'effort' in pdata and np.any(pdata['effort'] != 0)
                    has_f = 'flow' in pdata and np.any(pdata['flow'] != 0)
                    if has_e or has_f:
                        port_info.append(f"p{pid}({'e' if has_e else ''}{'f' if has_f else ''})")
                if port_info:
                    print(f"{prefix}│  端口: {port_info}")
            
            if comp.subcomponents:
                print(f"{prefix}│  子组件:")
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
    
    def plot_comparison(self, component_names: List[str], variable: str = 'energy', unit: str = '',
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
        plt.ylabel(variable.capitalize()+(f' ({unit})' if unit else ''))
        plt.title(f'{variable.capitalize()} Comparison')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

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
    print("BondGraphPost - 整合层次化与正确求解版本")
    print("\n核心特点:")
    print("1. 层次化处理：递归处理 mainmodel 和所有 submodel")
    print("2. 正确求解：使用 AX+F 方程从符号匹配提取端口变量")
    print("3. 完整映射：通过组件对象匹配，确保路径正确")
    print("\n使用:")
    print("post = BondGraphPost(mainmodel, (t, x))")
    print("post.list_components()")
    print("post.plot_component('subR.C2')")
    print("post.plot_component('subR.R1')")