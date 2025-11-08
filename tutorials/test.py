from BondGraphTools import reaction_builder
from BondGraphTools.reaction_builder import Reaction_Network
from BondGraphTools import simulate
from numpy import log, array, linspace, zeros, exp
import matplotlib.pyplot as plt
from sympy import init_printing, SparseMatrix, Eq, Symbol, lambdify, symbols
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy.integrate import solve_ivp
import seaborn as sns

init_printing()

class BiochemicalReactionSolver:
    """生化反应网络完整求解器 - 基于BondGraphTools"""
    
    def __init__(self):
        self.reaction_network = None
        self.bond_graph_model = None
        self.solution_data = {}
        self.full_model_equations = {}
        self.analysis_results = {}
        self.metabolites = []  # 代谢物列表
        
    def initialise_MM_model(self, temperature=310):
        """初始化 Michaelis-Menten 酶反应模型"""
        print("🧬 初始化 Michaelis-Menten 酶反应模型...")
        
        # 创建反应网络
        rn_MM = Reaction_Network(name='Michaelis-Menten enzyme', temperature=temperature)
        
        # 添加反应
        # E + S ⇌ C (酶-底物络合)
        rn_MM.add_reaction('E + S = C', name='R1')
        # C → E + P (产物形成)
        rn_MM.add_reaction('C = E + P', name='R2')
        
        # 添加化学恒定池 (chemostat) - 维持底物和产物浓度
        rn_MM.add_chemostat('S')  # 底物库
        rn_MM.add_chemostat('P')  # 产物库
        
        # 转换为键合图模型
        model = rn_MM.as_network_model()
        
        # 设置动力学参数
        (model/"C:E").set_param('k', 1)   # 酶的化学势参数
        (model/"C:C").set_param('k', 1)   # 络合物的化学势参数
        (model/"R:R1").set_param('r', 1)  # 反应R1的反应速率常数
        (model/"R:R2").set_param('r', 1)  # 反应R2的反应速率常数
        
        # 设置热力学参数
        R = reaction_builder.R  # 气体常数
        T = temperature
        K_S = 1   # 底物的平衡常数
        K_P = 1   # 产物的平衡常数
        x_S = 2   # 底物活度
        x_P = 1   # 产物活度
        
        # 设置化学势
        (model/"SS:S").set_param('e', R*T*log(K_S*x_S))  # 底物化学势
        (model/"SS:P").set_param('e', R*T*log(K_P*x_P))  # 产物化学势
        
        self.reaction_network = rn_MM
        self.bond_graph_model = model
        self.metabolites = list(rn_MM._chemostats.keys())  # 获取代谢物列表
        
        print(f"✓ 模型初始化完成")
        print(f"  - 反应数量: 2 (R1: E+S⇌C, R2: C→E+P)")
        print(f"  - 物种数量: {len(rn_MM.species)} ({rn_MM.species})")
        print(f"  - 代谢物库: {self.metabolites}")
        
        return rn_MM, model
    
    def full_equations(self, model):
        """获取模型的完整方程系统（基于文档中的方法）"""
        # Load full equations of model
        X, mapping, A, F, G = model.system_model()
        # AX + F(X) = 0
        # G(X) = 0
        AX = A * SparseMatrix(X) + F
        full_model_equations = {}
        for i in range(AX.rows):
            xi = X[i]
            eqn = xi - AX[i, 0]
            full_model_equations[str(xi)] = eqn
        
        self.full_model_equations = full_model_equations
        return full_model_equations
    
    def find_port(self, component, direction):
        """找到组件的端口（基于文档中的方法）"""
        if direction in ['f', 'forward']:
            index = 0
        elif direction in ['r', 'reverse']:
            index = 1
        return list(component.ports.keys())[index]
    
    def reaction_flux_expression(self, model, Re_comp, direction):
        """获取反应流量的数学表达式"""
        mapping = model.system_model()[1]
        port = self.find_port(Re_comp, direction)
        bond_index = mapping[1][port]
        
        if not self.full_model_equations:
            self.full_equations(model)
            
        V = self.full_model_equations[f'f_{bond_index}']
        return V
    
    def reaction_affinity_expression(self, model, Re_comp, direction):
        """获取反应亲和势的数学表达式"""
        mapping = model.system_model()[1]
        port = self.find_port(Re_comp, direction)
        bond_index = mapping[1][port]
        
        if not self.full_model_equations:
            self.full_equations(model)
            
        A = self.full_model_equations[f'e_{bond_index}']
        return A
    
    def reaction_affinity(self, model, Re_comp):
        """计算净反应亲和势"""
        Af = self.reaction_affinity_expression(model, Re_comp, 'f')
        Ar = self.reaction_affinity_expression(model, Re_comp, 'r')
        return Af - Ar
    
    def convert_to_function(self, expression, model):
        """将符号表达式转换为Python函数"""
        states = [Symbol(x) for x in model.state_vars.keys()]
        return lambdify(([states]), expression)
    
    def flux_function(self, model, Re_comp):
        """返回可用于计算流量的函数"""
        V = self.reaction_flux_expression(model, Re_comp, 'f')
        return self.convert_to_function(V, model)
    
    def reaction_affinity_function(self, model, Re_comp):
        """返回可用于计算反应亲和势的函数"""
        A = self.reaction_affinity(model, Re_comp)
        return self.convert_to_function(A, model)
    
    def find_species(self, model, species, metabolites):
        """找到对应物种的组件"""
        if species in metabolites:
            return model/f'SS:{species}'
        else:
            return model/f'C:{species}'
    
    def find_species_port(self, component, direction):
        """找到C组件对应的端口"""
        return list(component.ports.keys())[0]
    
    def species_potential(self, model, species):
        """返回物种化学势的符号表达式"""
        comp = self.find_species(model, species, self.metabolites)
        mapping = model.system_model()[1]
        bond_index = mapping[1][self.find_species_port(comp, 0)]
        
        if not self.full_model_equations:
            self.full_equations(model)
            
        potential = self.full_model_equations[f'e_{bond_index}']
        return potential
    
    def species_potential_func(self, model, species):
        """返回可用于计算物种化学势的函数"""
        potential = self.species_potential(model, species)
        return self.convert_to_function(potential, model)
    
    def solve_complete_system(self, timespan=(0., 3.), x0=None, n_points=1000):
        """完整系统求解"""
        if self.bond_graph_model is None:
            print("❌ 错误：请先初始化模型")
            return None
        
        print("🔄 开始求解系统...")
        
        # 默认初始条件
        if x0 is None:
            x0 = [1, 2]  # [E, C]
        
        # 使用BondGraphTools的simulate函数求解状态变量
        t, x = simulate(self.bond_graph_model, timespan=timespan, x0=x0)
        
        # 存储基本解
        self.solution_data = {
            'time': t,
            'states': x,
            'state_names': list(self.bond_graph_model.state_vars.keys())
        }
        
        print(f"✓ 状态变量求解完成")
        print(f"  - 时间点数: {len(t)}")
        print(f"  - 状态变量: {self.solution_data['state_names']}")
        
        return t, x
    
    def calculate_reaction_velocities(self):
        """计算反应速率"""
        if 'time' not in self.solution_data:
            print("❌ 错误：请先求解系统")
            return None
        
        print("🧪 计算反应速率...")
        
        model = self.bond_graph_model
        x = self.solution_data['states']
        
        # 获取反应组件
        R1 = model/"R:R1"
        R2 = model/"R:R2"
        
        # 创建速率函数
        V_R1_func = self.flux_function(model, R1)
        V_R2_func = self.flux_function(model, R2)
        
        # 计算速率时间序列
        V = [[V_R1_func(states), V_R2_func(states)] for states in x]
        V = array(V)
        
        self.analysis_results['reaction_velocities'] = {
            'R1': V[:, 0],
            'R2': V[:, 1],
            'total': V
        }
        
        print(f"✓ 反应速率计算完成")
        return V
    
    def calculate_reaction_affinities(self):
        """计算反应亲和势"""
        if 'time' not in self.solution_data:
            print("❌ 错误：请先求解系统")
            return None
        
        print("⚡ 计算反应亲和势...")
        
        model = self.bond_graph_model
        x = self.solution_data['states']
        
        # 获取反应组件
        R1 = model/"R:R1"
        R2 = model/"R:R2"
        
        # 创建亲和势函数
        A_R1_func = self.reaction_affinity_function(model, R1)
        A_R2_func = self.reaction_affinity_function(model, R2)
        
        # 计算亲和势时间序列
        A = [[A_R1_func(s), A_R2_func(s)] for s in x]
        A = array(A)
        
        self.analysis_results['reaction_affinities'] = {
            'R1': A[:, 0],
            'R2': A[:, 1],
            'total': A
        }
        
        print(f"✓ 反应亲和势计算完成")
        return A
    
    def calculate_power_consumption(self):
        """计算功率消耗"""
        if 'reaction_velocities' not in self.analysis_results:
            self.calculate_reaction_velocities()
        if 'reaction_affinities' not in self.analysis_results:
            self.calculate_reaction_affinities()
        
        print("🔋 计算功率消耗...")
        
        V = self.analysis_results['reaction_velocities']['total']
        A = self.analysis_results['reaction_affinities']['total']
        
        # 功率 = 反应速率 × 亲和势
        power = V * A
        total_power = [sum(p) for p in power]
        
        self.analysis_results['power'] = {
            'R1': power[:, 0],
            'R2': power[:, 1],
            'total': total_power,
            'individual': power
        }
        
        print(f"✓ 功率计算完成")
        return power, total_power
    
    def calculate_chemical_potentials(self):
        """计算化学势"""
        if 'time' not in self.solution_data:
            print("❌ 错误：请先求解系统")
            return None
        
        print("🧲 计算化学势...")
        
        model = self.bond_graph_model
        x = self.solution_data['states']
        
        chemical_potentials = {}
        
        # 对每个物种计算化学势
        for species in self.reaction_network.species:
            try:
                potential_func = self.species_potential_func(model, species)
                # 计算整个仿真过程中的化学势
                cp = [potential_func(s) for s in x]
                # 存储结果
                chemical_potentials[species] = cp
                print(f"  ✓ {species}: 化学势计算完成")
            except Exception as e:
                print(f"  ❌ {species}: 化学势计算失败 - {e}")
        
        self.analysis_results['chemical_potentials'] = chemical_potentials
        
        print(f"✓ 化学势计算完成")
        return chemical_potentials
    
    def generate_comprehensive_analysis(self):
        """生成综合分析"""
        print("📊 生成综合分析...")
        
        # 确保所有分析都已完成
        if 'reaction_velocities' not in self.analysis_results:
            self.calculate_reaction_velocities()
        if 'reaction_affinities' not in self.analysis_results:
            self.calculate_reaction_affinities()
        if 'power' not in self.analysis_results:
            self.calculate_power_consumption()
        if 'chemical_potentials' not in self.analysis_results:
            self.calculate_chemical_potentials()
        
        # 计算稳态值
        t = self.solution_data['time']
        steady_state_analysis = {}
        
        # 状态变量的稳态值
        for i, state_name in enumerate(self.solution_data['state_names']):
            final_value = self.solution_data['states'][-1][i]
            steady_state_analysis[f'steady_state_{state_name}'] = final_value
        
        # 反应速率的稳态值
        V_steady = {
            'R1': self.analysis_results['reaction_velocities']['R1'][-1],
            'R2': self.analysis_results['reaction_velocities']['R2'][-1]
        }
        steady_state_analysis['steady_state_velocities'] = V_steady
        
        # 化学势的稳态值
        cp_steady = {}
        for species, potentials in self.analysis_results['chemical_potentials'].items():
            cp_steady[species] = potentials[-1]
        steady_state_analysis['steady_state_chemical_potentials'] = cp_steady
        
        # 计算总能量变化
        total_power = self.analysis_results['power']['total']
        dt = t[1] - t[0] if len(t) > 1 else 0.01
        cumulative_energy = np.cumsum(total_power) * dt
        
        self.analysis_results['energy_analysis'] = {
            'cumulative_energy_dissipation': cumulative_energy,
            'total_energy_dissipated': cumulative_energy[-1],
            'average_power': np.mean(total_power),
            'peak_power': np.max(np.abs(total_power))
        }
        
        self.analysis_results['steady_state'] = steady_state_analysis
        
        print("✓ 综合分析完成")
        return self.analysis_results
    
    def plot_comprehensive_results(self, figsize=(16, 12)):
        """绘制综合分析结果"""
        if not self.analysis_results:
            self.generate_comprehensive_analysis()
        
        t = self.solution_data['time']
        x = self.solution_data['states']
        
        # 创建综合图表
        fig, axes = plt.subplots(3, 3, figsize=figsize)
        fig.suptitle('生化反应网络完整分析结果 - Michaelis-Menten 酶动力学', 
                    fontsize=16, fontweight='bold')
        
        # 1. 状态变量 (物种浓度)
        ax = axes[0, 0]
        ax.plot(t, x)
        ax.set_xlabel('时间')
        ax.set_ylabel('摩尔数')
        ax.set_title('物种动态 (状态变量)')
        ax.legend(self.solution_data['state_names'])
        ax.grid(True, alpha=0.3)
        
        # 2. 反应速率
        ax = axes[0, 1]
        V = self.analysis_results['reaction_velocities']['total']
        ax.plot(t, V)
        ax.set_xlabel('时间')
        ax.set_ylabel('反应速率')
        ax.set_title('反应速率')
        ax.legend(['R1 (E+S⇌C)', 'R2 (C→E+P)'])
        ax.grid(True, alpha=0.3)
        
        # 3. 反应亲和势
        ax = axes[0, 2]
        A = self.analysis_results['reaction_affinities']['total']
        ax.plot(t, A)
        ax.set_xlabel('时间')
        ax.set_ylabel('反应亲和势')
        ax.set_title('反应亲和势')
        ax.legend(['R1', 'R2'])
        ax.grid(True, alpha=0.3)
        
        # 4. 功率消耗
        ax = axes[1, 0]
        total_power = self.analysis_results['power']['total']
        ax.plot(t, total_power, 'k-', linewidth=2)
        ax.set_xlabel('时间')
        ax.set_ylabel('总功率')
        ax.set_title('系统功率消耗')
        ax.grid(True, alpha=0.3)
        
        # 5. 各反应功率贡献
        ax = axes[1, 1]
        power_individual = self.analysis_results['power']['individual']
        ax.plot(t, power_individual)
        ax.set_xlabel('时间')
        ax.set_ylabel('功率')
        ax.set_title('各反应功率贡献')
        ax.legend(['R1', 'R2'])
        ax.grid(True, alpha=0.3)
        
        # 6. 化学势
        ax = axes[1, 2]
        cp = self.analysis_results['chemical_potentials']
        for species, potentials in cp.items():
            ax.plot(t, potentials, linewidth=2, label=species)
        ax.set_xlabel('时间')
        ax.set_ylabel('化学势')
        ax.set_title('物种化学势')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 7. 累积能量耗散
        ax = axes[2, 0]
        cumulative_energy = self.analysis_results['energy_analysis']['cumulative_energy_dissipation']
        ax.plot(t, cumulative_energy, 'r-', linewidth=2)
        ax.set_xlabel('时间')
        ax.set_ylabel('累积能量')
        ax.set_title('累积能量耗散')
        ax.grid(True, alpha=0.3)
        
        # 8. 相位图 (如果有两个状态变量)
        ax = axes[2, 1]
        if len(self.solution_data['state_names']) >= 2:
            ax.plot(x[:, 0], x[:, 1], 'b-', linewidth=2)
            ax.set_xlabel(self.solution_data['state_names'][0])
            ax.set_ylabel(self.solution_data['state_names'][1])
            ax.set_title('相位图')
            ax.grid(True, alpha=0.3)
        
        # 9. 稳态分析
        ax = axes[2, 2]
        steady_state = self.analysis_results['steady_state']
        ss_velocities = steady_state['steady_state_velocities']
        
        reactions = list(ss_velocities.keys())
        velocities = list(ss_velocities.values())
        
        bars = ax.bar(reactions, velocities, alpha=0.7, color=['skyblue', 'orange'])
        ax.set_ylabel('稳态反应速率')
        ax.set_title('稳态反应速率')
        ax.grid(True, alpha=0.3)
        
        # 添加数值标签
        for bar, velocity in zip(bars, velocities):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{velocity:.4f}', ha='center', va='bottom')
        
        plt.tight_layout()
        plt.show()
        
        return fig
    
    def export_results(self, filename_prefix="michaelis_menten"):
        """导出分析结果"""
        if not self.analysis_results:
            self.generate_comprehensive_analysis()
        
        print("💾 导出分析结果...")
        
        # 创建综合数据DataFrame
        data_dict = {
            'time': self.solution_data['time']
        }
        
        # 添加状态变量
        for i, state_name in enumerate(self.solution_data['state_names']):
            data_dict[f'state_{state_name}'] = self.solution_data['states'][:, i]
        
        # 添加反应速率
        data_dict['velocity_R1'] = self.analysis_results['reaction_velocities']['R1']
        data_dict['velocity_R2'] = self.analysis_results['reaction_velocities']['R2']
        
        # 添加亲和势
        data_dict['affinity_R1'] = self.analysis_results['reaction_affinities']['R1']
        data_dict['affinity_R2'] = self.analysis_results['reaction_affinities']['R2']
        
        # 添加功率
        data_dict['power_R1'] = self.analysis_results['power']['R1']
        data_dict['power_R2'] = self.analysis_results['power']['R2']
        data_dict['total_power'] = self.analysis_results['power']['total']
        
        # 添加化学势
        for species, potentials in self.analysis_results['chemical_potentials'].items():
            data_dict[f'chemical_potential_{species}'] = potentials
        
        # 添加累积能量
        data_dict['cumulative_energy'] = self.analysis_results['energy_analysis']['cumulative_energy_dissipation']
        
        # 创建和导出DataFrame
        df = pd.DataFrame(data_dict)
        csv_filename = f"{filename_prefix}_complete_analysis.csv"
        df.to_csv(csv_filename, index=False)
        
        print(f"✓ 完整数据已导出到: {csv_filename}")
        
        # 导出摘要报告
        summary = {
            'model_info': {
                'name': 'Michaelis-Menten Enzyme Kinetics',
                'species': list(self.reaction_network.species),
                'reactions': ['R1: E + S ⇌ C', 'R2: C → E + P'],
                'metabolites': self.metabolites
            },
            'simulation_info': {
                'time_points': len(self.solution_data['time']),
                'time_range': [float(self.solution_data['time'][0]), 
                             float(self.solution_data['time'][-1])],
                'state_variables': self.solution_data['state_names']
            },
            'steady_state_analysis': self.analysis_results['steady_state'],
            'energy_analysis': self.analysis_results['energy_analysis']
        }
        
        import json
        json_filename = f"{filename_prefix}_analysis_summary.json"
        with open(json_filename, 'w', encoding='utf-8') as f:
            json.dump(summary, f, indent=2, ensure_ascii=False, default=str)
        
        print(f"✓ 分析摘要已导出到: {json_filename}")
        
        return df, summary

# 使用示例和演示
def run_complete_analysis():
    """运行完整的分析流程"""
    print("🧬 生化反应网络完整分析 - Michaelis-Menten 酶动力学")
    print("=" * 70)
    
    # 创建求解器实例
    solver = BiochemicalReactionSolver()
    
    # 初始化模型
    rn, model = solver.initialise_MM_model(temperature=310)
    
    # 求解系统
    t, x = solver.solve_complete_system(
        timespan=(0., 3.), 
        x0=[1, 2]  # 初始条件：E=1, C=2
    )
    
    # 生成综合分析
    analysis = solver.generate_comprehensive_analysis()
    
    # 打印关键结果
    print("\n📊 关键结果摘要:")
    print("-" * 50)
    
    # 稳态值
    ss = analysis['steady_state']
    print("🔄 稳态值:")
    for state in solver.solution_data['state_names']:
        print(f"  {state}: {ss[f'steady_state_{state}']:.6f}")
    
    print("\n⚡ 稳态反应速率:")
    for reaction, velocity in ss['steady_state_velocities'].items():
        print(f"  {reaction}: {velocity:.6f}")
    
    print("\n🧲 稳态化学势:")
    for species, potential in ss['steady_state_chemical_potentials'].items():
        print(f"  {species}: {potential:.6f}")
    
    print(f"\n🔋 能量分析:")
    energy = analysis['energy_analysis']
    print(f"  总能量耗散: {energy['total_energy_dissipated']:.6f}")
    print(f"  平均功率: {energy['average_power']:.6f}")
    print(f"  峰值功率: {energy['peak_power']:.6f}")
    
    # 可视化结果
    print("\n📈 生成可视化结果...")
    solver.plot_comprehensive_results()
    
    # 导出结果
    print("\n💾 导出分析数据...")
    df, summary = solver.export_results()
    
    print("\n✅ 完整分析完成！")
    print("📁 查看生成的 CSV 和 JSON 文件获取详细结果")
    
    return solver, analysis

# 高级分析示例
def advanced_parameter_study():
    """高级参数研究"""
    print("\n🔬 高级参数研究")
    print("=" * 40)
    
    # 研究不同初始条件的影响
    initial_conditions = [
        [0.5, 1.0], [1.0, 2.0], [1.5, 1.5], [2.0, 1.0]
    ]
    
    results_comparison = {}
    
    for i, x0 in enumerate(initial_conditions):
        print(f"\n分析初始条件 {i+1}: E={x0[0]}, C={x0[1]}")
        
        solver = BiochemicalReactionSolver()
        rn, model = solver.initialise_MM_model()
        
        t, x = solver.solve_complete_system(timespan=(0., 5.), x0=x0)
        analysis = solver.generate_comprehensive_analysis()
        
        # 存储关键指标
        results_comparison[f'condition_{i+1}'] = {
            'initial': x0,
            'final_E': analysis['steady_state']['steady_state_E'],
            'final_C': analysis['steady_state']['steady_state_C'],
            'total_energy_dissipated': analysis['energy_analysis']['total_energy_dissipated'],
            'average_power': analysis['energy_analysis']['average_power']
        }
    
    # 比较结果
    print("\n📊 参数研究结果比较:")
    comparison_df = pd.DataFrame(results_comparison).T
    print(comparison_df)
    
    return results_comparison

if __name__ == "__main__":
    # 运行完整分析
    solver, analysis = run_complete_analysis()
    
    # 取消注释以下行运行高级分析
    # advanced_parameter_study()