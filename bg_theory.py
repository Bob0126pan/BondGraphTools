"""
BondGraphTools 核心流程演示

展示完整流程：
1. 模块创建 (Component Creation)
2. 模块连接 (Connection)
3. 模型解析 (Model Analysis)
4. 方程转换 (Equation Generation)
5. 数值求解 (Numerical Solving)

示例：简单RC电路 (也可以理解为水箱-阀门系统)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import sympy as sp

print("="*70)
print("BondGraphTools 工作流程演示")
print("="*70)

# 尝试导入BondGraphTools
try:
    import BondGraphTools as bgt
    HAS_BGT = True
    print("✅ 使用真实的 BondGraphTools")
except ImportError:
    HAS_BGT = False
    print("⚠️  BondGraphTools 未安装，使用简化实现展示原理")
    print("   安装命令: pip install bondgraphtools\n")


# ==================== 简化实现（理解BGT原理） ====================

class SimpleBondGraph:
    """
    简化的BondGraph实现 - 展示核心概念
    
    核心要素：
    1. Component: 组件（R, C, I, Se, Sf等）
    2. Bond: 能量键，连接组件
    3. Port: 端口，有effort(e)和flow(f)两个变量
    4. Junction: 结点（0-junction, 1-junction）
    """
    
    def __init__(self, name):
        self.name = name
        self.components = []
        self.bonds = []
        self.state_vars = []  # 状态变量列表
        self.equations = []   # 方程列表（符号形式）
        
        print(f"\n📦 创建BondGraph模型: {name}")
    
    def add_C_element(self, name, C_value, initial_state=0.0):
        """
        添加容性元件 (C-element)
        
        本构关系：
        - q = C * e    (q: 广义位移, e: 势)
        - dq/dt = f    (f: 流)
        
        物理类比：
        - 电容：q=电荷, e=电压, f=电流
        - 水箱：q=体积, e=压力, f=流量
        """
        component = {
            'name': name,
            'type': 'C',
            'value': C_value,
            'state': sp.Symbol(f'q_{name}'),  # 状态变量
            'effort': sp.Symbol(f'e_{name}'),
            'flow': sp.Symbol(f'f_{name}'),
            'initial_state': initial_state
        }
        
        # 添加本构关系
        # e = q / C
        component['constitutive'] = sp.Eq(
            component['effort'],
            component['state'] / C_value
        )
        
        self.components.append(component)
        self.state_vars.append(component['state'])
        
        print(f"  ✓ 添加C元件: {name}, C={C_value}, q0={initial_state}")
        print(f"    本构关系: {component['constitutive']}")
        
        return component
    
    def add_R_element(self, name, R_value):
        """
        添加阻性元件 (R-element)
        
        本构关系：
        - e = R * f   (线性阻力)
        
        物理类比：
        - 电阻：e=电压, f=电流
        - 阻尼器：e=力, f=速度
        - 管道：e=压差, f=流量
        """
        component = {
            'name': name,
            'type': 'R',
            'value': R_value,
            'effort': sp.Symbol(f'e_{name}'),
            'flow': sp.Symbol(f'f_{name}')
        }
        
        # 本构关系: e = R * f
        component['constitutive'] = sp.Eq(
            component['effort'],
            R_value * component['flow']
        )
        
        self.components.append(component)
        
        print(f"  ✓ 添加R元件: {name}, R={R_value}")
        print(f"    本构关系: {component['constitutive']}")
        
        return component
    
    def add_I_element(self, name, I_value, initial_state=0.0):
        """
        添加惯性元件 (I-element)
        
        本构关系：
        - p = I * f   (p: 广义动量)
        - dp/dt = e
        
        物理类比：
        - 电感：p=磁链, e=电压, f=电流
        - 质量：p=动量, e=力, f=速度
        """
        component = {
            'name': name,
            'type': 'I',
            'value': I_value,
            'state': sp.Symbol(f'p_{name}'),
            'effort': sp.Symbol(f'e_{name}'),
            'flow': sp.Symbol(f'f_{name}'),
            'initial_state': initial_state
        }
        
        # f = p / I
        component['constitutive'] = sp.Eq(
            component['flow'],
            component['state'] / I_value
        )
        
        self.components.append(component)
        self.state_vars.append(component['state'])
        
        print(f"  ✓ 添加I元件: {name}, I={I_value}, p0={initial_state}")
        print(f"    本构关系: {component['constitutive']}")
        
        return component
    
    def add_0_junction(self, name):
        """
        添加0-结点 (0-junction)
        
        物理约束：
        - 所有连接的effort相等
        - 所有连接的flow之和为0
        
        类比：电路中的节点（KCL）
        """
        junction = {
            'name': name,
            'type': '0',
            'connected_components': []
        }
        self.components.append(junction)
        
        print(f"  ✓ 添加0-结点: {name}")
        print(f"    约束: 所有e相等, Σf=0")
        
        return junction
    
    def connect(self, comp1, comp2):
        """
        连接两个组件（创建Bond）
        
        Bond传递能量：Power = e * f
        """
        bond = {
            'from': comp1['name'],
            'to': comp2['name'],
            'effort': sp.Symbol(f"e_{comp1['name']}_{comp2['name']}"),
            'flow': sp.Symbol(f"f_{comp1['name']}_{comp2['name']}")
        }
        
        self.bonds.append(bond)
        
        print(f"  ✓ 连接: {comp1['name']} ←→ {comp2['name']}")
        
        return bond
    
    def analyze(self):
        """
        ===== 核心：模型解析 =====
        
        步骤：
        1. 因果关系分析 (Causality Assignment)
        2. 生成状态空间方程
        3. 符号化简
        """
        print(f"\n{'='*70}")
        print("🔍 模型解析过程")
        print(f"{'='*70}")
        
        # 第一步：识别状态变量
        print("\n1️⃣  识别状态变量:")
        for var in self.state_vars:
            print(f"   - {var}")
        
        # 第二步：应用本构关系
        print("\n2️⃣  本构关系:")
        for comp in self.components:
            if 'constitutive' in comp:
                print(f"   {comp['name']}: {comp['constitutive']}")
        
        # 第三步：应用结点约束（简化示例）
        print("\n3️⃣  结点约束:")
        for comp in self.components:
            if comp['type'] == '0':
                print(f"   {comp['name']}: 所有e相等, Σf=0")
        
        # 第四步：生成状态方程
        print("\n4️⃣  生成状态方程:")
        self._generate_state_equations()
        
    def _generate_state_equations(self):
        """
        生成状态空间方程: dx/dt = f(x, u, t)
        
        对于简单RC系统：
        - C元件：dq/dt = f
        - 通过结点约束和本构关系，表达f与其他变量的关系
        """
        # 这里手动构建简单RC电路的方程
        # 实际BGT会自动完成这个过程
        
        if len(self.state_vars) > 0:
            # 假设第一个状态是C元件的q
            q = self.state_vars[0]
            
            # 对于RC电路: dq/dt = -q/(R*C)
            # 这里简化表示
            print(f"   d{q}/dt = f({q})")
    
    def generate_numerical_model(self):
        """
        ===== 核心：转换为数值模型 =====
        
        将符号方程转换为可以数值求解的函数
        """
        print(f"\n{'='*70}")
        print("🔧 转换为数值模型")
        print(f"{'='*70}")
        
        # 收集参数
        params = {}
        for comp in self.components:
            if 'value' in comp:
                params[comp['name']] = comp['value']
        
        print(f"\n参数: {params}")
        
        # 为简单RC系统手动创建ODE函数
        # 实际BGT会自动生成
        
        C_comp = next((c for c in self.components if c['type'] == 'C'), None)
        R_comp = next((c for c in self.components if c['type'] == 'R'), None)
        
        if C_comp and R_comp:
            C_val = C_comp['value']
            R_val = R_comp['value']
            
            def ode_func(x, t):
                """状态方程: dx/dt = f(x, t)"""
                q = x[0]  # C元件的电荷
                
                # 通过本构关系
                e_C = q / C_val  # 电容电压
                
                # 通过结点约束: e_R = e_C
                e_R = e_C
                
                # 电阻定律: f_R = e_R / R
                f_R = e_R / R_val
                
                # 电容: dq/dt = -f_R (流出为负)
                dq_dt = -f_R
                
                return [dq_dt]
            
            print(f"\n生成的ODE函数:")
            print(f"  dq/dt = -q / (R*C)")
            print(f"  dq/dt = -q / ({R_val}*{C_val})")
            print(f"  dq/dt = -q / {R_val * C_val}")
            
            return ode_func
        
        return None
    
    def solve(self, timespan, initial_conditions):
        """
        ===== 核心：数值求解 =====
        
        使用ODE求解器求解状态方程
        """
        print(f"\n{'='*70}")
        print("🚀 数值求解")
        print(f"{'='*70}")
        
        # 生成数值模型
        ode_func = self.generate_numerical_model()
        
        if ode_func is None:
            print("❌ 无法生成数值模型")
            return None, None
        
        # 时间向量
        t = np.linspace(timespan[0], timespan[1], 200)
        
        # 初始条件
        x0 = initial_conditions
        
        print(f"\n求解参数:")
        print(f"  时间范围: {timespan}")
        print(f"  初始条件: {x0}")
        print(f"  时间步数: {len(t)}")
        
        # 求解ODE
        print(f"\n⏳ 求解中...")
        x = odeint(ode_func, x0, t)
        
        print(f"✅ 求解完成!")
        print(f"   解的形状: {x.shape}")
        
        return t, x


# ==================== 使用示例：RC电路 ====================

def example_RC_circuit():
    """
    示例：RC电路放电
    
    电路：
    C (电容, 初始电荷q0) ←→ 0-junction ←→ R (电阻)
    
    物理过程：
    - 电容通过电阻放电
    - 电压和电流随时间指数衰减
    
    理论解：
    - q(t) = q0 * exp(-t / (R*C))
    - V(t) = V0 * exp(-t / (R*C))
    """
    print("\n" + "="*70)
    print("示例：RC电路放电")
    print("="*70)
    
    # 创建模型
    model = SimpleBondGraph("RC_Circuit")
    
    # 添加组件
    print("\n📦 步骤1: 创建组件")
    C = model.add_C_element(name='C1', C_value=1.0, initial_state=10.0)  # 1F, q0=10C
    R = model.add_R_element(name='R1', R_value=2.0)  # 2Ω
    junction = model.add_0_junction(name='j0')
    
    # 连接组件
    print("\n🔗 步骤2: 连接组件")
    model.connect(C, junction)
    model.connect(junction, R)
    
    # 解析模型
    print("\n🔍 步骤3: 解析模型")
    model.analyze()
    
    # 求解
    print("\n🚀 步骤4: 数值求解")
    timespan = [0, 10]  # 0到10秒
    initial_conditions = [10.0]  # q0 = 10C
    
    t, x = model.solve(timespan, initial_conditions)
    
    # 后处理
    print("\n📊 步骤5: 结果分析")
    if x is not None:
        q = x[:, 0]  # 电荷
        V = q / C['value']  # 电压 V = q/C
        
        # 理论解
        R_val = R['value']
        C_val = C['value']
        tau = R_val * C_val
        q_theory = initial_conditions[0] * np.exp(-t / tau)
        
        print(f"\n时间常数 τ = R*C = {tau} 秒")
        print(f"最终电荷: q({t[-1]:.1f}s) = {q[-1]:.4f} C")
        print(f"理论值: q({t[-1]:.1f}s) = {q_theory[-1]:.4f} C")
        print(f"误差: {abs(q[-1] - q_theory[-1]):.2e}")
        
        # 可视化
        plt.figure(figsize=(12, 5))
        
        plt.subplot(1, 2, 1)
        plt.plot(t, q, 'b-', linewidth=2, label='数值解')
        plt.plot(t, q_theory, 'r--', linewidth=2, label='理论解')
        plt.xlabel('时间 (s)')
        plt.ylabel('电荷 q (C)')
        plt.title('RC电路放电 - 电荷')
        plt.legend()
        plt.grid(True)
        
        plt.subplot(1, 2, 2)
        plt.plot(t, V, 'b-', linewidth=2, label='电压')
        plt.xlabel('时间 (s)')
        plt.ylabel('电压 V (V)')
        plt.title('RC电路放电 - 电压')
        plt.legend()
        plt.grid(True)
        
        plt.tight_layout()
        # plt.savefig('/tmp/rc_circuit.png', dpi=100, bbox_inches='tight')
        print(f"\n📈 图表已保存")
        plt.show()


# ==================== 使用真实BGT（如果可用） ====================

def example_with_real_BGT():
    """使用真实的BondGraphTools"""
    if not HAS_BGT:
        print("\n⚠️  BondGraphTools未安装，跳过此示例")
        return
    
    print("\n" + "="*70)
    print("使用真实BondGraphTools: RC电路")
    print("="*70)
    
    # 创建模型
    model = bgt.new(name="RC_BGT")
    
    print("\n📦 创建组件:")
    # 创建组件
    C = bgt.new("C", value=1.0, name="C1")
    R = bgt.new("R", value=2.0, name="R1")
    junction = bgt.new("0", name="j0")
    
    # 添加到模型
    bgt.add(model, C, R, junction)
    print(f"  ✓ 添加了 C, R, 0-junction")
    
    print("\n🔗 建立连接:")
    # 连接
    bgt.connect(C, junction)
    bgt.connect(R, junction)
    print(f"  ✓ C ←→ 0-junction ←→ R")
    
    print("\n🔍 模型信息:")
    print(f"  组件数: {len(model.components)}")
    print(f"  状态变量: {model.state_vars}")
    
    print("\n🚀 运行仿真:")
    # 仿真
    timespan = [0, 10]
    x0 = [10.0]  # 初始电荷
    
    try:
        t, x = bgt.simulate(model, timespan=timespan, x0=x0)
        
        print(f"✅ 仿真完成")
        print(f"  时间点数: {len(t)}")
        print(f"  最终状态: {x[-1]}")
        
        # 绘图
        plt.figure(figsize=(8, 5))
        plt.plot(t, x[:, 0], 'b-', linewidth=2)
        plt.xlabel('时间 (s)')
        plt.ylabel('电荷 q (C)')
        plt.title('BondGraphTools: RC电路放电')
        plt.grid(True)
        plt.show()
        
    except Exception as e:
        print(f"❌ 仿真失败: {e}")


# ==================== 主程序 ====================

if __name__ == "__main__":
    # 运行简化版本（展示原理）
    example_RC_circuit()
    
    # 如果有BGT，运行真实版本
    example_with_real_BGT()
    
    print("\n" + "="*70)
    print("✅ 演示完成")
    print("="*70)
    print("\n关键流程回顾:")
    print("1️⃣  模块创建: 定义组件类型和参数")
    print("2️⃣  模块连接: 通过Bond连接组件")
    print("3️⃣  模型解析: 因果分析 + 方程生成")
    print("4️⃣  方程转换: 符号方程 → 数值函数")
    print("5️⃣  数值求解: ODE求解器求解状态方程")