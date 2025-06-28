

#### 下面描述生成通用的函数然后用来求解方程
import sympy as sp
import numpy as np
from scipy.interpolate import CubicSpline

def solve_dae_with_knowns(equations, known_vars, parameters):
    """
    使用已知变量和参数求解微分代数方程组
    :param equations: 方程列表
    :param known_vars: 已知变量字典 {var_name: values}
    :param parameters: 参数字典 {param_name: value}
    :return: 求解后的变量字典 {var_name: values}
    """
    solutions=sp.solve(equations, [e_0, f_0])
    dt=1
    for unknowvar in solutions:
        for func in solutions[unknowvar].atoms(sp.Derivative):
            if func.expr in known_vars:
                # 替换微分符号为对应的函数
                import  numpy as np
                value=np.gradient(known_vars[func.expr],dt)
                known_vars[func] = list(value)
        # 代入已知值和参数
        for i in range(3):
            # 代入已知值和参数
            substituted_dict=  {
                    k: v[i] if isinstance(v, list) else v
                    for k, v in known_vars.items()
                }
            value = solutions[unknowvar].subs(substituted_dict).subs(parameters)
            
            # 如果是简单表达式则计算数值
            if value.is_constant():
                value = value.evalf()
            varlist = value
            
        known_vars[unknowvar] = value


# === Step 1: 符号定义（与你之前的代码一致）===
t = sp.symbols('t')
C= sp.symbols('C')
q_0 = sp.Function('q_0')(t)
e_0 = sp.Function('e_0')(t)
f_0 = sp.Function('f_0')(t)
dq_0 = sp.diff(q_0,t)

# 创建方程
eq1 = sp.Eq(dq_0+q_0 , f_0)
eq2 = sp.Eq(q_0 , C * e_0)
# eq3 = sp.Eq(x + 2*y - z, a)

# 已知变量和参数
known_vars = {q_0: [2,2,3],t:[0,1,2]}
parameters = {C:1}
equations = [eq1, eq2]
solutions=sp.solve([eq1, eq2], [e_0, f_0])
known_vars[dq_0] = list(np.gradient(known_vars[q_0], 1))  # 假设 dt=1
dt=1
sub_solutions = [ v.subs(parameters) for k, v in solutions.items()]
f=sp.lambdify(tuple(known_vars.keys()), sub_solutions, modules='numpy')
arrayvalue=tuple(np.array(known_vars[var]) for var in known_vars.keys())
value=f(*arrayvalue)
