

#### 下面描述生成通用的函数然后用来求解方程
import sympy as sp
import numpy as np
from scipy.interpolate import CubicSpline

def solve_dae_with_knowns(equations, known_data, params=None):
    # 提取时间数据
    t_data = known_data['t']
    n_points = len(t_data)
    
    # 构建所有符号
    all_symbols = set()
    for eq in equations:
        all_symbols.update(eq.free_symbols)
    
    # 分离常数参数和变量
    params = params or {}
    param_syms = {sp.symbols(k) for k in params.keys()}
    variable_syms = list(all_symbols - param_syms)

    # 区分状态变量（可随时间变化）和其它变量
    time_varying_vars = []
    other_vars = []

    for sym in variable_syms:
        if any(str(sym).endswith(f"_{i}") for i in range(10)) or str(sym).startswith("d"):
            other_vars.append(sym)
        else:
            time_varying_vars.append(sym)

    # 构造插值函数
    interpolators = {}
    for var_name, values in known_data.items():
        if var_name == 't':
            continue
        var_sym = sp.symbols(var_name)
        if var_sym in all_symbols:
            interpolators[var_sym] = CubicSpline(t_data, values)

    # 替换参数
    subs_dict = {sp.symbols(k): v for k, v in params.items()}

    # 替换已知变量为插值函数
    numeric_equations = []
    for eq in equations:
        eq_subbed = eq.subs({
            var: interpolators[var](t_data) if var in interpolators else var
            for var in eq.free_symbols
        })
        numeric_equations.append(eq_subbed)

    # 解代数方程
    algebraic_solutions = {}
    for eq in numeric_equations:
        lhs, rhs = eq.lhs, eq.rhs
        if lhs.is_Function and rhs != 0:
            var = lhs.func(*lhs.args[:-1]) if lhs.is_Function else lhs
            algebraic_solutions[lhs] = rhs
        elif rhs.is_Function and lhs != 0:
            var = rhs.func(*rhs.args[:-1]) if rhs.is_Function else rhs
            algebraic_solutions[rhs] = lhs

    # 计算所有变量的值
    results = {'t': t_data.tolist()}
    
    # 添加已知变量
    for name, data in known_data.items():
        if name != 't':
            results[name] = data.tolist()

    # 求解代数变量
    for var, expr in algebraic_solutions.items():
        try:
            # 替换所有已知变量
            val = expr.evalf(subs=subs_dict)
            if isinstance(val, sp.Matrix):
                val = val.tolist()
            results[str(var)] = val
        except Exception as e:
            print(f"无法计算 {var}: {str(e)}")

    return results

# === Step 1: 符号定义（与你之前的代码一致）===
t_sym = sp.symbols('t')
q_0_sym = sp.Function('q_0')(t_sym)
dq_0_sym = sp.diff(q_0_sym, t_sym)
f_0_sym = sp.Function('f_0')(t_sym)
e_0_sym = sp.Function('e_0')(t_sym)
C = sp.symbols('C')

equations = [
   sp.Eq(dq_0_sym, f_0_sym),        # dq_0 = f_0
   sp.Eq(q_0_sym - C * e_0_sym, 0)    # q_0 = C * e_0
]

known_data = {
    't': np.array([0, 0.1, 0.2, 0.3, 0.4]),
    'q_0': np.array([0, 1, 2, 3, 4])
}

params = {
    'C': 2.0
}

results = solve_dae_with_knowns(equations, known_data, params)

# 打印结果
for var, values in results.items():
    print(f"{var}: {values}")