

import sympy as sp
import numpy as np
from scipy.interpolate import interp1d

def solve_state_equations_general(eq_strs, known_data, params, t_data):
    """
    通用状态方程求解器（已修复导数依赖）

    Parameters:
    - eq_strs: 状态方程（字符串列表），如 ['q_0 - C * e_0', 'dq_0 - f_0']
    - known_data: dict，已知变量的时间序列数据，变量名为键（如 'q_0'），值为 numpy 数组
    - params: dict，常数参数替换，如 {'C': 2.0}
    - t_data: 1D numpy 数组，时间序列

    Returns:
    - List of dicts：每个时刻对应的变量解（含已知 + 求解出）
    """
    t = sp.Symbol('t')

    # 替换 dq_0 为 Derivative(q_0(t), t)
    def convert_eq(eq_str):
        expr = eq_str
        for var in known_data:
            if var.startswith('d'):
                base = var[1:]
                expr = expr.replace(var, f"Derivative({base}(t), t)")
            else:
                expr = expr.replace(var, f"{var}(t)")
        for p in params:
            expr = expr.replace(p, p)
        return sp.sympify(expr)

    eqs = [convert_eq(eq_str) for eq_str in eq_strs]
    all_symbols = set().union(*[eq.free_symbols for eq in eqs])

    # 构造插值函数
    interp_funcs = {
        var: interp1d(t_data, data, kind='cubic', fill_value="extrapolate")
        for var, data in known_data.items()
        if not var.startswith('d')
    }

    # 数值导数函数（中心差分）
    def numeric_derivative(f, x, dx=1e-6):
        return (f(x + dx) - f(x - dx)) / (2 * dx)

    # 获取变量值
    def get_value(var, t_val):
        if var.startswith('d'):
            base = var[1:]
            f = interp_funcs[base]
            return float(numeric_derivative(f, t_val))
        else:
            return float(interp_funcs[var](t_val))

    results = []

    for t_val in t_data:
        substitutions = {sp.Symbol(p): v for p, v in params.items()}

        # 添加已知变量与导数
        for var in known_data:
            val = get_value(var, t_val)
            if var.startswith('d'):
                base = var[1:]
                substitutions[sp.Derivative(sp.Function(base)(t), t)] = val
            else:
                substitutions[sp.Function(var)(t)] = val

        # 解方程
        unknowns = list(all_symbols - set(substitutions.keys()) - {t})
        sol = sp.solve(eqs, unknowns, dict=True)

        if not sol:
            raise ValueError(f"No solution found at t={t_val}")

        row = {'t': t_val}
        for var in known_data:
            row[var] = get_value(var, t_val)
        for s in sol[0]:
            row[str(s)] = float(sol[0][s])

        results.append(row)

    return results

if __name__ == "__main__":
    # 示例数据
    eqs = ['q_0 - C * e_0', 'dq_0 - f_0']
    t_vals = np.linspace(0, 10, 100)
    q0_vals = np.sin(t_vals)

    data = {
        'q_0': q0_vals
    }

    params = {
        'C': 2.0
    }

    # 调用函数
    res = solve_state_equations_general(eqs, data, params, t_vals)

    # 查看前几项
    for r in res[:3]:
        print(r)