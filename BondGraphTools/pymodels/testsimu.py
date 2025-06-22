
import numpy as np
import sympy as sp
from scipy.interpolate import CubicSpline
from scipy.optimize import fsolve
import re
from scipy.integrate import solve_ivp
from scipy import interpolate

def solve_state_equations(equations, known_data, params=None):
    """
    基于DAE求解器的状态方程求解器
    
    参数:
    equations: 状态方程列表 (字符串表达式), 形式为['表达式1 = 0', '表达式2 = 0']
    known_data: 包含时间序列和已知变量数据的字典 
                {'t': [t0, t1, ...], 'var1': [y0, y1, ...], ...}
    params: 参数字典 {'param1': value1, 'param2': value2, ...}
    
    返回:
    包含所有变量时间序列数据的字典
    """
    # 默认参数处理
    params = params or {}
    
    # 符号化方程
    sym_equations = [sp.sympify(eq.replace('=', '-')) for eq in equations]
    
    # 提取所有符号
    all_symbols = set()
    for eq in sym_equations:
        all_symbols |= eq.free_symbols
    
    # 分离参数和变量
    param_symbols = set(sp.symbols(list(params.keys()))) if params else set()
    variable_symbols = all_symbols - param_symbols
    
    # 识别时间变量
    t_sym = sp.Symbol('t')
    if t_sym not in variable_symbols:
        t_sym = None
        for sym in variable_symbols:
            if str(sym) == 't':
                t_sym = sym
                break
    
    # 识别状态变量和导数变量
    state_vars = set()
    derivative_vars = {}
    for sym in variable_symbols:
        sym_str = str(sym)
        # 识别导数变量 (形式如 dq_0)
        if sym_str.startswith('d') and sym_str[1:] in [str(s) for s in variable_symbols]:
            base_var = sym_str[1:]
            derivative_vars[sym] = base_var
        else:
            state_vars.add(sym)
    
    # 确定已知和未知变量
    known_vars = set(known_data.keys()) - {'t'}
    if t_sym is not None:
        known_vars.add(str(t_sym))
    
    # 检查已知数据完整性
    for var in known_vars:
        if var not in known_data and var != 't':
            raise ValueError(f"缺少已知变量 '{var}' 的数据")
    
    # 准备时间序列
    t_data = np.array(known_data['t'])
    n_points = len(t_data)
    
    # 创建符号变量
    X = sp.IndexedBase('X')
    dX = sp.IndexedBase('dX')
    
    # 创建变量映射
    var_map = {}
    inv_var_map = {}
    state_var_list = []
    
    # 添加状态变量
    for i, var in enumerate(state_vars):
        var_map[var] = X[i]
        inv_var_map[X[i]] = var
        state_var_list.append(var)
    
    # 添加导数变量
    for deriv_sym, base_var in derivative_vars.items():
        base_sym = sp.Symbol(base_var)
        if base_sym in var_map:
            idx = list(var_map.keys()).index(sp.symbols(base_var))
            var_map[deriv_sym] = dX[idx]
    
    # 参数替换
    param_subs = {sp.Symbol(k): v for k, v in params.items()}
    
    # 构建残差函数
    residuals = []
    for eq in sym_equations:
        # 应用变量映射
        eq_mapped = eq.subs(var_map)
        # 应用参数替换
        eq_substituted = eq_mapped.subs(param_subs)
        residuals.append(eq_substituted)
    
    # 创建残差函数的lambda表达式
    t_sym = sp.Symbol('t')
    n_states = len(state_vars)
    
    # 创建残差函数
    F = sp.lambdify((t_sym, [X[i] for i in range(n_states)], [dX[i] for i in range(n_states)]), residuals, modules='numpy')
    
    # 定义DAE残差函数
    def residual_func(t, y, yp):
        # y是状态变量值, yp是状态变量的导数
        # 我们需要同时满足残差方程
        res = np.zeros_like(y)
        result = F(t, y, yp)
        for i in range(len(result)):
            res[i] = result[i]
        return res
    
    # 准备初始条件
    # 使用第一个时间点的值作为初始猜测
    y0 = np.zeros(n_states)
    yp0 = np.zeros(n_states)
    
    # 填充已知变量
    for i, var in enumerate(state_var_list):
        if str(var) in known_data:
            y0[i] = known_data[str(var)][0]
    
    # 对于导数变量，使用数值微分估计初始导数
    for deriv_sym, base_var in derivative_vars.items():
        base_sym = sp.Symbol(base_var)
        if base_sym in var_map:
            idx = list(var_map.keys()).index(sp.symbols(base_var))
            # 使用前向差分估计初始导数
            if n_points > 1:
                dt = t_data[1] - t_data[0]
                yp0[idx] = (known_data[base_var][1] - known_data[base_var][0]) / dt
    
    # 设置DAE求解器选项
    # 使用BDF方法，适合刚性问题
    sol = solve_ivp(
        lambda t, y: residual_func(t, y, np.zeros_like(y)),  # 对于DAE，我们需要特殊的处理
        [t_data[0], t_data[-1]],
        y0,
        method='BDF',
        t_eval=t_data,
        vectorized=True
    )
    
    # 检查求解是否成功
    if not sol.success:
        raise RuntimeError(f"DAE求解失败: {sol.message}")
    
    # 提取结果
    result = {'t': t_data.tolist()}
    
    # 添加状态变量结果
    for i, var in enumerate(state_var_list):
        result[str(var)] = sol.y[i].tolist()
    
    # 添加已知数据（可能包含不在状态变量中的量）
    for k, v in known_data.items():
        if k != 't' and k not in result:
            result[k] = v
    
    # 计算并添加导数变量
    for deriv_sym, base_var in derivative_vars.items():
        base_sym = sp.Symbol(base_var)
        if base_sym in var_map:
            idx = list(var_map.keys()).index(sp.symbols(base_var))
            # 使用样条插值计算导数
            cs = CubicSpline(t_data, sol.y[idx])
            deriv_vals = cs(t_data, 1)
            result[str(deriv_sym)] = deriv_vals.tolist()
    
    return result

if __name__ == "__main__":
    # 测试数据
    # 示例1：原始问题
    equations = ['q_0 - C * e_0', 'dq_0 - f_0']
    known_data = {
        't': [0, 1, 2, 3, 4],
        'q_0': [1, 0.5, 1.0, 0.5, 0]
    }
    params = {'C': 2.5}
    result = solve_state_equations(equations, known_data, params)

    # 示例2：更复杂的情况
    equations = [
        'm * D2(x) + c * D(x) + k * x - F',
        'I * D2(theta) + b * D(theta) + k_t * theta - T'
    ]
    known_data = {
        't': np.linspace(0, 10, 100).tolist(),
        'x': np.sin(np.linspace(0, 10, 100)).tolist(),
        'theta': np.cos(np.linspace(0, 10, 100)).tolist()
    }
    params = {'m': 1.0, 'c': 0.1, 'k': 2.0, 'I': 0.5, 'b': 0.05, 'k_t': 1.5}
    result = solve_state_equations(equations, known_data, params)

    # 示例3：带时间显式项
    equations = ['x*sin(t) - y*cos(t) - a', 'D(x)*t - D(y) - b']
    known_data = {
        't': np.linspace(0, 5, 50).tolist(),
        'x': np.exp(-0.2 * np.linspace(0, 5, 50)).tolist()
    }
    params = {'a': 0.5, 'b': 1.0}
    result = solve_general_state_equations(equations, known_data, params)