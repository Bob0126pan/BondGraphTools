
import numpy as np
import sympy as sp
from scipy.interpolate import CubicSpline
from scipy.optimize import fsolve
import re
from scipy.integrate import solve_ivp
from scipy import interpolate
import warnings

def solve_state_equations(equations, known_data, params=None):
    """
    改进的状态方程求解器，充分利用已知时间序列数据
    
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
    known_var_indices = {}  # 存储已知变量的索引
    
    # 添加状态变量 - 只添加未知的状态变量
    unknown_state_vars = []
    for var in state_vars:
        var_str = str(var)
        if var_str not in known_data or var_str == 't':
            # 只有未知变量才作为状态变量
            unknown_state_vars.append(var)
    
    # 为未知状态变量创建映射
    for i, var in enumerate(unknown_state_vars):
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
    n_states = len(unknown_state_vars)
    
    # 为已知变量创建插值函数
    interp_funcs = {}
    for var in known_vars:
        if var != 't' and var in known_data:
            # 创建线性插值函数
            interp_funcs[var] = interpolate.interp1d(
                t_data, known_data[var], 
                kind='linear', fill_value="extrapolate"
            )
    
    # 创建残差函数
    def residual_func(t, y, yp):
        # 准备变量值字典
        var_values = {}
        
        # 添加时间
        var_values[t_sym] = t
        
        # 添加已知变量的值（通过插值）
        for var, func in interp_funcs.items():
            var_values[sp.Symbol(var)] = func(t)
        
        # 添加未知状态变量及其导数
        for i in range(n_states):
            var_values[X[i]] = y[i]
            var_values[dX[i]] = yp[i]
        
        # 计算残差
        res = np.zeros(n_states)
        for i, eq in enumerate(residuals):
            try:
                # 使用数值计算残差
                res[i] = float(eq.evalf(subs=var_values))
            except Exception as e:
                warnings.warn(f"Error evaluating equation {i} at t={t}: {str(e)}")
                res[i] = 0.0
        
        return res
    
    # 准备初始条件
    y0 = np.zeros(n_states)
    yp0 = np.zeros(n_states)
    
    # 填充未知状态变量的初始值
    for i, var in enumerate(unknown_state_vars):
        var_str = str(var)
        # 如果没有提供初始值，使用0
        if var_str in known_data:
            y0[i] = known_data[var_str][0]
        else:
            y0[i] = 0.0
    
    # 设置DAE求解器选项
    # 使用BDF方法，适合刚性问题
    if n_states > 0:
        # 只有存在未知状态变量时才需要求解
        sol = solve_ivp(
            lambda t, y: residual_func(t, y, np.zeros_like(y)),  # Dummy derivative
            [t_data[0], t_data[-1]],
            y0,
            method='BDF',
            t_eval=t_data,
            vectorized=False
        )
        
        # 检查求解是否成功
        if not sol.success:
            warnings.warn(f"DAE求解警告: {sol.message}")
    else:
        # 没有未知状态变量，创建空解
        sol = type('', (), {})()
        sol.y = np.zeros((0, len(t_data)))
        sol.t = t_data
    
    # 提取结果
    result = {'t': t_data.tolist()}
    
    # 添加未知状态变量结果
    for i, var in enumerate(unknown_state_vars):
        result[str(var)] = sol.y[i].tolist()
    
    # 添加已知数据（覆盖插值值）
    for k, v in known_data.items():
        if k != 't':
            result[k] = v
    
    # 计算并添加导数变量
    for deriv_sym, base_var in derivative_vars.items():
        base_sym = sp.Symbol(base_var)
        base_str = str(base_sym)
        
        if base_str in result:
            # 使用样条插值计算导数
            cs = CubicSpline(t_data, result[base_str])
            deriv_vals = cs(t_data, 1)
            result[str(deriv_sym)] = deriv_vals.tolist()
        elif base_str in known_data:
            # 使用已知数据计算导数
            cs = CubicSpline(t_data, known_data[base_str])
            deriv_vals = cs(t_data, 1)
            result[str(deriv_sym)] = deriv_vals.tolist()
    
    # 对于非状态变量和导数变量的方程，使用代数求解
    algebraic_vars = state_vars - set(unknown_state_vars) - set(derivative_vars.keys())
    
    if algebraic_vars:
        # 创建代数求解的符号表达式
        algebraic_eqs = [eq.subs(param_subs) for eq in sym_equations]
        
        # 为每个时间点求解代数方程
        for var in algebraic_vars:
            result[str(var)] = np.zeros(n_points)
        
        for i in range(n_points):
            # 当前时间点的所有已知值
            current_vals = {
                't': t_data[i]
            }
            for k in result:
                if k != 't' and k in result:
                    current_vals[k] = result[k][i]
            
            # 尝试符号求解
            try:
                sol_dict = sp.solve(algebraic_eqs, list(algebraic_vars), dict=True)
                if sol_dict:
                    for var in algebraic_vars:
                        expr = sol_dict[0].get(var, 0)
                        if expr:
                            # 代入当前值
                            value = expr.subs(current_vals)
                            try:
                                result[str(var)][i] = float(value)
                            except TypeError:
                                result[str(var)][i] = float(sp.re(value))
            except Exception:
                # 符号求解失败，使用数值求解
                def algebraic_func(x):
                    subs_dict = dict(zip(algebraic_vars, x))
                    subs_dict.update(current_vals)
                    return [float(eq.subs(subs_dict)) for eq in algebraic_eqs]
                
                x0 = [result[str(var)][i-1] if i > 0 else 0 for var in algebraic_vars]
                sol = fsolve(algebraic_func, x0)
                
                for j, var in enumerate(algebraic_vars):
                    result[str(var)][i] = sol[j]
    
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