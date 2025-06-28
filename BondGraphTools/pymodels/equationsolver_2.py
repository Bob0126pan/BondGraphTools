import sympy as sp
from collections import deque, defaultdict
from typing import Dict, List, Set, Tuple, Union, Optional

## 这个方程是可能比较复杂的方程求解器 本项目暂时不用，其他的用

class EquationSolver:
    def __init__(self, equations: List[sp.Eq], known_vars: Dict[sp.Symbol, float], 
                 parameters: Dict[sp.Symbol, float], independent_var: sp.Symbol = sp.symbols('t')):
        """
        初始化方程求解器
        
        :param equations: 原始方程列表
        :param known_vars: 已知变量字典 {符号: 值}
        :param parameters: 参数字典 {符号: 值}
        :param independent_var: 自变量符号 (默认为时间 t)
        """
        self.raw_equations = equations
        self.known_vars = known_vars
        self.parameters = parameters
        self.independent_var = independent_var
        
        # 预处理后的方程组
        self.equations = []
        # 方程状态数据库
        self.equation_db = []
        # 求解结果
        self.solutions = {}
        # 需要数值求解的方程组
        self.numeric_groups = []
        # 符号映射
        self.symbol_map = {}
        
        # 初始化求解状态
        self._initialize_solver()


    def _initialize_solver(self):
        """初始化求解器状态"""
        self.solutions = self.known_vars.copy()
        self.solutionsNum=len(next(iter(self.known_vars.values()))) if  isinstance(next(iter(self.known_vars.values())), (list, tuple)) else 1
        self.singleresult = False if self.solutionsNum > 1 else True
        # self._preprocess_equations()
        self._classify_symbols()
        self._build_equation_database()
        self._create_symbol_mappings()

    def _preprocess_equations(self):
        """预处理方程：添加微分关系定义"""
        self.equations = self.raw_equations.copy()
        
        # 收集所有符号
        all_symbols = set()
        for eq in self.equations:
            all_symbols |= eq.free_symbols
        
        # 识别微分符号并添加定义方程
        for symbol in all_symbols:
            if symbol.name.startswith('d') and symbol.name[1:].isidentifier():
                base_name = symbol.name[1:]
                base_var = sp.Function(base_name)(self.independent_var)
                
                # 创建微分定义方程
                diff_eq = sp.Eq(symbol, sp.Derivative(base_var, self.independent_var))
                self.equations.append(diff_eq)
                
                # 添加到符号映射
                self.symbol_map[symbol] = base_var

    def _classify_symbols(self):
        """分类符号：参数、代数变量、函数变量、微分符号"""
        self.equations = self.raw_equations.copy()
        # 收集所有符号
        self.all_symbols = set()
        for eq in self.equations:
            self.all_symbols |= eq.free_symbols | eq.atoms(sp.Function)
        
        # 排除参数和已知变量
        self.parameters_set = set(self.parameters.keys())
        self.known_set = set(self.known_vars.keys())
        
        # 分类剩余符号
        self.algebraic_vars = set()
        self.function_vars = set()
        self.diff_symbols = set()
        
        for symbol in self.all_symbols- self.parameters_set - self.known_set:
            if isinstance(symbol, sp.Function):
                self.function_vars.add(symbol)
            elif symbol in self.symbol_map.values():
                self.function_vars.add(symbol)
            elif symbol in self.symbol_map:
                self.diff_symbols.add(symbol)
            else:
                self.algebraic_vars.add(symbol)
        
        # 完整未知量集合
        self.unknowns = self.algebraic_vars | self.function_vars

    def _build_equation_database(self):
        """构建方程数据库"""
        self.equation_db = []
        for i, eq in enumerate(self.equations):
            # 提取方程中的所有符号
            all_symbols = eq.free_symbols | eq.atoms(sp.Function)
            # 分类方程类型
            if any(isinstance(term, sp.Derivative) for term in eq.atoms(sp.Derivative)):
                eq_type = 'DIFFERENTIAL'
            elif any(sym in self.diff_symbols for sym in all_symbols):
                eq_type = 'DIFFERENTIAL'
            else:
                eq_type = 'ALGEBRAIC'
            
            # 计算未知量
            unknowns = [sym for sym in all_symbols if sym in self.unknowns]
            
            self.equation_db.append({
                'id': i,
                'equation': eq,
                'type': eq_type,
                'symbols': all_symbols,
                'unknowns': set(unknowns),
                'solved': False
            })

    def _create_symbol_mappings(self):
        """创建符号映射关系"""
        # 方程到未知量映射
        self.eq_to_unknowns = {}
        # 未知量到方程映射
        self.unknown_to_eqs = defaultdict(list)
        
        for eq_data in self.equation_db:
            eq_id = eq_data['id']
            unknowns = eq_data['unknowns']
            self.eq_to_unknowns[eq_id] = unknowns
            
            for unknown in unknowns:
                self.unknown_to_eqs[unknown].append(eq_id)

    def _update_equation_unknowns(self, solved_var: sp.Symbol):
        """更新方程未知量状态"""
        # 从所有相关方程中移除已解变量
        for eq_id in self.unknown_to_eqs.get(solved_var, []):
            if solved_var in self.eq_to_unknowns[eq_id]:
                self.eq_to_unknowns[eq_id].remove(solved_var)
                
                # 更新方程数据库
                for eq_data in self.equation_db:
                    if eq_data['id'] == eq_id:
                        eq_data['unknowns'].discard(solved_var)
                        break

    def _solve_algebraic_equation(self, eq_data: dict) -> bool:
        """尝试求解代数方程"""
        eq = eq_data['equation']
        unknowns = list(eq_data['unknowns'])
        
        # 无未知数 - 验证方程
        if len(unknowns) == 0:
            # 代入已知值和参数
            substituted = eq.subs(self.solutions).subs(self.parameters)
            if not substituted:
                raise ValueError(f"Equation {eq} is not satisfied with given values")
            eq_data['solved'] = True
            return True
        
        # 单变量方程 - 直接求解
        if len(unknowns) == 1:
            var = unknowns[0]
            try:
                # 尝试求解方程
                solution = sp.solve(eq, var, dict=True)
                if not solution:
                    return False
                
                # 处理多解情况（取第一个有效解）
                for sol in solution:
                    if var in sol:
                        # 代入已知值和参数
                        self.solve_equation_var(var, sol[var])
                        eq_data['solved'] = True
                        return True
                return False
            except NotImplementedError:
                return False
        
        return False

    def _solve_differential_equation(self, eq_data: dict) -> bool:
        """尝试求解微分方程"""
        eq = eq_data['equation']
        fun_var=next(iter((eq_data['unknowns'])))
        
        try:
            # 尝试符号求解
            solution = sp.solve(eq, fun_var, dict=True)
            if not solution:
                return False
            
            for sol in solution:
                ## 对有微分的方程，先替换微分项
                for func in sol[fun_var].atoms(sp.Derivative):
                    if func.expr in self.solutions:
                        # 替换微分符号为对应的函数
                        import  numpy as np
                        value=np.gradient(self.solutions[func.expr],self.dt)
                        self.solutions[func]=list(value) 
                    # 代入已知值和参数

                    
                    self.solve_equation_var(fun_var, sol[fun_var])
                    eq_data['solved'] = True
                    return True
        except (NotImplementedError, ValueError):
            # 符号求解失败
            return False

    def _solve_single_variable_equations(self):
        """迭代求解单变量方程"""
        queue = deque([eq['id'] for eq in self.equation_db])
        max_iter = 2 * len(self.equation_db)  # 防止无限循环
        solved_count = 0
        
        while queue and solved_count < max_iter:
            eq_id = queue.popleft()
            eq_data = next(eq for eq in self.equation_db if eq['id'] == eq_id)
            
            # 跳过已解方程
            if eq_data['solved']:
                continue
                
            # 更新未知量计数（排除已求解变量）
            current_unknowns = [u for u in eq_data['unknowns'] if u not in self.solutions]
            eq_data['unknowns'] = set(current_unknowns)
            
            success = False
            if len(current_unknowns) == 0:
                # 验证无未知数方程
                success = self._verify_equation(eq_data)
            elif eq_data['type'] == 'ALGEBRAIC' and len(current_unknowns) == 1:
                # 求解单变量代数方程
                success = self._solve_algebraic_equation(eq_data)
            elif eq_data['type'] == 'DIFFERENTIAL' and len(current_unknowns) == 1:
                # 尝试求解微分方程
                success = self._solve_differential_equation(eq_data)
            
            if success:
                # 成功求解，重新检查所有相关方程
                for u in current_unknowns:
                    for related_eq_id in self.unknown_to_eqs.get(u, []):
                        if related_eq_id not in queue:
                            queue.append(related_eq_id)
            else:
                # 放回队列等待后续处理
                queue.append(eq_id)
            
            solved_count += 1

    def _verify_equation(self, eq_data: dict) -> bool:
        """验证无未知数的方程"""
        eq = eq_data['equation']
        # 代入所有已知值和参数
        substituted = eq.subs(self.solutions).subs(self.parameters)
        
        # 简化并验证
        simplified = sp.simplify(substituted)
        if simplified == True:
            eq_data['solved'] = True
            return True
        elif simplified.is_Relational:
            # 处理不等式
            if simplified == False:
                raise ValueError(f"Equation {eq} is not satisfied")
            eq_data['solved'] = True
            return True
        else:
            # 尝试数值验证
            try:
                if abs(simplified.evalf()) < 1e-10:  # 数值容差
                    eq_data['solved'] = True
                    return True
            except TypeError:
                pass
            
            raise ValueError(f"Equation {eq} could not be verified")

    def _cluster_equations(self) -> List[List[dict]]:
        """聚类未解方程（按共享变量分组）"""
        # 获取未解方程
        unsolved_eqs = [eq for eq in self.equation_db if not eq['solved']]
        
        # 创建图结构用于聚类
        graph = {}
        for eq in unsolved_eqs:
            graph[eq['id']] = set()
            for other in unsolved_eqs:
                if eq['id'] != other['id'] and eq['unknowns'] & other['unknowns']:
                    graph[eq['id']].add(other['id'])
        
        # 连通分量分析
        visited = set()
        clusters = []
        
        for eq in unsolved_eqs:
            if eq['id'] not in visited:
                cluster = []
                stack = [eq['id']]
                while stack:
                    node = stack.pop()
                    if node not in visited:
                        visited.add(node)
                        cluster.append(next(eq for eq in unsolved_eqs if eq['id'] == node))
                        stack.extend(graph[node] - visited)
                clusters.append(cluster)
        
        return clusters

    def _solve_equation_cluster(self, cluster: List[dict]):
        """尝试求解方程簇"""
        # 提取所有方程和变量
        equations = [eq['equation'] for eq in cluster]
        variables = set()
        for eq in cluster:
            variables |= eq['unknowns']
            symbols = eq['symbols']
        
        try:
            # 尝试求解方程组
            solutions = sp.solve(equations, *variables, dict=True)
            if not solutions:
                raise ValueError("No solution found")
            
            
            # 处理多解情况
            for sol in solutions:
                for var, expr in sol.items():
                    # 只处理未求解的变量
                    if var not in self.solutions:
                        self.solve_equation_var(var,expr)
            
            # 标记所有方程为已解
            for eq in cluster:
                eq['solved'] = True
                
        except (NotImplementedError, ValueError):
            # 符号求解失败，标记为数值求解
            self.numeric_groups.append({
                'equations': equations,
                'variables': variables,
                'type': 'ALGEBRAIC' if all(eq['type'] == 'ALGEBRAIC' for eq in cluster) else 'DIFFERENTIAL'
            })

    def _solve_remaining_systems(self):
        """求解剩余的方程组"""
        # 聚类未解方程
        clusters = self._cluster_equations()
        
        # 尝试求解每个簇
        for cluster in clusters:
            self._solve_equation_cluster(cluster)

    def solve(self, numeric_solver: Optional[callable] = None,dt=0.1) -> Dict[sp.Symbol, Union[float, sp.Expr]]:
        """
        主求解方法
        
        :param numeric_solver: 可选的数值求解器函数
        :return: 求解结果字典
        """
        # 步骤1: 求解单变量方程
        self.dt=dt
        self._solve_single_variable_equations()
        
        # 步骤2: 求解剩余方程组
        self._solve_remaining_systems()
        
        # 步骤3: 数值求解（如果提供求解器）
        if numeric_solver and self.numeric_groups:
            numeric_solutions = {}
            for group in self.numeric_groups:
                solutions = numeric_solver(group['equations'], group['variables'])
                numeric_solutions.update(solutions)
            
            # 更新数值解
            for var, value in numeric_solutions.items():
                self.solutions[var] = value
        
        # 验证所有方程
        self._verify_all_equations()
        
        return self.solutions

    def _verify_all_equations(self):
        """验证所有方程是否满足"""
        for eq_data in self.equation_db:
            eq = eq_data['equation']
            try:
                # 代入所有已知解和参数
                substituted = eq.subs(self.solutions).subs(self.parameters)
                simplified = sp.simplify(substituted)
                
                if simplified != True and simplified.is_Relational:
                    if simplified == False:
                        raise ValueError(f"Equation {eq} not satisfied after solving")
            except Exception as e:
                print(f"Warning: Verification failed for equation {eq}: {str(e)}")

    def get_unsolved_equations(self) -> List[sp.Eq]:
        """获取未解方程列表"""
        return [eq['equation'] for eq in self.equation_db if not eq['solved']]

    def get_numeric_groups(self) -> List[dict]:
        """获取需要数值求解的方程组"""
        return self.numeric_groups

    def solve_equation_var(self, var: sp.Symbol,expr: Optional[sp.Expr] = None):
        varlist=[] if not self.singleresult else None
                            
        for i in range(self.solutionsNum):
            # 代入已知值和参数
            substituted_dict=  {
                    k: v[i] if isinstance(v, list) else v
                    for k, v in self.solutions.items()
                }
            value = expr.subs(substituted_dict).subs(self.parameters)
            
            # 如果是简单表达式则计算数值
            if value.is_constant():
                value = value.evalf()
            varlist = value if self.singleresult else varlist + [value]
        
        self.solutions[var]=varlist
        self._update_equation_unknowns(var)    


# 示例1: 纯代数方程系统
def test_algebraic_system():
    # 定义符号
    x, y, z = sp.symbols('x y z')
    a, b = sp.symbols('a b')
    
    # 创建方程
    eq1 = sp.Eq(x + y + z, 10)
    eq2 = sp.Eq(2*x - y + 3*z, 5)
    # eq3 = sp.Eq(x + 2*y - z, a)
    
    # 已知变量和参数
    known_vars = {z: [2,2,3]}
    parameters = {a: 1, b: 3}
    
    # 创建求解器
    solver = EquationSolver([eq1, eq2], known_vars, parameters)
    
    # 求解
    solutions = solver.solve()
    print("Algebraic System Solutions:")
    for var, val in solutions.items():
        print(f"{var}: {val}")

def test_step_solve_system():
    # 定义符号
    t = sp.symbols('t')
    C= sp.symbols('C')
    q_0 = sp.Function('q_0')(t)
    e_0 = sp.Function('e_0')(t)
    f_0 = sp.Function('f_0')(t)
    dq_0 = sp.diff(q_0,t)
    
    # 创建方程
    eq1 = sp.Eq(dq_0+q_0+e_0 , f_0)
    eq2 = sp.Eq(q_0 , C * e_0)
    # eq3 = sp.Eq(x + 2*y - z, a)
    
    # 已知变量和参数
    known_vars = {q_0: [2,2,3],t:[0,1,2]}
    parameters = {C:1}
    
    # 创建求解器
    sp.dsolve(eq1, f_0)
    print(eq1.free_symbols)
    print(eq1.atoms(sp.Function))
    solver = EquationSolver([eq1, eq2], known_vars, parameters,independent_var=t)
    
    # 求解
    solutions = solver.solve(dt=0.1)
    print("Algebraic System Solutions:")
    for var, val in solutions.items():
        print(f"{var}: {val}")

# 示例2: 混合代数微分方程系统
def test_mixed_system():

    # 定义符号
    t = sp.symbols('t')
    x = sp.Function('x')(t)
    y = sp.Function('y')(t)
    dx = sp.symbols('dx')
    dy = sp.symbols('dy')
    k = sp.symbols('k')
    
    # 创建方程
    eq1 = sp.Eq(dx, -k * x)
    eq2 = sp.Eq(dy, k * x - 0.5 * y)
    eq3 = sp.Eq(x + y, 100)
    
    # 初始条件
    known_vars = {x.subs(t, 0): 100, y.subs(t, 0): 0}
    parameters = {k: 0.1}
    
    # 创建求解器
    solver = EquationSolver([eq1, eq2, eq3], known_vars, parameters, independent_var=t)
    
    # 求解
    solutions = solver.solve()
    print("\nMixed System Solutions:")
    for var, val in solutions.items():
        print(f"{var}: {val}")

# 运行测试
if __name__ == "__main__":
    test_algebraic_system()
    test_step_solve_system()
    # test_mixed_system()