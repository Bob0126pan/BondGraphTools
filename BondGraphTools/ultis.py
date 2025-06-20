def set_parameters(spec, params, prefix=""):
    """递归设置参数到组件规范"""
    if not params:
        return
        
    # 处理当前层级的子组件
    for sub_name, sub_spec in spec.get("subcomponents", {}).items():
        sub_prefix = f"{prefix}{sub_name}."
        
        # 收集属于当前子组件的参数
        sub_params = {}
        for key in list(params.keys()):
            if key.startswith(sub_prefix):
                # 提取相对参数名（去掉子组件前缀）
                rel_key = key[len(sub_prefix):]
                sub_params[rel_key] = params.pop(key)
        
        # 如果有属于此子组件的参数，递归设置
        if sub_params:
            set_parameters(sub_spec, sub_params, prefix="")
    
    # 设置当前组件的直接参数
    if prefix == "":  # 只在最顶层处理直接参数
        for key, val in params.items():
            # if '.' not in key:  # 直接参数（无点号）
            #     if "value" not in spec:
            #         spec["value"] = {}
            spec["value"][key] = val