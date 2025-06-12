import json
import BondGraphTools as bgt

def build_model_from_json(json_data):
    """
    根据 JSON 配置构建复合元件模型
    json_data 格式示例：
    {
      "name": "ThrottleValve",
      "components": [
        {"type": "R", "name": "FlowR", "params": {"value": 0.6}},
        {"type": "R", "name": "LeakR", "params": {"value": 1e-8}},
        {"type": "0", "name": "Junction"}
      ],
      "connections": [
        ["FlowR", "Junction"],
        ["LeakR", "Junction"]
      ],
      "exposed_ports": {
        "FlowR": "in",
        "LeakR": "out"
      }
    }
    """
    # 创建复合模型
    model_name = json_data.get("name", "CompositeModel")
    model = bgt.new(name=model_name)

    # 先创建所有元件，存进字典便于引用
    components = {}
    for comp in json_data.get("components", []):
        ctype = comp.get("type")
        cname = comp.get("name")
        params = comp.get("params", {})
        # bgt.new 的关键字参数，value等可直接传
        components[cname] = bgt.new(ctype, name=cname, **params)
        bgt.add(model, components[cname])

    # 建立连接
    for conn in json_data.get("connections", []):
        # 连接两个或多个元件端口，简化默认连接全部端口
        elems = [components[name] for name in conn]
        bgt.connect(*elems)

    # 暴露端口
    for internal_name, alias in json_data.get("exposed_ports", {}).items():
        if internal_name in components:
            internal_component = components[internal_name]
            # 暴露该元件的第一个端口
            # 如果需要暴露具体端口，可在 JSON 中扩展字段
            bgt.expose(model, internal_component.ports[0], alias=alias)

    return model


# file: subsystem_builder.py

import BondGraphTools as bgt

def build_composite_component(spec: dict):
    system = bgt.new(name="Composite")

    components = {}
    for name, desc in spec.get("subcomponents", {}).items():
        comp = bgt.new(desc["type"], name=name, value=desc.get("value", {}))
        components[name] = comp
        bgt.add(system, comp)

    for connection in spec.get("connections", []):
        elems = [components[name] for name in connection]
        bgt.connect(*elems)

    for comp_name, port_alias in spec.get("expose_ports", {}).items():
        comp = components[comp_name]
        from BondGraphTools.actions import _unpack_port_arg
        component, port=_unpack_port_arg(comp)  # 确保端口被正确解析
        port = comp.ports.get("e") or next(iter(comp.ports.values()))
        bgt.expose(system, port, alias=port_alias)

    return system

# 用法示例
if __name__ == "__main__":
    # 从文件读配置（可根据需要替换为你的json路径）
    import os
    with open(os.path.join(os.path.dirname(__file__), "RC.json"), "r") as f:
        data = json.load(f)

    # model = build_model_from_json(data)
    # print(model)

    
    rc_block = build_composite_component(data["components"]["RCBlock"])

    # 用法示例
    model = bgt.new()
    source = bgt.new("Se", value=1.0)
    bgt.add(model, rc_block, source)
    bgt.connect(source, rc_block)
    bgt.draw(model)
