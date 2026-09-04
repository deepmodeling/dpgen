# DP-GEN Run 工作流说明

在需要了解范围、准备策略、项目布局或结果汇报时加载本 reference。

## 循环

`dpgen run` 会重复以下五个活动：

1. 训练一组 Deep Potential 模型；
1. 通常使用 LAMMPS MD 探索构型空间；
1. 选取处于 model-deviation 候选区间的构型；
1. 使用配置的第一性原理后端标注选中构型；
1. 将标注加入下一轮训练。

`param.json` 保存科学/工作流设置，`machine.json` 保存命令、上下文和资源，
启动命令为 `dpgen run param.json machine.json`。准备阶段需要 Python 和 DP-GEN；
执行阶段还需要 DeePMD-kit、探索引擎、选定的 FP 软件和调度器运行环境。

## 准备策略

先检查已有配置、训练输入、数据集元数据和 machine 模板。只询问无法发现的值，
修改现有文件而不是无必要地重建。

保留 descriptor 与 fitting 设置、训练和 FP 后端、阈值、`type_map` 顺序、
`model_devi_jobs`、ensemble 数量、温度、压力和 MD ensemble。发现问题时先说明并
请求方向，不要静默修改。

原有的激活命令、可执行文件、队列、分区、账号、路径和调度器参数必须原样复用。
不得猜测站点策略、conda 环境或 module 配置。

## 通用布局

```text
project/
|-- param.json
|-- machine.json
|-- init_data/
|   `-- raw_xxx/          # DeepMD NumPy system(s)
|-- assets/               # 构型和运行所需支持文件
`-- iter.*/               # 由 DP-GEN 创建
```

重复实验应放在独立且命名清晰的目录中，并从审阅过的基础配置派生。后端专用文件
放在 `assets/` 或选定后端明确要求的路径下。

## 汇报契约

执行前汇报两个 JSON 的绝对路径、精确命令、验证结果、未解决输入和计费阶段。只有用户
确认该精确命令后才执行。

执行后汇报当前迭代和阶段状态、失败/等待任务、主要日志和输出、候选与标注数量、
准确率趋势以及下一步要检查的文件。迭代停滞或准确率不提升时，使用
[监控说明](monitoring.md)。

## 防护规则

- 两个 JSON 都存在且通过验证前不得运行。
- 数据、训练和标注过程中的 `type_map` 必须一致。
- 不得盲目覆盖可用模板，也不能假设外层激活会传递到任务环境。
- 无法发现缺失的科学或站点值时，停止而不是猜测。

官方概览：https://docs.deepmodeling.com/projects/dpgen/en/latest/run/index.html
