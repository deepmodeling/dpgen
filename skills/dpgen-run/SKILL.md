---
name: dpgen-run
description: 使用 param.json 和 machine.json 准备、解释、验证并运行 DP-GEN 并行学习流程。
license: LGPL-3.0-or-later
metadata:
  author: MatMaster
  version: 0.4.0
  repository: https://github.com/deepmodeling/dpgen
---

# DP-GEN Run

用于迭代训练、探索、标注和再训练流程。

## 使用契约

- 科学参数放在 `param.json`，执行参数放在 `machine.json`。
- 精确运行命令为 `dpgen run param.json machine.json`。
- 默认只准备和验证；只有用户明确确认后才执行。

## 工作流程

1. 先检查已有数据和配置，再询问缺失输入。
1. 修改现有可用文件，并保留科学及站点选择。
1. 只加载当前任务需要的 reference。
1. 验证 schema、生成的任务输入、阶段状态和监控信号。
1. 汇报路径、命令、分体系证据、风险和后续检查。

## Reference

- [工作流说明](references/workflow.md)：范围、布局、策略和汇报。
- [配置、验证与执行](references/validation-and-run.md)：两个 JSON、生成输入检查及安全启动/恢复。
- [监控与排障](references/monitoring.md)：准确率趋势、阶段进度和证据化诊断。

## 防护规则

- 不得合并两个 JSON，也不得臆造路径、队列、模块或命令。
- 保留 descriptor、backend、阈值、计划、ensemble 和 `type_map` 选择。
- 必须准确写成 `se_atten_v2`；当前 `dpgen run` 不接受 `fp_style: "none"`。
- 无法安全发现必需的科学或执行输入时，停止并询问。
