# DP-GEN 配置、验证与执行

在创建、审阅、验证、启动或恢复 `dpgen run` 配置时加载本 reference。

## 1. `param.json` 要点

应保留或收集：

- `type_map`（有序元素符号）、可选的 `mass_map`、`init_data_prefix` 和 `init_data_sys`；
- `sys_configs_prefix` 与 `sys_configs`（系统列表，每个系统又是结构列表）；
- `numb_models` 以及完整的 `default_training_param.model`、`learning_rate`、`loss` 和 `training`；
- `model_devi_dt`、`model_devi_skip`、信任阈值、清理标志和 `model_devi_jobs`；
- `fp_style`、`fp_task_min`、`fp_task_max` 及后端专用 FP 输入。

DeePMD-kit 1.x 使用 `training.systems`，2.x/3.x 使用
`training.training_data.systems`。将对应版本的 systems 留空，由 DP-GEN 根据
`init_data_sys` 填充。使用版本相关 descriptor 或功能前，核对 `deepmd_version`
和已安装的软件栈。

`train_backend` 支持 `tensorflow` 和 `pytorch`，但实际兼容性仍取决于已安装的
DeePMD-kit 与 descriptor。核对组合并准确拼写 `se_atten_v2`。普通 LAMMPS 探索使用
默认 `model_devi_engine`；替代引擎只有在验证后才能设置。当前 `dpgen run` 不接受
`fp_style: "none"`。

每个 `init_data_sys` 目录必须包含 `type_map.raw`、`type.raw` 和 `set.000/`；
每个 `type_map.raw` 必须与 `param.json.type_map` 完全一致。

## 2. `machine.json` 要点

保持独立的 `train`、`model_devi` 和 `fp` block。每个 block 提供 `command`、
`machine`（`batch_type`、`context_type`、`local_root` 以及需要时的
`remote_root`/profile）和 `resources`（节点数、CPU/GPU 数、`group_size`、
队列/分区参数及调度器 `source_list`）。通用 local、scheduler 和 SSH 结构见
[machine schema 说明](../../../doc/run/example-of-machine.md)。

启动 shell 与派发任务处于不同环境。调度器任务必须通过 `source_list` 激活自己的
运行环境，外层激活不会可靠继承。训练命令必须选择与 `train_backend` 相同的后端
（DeePMD-kit 3.x 的 PyTorch 可能需要 `dp --pt`）。应根据已安装的 DPDispatcher
验证命令、context、batch 别名、路径和资源，不要凭经验统一改写。

## 3. 启动前验证生成输入

```bash
dpgen -h
python -m json.tool param.json
python -m json.tool machine.json
```

用 `run_jdata_arginfo()` 规范化 `param.json`，用 `convert_mdata()` 转换
`machine.json`。解析所有数据/结构路径，比较 type map，检查可执行文件和调度器权限，
并核对 FP 输入及成本上限。

提交前逐个检查每个阶段/任务实际生成的 `input.json`。不能从根目录 JSON 推断最终配置：
核对模型数、seed、checkpoint、数据系统、`sys_idx`、温度/压力和最终训练步数。存在
reuse 或迭代 override 时，同时断言默认值和生效值。

将目标迭代与 `record.dpgen`、阶段目录、DPDispatcher work base/提交元数据及日志交叉核对。
仅有目录或 job 并不代表阶段已完成。

## 4. 恢复与重启

从 `record.dpgen` 恢复属于可能产生费用的操作。同一工作目录只能保留一个 controller。
编辑状态前先停止 controller、备份状态文件、记录阶段/job ID 及新旧设置，并查询调度器。

活动阶段的 command、context、remote root、队列、资源、`group_size`、flags、激活方式、
任务列表或 common files 变化，都可能产生新的 DPDispatcher submission identity。新 identity
不会导入旧 identity 的完成状态。要区分旧提交恢复和使用新 remote root 的有意重跑；连接超时
本身不是重新提交理由。

恢复的正证据包括旧 submission hash、未变化的 job ID 和恢复日志；同时不能出现新的提交日志。
`record.dpgen` 是工作流指针，不是调度器完成记录。当前计费阶段推进后，才能应用后续迭代的
machine 修改。

## 5. 确认与检查

展示文件、验证摘要、未解决风险和 `dpgen run param.json machine.json`。只有明确确认后才执行。
启动后检查 `iter.*`、阶段日志、失败/等待状态、选中和标注数量，以及
[监控输出](monitoring.md)。

外部参考：https://docs.deepmodeling.com/projects/dpgen/en/latest/run/index.html
