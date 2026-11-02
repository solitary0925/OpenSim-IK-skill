# OpenSim Vicon IK Skill

**v1.4.0**：从动态 Vicon C3D 和已缩放的原膝 CGM24 模型运行 OpenSim 逆运动学（IK），保留配置、原始帧号、误差和复核记录。可通过 skill 调用，也可直接运行 Python 脚本。

## 功能

- 统一 JSON 配置和运行前检查：路径、标记对应、单位、完整性、时间窗、权重、坐标旋转和排除清单。
- 输入受阻时单独记录该试次，其他正常试次继续运行；环境、模型或全局配置错误时停止整批。
- 输出缺失点名、所有受影响标记的缺失帧数和连续缺失范围；分别记录整条数据及所选窗口。
- 在工作模型副本中锁定无标记观测的 MTP/腰椎坐标，保留原膝单旋转自由度和可观测的踝/距下关节坐标。
- 从保存的 IK MOT 独立重建标记误差，与 OpenSim 原生误差比较。
- 定位异常在整条试次中的位置，例如“第10–20帧 / 共200帧”，另存原始 C3D 帧号和时间映射。
- 为待复核试次生成整条数据的误差时间线及模型点/实验点叠加图，默认最多3个叠加帧。
- 可配置后续分析保留区间，分别统计区间内外异常及区间内 RMS/Max，保留全段证据。
- 分开记录输入受阻、用户确认排除、求解失败、已求解但待复核及未触发数值筛查的试次。

## 使用

需要可导入 `opensim`、`numpy` 的 Python 环境；开启绘图时还需要 `matplotlib`。已验证的环境为 Python 3.11.14、OpenSim Python API 4.5.2、NumPy 1.25.2。请使用实际安装了这些依赖的 Python 解释器。

1. 下载本仓库，或克隆为 `opensim-vicon-ik` 文件夹。
2. 准备自己的已缩放兼容模型、动态 C3D 和可选的显示几何文件。
3. 复制 [配置示例](assets/config.example.json) 为 `my_study.json`，填写路径和设置。相对路径以配置文件所在目录为基准。
4. 先检查，再运行：

```sh
python scripts/run_config.py --config my_study.json --check-only --report-dir preflight_report
python scripts/run_config.py --config my_study.json
```

批处理使用 `data` 指定动态 C3D 所在目录；单试次使用 `c3d` 指定文件，两者选一个。批处理只扫描目录顶层，跳过文件名包含 `static` 的 C3D。标记 IK 不要求 VSK 或测力数据。

检查中 `ready=false` 表示存在受阻项；若 `can_run=true`，正常运行命令仍会处理通过检查的试次。混合批次完成后返回码2表示有受阻或求解失败试次，并不表示全部失败。返回码0也可能包含需要数值复核的结果。

已有结果可只更新复核报告，无需重跑 IK：

```sh
python scripts/run_config.py --config my_study.json --review-only
```

完整工作流见 [SKILL.md](SKILL.md)，配置字段、状态定义、时间映射和输出格式见 [配置与输出说明](references/configuration.md)。

## 主要输出

- 每个成功试次：工作模型、`markers.trc`、`ik.mot`、IK setup、日志、原生/独立误差、`metrics.json`、`review.json` 和中文报告。
- 时间对应：`timing.json`、`time_mapping.csv`；所选窗口的 adapter 时间不重新归零。
- 异常复核：`异常时间段.csv`、`异常帧可视化/`，包含原始帧号与整条数据中的位置。
- 输入受阻：`blocked_inputs.json`、`输入受阻明细.md`、`输入受阻帧范围.csv`。
- 批次记录：`batch_metrics.json`、`review_status.json`、`trial_selection.json`、`核查汇总.md`。

## 适用范围

适用于当前已验证的原膝 CGM24 模型及兼容标记定义；不进行缩放、ID、RRA、SO、CMC 或 JRA。默认坐标转换为模型 XYZ = 实验室(-X, Z, Y)，更换实验室或会话时须确认坐标约定。

程序不自动插值、交换标签、移动模型点位或裁剪数据。标记缺失不会自动变成用户废弃决定；历史结果保留，但当前受阻试次不能凭旧MOT被认定为本次有效结果。缺失点在已明确选择的求解窗口之外时，不阻止该有效窗口运行。

RMS>20 mm、Max>40 mm 和邻帧角度变化>20°用于数值筛查，仍需结合动作、点位、图像和有效分析区间复核。未触发数值筛查不代表已经完成解剖验证。帧号时间不等于墙钟或已验证的外部测力同步；程序不推断力台接触事件或完整步态周期。

仓库仅包含通用 skill、脚本和示例配置，不包含受试者数据、人体模型、OpenSim 软件、本机配置或实验输出。`PACKAGE_MANIFEST.json` 保存发布的 skill 文件 SHA256；README 和仓库忽略规则属于仓库附加文件。
