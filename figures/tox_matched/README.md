# 与Figure 1g/h统一格式的TOX图

TOX_matched为双抗第4周（7,445细胞），TOX_baseline_matched为基线。点面积从旧TOX图的2 pt²改为与g/h相同的3 pt²，并统一实际坐标轴尺寸、颜色、字体、点排序和色条版式。移除原TOX图的cluster文字以匹配g/h。TOX数值及原色阶上下限不变，基线/第4周TOX共用0–2.497504670715671。色条名称为Normalized expression。

PDCD1_HAVCR2_TOX_matched.pdf提供三个同尺寸面板便于对照。插入论文时请保持三个面板相同宽高，避免单独缩放导致点大小看起来再次不同。

复现：使用主绘图包的Python依赖，运行 `python plot_tox_matched.py`。该脚本直接读取同目录的TOX_matched_source.csv。图像保存到outputs目录。参数见Plot_parameters.json。未重新聚类或改变任何表达值。
