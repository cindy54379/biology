# Figure 1f：Module 8细胞毒性评分备选版本

- Figure1f_cytotoxic_Module8：单独的f面板（PNG/PDF/SVG）。
- Figure_1_CD8_atlas_and_programs：完整Figure 1，仅f换为Module 8；其余面板数据不变，页脚更新对应显示规则。
- 细胞范围：双抗Nivo+Rela第4周肿瘤CD8细胞，n=7,445，11位患者。
- 使用作者Module_8，不是Novershtern effector参考程序。分数未修改，颜色范围为该群体第2–98百分位（−0.08229235256210185到1.2220180902626834）。点面积3 pt²，与g/h一致；颜色表示模块得分，不是细胞数量或单基因表达。
- 单面板复现：安装原包依赖，运行 `python plot_cytotoxic_panel.py`，读取同目录Cytotoxic_module8_source.csv。
- 完整Figure 1复现：将replot_Figure1_cytotoxic.py放进原ZIP解压数据目录，运行 `python replot_Figure1_cytotoxic.py`，输出到cytotoxic_option_figures。
- 这是备选图，原ZIP中默认的Novershtern effector版本未被替换。
