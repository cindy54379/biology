# 黑色素瘤 Nivo + Rela 单细胞再分析绘图包

[下载完整ZIP：NivoRela_Replot.zip](https://github.com/cindy54379/biology/raw/refs/heads/main/downloads/NivoRela_Replot.zip)

ZIP约12 MB，包含重绘所需的全部CSV数据、JSON参数、Python代码、固定依赖版本、中英文说明及复现验证记录。无需额外下载GEO对象，无需原工作区的中间文件。

## 复现

建议Python 3.12。解压全部文件至同一文件夹，保持原文件名，在该文件夹打开终端：

```bash
python -m pip install -r requirements_replot.txt
python replot_from_csv.py
```

输出位于 `reproduced_figures/`，每套图包含PNG、PDF、SVG。

[详细中文操作说明](REPRODUCE_CN.txt)

## 图件范围

- Figure 1：CD8图谱、程序及标志基因，含TOX。
- Figure 2：拟时序第10箱中的耗竭/细胞毒性共表达。
- Figure 3：独立临床研究的发表汇总结果。
- Figure S1：患者和阈值敏感性。
- Figure S2：全部肿瘤CD8细胞。
- TOX分布及Effector显示方式核查图。

已从ZIP独立解压并运行，当前环境中7套PNG均与交付图一致，记录见包内Validation.csv。其他系统的字体及软件差异可能造成轻微排版变化。

这是使用作者处理后坐标、cluster和评分的再绘图，不是从FASTQ独立重新计算所有单细胞分析。Figure 2主图仍为第10箱，全体细胞版本为Figure S2。临床文件为发表汇总值，非患者级随访记录。

## 来源

- Cell 2024：doi:10.1016/j.cell.2024.06.036；GSE225063；作者Zenodo记录11222247。
- Nature Medicine 2026：doi:10.1038/s41591-026-04756-1，RELATIVITY-047。
- JCO 2025：doi:10.1200/JCO-25-00494，NCT02519322。

逐文件来源及定义见ZIP中的Source_manifest.csv、Field_dictionary.csv及图注。

文件完整性校验见[SHA256SUMS](downloads/SHA256SUMS.txt)。
