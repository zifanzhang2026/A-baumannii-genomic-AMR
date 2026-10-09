import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --- 设置全局字体为 Times New Roman ---
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['axes.unicode_minus'] = False 

# 1. 读取数据
df = pd.read_csv(r'D:\A.baumannii\data\metadata.csv')

# 2. 数据预处理
years = df['year'].dropna()
years = pd.to_numeric(years, errors='coerce').dropna()
years = years[(years >= 1900) & (years <= 2025)].astype(int)

# 3. 统计每年的频数 (Frequency)
yearly_counts = years.value_counts().sort_index()

# 4. 计算 Y 轴的值: log10(Frequency + 1)
y_values = np.log10(yearly_counts + 1)
x_values = yearly_counts.index

# 5. 开始绘图
plt.figure(figsize=(9, 6))

# 绘制散点图
plt.scatter(x_values, y_values, color='#ED7D31', edgecolor='#F8CBAD', s=50, zorder=3)

# 设置坐标轴标签 (这里将 fontsize 调大到了 14)
plt.xlabel('Period', fontsize=14)
plt.ylabel('Genomes deposited(log(Frequency+1))', fontsize=14)

# ==========================================
# 在这里修改坐标轴数字的大小！(fontsize=12)
# ==========================================
plt.xticks(np.arange(1900, 2026, 20), fontsize=12)
plt.yticks(np.arange(0, 5, 1), fontsize=12)

# 获取当前坐标轴
ax = plt.gca()

# 去除右侧和顶部的边框线
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# 去除坐标轴上的刻度短线 (仅保留数字)
ax.tick_params(axis='both', bottom=False, left=False)

# 调整边距并显示图表
plt.tight_layout()
plt.savefig('散点图.png', dpi=300, bbox_inches='tight')
plt.show()
