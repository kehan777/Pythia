#####################
import os
import pandas as pd
import plotly.express as px
import numpy as np
import plotly.colors

# 读取数据并保持顺序
positions = []
data = []
#>>>>>>>>>>>>>>>>>>>>>>>>>>
with open('D:/conda/jupyter/Pythia/inputs/obj01_pred_mask.txt', 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) == 2:
            mutation = parts[0]
            value = float(parts[1])
            # 解析突变信息
            wt = mutation[0]  # 野生型氨基酸
            pos = mutation[1:-1]  # 位置（保持字符串形式）
            mut = mutation[-1]  # 突变型氨基酸
            pos_label = f"{wt}{pos}"  # 位置标签
            
            #>>> 跳过 C、G、P 的突变
            if mut in ['C', 'G', 'P']:
                continue
                
            # 记录位置顺序
            if pos_label not in positions:
                positions.append(pos_label)
                
            data.append([pos_label, mut, value])

# 转换为DataFrame
df = pd.DataFrame(data, columns=['Position', 'Mutation', 'Energy'])

# >>>> 定义要显示的位置区间
selected_positions = []
for pos in positions:
    pos_num = int(''.join(filter(str.isdigit, pos)))  # 提取位置数字
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if (28 <= pos_num <= 35) or (49 <= pos_num <= 67) or (100 <= pos_num <= 111):
        selected_positions.append(pos)

# 过滤数据，只保留选定的位置
df = df[df['Position'].isin(selected_positions)]

# 创建矩阵（保持原始顺序）
amino_acids = sorted(df['Mutation'].unique())  # 突变氨基酸按字母顺序
matrix = pd.DataFrame(index=amino_acids, columns=selected_positions)  # 保持原始位置顺序

# 填充矩阵，并保留2位有效数字
for _, row in df.iterrows():
    matrix.at[row['Mutation'], row['Position']] = round(row['Energy'], 2)

# 获取数据的最小值和最大值
zmin = np.nanmin(matrix.values)
zmax = np.nanmax(matrix.values)

# 自定义颜色比例
# 计算0在颜色比例中的位置
zero_pos = (0 - zmin) / (zmax - zmin)

# 创建渐变色
# 蓝色渐变：绝对值越大，颜色越深
blue_gradient = plotly.colors.sample_colorscale('Blues', np.linspace(1, 0, 5))  # 反转渐变方向
# 红色渐变：值越大，颜色越深
red_gradient = plotly.colors.sample_colorscale('Reds', np.linspace(0.01, 1, 5))  # 红色渐变

# 构建颜色比例
custom_colorscale = []
# 蓝色渐变部分（从最深到最浅）
for i, color in enumerate(blue_gradient):
    pos = i / (len(blue_gradient) - 1) * zero_pos
    custom_colorscale.append([pos, color])
# 红色渐变部分
for i, color in enumerate(red_gradient):
    pos = zero_pos + i / (len(red_gradient) - 1) * (1 - zero_pos)
    custom_colorscale.append([pos, color])

# 使用plotly创建交互式热图
fig = px.imshow(matrix.astype(float),
                labels=dict(x="Wildtype Amino Acid", 
                           y="Mutated Amino Acid", 
                           color="Energy Change"),
                color_continuous_scale=custom_colorscale,
                aspect="auto",
                zmin=zmin,
                zmax=zmax,
                text_auto=True)  # 自动显示数值

# 更新布局
fig.update_layout(
    title='Protein Mutation Energy Changes',
    xaxis_nticks=len(selected_positions),  # 显示所有x轴刻度
    yaxis_nticks=len(amino_acids)  # 显示所有y轴刻度
)

# 保存为HTML文件
output_path = os.path.join(os.getcwd(), 'mutation_energy_heatmap_selected.html')
fig.write_html(output_path)

print(f"图表已保存到：{output_path}")

# 显示图表（可选）
fig.show()
