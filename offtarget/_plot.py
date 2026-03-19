import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import sys
# 读取数据文件
def readdata(filename):
    genes = []
    values = []
    lastfourcols = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # 分割行
            parts = line.split(':')
            
            # 提取基因名称 (第一部分按分割的第一部分)
            genename = parts[0].split('_')[0]

            
            # 提取第8列值 (索引7)
            try:
                value = float(parts[7])
            except (IndexError, ValueError):
                continue
                
            # 提取最后四列
            if len(parts) >= 4:
                lastfour = parts[-4:]
            else:
                lastfour = [''] * 4
                
            genes.append(genename)
            values.append(value)
            lastfourcols.append(lastfour)
    
    return genes, values, lastfourcols

# 创建PDF报告
def create_report(genes, values, lastfourcols, output_pdf):
    # 转换为numpy数组便于排序
    values = np.array(values)
    
    # 获取排序索引（升序，最小的在前）
    sorted_indices = np.argsort(values)
    
    # 获取前15个最小值的基因
    top15_indices = sorted_indices[:15]
    top15_genes = [genes[i] for i in top15_indices]
    top15_values = [values[i] for i in top15_indices]
    top15_lastfour = [lastfourcols[i] for i in top15_indices]
    
    # 创建PDF - 显著增加高度到20英寸
    with PdfPages(output_pdf) as pdf:
        # 创建图形和左右布局 - 大幅增加高度
        fig = plt.figure(figsize=(16, 20))  # 高度从12增加到20英寸
        gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1], height_ratios=[1])  # 明确指定高度比例
        
        # 左侧：Rank Plot
        ax1 = plt.subplot(gs[0])
        
        # 绘制所有点的rank plot
        ranks = np.arange(1, len(values) + 1)
        ax1.plot(ranks, values[sorted_indices], 'b-', alpha=0.5, linewidth=1.5)
        
        # 标记前15个最小的点
        ax1.scatter(ranks[:15], values[sorted_indices][:15], color='red', s=80, zorder=5)
        
        # 添加前15个基因的标签 - 增加字体大小和偏移量
        for i in range(15):
            ax1.annotate(top15_genes[i], 
                         (ranks[i], values[sorted_indices][i]),
                         textcoords="offset points", 
                         xytext=(10, 10),  # 增加偏移量
                         ha='left',
                         fontsize=12,     # 增加字体大小
                         bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
        
        # 设置图表标题和标签 - 增加字体大小和间距
        ax1.set_title('Off-target Ranking Plot', fontsize=18, pad=30)
        ax1.set_xlabel('Genes Rnak', fontsize=16, labelpad=15)
        ax1.set_ylabel('MFE(Lower values, higher off-target risk)', fontsize=16, labelpad=15)
        ax1.grid(True, linestyle='--', alpha=0.7)
        
        # 调整坐标轴刻度标签大小
        ax1.tick_params(axis='both', which='major', labelsize=12)
        
        ax2 = plt.subplot(gs[1])
        ax2.axis('off')
        
        # 准备表格数据 - 保留所有空格
        # 准备表格数据 - 保留所有空格
        cell_text = []
        for i in range(15):
            # 将最后四列数据组合成一个单元格，每行一个值
            # 使用非断空格(\u00A0)替换普通空格，确保空格被保留
            processed_data = [s.replace(' ', '\u00A0') for s in top15_lastfour[i]]
            
            # 在不同行前面添加不同的前缀
            for j in range(len(processed_data)):
                if j == 1:  # 第二行（索引从0开始）
                    processed_data[j] = "RNA:" + processed_data[j]
                elif j == 2:  # 第三行（索引从0开始）
                    processed_data[j] = "ASO:" + processed_data[j].replace('U', 'T')
                else:
                    processed_data[j] = "\u00A0\u00A0\u00A0\u00A0" + processed_data[j]  # 3个非断空格
            
            cell_data = '\n'.join(processed_data)
            cell_text.append([cell_data])

        # 创建表格
        table = ax2.table(
            cellText=cell_text,
            rowLabels=top15_genes,
            colLabels=['Last Four Columns'],
            loc='center',
            cellLoc='center'
        )
        
        # 设置表格样式 - 显著增加行高，保持字体大小不变
        table.auto_set_font_size(False)
        
        # 设置等宽字体 - 这是关键部分
        monospace_font = {'family': 'monospace', 'size': 12}
        table.set_fontsize(12)
        
        # 应用等宽字体到所有单元格
        for (i, j), cell in table.get_celld().items():
            cell.set_text_props(fontproperties=monospace_font)
            # 设置单元格对齐方式为居中
            cell.set_text_props(ha='center', va='center')
        
        # 设置表头字体
        for (i, j), cell in table.get_celld().items():
            if i == 0:  # 表头行
                cell.set_text_props(weight='bold')
        
        # 调整表格缩放
        table.scale(1, 6)         # 显著增加行高
        
        # 设置表格标题 - 增加字体大小和间距
        ax2.set_title('Top 15 Off-target Genes', fontsize=18, pad=40)  # 增加pad值
        
        # 调整布局 - 增加边距
        plt.tight_layout(rect=[0, 0.05, 1, 0.95])  # 添加边距
        
        # 保存到PDF
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()


# 主程序
if __name__ == "__main__":
    input_file = sys.argv[1]
    output_pdf = sys.argv[2]
    
    genes, values, lastfourcols = readdata(input_file)
    create_report(genes, values, lastfourcols, output_pdf)
    
    print(f"报告已生成: {output_pdf}")
