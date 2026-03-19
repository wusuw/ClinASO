#!/usr/bin/env python3
import sys
import os
from reportlab.lib.pagesizes import letter, A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch, cm
from reportlab.lib import colors
from reportlab.graphics.shapes import Drawing, Rect, String, Line, Polygon
from reportlab.graphics import renderPDF
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY
from reportlab.pdfgen import canvas
import math

# 定义网页配色方案
WEB_COLORS = {
    'primary_dark': '#1a3a6c',      # 深蓝
    'primary_medium': '#2c5e92',    # 中蓝
    'primary_light': '#3a7bb8',     # 浅蓝
    'background_light': '#f5f7fa',  # 背景浅色
    'background_medium': '#e4edf5', # 背景中色
    'accent_green': '#2e7d32',      # 成功绿色
    'accent_red': '#c62828',        # 错误红色
    'text_dark': '#333333',         # 深色文字
    'text_light': '#666666',        # 浅色文字
    'white': '#ffffff',             # 白色
}

def reverse_complement(seq):
    """生成反向互补序列"""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
            'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
            'N': 'N', 'n': 'n'}
    return ''.join(comp.get(base, base) for base in reversed(seq))

def parse_vcf_file(filename):
    """解析VCF文件"""
    variants = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 8:
                variant = {
                    'CHROM': fields[0],
                    'POS': int(fields[1]),
                    'ID': fields[2],
                    'REF': fields[3],
                    'ALT': fields[4],
                    'QUAL': fields[5],
                    'FILTER': fields[6],
                    'INFO': fields[7]
                }
                if len(fields) > 8:
                    variant['FORMAT'] = fields[8]
                    variant['SAMPLES'] = fields[9:]
                variants.append(variant)
    return variants

def parse_frequency_vcf(filename):
    """解析频率VCF文件"""
    variants = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 9:
                variant = {
                    'CHROM': fields[0],
                    'POS': int(fields[1]),
                    'ID': fields[2],
                    'REF': fields[3],
                    'ALT': fields[4],
                    'QUAL': fields[5],
                    'FILTER': fields[6],
                    'INFO': fields[7],
                    'FORMAT': fields[8],
                    'SAMPLES': fields[9:]
                }
                variants.append(variant)
    return variants

def hex_to_color(hex_color):
    """将十六进制颜色转换为ReportLab颜色"""
    hex_color = hex_color.lstrip('#')
    r = int(hex_color[0:2], 16) / 255.0
    g = int(hex_color[2:4], 16) / 255.0
    b = int(hex_color[4:6], 16) / 255.0
    return colors.Color(r, g, b)

def create_sequence_drawing(aseq, chrom, start_pos, end_pos, snps, page_width, strand):
    """创建序列图 - 直接在PDF中绘制"""
    seq_length = len(aseq)
    
    # 根据页面宽度计算每个碱基的宽度
    available_width = page_width - 2 * cm
    base_width = min(15, available_width / seq_length)
    
    # 计算总宽度
    total_width = seq_length * base_width + 2 * cm
    
    # 创建绘图对象
    drawing = Drawing(total_width, 150)  # 调整高度
    
    # 颜色定义
    primary_dark = hex_to_color(WEB_COLORS['primary_dark'])
    primary_medium = hex_to_color(WEB_COLORS['primary_medium'])
    accent_red = hex_to_color(WEB_COLORS['accent_red'])
    
    # 添加链方向指示器
    strand_label = String(total_width/2, 130, f"Target Strand: {strand}", fontSize=10, 
                         textAnchor="middle", fillColor=primary_dark)
    drawing.add(strand_label)
    
    # 添加起始位置标签（只在第一个碱基位置）
    start_label = String(cm, 115, f"{chrom}:{start_pos}", fontSize=9, 
                        textAnchor="start", fillColor=primary_dark)
    drawing.add(start_label)
    
    # 添加结束位置标签（只在最后一个碱基位置）
    end_label = String(total_width - cm, 115, f"{chrom}:{end_pos}", fontSize=9, 
                      textAnchor="end", fillColor=primary_dark)
    drawing.add(end_label)
    
    # 绘制序列线
    line = Line(cm, 110, total_width - cm, 110, strokeColor=primary_dark, strokeWidth=1)
    drawing.add(line)
    
    # 绘制每个碱基的方块和标签
    for i, base in enumerate(aseq):
        x = cm + i * base_width
        y = 80  # 碱基方块的位置
        
        # 绘制碱基方块 - 使用网页配色
        rect = Rect(x, y, base_width-2, 20, 
                   fillColor=hex_to_color(WEB_COLORS['background_medium']), 
                   strokeColor=primary_dark, strokeWidth=0.5)
        drawing.add(rect)
        
        # 添加碱基文本（居中）
        text_x = x + (base_width-2)/2
        text_y = y + 10
        base_text = String(text_x, text_y, base, fontSize=8, textAnchor="middle",
                          fillColor=primary_dark)
        drawing.add(base_text)
        
        # 检查是否有SNP在这个位置
        pos_in_genome = start_pos + i
        for snp in snps:
            if snp['POS'] == pos_in_genome:
                # 在碱基下方标记SNP（红色三角形）
                marker_x = text_x
                marker_y = y - 15
                
                # 绘制红色三角形标记
                marker_size = 5
                triangle = Polygon(
                    points=[marker_x, marker_y, 
                           marker_x - marker_size, marker_y - marker_size,
                           marker_x + marker_size, marker_y - marker_size],
                    fillColor=accent_red,
                    strokeColor=accent_red,
                    strokeWidth=0.5
                )
                drawing.add(triangle)
                
                # 在标记下方添加SNP ID（只显示ID，不显示其他信息）
                snp_id = snp['ID']
                if len(snp_id) > 12:
                    snp_id = snp_id[:10] + "..."
                
                snp_label = String(text_x, y - 25, snp_id, fontSize=6, 
                                  textAnchor="middle", fillColor=accent_red)
                drawing.add(snp_label)
                break
    
    return drawing

def format_frequency_table(variants):
    """格式化频率数据表格"""
    population_map = [
        {"name": "European", "group": "Sub"},
        {"name": "African Others", "group": "Sub"},
        {"name": "East Asian", "group": "Sub"},
        {"name": "African American", "group": "Sub"},
        {"name": "Latin American 1", "group": "Sub"},
        {"name": "Latin American 2", "group": "Sub"},
        {"name": "Other Asian", "group": "Sub"},
        {"name": "South Asian", "group": "Sub"},
        {"name": "Other", "group": "Sub"},
        {"name": "African", "group": "Sub"},
        {"name": "Asian", "group": "Sub"},
        {"name": "Total", "group": "Global"}
    ]
    
    tables = []
    
    for variant in variants:
        title = f"{variant['ID']} - {variant['CHROM']}:{variant['POS']} (Ref: {variant['REF']}/Alt: {variant['ALT']})"
        tables.append([title])
        
        header = ["Population", "Group", "Sample Size", "Ref Allele", "Alt Allele", 
                 "Ref HOMOZ", "Alt HOMOZ", "HETEROZ", "HWEP"]
        
        table_data = [header]
        
        for i, sample_data in enumerate(variant.get('SAMPLES', [])):
            if i >= len(population_map):
                break
                
            population = population_map[i]
            fields = sample_data.split(':')
            
            if len(fields) >= 6:
                AN, AC, HWEP, GR, GV, GA = fields
                sample_size = int(AN) // 2 if AN and AN != '.' else 0
                hwep_value = 'N/A' if HWEP == '-1' else f"{float(HWEP):.4f}"
                
                row = [
                    population['name'],
                    population['group'],
                    str(sample_size),
                    variant['REF'],
                    variant['ALT'],
                    GR if GR else '0',
                    GV if GV else '0',
                    GA if GA else '0',
                    hwep_value
                ]
                table_data.append(row)
        
        tables.append(table_data)
    
    return tables

def create_pdf(gene_name, aseq, chrom, genome_start, genome_end, strand, output_filename,dir_o):
    """创建PDF文件"""
    
    # 读取VCF文件
    snps_variants = parse_vcf_file(f"{dir_o}/snps_output.vcf")
    freq_variants = parse_frequency_vcf(f"{dir_o}/snpsfreq_output.vcf")
    
    # 创建PDF文档
    doc = SimpleDocTemplate(output_filename, pagesize=A4, 
                           topMargin=0.5*inch, bottomMargin=0.5*inch,
                           leftMargin=0.5*inch, rightMargin=0.5*inch)
    styles = getSampleStyleSheet()
    
    # 定义网页配色
    primary_dark = hex_to_color(WEB_COLORS['primary_dark'])
    primary_medium = hex_to_color(WEB_COLORS['primary_medium'])
    background_light = hex_to_color(WEB_COLORS['background_light'])
    background_medium = hex_to_color(WEB_COLORS['background_medium'])
    
    # 创建自定义样式 - 使用网页配色
    title_style = ParagraphStyle(
        'CustomTitle',
        parent=styles['Heading1'],
        fontSize=16,
        textColor=primary_dark,
        spaceAfter=12,
        alignment=TA_CENTER
    )
    
    heading_style = ParagraphStyle(
        'CustomHeading',
        parent=styles['Heading2'],
        fontSize=12,
        textColor=primary_medium,
        spaceAfter=6
    )
    
    # 内容列表
    story = []
    
    # 第一部分：显示接收的信息
    story.append(Paragraph("ASO SNP Analysis Report", title_style))
    story.append(Spacer(1, 12))
    
    info_data = [
        ["Parameter", "Value"],
        ["Gene Name", gene_name],
        ["ASO 5'->3'", reverse_complement(aseq)],
        ["Target Sequence 5'->3'", aseq],
        ["Chromosome", chrom],
        ["Genome Start", str(genome_start)],
        ["Genome End", str(genome_end)],
        ["Strand", strand],
        ["Sequence Length", str(len(aseq))],
        ["Number of SNPs", str(len(snps_variants))]
    ]
    
    info_table = Table(info_data, colWidths=[2*inch, 3*inch])
    info_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), primary_dark),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 10),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), background_medium),
        ('GRID', (0, 0), (-1, -1), 1, primary_dark),
        ('FONTSIZE', (0, 1), (-1, -1), 9),
    ]))
    
    story.append(Paragraph("Part 1: Input Parameters", heading_style))
    story.append(info_table)
    story.append(Spacer(1, 20))
    
    # 第二部分：序列图
    story.append(Paragraph("Part 2: Sequence with SNP Annotations", heading_style))
    
    # 创建序列图
    page_width = A4[0] - inch
    seq_drawing = create_sequence_drawing(aseq, chrom, genome_start, genome_end, 
                                         snps_variants, page_width, strand)
    
    # 添加序列图到PDF
    story.append(seq_drawing)
    story.append(Spacer(1, 15))
    
    # 添加图例 - 使用网页配色
    legend_data = [
        ["Symbol", "Description"],
        ["Blue square", "Base nucleotide"],
        ["Red triangle below square", "SNP position marker"],
        ["Red text below marker", "SNP ID"],
        ["Start/End labels", f"{chrom}:{genome_start} and {chrom}:{genome_end}"]
    ]
    
    legend_table = Table(legend_data, colWidths=[1.5*inch, 3*inch])
    legend_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), primary_medium),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('GRID', (0, 0), (-1, -1), 1, primary_dark),
        ('FONTSIZE', (0, 0), (-1, -1), 8),
        ('BACKGROUND', (0, 1), (-1, -1), background_light),
    ]))
    
    story.append(Paragraph("Legend:", styles['Heading3']))
    story.append(legend_table)
    story.append(Spacer(1, 20))
    
    # 第三部分：SNP信息表格
    story.append(Paragraph("Part 3: SNP Variants", heading_style))
    
    if snps_variants:
        snp_data = [["CHROM", "POS", "ID", "REF", "ALT"]]
        
        for variant in snps_variants:
            row = [
                variant['CHROM'],
                str(variant['POS']),
                variant['ID'][:15] + "..." if len(variant['ID']) > 15 else variant['ID'],
                variant['REF'],
                variant['ALT'],
            ]
            snp_data.append(row)
        
        snp_table = Table(snp_data, repeatRows=1)
        snp_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), primary_dark),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 8),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), background_medium),
            ('GRID', (0, 0), (-1, -1), 1, primary_dark),
            ('FONTSIZE', (0, 1), (-1, -1), 7),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, background_light])
        ]))
        
        story.append(snp_table)
    else:
        story.append(Paragraph("No SNP variants found", styles['Normal']))
    
    story.append(Spacer(1, 20))
    
    # 第四部分：频率信息表格
    story.append(Paragraph("Part 4: SNP Frequency Data", heading_style))
    
    if freq_variants:
        frequency_tables = format_frequency_table(freq_variants)
        
        for i, table_data in enumerate(frequency_tables):
            if len(table_data) == 1:
                if i > 0:
                    story.append(Spacer(1, 10))
                story.append(Paragraph(table_data[0], styles['Heading3']))
            else:
                freq_table = Table(table_data, repeatRows=1)
                freq_table.setStyle(TableStyle([
                    ('BACKGROUND', (0, 0), (-1, 0), primary_medium),
                    ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
                    ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                    ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                    ('FONTSIZE', (0, 0), (-1, 0), 6),
                    ('BACKGROUND', (0, 1), (-1, -1), background_light),
                    ('GRID', (0, 0), (-1, -1), 1, primary_dark),
                    ('FONTSIZE', (0, 1), (-1, -1), 5),
                ]))
                story.append(freq_table)
    else:
        story.append(Paragraph("No frequency data found", styles['Normal']))
    
    story.append(Spacer(1, 20))
    
    # 生成PDF
    doc.build(story)
    print(f"PDF generated successfully: {output_filename}")

def main():
    """主函数"""
    if len(sys.argv) != 8:
        print("Usage: python script.py <GENE_NAME> <ASEQ> <CHROM> <GENOME_START> <GENOME_END> <STRAND>")
        print("Example: python script.py DGAT2 ATGCGACT 11 75785889 75785903 +")
        print("Strand: + for positive strand, - for negative strand")
        sys.exit(1)
    
    gene_name = sys.argv[1]
    aseq = sys.argv[2]
    chrom = sys.argv[3]
    genome_start = int(sys.argv[4])
    genome_end = int(sys.argv[5])
    strand = sys.argv[6]
    dir_o = sys.argv[7]
    print(dir_o)
    # 验证链方向参数
    if strand not in ['+', '-']:
        print("Error: Strand must be either '+' or '-'")
        sys.exit(1)
    
    # 处理序列：如果是正链(+)，需要反向互补；如果是负链(-)，保持原样
    original_aseq = aseq
    if strand == '+':
        rna = reverse_complement(aseq)
        print(f"Positive strand detected. Converted sequence: {original_aseq} -> {rna}")
    else:
        rna = aseq
        print(f"Negative strand detected. Using original sequence: {rna}")
    
    output_filename = f"{dir_o}/SNP_Analysis_Report.pdf"
    
    # 检查必要的文件是否存在
    if not os.path.exists(f'{dir_o}/snps_output.vcf'):
        print("Error: snps_output.vcf file not found")
        print("Please make sure snps_output.vcf exists in the current directory")
        sys.exit(1)
    
    if not os.path.exists(f"{dir_o}/snpsfreq_output.vcf"):
        print("Error: snpsfreq_output.vcf file not found")
        print("Please make sure snpsfreq_output.vcf exists in the current directory")
        sys.exit(1)
    
    try:
        create_pdf(gene_name, rna, chrom, genome_start, genome_end, strand, output_filename,dir_o)
    except Exception as e:
        print(f"Error generating PDF: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()