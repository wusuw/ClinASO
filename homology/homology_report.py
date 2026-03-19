#!/usr/bin/env python3
"""Generate a PDF report for Homology Analysis results - preserves original RNAhybrid formatting."""

import sys
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_CENTER
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer,
                                 Table, TableStyle, HRFlowable)
from reportlab.lib.colors import HexColor

PRIMARY_DARK = HexColor('#1a365d')
PRIMARY_MEDIUM = HexColor('#2b6cb0')
TEXT_DARK = HexColor('#2d3748')
TEXT_LIGHT = HexColor('#718096')
BG_LIGHT = HexColor('#f7fafc')

def get_species_display(raw):
    mapping = {
        'human': 'Human', 'mouse': 'Mouse', 'rat': 'Rat',
        'pig': 'Pig', 'crab_eating_macaque': 'Crab-eating Macaque',
        'rabbit': 'Rabbit', 'guinea_pig': 'Guinea Pig',
    }
    for key, val in mapping.items():
        if key in raw.lower():
            return val
    return raw.split('_')[0].capitalize()

def create_pdf(gene_name, aso_sequence, analysis_type, input_file, output_file):
    with open(input_file, 'r') as f:
        raw_text = f.read()

    # Replace miRNA -> ASO
    raw_text = raw_text.replace('miRNA', 'ASO')

    # Split into sections by "target:" prefix
    sections = []
    current_section = None
    for line in raw_text.split('\n'):
        if line.startswith('target:'):
            if current_section is not None:
                sections.append(current_section)
            current_section = {
                'species': get_species_display(line),
                'lines': [line],
            }
        elif current_section is not None:
            current_section['lines'].append(line)
    if current_section is not None:
        sections.append(current_section)

    doc = SimpleDocTemplate(output_file, pagesize=A4,
                           topMargin=0.6*inch, bottomMargin=0.6*inch,
                           leftMargin=0.6*inch, rightMargin=0.6*inch)

    story = []

    # Header
    story.append(Paragraph("ASO Homology Analysis Report", ParagraphStyle(
        'Title', fontSize=16, textColor=PRIMARY_DARK,
        spaceAfter=4, alignment=TA_CENTER, fontName='Helvetica-Bold')))
    story.append(Paragraph("ClinASO Platform | Yunnan University", ParagraphStyle(
        'Sub', fontSize=9, textColor=TEXT_LIGHT,
        spaceAfter=10, alignment=TA_CENTER)))
    story.append(HRFlowable(width="100%", thickness=1.5, color=PRIMARY_MEDIUM, spaceAfter=10))

    # Parameters
    story.append(Paragraph("Analysis Parameters", ParagraphStyle(
        'H', fontSize=12, textColor=PRIMARY_MEDIUM,
        spaceAfter=8, fontName='Helvetica-Bold')))
    species_label = analysis_type.replace('_', ' ').title() if analysis_type != 'all' else 'All Species (6)'
    param_data = [
        ['Gene Name', gene_name],
        ['ASO Sequence', aso_sequence],
        ['Species', species_label],
    ]
    param_table = Table(param_data, colWidths=[1.5*inch, 4*inch])
    param_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (0, -1), BG_LIGHT),
        ('GRID', (0, 0), (-1, -1), 0.5, colors.Color(0.85, 0.85, 0.85)),
        ('FONTSIZE', (0, 0), (-1, -1), 9),
        ('TOPPADDING', (0, 0), (-1, -1), 4),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
        ('LEFTPADDING', (0, 0), (-1, -1), 8),
    ]))
    story.append(param_table)
    story.append(Spacer(1, 12))

    # Results
    story.append(Paragraph("RNAhybrid Alignment Results", ParagraphStyle(
        'H2', fontSize=12, textColor=PRIMARY_MEDIUM,
        spaceAfter=8, fontName='Helvetica-Bold')))

    mono = ParagraphStyle('Mono', fontSize=9, fontName='Courier',
                           textColor=TEXT_DARK, leading=13, leftIndent=10)

    for sec in sections:
        # Species label
        sp_table = Table(
            [[Paragraph(f"<b>{sec['species']}</b>", ParagraphStyle(
                'SH', fontSize=10, textColor=PRIMARY_DARK, fontName='Helvetica-Bold'))]],
            colWidths=[5.5*inch])
        sp_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (0, 0), BG_LIGHT),
            ('LEFTPADDING', (0, 0), (-1, -1), 10),
            ('TOPPADDING', (0, 0), (-1, -1), 3),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 3),
        ]))
        story.append(sp_table)

        # Raw alignment lines with spaces preserved (use &nbsp; for non-breaking spaces)
        for line in sec['lines']:
            # Escape HTML special chars but preserve spaces
            safe = (line.replace('&', '&amp;')
                       .replace('<', '&lt;')
                       .replace('>', '&gt;')
                       .replace(' ', '&nbsp;'))
            story.append(Paragraph(safe, mono))
        story.append(Spacer(1, 10))

    # Footer
    story.append(HRFlowable(width="100%", thickness=0.5, color=TEXT_LIGHT, spaceBefore=6, spaceAfter=4))
    story.append(Paragraph(
        "<i>ClinASO Platform - School of Life Sciences, Yunnan University</i>",
        ParagraphStyle('Foot', fontSize=7, textColor=TEXT_LIGHT, alignment=TA_CENTER)))

    doc.build(story)
    print(f"PDF generated: {output_file}")


if __name__ == '__main__':
    if len(sys.argv) < 6:
        print("Usage: python homology_report.py <gene_name> <aso_sequence> <analysis_type> <input_file> <output_file>")
        sys.exit(1)
    create_pdf(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
