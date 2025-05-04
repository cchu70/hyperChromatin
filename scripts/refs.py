import seaborn as sns

# Seaborn colorblind palette
cell_types = ['TAC-1', 'TAC-2', 'IRS', 'Medulla', 'Hair Shaft-cuticle.cortex']
pal = sns.color_palette('colorblind', len(cell_types)).as_hex()

celltype_colors = dict(zip(cell_types, pal))
