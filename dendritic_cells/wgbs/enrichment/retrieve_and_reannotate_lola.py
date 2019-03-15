"""Create annotation files for LOLA codex and encode experiments

The idea is to leave the index.txt files untouched and exchange the annotations
after LOLA has run in order to comply with my downstream requirements

This way, LOLA does not need to be rerun when annotations change

Output annotation files comply to my standard annotation format, ie. they
have the columns

- db: database, e.g. tf-chipseq
- sub_db: biologically meaningful subset, e.g. encode, codex
- short_exp_name: unique ID; short descriptive name for the region set used in plots,
                  Merging is done based on short_exp_name
- long_exp_name: unique ID; long descriptive name, e.g. for hover text
- biol_entity: used when pooling experiments and when checking whether there is only one experiment per biological entity (e.g. per TF)

Note that the subdatabase (sub_db) annotation used to distinguish biologically
meaningful sets is different from the collections, used to create manageable
chunks of test regions with respect to memory consumption
"""

import subprocess
import pandas as pd

from dendritic_cells.base_paths import enrichment_database_dir

lola_chipseq_dirname = 'lola_chipseq_2018-04-12'
lola_chipse_dir = enrichment_database_dir + '/lola_chipseq_dirname'

subprocess.run(f"""
cd {enrichment_database_dir}
mkdir {lola_chipseq_dirname}
cd {lola_chipseq_dirname}
wget http://big.databio.org/regiondb/LOLACore_180412.tgz
tar -xzf LOLACore_180412.tgz
mv nm/t1/resources/regions/LOLACore/mm10 .
rm -r nm
""")





# Resources
# =============================================================================

# Codex cell super type mapping
# -----------------------------------------------------------------------------
# Note that there are many name and spelling variants of the same super type

# codex_anno.super_cell_type.unique()
orig_names = ['Mast', 'Multipotent myeloid progenitor',
              'Haematopoietic progenitor', 'Leukemogenic',
              'Mouse ErythroLeukaemic', 'Macrophages', 'Adipocytes', 'T-Cells',
              'B-Cells', 'Pre-adipocyte', 'Erythroid',
              'Megakaryocyte Progenitors', 'AML pre-leukaemic',
              'Embryonic stem cells', 'Embryonic Stem Cell',
              'Embryonic Stem Cells', 'Fibroblast', 'Osteoblast precursor',
              'Osteoblast', 'Myeloid progenitor cells', 'Embryonic Stem cell',
              'neural progenitor cells', 'Mesoderm', 'Plasmablasts',
              'Embryoid bodies', 'Megakaryocyte', 'Leukaemia',
              'Hematopoietic Progenitor Cells', 'Pancreas Cells',
              'Neural Tube Cells', 'Embryonic stem cell', 'NK cells',
              'Macrophage', 'Lymphoid cells', 'Dendritic',
              'Primary murine bone marrow macrophage cells',
              'Mouse Embryonic fibroblasts', 'Hematopoietic Stem Cells',
              'Erythroid Progenitors', 'Erythroid progenitor',
              'Myeloid Progenitors', 'Cardiac muscle', 'Neuron',
              'Embyonic stem cell', 'Haemogenic endothelium',
              'Haematopoietic progenitors', 'Megakaryoblastic cells',
              'Haematopoietic precursor and progenitor cells',
              'Haematopoietic stem and progenitor cells', 'Myeloblastic',
              'Neural progenitor cells', 'embryonic stem cell', 'Thymus',
              'Mouse ErythroLeukaemic ', 'striatal cells', 'embryonic stem cells',
              'Embryonic fibroblast cells', 'Embryonic Stem cells',
              'Motor neuron progenitors', 'Embryonic Trunk']

new_names = ['Mast', 'MMP',
             'HPC', 'Leukemogenic',
             'MEL', 'Macrophage', 'Adipocyte', 'T-cell',
             'B-cell', 'Pre-adipocyte', 'Erythroid',
             'MkPs', 'AML pre-leukaemic',
             'ESC', 'ESC',
             'ESC', 'Fibroblast', 'Osteoblast precursor',
             'Osteoblast', 'MeyloP', 'ESC',
             'neural progenitor cells', 'Mesoderm', 'Plasmablasts',
             'Embryoid bodies', 'Megakaryocyte', 'Leukaemia',
             'HPC', 'Pancreas Cells',
             'Neural Tube Cells', 'ESC', 'NK cells',
             'Macrophage', 'Lymphoid cells', 'Dendritic',
             'BMDM',
             'MEF', 'HSC',
             'EryP', 'EryP',
             'MeyloP', 'Cardiac muscle', 'Neuron',
             'ESC', 'Haemogenic endothelium',
             'HPC', 'Megakaryoblastic cells',
             'HPC',
             'HSPC', 'Myeloblastic',
             'Neural progenitor cells', 'ESC', 'Thymus',
             'MEL', 'striatal cells', 'ESC',
             'Embryonic fibroblast cells', 'ESC',
             'Motor neuron progenitors', 'Embryonic Trunk']

codex_cell_super_type_mapping = dict(zip(orig_names, new_names))

# Curate encode annotations
# =============================================================================

print('Working on encode data')

encode_index_df = pd.read_csv(snakemake.input.mouse_encode_orig_index_file[0],
                              sep = '\t')

encode_anno_df = encode_index_df.copy()

# remove antibody annotations from antibody name
# so that name is only TF
# (all experiments for the same TF have the same antibody annotation, if present)
encode_anno_df.antibody = (encode_anno_df.antibody.str.replace('_\(.*\)$', '')
                           # manually curate some names
                           # GATA-1 and GATA1
                           .replace({'GATA-1': 'GATA1',
                                     'Pol2(phosphoS2)': 'Pol2_phosphoS2'}))

# running experiment ID, will be used for merging
encode_anno_df['exp_id'] = 'MEC' + pd.Series(range(1, encode_anno_df.shape[0] + 1)).astype(str)

# Unique experiment name comprising db name, antibody, cell type, id
# For display in plots
encode_anno_df['exp_name_short'] = (encode_anno_df.antibody
                                    + '_' + encode_anno_df.cellType
                                    + '_' + encode_anno_df['exp_id'])

# Long description, e.g. for hover text
encode_anno_df['exp_name_long'] = (encode_anno_df['exp_name_short']
                                   + ' treatment: ' + encode_anno_df.treatment
                                   + ' description: ' + encode_anno_df.description)

# Add missing fields required be the enrichment result format of the workflow
encode_anno_df['db'] = 'dna_binding_chipseq'
encode_anno_df['sub_db'] = 'encode'
encode_anno_df['biol_entity'] = encode_anno_df['antibody']

# reorder columns - not enforced by file format, but nicer
encode_anno_df = encode_anno_df[[
    'db', 'sub_db',
    'biol_entity', 'exp_name_short',
    'antibody', 'cellType', 'treatment', 'description',
    'exp_id', 'filename', 'exp_name_long',
]]

encode_anno_df.rename(index={'cellType': 'cell_type'}, inplace=True)

# lola discards unknown columns, so have to save expid in standard field
encode_index_df['description'] = encode_anno_df['exp_id']

encode_index_df.to_csv(snakemake.output.mouse_encode_curated_index_file[0],
                       sep='\t', header=True, index=False)
encode_anno_df.to_csv(snakemake.output.mouse_encode_annotations[0],
                      sep='\t', header=True, index=False)


# Curate codex annotations
# =============================================================================

print('Working on codex data')

codex_index = pd.read_csv(snakemake.input.codex_orig_index_file[0], sep='\t')

codex_anno = codex_index.copy()

# codex_anno has a super and sub cell type, named cellType and cellTypeSubtype
# cellTypeSubtype is actually composed of sample type marker (e.g. '^[CL]'
# for cell line and a cell subtype
# rename
codex_anno = codex_anno.rename(columns={'cellType': 'super_cell_type',
                                        'cellTypeSubtype': 'mixed_type_column'})

codex_anno[['sample_type', 'cell_subtype']] = (codex_anno
    .mixed_type_column.str.extract(
        '\[(?P<sample_type>.*)\] (?P<cell_subtype>.*)$', expand=True))

# cell super and subtype names are very long. Get under some control
# manual curation later on will still be necessary
codex_anno.super_cell_type.replace(codex_cell_super_type_mapping, inplace=True)

# create celltype column by combining super and subtype
codex_anno['cell_type'] = codex_anno.super_cell_type + '_' + codex_anno.cell_subtype

# running experiment ID
codex_anno['exp_id'] = 'CDX' + pd.Series(range(1, codex_anno.shape[0] + 1)).astype(str)

# Unique experiment name comprising db name, antibody, cell type, id
codex_anno['exp_name_short'] = (codex_anno.antibody
                                + '_' + codex_anno.cell_type
                                + '_' + codex_anno['exp_id'])

# Long description, e.g. for hover text
codex_anno['exp_name_long'] = (codex_anno['exp_name_short']
                               + ' description: ' + codex_anno.description)

# Add missing fields required be the enrichment result format of the workflow
codex_anno['db'] = 'dna_binding_chipseq'
codex_anno['sub_db'] = 'codex'
codex_anno['biol_entity'] = codex_anno['antibody']

# must assign these columns to index before dropping them from anno df
# cellType overwritten since there are two cell type columns originally
codex_index['cellType'] = codex_anno['cell_type']
# lola discards unknown columns, so have to save expid in standard field
codex_index['description'] = codex_anno['exp_id']

# reorder and select columns - not enforced by file format, but nicer
codex_anno = codex_anno[[
    'db', 'sub_db',
    'biol_entity', 'exp_name_short',
    'antibody', 'description',
    'exp_id', 'filename', 'exp_name_long',
]]

codex_index.to_csv(snakemake.output.codex_curated_index_file[0],
                   sep='\t', header=True, index=False)
codex_anno.to_csv(snakemake.output.codex_annotations[0],
                  sep='\t', header=True, index=False)


# common annotation
common_annotation = pd.concat([encode_anno_df, codex_anno], axis=0)

common_annotation.to_csv(snakemake.output.common_annotation,
                         sep='\t', header=True, index=False)


