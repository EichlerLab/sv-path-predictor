import pandas as pd
import pandas.errors
import numpy as np
import os
from pybedtools import BedTool
from io import StringIO

configfile: 'config/config.yaml'

thresholds = config.get('thresholds')

CADD_CUTOFF = thresholds.get('cadd', 15)
TADA_CUTOFF = thresholds.get('tada', 0.8)
STRV_CUTOFF = thresholds.get('strv', 0.75)
MANIFEST = config.get('manifest', 'config/manifest.tab')
SV_PATH_SNAKE_DIR = os.path.dirname(workflow.snakefile)

manifest_df = pd.read_csv(MANIFEST, sep='\t', index_col='sample')

def find_vcf(wildcards):
	return manifest_df.at[wildcards.sample, 'vcf_file']

def get_tada(wildcards):
	"""Get non-empty input files for TADA"""
	files = {
		'DEL': checkpoints.split_by_type.get(sample=wildcards.sample).output.del_bed,
		'DUP': checkpoints.split_by_type.get(sample=wildcards.sample).output.dup_bed
	}
	outputs=[]
	for sv_key in files:
		if os.path.getsize(files[sv_key]):
			outputs.append(f'results/{wildcards.sample}/tada/{sv_key}/Annotated_Predicted_TEST.csv')
	return outputs

localrules: all, append_cadd_output, clean_up

rule all:
	input:
		expand(['results/{sample}/{sample}_annotated.bed', 'cleaned_{sample}.done'], sample=manifest_df.index)

rule convert_vcf:
	input:
		vcf = find_vcf
	output:
		tab = temp('results/{sample}/sv_all.tab')
	threads: 1
	resources:
		mem = 8,
		hrs = 24
	envmodules:
		'modules',
		'modules-init',
		'modules-gs/prod',
		'modules-eichler/prod',
		'bcftools/1.9'
	shell:
		'''
		bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%SVLEN\\t%SVTYPE\\n' {input.vcf} > {output.tab}
		'''

rule convert_tab:
	input:
		tab = rules.convert_vcf.output.tab
	output:
		bed = 'results/{sample}/sv-all.bed'
	threads: 1
	resources:
		mem = 8,
		hrs = 24
	run:
		df = pd.read_csv(input.tab, sep='\t', header=None, keep_default_na=False, names=['#CHROM', 'POS', 'REF', 'ALT', 'SVLEN', 'SVTYPE'])
		df['SVLEN'] = df.apply(lambda row: len(row['ALT']) - len(row['REF']), axis=1)
		df['SVLEN'] = np.abs(df['SVLEN'])
		df['END'] = df.apply(lambda row: row['POS']+row['SVLEN'] if row['SVTYPE'] != 'INS' else row['POS']+1, axis=1)
		df = df.loc[df['SVLEN'] >= 50]
		df['ID'] = df['#CHROM']+'-'+df['POS'].astype(str)+'-'+df['SVTYPE']+'-'+df['SVLEN'].astype(str)
		df[['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN']].drop_duplicates().sort_values(['#CHROM', 'POS']).to_csv(output.bed, sep='\t', index=False)


checkpoint split_by_type:
	input:
		bed = rules.convert_tab.output.bed
	output:
		dup_bed = 'results/{sample}/sv_DUP.bed',
		del_bed = 'results/{sample}/sv_DEL.bed',
		ins_bed = 'results/{sample}/sv_INS.bed',
		cadd_input_bed = 'results/{sample}/cadd-sv/id_{sample}.bed',
		temp_cadd_input = 'input/id_{sample}.bed'
	threads: 1
	resources:
		mem = 8,
		hrs = 24
	run:
		df = pd.read_csv(input.bed, sep='\t')
		df.loc[(df['SVTYPE'] == 'DUP')][['#CHROM', 'POS', 'END', 'SVTYPE', 'ID']].to_csv(output.dup_bed, sep='\t', index=False, header=False)
		df.loc[df['SVTYPE'] == 'DEL'][['#CHROM', 'POS', 'END', 'SVTYPE', 'ID']].to_csv(output.del_bed, sep='\t', index=False, header=False)
		df.loc[df['SVTYPE'] == 'INS'][['#CHROM', 'POS', 'END', 'SVTYPE', 'ID']].to_csv(output.ins_bed, sep='\t', index=False, header=False)

		# write input for cadd
		shell('cat {output.dup_bed} {output.del_bed} {output.ins_bed} | sort -k 1,1 -k2,2n > {output.cadd_input_bed} && ln -sf $( readlink -f {output.cadd_input_bed} ) {output.temp_cadd_input}')

# SDIR must be present so it's exposed to the module's snakefile
SDIR = os.path.dirname(config['cadd_config']['snakefile'])
module cadd_sv:
    snakefile: config['cadd_config']['snakefile']
    config: config['cadd_config']

use rule * from cadd_sv as cadd_sv_*

use rule prep_chr1 from cadd_sv as cadd_sv_prep_chr1 with:
	input:
		'results/{set}/cadd-sv/id_{set}.bed'

rule append_cadd_output:
	input:
		bed = 'output/{sample}_score.bed'
	output:
		transformed = 'results/{sample}/cadd-sv/cadd-sv_annotated-{sample}.bed'
	threads: 1
	resources:
		mem = 8,
		hrs = 24
	run:
		# fix header
		with open(input.bed) as file:
			data = file.read()
		data = data[:data.index('\n')].replace(' ','\t') + data[data.index('\n'):]
		stream_data = StringIO(data)

		cadd = pd.read_csv(stream_data,sep='\t',index_col=False)

		cadd.rename(columns={'name': 'ID'}, inplace=True)
		cadd.to_csv(output.transformed, sep='\t', index=False)


rule tada_per_type:
	input:
		bed = 'results/{sample}/sv_{svtype}.bed'
	output:
		annotated = temp('results/{sample}/tada/{svtype}/Annotated_Predicted_TEST.csv'),
	wildcard_constraints:
		svtype='DEL|DUP'
	threads: 1
	resources:
		mem = 8,
		hrs = 24
	envmodules:
		'modules',
		'modules-init',
		'modules-gs/prod',
		'modules-eichler/prod',
		'TADA/1.0.2'
	shell:
		'''
		predict_variants -d -t {wildcards.svtype} -v {input.bed} -o $( dirname {output.annotated} )
		'''

rule tada_combine:
	input:
		annotated = get_tada
	output:
		combined = 'results/{sample}/tada/tada_annotated-{sample}.bed'
	threads: 1
	resources:
		mem = 8,
		hrs = 24
	run:
		df = pd.DataFrame()
		for file in input.annotated:
			if 'DEL' in file:
				sv_df = pd.read_csv(file,sep='\t')
				sv_df['ID'] = sv_df['CHR'] + '-' + sv_df['START'].astype(str) + '-' + 'DEL' + '-' + (sv_df['END'] - sv_df['START']).astype(str)
				df = pd.concat([df,sv_df])
			else:
				sv_df = pd.read_csv(file,sep='\t')
				sv_df['ID'] = sv_df['CHR'] + '-' + sv_df['START'].astype(str) + '-' + 'DUP' + '-' + (
							sv_df['END'] - sv_df['START']).astype(str)
				df = pd.concat([df, sv_df])
		df.to_csv(output.combined, sep='\t', index=False)


rule strvctvre:
	input:
		bed = expand('results/{{sample}}/sv_{svtype}.bed', svtype=['DEL','DUP'])
	output:
		temp_input = temp('results/{sample}/strvctvre/temp_input.bed'),
		annotated = 'results/{sample}/strvctvre/strvctvre_annotated-{sample}.bed'
	threads: 1
	resources:
		mem = 8,
		hrs = 24
	envmodules:
		'modules',
		'modules-init',
		'modules-gs/prod',
		'modules-eichler/prod',
		'bedtools/2.29.2',
		'strvctvre/1.7'
	shell:
		'''
		cat {input.bed} > {output.temp_input}
		StrVCTVRE.py -i temp_input.bed -o {output.annotated} -f bed && \
		sed -ibackup '1 i\CHR\tSTART\tEND\tSVTYPE\tID\tSCORE' {output.annotated}
		cat {output.annotated} | grep -v 'not' > {output.annotated}backup
		mv {output.annotated}backup {output.annotated}
		'''

rule combine_scores:
	input:
		cadd = rules.append_cadd_output.output.transformed,
		tada = rules.tada_combine.output.combined,
		strv = rules.strvctvre.output.annotated
	output:
		combined = 'results/{sample}/{sample}_scored.bed'
	params:
		c_cutoff = CADD_CUTOFF,
		s_cutoff = STRV_CUTOFF,
		t_cutoff = TADA_CUTOFF
	threads: 1
	resources:
		mem = 8,
		hrs = 24
	run:
		try:
			cadd_df = pd.read_csv(input.cadd, sep='\t', usecols=['CADDSV-score','ID'], index_col='ID')
			cadd_df['CADDSV-score'] = cadd_df.apply(lambda row: 1 if (row['CADDSV-score'] >= float(params.c_cutoff)) else 0,axis=1)
			cadd_df = cadd_df[~cadd_df.index.duplicated(keep='first')]
		except pandas.errors.EmptyDataError:
			cadd_df = pd.DataFrame(columns=['CADDSV-score','ID']).set_index('ID')
			print(f'Note: {input.cadd} was empty. Skipping.')

		try:
			tada_df = pd.read_csv(input.tada, sep='\t', usecols=['Pathogenicity Score','Pathogenicity Label', 'ID'], index_col='ID')
			tada_df['Pathogenicity Score'] = tada_df.apply(lambda row: 1 if (row['Pathogenicity Score'] >= float(params.t_cutoff) and row['Pathogenicity Label'] == 1) else 0,axis=1)
		except pandas.errors.EmptyDataError:
			tada_df = pd.DataFrame(columns=['Pathogenicity Score','Pathogenicity Label','ID']).set_index('ID')
			print(f'Note: {input.tada} was empty. Skipping.')

		try:
			strv_df = pd.read_csv(input.strv,sep='\t',usecols=['SCORE', 'ID'],index_col='ID')
			strv_df['SCORE'] = strv_df.apply(lambda row: 1 if (row.SCORE >= float(params.s_cutoff)) else 0,axis=1)
		except pandas.errors.EmptyDataError:
			strv_df = pd.DataFrame(columns=['SCORE', 'ID']).set_index('ID')
			print(f'Note: {input.strv} was empty. Skipping.')

		df = pd.concat([strv_df,tada_df,cadd_df], axis=1)
		df.fillna(0, inplace=True)

		df['AGG_SCORE'] = df.SCORE.astype(int) + df['Pathogenicity Score'].astype(int) + df['CADDSV-score'].astype(int)
		df = df.loc[df['AGG_SCORE'] >= 1]

		df['CHR'] = [x.split('-')[0] for x in df.index]
		df['START'] = [x.split('-')[1] for x in df.index]
		df['END'] = [int(x.split('-')[1]) + int(x.split('-')[3]) if x.split('-')[2] != 'INS' else int(x.split('-')[1]) + 1 for x in df.index]

		df.reset_index(inplace=True)
		df.sort_values(['CHR', 'START']).to_csv(output.combined, sep='\t', index=False, header=False, columns=['CHR','START','END','ID','AGG_SCORE'])


rule annotate_sv:
	input:
		scored = rules.combine_scores.output.combined,
		refseq = config['anno']['hg38']['ncbi_refseq'],
	output:
		scored_and_annotated = 'results/{sample}/{sample}_annotated.bed'
	params:
		window = config.get(config['anno']['window'], 2500)
	threads: 1
	resources:
		mem = 8,
		hrs = 24
	run:
		anno = pd.read_csv(input.refseq,sep='\t',names=['chrom', 'txStart', 'txEnd', 'name', 'score', 'strand'])

		a=BedTool(input.scored)
		b=BedTool.from_dataframe(anno)

		hits=a.window(b, w=params.window)
		target_cols = ['chr','start','end','id','predictor_concordance','refseq-exons']
		try:
			assert hits[0] >= 1
			bed = hits.groupby(g=[1, 2, 3, 4, 5], c=9, o=['collapse']).to_dataframe(names=target_cols)
		except IndexError:
			print(f'Note: No intersection happened for {wildcards.sample}, making empty file.')
			bed = pd.DataFrame(columns=target_cols)

		bed.to_csv(output.scored_and_annotated, sep='\t', index=False)

rule clean_up:
	input:
		requiremnt = rules.annotate_sv.output.scored_and_annotated,
	output:
		done = "cleaned_{sample}.done",
	threads: 1
	resources:
		mem = 1,
		hrs = 24
	shell:
		'''
		rm -rf {wildcards.sample}/ && touch {output.done}
		'''
onsuccess:
	# rm cadd-sv outputs and softlinks
	shell("rm -rf input/ beds/ output/; rm annotations models data; rm *.done")