import pandas as pd
import json

for idx in range(len(snakemake.input)):
	samp_name = snakemake.input[idx].split('/')[-2]
	with open(snakemake.input[idx]) as f:
		data = json.load(f)
	frame = pd.DataFrame.from_dict(data['coverageUniq']['total'],orient='index')
	sample = frame.loc[['intron','exon','intergenic','total']]
	sample.columns = [samp_name]
	if idx==0:
		dataframe=sample
	else:
		dataframe=pd.concat([dataframe,sample],axis=1)
dataframe.to_csv(snakemake.output[0], sep='\t')

