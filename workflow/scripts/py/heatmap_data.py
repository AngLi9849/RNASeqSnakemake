import pandas as pd
import numpy as np
from numpy import log2

scale = pd.read_csv(snakemake.input.size_table,sep='\t',header=0,index_col=0)

control = snakemake.params.control
treat = snakemake.params.treat

samples = snakemake.params.sample_table

control_samples = samples[samples.condition==control].sample_name.tolist()
treat_samples = samples[samples.condition==treat].sample_name.tolist()

control_ls = []
treat_ls = []

bin_n = snakemake.config['heatmap']['bin_number']

main_bin = snakemake.params.main_bin
bef_bin = snakemake.params.bef_bin
plotbef_bin = snakemake.params.plotbef_bin
plotaft_bin = snakemake.params.plotaft_bin
total_bin = main_bin + plotbef_bin + plotaft_bin

bin_size = total_bin/bin_n
start_bin = plotbef_bin + bef_bin
start_pos = int(start_bin/bin_size)
sense_ls = snakemake.input.sense_mx
antisense_ls=snakemake.input.antisense_mx

heat_long=[]
for s in ['sense','antisense'] :
  ls = str(s) + '_ls'
  l = globals()[ls]
  if len(l)>0:
    control_s={}
    treat_s={}
    for i in l:
      mx = pd.read_csv(i,sep='\t',header=0,index_col=0)
      mx = mx.groupby(np.arange(len(mx.columns))//bin_size, axis=1).mean()
      mx = mx*scale.loc[mx.index.name,"scale_factor"]
      if mx.index.name in control_samples:
        control_s[str(mx.index.name)] = mx
      elif mx.index.name in treat_samples:
        treat_s[str(mx.index.name)] = mx
      else:   
        continue
    control_mx = sum(control_s)/len(control_s)
#    control_s['All'] = control_mx
    treat_mx = sum(treat_s)/len(treat_s)
    treat_s['All'] = treat_mx
    for exp in treat_s : 
      if snakemake.wildcards.pair=="paired" :
        control_mx = control_s[samples.sample_name[(samples.condition==control) & (samples.protocol==samples.loc[exp,"protocol"][0]) & (samples.replicate==samples.loc[exp,"replicate"][0])][0]]
      treat_mx = exp_mx
      mean_mx = (control_mx + treat_mx)/2
#      mask_mx = (mean_mx > snakemake.config['heatmap']['min_fc_cov'])*1
      fc_mx = (treat_mx / mean_mx).fillna(1)
      fc_mx = fc_mx-1
#      fc_mx = fc_mx * mask_mx
      mx_bin = len(fc_mx.columns)
      fc_mx.columns = range(0-start_pos,mx_bin - start_pos,1) if (s == "sense") else range(mx_bin-start_pos,0 - start_pos,-1)
      fc_mx.index.name="featureID"
      fc_long = pd.melt(fc_mx.reset_index(), id_vars=fc_mx.index.name, value_vars=fc_mx.columns.tolist())
      fc_long.columns=['featureID','Position','heat']
      fc_long['Sense']=s.capitalize()
      fc_long['Sample'] = 
      mean_mx.columns = range(0-start_pos,mx_bin - start_pos,1) if (s == "sense") else range(mx_bin-start_pos,0 - start_pos,-1)
      mean_mx.index.name="featureID"
      mean_long = pd.melt(mean_mx.reset_index(), id_vars=mean_mx.index.name, value_vars=mean_mx.columns.tolist())
      mean_long.columns=['featureID','Position','Cov']
      fc_long['coverage']=mean_long['Cov']
      heat_long.append(fc_long)
  else:
    continue

heat_long = pd.concat(heat_long)
heat_long.to_csv(snakemake.output.heat_data, sep='\t', header=True, index=False)

