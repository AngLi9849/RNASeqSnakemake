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

for l in [snakemake.input.sense_mx,snakemake.input.antisense_mx]:
  if len(l)>0:
    control_s=[]
    treat_s=[]
    for i in l:
      mx = pd.read_csv(i,sep='\t',header=0,index_col=0)
      mx = mx.groupby(np.arange(len(mx.columns))//bin_size, axis=1).mean()
      mx = mx*scale.loc[mx.index.name,"scale_factor"]
      if mx.index.name in control_samples:
        control_s.append(mx)
      elif mx.index.name in treat_samples:
        treat_s.append(mx)
      else:   
        continue
    control_sum = sum(control_s)
    treat_sum = sum(treat_s)
    control_ls.append(control_sum)
    treat_ls.append(treat_sum) 
  else:
    continue



control_mx = pd.concat(control_ls)
treat_mx = pd.concat(treat_ls)
mean_mx = (control_mx + treat_mx)/2
mask_mx = (mean_mx > snakemake.config['heatmap']['min_fc_cov'])*1
fc_mx = (treat_mx / mean_mx).fillna(1)
fc_mx = fc_mx-1
fc_mx = fc_mx * mask_mx

mx_bin = len(fc_mx.columns)
fc_mx.columns = range(0-start_pos,mx_bin - start_pos)

fc_mx.index.name="featureID"

fc_long = pd.melt(fc_mx.reset_index(), id_vars=fc_mx.index.name, value_vars=fc_mx.columns.tolist())
fc_long.columns=['featureID','Position','heat']

fc_long.to_csv(snakemake.output.heat_data, sep='\t', header=True, index=False)

