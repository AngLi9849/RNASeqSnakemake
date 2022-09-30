import pandas as pd
import numpy as np

scale = pd.read_csv(snakemake.input.size_table,sep='\t',compression='gzip',header=0,index_col=0)

control = snakemake.params.control
treat = snakemake.params.treat

samples = snakemake.params.sample_table

control_samples = samples[samples.condition==control].sample_name.tolist()
treat_samples = samples[samples.condition==treat].sample_name.tolist()

control_ls = []
treat_ls = []

main_bin = snakemake.params.main_bin
bef_bin = snakemake.params.bef_bin
plotbef_bin = snakemake.params.plotbef_bin
plotaft_bin = snakmake.params.plotaft_bin
total_bin = main_bin + plotbef_bin + plotaft_bin

bin_size = total_bin//100
start_bin = plotbef_bin + bef_bin
start_pos = int(start_bin//bin_size)

for l in [snakemake.input.sense_mx,snakemake.input.antisense_mx]:
  if len(l)>0:
    for i in snakemake.input.sense_mx:
      mx = pd.read_csv(i,sep='\t',header=0,index_col=0)
      mx = mx.groupby(np.arange(len(mx.columns))//bin_size, axis=1).mean()
      mx = mx*scale.loc[mx.index.name].scale_factor
      if mx.index.name in control_samples:
        control_ls.append(mx)
      elif mx.index.name in treat_sampels:
        treat_ls.append(mx)
      else:
        continue 
   else:
     continue

control_mx = sum(control_ls)
treat_mx = sum(treat_ls)
mean_mx = (control_mx + treat_mx)//2
fc_mx = (treat_mx // mean_mx).fillna(1)

