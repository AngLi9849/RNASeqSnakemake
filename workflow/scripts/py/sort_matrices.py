import sys
import pandas as pd

bef_ls=[]

for i in snakemake.input.bef:
    bef_ls.append( pd.read_csv(i,sep='\t',header=None,index_col=0, compression='gzip') )

bef = None if bef_ls==[] else pd.concat(bef_ls)

main_ls=[]

for i in snakemake.input.main:
    main_ls.append( pd.read_csv(i,sep='\t',header=None,index_col=0, compression='gzip') )

main = None if main_ls==[] else pd.concat(main_ls)

aft_ls=[]

for i in snakemake.input.aft:
    aft_ls.append( pd.read_csv(i,sep='\t',header=None,index_col=0, compression='gzip') )

aft = None if aft_ls==[] else pd.concat(aft_ls)

ls=[bef,main,aft]

start = 0 - snakemake.params.plotbef_bin - snakemake.params.bef_bin

end = 0 + snakemake.params.main_bin + snakemake.params.plotaft_bin - snakemake.params.bef_bin

mx=pd.concat(ls,axis=1)

mx.columns = range(start,end)

mx.index.name = "id"

mx.to_csv(snakemake.output.sum_mx, sep='\t', header=True, index=True, compression='gzip')

norm_mx = mx.div(mx.sum(axis=1),axis=0)*(mx.sum().sum()/len(mx))

norm_mx.to_csv(snakemake.output.norm_mx, sep='\t', header=True, index=True, compression='gzip')

