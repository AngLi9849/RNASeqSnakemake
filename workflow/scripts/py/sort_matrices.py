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

start = 0 - snakemake.params.bef_bin

end = 0 + snakemake.params.main_bin + snakemake.params.aft_bin

mx=pd.concat(ls,axis=1)

mx.columns = range(start,end)

mx.index.name = "id"

mx.to_csv(snakemake.output.matrix, sep='\t', header=True, index=True, compression='gzip')



