import sys
import pandas as pd

lengths = pd.read_csv(snakemake.input.bed,sep='\t',header=None,index_col=7)[4]


bef_ls=[]

for i in snakemake.input.bef:
    bef_ls.append( pd.read_csv(i,sep='\t',header=None,index_col=0, compression='gzip') )

bef = None if bef_ls==[] else pd.concat(bef_ls)

main_ls=[]

for i in snakemake.input.main:
    main_ls.append( pd.read_csv(i,sep='\t',header=None,index_col=0, compression='gzip') )

main = None if main_ls==[] else pd.concat(main_ls)
main_bin = snakemake.params.main_bin

if snakemake.params.main_int:
    main_x_sum=main.sum(axis=1)
else:
    lengths = pd.read_csv(snakemake.input.bed,sep='\t',header=None,index_col=7)[4]
    main_length=main.apply(lambda row:
        lengths.loc[row.name].sum() if (row.name in lengths.index) else main_bin,
        axis=1)
    main_x_sum=main.sum(axis=1)*main_length


#main_mean = main_x_sum[main_x_sum!=0].median() if snakemake.wildcards.mean=="median" else main_x_sum[main_x_sum!=0].mean()


aft_ls=[]

for i in snakemake.input.aft:
    aft_ls.append( pd.read_csv(i,sep='\t',header=None,index_col=0, compression='gzip') )

aft = None if aft_ls==[] else pd.concat(aft_ls)

ls=[bef,main,aft] if snakemake.wildcards.sense=="sense" else [aft,main,bef]

start = 0 - snakemake.params.plotbef_bin - snakemake.params.bef_bin

end = 0 + snakemake.params.main_bin + snakemake.params.plotaft_bin - snakemake.params.bef_bin

mx=pd.concat(ls,axis=1)

mx.columns = range(start,end)

mx.index.name = snakemake.wildcards.sample

mx.to_csv(snakemake.output.mx, sep='\t', header=True, index=True, compression='gzip')
#mx.to_csv(snakemake.output.sum_mx, sep='\t', header=True, index=True, compression='gzip')

#norm_mx = mx.div(main_x_sum,axis=0)*(main_mean)
#norm_mx = norm_mx.fillna(0)

#norm_mx.to_csv(snakemake.output.norm_mx, sep='\t', header=True, index=True, compression='gzip')

