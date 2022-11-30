# All Python variables and functions to be used in the whole pipeline

import glob

import pandas as pd
import hashlib as hashlib
import re
from os.path import dirname
from os.path import basename
from snakemake.remote import FTP
from snakemake.utils import validate


ftp = FTP.RemoteProvider()

# Generic Utility and helper functions

def capitalize(line):
    return ' '.join(s[:1].upper() + s[1:] for s in line.split(' '))

#validate(config, schema="../schemas/config.schema.yaml")

# Read references_info config table into pandas dataframe
references = (pd.read_csv(config["references_info"], sep="\t", dtype={"species": str, "ensembl_build": str },comment="#"))
references.columns=references.columns.str.strip()
references = references.applymap(lambda x: x.strip() if isinstance(x, str) else x)
references = references.mask(references=='')
references = (references.set_index(["species"], drop=False).sort_index())

# Helper functions for references 

def get_ref_source(species):
    if not pd.isna(references.loc[species,"genome_dir"]):
        return str("local_" + species)
    else: 
        if (species=="homo_sapiens") & (references.loc[species,"ensembl_build"]=="GRCh38") & (config["MANE_gene_range"]):
            return "MANE_{S}.{B}.{R}".format(S=species,B=references.loc[species,"ensembl_build"],R=references.loc[species,"ensembl_release"])
        else:
            return "ensembl_{S}.{B}.{R}".format(S=species,B=references.loc[species,"ensembl_build"],R=references.loc[species,"ensembl_release"])


def get_genome(species):
    if pd.notna(references.loc[species,"genome_dir"]):
        return "resources/genomes/local_{S}_genome.fasta".format(S=species)
    else :
        return "resources/genomes/ensembl_{S}.{B}.{R}_genome.fasta".format(S=species,B=references.loc[species,"ensembl_build"],R=references.loc[species,"ensembl_release"])


def get_annotation(species):
    if not pd.isna(references.loc[species,"annotation_dir"]):
        return "resources/annotations/local_{S}_genome.gtf".format(S=species)
    else :
        return "resources/annotations/ensembl_{S}.{B}.{R}_genome.gtf".format(S=species,B=references.loc[species,"ensembl_build"],R=references.loc[species,"ensembl_release"])

# Read feature_parameters config table into pandas dataframe
features = (pd.read_csv(config["feature_parameters"], sep="\t", dtype={"feature_name": str, "feature": str, "subfeat": str , "annotation_bed6": str}, comment="#"))
features.columns = features.columns.str.strip()
features = features.applymap(lambda x: x.strip() if isinstance(x, str) else x)
features = features.mask( features == '')
features = (features.set_index(["feature_name"], drop=False).sort_index())
feat_prefix = features.iloc[:, features.columns.get_loc("region"): (features.columns.get_loc("len_aft")+1)]
feat_plot = features.iloc[:,features.columns.get_loc("region"): features.columns.get_loc("s2b_min")]
features["prefix"] = feat_prefix.apply(lambda row: 
    ''.join([str(a)+ "_" + str(b) + "." for a,b in zip(row.index.tolist(), row.tolist())]) 
    , axis=1) 
features["plot_prefix"] = feat_plot.apply(lambda row:
    ''.join([str(a)+ "_" + str(b) + "." for a,b in zip(row.index.tolist(), row.tolist())])
    , axis=1)

features['prefix_md5'] = [hashlib.md5(i.encode('utf-8')).hexdigest() for i in features.prefix.tolist()]
features['plot_md5'] = [hashlib.md5(i.encode('utf-8')).hexdigest() for i in features.plot_prefix.tolist()]

# Helper functions for features
def get_feature_type(feature):
    if feature in features["feature_name"].tolist():
        return "custom-" + str(features.loc[feature,"prefix_md5"])
    else:
        return "gtf"


def get_root_feature(feature):
    base=features.loc[feature].squeeze(axis=0)["feature"]
    if not pd.isna(features.loc[feature].squeeze(axis=0)["annotation_bed"]):
        return str(feature)
    else:
        if base in features["feature_name"].tolist():
            return get_root_feature(base)
        else:
            return str(features.loc[feature].squeeze(axis=0)["feature"])

def get_feature_validity(feature):
    feat=features.loc[feature].squeeze(axis=0)["feature"]
    if not pd.isna(features.loc[feature].squeeze(axis=0)["annotation_bed"]):
        return "provided"
    else:
        if feat in features["feature_name"].tolist():
            return get_feature_validity(feat)
        else:
            if not config["features"]["validate_features"]:
                return "annotated"
            else:
                return "validated"

features["sense_dir"]=features.apply(lambda row:
    -1 if str(row.sense)[0]=="-" else 1
    , axis=1)

def get_part_bin_number(feature,part):
    row=features.loc[feature]
    if part == "main" :
        return int(row.bin_n)
    else :
        opposite="plotbef" if part=="plotaft" else "plotaft"
        main_bin=int(row.bin_n)
        main_is_int=not ( ("x" in str(row.len_bef)) or ("x" in str(row.len_aft)) or str(row.section)=="body" )
        r = 1 if pd.isna(row.befaftr) else (float(row.befaftr)/(1-float(row.befaftr)))
        if ("x" in str(row.plotbef)) and ("x" in str(row.plotaft)) :
            bef = float(re.search(r".*(?<!x)",str(row.plotbef)).group())
            aft = float(re.search(r".*(?<!x)",str(row.plotaft)).group())
            part_bin = bef if part=="plotbef" else aft
            part_bin = part_num * main_bin * r / (bef + aft)
        elif "x" in str(row[part]) :
            part_r=float(re.search(r".*(?<!x)",str(row[part])).group())
            part_bin= (part_r * main_bin * r)
        elif "x" in str(row[opposite]) :
            opp_r=float(re.search(r".*(?<!x)",str(row[opposite])).group())
            if main_is_int:
                main_len=int(row.len_bef) + int(row.len_aft)
                part_bin= (main_bin * r * int(row[part]) / main_len)
            else :                
                part_bin= opp_r * main_bin * r
        elif main_is_int :
            main_len= int(row.len_bef) + int(row.len_aft)
            part_bin = (main_bin * r * int(row[part]) / main_len)
        else :
            part_bin = ( main_bin * r * int(row[part]) / ( int(row[part]) + int(row[opposite]) ) )
    return int(part_bin)

def get_bef_bin_number(feature):
    row=features.loc[feature]
    if row.section == "body" :
        return 0
    else :
        opposite="len_aft"
        main_bin=int(row.bin_n)
        main_is_int=not ( ("x" in str(row.len_bef)) or ("x" in str(row.len_aft)) or str(row.section)=="body" )
        if ("x" in str(row.len_bef)) and ("x" in str(row.len_aft)) :
            bef = float(re.search(r".*(?<!x)",str(row.len_bef)).group())
            aft = float(re.search(r".*(?<!x)",str(row.len_aft)).group())
            part_num = bef if row.sense_dir==1 else aft
            part_bin = part_num * main_bin / (bef + aft)
        elif main_is_int :
            main_len= int(row.len_bef) + int(row.len_aft)
            bef = int(row.len_bef) if row.sense_dir==1 else int(row.len_aft)
            part_bin = (main_bin * bef / main_len)
    return int(part_bin)

features["valid"] = features.apply(lambda row:
    get_feature_validity(row.feature_name)
    , axis=1)
features["type"] = features.apply(lambda row:
    get_feature_type(row.feature_name)
    , axis=1)
features["plotbef_bin"]=features.apply(lambda row: 
    get_part_bin_number(row.feature_name,"plotbef")
    , axis=1)
features["plotaft_bin"]=features.apply(lambda row: 
    get_part_bin_number(row.feature_name,"plotaft")
    , axis=1)
features["bef_bin"]=features.apply(lambda row:
    get_bef_bin_number(row.feature_name)
    , axis=1)

features["antisense"]=features.apply(lambda row:
    not (pd.isna(row.sense) or len(str(row.sense))==1) 
    , axis=1)
features["is_main_int"]=features.apply(lambda row:
    not ( ("x" in str(row.len_bef)) or ("x" in str(row.len_aft)) or str(row.section)=="body" )
    , axis=1)


# Read protocol configurations into pandas data.frame
protocols = (pd.read_csv(config["protocols"], sep="\t", dtype={"protocol": str,"dedup": str, "demulti": str}, comment="#"))
protocols.columns=protocols.columns.str.strip()
protocols = protocols.applymap(lambda x: x.strip() if isinstance(x, str) else x)
protocols = protocols.mask(protocols == '')
protocols = (protocols.set_index(["protocol"], drop=False).sort_index())

# Read gene_sets
gene_sets = (pd.read_csv(config["gene_sets"], sep="\t", dtype={"set_name": str, "genes": str}, comment="#"))
gene_sets.columns=gene_sets.columns.str.strip()
gene_sets = gene_sets.applymap(lambda x: x.strip() if isinstance(x, str) else x)
gene_sets = gene_sets.mask(gene_sets == '')
gene_sets = (gene_sets.set_index(["set_name"], drop=False).sort_index())

# Read read processing configs
reads = (pd.read_csv(config["read_alignments"], sep="\t", dtype={ "name": str}, comment="#"))
reads.columns=reads.columns.str.strip()
reads = reads.applymap(lambda x: x.strip() if isinstance(x, str) else x)
reads = reads.mask(reads == '')
reads = (reads.set_index(["name"], drop=False).sort_index())
reads["single_nuc"] = reads.apply(lambda row:\
    row.sn_pos > 0 or row.sn_pos < 0,
    axis=1
)

   
# Read experimentss config table into pandas dataframe
experiments = (pd.read_csv(config["experiments"], sep="\t", dtype={"protocol": str, "sample_lineage" : str}, comment="#"))
experiments.columns=experiments.columns.str.strip()
experiments = experiments.applymap(lambda x: x.strip() if isinstance(x, str) else x)
experiments = experiments.mask(experiments == '')
experiments["experiment"] = experiments.apply( lambda row: \
    ( str( row.treatment ) + "_vs_" + str( row.control ) + "_" + str( row.protocol ) ), axis = 1
)
experiments = (experiments.set_index(["experiment"], drop=False).sort_index())

experiments["trs_val"]=experiments.apply(
    lambda row: row.trs_val & protocols.loc[row.protocol,"trs_val"], axis=1)
experiments["splice"]=experiments.apply(
    lambda row: protocols.loc[row.protocol,"splice"], axis=1)
experiments["norm_feat"]=experiments.apply(
    lambda row: protocols.loc[row.protocol,"norm_feat"], axis=1)
experiments["demultimap"]=experiments.apply(
    lambda row: protocols.loc[row.protocol,"demulti"], axis=1)
experiments["deduplicate"]=experiments.apply(
    lambda row: protocols.loc[row.protocol,"dedup"], axis=1)
experiments["sep_spl"]=experiments.apply(
    lambda row: protocols.loc[row.protocol,"sep_spl"], axis=1)

experiments["GOI"]=experiments.apply(lambda row: \
     gene_sets.loc[list(set( [x.strip(' ') for x in str(row.gene_sets).split(",")] ).intersection(gene_sets.set_name.tolist())),"genes"].str.cat() if not pd.isna(row.gene_sets) else "" , axis=1)

experiments["norm_group"]=experiments.apply(
    lambda row: sorted(
        pd.unique(
            [row.control,row.treatment] + \
            experiments[experiments.sample_species == row.sample_species][experiments.control==row.control][ experiments.protocol==row.protocol]["treatment"].tolist()
        )
    ), axis=1)

experiments["group_name"] = experiments.apply(
    lambda row: "_and_".join(row.norm_group) + "_" + str(row.protocol), 
    axis = 1 )

# Helper Functions for experiments

def get_source(experiment):
    if not pd.isna(experiments.loc[experiment,"spikein_species"]):
        return "combined_{sample_species}_and_{spikein_species}".format(
            sample_species=get_ref_source(experiments.loc[experiment,"sample_species"]),
            spikein_species=get_ref_source(experiments.loc[experiment,"spikein_species"])
        )
    else :
        return get_ref_source(experiments.loc[experiment,"sample_species"])

def get_sample_source(experiment):
    return get_ref_source(experiments.loc[experiment,"sample_species"])

# Assign reference to each experiment
experiments["reference"] = experiments.apply( lambda row: get_source(row.experiment), axis = 1 )
experiments["normaliser"] = experiments.apply( lambda row: "spikein" if row.spikein else "internal" , axis = 1 )

# Read samples config table into pandas dataframe
samples = (pd.read_csv(config["samples"], sep="\t", dtype={"protocol": str, "replicate": str, "unit_name": str}, comment="#"))
samples["sample_name"]=samples.apply(lambda row: str(row.condition) + "_" + str(row.protocol) + "_Replicate_" + str(row.replicate), axis=1)
samples=samples.set_index(["sample_name","unit_name"], drop=False).sort_index()
samples = samples.mask(samples == '')

# Helper functions for samples

def is_paired_end(sample):
    sample_units = samples.loc[sample]
    fq2_null = sample_units["fq2"].isnull()
    sra_null = sample_units["sra"].isnull()
    paired = ~fq2_null | ~sra_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "invalid read files for sample {}, must be all paired end or all single end".format(
        sample
    )
    return all_paired

def get_experiment_samples(experiment):
    exp = experiments.loc[experiment].squeeze()
    cond = [exp.control,exp.treatment]
    sample = samples[samples.protocol == exp.protocol][samples.condition.isin(cond)]
    return sample

def get_experiment_qc_samples(experiment):
    exp = experiments.loc[experiment].squeeze()
    cond = [exp.control,exp.treatment]
    sample = samples[samples.protocol == exp.protocol][samples.condition.isin(cond)]
    sample_fq1 = sample[['sample_name','unit_name']][pd.notna(sample.fq1)]
    sample_fq1['fq']='fq1'
    sample_fq2 = sample[['sample_name','unit_name']][pd.notna(sample.fq2)]
    sample_fq2['fq']='fq2'
    sample=pd.concat([sample_fq1,sample_fq2])
    return sample

def get_experiment_controls(experiment):
    exp = experiments.loc[experiment].squeeze()
    cond = [exp.control]
    sample = samples[samples.protocol == exp.protocol][samples.condition.isin(cond)]
    return sample

def get_experiment_treatments(experiment):
    exp = experiments.loc[experiment].squeeze()
    cond = [exp.treatment]
    sample = samples[samples.protocol == exp.protocol][samples.condition.isin(cond)]
    return sample

def get_experiment_readlen(experiment):
    sample = get_experiment_samples(experiment)
    readlen=int(sample.readlen.mean())
    return readlen

def is_experiment_readpaired(experiment):
    sample = get_experiment_samples(experiment)
    paired=sample.apply(lambda row: is_paired_end(row.sample_name), axis=1).all()
    return paired

def get_norm_group_samples(norm_group):
    exp = experiments[experiments.group_name==norm_group]
    cond = exp.norm_group[0]
    protocol = exp.protocol[0]
    sample = samples[samples.protocol == protocol][samples.condition.isin(cond)]
    return sample

experiments["paired"]=experiments.apply( lambda row: "paired" if row.pairRep else "unpaired" , axis = 1 )
experiments["diff_lineage"]=experiments.sample_lineage if config['features']['validate_features'] else "genome"
experiments["demulti"]=experiments.apply(lambda row: "Demultimapped" if row.demultimap else "", axis=1)
experiments["dedup"]=experiments.apply(lambda row: "Deduplicated" if row.deduplicate else "", axis=1)


# Create a summary data frame to assign samples to lineages
lineage=[]

for i in range(0,len(experiments)):
    samps = get_experiment_samples(experiments.experiment[i])[["sample_name","unit_name"]]
#    samps["experiment"] = experiments.experiment[i]
    samps["reference"] = get_source(experiments.experiment[i])
    samps["lineage"] = experiments.sample_lineage[i]
    samps["species"] = experiments.sample_species[i]
    samps["group"] = experiments.group_name[i]
    samps["trs_val"] = experiments.trs_val[i]
    samps["splice"] = experiments.splice[i]
    lineage.append(samps)

lineage = pd.concat(lineage).drop_duplicates()
lineage=lineage.set_index(["species","lineage"], drop=False).sort_index()

def get_reference_species(reference):
    species=lineage[lineage.reference==reference].species[0]
    return species

# Define keywords for obtaining sample-specific bigwigs
bw_samples=[]

for i in range(0,len(experiments)):
    samps = get_experiment_samples(experiments.experiment[i])[["sample_name","unit_name","strandedness"]]
    samps["experiment"]=experiments.experiment[i]
    samps["reference"] = get_source(experiments.experiment[i])
    samps["lineage"] = experiments.sample_lineage[i]
    samps["group"] = experiments.group_name[i]
    samps["spikein"]=experiments.normaliser[i]
    samps["deduplicate"]=experiments.deduplicate[i]
    samps["demultimap"]=experiments.demultimap[i]
    samps["sep_spl"]=experiments.sep_spl[i]
    samps["norm_feat"]=experiments.norm_feat[i]
    samps["pairRep"]="paired" if experiments.pairRep[i] else "unpaired"
    samps["stranded"]=samps.apply(lambda row: "unstranded" if (row.strandedness=="no" or pd.isna(row.strandedness))  else "stranded", axis=1)
    unstrand=samps[samps.stranded=="unstranded"]
    unstrand["strand"]="unstranded"
    fwd=samps[samps.stranded=="stranded"]
    fwd["strand"]="fwd"
    rev=fwd.copy(deep=True)
    rev["strand"]="rev"
    bw_samples.append(unstrand)
    bw_samples.append(fwd)
    bw_samples.append(rev)

bw_samples = pd.concat(bw_samples).drop_duplicates()
bw_samples=bw_samples.set_index(["sample_name","unit_name"], drop=False).sort_index()

bw_samples["demulti"]=bw_samples.apply(
    lambda row: ['','Demultimapped'] if row.demultimap=="BOTH" else ['Demultimapped'] if row.demultimap=="TRUE" else [''], axis=1)
bw_samples["dedup"]=bw_samples.apply(
    lambda row: ['','Deduplicated'] if row.deduplicate=="BOTH" else ['Deduplicated'] if row.deduplicate=="TRUE" else [''], axis=1)
bw_samples["splice_prefix"]=bw_samples.apply(
    lambda row: ['All','Spliced','Unspliced'] if row.sep_spl else ['All'], axis=1)
bw_samples["demulti"]=bw_samples["demulti"].astype('object')
bw_samples["dedup"]=bw_samples["dedup"].astype('object')
bw_samples["splice_prefix"]=bw_samples["splice_prefix"].astype('object')

results=[]

for i in range(0,len(experiments)):
    samps = get_experiment_samples(experiments.experiment[i])[["sample_name","unit_name","strandedness"]]
    samps["experiment"] = experiments.experiment[i]
    samps["reference"] = get_source(experiments.experiment[i])
    samps["lineage"] = experiments.sample_lineage[i]
    samps["diff_lineage"] = experiments.diff_lineage[i]
    samps["group"] = experiments.group_name[i]
    samps["spikein"] = experiments.normaliser[i]
    samps["deduplicate"] = experiments.deduplicate[i]
    samps["demultimap"] = experiments.demultimap[i]
    samps["sep_spl"] = experiments.sep_spl[i]
    samps["norm_feat"] = experiments.norm_feat[i]
    samps["pairRep"] = "paired" if experiments.pairRep[i] else "unpaired"
    samps["stranded"] = samps.apply(lambda row: "unstranded" if (row.strandedness=="no" or pd.isna(row.strandedness))  else "stranded", axis=1)
    samps["demulti"] = samps.apply(lambda row: "Demultimapped" if (row.demultimap=="TRUE")  else "", axis=1)
    samps["dedup"] = samps.apply(lambda row: "Deduplicated" if (row.deduplicate=="TRUE")  else "", axis=1)
    samps["splice_prefix"] = "All" 
    samp_ls=None
    samp_ls=[]
    samp_ls.append(samps)
    dedup_both=samps[samps.deduplicate=="BOTH"]
    dedup_both["dedup"]="Deduplicated"
    samp_ls.append(dedup_both)
    samps=pd.concat(samp_ls).drop_duplicates()
    demulti_both=samps[(samps.demultimap=="BOTH")]
    demulti_both["demulti"]="Demultimapped"
    samp_ls.append(demulti_both)
    samps=pd.concat(samp_ls).drop_duplicates()
    spliced=samps[samps.sep_spl==True]
    spliced["splice_prefix"]="Spliced"
    unspliced=spliced.copy(deep=True)
    unspliced["splice_prefix"]="Unspliced"
    samp_ls.append(spliced)
    samp_ls.append(unspliced)
    samps=pd.concat(samp_ls).drop_duplicates()
    unstrand=samps[samps.stranded=="unstranded"]
    unstrand["strand"]="unstranded"
    fwd=samps[samps.stranded=="stranded"]
    fwd["strand"]="fwd"
    rev=fwd.copy(deep=True)
    rev["strand"]="rev"
    results.append(unstrand)
    results.append(fwd)
    results.append(rev)

results = pd.concat(results).drop_duplicates()
results = results.set_index(["sample_name","unit_name"], drop=False).sort_index()

feat_res=[]

dif_exp=features[features.dif_exp==True]
dif_exp["diff"]="expression"
dif_spl=features[features.dif_spl==True]
dif_spl["diff"]="splicing_ratio"

feat_res.append(dif_exp)
feat_res.append(dif_spl)

feat_res=pd.concat(feat_res).drop_duplicates()
feat_res = feat_res.set_index(["feature_name"], drop=False).sort_index()


# Define groups for group analysis
groups = (pd.read_csv(config["groups"], sep="\t", dtype={"protocols": str, "treatment" : str, "splice" : str}, comment="#"))
groups.columns=groups.columns.str.strip()
groups = groups.applymap(lambda x: x.strip() if isinstance(x, str) else x)
groups = groups.mask(groups == '')

groups["experiments"]=groups.apply(lambda row:
    ",".join(
        experiments.experiment[
            (experiments.treatment==row.treatment) & \
            (experiments.protocol.isin([x.strip() for x in row.protocols.split(",")]))
        ].tolist()
    ),
    axis=1
)
groups["splice"]=groups.apply(lambda row:
    "Spliced" if str(row.splice)=="TRUE" else "Unspliced" if str(row.splice)=="FALSE" else "All",
    axis=1
)
groups["splice"]=groups.apply(lambda row:
    ",".join([row.splice] * len(row.experiments.split(","))),
    axis=1
)
groups["feature"]=groups.apply(lambda row:
    ",".join([row.feature] * len(row.experiments.split(","))),
    axis=1
)
groups["group"]=groups.apply(lambda row:
    ",".join([str(row.group)] * len(row.experiments.split(","))),
    axis=1
)
groups["difference"]=groups.apply(lambda row:
    ",".join([row.difference] * len(row.experiments.split(","))),
    axis=1
)
groups=pd.DataFrame(
    list(
        zip(
            (",".join(groups.group.tolist())).split(","),
            (",".join(groups.experiments.tolist())).split(","),
            (",".join(groups.splice.tolist())).split(","),
            (",".join(groups.difference.tolist())).split(","),
            (",".join(groups.feature.tolist())).split(","),
        )
    ),
    columns=["group","experiment","splice","difference","feature"]
)

groups[["control","treatment","protocol","reference","norm_feat","normaliser","paired","diff_lineage","demulti","dedup"]]=groups.apply(lambda row:
    experiments.loc[row.experiment,["control","treatment","protocol","reference","norm_feat","normaliser","paired","diff_lineage","demulti","dedup"]],
    axis=1
)

groups["group_title"]=groups.apply(lambda row: 
    "_".join(list(filter(lambda i: i!=None,groups[groups.group==str(row.group)][["treatment","protocol","feature","difference"]].apply(lambda col: pd.unique(col).item() if len(pd.unique(col))==1 else None)))),
    axis=1
)

groups["title"]=groups.apply(lambda row:
    "_".join(list(filter(lambda i: i!=None,groups[groups.group==str(row.group)][["treatment","protocol","feature","difference"]].apply(lambda col: row[col.name] if len(pd.unique(col))>1 else None)))),
    axis=1
)
groups["type"]=groups.apply(lambda row:
    features.loc[row.feature,"type"] if (row.feature in features["feature_name"].tolist()) else "gtf",
    axis=1
)
groups["sample_source"]=groups.apply(lambda row: get_sample_source(row.experiment),axis=1)
groups["valid"]=groups.apply(lambda row: "annotated" if str(row.diff_lineage)=="genome" else "validated",axis=1)
groups["group_lineage"]=groups.apply(lambda row: row.diff_lineage if len(pd.unique(groups[groups.group==row.group]["diff_lineage"]))==1 else "genome", axis=1)

groups["base"]=groups.apply(lambda row: features.loc[row.feature,"feature"],axis=1)
groups["base_length"]=groups.apply(lambda row: (features.loc[row.feature,"section"]!="body") & (features.loc[row.feature,"is_main_int"]),axis=1)


groups["descript"]=groups.apply(lambda row: 
    "Differences in " + re.sub("[_\.]"," ",str(row.experiment)) + " " + str(row.splice) + " " + str(row.demulti) + " " + str(row.dedup) + " " + str(row.normaliser) + " " + str(row.norm_feat) + " normalised " + str(row.difference) + " of " + str(row.feature) + ", which are based on genes annotated in "  + capitalize(re.sub("[_\.]"," ",str(row.reference))) + " genome" + ( (", and validated in " + str(row.diff_lineage) + " data.") if row.valid=="validated" else ".") ,axis=1)

groups["md5"]=groups.apply(lambda row:
    hashlib.md5((",".join(groups[groups.group==str(row.group)].apply(lambda i: ",".join(str(i))))).encode('utf-8')).hexdigest(),
    axis=1
)
groups["gene_sets"]=groups.apply(lambda row: experiments.loc[row.experiment,"gene_sets"],axis=1)
groups=groups.set_index(["md5"],drop=False)

group_config=groups[["group","md5","group_lineage","group_title"]].drop_duplicates()
group_config["sample_source"]=group_config.apply(lambda row: pd.unique(groups.loc[row.md5,"sample_source"]).item() if len(pd.unique(groups.loc[row.md5,"sample_source"]))==1 else None, axis=1)
group_config["valid"]=group_config.apply(lambda row: "annotated" if row.sample_source=="genome" else "validated", axis=1)
group_config["gene_sets"]=group_config.apply(lambda row: 
    ",".join(list(filter(lambda i: i!='',pd.unique(",".join(groups.loc[row.md5,"gene_sets"].tolist()).split(",")).tolist()))),
    axis=1
)

groups=groups[groups.md5.isin(group_config[group_config.sample_source!=None]["md5"])]

#Set variables for result filenames
DEMULTI=['Demultimapped',''] if config["remove_multimappers"]=='BOTH' else 'Demultimapped' if config["remove_multimappers"] else ''
DEDUP=['Deduplicated',''] if config["remove_duplicate_reads"]=='BOTH' else 'Deduplicated' if config["remove_duplicate_reads"] else ''
SPLICE=['All','Spliced','Unspliced'] if config['seperate_spliced_reads'] else 'All'

VALID='validated' if config['features']['validate_features'] else 'annotated'
TAG=list(config['ensembl_tags'])
STRAND_BIGWIG=""

#STRAND_BIGWIG=['unstranded','fwd','rev'] if config['strand_specific_bigwigs'] else ['unstranded']
#STRAND_META=['stranded']

USAGE=[""] 
#USAGE=["Exon","Intron","ExonVariants"]

BIOTYPE=list(config["biotypes"])

SCORE=["sum","per_gene"] if config["metagene"]["norm_per_gene"] else "per_gene"
MX_NORM="norm" if config["metagene"]["norm_per_gene"] else "sum"
MX_MEAN=["median","mean"] if config["metagene"]["norm_to_median"]=="BOTH" else ["median"] if config["metagene"]["norm_to_median"] else ["mean"]

#Functions for generating results
def get_multiqc():
    qc = expand(
        "qc/multiqc/{exp.experiment}.multiqc.html",
        exp=experiments.itertuples()   
    )
    return qc

def get_bams():
    bams = expand(
        "star/{sample.sample_name}-{sample.unit_name}/{splice}Aligned{demulti}{dedup}.sortedByCoord.out.bam",
        sample=samples.itertuples(), counts=COUNTS_BIGWIG, demulti=DEMULTI, dedup=DEDUP,strand=STRAND_BIGWIG, splice=SPLICING
    )
    return bams

def get_norm_bigwigs():
    bigwigs = expand(
        "norm_bw/{sample.group}/{sample.reference}/bigwigs/{sample.pairRep}.{sample.spikein}_{sample.norm_feat}ReadCount_normalised/{sample.splice_prefix}Aligned{sample.demulti}{sample.dedup}/{sample.sample_name}_{sample.unit_name}.{sample.strand}_{sample.splice_prefix}.normalised.bigwig",
        sample=results.itertuples(), splice=SPLICE
    )
    return bigwigs

def get_feature_counts():
    counts = expand(
        "featurecounts/{experiment}/{reference}/{sample.sample_name}-{sample.unit_name}/{splice}Aligned{demulti}{dedup}.{valid}_{tag}.gene.counts.tab",
        sample=samples.itertuples(), valid=VALID, tag=TAG, demulti=DEMULTI, dedup=DEDUP,strand=STRAND_BIGWIG, splice=SPLICING
    ),
    return counts

def get_diffexp_docx():
    counts = expand(
        "diff_plots/{exp.experiment}/{exp.reference}/differential_expression/{exp.pairRep}.{exp.spikein}_{exp.norm_feat}ReadCount_normalised.{mean}_{norm}/{exp.experiment}.{exp.splice_prefix}_Aligned{exp.demulti}{exp.dedup}.{exp.diff_lineage}_{valid}.custom-{feature.prefix_md5}.{tag}.{feature.feature_name}.docx",
        exp=results.itertuples(), valid=VALID, tag=TAG,strand=STRAND_BIGWIG, splice=SPLICE, feature=features[features.dif_exp.tolist()].itertuples(), mean=MX_MEAN, norm=MX_NORM,
    ),
    return counts

def get_diffsplice_docx():
    counts = expand(
        "diff_plots/{exp.experiment}/{exp.reference}/differential_splicing_ratio/{exp.pairRep}.{exp.spikein}_{exp.norm_feat}ReadCount_normalised.{mean}_{norm}/{exp.experiment}.{exp.splice_prefix}_Aligned{exp.demulti}{exp.dedup}.{exp.diff_lineage}_{valid}.custom-{feature.prefix_md5}.{tag}.{feature.feature_name}.docx",
        exp=results.itertuples(), valid=VALID, tag=TAG,strand=STRAND_BIGWIG, splice=SPLICE, feature=features[features.dif_spl.tolist()].itertuples(), mean=MX_MEAN, norm=MX_NORM,
    ),
    return counts

def get_differential_reports():
    docx = expand(
        "diff_reports/experiment_reports/{exp.experiment}.{exp.diff_lineage}.{tag}.{exp.pairRep}.{exp.spikein}.{exp.norm_feat}_normalised.{mean}_{norm}/{exp.experiment}.{exp.splice_prefix}_Aligned{exp.demulti}{exp.dedup}.differential_report.docx",
        exp=results.itertuples(), valid=VALID, tag=TAG, splice=SPLICE, mean=MX_MEAN, norm=MX_NORM
    ),
    return docx

def get_group_reports():
    docx = expand(
        "group_reports/Group.{group.group}.{group.group_title}.{group.md5}.{tag}.docx",
         group=group_config.itertuples(),
         tag=TAG,
    ),
    return docx

def get_meta_data():
    counts = expand(
        "meta_data/{exp.experiment}/{exp.reference}/differential_expression/{exp.pairRep}.{exp.spikein}_{exp.norm_feat}ReadCount_normalised/{exp.experiment}.{exp.splice_prefix}_Aligned{exp.demulti}{exp.dedup}.{exp.diff_lineage}_{valid}.custom-{feature.prefix_md5}.{tag}.{feature.feature_name}.{mean}_{norm}.mx_data.tab",
        exp=results.itertuples(), valid=VALID, tag=TAG,strand=STRAND_BIGWIG, splice=SPLICE, feature=features[features.dif_exp.tolist()].itertuples(), mean=MX_MEAN, norm=MX_NORM
    ),
    return counts

def get_heat_data():
    counts = expand(
        "heat_data/{exp.experiment}/{exp.reference}/differential_expression/{exp.pairRep}.{exp.spikein}_{exp.norm_feat}ReadCount_normalised/{exp.experiment}.{exp.splice_prefix}_Aligned{exp.demulti}{exp.dedup}.{exp.diff_lineage}_{valid}.custom-{feature.prefix_md5}.{tag}.{feature.feature_name}.heat_data.tab",
        exp=results.itertuples(), valid=VALID, tag=TAG,strand=STRAND_BIGWIG, splice=SPLICE, feature=features[features.dif_exp.tolist()].itertuples()
    ),
    return counts

def get_differential_exp():
    diff_exp = expand(
        "results/{contrast.experiment}/differential_expression/{contrast.contrast}/{contrast.normaliser}_normalised.{splice}Aligned{demulti}{dedup.{feature}",
        contrast=contrasts.itertuples(), dedup=DEDUP, feature=features[features.diffexp==True].itertuples(), demulti=DEMULTI, splice=SPLICING,
    ),
    return diff_exp

def get_differential_usage():
    diff_use=expand(
        "results/{contrast.experiment}/differential_usage/{usage}_Usage/{splice}Aligned{demulti}{dedup}",
        contrast=contrasts.itertuples(), dedup=DEDUP, demulti=DEMULTI, splice=SPLICING, usage=USAGE,
    ),   
    return diff_use

def get_differential_splicing():
    diffsplice = expand(
        "results/{sample.experiment}/differential_splicing/AllAligned{demulti}{dedup}",
        sample=samples.itertuples(), dedup=DEDUP, demulti=DEMULTI, splice=SPLICING,
    )
    return diffsplice if config["differential_analysis"]["splicing"] else ""

def get_meta_profiles():
    pdf = expand(
        "results/{sample.experiment}/meta_profiles/{splice}Aligned{demulti}{dedup}_normalised_by_{sample.normaliser}_{counts}/{biotype}_{feature.feature}/{splice}_{biotype}_{feature}_normalised_{score}.{strand}.meta_profile.pdf",
        sample=samples.itertuples(), counts=COUNTS_BIGWIG, demulti=DEMULTI, dedup=DEDUP,strand=STRAND_META,biotype=BIOTYPE,feature=features.itertuples,score=SCORE, splice=SPLICING,
    )
    return pdf

#Functions for processing
def get_part_bin_number(feature,part):     
    if part == "main" :
        return int(features.loc[feature,"bin_n"])
    else :
        if "x" in str(features.loc[feature,part]) :
            part_bin = float(re.search(r".*(?<!x)",str(features.loc[feature,part])).group()) * features.loc[feature,"bin_n"]
        else :
            part_bin = int(features.loc[feature,part])/int
        if pd.isna(features.loc[feature,"befaftr"]) :
            return int(part_bin)
        else : 
           bef=float(re.search(r".*(?<!x)",str(features.loc[feature,"plotbef"])).group()) if "x" in str(features.loc[feature,"plotbef"]) else int(features.loc[feature,"plotbef"])
           aft=float(re.search(r".*(?<!x)",str(features.loc[feature,"plotaft"])).group()) if "x" in str(features.loc[feature,"plotaft"]) else int(features.loc[feature,"plotaft"])
           
           part_bin = float(features.loc[feature,part]) * int(features.loc[feature,"bin_n"]) * float(features.loc[feature,"befaftr"]) / (float(features.loc[feature,"plotaft"]) + float(features.loc[feature,"plotbef"]))
        return int(part_bin)

def feature_descript(wildcards):
    descript = ( str(wildcards.feature) + " is based on " + "{v}"  + " " + str(wildcards.tag) + " {root}s" + "{ref}." ).format(
        v="annotated" if wildcards.valid=="annotated" else ( str( experiments.loc[wildcards.experiment].squeeze(axis=0)["sample_lineage"] ) + " validated" ) if wildcards.valid=="validated" else "provided",
        root=get_root_feature(wildcards.feature),
        ref="" if wildcards.valid=="provided" \
            else \
            (" in local " + \
            str(experiments.loc[wildcards.experiment].squeeze(axis=0)["sample_species"]).capitalize().replace("_"," ") + \
            " genome")  if not \
            pd.isna(references.loc[experiments.loc[wildcards.experiment].squeeze(axis=0)["sample_species"]].squeeze(axis=0)["annotation_dir"]) \
            else \
            (" in Ensembl " +
            str(experiments.loc[wildcards.experiment].squeeze(axis=0)["sample_species"]).capitalize().replace("_"," ") + " " + \
            str( references.loc[experiments.loc[wildcards.experiment].squeeze(axis=0)["sample_species"]].squeeze(axis=0)["ensembl_build"] )  + " release " + \
            str(references.loc[experiments.loc[wildcards.experiment].squeeze(axis=0)["sample_species"]].squeeze(axis=0)["ensembl_release"]) + \
            " genome" + \
            ("" if wildcards.valid=="validated" \
                else \
                " with transcript support level " + \
                str(features.loc[wildcards.feature].squeeze(axis=0)["tsl"]) + \
                " or better" )
            ),
    )
    return descript
 
def get_cutadapt_input(wildcards):
    sample = samples.loc[wildcards.sample].loc[wildcards.unit].squeeze(axis=0)

    if pd.isna(sample["fq1"]):
        # SRA sample (always paired-end for now)
        accession = sample["sra"]
        return expand("sra/{accession}_{read}.fastq", accession=accession, read=[1, 2])

    if sample["fq1"].endswith("gz"):
        ending = ".gz"
    else :
        ending = ""

    if pd.isna(sample["fq2"]):
        # single end local sample
        return "reads/raw/{S}_{U}_fq1_raw.fastq{E}".format(
            S=sample.sample_name, U=sample.unit_name, E=ending
        )
    else :
         # paired end local sample
        return expand(
            "reads/raw/{S}_{U}_{{read}}_raw.fastq{E}".format(
                 S=sample.sample_name, U=sample.unit_name, E=ending
            ),
            read=["fq1", "fq2"],
        )

def get_rmats_sample_libtype(sample,unit):
    strand=samples.loc[sample].loc[unit,"strandedness"]
    lib="fr-secondstrand" if strand=="reverse" else "fr-firststrand" if strand=="yes" else "fr-unstranded"
    return lib

def get_rmats_exp_libtype(experiment):
    sample=get_experiment_samples(experiment)
    strands=pd.unique(sample.strandedness.tolist())
    strand=strands.item() if len(strands)==1 else "no"
    lib="fr-secondstrand" if strand=="reverse" else "fr-firststrand" if strand=="yes" else "fr-unstranded"
    return lib

def sort_raw_reads(wildcards):
    files = list(
        sorted(glob.glob(samples.loc[wildcards.sample].loc[wildcards.unit].squeeze(axis=0)[wildcards.fq]))
    )
    assert len(files) > 0
    return files

def get_raw_fastqc_input(wildcards):
    sample = samples.loc[wildcards.sample].loc[wildcards.unit]
    
    if sample["fq1"].endswith("gz"):
        ending = ".gz"
    else :
        ending = ""

    if pd.isna(sample["fq2"]):
        # single end sample
        return "reads/raw/{S}_{U}_fq1_raw.fastq{E}".format(
            S=sample.sample_name, U=sample.unit_name, E=ending
        )
    else :
        # paired end sample
        return "reads/raw/{S}_{U}_{fq}_raw.fastq{E}".format(
                S=sample.sample_name, U=sample.unit_name, E=ending, fq=wildcards.fq,
        )



def get_fq(wildcards):
    if not pd.isna(samples.loc[(wildcards.sample,wildcards.unit),"adapters"]):
        # activated trimming, use trimmed data
        if is_paired_end(wildcards.sample):
            # paired-end sample
            return dict(
                zip(
                    ["fq1", "fq2"],
                    expand(
                        "reads/trimmed/{sample}_{unit}_{group}_trimmed.fastq.gz",
                        group=["fq1", "fq2"],
                        **wildcards,
                    )
                )
            )
        # single end sample
        return {"fq1": "reads/trimmed/{sample}_{unit}_single.fastq.gz".format(**wildcards)}
    else :
        # no trimming, use raw reads
        u = samples.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
        if not is_paired_end(wildcards.sample):
            return {"fq1": f"{u.fq1}"}
        else :
            return {"fq1": f"{u.fq1}", "fq2": f"{u.fq2}"}

def get_salmon_fq(wildcards):
    sample = samples.loc[wildcards.sample].loc[wildcards.unit].squeeze(axis=0)

    if pd.isna(sample["fq1"]):
        accession = sample["sra"]
        return expand("sra/{accession}_{read}.fastq", accession=accession, read=[1, 2])

    if sample["fq1"].endswith("gz"):
        ending = ".gz"
    else :
        ending = ""
    if not pd.isna(samples.loc[(wildcards.sample,wildcards.unit),"adapters"]):
        trim="trimmed"
        if is_paired_end(wildcards.sample):
            single="trimmed"
        else:
            single="single"
    else:
        trim="raw"
        single="raw"  
    if pd.isna(sample["fq2"]):
        return "reads/{Trim}/{S}_{U}_fq1_{Single}.fastq{E}".format(
            S=sample.sample_name, U=sample.unit_name, Trim=trim, Single=single, E=ending
        )
    else :
        return expand(
            "reads/{Trim}/{S}_{U}_{{read}}_{Single}.fastq{E}".format(
                 S=sample.sample_name, U=sample.unit_name, Trim=trim, Single=single, E=ending
            ),
            read=["fq1", "fq2"],
        )
    
def get_salmon_input(wildcards, input):
    def add_gunzip(a):
       return "<(gunzip -c " + str(a) + ")"
    if input.fq[0].endswith("gz"):
       fqs=[add_gunzip(i) for i in input.fq]
    if is_paired_end(wildcards.sample):
       return "-1 " + fqs[0] + " -2 " + fqs[1]
    else :
       return "-r " + fqs[0]


def get_strandedness(samples):
    if "strandedness" in samples.columns:
        return samples["strandedness"].tolist()
    else :
        strand_list = ["none"]
        return strand_list * samples.shape[0]

def get_sample_strandedness(wildcards):
      if samples.loc[wildcards.sample].loc[wildcards.unit].squeeze(axis=0)['strandedness'] == "yes" :
          return 1
      else :
          if samples.loc[wildcards.sample].loc[wildcards.unit].squeeze(axis=0)['strandedness'] == "reverse":
              return 2
          else :
              return 0

def get_read_flag(w):
    strand=get_sample_strandedness(w)
    if reads.loc[w.read,"sn_pos"] > 0 :
        flag=128 if strand==2 else 64
    elif reads.loc[w.read,"sn_pos"] < 0 :
        flag=128 if strand!=2 else 64
    return flag
        

def get_bamcov_options(w):
    strand = get_sample_strandedness(w)
    paired = 1 if is_paired_end(w.sample) else 0
    offset_mod = -1 if ((strand==2 and paired==0) or (paired==1 and reads.loc[w.read,"sn_pos"] < 0) else 1
    offset = ("--Offset " + str(reads.loc[w.read,"sn_pos"]*mod) + " " + str(reads.loc[w.read,"sn_pos"]*mod) )if reads.loc[w.read,"single_nuc"] else ""
    saminc = ("--samFlagInclude " + str(get_read_flag(w))) + if (paired==1 and reads.loc[w.read,"single_nuc"]) else ""
    options = offset + " " + saminc
    return options
    
def get_fc_sn_opt(w):
    strand = get_sample_strandedness(w)
    paired = 1 if is_paired_end(w.sample) else 0
    offset_mod = -1 if ((strand==2 and paired==0) or (paired==1 and reads.loc[w.read,"sn_pos"] < 0) else 1
    offset = reads.loc[w.read,"sn_pos"]*mod
    offset_val = abs(offset)
    pos = "--read2pos " + ("5" if offset > 0 else "3")
    shift_type = "--readShiftType " + ("downstream" if offset > 0 else "upstream")
    options = end + " " + shift_type
    return options 

def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))


def get_fastqs(wc):
    if config["trimming"]["activate"]:
        return expand(
            "reads/trimmed/{sample}/{unit}_{read}.fastq.gz",
            unit=samples.loc[wc.sample, "unit_name"],
            sample=wc.sample,
            read=wc.read,
        )
    sample = samples.loc[wc.sample]
    if all(pd.isna(unit["fq1"])):
        # SRA sample (always paired-end for now)
        accession = sample["sra"]
        return expand(
            "sra/{accession}_{read}.fastq", accession=accession, read=wc.read[-1]
        )
    fq = "fq{}".format(wc.read[-1])
    return samples.loc[wc.sample, fq].tolist()


