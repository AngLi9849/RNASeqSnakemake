# All Python variables and functions to be used in the whole pipeline

import glob

import pandas as pd
import hashlib as hashlib
from snakemake.remote import FTP
from snakemake.utils import validate

ftp = FTP.RemoteProvider()

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
    else :
        return "ensembl_{S}.{B}.{R}".format(S=species,B=references.loc[species,"ensembl_build"],R=references.loc[species,"ensembl_release"])


def get_genome(species):
    if not pd.isna(references.loc[species,"genome_dir"]):
        return "resources/genomes/local_{S}_genome.fasta".format(S=species)
    else :
        return "resources/genomes/ensembl_{S}.{B}.{R}_genome.fasta".format(S=species,B=references.loc[species,"ensembl_build"],R=references.loc[species,"ensembl_release"])


def get_annotation(species):
    if not pd.isna(references.loc[species,"annotation_dir"]):
        return "resources/annotations/local_{S}_genome.gtf".format(S=species)
    else :
        return "resources/annotations/ensembl_{S}-{B}-{R}_genome.gtf".format(S=species,B=references.loc[species,"ensembl_build"],R=references.loc[species,"ensembl_release"])

# Read feature_parameters config table into pandas dataframe
features = (pd.read_csv(config["feature_parameters"], sep="\t", dtype={"feature_name": str, "feature": str, "subfeat": str , "annotation_bed6": str}, comment="#"))
features.columns = features.columns.str.strip()
features = features.applymap(lambda x: x.strip() if isinstance(x, str) else x)
features = features.mask( features == '')
features = (features.set_index(["feature_name"], drop=False).sort_index())
feat_prefix = features.iloc[:, features.columns.get_loc("region"): features.columns.get_loc("plotbef")]
feat_plot = features.iloc[:,features.columns.get_loc("region"): features.columns.get_loc("annotation_bed")]
features["prefix"] = feat_prefix.apply(lambda row: 
    ''.join([str(a)+ "_" + str(b) + "." for a,b in zip(row.index.tolist(), row.tolist())]) 
    , axis=1) 
features["plot_prefix"] = feat_plot.apply(lambda row:
    ''.join([str(a)+ "_" + str(b) + "." for a,b in zip(row.index.tolist(), row.tolist())])
    , axis=1)

features['prefix_md5'] = [hashlib.md5(i.encode('utf-8')).hexdigest() for i in features.prefix.tolist()]
features['plot_md5'] = [hashlib.md5(i.encode('utf-8')).hexdigest() for i in features.plot_prefix.tolist()]

# Helper functions for features

def get_feature_validity(feature):
    feat=features.loc[feature].squeeze(axis=0)["feature"]
    if not pd.isna(features.loc[feature].squeeze(axis=0)["annotation_bed"]):
        return "provided"
    else:
        if feat in features["feature_name"].tolist():
            return get_feature_validity(feat)
        else:
            if config["features"]["validate_features"]:
                return "validated"
            else:
                return "annotated"

def get_root_feature(feature):
    base=features.loc[feature].squeeze(axis=0)["feature"]
    if not pd.isna(features.loc[feature].squeeze(axis=0)["annotation_bed"]):
        return str(feature)
    else:
        if base in features["feature_name"].tolist():
            return get_root_feature(base)
        else:
            return str(features.loc[feature].squeeze(axis=0)["feature"])

def get_prefixed_feature(f):
    

# Read experimentss config table into pandas dataframe
experiments = (pd.read_csv(config["experiments"], sep="\t", dtype={"protocol": str, "sample_lineage" : str}, comment="#"))
experiments.columns=experiments.columns.str.strip()
experiments = experiments.applymap(lambda x: x.strip() if isinstance(x, str) else x)
experiments = experiments.mask(experiments == '')
experiments["experiment"] = experiments.apply( lambda row: \
    ( str( row.treatment ) + "_vs_" + str( row.control ) + "_" + str( row.protocol ) ), axis = 1
)
experiments = (experiments.set_index(["experiment"], drop=False).sort_index())


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

def get_experiment_samples(wildcards):
    exp = experiments.loc[wildcards.experiment].squeeze()
    cond = [exp.control,exp.treatment]
    sample = samples[samples.protocol == exp.protocol][samples.condition.isin(cond)]
    return sample

def get_lineage_sj_samples(wildcards):
    exp = experiments[ (experiments.sample_lineage==wildcards.lineage).tolist() & experiments.splice & (experiments.sample_species == wildcards.species).tolist() ]
    cond = [exp.control,exp.treatment]
    sample=samples[(samples.condition + "_" + samples.protocol).isin(
        (exp["control"] + "_" + exp["protocol"]).tolist() + (exp["treatment"] + "_" + exp["protocol"]).tolist()
    )]
    return sample

#Set variables for result filenames
DEMULTI=['Demultimapped',''] if config["remove_multimappers"]=='BOTH' else 'Demultimapped' if config["remove_multimappers"] else ''
DEDUP=['Deduplicated',''] if config["remove_duplicate_reads"]=='BOTH' else 'Deduplicated' if config["remove_duplicate_reads"] else ''
SPLICE=['All','Spliced','Unspliced'] if config['seperate_spliced_reads'] else 'All'

VALID=['validated'] if config['features']['validate_features'] else ['annotated']
TAG=list(config['ensembl_tags'])

STRAND_BIGWIG=['unstranded','fwd','rev'] if config['strand_specific_bigwigs'] else ['unstranded']
STRAND_META=['stranded']

USAGE=[""] 
#USAGE=["Exon","Intron","ExonVariants"]

BIOTYPE=list(config["biotypes"])

SCORE=["sum","per_gene"] if config["metagene"]["plot_sum"] else "per_gene"

#Functions for generating results
def get_bams():
    bams = expand(
        "star/{sample.sample_name}-{sample.unit_name}/{splice}Aligned{demulti}{dedup}.sortedByCoord.out.bam",
        sample=samples.itertuples(), counts=COUNTS_BIGWIG, demulti=DEMULTI, dedup=DEDUP,strand=STRAND_BIGWIG, splice=SPLICING
    )
    return bams

def get_bigwigs():
    bigwigs = expand(
        "results/{sample.experiment}/bigwig/{splice}Aligned{demulti}{dedup}_normalised_by_{sample.normaliser}_{normaliser}/{sample.sample_name}_{sample.unit_name}.{strand}_{splice}.bigwig",
        sample=samples.itertuples(), demulti=DEMULTI, dedup=DEDUP,strand=STRAND_BIGWIG, splice=SPLICING
    )
    return bigwigs

def get_feature_counts():
    counts = expand(
        "featurecounts/{experiment}/{reference}/{sample.sample_name}-{sample.unit_name}/{splice}Aligned{demulti}{dedup}.{valid}_{tag}.gene.counts.tab",
        sample=samples.itertuples(), valid=VALID, tag=TAG, demulti=DEMULTI, dedup=DEDUP,strand=STRAND_BIGWIG, splice=SPLICING
    ),
    return counts

def get_genebody_diffexp_docx():
    counts = expand(
        "results/{exp.experiment}/{exp.reference}/differential_expression/{exp.normaliser}_{exp.norm_feat}ReadCount_normalised.{splice}_Aligned{demulti}{dedup}.genome_{valid}.{tag}.GeneBody.docx",
        exp=experiments.itertuples(), valid=VALID, tag=TAG, demulti=DEMULTI, dedup=DEDUP,strand=STRAND_BIGWIG, splice=SPLICE,  
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
def feature_descript(wildcards):
    descript = ( str(wildcards.feature) + " is based on " + "{v}"  + " " + str(wildcards.tag) + " {root}s" + "{ref}." ).format(
        v="annotated" if wildcards.valid=="annotated" else ( str( experiments.loc[wildcards.experiment].squeeze(axis=0)["sample_lineage"] ) + "validated" ) if wildcards.valid=="validated" else "provided",
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
            " genome with transcript support level " + \
            str(features.loc[wildcards.feature].squeeze(axis=0)["tsl"]) + \
            " or better" ),
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


