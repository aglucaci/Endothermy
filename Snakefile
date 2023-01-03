"""
@Description: Fish endothermy examination
@Authors: Alexander Lucaci, Avery Selberg
@Version: 2023.1
"""

#----------------------------------------------------------------------------
# Imports
#----------------------------------------------------------------------------
import os
import sys
import json
import csv
from pathlib import Path
import glob

#----------------------------------------------------------------------------
# Declares
#----------------------------------------------------------------------------
with open("config.json", "r") as input_sc:
  config = json.load(input_sc)
#end with

with open("cluster.json", "r") as input_c:
  cluster = json.load(input_c)
#end with

#----------------------------------------------------------------------------
# User Settings
#----------------------------------------------------------------------------
BASEDIR = os.getcwd()
print("# We are operating out of base directory:", BASEDIR)

# Which project are we analyzing?
LABEL = config["Label"]

DATA_DIRECTORY = os.path.join(BASEDIR, "data", "CleanCandidateGenes")
EXTENSION = "fasta"

PARTITION_LIST = config["partition"].split(",")
print("# Testing Partition(s):", PARTITION_LIST)

CANDIDATE_GENES = [os.path.basename(x) for x in glob.glob(os.path.join(DATA_DIRECTORY, "*." + EXTENSION))]
print("# We will process", len(CANDIDATE_GENES), "candidate gene files")

#----------------------------------------------------------------------------
# Software Settings
#----------------------------------------------------------------------------
HYPHY    = "hyphy"
HYPHYMPI = "HYPHYMPI"
HYPHY_ANALYSES = os.path.join(BASEDIR, "hyphy-analyses")
STRIKE_AMBIGS_BF = os.path.join(BASEDIR, "scripts", "strike-ambigs.bf")
OUTDIR = os.path.join(BASEDIR, "results", LABEL)

#hyphy-analyses/remove-duplicates/remove-duplicates.bf
REMOVE_DUPS_BF = os.path.join(HYPHY_ANALYSES, "remove-duplicates", "remove-duplicates.bf")
BUSTED_PH_BF = os.path.join(HYPHY_ANALYSES, "BUSTED-PH", "BUSTED-PH.bf")
LABEL_TREE_BF = os.path.join(BASEDIR, "hyphy-analyses", "LabelTrees", "label-tree.bf")

# Create output directories
Path(os.path.join(BASEDIR, "results")).mkdir(parents=True, exist_ok=True)
Path(OUTDIR).mkdir(parents=True, exist_ok=True)

# Settings, these can be passed in or set in a config.json type file
PPN = cluster["__default__"]["ppn"] 

#----------------------------------------------------------------------------
# rule all
#----------------------------------------------------------------------------
rule all:
    input:
        expand(os.path.join(OUTDIR, "{GENE}.SA"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.SA.noDups"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.SA.noDups.dst"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.SA.noDups.raxml.bestTree"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.SA.noDups.BUSTEDS.json"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.SA.noDups.BUSTED.json"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.SA.noDups.BUSTED+MH.json"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.SA.noDups.BUSTEDS+MH.json"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.SA.noDups.Labeled_{P}.nwk"), GENE=CANDIDATE_GENES, P=PARTITION_LIST),
        expand(os.path.join(OUTDIR, "{GENE}.SA.noDups.Labeled_{P}.BUSTEDS+PH.json"), GENE=CANDIDATE_GENES, P=PARTITION_LIST)
        #expand(os.path.join(OUTDIR, "{GENE}.SA.noDups.GARD.json"), GENE=CANDIDATE_GENES),
        #expand(os.path.join(OUTDIR, "{GENE}.SA.noDups.meme-all.json"), GENE=CANDIDATE_GENES),
        #expand(os.path.join(OUTDIR, "{GENE}.SA.noDups.meme_{P}.json"), GENE=CANDIDATE_GENES, P=PARTITION_LIST),
        #expand(os.path.join(OUTDIR, "{GENE}.SA.noDups.cfel-{P}.json"), GENE=CANDIDATE_GENES, P=PARTITION_LIST),
        #expand(os.path.join(OUTDIR, "{GENE}.SA.noDups.RELAX-{P}.json"), GENE=CANDIDATE_GENES, P=PARTITION_LIST),
        #expand(os.path.join(OUTDIR, "{GENE}.SA.noDups.aBSREL.json"), GENE=CANDIDATE_GENES),
        #expand(os.path.join(OUTDIR, "{GENE}.SA.noDups.aBSREL-FG.json"), GENE=CANDIDATE_GENES)
 #end input
#end rule

#----------------------------------------------------------------------------
# Individual rules
#----------------------------------------------------------------------------
rule strike_ambigs:
   input:
       input_msa = os.path.join(DATA_DIRECTORY, "{GENE}")
   output:
       output = os.path.join(OUTDIR, "{GENE}.SA")
   shell:
      "{HYPHY} {STRIKE_AMBIGS_BF} --alignment {input.input_msa} --output {output.output}"
#end rule

rule remove_duplicates:
    input:
        input = rules.strike_ambigs.output.output
    output:
        output = os.path.join(OUTDIR, "{GENE}.SA.noDups")
    shell:
       "{HYPHY} {REMOVE_DUPS_BF} --msa {input.input} --output {output.output} ENV='DATA_FILE_PRINT_FORMAT=9'"
#end rule

rule tn93:
    input:
       input = rules.remove_duplicates.output.output
    output:
       output = os.path.join(OUTDIR, "{GENE}.SA.noDups.dst")
    shell:
       "tn93 -t 1 -o {output.output} {input.input}"
#end rule

rule raxml:
    params:
        THREADS = PPN
    input:
        input = rules.remove_duplicates.output.output
    output:
        output = os.path.join(OUTDIR, "{GENE}.SA.noDups.raxml.bestTree")
    shell:
        "raxml-ng --model GTR --msa {input.input} --threads {params.THREADS} --tree pars{{5}} --force"
#end rule 

#----------------------------------------------------------------------------
# Recombination detection
#----------------------------------------------------------------------------

#rule recombination:
#    input: 
#        input = rules.remove_duplicates.output.output
#    output: 
#        output = os.path.join(OUTDIR, "{GENE}.SA.noDups.GARD.json")
#    shell: 
#        "mpirun -np {PPN} {HYPHYMPI} LIBPATH={RES} GARD --alignment {input.input} --rv GDD --output {output.output}"
#end rule

#----------------------------------------------------------------------------
# Label tree partitions for BUSTED-PH
#----------------------------------------------------------------------------

rule label_tree:
    input:
        tree = rules.raxml.output.output,
        partition_file = os.path.join(BASEDIR, "data", "Partitions", "{P}.txt")    
    output:
        output = os.path.join(OUTDIR, "{GENE}.SA.noDups.Labeled_{P}.nwk")
    shell:
        "{HYPHY} {LABEL_TREE_BF} --tree {input.tree} --list {input.partition_file} --output {output.output} --internal-nodes 'Parsimony' --label FOREGROUND" 

#----------------------------------------------------------------------------
# Run BUSTED-PH
#----------------------------------------------------------------------------
rule busted_ph:
    input:
        msa = rules.remove_duplicates.output.output,
        tree = rules.label_tree.output.output
    output:
        output = os.path.join(OUTDIR, "{GENE}.SA.noDups.Labeled_{P}.BUSTEDS+PH.json")
    shell:
        "{HYPHY} {BUSTED_PH_BF} --alignment {input.msa} --tree {input.tree} --srv Yes --starting-points 10 --output {output.output} --branches FOREGROUND"
        
#----------------------------------------------------------------------------
# Run BUSTEDS
#----------------------------------------------------------------------------

rule busteds:
    input:
        msa = rules.remove_duplicates.output.output,
        tree = rules.raxml.output.output
    output:
        output = os.path.join(OUTDIR, "{GENE}.SA.noDups.BUSTEDS.json")
    shell:
        "mpirun -np {PPN} {HYPHYMPI} BUSTED --alignment {input.msa} --tree {input.tree} --srv Yes --starting-points 10 --output {output.output}" 

#----------------------------------------------------------------------------
# Run BUSTED
#----------------------------------------------------------------------------

rule busted:
    input:
        msa = rules.remove_duplicates.output.output,
        tree = rules.raxml.output.output
    output:
        output = os.path.join(OUTDIR, "{GENE}.SA.noDups.BUSTED.json")
    shell:
        "mpirun -np {PPN} {HYPHYMPI} BUSTED --alignment {input.msa} --tree {input.tree} --srv No --starting-points 10 --output {output.output}" 

#----------------------------------------------------------------------------
# Run BUSTED+MH
#----------------------------------------------------------------------------

rule bustedmh:
    input:
        msa = rules.remove_duplicates.output.output,
        tree = rules.raxml.output.output
    output:
        output = os.path.join(OUTDIR, "{GENE}.SA.noDups.BUSTED+MH.json")
    shell:
        "mpirun -np {PPN} {HYPHYMPI} BUSTED --alignment {input.msa} --tree {input.tree} --srv No --starting-points 10 --output {output.output} --multiple-hits Double+Triple" 

#----------------------------------------------------------------------------
# Run BUSTED+S+MH
#----------------------------------------------------------------------------

rule bustedsmh:
    input:
        msa = rules.remove_duplicates.output.output,
        tree = rules.raxml.output.output
    output:
        output = os.path.join(OUTDIR, "{GENE}.SA.noDups.BUSTEDS+MH.json")
    shell:
        "mpirun -np {PPN} {HYPHYMPI} BUSTED --alignment {input.msa} --tree {input.tree} --srv Yes --starting-points 10 --output {output.output} --multiple-hits Double+Triple" 

#----------------------------------------------------------------------------
# End of file
#----------------------------------------------------------------------------
















