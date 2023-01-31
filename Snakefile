"""
@Description: Fish endothermy examination
@Authors: Alexander Lucaci, Avery Selberg
@Version: 2023.1
"""

# =============================================================================
# Imports
# =============================================================================

import os
import sys
import json
import csv
from pathlib import Path
import glob
from Bio import SeqIO

# =============================================================================
# Declares
# =============================================================================

with open("config.json", "r") as input_sc:
  config = json.load(input_sc)
#end with

with open("cluster.json", "r") as input_c:
  cluster = json.load(input_c)
#end with

# =============================================================================
# User Settings
# =============================================================================

BASEDIR = os.getcwd()
print("# We are operating out of base directory:", BASEDIR)

# Which project are we analyzing?
LABEL = config["Label"]

DATA_DIRECTORY = os.path.join(BASEDIR, "data", config["DataDirectory"])
EXTENSION = "fasta"

PARTITION_LIST = config["Partition"].split(",")
print("# Testing Partition(s):", PARTITION_LIST)

CANDIDATE_GENES = [os.path.basename(x) for x in glob.glob(os.path.join(DATA_DIRECTORY, "*." + EXTENSION))]
print("# We will process", len(CANDIDATE_GENES), "candidate gene files")

# =============================================================================
# Software Settings
# =============================================================================

HYPHY    = "hyphy"
HYPHYMPI = "HYPHYMPI"
HYPHY_ANALYSES = os.path.join(BASEDIR, "hyphy-analyses")
STRIKE_AMBIGS_BF = os.path.join(BASEDIR, "scripts", "strike-ambigs.bf")
OUTDIR = os.path.join(BASEDIR, "results", LABEL)

#hyphy-analyses/remove-duplicates/remove-duplicates.bf
REMOVE_DUPS_BF = os.path.join(HYPHY_ANALYSES, "remove-duplicates", "remove-duplicates.bf")
BUSTED_PH_BF = os.path.join(HYPHY_ANALYSES, "BUSTED-PH", "BUSTED-PH.bf")
LABEL_TREE_BF = os.path.join(HYPHY_ANALYSES, "LabelTrees", "label-tree.bf")
PREMSA = os.path.join(HYPHY_ANALYSES, "codon-msa", "pre-msa.bf")
POSTMSA = os.path.join(HYPHY_ANALYSES, "codon-msa", "post-msa.bf")

# Create output directories
Path(os.path.join(BASEDIR, "results")).mkdir(parents=True, exist_ok=True)
Path(OUTDIR).mkdir(parents=True, exist_ok=True)

# Settings, these can be passed in or set in a config.json type file
PPN = cluster["__default__"]["ppn"] 

# =============================================================================
# rule all
# =============================================================================
rule all:
    input:
        expand(os.path.join(OUTDIR, "{GENE}.fa"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.fa_protein.fas"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.fa_nuc.fas"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.fa_protein.aln"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.fa_codons.fasta"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.fa_codons_duplicates.json"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.fa_codons.ID.fasta"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.fa_codons.ID.SA.fasta"), GENE=CANDIDATE_GENES),
        expand(os.path.join(OUTDIR, "{GENE}.fa_codons.ID.SA.fasta.raxml.bestTree"), GENE=CANDIDATE_GENES),
        expand(os.path.join(BASEDIR, "data", "Partitions", "{P}_BG.txt"), P=PARTITION_LIST),
        expand(os.path.join(OUTDIR, "{GENE}.fa_codons.ID.SA.fasta.raxml.bestTree.labeled_fgOnly_{P}.nwk"), GENE=CANDIDATE_GENES, P=PARTITION_LIST),
        expand(os.path.join(OUTDIR, "{GENE}.fa_codons.ID.SA.fasta.raxml.bestTree.labeled_{P}.nwk"), GENE=CANDIDATE_GENES, P=PARTITION_LIST),
        expand(os.path.join(OUTDIR, "{GENE}.fa_codons.ID.SA.fasta.{P}.BUSTEDPH.json"), GENE=CANDIDATE_GENES, P=PARTITION_LIST)
    #end input
#end rule

# =============================================================================
# Individual rules
# =============================================================================

rule clean:
    input:
        input = os.path.join(DATA_DIRECTORY, "{GENE}")
    output:
        output = os.path.join(OUTDIR, "{GENE}.fa")
    shell:
       "bash scripts/cleaner.sh {input.input} {output.output}"
#end rule

#----------------------------------------------------------------------------
# Alignment
#----------------------------------------------------------------------------

rule pre_msa:
    input: 
        codons = rules.clean.output.output
    output: 
        protein_fas = os.path.join(OUTDIR, "{GENE}.fa_protein.fas"),
        nucleotide_fas = os.path.join(OUTDIR, "{GENE}.fa_nuc.fas")
    shell: 
        "mpirun -np {PPN} {HYPHYMPI} {PREMSA} --input {input.codons}"
#end rule 

rule mafft:
    input:
        protein = rules.pre_msa.output.protein_fas
    output:
        protein_aln = os.path.join(OUTDIR, "{GENE}.fa_protein.aln")
    shell:
        "mafft --auto {input.protein} > {output.protein_aln}"
#end rule

rule post_msa:
    input: 
        protein_aln = rules.mafft.output.protein_aln,
        nucleotide_seqs = rules.pre_msa.output.nucleotide_fas  
    output: 
        codons_fas = os.path.join(OUTDIR, "{GENE}.fa_codons.fasta"),
        duplicates_json = os.path.join(OUTDIR, "{GENE}.fa_codons_duplicates.json")
    shell: 
        "mpirun -np {PPN} {HYPHYMPI} {POSTMSA} --protein-msa {input.protein_aln} --nucleotide-sequences {input.nucleotide_seqs} --output {output.codons_fas} --duplicates {output.duplicates_json}"
#end rule 

#----------------------------------------------------------------------------
# Alignment Quality Control
#----------------------------------------------------------------------------

# Trim the end of the FASTA ID that pre-msa adds
rule trim_fasta_id:
    input:
        input = rules.post_msa.output.codons_fas
    output:
        output = os.path.join(OUTDIR, "{GENE}.fa_codons.ID.fasta")
    run:
        records = []
        with open(input.input, "r") as handle:
            for m, record in enumerate(SeqIO.parse(handle, "fasta")):
                _id = record.id.split("_")
                _id = "_".join(_id[:-1])
                print(record.id, _id)
                record.id = _id
                record.description = _id
                records.append(record)
            #end for
        #end with
        # Write to file
        SeqIO.write(records, output.output, "fasta")
    #end run
#end rule

rule strike_ambigs:
   input:
       input_msa = rules.trim_fasta_id.output.output
   output:
       output = os.path.join(OUTDIR, "{GENE}.fa_codons.ID.SA.fasta")
   shell:
      "{HYPHY} {STRIKE_AMBIGS_BF} --alignment {input.input_msa} --output {output.output}"
#end rule

#----------------------------------------------------------------------------
# Phylogenetic inference
#----------------------------------------------------------------------------

rule raxml_ng:
    params:
        THREADS = PPN
    input:
        input = rules.strike_ambigs.output.output
    output:
        output = os.path.join(OUTDIR, "{GENE}.fa_codons.ID.SA.fasta.raxml.bestTree")
    shell:
        "raxml-ng --model GTR+G --msa {input.input} --threads {params.THREADS} --force --redo"
#end rule 

#----------------------------------------------------------------------------
# Prepare Background Species for each Scenario
#----------------------------------------------------------------------------

rule prepare_bg:
    input:
        foreground_file = os.path.join(BASEDIR, "data", "Partitions", "{P}_FG.txt"),
        nuisance_file = os.path.join(BASEDIR, "data", "Partitions", "{P}_Nuisance.txt"),
        all_species = os.path.join(BASEDIR, "data", "Partitions", "All_Species.txt")
    output:
        background_file = os.path.join(BASEDIR, "data", "Partitions", "{P}_BG.txt")
    run:
        # --- Read in data
        with open(input.all_species, "r") as f:
            all_species = [line for line in f]
        #end with
        with open(input.nuisance_file, "r") as f:
            nu_species = [line for line in f]
        #end with
        with open(input.foreground_file) as f:
            fg_species = [line for line in f]
        #end with
        # --- Create BG partition
        species_no_nu = [item for item in all_species if item not in nu_species]
        bg_species = [item for item in species_no_nu if item not in fg_species]
        # --- Write to file
        with open(output.background_file, "w") as f:
            for line in bg_species:
                f.write(line)
            #end for
        #end with
#end rule

#----------------------------------------------------------------------------
# Label tree partitions for BUSTED-PH
#----------------------------------------------------------------------------

rule label_tree_fg:
    input:
        tree = rules.raxml_ng.output.output,
        partition_file = os.path.join(BASEDIR, "data", "Partitions", "{P}_FG.txt"),   
    output:
        output = os.path.join(OUTDIR, "{GENE}.fa_codons.ID.SA.fasta.raxml.bestTree.labeled_fgOnly_{P}.nwk")
    shell:
        "{HYPHY} {LABEL_TREE_BF} --tree {input.tree} --list {input.partition_file} --output {output.output} --internal-nodes 'Parsimony' --label FOREGROUND" 
#end rule

rule label_tree_bg:
    input:
        tree = rules.label_tree_fg.output.output,
        partition_file = os.path.join(BASEDIR, "data", "Partitions", "{P}_BG.txt")
    output:
        output = os.path.join(OUTDIR, "{GENE}.fa_codons.ID.SA.fasta.raxml.bestTree.labeled_{P}.nwk")
    shell:
        "{HYPHY} {LABEL_TREE_BF} --tree {input.tree} --list {input.partition_file} --output {output.output} --internal-nodes 'Parsimony' --label BACKGROUND" 
#end rule


#----------------------------------------------------------------------------
# Selection analyses
#----------------------------------------------------------------------------

rule busted_ph:
    input:
        msa =  rules.strike_ambigs.output.output,
        tree = rules.label_tree_bg.output.output
    output:
        output = os.path.join(OUTDIR, "{GENE}.fa_codons.ID.SA.fasta.{P}.BUSTEDPH.json")
    shell:
        "{HYPHY} {BUSTED_PH_BF} --alignment {input.msa} --tree {input.tree} --srv Yes --starting-points 10 --output {output.output} --branches FOREGROUND"
#end rule


#----------------------------------------------------------------------------
# End of file
#----------------------------------------------------------------------------















