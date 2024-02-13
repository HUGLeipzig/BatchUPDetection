
##################################################
# UPD-FINDER
# Input: 
#       - a folder full of ROH-RG result files
#       - a folder full of isec files
#       - a list determining the family relations of the roh and isec files
# Process:
#       - all files are read, processed and interpreted
#       - roh coverage per chromosome is calculated
#       - inheritance ratio per chromosome is calculated
# Output:
#       - a table containing the roh coverage and inheritance ratios per sample per chromosome
#       - interactive visualization of the roh coverage and inheritance ratios
#       - tagging of samples/chromosmes as ["likely_consanguinous", "upd_likely", "isodisomy", "heterodisomy"]
# Example usage:
#   python upd_finder.py roh_folder isec_folder family_files vcf_folder result_folder
##################################################

##################################################
# imports
import pandas as pd
import altair as alt
from cyvcf2 import VCF
import os.path
import sys
from multiprocessing import Pool

##################################################
# input arguments

upd_table_out = "upd_finder_overview.xlsx"

roh_folder = sys.argv[1]
isec_folder = sys.argv[2]
family_files = sys.argv[3]
vcf_folder = sys.argv[4]
out_dir = sys.argv[5]

##################################################
# settings
roh_suffix = "_roh.txt"
isec_suffix = "_isec.txt"
roh_cols = ["RG", "sample", "chr", "start", "end", "length", "markers", "qual"]
isec_cols = ["chr", "start", "ref", "alt", "snv_occurence"]
chr_list = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
assembly = "hg38" #"hg19"
genome_file = "data/" + assembly + ".genome"
inheritance_dict_trio = {
                "001"   :   "unknown",
                "010"   :   "unknown",
                "011"   :   "unknown",
                "100"   :   "de novo",
                "101"   :   "paternal",
                "110"   :   "maternal",
                "111"   :   "biparental",
            }

inheritance_dict_duo_mother = {
                "00"   :   "unknown",
                "01"   :   "unknown",
                "10"   :   "not maternal",
                "11"   :   "maternal",
            }

inheritance_dict_duo_father = {
                "00"   :   "unknown",
                "01"   :   "unknown",
                "10"   :   "not paternal",
                "11"   :   "paternal",
            }

# flagging cutoffs:
roh_high_cutoff = 0.7
roh_high_tag = "roh_high"
roh_high_mixed_start = 0.2
roh_high_mixed_end = 0.7
roh_high_mixed_tag = "roh_high_mixed"

inh_ratio_high_trio_cutoff = 2
inh_ratio_high_duo_cutoff = 5
inh_ratio_high_tag = "inh_ratio_high"

consanguin_min_chr_count = 3
consanguin_roh_cutoff = 0.1
consanguin_tag = "likely_consanguinous"
not_consanguin_tag = "unlikely_consanguinous"

min_snvs_per_chr = 200
snv_per_chr_warning = "insufficient_snvs"

errors = 0


##################################################
# functions
def set_up_results_table():
    df_upd_finder = pd.DataFrame(
        columns=["index_file",
            "mother_file",
            "father_file",
            "setup",
            "chr",
            "perc_roh",
            "mat_over_pat",
            "mat_over_notmat",
            "pat_over_mat",
            "pat_over_notpat"
            ]
    )
    return df_upd_finder

def read_vcf_file(vcf_file):
    df_vcf_variants = pd.DataFrame(columns = ["chr", "start", "end"])

    li_chr = []
    li_start = []
    li_end = []

    for variant in VCF(vcf_file):
        
        # some variants have no alt allele, these return "[]" as ALT 
        if not variant.ALT == []:

            li_chr.append(variant.CHROM)
            li_start.append(variant.POS)
            li_end.append(variant.POS + len(variant.REF))
    
    df_vcf_variants["chr"] = li_chr
    df_vcf_variants["start"] = li_start
    df_vcf_variants["end"] = li_end

    df_vcf_variants["start"] = df_vcf_variants["start"].astype(int)
    df_vcf_variants["end"] = df_vcf_variants["end"].astype(int)

    return df_vcf_variants

def get_chromosome_lengths(vcf_file):
    global min_snvs_per_chr
    
    li_chr = []
    li_length=[]

    df_vcf_variants = read_vcf_file(vcf_file)
    chr_list = set(list(df_vcf_variants["chr"]))

    for chr in chr_list:
        df = df_vcf_variants[df_vcf_variants["chr"] == chr]
        
        snvs_per_chr = df.shape[0]
        if snvs_per_chr < min_snvs_per_chr:
            chr_length = 0
        else:
            first_variant = df.iloc[0]["start"]
            last_variant = df.iloc[-1]["start"]
            chr_length = last_variant-first_variant
            
        if not "chr" in str(chr):
            chr = "chr"+str(chr)

        li_chr.append(chr)
        li_length.append(chr_length)

    df_chromosomes = pd.DataFrame()
    df_chromosomes["chr"] = li_chr
    df_chromosomes["length"] = li_length

    df_chromosomes = df_chromosomes.set_index('chr')

    return df_chromosomes

def roh_inh_scatter(df, trait, df_cutoffs):
    
    points = alt.Chart(df).mark_point().encode(
        alt.X('perc_roh:Q', 
                title='perc_roh',
                scale=alt.Scale(domain=[0, 1],
                clamp=True,)
                ),
        alt.Y(trait+":Q",
            title=trait,
            scale=alt.Scale(domain=[0, 20],
            clamp=True)
            ),
        tooltip=["index_file", "chr", "perc_roh", "tags", "consanguinity"],
        color="consanguinity"
    ).properties(
        width=800,
        height=600,
        title=trait
    ).interactive(
    )

    roh_high_rect = alt.Chart(df_cutoffs).mark_rect(opacity=0.1, color="orange").encode(
                x = "roh_high_start",
                x2 = "roh_high_end",
                y = alt.value(0),
                y2 = alt.value(10000)
            )
    roh_mixed_rect = alt.Chart(df_cutoffs).mark_rect(opacity=0.1, color="yellow").encode(
                    x = "roh_mixed_start",
                    x2 = "roh_high_mixed_end",
                    y = alt.value(0),
                    y2 = alt.value(10000)
                )
    inh_duo_rect = alt.Chart(df_cutoffs).mark_rect(opacity=0.1, color="yellow").encode(
                    x = alt.value(0),
                    x2 = alt.value(10000),
                    y = "inh_ratio_high_cutoff_start",
                    y2 = "inh_ratio_high_cutoff_end"
                )    
    inh_trio_rect = alt.Chart(df_cutoffs).mark_rect(opacity=0.1, color="orange").encode(
                    x = alt.value(0),
                    x2 = alt.value(10000),
                    y = "inh_ratio_high_cutoff_end",
                    y2 = "inh_ratio_high_cutoff_end_duos"
                )    

    plot = points + roh_high_rect + roh_mixed_rect + inh_duo_rect + inh_trio_rect

    return plot

def inh_plot(df, df_cutoffs):
    #TODO: color by tag
    roh_plot = alt.Chart(df).mark_point().encode(
        alt.X('chr:N', title='chromosome'),
        alt.Y("perc_roh:Q",
            scale=alt.Scale(domain=[0, 1],
            clamp=True)),
        tooltip=["index_file", "chr", "perc_roh", "tags", "consanguinity"],
        color="consanguinity"
    ).properties(
        width=800,
        height=600,
        title="perc_roh"
    ).interactive(
    )

    roh_high_rect = alt.Chart(df_cutoffs).mark_rect(opacity=0.1, color="orange").encode(
            y = "roh_high_start",
            y2 = "roh_high_end"
        )
    roh_mixed_rect = alt.Chart(df_cutoffs).mark_rect(opacity=0.1, color="yellow").encode(
            #x = "roh_high_start",
            #x2 = "roh_high_end",
            y = "roh_mixed_start",
            y2 = "roh_high_mixed_end"
        )
    plot = roh_plot + roh_high_rect + roh_mixed_rect

    return plot

def collect_roh_inh(row):

    global roh_folder, \
        roh_suffix, \
        isec_suffix, \
        vcf_folder, \
        roh_cols, \
        errors

    li_index_file = []
    li_mother_file = []
    li_father_file = []
    li_setup = []
    li_chr = []
    li_perc_roh = []
    li_mat_over_pat = []
    li_mat_over_notmat= []
    li_pat_over_mat = []
    li_pat_over_notpat = []

    # check if file exists
    if (os.path.isfile(roh_folder+row["index_file"]+roh_suffix)) \
        and (os.path.isfile(isec_folder+row["index_file"]+isec_suffix) \
        and (os.path.isfile(vcf_folder+row["index_file"]))):
        
        # collect chromosome lenghts from vcf file
        df_chromosomes = get_chromosome_lengths(vcf_folder+row["index_file"])

        # determine setup and check for files
        if row["mother_file"] and row["father_file"]:
            setup = "trio"
        elif row["mother_file"] and not row["father_file"]:
            setup = "duo_mother"
        elif not row["mother_file"] and row["father_file"]:
            setup = "duo_father"
        elif not row["mother_file"] and not row["father_file"]:
            setup = "single"
        
        # get roh coverage
        df_roh = pd.read_csv(roh_folder+row["index_file"]+roh_suffix, sep="\t", names=roh_cols)
        if df_roh.shape[0] == 0:
            print("ERROR: ROHs empty: " + roh_folder+row["index_file"]+roh_suffix)
            errors += 1
        else:
            # add "chr" if not present
            if "chr" not in str(df_roh.iloc[0]["chr"]):
                df_roh["chr"] = "chr" + df_roh["chr"].astype(str)

            for chr in chr_list:
                if not chr in df_chromosomes.index or df_chromosomes.at[chr,"length"] == 0:
                    perc_covered_by_rohs = pd.NA
                else:
                    df_roh_chr = df_roh[df_roh["chr"] == chr]
                    chr_length = df_chromosomes.at[chr,"length"]
                    total_lengths_of_rohs = df_roh_chr["length"].sum()
                    perc_covered_by_rohs = total_lengths_of_rohs / chr_length

                li_index_file.append(row["index_file"])
                li_mother_file.append(row["mother_file"])
                li_father_file.append(row["father_file"])
                li_setup.append(setup)
                li_chr.append(chr)
                li_perc_roh.append(perc_covered_by_rohs)
            
            # get inheritance
            if setup == "single":
                mop=""
                monm="" 
                pom="" 
                ponp=""
                for chr in chr_list:
                    li_mat_over_pat.append(mop)
                    li_mat_over_notmat.append(monm)
                    li_pat_over_mat.append(pom)
                    li_pat_over_notpat.append(ponp)
            else:
                df_isec = pd.read_csv(isec_folder+row["index_file"]+isec_suffix, sep="\t", names=isec_cols, dtype=str)

                if "chr" not in str(df_isec.iloc[0]["chr"]):
                    df_isec["chr"] = "chr" + df_isec["chr"].astype(str)

                if df_isec.shape[0] > 10:
                    if setup == "trio":
                        df_isec["snv_occurence"] = df_isec["snv_occurence"].replace(inheritance_dict_trio, regex=True)
                    if setup == "duo_mother":
                        df_isec["snv_occurence"] = df_isec["snv_occurence"].replace(inheritance_dict_duo_mother, regex=True)
                    if setup == "duo_father":
                        df_isec["snv_occurence"] = df_isec["snv_occurence"].replace(inheritance_dict_duo_father, regex=True)
                    

                    for chr in chr_list:
                        df_chr = df_isec[df_isec["chr"] == chr]
                        mat_variants = df_chr[df_chr["snv_occurence"] == "maternal"].shape[0]
                        pat_variants = df_chr[df_chr["snv_occurence"] == "paternal"].shape[0]
                        notmat_variants = df_chr[df_chr["snv_occurence"] == "not maternal"].shape[0]
                        notpat_variants = df_chr[df_chr["snv_occurence"] == "not paternal"].shape[0]

                        if (notmat_variants > 0) and (mat_variants > 0):
                            monm = mat_variants / notmat_variants
                        else:
                            monm = 0
                        if (notpat_variants > 0) and (pat_variants > 0):
                            ponp = pat_variants / notpat_variants
                        else:
                            ponp = 0
                        if (mat_variants>0) and (pat_variants>0):
                            mop = mat_variants / pat_variants
                            pom = pat_variants / mat_variants
                        else:
                            mop = 0
                            pom = 0

                        li_mat_over_pat.append(mop)
                        li_mat_over_notmat.append(monm)
                        li_pat_over_mat.append(pom)
                        li_pat_over_notpat.append(ponp)
                else:
                    for chr in chr_list:
                        li_mat_over_pat.append("")
                        li_mat_over_notmat.append("")
                        li_pat_over_mat.append("")
                        li_pat_over_notpat.append("")

            # merge all results in one table
            df_upd_finder = set_up_results_table()
            df_upd_finder["index_file"] = li_index_file
            df_upd_finder["mother_file"] = li_mother_file
            df_upd_finder["father_file"] = li_father_file
            df_upd_finder["setup"] = li_setup
            df_upd_finder["chr"] = li_chr
            df_upd_finder["perc_roh"] = li_perc_roh
            df_upd_finder["mat_over_pat"] = li_mat_over_pat
            df_upd_finder["mat_over_notmat"] = li_mat_over_notmat
            df_upd_finder["pat_over_mat"] = li_pat_over_mat
            df_upd_finder["pat_over_notpat"] = li_pat_over_notpat

            df_upd_finder = df_upd_finder.drop_duplicates()

            return df_upd_finder

    else:
        print("ERROR: files not found: " + roh_folder+row["index_file"]+roh_suffix + " / " + isec_folder+row["index_file"]+isec_suffix + " / " + vcf_folder+row["index_file"])
        errors += 1

##################################################
# main
#####################
# Setup
# setup final results table
df_upd_finder = set_up_results_table()

# read family_files.csv
df_family_files = pd.read_csv(family_files, sep=",")
df_family_files = df_family_files.fillna("")

n=0
param_list = []
df_result_list = []

#####################
# collect roh and isec files
for idx, row in df_family_files.iterrows():
    #TODO: find better way to create list of rows from df:
    param_list.append(row)
    
    n+=1
    #if n > 100: break
print("starting data collection...")

with Pool(46) as p:
    df_result_list = p.map(collect_roh_inh, param_list)

df_upd_finder = pd.concat(df_result_list, axis=0)

#####################
# Flagging
print("flagging...")

df_cutoffs = pd.DataFrame()
df_cutoffs.at[0,"roh_high_start"] = float(roh_high_cutoff)
df_cutoffs.at[0,"roh_high_end"] = float(10000)
df_cutoffs.at[0,"roh_mixed_start"] = float(roh_high_mixed_start)
df_cutoffs.at[0,"roh_high_mixed_end"] = float(roh_high_mixed_end)
df_cutoffs.at[0,"inh_ratio_high_cutoff_start"] = float(inh_ratio_high_trio_cutoff)
df_cutoffs.at[0,"inh_ratio_high_cutoff_end"] = float(inh_ratio_high_duo_cutoff)
df_cutoffs.at[0,"inh_ratio_high_cutoff_end_duos"] = float(10000)


samples = set(list(df_upd_finder["index_file"]))

li_tags = [[]] * 22
li_chr = []
li_sample = []
li_tags = []
li_con= []


c = 0
n=0
for sample in samples:
    if sample:
        df_sample = df_upd_finder[df_upd_finder["index_file"] == sample]
        
        # cansanguinity flags
        consanguin = False
        df_roh_cut = df_sample[df_sample["perc_roh"] >= consanguin_roh_cutoff]
        if df_roh_cut.shape[0] >= consanguin_min_chr_count:
            consanguin = True

        # chr specific flags
        for idx, row in df_sample.iterrows():
            chr_tags = []

            # ROH
            if not pd.isna(row["perc_roh"]):
                if row["perc_roh"] >= roh_high_cutoff:
                    chr_tags.append(roh_high_tag)
                if (row["perc_roh"] >= roh_high_mixed_start) and (row["perc_roh"] <= roh_high_mixed_end):
                    chr_tags.append(roh_high_mixed_tag)
            else:
                chr_tags.append(snv_per_chr_warning)
            # Inheritance
            if row["mat_over_pat"] and row["pat_over_mat"]:
                if (row["mat_over_pat"] >= inh_ratio_high_trio_cutoff) or (row["pat_over_mat"] >= inh_ratio_high_trio_cutoff):
                    chr_tags.append(inh_ratio_high_tag)
            if row["mat_over_notmat"] or row["pat_over_notpat"]:
                if (row["mat_over_notmat"] >= inh_ratio_high_duo_cutoff) or (row["pat_over_notpat"] >= inh_ratio_high_duo_cutoff):
                    chr_tags.append(inh_ratio_high_tag)
            
            if consanguin:
                li_con.append(consanguin_tag)
            else:
                li_con.append(not_consanguin_tag)
                
            li_sample.append(sample)
            li_chr.append(row["chr"])
            li_tags.append(chr_tags)
        n += 1
    c += 1
    #if c > 1000: break

df_tags = pd.DataFrame()
df_tags["index_file"] = li_sample
df_tags["chr"] = li_chr
df_tags["tags"] = li_tags
df_tags["consanguinity"] = li_con

# merge flag df and sample df
df_upd_finder_tagged = pd.merge(df_upd_finder, df_tags, how = "left", on=["index_file", "chr"])
df_upd_finder_tagged.to_excel(out_dir+upd_table_out, index=False)

#####################
# Plotting
# plot roh and ratios for trios and duos
print("plotting...")
# mat over pat:
trait = "mat_over_pat"
df_upd_finder_trio = df_upd_finder_tagged[df_upd_finder_tagged["setup"]=="trio"]
mop_plot = roh_inh_scatter(df_upd_finder_trio, trait, df_cutoffs)
mop_plot.save(out_dir+trait+"_scatter.html")

# pat over mat:
trait = "pat_over_mat"
df_upd_finder_trio = df_upd_finder_tagged[df_upd_finder_tagged["setup"]=="trio"]
pom_plot = roh_inh_scatter(df_upd_finder_trio, trait, df_cutoffs)
pom_plot.save(out_dir+trait+"_scatter.html")

# mat over not mat:
trait = "mat_over_notmat"
df_upd_finder_trio = df_upd_finder_tagged[df_upd_finder_tagged["setup"]=="duo_mother"]
monm_plot = roh_inh_scatter(df_upd_finder_trio, trait, df_cutoffs)
monm_plot.save(out_dir+trait+"_scatter.html")

# pat over not pat:
trait = "pat_over_notpat"
df_upd_finder_trio = df_upd_finder_tagged[df_upd_finder_tagged["setup"]=="duo_father"]
ponp_plot = roh_inh_scatter(df_upd_finder_trio, trait, df_cutoffs)
ponp_plot.save(out_dir+trait+"_scatter.html")

# plot ROHs
roh_plot = inh_plot(df_upd_finder_tagged, df_cutoffs)
roh_plot.save(out_dir+"roh_per_chr_single.html")

print(errors, " files not found")
