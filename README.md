# BatchUPDetection

workflow to detect UPDs in large cohorts

## Requirements

Unix packages:

- nextflow
- bcftools

python packages:

- pandas
- altair
- cyvcf2
- multiprocess
- openpyxl

## Usage

1. Install all required packages and python libraries.
2. Collect all vcf-files to analyse in a specific folder (=vcf_folder)
3. Create a family file according to family_file_template.csv
4. Run nextflow pipeline as:

```
nextflow run scripts/batchUPdetection.nf \
--family_files path/to/family_file.csv \
--vcf_folder path/to/vcf_folder/ \
--outdir_roh path/to/roh_output/ \
--outdir_isec path/to/isec_output/ \
--outdir_results path/to/result_output/
```

Folders roh_output, isec_output and result_output will be used to write corresponding result files

## Results

The workflow will create the following files:

- ROH files: output of bcftools roh
- isec files: output of bcftools isec
- a series of html files with interactive plots:
  - ...scatter.html: these plots show the ratio of maternal over paternal over maternal, maternal over non-maternal and paternal over non-paternal variants per chromosome vs ROH coverage. These plots can be used to identify UPDs.
- upd_finder_overview:
  - this sheet contains all chromosomes of all samples with:

    - file names (index, mother, father)
    - setup (single, duo, trio)
    - chromosome
    - perc_roh: fraction of the chromosome that is covered by ROHs
    - mat_over_pat: ratio of maternal over paternal variants
    - mat_over_notmat: ratio of maternal over non-maternal variants (for duos)
    - pat_over_mat: ratio of paternal over maternal variants
    - pat_over_notpat: ratio of paternal over non-paternal variants (for duos)
    - tags: all tags assigned to the respective chromosome
    - consanguinity: wether or not a consanguinity is suspected

### Tagging and cutoffs

In order to identify UPDs and consanguinity, several tags with a variety of cutoffs are used:


| Flag                   | Cut-Offs                                    | Status/Interpretation                                          |
| ---------------------- | ------------------------------------------- | -------------------------------------------------------------- |
| consanguinity unlikely | <3 chromosomes with >10% ROHs               | -                                                              |
| consanguinity likely   | >=3 chromosomes >10% ROHs                   | handle potential UPD flags with extra care                     |
| roh_high               | > 70% covered by ROH                        | potential Isodisomy or deletions                               |
| roh_high_mixed         | 20-70% covered by ROH                       | potential Isodisomy or mixed UPD                               |
| inh_ratio_high         | in duos >2; in trios > 5; inheritance ratio | potential heterodisomy                                         |
| insufficient_snv       | <200 SNVs/chr                               | number of SNVs is insufficient to reliably detect UPD features |
