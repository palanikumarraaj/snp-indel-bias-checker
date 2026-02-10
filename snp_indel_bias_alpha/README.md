# snp-indel-bias checker
## ALPHA VERSION 1.0

## Why Bias occurs in Variant calling ?

- Issues due to algorithms potential limitations.
- Read depth and coverage promotes false calls in short reads.
- WES + Short read tech + Repeats or Sequence similarity manifest bias.
- Orthologous or Pseudogenes or Longer repeats on genome promotes alignment limiations.
- Issues and Limiation starts from alignment in short reads data that amplified in downstream analysis.

## Impact of Bias in Downstream analysis

Among all the highly cited variant calling algorithms and tools, the issue of false calling or missed calls over a short window region of high polymorphism or repeats or high similar reference is not yet solved fully.

I have tested with GATK Haplotypecaller, DeepVariant and Freebayes to test the overall consistency and possible bias calls.
As expected, the false call rate is not yet nullifed or solved. It is noted that 1 to 2 % of false call found across multiple WES samples, when I tested. Depends on the target region similarity, read depth and coverage, repeats and NGS tech used, the variation can be from 2 to 5% also.

- For Automated pipelines and Submission of public dataset, these bias can affect the productivity or utility.
- Automation in reports, file handlers and resource distribution can have better control if bias checked.
- Enterprise can point out false calls and perform better RnD if noted early and removed accordingly.
- Bulk batch with combination of ML/DL can point out anomalies better and faster. 

## Usage Options
2 scripts have provided for usage and outputs flexibility ; Both has almost same filter options.

**2 key options used are `-w` and `-m` for window size and min-neighbors respectively**

#### AlphaA1
> indel_snp_bias_finder_alphaA1.py

To know better about the options :

```
python3 indel_snp_bias_finder_alphaA1.py --help
````
- Generates a single output file for biased or potential false call regions
- Takes user option `-o` for output file name
#### AlphaB1
> indel_snp_bias_finder_alphaB1.py
- Generates two VCF output files as `_biased.vcf` and `_unbiased.vcf` from the given input file name
- As name suggest, potential `biased or false calls` are filtered in the output file `_biased.vcf`
- Ideally you may expect 1 to 2% of false call rates from your given input VCF files
- Prefix option is required for output file generation, given with option `-p`

To know better about the options :

```
python3 indel_snp_bias_finder_alphaB1.py --help
````
