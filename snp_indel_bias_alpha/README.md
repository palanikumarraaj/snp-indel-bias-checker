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
