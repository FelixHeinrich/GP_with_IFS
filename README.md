# GP_with_IFS
Genomic prediction with incremental feature selection 

The script determines the optimal number of SNPs to use for genomic prediction by ranking SNPs according to a GWAS and adding them incrementally to a random forest prediction model.

Genotype and phenotype data needs to be in the binary format used by [PLINK](https://www.cog-genomics.org/plink/1.9/formats).

The script will save the results in a table as well as a figure, which show the prediction accuracy when using different numbers of SNPs as input.

## Requirements

R 3.6.3 or higher

R packages data.table, ggplot2 and ranger

[PLINK](https://www.cog-genomics.org/plink/1.9/) 1.9 needs to be installed

## Usage

```
Rscript GenomicPrediction_with_IncrementalFeatureSelection.R plinkBinaryPrefix threadCount
```
## Examples
The datasets we analyzed in our publication can be found as compressed PLINK binary files under [Datasets](https://github.com/FelixHeinrich/GP_with_IFS/tree/main/Datasets) (see also their corresponding publications [https://doi.org/10.3389/fgene.2018.00128](https://doi.org/10.3389/fgene.2018.00128) (Chicken), [https://doi.org/10.1534/g3.111.001453](https://doi.org/10.1534/g3.111.001453) (Pig) and [https://doi.org/10.1534/g3.119.400498](https://doi.org/10.1534/g3.119.400498) (Plants)).

## License

This project is licensed under the **GPL-3.0 License** - see [LICENSE](LICENSE) for more information.
