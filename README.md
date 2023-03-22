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
Description of arguments:

```
plinkBinaryPrefix     name of PLINK files without filename extension (.bed, .bim and .fam)
threadCount           number of threads to use
```

## Examples
The datasets we analyzed in our publication can be found as compressed PLINK binary files under [Datasets](https://github.com/FelixHeinrich/GP_with_IFS/tree/main/Datasets) (see also their corresponding publications [https://doi.org/10.3389/fgene.2018.00128](https://doi.org/10.3389/fgene.2018.00128) (Chicken), [https://doi.org/10.1534/g3.111.001453](https://doi.org/10.1534/g3.111.001453) (Pig) and [https://doi.org/10.1534/g3.119.400498](https://doi.org/10.1534/g3.119.400498) (Plants)).

#### Exemplary usage
1. Download the [script](https://github.com/FelixHeinrich/GP_with_IFS/blob/main/GenomicPrediction_with_IncrementalFeatureSelection.R) as well as the [soy_HT](https://github.com/FelixHeinrich/GP_with_IFS/blob/main/Datasets/soy_HT.tar.gz) dataset and put them both in the same directory
2. Open a terminal and go to the directory where the downloaded files are located
3. Unpack the dataset with 
```
tar -xzf soy_HT.tar.gz
```
4. Run the script using 10 threads with  
```
Rscript GenomicPrediction_with_IncrementalFeatureSelection.R soy_HT 10
```

## License

This project is licensed under the **GPL-3.0 License** - see [LICENSE](LICENSE) for more information.
