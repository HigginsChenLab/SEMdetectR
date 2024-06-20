## TO-DO
1) remove probes that have stats calculated as NAs if using clustering (it fails otherwise)
2) perform imputation on DNAme and remove probes w/ non-random NAs
3) Report individual SEMs


## Installation

You can install the development version of SEMdetectR from [this folder](https://github.com/yaromar/SEMdetectR) with:

``` r
install.packages("devtools")
devtools::install_github("HigginsChenLab/SEMdetectR")
```


## Usage

Here is a basic example of how to use SEMdetectR to detect Stochastic Epigenetic Mutations in DNA methylation beta values:

```r
library(SEMdetectR)

#Assuming `DNA_methylation_betas` is your input dataframe with rows = samples and columns = probes
results <- detectSEM(DNA_methylation_betas, num_cores=4, cluster=TRUE, rf=FALSE) # detects SEMs with the original IQR-based method, separately for unmethylated, intermediate, and methylated probes
results_rf <- detectSEM(DNA_methylation_betas, num_cores=4, cluster=TRUE, rf=TRUE, array="450k") # detects SEMs with the RF-based method, separately for unmethylated, intermediate, and methylated probes found on Illumina 450k array. Even if the parameter "cluster" is set to FALSE, it will still need to cluster probes, because this information is used by the RF models.
```

## Contributing

We welcome contributions to SEMdetectR! If you have suggestions for improvements or encounter any issues, please open an issue or submit a pull request on GitHub.


## License

SEMdetectR is distributed under the terms of the [MIT License](LICENSE.md).


## Contact

For questions or support, please open an issue on the GitHub repository, or contact Yaroslav Markov at yaroslav.markov@yale.edu


## Citation

We kindly request that you cite the following paper if you use SEMdetectR in your research:

> Markov, Y., Levine, M., & Higgins-Chen, A. T. (2024). Reliable detection of stochastic epigenetic mutations and associations with cardiovascular aging. GeroScience. 10.1007/s11357-024-01191-3

A BibTeX entry for LaTeX users is:

```bibtex
@article{Markov2024,
  title={Reliable detection of stochastic epigenetic mutations and associations with cardiovascular aging},
  author={Markov, Y., Levine, M., & Higgins-Chen, A. T.},
  journal={GeroScience},
  volume={},
  number={},
  pages={},
  year={2024},
  publisher={},
  doi={10.1007/s11357-024-01191-3}
}
```

Additionally, please cite Gentilini et al. (2015) if using the original IQR-based method and Houseman et al. (2012) if using the RF-based method (it uses cell count estimation described in the publication):

> Gentilini, D., Garagnani, P., Pisoni, S., et al. (2015). Stochastic epigenetic mutations (DNA methylation) increase exponentially in human aging and correlate with X chromosome inactivation skewing in females. Aging (Albany NY). 7, 8, 568-78. 10.18632/aging.100792

```bibtex
@article{Gentilini2015,
  title={Stochastic epigenetic mutations (DNA methylation) increase exponentially in human aging and correlate with X chromosome inactivation skewing in females},
  author={Gentilini, D., Garagnani, P., Pisoni, S., et al.},
  journal={Aging (Albany NY)},
  volume={7},
  number={8},
  pages={568-78},
  year={2015},
  publisher={},
  doi={10.18632/aging.100792}
}
```

> Houseman, E.A., Accomando, W.P., Koestler, D.C. et al. (2012). DNA methylation arrays as surrogate measures of cell mixture distribution. BMC Bioinformatics. 13, 86. 10.1186/1471-2105-13-8

```bibtex
@article{Houseman2012,
  title={DNA methylation arrays as surrogate measures of cell mixture distribution},
  author={Houseman, E.A., Accomando, W.P., Koestler, D.C. et al.},
  journal={BMC Bioinformatics},
  volume={13},
  number={86},
  pages={},
  year={2012},
  publisher={},
  doi={10.1186/1471-2105-13-8}
```
