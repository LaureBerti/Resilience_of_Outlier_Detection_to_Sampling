# Resilience of Outlier Methods to Sampling

> We  introduce the notion of  *resilience  to sampling* for outlier detection methods. Orthogonal to  traditional quality performance metrics such as precision/recall, resilience represents the extent up to which the number of outliers detected by a method applied to the whole dataset is preserved when the method is actually applied to a set of samples from a given sampling scheme. 
        
> We propose a novel approach for estimating the resilience to sampling of both individual outlier methods and their ensembles. 

> We performed an extensive experimental study on synthetic and real-world datasets where we study seven diverse and representative outlier detection methods, compare their results obtained from samples versus the ones obtained from the whole datasets and evaluate the accuracy of our resilience estimates. 

> We observed that the methods are not equally resilient to a given sampling scheme and it is often the case that careful joint selection of both the sampling scheme and the outlier detection method is necessary. Our findings have practical impact for detecting outliers from Big Data and we hope that this will initiate research on designing new outlier detection algorithms that are resilient to sampling.

# For More detail

Please find the description of our approach in [our paper on Axiv:](https://arxiv.org/abs/1907.13276)
* Laure Berti-Équille, Ji Meng Loh, Saravanan Thirumuruganathan: Are Outlier Detection Methods Resilient to Sampling? CoRR abs/1907.13276 (2019)


***

# Projet Structure
        
The project is organized into 3 folders:
* *Real-worldData* folder: Eight real-world datasets have been used in the study with the ground truth (when available); they have been originally extracted from UCI (https://archive.ics.uci.edu/ml/datasets.html),  ODDS (http://odds.cs.stonybrook.edu/), and from http://www.dbs.ifi.lmu.de/research/outlier-evaluation/ and cleaned and they are available in the *RealWorldData* folder with outlier detection as 0 or 1 flags in the last columns for each outlier detection methods considered in the study;
* *SyntheticData* folder includes synthetic datasets with 1,000, 5,000, and 10,000 records that have been generated from independent bivariate normal distribution with mean (0, 0) and standard deviations 1 and 2 (with the respective R script) and they can be found as two separate archives.
* *Code* folder includes the Python script for outlier ensembling and 4 R scripts for generating the synthetic data with controlled distribution of outliers, for sampling (random, blocking and partitioning) and for outlier detection.


# Reproducing Results

After downloading the datasets, you make sure to change the path to access (read/write) the input/output datasets and you can re-run the results presented in our paper simply by executing corresponding programs. 

>  For the dataset *Arrhythmia* only (due to space limitation): we share the files obtained after sampling and re-detection of outliers in three separate folders: random, blocking and partitioning. the last columns of each file indicates the outlier detection flag from the sample for each method (prefixed with "s"). The results for the other datasets can be obtained by running the R script for sampling and re-running the outlier detection over the samples and averaging the results of the 100 generated experimental runs.

Our experiments are run in Python 2.7 on a MacOS Sierra with Intel 2.8 GHz Core i7 processor processor and 16 GB memory. 
