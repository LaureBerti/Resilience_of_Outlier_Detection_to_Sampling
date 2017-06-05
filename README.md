# Resilience
Estimating the Resilience to Sampling of Outlier Detection Methods

> We  introduce the notion of  *resilience  to sampling* for outlier detection methods. Orthogonal to  traditional quality performance metrics such as precision/recall, resilience represents the extent up to which the number of outliers detected by a method applied to the whole dataset is preserved when the method is actually applied to a set of samples from a given sampling scheme. 
        
> We propose a novel approach for estimating the resilience to sampling of both individual outlier methods and their ensembles. 

> We performed an extensive experimental study on synthetic and real-world datasets where we study seven diverse and representative outlier detection methods, compare their results obtained from samples versus the ones obtained from the whole datasets and evaluate the accuracy of our resilience estimates. 

> We observed that the methods are not equally resilient to a given sampling scheme and it is often the case that careful joint selection of both the sampling scheme and the outlier detection method is necessary. Our findings have practical impact for detecting outliers from Big Data and we hope that this will initiate research on designing new outlier detection algorithms that are resilient to sampling.

***

# Projet Structure
        
The project is organized into 3 folders:
* Real-world Data sets: Eight real-world datasets have been used in the study with the ground truth (when available); they have been originally extracted from UCI (https://archive.ics.uci.edu/ml/datasets.html),  ODDS (http://odds.cs.stonybrook.edu/), and from http://www.dbs.ifi.lmu.de/research/outlier-evaluation/ and cleaned and they are available in the *RealWorldData* folder with outlier detection as 0 or 1 flags in the last columns for each outlier detection methods considered in the study;
* Synthetic Datasets with 1,000, 5,000, and 10,000 records have been generated from  independent bivariate normal distribution with mean (0, 0) and standard deviations 1 and 2 and they can be found in the *SyntheticData* folder as two separate archives.
* *Code* folder includes Python script for outlier ensembling and R scripts for generating the synthetic data with controlled distribution of outliers, for sampling (random, blocking and partitioning) and for outlier detection.


# Reproducing Results

After downloading the datasets, you make sure to change the path to access (read/write) the dataset and you can re-run the results presented in our paper simply by executing corresponding programs. 

Our experiments are run in Python 2.7 on a MacOS Sierra with Intel 2.8 GHz Core i7 processor processor and 16 GB memory. 
