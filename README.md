# PRIDE_Code
[ReadMe.txt file]

There are three files (.R) for illustrating our proposed method "PRIDE".

[PRIDE_Function.R file]
The "PRIDE_Function.R" file is the function for applying PRIDE to estimate probabilistic differential edges. 
This function can generate posterior samples for "BETA" and "Gamma" in the PRIDE model, which can be applied to estimate probabilistic differential network
and the relative strength of differential edges. 
The output of PRIDE_Function is a "list", which is the default setting by using the "coda" package to store the posterior samples.

[Example.R file]
The "Example. R" file consists of a toy example to illustrate applying PRIDE to conduct differential network analysis. 
This example will use the comparison between AR(2) network and AR(1) network to illustrate how PRIDE can provide probabilistic inference about differential edges. 

[Performance_Criteria.R file]
The "Performance_Criteria.R" file is the function for calculating the performance of estimating differential network by our proposed methods.  
By utilizing this function, one can generate those commonly used criteria when presenting estimation accuracy, such as F1-score, MCC score, etc.
The true positive (TP) here means: correctly identifying the true differential edge as positive.
