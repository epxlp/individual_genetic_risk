# Going from association to prediction
During unit 2 you identified several genetic variants associated with disease or continuous traits.
In this session we will see if we can use these findings to predict whether an individual will go on to suffer from disease (or their likely value of a continuous trait).
Here are the results for a GWAS of BMI:


| SNP	| CHR	| BP	| A1	| BETA	| P	| MAF | var_exp |
| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	|
| rs12748679	| 1	| 72644585	| C	| -0.598	| 1.3E-10	| 0.205	| 0.0038 |
| rs630372	| 1	| 177885762	| A	| 0.687	| 5.5E-15	| 0.236	| 0.0055 |
| rs2867125	| 2	| 622827	| T	| -0.606	| 1.1E-09	| 0.168	| 0.0033 |
| rs13130484	| 4	| 45175691	| T	| 0.482	| 2.2E-10	| 0.434	| 0.0037 |
| rs2920930	| 8	| 76733973	| T	| -0.470	| 1.1E-08	| 0.296	| 0.0030 |
| rs8050136	| 16	| 53816275	| A	| 0.594	| 9.5E-15	| 0.395 | 0.0055 |
| rs12970134	| 18	| 57884750	| A	| 0.950	| 4.0E-29	| 0.267	| 0.0114 |


<br><br>
> **Task: Examine the strength of the evidence for the association (P-value), the magnitude of the association (beta) and the variance in BMI explained for each of these SNPs.  
How convinced are you that these variants are associated with BMI?**



<br><br>
> **Question: These variants are all associated with BMI, but what determines how good a variant is for prediction?**

<br><br>
Now imagine we also have identified another variant (VAR1) associated with BMI. The results for this variant are shown in the table below, alongside one of the variants from the table above.



Compare the r<sup>2</sup>, beta and MAF between the following two variants:

| SNP	| CHR	| BP	| A1	| BETA	| P	| MAF | var_exp |
| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	|
| rs12970134	| 18	| 57884750	| A	| 0.950	| 4.0E-29	| 0.267	| 0.0114 |
| VAR1	| 18	| 57884751	| G	| 22.546	| 1.0E-19	| 0.0003	| 0.0100 |

<br><br>
> **Question: Which is the better predictor of BMI in the population?**
<br><br>
> **Question: If someone carries the risk allele for the first variant - how informative is this for this individual?**
<br><br>
> **Question: If someone carries the risk allele for the second variant - how informative is this for this individual?**
<br><br>
> **Question: If someone does not carry the risk allele for the second variant - how informative is this for this individual?**

# Combining SNPs into a Polygenic Risk Score (PRS)

First we are going to generate a file that details the SNPs we want to use in the score. 

Use the following bash command to generate a snps_for_score.txt file.

```
echo 'rs571312        A       1
rs630372        A       1
rs8050136       A       1
rs12748679      T       1
rs13130484      T       1
rs2867125       C       1
rs2920930       G       1' >  ~/Documents/ibsc_unit2/data/snps_for_score.txt
```

The third column in this file denotes the weight that should be applied to each SNP

<br><br>
> **Question: Is this a weighted or an unweighted score?**

<br><br>
> **Question: What values would you include in column 3 to make the other type of score?**


Use this Plink command to generate a risk score for each person in the sample:

```
plink --bfile ~/Documents/ibsc_unit2/data/geno_qc --score ~/Documents/ibsc_unit2/data/snps_for_score.txt --out ~/Documents/ibsc_unit2/data/BMI_score
```

Now start up R and then run the following command to load the data into R:

```
score <- read.table(“~/Documents/ibsc_unit2/data/BMI_score.profile”, header=T)
phen <- read.table(“~/Documents/ibsc_unit2/data/phen_clean.txt”)
```
Run the following commands to investigate the distribution of the genetic risk score
```
BMI_score <- score$CNT2
hist(BMI_score)
hist(BMI_score[phen$BMIcat=="underweight" | phen$BMIcat=="healthy" | phen$BMIcat=="overweight"], col=rgb(1,0,0,0.5), xlim=c(0,14), main = "BMI categories", xlab="genetic_score")
hist(BMI_score[phen$BMIcat=="obese"], col=rgb(0,0,1,0.5), add=T)
```

<br><br>
> **Question: Is there a difference in the PRS between obese and non-obese?**


# Population-level predictions

What is the difference in obesity risk between those in the top 5% and bottom 5% of the score?
```
quantile(BMI_score, c(.05, .95))
sum(phen$BMIcat=="obese" & BMI_score<5)/sum(BMI_score<5)
sum(phen$BMIcat=="obese" & BMI_score>10)/sum(BMI_score>10)
```

There seems to be quite a large difference in risk. Those in the bottom 5% of the risk score have only a 30% risk of being obese, whilst those in the top 5% of the risk score have a 60% risk of being obese.

<br><br>
> **Question: Would an intervention targeted at children in the top 5% of the risk score be more effective than one used for all children?**
<br><br>
> **Question: Would an intervention targeted at children in the top 5% of the risk score be more effective than one used for a random 5% of children?**


You can calculate the sensitivity & specificity for this ‘top 5%’ threshold using the following commands:

Sensitivity (or True Positive Rate) 	= True Positives / Positives 
=  Estimated Obese / Obese

Specificity (or True Negative Rate) = True Negatives / Negatives
=  Estimated Not-obese / Not-obese

```
obese <- ifelse(phen$BMIcat==”obese”, 1, 0)
pred_obese <- ifelse(BMI_score>10, 1, 0)
table(obese, pred_obese)

P <- length(which(obese==1))
N <- length(which(obese==0))
TP <- length(which(obese==1 & pred_obese==1))
TN <- length(which(obese==0 & pred_obese==0))

sens <- TP/P
sens
spec <- TN/N
spec

length(which(pred_obese==1))
```
<br><br>
>**Question: How sensitive and specific is this threshold? How many individuals would be selected? Comment on these values.**

Now try a different threshold for predicting obesity. Try the mean score of obese individuals (7.54).
Adapt the code above to answer the same questions:

<br><br>
>**Question: How sensitive and specific is this threshold? How many individuals would be selected? Comment on these values.**




You can generate a ROC curve and calculate AUC using the following commands:

```
install.packages("pROC")
library(pROC)
plot(roc(obese, BMI_score), print.auc=TRUE)
```

<br><br>
>**Question: Do you think an intervention for obesity should be targeted at children in the top 5% of the risk score? What factors affect your decision?**

<br><br>
> **If there was an intervention that was either very costly or had objectionable side effects or consequences, how would that impact what threshold you might use to classify people as 'at risk'**


# Individual level prediction

There’s clearly a difference in the genetic_score distribution between healthy and underweight versus obese individuals. However, if you knew someone had a genetic score of 4, how useful would that be?
```
table (BMI_score, phen$BMIcat)
summary(phen$BMI[BMI_score==4])
boxplot(phen$BMI~BMI_score)
```

Of 335 people in this category: 7% are underweight, 27% are healthy, 35% are overweight, 32% are obese. So you might say that they are quite likely to be overweight or obese, but the range of BMIs for people in this category is very wide (11.5 to 43.3), so it would be difficult to predict their BMI with any accuracy.


Compare the mean BMIs of people at the top and bottom of the central 90%
```
quantile(BMI_score, c(.05, .95))
summary(phen$BMI[BMI_score==5])
summary(phen$BMI[BMI_score==10])
```

90% of people fall between BMI_scores of 5 and 10, the mean BMI difference between the top and bottom of this range is: 2.81kg/m2 BMI=30 versus BMI=27. 

<br><br>
> **For the majority of people knowing your BMI score is very useful / somewhat useful / not very useful (delete as applicable)**

Another way ti think about the sensitivity and specificity values that we calculate earlier is to work out how often we would be right or wrong if we used BMIscore>10 as a prediction for obesity.

```
obese <- ifelse(phen$BMIcat=="obese", 1, 0)
pred_obese <- ifelse(BMI_score>10, 1, 0)
table(obese, pred_obese)
```
<br><br>
> **Give the numbers in each category:  
given low risk prediction & did not become obese:  
given low risk prediction & did become obese:  
given high risk prediction & did become obese:  
given high risk prediction & did not become obese: 
<br><br> 
Question: what proportion of the time would you be right?**


What about if you took BMI_score>=13 as the cut-off?

```
pred_obese <- ifelse(BMI_score>13, 1, 0)
table(obese, pred_obese)
```
<br><br>
> **Give the numbers in each category:  
given low risk prediction & did not become obese:  
given low risk prediction & did become obese:  
given high risk prediction & did become obese:  
given high risk prediction & did not become obese:**

<br><br>
In this dataset for the 2 people you predicted would be obese, you would be right, but you would have incorrectly told 3524 people that they were not in the high risk group.

<br><br>
> **Question: What proportion of people have a genetic risk score>=13?**

This is a similar situation to the rare highly penetrant mutation you considered at the beginning of this session.

<br><br>
> **Sum up in your own words how useful calculating a polygenic risk score for BMI is. Consider what would affect your assessment?**
<br><br>


# Published example

Here’s some results from a recent publication on the genetics of male-patterned baldness. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5308812/


Other studies you might want to read that include genetic prediction of complex traits:
<br><br>


<br><br>
> **Question: Can you say anything about the genetic architecture of male-patterned baldness from the table above? What implication does this have for genetic prediction?**

They don’t report the variance explained for the genome-wide significant SNPs, but do estimate the total variance explained by all common SNPs to be ~47%.

<br><br>
> **Question: What is important about this type of estimate?**

Speliotes 2010, Nature Genetics 42(11):937-948: predict of BMI from genetic data
Morandi PloS ONE 2012, 7(11):e49919: compared genetic prediction to prediction using traditional risk factors
Timmers BioRxiv 2018, doi.org/10.1101/363036: Peter Joshi's work on genetic prediction of lifespan
