There are many ways in which you can estimate an individual’s risk of a particular disease or predict their likely trait value. How you do this will depend on the trait, the circumstances of the individual and the reason for wanting to make such a prediction.

In today’s session we will demonstrate estimating risk and predicting trait values in a variety of settings using different types of information.


# The data
Because we have our genes from birth, we could use an individual’s genetic data to estimate their risk having a certain disease later in life. In practice you might be more interested in doing this for cancer or heart disease, but we’ll do this for a trait that you have data for: **BMI/obesity**.

# Recapping association results
In the GWAS exercise in unit 2 you found a number of genetic variants associated with BMI. 
Let’s start by reminded ourselves of the association evidence between being BMI and the associated SNPs.


| CHR	| SNP	| A1	| BETA	| A1_freq	| var_exp	| OR	| STAT	| P |
| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| --- |
| 1	| rs12748679	| C	| -0.5972	| 	| 0.0038	| 	| -6.42	| 1.4e-10 |
| 1	| rs630372	| A	| 0.6878	| 	| 	| 	| 7.831	| 5.4e-15 |
| 2	| rs2867125	| T	| -0.6013	| 	| 	| 	| -6.04	| 2.2e-8 |
| 4	| rs13130484	| T	| 0.4853	| 	| 	| 	| 6.391	| 1.7e-10 |
| 8	| rs2920930	| T	| -0.4611	| 	| 	| 	| -5.604	| 2.2e-8 |
| 16	| rs8050136	| A	| 0.5906	| 	| 	| 	| 7.702	| 1.5e-14 |
| 18	| rs571312	| A	| 0.9928	| 	| 	| 	| 11.26	| 3.3e-29 |  

<br><br>
> **Task: Examine the strength of the evidence for the association (P-value), the magnitude of the association (beta/OR) and the variance in BMI explained for each of these SNPs, that you estimated in unit 2.  
How convinced are you that these variants are associated with BMI/being overweight?**

<br><br>
Some of these are very strong associations. p=3e-29 is extremely strong evidence of an association. Despite this the effect sizes and the variance explained are quite small {elaborate}. In unit 2 you estimated that the total variance explained for these 7 SNPs was 4%. 

<br><br>
> **Question: These variants are all associated with being overweight, but what do you think determines how good a variant is for prediction?**

<br><br>
**BUT – these only tell us how useful a genotype is ‘on average’ in the population. How useful a genotype is for prediction for an individual also depends on the individual’s genotype.
Some variants might have very large effect, but be extraordinarily rare. This variant would have a low variance explained, but for someone who had this genotype, this would be very informative. Discuss this with a neighbour**

# Generating a genetic risk score

Use these commands to extract the 7 SNPs and generate a risk score for each person in the sample:

```
plink --bfile geno_raw --extract snps_for_score.txt –recode A --out BMI_snps
plink --bfile geno_raw --score snps_for_score.txt --out BMI_score
```

Now we will load the data into R

```
snps <- read.table("~/ibsc_unit2/data/BMI_snps.raw”, header=T)
score <- read.table(“~/ibsc_unit2/data/BMI_score.profile”, header=T)
phen <- read.table(“~/ibsc_unit2/data/phen_clean.txt”)
```
Run the following commands to investigate the distribution of the genetic risk score
```
BMI_score <- score$CNT2
hist(BMI_score)
hist(BMI_score[phen$BMIcat=="underweight" | phen$BMIcat=="healthy"], col=rgb(1,0,0,0.5), xlim=c(0,14), ylim=c(0,600), main = "BMI categories", xlab="genetic_score")
hist(BMI_score[phen$BMIcat=="obese"], col=rgb(0,0,1,0.5), add=T)
```

# Population-level predictions

What is the difference in obesity risk between those in the top 5% and bottom 5% of the score?
```
quantile(BMI_score, c(.05, .95))
sum(phen$BMIcat=="obese" & BMI_score<5)/sum(BMI_score<5)
sum(phen$BMIcat=="obese" & BMI_score>10)/sum(BMI_score>10)
```

There seems to be quite a significant difference in risk. Those in the bottom 5% of the risk score have only a 30% risk of being obese, whilst those in the top 5% of the risk score have a 60% risk of being obese.

<br><br>
> **Question: Do you think it would be an economically effective strategy to target an intervention only at children that were in the top 5% of the risk score?  
Although this is unlikely to happen for an outcome such as BMI, think about if the outcome were cancer or Alzheimer’s disease.  
Discuss this with a neighbour**

<br><br>
**As a genetic risk score increases there is a sliding scale of increased risk, but any health policies would likely have to define a threshold and this would have sensitivity and specificity implications.**

# Individual level prediction

There’s clearly a difference in the genetic_score distribution between healthy and underweight versus obese individuals. However, if I told you someone had a genetic score of 4, what would that tell you?

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
> **For the majority of people knowing your BMI score is very useful / somewhat useful / pretty useless (delete as applicable)**

If we classed everyone with BMI_score >10 as at risk of obesity. How often would we be right and wrong?

```
obese <- phen$BMIcat
levels(obese) <- c("control", "obese", "control", "control")
risk <- factor(BMI_score)
levels(risk) <- c("<=10", "<=10", "<=10", "<=10", "<=10", "<=10", "<=10", "<=10", "<=10", ">10", ">10", ">10")
table(obese, risk)
```
<br><br>
> **Give the numbers in each category:  
given low risk prediction & did not become obese:  
given low risk prediction & did become obese:  
given high risk prediction & did become obese:  
given high risk prediction & did not become obese: 
<br><br> 
Question: what proportion of the time would you be right?**

<br><br>
You are only going to be right 57% of the time. So for an individual you can make a prediction, but it’s not much better than flipping a coin

What about if you took BMI_score>=13 as the cut-off?

```
risk <- factor(BMI_score)
levels(risk) <- c("<=10", "<=10", "<=10", "<=10", "<=10", "<=10", "<=10", "<=10", "<=10", "<=10", "<=10", ">10")
table(obese, risk)
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
> **Question: What proportion of people have a genetic risk score>=13?  
Can we consider this a similar situation as the rare high penetrant mutations you were considering in the other part of this practical?**
