## Background Information

## Multi-Trait Score
 
We wanna combine several meaningful traits into a single numerical value that reflects the overall performance of each genotype.
This score is based only on traits that show clear and significant! genetic variation.
 
### BLUP = Best Linear Unbiased Prediction
A statistical estimate of the genetic performance of a genotype. BLUPs represent the genetic value of a genotype for a trait!

Corrects for:
- environmental effects
- replicate/block effects
- random error
 
### EBLUP = Estimated BLUP
- This is simply: EBLUP = BLUP + grand mean
- BLUPs alone are centered around zero so they are sucky to read
- EBLUPs shift the predictions back to the original trait scale, making them comparable and interpretable.
- For PCA and multi-trait scoring, EBLUPs are typically used so we are doing that, too.
 

## Traits Included in the Multi-Trait Score
Only traits with a significant genotype effect (p < 0.05)
 
Included traits (BLUP/EBLUP used):
- SQI (Seed Quality Index) this I have calculated using the measured traits from the seed quality sheet they provided [Updated with Sarah's input!]
- CCM (Chlorophyll Content)
- HGW (Hundred Grain Weight)
- Fresh Mass (FM) 
  - "After sampling" was incomplete and I left it out because I didn't see the use of it.
- Dry Weight (DW)
- Plant Height (PH)
- BBCH

## Traits Excluded from the Multi-Trait Score
I excluded them because the genotype had no significant effect on the trait (Flowering Time) or because I was uncertain.
 
Flowering Type (FT): Genotype had no significant effect on flowering type


