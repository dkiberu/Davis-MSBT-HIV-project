# Simulation of clonal expansion and latent HIV reservoir size prediction using machine learning in virally suppressed individuals in Uganda
This repository includes scripts used in mathematical modelling of the latent viral reservoir (LVR) and machine learning prediction of the LVR size using clinical data from a cohort of Ugandans living with HIV.\
The LVR is made up of resting cells infected by HIV in such a way that viral particles are not released and as such the immune system cannot target them. Anti-retroviral drugs (ARVs) are also unable to interrupt the LVR and upon ceasing therapy, the reservoir can be reactivated leading to rebound HIV infection. During effective therapy, the reservoir is maintained through clonal expansion, which can be homeostatic in order to preserve the number of resting cells or antigen induced so as to combat infection. Given the multiple endemic diseases in countries such as Uganda, immune activation and clonal expansion of the LVR should play a major role in the dynamics of this reservoir.\ 
Our results from modeling this reservoir showed that most individuals were having an increase in the latent viral reservoir size which could imply that frequent immune challenge may negatively impact the dynamics of the reservoir causing it to increase in size thus hindering clearance of HIV infection. For individuals with a decreasing reservoir, the dynamics were similar to those reported in countries with less endemic diseases, implying that immune activation, if it exists, may have minimal effect on the dynamics of the latent HIV reservoir.\
We also explored the potential of machine learning in predicting the size of the LVR since measurement of this reservoir is expensive, time consuming and requires expertise. We used clinical data collected over time including; CD4 count, CD8 count, gender, time on anti-retroviral therapy (ART), pol and gp41 phylogenetic distances, to predict the LVR size in terms of infectious units per million CD4+ T cells (IUPM). Our findings showed that machine learning could indeed be used to predict the LVR size using other clinical measurements though larger datasets would be needed and other metrics such as how early individuals initiated ART after infection should be included for more robust predictions. The best performing algorithm was the light gradient boost model (LGBM). 
