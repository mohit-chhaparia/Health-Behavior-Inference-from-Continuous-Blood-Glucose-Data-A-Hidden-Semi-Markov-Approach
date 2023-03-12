# Health Behavior Inference from Continuous Blood Glucose Data: A Hidden Semi Markov Approach

#### Background
Diabetes is a condition when the body doesn't produce enough insulin or fails to use it as efficiently as it should. As per American Diabetes Association, in 2019 about 12.84\% Americans were children and adolescents who had type-I diabetes. For patients with diabetes, hypoglycemia is a condition in which blood sugar (glucose) is lower than the standard range whereas hyperglycemia is a condition in which blood sugar (glucose) is higher than the standard range. Hypoglycemic events can lead to serious life threatening consequences whereas hyperglycemic events can lead to slow and permanent damage to internal organs for patients with type I diabetes.
#### Objective
This research is aimed to develop a model to predict the probability of hypoglycemia in the next 1 hour at each 5 minute intervals. The model is expected to have accuracy comparable to the ML models, better interpretability, and ability to forecast events like hypoglycemia, hyperglycemia, glucose values, etc. 
#### Methods
The research implements the Hidden semi-Markov model with the help of the R package mhsmm, and custom user defined distributions and applies a Monte Carlo approach for forecasting. 
#### Results
Patient-specific and Population-level models are developed and the results are explained by comparing the predicted probability of hypoglycemia with the observed glucose values. For a specific threshold on the population-level model, the sensitivity, and specificity for 30 minute ahead forecast are 93.03\% and 72.50\% and for 60 minute ahead forecast are 89.76\% and 65.27\%, respectively. The 30 minute ahead forecast and 60 minute ahead forecast ROC-AUC for the population level model are 0.9035 and 0.8214, respectively. In literature the 30 minute ahead prediction sensitivity, specificity, and ROC-AUC are generally in the range 74\%-95\%, 79\%-96\%, and 0.73-0.93, respectively.
#### Conclusions
For hypoglycemia prediction, the HSMM model provides better explainability of a patient's physiological latent states compared to the ML models and comparable sensitivity and specificity. The prediction accuracy can be further improved by introducing other parameters like carbohydrates and insulin as covariates directly into the model.
