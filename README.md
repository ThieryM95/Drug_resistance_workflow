# Drug_resistance_workflow
This is the first release of the analysis code for:

The influence of biological, epidemiological, and treatment factors on the establishment and spread of drug-resistant Plasmodium falciparum Thiery Masserey1,2, Tamsin Lee1,2, Monica Golumbeanu1,2, Andrew J Shattock1,2, Sherrie L Kelly1,2, Ian M Hastings3, Melissa A Penny1,2*

1Swiss Tropical and Public Health Institute, Basel, Switzerland

2 University of Basel, Basel, Switzerland

3 Liverpool School of Tropical Medicine, Liverpool, UK

*corresponding author: melissa.penny@unibas.ch

In the study mentioned above, we have developed a disease modelling approach with emulator-based global sensitivity analysis to systematically quantify which factors drive the establishment and spread of drug-resistant Plasmodium falciparum parasites. This analysis was performed using an individual-based malaria model (https://github.com/SwissTPH/openmalaria/wiki).

The folder "WF_spread" contains the code used to assess via global sensitivity analysis the impact of factors on the rate of spread parasites having different degrees of resistance to artemisinin derivatives or their partner drugs when used in monotherapy and combination (as ACTs). The analysis need the following steps:(i) randomly sampling combinations of parameters of the factors of interest; (ii) simulating and estimating the rate of spread of the resistant genotype for each parameter combination in OpenMalaria; (iii) training an HGP to learn the relationship between the input and output with iterative improvements to fitting through adaptive sampling, and (iv) performing a global sensitivity analysis based on the Sobol variance decomposition. The code asks to specify the analysis details in the script "Option.R" (such as which treatment profile is modelled, how many settings are modelled, etc.). The script "launch.R" calls all the necessary functions to perform all the analysis steps specified in the script "Option.R".

The folder "wf_establishment" contains the code that aimed to assess the probability of establishment of mutations with known selection coefficients (rate of spread). The probability of establishment is assessed by simulating the emergence of the mutations conferring drug resistance in our model and monitoring if they esthablished or not. As above, the code asks to specify the analysis details in the script "Option.R" (such as which treatment profile is modelled, how many settings are modelled, etc.). Then, the analysis can be performed by running the script "launch.R". The script "launch.R" call all the necessary functions to perform the analysis specified in the script "Option.R".

Note that the plotting code uses the folder structure and working directories used by the researchers. The file paths will have to be adjusted to replicate figures.
