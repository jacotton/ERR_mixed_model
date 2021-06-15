# mixed-effects model to estimate individual level egg-reduction rates

Code to implement a mixed-effects model for egg reduction rates. These scripts fit a poisson GLM to count data taken pre- and post- anthelminthic treatment, 
with pre/post and some other potential confounders as main effects and individual and school-level nested random effects in a Bayesian setting, using MCMCglmm.

Model is from Martin Walker (Martin Walker (mwalker@rvc.ac.uk).

Code by Tom Crellen (thomas.crellen@bdi.ox.ac.uk), and Martin, plus a tiny bit by James Cotton (jc17@sanger.ac.uk).

See:

Martin Walker, Thomas S. Churcher, María-Gloria Basáñez,
Models for measuring anthelmintic drug efficacy for parasitologists,
Trends in Parasitology 30:528-537 (2014).

and:

Walker, M., Mabud, T.S., Olliaro, P.L. et al. 
New approaches to measuring anthelminthic drug efficacy: parasitological responses of childhood schistosome infections to treatment with praziquantel. Parasites Vectors 9, 41 (2016).

for an explanation of the model
This code was originally written and used in:

Crellen T, Walker M, Lamberton PH, et al. 
Reduced Efficacy of Praziquantel Against Schistosoma mansoni Is Associated With Multiple Rounds of Mass Drug Administration. 
Clin Infect Dis. 63(9):1151-1159. (2016).

Please cite this paper, and the two Walker et al. papers above if you use these in published work.
The script called Crellen_CID.trimmed_script.R should recreate this analysis exactly.

The individual ERR results were recalculated (identically up to Markov error)  then also used, and marginal distributions calculated for some different subsets of individuals in:

Berger D, Crellen T, Lamberton PH, et al.
Extensive population genomics diversity with limited selection of Schistosoma mansoni despite repeated mass drug administraiotn. 
Nature Communications (in press, 2021).

The script ERR.R will recreate exactly this analysis.
