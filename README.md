# Multi state hazard models (MSH), for interval censored data
Multi-state hazard models (MSH) have been a key part of the R survival 
package, and play a significant role in our group's medical research work.
This is a library extends these methods to interval censored data. One subset
or extension of interval censored MSH are hidden Markov models (HMM), this the
root of our combined acronym for the package name: MSH + HMM = mshmm.

Three things have spurred the development of the package.
The primary one has been the Mayo Clinic study of aging (MCSA), 
and the work therin to understand the patterns and 
components of cognitive decline, including the impacts thereon of amyloid and
tau depostion (an important driver of Alzheimer's disease, a dementia subtype), 
of accumulated cerbral vascular compromise (itself multifactorial), neuron and
neural function loss,clinical risk factors such as age, sex, hypertension, etc.,
and measurement of the cognitive decline itself.
This dictates many of types of model we want to fit, in particular that direct
measurement of many of these underlying processes is only available during
postmortem pathologic examination of the brain,
so we must rely on indirect ones. (Unlike a potential breast or prostate
cancer, say, a needle biospy of the living brain is not a diagnostic option.)
A second is that our more complex models have required more sophisticated 
numerical methods, e.g., Levenberg-Marquardt and trust region additions to the
iteration, imposition of constraints, and for some, an MCMC approach. 
The optim routine underlying the msm package was simply not up to the task.
The third was to make these often quite complex models easier for the user to
specify by extending user interface ideas worked out for the semi-parametric
MSH methods found in the survival package, and in the process make the
similarities of the two MSH approaches clear.

This current version is intended for sharing within a small group of developers,
until we have worked all the test cases.


