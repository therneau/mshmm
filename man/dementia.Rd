\name{dementia}
\alias{dementia}
\docType{data}
\title{
Cognitive change over time.
}
\description{
  This is a manufactured data set based on a subset of subjects
  from the Mayo Clinic Study of Aging.  It replicates many important
  aspects of that data such as age distribution, covariate distrbutions
  and correlations, and number of visits.
  }
\usage{data("dementia")}
\format{
  A data frame with 24983 observations on 4200 subjects,
  with the following 10 variables.
  \describe{
    \item{\code{ptnum}}{patient identifier of 1, 2, \ldots (not in
      order of enrollment)}
    \item{\code{age}}{age at this observation}
    \item{\code{male}}{0= female, 1= male}
    \item{\code{educ}}{years of education, less than 8 coded as 8}
    \item{\code{apoepos}}{1= one or more e4 alleles in the
      apolipoprotein E (APOE) gene, 0= no e4 alleles}
    \item{\code{cmc}}{cardiovascular and metabolic conditions (CMC)
      score, the sum of seven conditions proposed by the US Department of
      Heath and Human services as indicators of vacular health:
      hypertension, hyperlipidemia, cardiac arrhythmias, coronary artery
      disease, congestive heart failure, diabetes mellitus, and stroke}
     \item{\code{centiloid}}{PET measurment of amyloid level in the brain using
       Pittsburgh compound B (PiB), expressed in centiloid units}
     \item{\code{tau}}{PET measurement of tau level in the brain using
       flortaucipir, expressed in SUVR units}
     \item{\code{p.tau181}}{assay value of phosphorolated tau-181}
     \item{\code{dx}}{diagnosis at this visit, one of CU, MCI, ND or dementia;
       a factor variable; cognitively unimpaired, mild cognitive
       impairment, not demented (CU or MCI) and dementia}
     \item{\code{memory}}{a scaled memory score, based on a battery of
       cognitive tests}
     \item{\code{studydx}}{diagnosis at a study visit, one of CU
       (cognitively unimpaired), MCI (mild cognitive impairment) or dementia}
       a factor variable; cognitively unimpaired, mild cognitive
       \item{\code{offdx}}{diagnosis from a off-study office visit, one of
	 ND (not demented) or dementia}
     \item{\code{death}}{0=alive, 1=dead}
  }
}
\details{
The data set is modeled on a subset of the MCSA data used in Jack (2022).
Covariates were further modified by the addition of random noise to the
visit intervals, sex, and all continuous variables.
As such this data will not reproduce the results in the reference,
nor should it be used for clinical descison making.  It is useful for
illustration of statistical methods, however.

The Mayo Clinic study of aging is based on an age and sex stratified
random sample of subjects from Olmsted County, Minnesota.
Enrolled subjects are scheduled for return visits every 15 months.
This subset was based on those who were enrolled at age 50 or greater,
not demented at enrollment, who had at least one post-enrollment visit.
Patients who cease active particpation are followed passively through
their medical records, this further followup is available on all
subjects through the infrastructure provided by the Rochester
Epidemiology project; thus with respect to the death and dementia endpoints
informative censoring is not be an issue. For visits that are not in
person, which includes phone or spouse interview and off study record
review, the cognitively unimpaired (CU) vs. mild cognitive impairment
(MCI) states cannot be distinguished, but non-demented (ND) vs. dementia can
be ascertained.

Participation in the imaging portion of the study (MRI, amyloid PET,
tau PET, optional FDG PET) requires a second patient visit of
2-3 days, partly as a consequence of this less than half the MCSA
subjects chose to partcicpate in the imaging portion of the study.
A small further number are excluded due to technical issues, e.g., a
metallic implant.
Cumulative radiation exposure limits the total number of PET scans;
for this reason subjects with no current evidence of amyloid accumulation
were often scheduled for rescan at every other study visit rather than
at each return.
Since the PiB ligand did not become available until year 6 of the study
a subset of those who elected to participate in imagine will be missing
PET for their first visits.
For further details see Jack (2022) and the references therein.

The data set is laid out in \emph{timeline} form.
Each row of data contains an identifier, a time, and covariate values
that were measured at that time point, with missing values for those that
were not. Rows where death=1, for
instance, will normally be missing all the continuous variables.
The exception is fixed covariates
such as sex and APOE genetic status, which are copied forward across all
data rows. (Although the variable will not have been formally re-measured
at each visit its value at that time is nevertheless known.)
The seven carido-metabolic conditions are coded as
"any history of" and thus the CMC count can only increase.
}
\references{
 Clifford R. Jack Jr., Terry M. Therneau, et.al (2022)
    Long-term associations between amyloid positron emission tomography,
         sex, apolipoprotein E and incident dementia and mortality among 
         individuals without dementia: hazard ratios and absolute risk,
    Brain Communications, doi:10.1093/braincomms/fcac017
}
\examples{
data(dementia)
table(table(dementia$ptnum))  # number of rows per subject
}
\keyword{datasets}
