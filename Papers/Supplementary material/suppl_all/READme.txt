
This zipped folder contains R-scripts and an R-package that, together, carry out
the computation presented in the paper

"Dependent Modeling of Temporal Sequences of Random Partitions."


This READme.txt file contains a brief description of each file.  Note that Each 
of the R-scripts write output to external files.  The lines of R-code that
carry this out will write files to the local directory.  Additionally, code that 
organizes results from simulation studies into tables and figures reads from 
external files the contain model fits and simulation study results.  These are 
provided for your convenience. As you run the code and write to the local directory, 
it will the also read output from the local directory.   


drpm R-package - 

  The drpm R-package contains C code and R-wrappers that permit fitting the 
  temporally dependent random partition models (5)and (9) described in the 
  main article. The spatio-temporal model described in Section 4.2 is also 
  available. In addition, the package contains routines to fitting the wddp 
  and the lddp that are included as competitors in the third simulation study.

  This package can be installed locally by unzipping it and in the command 
  line moving to the directory where the package is found and then typing the 
  following command.
    
  > R CMD INSTALL drpm

  Alternatively, the package can also be found in the first authors' Git hub 
  repository.  This can be installed into R using the following code in an
  R-session.

  > devtools::install_github("https://github.com/gpage2990/drpm.git")


MonteCarloExperiments_laggedARIandCorrelations.R - 

  Contains R-code that carries out the monte carlo simulation that is used
  to produce Figures 1 and 2 in the main document. Also, the code used to 
  carry out the  monte carlo simulation that produced Figure 3 from simulation 
  study 2 that is detailed in Section 3.2 is contained in this file. 



simstudy1_rev.R - 

  This file contains the R-code that runs the 1st simulation study described in 
  Section 3.1 of the main article. It also contains code that organizes results.


simstudy3_rev.R - 

  This file contains the R-code that runs the 3rd simulation study described in 
  Section 3.3 of the main article. Note that all output from the simulation study
  is read to external files.  In addition, plots based on the simulation study
  output are saved to external directories.  Directories that contain simulation 
  study output and those to which plots are saved need to be customized. 


PM10DependentPartitionFullModel_rev - 

  This file fits the 16 models detailed in Sections 4.1 and 4.2 of the main document
  and produces Table 3 and Figures 5 - 8. 


Functions.R - 
  
  This file contains utility functions that are used in the files described above to
  carry out simulation studies and applications


execute_simulation.sh - 
  
  Bash script that permits running the 3rd simulation study parallelized.  Run the
  simulation can be done in parallel by executing the following bash command.  Note
  that GNU Parallel must be installed.  
  
  > ./execute_simulation.sh | parallel 



The last four objects in the JCGS_Codes zipped folder are folders that contain
results from the simulation studies and the .RData object that contains the all
16 models fits associated with the PM10 data.  These are provided to facilitate
running code or corroborating results in the paper.  Specifically each file is 

simstudy1results - 

  Folder that contains all the output from the simulation study in Section 3.1

marginal_correlations - 

  Folder containing all output that displayed in Section 3.2 from the second
  simulation study

simstudy3results - 

  Folder that contains all the output from the simulation study in Section 3.1


PM10_ALL_SpCo_4_rev.RData - 

  .RData object created to save R-session where all 16 models were fit to the 
   PM10 data.    





  