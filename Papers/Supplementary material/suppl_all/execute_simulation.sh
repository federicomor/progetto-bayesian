#!/usr/bin/env Rscript                                                                                                                                                                                                                                                                                                                                                                                                                                         
jobnum <- 1
for ( datatype  in 1:2 ) {
  for ( stdev in c(0.5, 1.0) ) {
    for ( correlation in c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99) ) {
      for( Tm in c(5, 10) ) {
#        cat(sprintf("nohup ./SimStudyJCGSrev.R 100 %s 100 %s %s %s > out%s.txt & \n",datatype, Tm, correlation, stdev, jobnum))
         cat(sprintf("./SimStudyJCGSrev.R 100 %s 100 %s %s %s  \n",datatype, Tm, correlation, stdev))                                                                                                                                                                                                                                                                                                                                                         
         jobnum  <- jobnum + 1
      }
    }
  }
}
