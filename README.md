# Project A2

Goal: Clustering weekly data of one year of PM10 (plus covariates - see AGRIMONIA project)   
Models & R packages:  
- drpm (Page, Quintana, Dahl (2022) "Dependent Modeling of Temporal Sequences of Random Partitions", JCGS, R-package on Github)    
- ppmSuite (various models implemented, also PPMs) on https://cran.r-project.org/web/packages/ppmSuite/index.html   
Data: old https://zenodo.org/records/7563265 new https://zenodo.org/records/7956006   


# Models order
**Important**: for models fitting let's use `df_weekly_scaled_centered` loaded from `include.R`, which has (now it's done):

- numerical covariates scaled
- spatial coordinates also scaled (fitting appears to be better)
- PM10 just centered, not scaled

And for models evaluation see file `Metrics and workflow.Rmd` (now it's done).
See also that for understand how to create the plots of the clusters.

**Models recap**:

- DRPM (Federico): with time and space
- Gaussian PPMx (Federica): just covariates
- SPPM (Oswaldo): just space
- Curve ppmx (Aby): time and covariates
- Linear (Giulia): time and covariates, simple baseline model

| Model                    | Time     | Space    | Covariates | Clustering|
|--------------------------|----------|----------|------------|-----------|
| DRPM (Federico)          |  ✓       | ✓          |  X      | ✓      | 
| Gaussian PPMx (Federica) |  X       |   X        |  ✓      | ✓      | 
| SPPM (Oswaldo)           |  X       |    ✓       | X      | ✓      |
| Curve PPMx (Aby)         |   ✓     | X           |  ✓     | ✓      |
| CarBayesST               |   ✓     | ✓           |  ✓      |  X    |
| Linear                   |    ✓    |    X        |    ✓    |   X    |


The models without time means that time is not inside the model. Then we can "force" time to be included by looping or doing some other tricks, but time will never be inside the model.

The models without space could consider to use space as a covariate? like the locations longitude and latitude as numerical covariates.

# Report chapters idea
Report here at https://it.overleaf.com/6577451481wwsqvrgxswjy#c9ce5d

1. Data inspection + problem presentation
	- 1.0 Introduction (Ettore)
	- 1.1. NAs treatment, stations descarded (Giulia)
	- 1.2. trends, choice of the year (Federico, but help may be needed)
	- ?
2. Models 
	- 2.1. SPPM, as it is just spatial, like the start point (Oswaldo)
	- 2.2. DRPM, space plus time (Federico)
	- 2.3. Gaussian PPMx, only covariates (Federica) 
	- 2.4. Curve PPMx, covariates and time (Aby)
        - 2.5. Giulia's linear model (Giulia)
 -  2.6  CARBayesST + covariate selection, so Ettore and Giulia work, to prepare for the models with covariates  (only cited as a not determinant attempt) + variable selection, different for every model but summerized in 1 

3. Model choice
	- 3.1. model comparison with the metrics used (Federico)
	- 3.2. why model X is the "best" 
	- CV (?), MSE, ARI  and other metrics used to compare the models (Federico)
4. Analysis of the results
	- 4.1. How do we interpret the model? Why some covariates seem more important than others? (Ettore)
5. Conclusion (Ettore)
6. Further Developments (Ettore)
	- 6.1. model averaging Package‘AICcmodavg’
	- 6.2. weekend-weekday division
	- 6.3. prior dall'anno precedente
7. Appendix
	- plot librairies, animations, html... (Ettore)
	- link (Federico)


## Comments on the models
Write them here so are easily available to everyone, and to be able to build a general frame of the information we can extract.
Write them below each mode model plot.  
Use the visualize html page to see peculiar features of your model, as we are ok interested in common comments, on which all models agree, but also some patterns that a certain model may spot and others models not. So especially use the "Covariates visualization" section in the html page.

1. General comments
	- read page 49 of report/allegato pria ecc for understading the causes of pollution "l'influenza delle condizioni meteorologiche sulle concentrazioni degli inquinanti".

![](./src/figures/DRPM/Time%20Series/plt_modemap.png)

1. Comments:
	- DRPM shows as the other models a clear division in mountain (alpi), pre-mountain (prealpi), milan area, the south-east area (the area of mantova) and the one on the south part (around pavia), but also defines another area between milan and the prealpi, passing through monza and brescia (which maybe is ok as they are both highly industrial areas? brescia for sure).
	- the most differentiating covariates are
		- altitude (as all the models)
		- EM_nh3_agr_soils (like in all models, but with different patterns) toghether with the other antropogenic (nox sum, nh3, ecc)
	- the EM_nox_traffic shows a strange pattern, where the yellow and and red curves swap for intensities, as also Gaussian PPMx model will show. This is because "Nel luglio 2018, la Regione Lombardia ha adottato il Piano Aria Integrato Regionale (PAIR) 2018-2020, che contiene le misure per ridurre le emissioni inquinanti da diverse fonti, tra cui il traffico veicolare, i generatori di calore a biomassa legnosa, il settore industriale e agricolo". Go in report folder and check the file named "Allegato pria ecc"
	- the selection of a single cluster of two units in the south zone (the light blue ones) is done also by all the other models but Gaussian PPMx

![](./src/figures/sPPM/Time%20Series/plt_modemap.png)

1. Comments:
	- ecc
	- ecc


![](./src/figures/Gaussian%20PPMx/Time%20Series/plt_modemap.png)

1. Comments:
	- Altitude: altitude plays an important role in the clustering of pm10 for the year 2018. The clusters are strongly influenced by the altitude of the stations. In particular, the blue cluster is characterized by stations that have a very high altitude compared to the average, on the contrary the red cluster includes stations with an altitude below the average.
- PM10: PM10 levels are very different for the blue cluster (below average) and the red cluster (above average) particularly in the first and last part of the year. The trends are different for the different clusters.
- Temp_2m: the temperature values for the different clusters are different (obviously, different latitude...) but we have the same trend.
- Wind_speed_10m_max : red and yellow clusters have the same trend and wind speed values (above average). blue and green clusters have the same trend and wind speed values (below the average).
- Precipitation: we have high precipitation for the blue cluster, particularly at the beginning and end of the year, when the pm10 level is at its minimum. Precipitation most likely reduces pm10.
- NH3 livestock: the red cluster highlights a strong presence of NH3 due to livestock, particularly in the central part of the year. The other clusters have smaller values. Probably the presence of NH3 due to livestock contributes to raising PM10.
Obviously for this cluster I have a large number of animals.
- Nh3 agr_soil: see point before. Even more evidence for the red cluster (strong presence of NH3).
- Nh3 _agr_waste_burn: see point before.
- Nox: the yellow cluster has higher nox values than all the other clusters, but it is not the most "polluted" cluster. PM10 is therefore more induced by NH3 rather than NOX
NH3 in all its origins (waste, pastures,..) negatively influences the PM10 level, more than NOX.

- Li: the number of animals is higher in the red cluster. (see comment before)
- Land: from the land use plot it can be seen that the blue cluster contains less exploited stations at ground level, therefore a positive factor for the pm10 level.
In short, nothing that wasn't known, but at least we have statistical and scientific evidence.



![](./src/figures/Curve%20PPMx/Time%20Series/plt_modemap.png)

1. Comments:
	- ecc
	- ecc


