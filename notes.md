# WORKFLOW SUGGESTIONS
Test with only ppmsuit simple models (spatial ppm with functional data also)
and variations there to include covariates.
Ie start with the "incomplete" model to see how they behave even if not being perfect, then we will experience the more complex ones.

The request however is only to understand the papers and to know how to implement them; it is an application project in fact, there is nothing new or surprising to propose.
That is, modifications to basic algorithms are appreciated, but they are not "necessary", and must not be exceptional, just small improvements reasoned.
There is therefore to understand well the tools that others have created (packages pmmsuite, drpm, etc.).
So the question is much easier (for the basic level): understand the papers and use the tools to implement them. Others have templates to write by themselves, but because theirs are easier, here the papers deal with difficult things
knowing how to understand the effects of parameter changes

# CODE SUGGESTIONS
When we have more clusters, use the adjusted random index to measure how similar the clusters are. That is, it is not enough just plotting with colors, it is necessary to be quantitative.

From the packages we get MCMC with cluster labels, that is, each line has its units with their assigned labels. From the matrix we then have to obtain a single cluster for all the units: for that it uses the binder loss function, with the packet salso, that it chooses the best clustering.
For sauce serves R>4.2