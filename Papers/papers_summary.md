### paper BayesianPoker
They say they didnt actually did a good work, so there is room for us.    
Their model just worked on two players, the bot and a second player, losed basically all matches against humans, while ok being a bit better than basic models (like a random player).

So this idea may be ok for us but we need to think really well about designing the model. Like setting priors, something like   

- for each player, to learn their strategies (aggressive/conservative) during the match
- for our style of play in the match
- for each player about their winning expecation in each round
- for our winning expectation in each round
- etc, so may be difficult to model