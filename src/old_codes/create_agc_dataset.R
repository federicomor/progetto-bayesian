
load("../data/data_agc_lomb.Rdata")


# Dividi il dataset in due parti
parte1 <- AGC_Dataset[1:961915, ]  
parte2 <- AGC_Dataset[961916:nrow(AGC_Dataset), ]

# Salva le due parti in file RData separati
save(parte1, file = "../data/data_agc_lomb_part1.RData")
save(parte2, file = "../data/data_agc_lomb_part2.RData")
