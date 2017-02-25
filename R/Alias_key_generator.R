chick <- 
  read.table("/Users/Luke/Dropbox/Luke/R_projects/Ceuta_ASR_Matrix_Modeling/data/chick_survival_data.txt",
             header = TRUE, colClasses = c("factor", "character","factor",
                                           "numeric","factor","factor", "numeric"))

fledgling_adult <- 
  read.table("/Users/Luke/Dropbox/Luke/R_projects/Ceuta_ASR_Matrix_Modeling/data/fledgling_adult_survival_data.txt",
             header = TRUE, colClasses = c("factor","character","factor","factor"))

all_bird_rings <- as.factor(levels(as.factor(c(as.character(chick$ring), as.character(fledgling_adult$ring)))))
numeric_rings <- paste("bird", as.numeric(all_bird_rings), sep = "_")
ring_table <- data.frame(all_bird_rings, numeric_rings)

write.table(ring_table, file = "/Users/Luke/Dropbox/Luke/R_projects/Ceuta_ASR_Matrix_Modeling/data/Ceuta_ring_alias_key.txt")
colnames(ring_table) <- c("ring", "alias")

chick2 <- left_join(chick, ring_table, by = "ring")[c(8,2:7)]
chick3 <- left_join(chick2, brood_table, by = "brood_ID")[c(1:5,8,7)]
colnames(chick3) <- c("bird_ID", "ch", "year", "day_of_season", "sex", "brood_ID", "clutch_size")
write.table(chick3, file = "/Users/Luke/Dropbox/Luke/R_projects/Ceuta_ASR_Matrix_Modeling/data/chick_mark-recapture_data.txt", row.names = FALSE, col.names = TRUE)

fledgling_adult2 <- left_join(fledgling_adult, ring_table, by = "ring")[c(5,2:4)]
colnames(fledgling_adult2) <- c("bird_ID", "ch", "sex", "age")
write.table(fledgling_adult2, file = "/Users/Luke/Dropbox/Luke/R_projects/Ceuta_ASR_Matrix_Modeling/data/fledgling_adult_mark-recapture_data.txt", row.names = FALSE, col.names = TRUE)


length(levels(as.factor(as.character(chick2$alias))))
length(levels(as.factor(as.character(fledgling_adult2$alias))))