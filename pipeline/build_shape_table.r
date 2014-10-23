
mgw <- read.delim("MGW.csv", stringsAsFactors=FALSE, header=FALSE, sep=",")

shape.df <- data.frame(pentamer=gsub(">", "", mgw$V1[seq(from=1, to=nrow(mgw), by=2)]),
                       mgw=mgw$V3[seq(from=2, to=nrow(mgw), by=2)],
                       stringsAsFactors=FALSE)

propt <- read.delim("ProT.csv", stringsAsFactors=FALSE, header=FALSE, sep=",")

shape.df <- merge(shape.df, data.frame(pentamer=gsub(">", "", propt$V1[seq(from=1, to=nrow(propt), by=2)]),
                       prop_twist=propt$V3[seq(from=2, to=nrow(propt), by=2)],
                       stringsAsFactors=FALSE))

helt <- read.delim("HelT.csv", stringsAsFactors=FALSE, header=FALSE, sep=",")

shape.df <- merge(shape.df, data.frame(pentamer=gsub(">", "", helt$V1[seq(from=1, to=nrow(helt), by=2)]),
                            hel_twist_left=helt$V2[seq(from=2, to=nrow(helt), by=2)],
                            hel_twist_right=helt$V3[seq(from=2, to=nrow(helt), by=2)],
                            stringsAsFactors=FALSE))


roll <- read.delim("Roll.csv", stringsAsFactors=FALSE, header=FALSE, sep=",")

shape.df <- merge(shape.df, data.frame(pentamer=gsub(">", "", roll$V1[seq(from=1, to=nrow(roll), by=2)]),
                            roll_left=roll$V2[seq(from=2, to=nrow(roll), by=2)],
                            roll_right=roll$V3[seq(from=2, to=nrow(roll), by=2)],
                            stringsAsFactors=FALSE))

saveRDS(shape.df, file="shape.df.rds")
