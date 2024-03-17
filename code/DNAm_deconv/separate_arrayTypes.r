library(minfi)

args<-commandArgs(trailingOnly=TRUE)

cohort = args[1]
savedir <- args[2]

ls.files <- list.files(cohort)


i <- 1
files <- c()
arrayTypes <- c()

for (file in ls.files){
    if (grepl("idat", file, fixed=TRUE)){
        print(i)
        filename <- paste(cohort, file, sep="")
        data <- read.metharray(filename)
        files[i] <- filename
        arrayTypes[i] <- annotation(data)[["array"]]
        i <- i+1
        }
    }
    
files <- as.data.frame(files)
arrayTypes <- as.data.frame(arrayTypes)

df <- cbind(files, arrayTypes)

unique(df$arrayTypes)
write.csv(df, paste0(savedir, "/arrayTypes.csv"))
