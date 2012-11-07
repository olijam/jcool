THREADS <- 5
library(doMC); registerDoMC(THREADS)
library(Ripf)

vars <- c("keywords","")
day <- read.csv("../data/pday4655.csv.gz",as.is=TRUE)

keywords <- day$keywords
keep <- names(sort(table(keywords),dec=TRUE)[1:2000])
keywords[!(keywords %in% keep)] <- "keywords"

pres <- foreach(i=1:200) %dopar% {
  print(sprintf("Thread %d:\n",i))
  r <- sample(1:length(keywords),200000)
  f.keywords <- factor(keywords[r])
  m.keywords <- list()
  for(j in 1:nlevels(f.keywords)) {
    m.keywords[[levels(f.keywords)[j]]] <- as.integer(f.keywords) == j
  }
  m.keywords <- as.data.frame(m.keywords)

  names(minfo(day$clicks[r],m.keywords,iterate=TRUE,progress=TRUE,count=30))
}
