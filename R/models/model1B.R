#!/usr/local/bin/R -q --no-save
args <- commandArgs(trailingOnly==TRUE)

print(args)
stop()
THREADS <- 10
ITERATIONS <- 400
#library(bigmemory)
#options(bigmemory.typecast.warning=FALSE)
#library(biganalytics)
#library(biglm)
library(doMC); registerDoMC(THREADS)
library(Rminfo)

var_cols <- c(
  "hour_id",
#  "date",
#  "wday",
#  "tday",
  "campaign_name",
  "daily_budget",
  "company_name",
  "advertiser_name",
  "cost_type",
  "adx_category",
  "adx_product",
  "campaign_type",
  "ad_size",
  "ad_position",
  "ad_service_domain",
  "ad_file_type",
  "ad_content_type1",
  "ad_content_type2",
  "recency",
  "domain",
  "metro_name",
  "region_name",
  "country_name",
  "dma_code",
  "keywords",
  "context_l1",
  "context_l2",
  "exchange")

recs <- read.csv("~/data/day4689.csv.gz",as.is=TRUE)
gc()

pres <- foreach(i=1:ITERATIONS) %dopar% {
  r <- sample(1:nrow(recs), nrow(recs) * 0.20)
  t <- minfo(recs$clicks[r],recs[r,var_cols],count=15,progress=TRUE)
  gc()
  return(t)
}

t <- data.frame()

for(i in 1:length(pres)) t <- rbind(t,cbind(pos=rownames(pres[[i]]),pres[[i]]))

t2 <- aggregate(cbind(1,pos,minfo,levels) ~ factor + level, data=t, sum)
t2$pos <-  t2$pos/t2$"1"
t2$minfo <- t2$minfo/t2$"1"
t2$levels <- t2$levels/t2$"1"

t2$impressions <- 1
t2$clicks <- 1
t2$ctr <- 1

for(i in 1:nrow(t2)) {
  r <- which(recs[,t2$factor[i]] == t2$level[i])
  t2$impressions[i] <-  sum(recs$impressions[r])
  t2$clicks[i]      <-  sum(recs$clicks[r])
  t2$ctr[i]   <- t2$clicks[i] / t2$impressions[i]
}
t2 <- rbind(t2,
  cbind(
    factor="ALL",
    level="ALL",
    "1"=ITERATIONS,
    pos=ITERATIONS,
    minfo=1,
    levels=1,
    impressions=sum(recs$impressions),
    clicks=sum(recs$clicks),
    ctr=sum(recs$clicks)/sum(recs$impressions)
  ))

  write.csv(t2,file="day4689_report.csv")
