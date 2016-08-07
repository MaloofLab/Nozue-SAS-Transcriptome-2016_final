library(reshape2)
library(lmerTest)
library(plyr)

categories <- read.csv("~/Downloads/test.wozero.melt.hormone.signifiant.s.csv",row.names=1,as.is=TRUE)

summary(categories)
head(categories)

load("~/Downloads/summary.vst.response.kazu.Rdata")

head(summary.vst.response.kazu[1:10])

summary.vst.response.kazu$AGI <- row.names(summary.vst.response.kazu)

expr.melt <- melt(summary.vst.response.kazu,id.vars= "AGI",variable.name = "sample", value.name = "log2.dif")

head(expr.melt)

expr.melt$gt <- as.factor(regmatches(expr.melt$sample,regexpr("^[0-9]{1,2}",expr.melt$sample)))

expr.melt$time <- as.factor(regmatches(expr.melt$sample,regexpr("[0-9]{1,2}hr",expr.melt$sample)))

expr.melt$rep <- as.factor(regmatches(expr.melt$sample,regexpr("[A-Z]$",expr.melt$sample)))

head(expr.melt)
tail(expr.melt)

summary(expr.melt)

summary(categories)
table(categories$my.category)

expr.melt.merge <- merge(expr.melt,categories,by="AGI",all=FALSE)

expr.melt.merge$gt <- relevel(expr.melt.merge$gt,ref="10")

table(expr.melt.merge$gt,expr.melt.merge$time)

mutants <- levels(expr.melt.merge$gt)[-1] #1 is Col (Gt 10)

category.names<- unique(expr.melt.merge$my.category)

results <- lapply(mutants,function(gt) { 
  sapply(category.names,function(cat.name) {
    tmp <- expr.melt.merge[(expr.melt.merge$gt=="10" | expr.melt.merge$gt==gt) & expr.melt.merge$my.category==cat.name,]
    lmer1 <- lmer(log2.dif ~ gt*time + (1|AGI), data=tmp)
    tmp.results <- anova(lmer1)[, "Pr(>F)"]
    names(tmp.results) <- row.names(anova(lmer1))
    tmp.results[-2] # we do not care about time
    })
} )

names(results) <- mutants

results

results.df <- ldply(results,.id="gt")

results.df <- cbind(test=rep(c("gt","gt:time")), results.df)

results.df

results.df.adj <- matrix(
  p.adjust(unlist(results.df[,-1:-2]),"holm"),
  nrow=nrow(results.df))

results.df.adj <- cbind(results.df[,1:2],results.df.adj)

colnames(results.df.adj) <- colnames(results.df)

results.df.adj


