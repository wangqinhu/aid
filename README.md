Installation
============

```R
library(devtools)
install_github(repo = "wangqinhu/aid")
```

Syntax
======

```R
library(aid)
# grade data
demo1 <- system.file("extdata", "demo1.tsv", package="aid")
aid(demo1, type = "grade", alternative="less", paired=FALSE)
# lesion data
demo2 <- system.file("extdata", "demo2.tsv", package="aid")
aid(demo2, type = "lesion", alternative="two.sided", paired=FALSE)
# biomass data
demo3 <- system.file("extdata", "demo3.tsv", package="aid")
aid(demo3, type = "biomass", alternative="greater", paired=FALSE)
```
