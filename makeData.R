
# script to "create" all the datasets included in OCNet
require(OCNet)
load(file="data/OCN_4.rda")
class(OCN_4) <- "OCN"
save(OCN_4, file="data/OCN_4.rda")
tools::resaveRdaFiles("data/OCN_4.rda")


set.seed(1);
OCN_20 <- create_OCN(20,20)
save(OCN_20, file="data/OCN_20.rda")
tools::resaveRdaFiles("data/OCN_20.rda")




## To not reconstruct:
setwd("data/")
files <- system("ls *.rda", intern=T)
for (i in files) {
    load(file=i)
    vname <- stringr::str_remove_all(i, ".rda")
    tmp <- get(vname)
    class(tmp) <- "OCN"
    assign( vname, tmp)
    save(list=vname, file=i)
    tools::resaveRdaFiles(i)    
}


