setwd("")

files1 = file.info(list.files(pattern="*.raw"))

files1 = files1[with(files1, order(as.POSIXct(mtime))), ]

files1 <- rownames(files1)

files1 <- files1[-grep("blk", files1)]

files1 <- files1[-grep("Mix", files1)]

files1 <- files1[-grep("NCE", files1)]

files1 <- files1[-grep("dln", files1)]

files1 <- gsub(".raw","",files1)



setwd("e://Kevin/170121 iPOP plasma batch2/RPLC/pos/raw/")

files2 = file.info(list.files(pattern="*.raw"))

files2 = files2[with(files2, order(as.POSIXct(mtime))), ]

files2 <- rownames(files2)

files2 <- files2[-grep("blk", files2)]

files2 <- files2[-grep("Mix", files2)]

files2 <- files2[-grep("NCE", files2)]

files2 <- files2[-grep("dln", files2)]

files2 <- gsub(".raw","",files2)  



setwd("e://Kevin/170322 iPOP plasma batch3/RPLC/pos/raw/")

files3 = file.info(list.files(pattern="*.raw"))

files3 = files3[with(files3, order(as.POSIXct(mtime))), ]

files3 <- rownames(files3)

files3 <- files3[-grep("blk", files3)]

files3 <- files3[-grep("NCE", files3)]

files3 <- files3[-grep("dln", files3)]

files3 <- gsub(".raw","",files3)  



files <- c(files1, files2, files3)