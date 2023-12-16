library(readVCF)

filepath <- system.file("extdata", "LCT.vcf", package = "readVCF")
A <- readVCFgenotypes(filepath)
stopifnot(sum(A) == 107223)

filepathgz <- system.file("extdata", "LCT.vcf.gz", package = "readVCF")
B <- readVCFgenotypes2(filepathgz, character(0))
stopifnot(all(A == B))

C <- readVCFgenotypes2(filepathgz, c("2:136402646-136402781", "2:136408821-136408831"))
stopifnot(all( C == A[,match(colnames(C), colnames(A))]))

filepathi <- system.file("extdata", "index.vcf.gz", package = "readVCF")
A <- readVCFgenotypes2(filepathi, c("1:9999943-9999950", "10:3000100-3000300"))


