filepath <- system.file("extdata", "LCT.vcf", package = "readVCF")
A <- readVCFgenotypes(filepath)

filepathgz <- system.file("extdata", "LCT.vcf.gz", package = "readVCF")
B <- readVCFgenotypes2(filepathgz, character(0))
C <- readVCFgenotypes2(filepathgz, c("2:136402646-136402781", "2:136408821-136408831"))

