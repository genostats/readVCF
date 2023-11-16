filepath <- system.file("extdata", "LCT.vcf.gz", package = "readVCF")
test_htsVCF(filepath, c("2:136402646-136402781", "2:136408821-136408831"))
test_pointer(filepath, c("2:136402646-136402781", "2:136408821-136408831"))


filepath2 <- system.file("extdata", "index.vcf.gz", package = "readVCF")
test_htsVCF("inst/extdata/index.vcf.gz", "10:-")

