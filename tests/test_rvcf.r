library(readVCF)

filepath <- system.file("extdata", "LCT.vcf", package = "readVCF")
A <- readVCF:::readVCFgenotypes(filepath)
stopifnot(sum(A) == 107223)

filepathgz <- system.file("extdata", "LCT.vcf.gz", package = "readVCF")
B <- readVCF:::readVCFgenotypes2(filepathgz, character(0))
stopifnot(all(A == B))

C <- readVCF:::readVCFgenotypes2(filepathgz, c("2:136402646-136402781", "2:136408821-136408831"))
stopifnot(all( C == A[,match(colnames(C), colnames(A))]))


# -------------- opening without regions ----------------
object <- openVCF(filepathgz)

s <- getSamples(object)
stopifnot(all(s == rownames(A)))

# lire une ligne
stopifnot(VCFnext(object))
a <- getLine(object, "GT")
stopifnot(all(a$genotypes == A[,1]))

# lire la suivante
stopifnot(VCFnext(object))
b <- getLine(object, "GT")
stopifnot(all(b$genotypes == A[,2]))

# ------------- opening with regions --------------

object <- openVCF(filepathgz, regions = c("2:136402646-136402781", "2:136408821-136408831"))

s <- getSamples(object)
stopifnot(all(s == rownames(C)))

# lire une ligne
stopifnot(VCFnext(object))
a <- getLine(object, "GT")
stopifnot(all(a$genotypes == C[,1]))

# lire la suivante
stopifnot(VCFnext(object))
b <- getLine(object, "GT")
stopifnot(all(b$genotypes == C[,2]))
