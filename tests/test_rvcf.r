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

# juste pour vérifier que pas de leaks
regs <- getRegions(object)
chroms <- getChroms(object)

#object <- openVCF("~/Stage/big.vcf.gz")

# ------------- reading and writing dosage file --------------

#to change with your local file that has dosages

#to note : writeDosage will give back samples
#writeDosage("inst/extdata/dosages.vcf.gz", "truc", character(0))
#dos_object <- openVCF("inst/extdata/dosages.vcf.gz")
#for potentially reading the file in R 
#written <- readBin("file_to_read","double",1000000000,size=4) #cos float 4 bytes
#VCFnext(dos_object)
#val <- getLine(dos_object, "DS")$values
