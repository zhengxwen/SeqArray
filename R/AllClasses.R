# register old-style (S3) classes and inheritance
setOldClass("gds.class")
setOldClass("gdsn.class")


# create class definitions
setClass("SeqVarGDSClass", contains="gds.class")
# setClass("SeqVarNodeNewInfoClass", slots = c(
#         gdsn = "gdsn.class", number = "integer", type = "character"
#     ))


# test the validity of objects
setValidity("SeqVarGDSClass",
    function(object)
    {
        if (!inherits(object, "gds.class"))
            return("object should inherited from `gds.class'.")

        n <- index.gdsn(object, "description", silent=TRUE)
        if (is.null(n))
            return("Description variable must exist!")

        var.names <- ls.gdsn(object)
        if (!all(c("sample.id", "variant.id", "position",
            "chromosome", "allele", "genotype") %in% var.names))
        {
            return("sample.id, variant.id, position, chromosome, allele, and genotype are required variables.")
        }
        TRUE
    }
)
