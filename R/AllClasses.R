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
            return("object should inherited from 'gds.class'.")

        # validObject() may validate this superclass against a hollow S4
        # coercion shell (e.g. as(, "SeqVarGDSClass")) that carries
        # no live gds connection; only a live handle is an S3 list. Skip the
        # shell so it does not error on ls.gdsn()/`$root`.
        if (typeof(object) != "list")
            return(TRUE)

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
