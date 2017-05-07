if (Sys.info()[['sysname']] != "Windows")
{
	# according to the limit of 32-bit Windows
	BiocGenerics:::testPackage("SeqArray")
}
