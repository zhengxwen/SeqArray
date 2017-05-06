if (Sys.info()[['sysname']] != "Windows")
{
	# memory issue on 32-bit Windows
	BiocGenerics:::testPackage("SeqArray")
}
