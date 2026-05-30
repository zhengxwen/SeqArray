# SeqArray C/C++ Code Audit Results

I scanned all 17 C/C++ source files in `src/`. Here's a summary of the bugs and issues found, grouped by severity:

---

## CRITICAL (4 issues)

| # | File | Lines | Issue |
|---|------|-------|-------|
| 1 | ConvToGDS.cpp | 153–160, 196–198, 237–252 | **Mutating R immutable string objects.** Code obtains `const char*` from `CHAR(STRING_ELT(...))` then casts away const and writes via `memcpy((void*)base, ...)` and `((char*)base)[sign_pos] = ','`. This corrupts R's internal shared string pool — any other R object referencing the same `CHARSXP` will see corrupted data. |
| 2 | ConvGDS2VCF.cpp | 91–92 | ~~**Signed integer overflow in `fast_itoa()`.**~~ **FIXED.** `fast_itoa()` now checks `val == INT32_MIN` first and outputs `"NA"` directly; also rewritten with division-free magic-multiply approach. |
| 3 | ConvGDS2VCF.cpp | 81 | **Size overflow in `LineBuf_NeedSize()`.** `size_t n = p + st` can wrap to a small value, causing an undersized buffer allocation and subsequent heap overflow. |
| 4 | FileMerge.cpp | 131, 163, 231, 252 | **Integer overflow in buffer size calculations.** `nsamp * ploidy` uses `int` arithmetic with no overflow guard — if the product exceeds `INT_MAX`, the allocated vector is undersized and subsequent writes corrupt the heap. |

---

## HIGH (10 issues)

| # | File | Lines | Issue |
|---|------|-------|-------|
| 5 | ConvVCF2GDS.cpp | 989–1001, 1110–1120 | **Integer overflow in VCF integer parsing.** `val = val*10 + (ch-'0')` overflows *before* the subsequent bounds check executes — UB in signed arithmetic. |
| 6 | ConvGDS2VCF.cpp | 297, 324, 340 | **Overflow in `FORMAT_Write()`.** `(n-1)*Step` pointer arithmetic can overflow, causing out-of-bounds reads. |
| 7 | GetData.cpp | 638–680 | **Unchecked array index in `get_list()`.** `STRING_ELT(val, pt+i)` without verifying `pt+i < XLENGTH(val)`. |
| 8 | Index.cpp | 290 | **Stack buffer overflow in `AddChrom()`.** Fixed-size `string txt[65536]` — if data exceeds this, stack is corrupted. |
| 9 | Index.cpp | 95 | **Out-of-bounds access in `CIndex::GetInfo()`.** `Values[AccIndex]` after increment without bounds check. |
| 10 | Methods.cpp | 235, 361 | **Division by zero.** `FC_AF_DS_Ref` divides by `AFreq_Ploidy` (could be 0); `FC_AF_DS_Allele` does `n/m` without checking `m != 0`. |
| 11 | ReadBySample.cpp | 172–174 | **Integer overflow in buffer resize.** `CellCount *= DLen[2]` without overflow check before `resize()`. |
| 12 | ReadByVariant.cpp | 247–248 | **Unsafe 64→32 bit cast.** `IndexRaw` (C_Int64) cast to C_Int32 without range validation. |
| 13 | FileMerge.cpp | 278–316, 343–408 | **Memory leak on exception.** `new CApply_Variant_*` objects in `NodeList`/`Files` are never deleted if an exception is thrown. |
| 14 | Methods.cpp | 69 | **Null pointer dereference.** `INTEGER(GET_DIM(Geno))` without null check — if `Geno` has no dim attribute, crash. |

---

## MEDIUM (8 issues)

| # | File | Lines | Issue |
|---|------|-------|-------|
| 15 | ConvToGDS.cpp | 250–262 | **Unprotected SEXP.** Old `dosage` SEXP could be GC'd while `src` pointer still references its memory. |
| 16 | GetData.cpp | 301–302 | **Truncation with long chromosome names.** `snprintf(buf, 1024, "%s:%d", Chrom[i].c_str(), ...)` silently truncates. |
| 17 | Index.cpp | 246 | **Integer overflow in sum.** `ntot += len[i]` without overflow check. |
| 18 | Index.cpp | 549 | **Unchecked multiplication.** `(pE - pSt) * numPloidy` can overflow. |
| 19 | vectorization.cpp | 1650+ | **SIMD overread.** SSE/AVX loads process 16/32 bytes without verifying remaining buffer length matches vector width. |
| 20 | samtools_ext.c | 57–60 | **Integer truncation.** `(unsigned int)(size*nitems)` truncates on 64-bit systems. |
| 21 | ConvVCF2GDS.cpp | 977, 1014 | **Shared buffer modification.** Setting `*end = 0` for `strtod()` relies on undocumented extra buffer space. |
| 22 | Methods.cpp | 125 | **Overflow before comparison.** `m * num_samp` overflows in `int` before being compared to `n`. |

---

## Key Recommendations

1. **Issue #1 is the most impactful real-world bug** — mutating R's CHARSXP objects causes silent data corruption. Fix by `PROTECT`ing a new `mkChar()` string instead of modifying in place.
2. **Issue #2** — handle INT32_MIN explicitly: `if (val == INT32_MIN) { /* write "-2147483648" directly */ }`.
3. **Issue #3** — add overflow check: `if (st > SIZE_MAX - p) Rf_error("buffer overflow")`.
4. **Issues #4, #11, #12** — use `size_t` or add overflow guards before arithmetic used for allocation sizes.
5. **Issue #5** — perform bounds check *before* the multiplication: `if (val > (INT_MAX - (ch-'0')) / 10)`.
6. **Issue #13** — use RAII (e.g., `std::unique_ptr`) or a cleanup scope guard for dynamic allocations.
