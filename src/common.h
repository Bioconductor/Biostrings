#ifndef _R_COMMON_H_
#define _R_COMMON_H_

#define TEST_CALL_DEF

#ifdef TEST_CALL_DEF
#if defined __GNUC__ && __GNUC__ >= 2
/*
 * These macros may work with other compilers. However, they are
 * just for compile time error checking for wrong .Call and
 * .External function registration and so, need not be used except
 * with gcc (which is used by the developers).
 */

#define PREPROC_NULL_ARGS(fname, n) (fname)(PREPROC_SEQ_ARGS(n, R_NilValue))

#define PREPROC_SEQ_ARGS(n, v) PREPROC_SEQ_ ## n (v)

#define PREPROC_SEQ_0(v) /* nothing */
#define PREPROC_SEQ_1(v) (v)
#define PREPROC_SEQ_2(v) (v), (v)
#define PREPROC_SEQ_3(v) (v), (v), (v)
#define PREPROC_SEQ_4(v) PREPROC_SEQ_2(v), PREPROC_SEQ_2(v)
#define PREPROC_SEQ_5(v) PREPROC_SEQ_3(v), PREPROC_SEQ_2(v)
#define PREPROC_SEQ_6(v) PREPROC_SEQ_3(v), PREPROC_SEQ_3(v)
#define PREPROC_SEQ_7(v) PREPROC_SEQ_4(v), PREPROC_SEQ_3(v)
#define PREPROC_SEQ_8(v) PREPROC_SEQ_4(v), PREPROC_SEQ_4(v)
#define PREPROC_SEQ_9(v) PREPROC_SEQ_5(v), PREPROC_SEQ_4(v)
#define PREPROC_SEQ_10(v) PREPROC_SEQ_5(v), PREPROC_SEQ_5(v)
#define PREPROC_SEQ_11(v) PREPROC_SEQ_6(v), PREPROC_SEQ_5(v)
#define PREPROC_SEQ_12(v) PREPROC_SEQ_6(v), PREPROC_SEQ_6(v)
#define PREPROC_SEQ_13(v) PREPROC_SEQ_7(v), PREPROC_SEQ_6(v)
#define PREPROC_SEQ_14(v) PREPROC_SEQ_7(v), PREPROC_SEQ_7(v)
#define PREPROC_SEQ_15(v) PREPROC_SEQ_8(v), PREPROC_SEQ_7(v)
#define PREPROC_SEQ_16(v) PREPROC_SEQ_8(v), PREPROC_SEQ_8(v)
#define PREPROC_SEQ_17(v) PREPROC_SEQ_9(v), PREPROC_SEQ_8(v)
#define PREPROC_SEQ_18(v) PREPROC_SEQ_9(v), PREPROC_SEQ_9(v)
#define PREPROC_SEQ_19(v) PREPROC_SEQ_10(v), PREPROC_SEQ_9(v)
#define PREPROC_SEQ_20(v) PREPROC_SEQ_10(v), PREPROC_SEQ_10(v)
#define PREPROC_SEQ_21(v) PREPROC_SEQ_11(v), PREPROC_SEQ_10(v)
#define PREPROC_SEQ_22(v) PREPROC_SEQ_11(v), PREPROC_SEQ_11(v)
#define PREPROC_SEQ_23(v) PREPROC_SEQ_12(v), PREPROC_SEQ_11(v)
#define PREPROC_SEQ_24(v) PREPROC_SEQ_12(v), PREPROC_SEQ_12(v)
#define PREPROC_SEQ_25(v) PREPROC_SEQ_13(v), PREPROC_SEQ_12(v)
#define PREPROC_SEQ_26(v) PREPROC_SEQ_13(v), PREPROC_SEQ_13(v)
#define PREPROC_SEQ_27(v) PREPROC_SEQ_14(v), PREPROC_SEQ_13(v)
#define PREPROC_SEQ_28(v) PREPROC_SEQ_14(v), PREPROC_SEQ_14(v)
#define PREPROC_SEQ_29(v) PREPROC_SEQ_15(v), PREPROC_SEQ_14(v)
#define PREPROC_SEQ_30(v) PREPROC_SEQ_15(v), PREPROC_SEQ_15(v)
#define PREPROC_SEQ_31(v) PREPROC_SEQ_16(v), PREPROC_SEQ_15(v)
#define PREPROC_SEQ_32(v) PREPROC_SEQ_16(v), PREPROC_SEQ_16(v)
#define PREPROC_SEQ_33(v) PREPROC_SEQ_17(v), PREPROC_SEQ_16(v)
#define PREPROC_SEQ_34(v) PREPROC_SEQ_17(v), PREPROC_SEQ_17(v)
#define PREPROC_SEQ_35(v) PREPROC_SEQ_18(v), PREPROC_SEQ_17(v)
#define PREPROC_SEQ_36(v) PREPROC_SEQ_18(v), PREPROC_SEQ_18(v)
#define PREPROC_SEQ_37(v) PREPROC_SEQ_19(v), PREPROC_SEQ_18(v)
#define PREPROC_SEQ_38(v) PREPROC_SEQ_19(v), PREPROC_SEQ_19(v)
#define PREPROC_SEQ_39(v) PREPROC_SEQ_20(v), PREPROC_SEQ_19(v)
#define PREPROC_SEQ_40(v) PREPROC_SEQ_20(v), PREPROC_SEQ_20(v)
#define PREPROC_SEQ_41(v) PREPROC_SEQ_21(v), PREPROC_SEQ_20(v)
#define PREPROC_SEQ_42(v) PREPROC_SEQ_21(v), PREPROC_SEQ_21(v)
#define PREPROC_SEQ_43(v) PREPROC_SEQ_22(v), PREPROC_SEQ_21(v)
#define PREPROC_SEQ_44(v) PREPROC_SEQ_22(v), PREPROC_SEQ_22(v)
#define PREPROC_SEQ_45(v) PREPROC_SEQ_23(v), PREPROC_SEQ_22(v)
#define PREPROC_SEQ_46(v) PREPROC_SEQ_23(v), PREPROC_SEQ_23(v)
#define PREPROC_SEQ_47(v) PREPROC_SEQ_24(v), PREPROC_SEQ_23(v)
#define PREPROC_SEQ_48(v) PREPROC_SEQ_24(v), PREPROC_SEQ_24(v)
#define PREPROC_SEQ_49(v) PREPROC_SEQ_25(v), PREPROC_SEQ_24(v)
#define PREPROC_SEQ_50(v) PREPROC_SEQ_25(v), PREPROC_SEQ_25(v)
#define PREPROC_SEQ_51(v) PREPROC_SEQ_26(v), PREPROC_SEQ_25(v)
#define PREPROC_SEQ_52(v) PREPROC_SEQ_26(v), PREPROC_SEQ_26(v)
#define PREPROC_SEQ_53(v) PREPROC_SEQ_27(v), PREPROC_SEQ_26(v)
#define PREPROC_SEQ_54(v) PREPROC_SEQ_27(v), PREPROC_SEQ_27(v)
#define PREPROC_SEQ_55(v) PREPROC_SEQ_28(v), PREPROC_SEQ_27(v)
#define PREPROC_SEQ_56(v) PREPROC_SEQ_28(v), PREPROC_SEQ_28(v)
#define PREPROC_SEQ_57(v) PREPROC_SEQ_29(v), PREPROC_SEQ_28(v)
#define PREPROC_SEQ_58(v) PREPROC_SEQ_29(v), PREPROC_SEQ_29(v)
#define PREPROC_SEQ_59(v) PREPROC_SEQ_30(v), PREPROC_SEQ_29(v)
#define PREPROC_SEQ_60(v) PREPROC_SEQ_30(v), PREPROC_SEQ_30(v)
#define PREPROC_SEQ_61(v) PREPROC_SEQ_31(v), PREPROC_SEQ_30(v)
#define PREPROC_SEQ_62(v) PREPROC_SEQ_31(v), PREPROC_SEQ_31(v)
#define PREPROC_SEQ_63(v) PREPROC_SEQ_32(v), PREPROC_SEQ_31(v)
#define PREPROC_SEQ_64(v) PREPROC_SEQ_32(v), PREPROC_SEQ_32(v)
#define PREPROC_SEQ_65(v) PREPROC_SEQ_33(v), PREPROC_SEQ_32(v)


/*
 * Create the .Call registration definition and do a compile time
 * check to make sure fname is a function taking nargs SEXP's and
 * returning a SEXP as value. If there is a type mismatch, gcc would
 * produce a warning. If the number of arguments do not match, gcc
 * will produce an error.
 * 
 * This only works if nargs is an integer >= 0 and <= 65. Currently
 * these are the only valid nargs values for R .Call functions.
 *
 * This depends on sizeof() returning a compile time constant without
 * evaluating its argument.
 * 
 */
#define CALL_DEF(fname, nargs) \
    { #fname, (DL_FUNC)&(fname), nargs+sizeof(TYPEOF(PREPROC_NULL_ARGS(fname, nargs)))-sizeof(TYPEOF(R_NilValue))}

#define EXTERNAL_DEF(fname) CALL_DEF(fname, 1)

#endif
#endif

#ifndef CALL_DEF
#define CALL_DEF(fname, nargs) { #fname, (DL_FUNC)&fname, nargs}
#endif

#endif
