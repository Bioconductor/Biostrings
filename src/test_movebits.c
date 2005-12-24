#include <stdio.h>
#include <limits.h>

#define CHAR_SIZE               (sizeof(char))
#define LONG_SIZE               (sizeof(long))
#define BITS_PER_CHAR           CHAR_BIT
#define BITS_PER_SIZEOF_UNIT    (BITS_PER_CHAR / (int) CHAR_SIZE)
#define BITS_PER_LONG           ((int) LONG_SIZE * BITS_PER_SIZEOF_UNIT)

typedef unsigned long ShiftOrWord_t; /* hopefully this will be 32-bit
                                      * on 32-bit machines and 64-bit
                                      * on 64-bit machines */

static void printULBits(unsigned long bits)
{
    unsigned long current_bit = 1UL << (BITS_PER_LONG-1);
    int i;

    for (i = 0; i < BITS_PER_LONG; i++) {
        printf("%d", (bits & current_bit) != 0UL);
        if ((i % 8) == 7) {
            printf(" ");
        }
        current_bit >>= 1;
    }
    printf("-> %uL\n", bits);
    return;
}

/*
 * M_k_len MUST be >= 1
 */
static void debug_ShowBits(int M_k_len, ShiftOrWord_t *M_k)
{
    static int i;

    for (i = 0; i < M_k_len; i++) {
        printf("M_k[%d] = ", i);
        printULBits(M_k[i]);
    }
    return;
}

/*
 * M_k_len MUST be >= 1
 */
static void debug_InitBits(int M_k_len, ShiftOrWord_t *M_k)
{
    static int i;

    M_k[0] = ~0UL;
    for (i = 1; i < M_k_len; i++) {
        M_k[i] = M_k[i-1] << 1;
    }
    return;
}

/*
 * M_k_len MUST be >= 1
 */
static void debug_MoveBits(int M_k_len, ShiftOrWord_t *M_k, ShiftOrWord_t U_xptr_k)
{
    static ShiftOrWord_t tmpA, tmpB;
    static int i;

    tmpA = M_k[0];
    M_k[0] = (tmpA << 1) | U_xptr_k;
    for (i = 1; i < M_k_len; i++) {
        tmpB = tmpA;
        tmpA = M_k[i];
        M_k[i] = ((tmpA << 1) | U_xptr_k) &
                 tmpB & (tmpB << 1) & (M_k[i-1] << 1);
    }
    return;
}

int main()
{
    ShiftOrWord_t M_k[4], mover;

    printULBits(5UL);
    printULBits(~0UL);
    printULBits(1UL << (BITS_PER_LONG-1));
    debug_InitBits(4, M_k);
    debug_ShowBits(4, M_k);

    debug_InitBits(4, M_k);
    mover = 0UL;
    printf("With mover = ");
    printULBits(mover);
    debug_MoveBits(4, M_k, mover);
    debug_ShowBits(4, M_k);

    debug_InitBits(4, M_k);
    mover = 1UL;
    printf("With mover = ");
    printULBits(mover);
    debug_MoveBits(4, M_k, mover);
    debug_ShowBits(4, M_k);

    debug_InitBits(4, M_k);
    mover = 2UL;
    printf("With mover = ");
    printULBits(mover);
    debug_MoveBits(4, M_k, mover);
    debug_ShowBits(4, M_k);

    debug_InitBits(4, M_k);
    mover = 3UL;
    printf("With mover = ");
    printULBits(mover);
    debug_MoveBits(4, M_k, mover);
    debug_ShowBits(4, M_k);

    debug_InitBits(4, M_k);
    mover = 4UL;
    printf("With mover = ");
    printULBits(mover);
    debug_MoveBits(4, M_k, mover);
    debug_ShowBits(4, M_k);

    debug_InitBits(4, M_k);
    mover = 5UL;
    printf("With mover = ");
    printULBits(mover);
    debug_MoveBits(4, M_k, mover);
    debug_ShowBits(4, M_k);

    debug_InitBits(4, M_k);
    mover = 2UL + 8UL + 32UL + 128UL + 512UL + 2048UL + 8192UL + 32768UL;
    printf("With mover = ");
    printULBits(mover);
    debug_MoveBits(4, M_k, mover);
    debug_ShowBits(4, M_k);

    return 0;
}
