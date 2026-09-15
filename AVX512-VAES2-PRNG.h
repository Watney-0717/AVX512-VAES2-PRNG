#ifndef ULTRAFAST_RNG32X_AVX512_8WAY_CHEAPSPLIT_HMIX1_FINALDIAG_VAES2_SANDWICH_OPT_LUTFIXED_FINALCANDIDATE_6_2_H
#define ULTRAFAST_RNG32X_AVX512_8WAY_CHEAPSPLIT_HMIX1_FINALDIAG_VAES2_SANDWICH_OPT_LUTFIXED_FINALCANDIDATE_6_2_H

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <immintrin.h>

#ifndef WEYL_INC_32
#define WEYL_INC_32 0x9E3779B9u
#endif

typedef struct {
    __m512i gA, gB;
    __m512i vaes_rk;
} UFR32_8Way_CheapSplit_HMix1_VAES1;

static inline __m512i ufr32_rotl32(__m512i x, int r) {
    return _mm512_or_si512(_mm512_slli_epi32(x, r),
                           _mm512_srli_epi32(x, 32-r));
}

static inline __m512i ufr32_gfni(__m512i x) {
    const __m512i gm = _mm512_set1_epi64(
        (long long)UINT64_C(0x1F3E7CFCF9F3E7CE));
    return _mm512_gf2p8affine_epi64_epi8(x, gm, 0);
}

/* Lane 2: VPSHUFB + VPTERNLOGD */
static inline __m512i ufr32_pshufb_ternlog(__m512i x) {
    const __m512i mask = _mm512_set1_epi8(0x0F);
    static const uint8_t lut_lo_bytes[64] __attribute__((aligned(64))) = {
        0x00,0x07,0x0D,0x0B,0x06,0x0C,0x05,0x09,
        0x0E,0x03,0x0F,0x01,0x08,0x04,0x02,0x0A,
        0x00,0x07,0x0D,0x0B,0x06,0x0C,0x05,0x09,
        0x0E,0x03,0x0F,0x01,0x08,0x04,0x02,0x0A,
        0x00,0x07,0x0D,0x0B,0x06,0x0C,0x05,0x09,
        0x0E,0x03,0x0F,0x01,0x08,0x04,0x02,0x0A,
        0x00,0x07,0x0D,0x0B,0x06,0x0C,0x05,0x09,
        0x0E,0x03,0x0F,0x01,0x08,0x04,0x02,0x0A};
    static const uint8_t lut_hi_bytes[64] __attribute__((aligned(64))) = {
        0x00,0x0E,0x0B,0x06,0x09,0x01,0x0D,0x0C,
        0x03,0x0F,0x05,0x08,0x07,0x0A,0x04,0x02,
        0x00,0x0E,0x0B,0x06,0x09,0x01,0x0D,0x0C,
        0x03,0x0F,0x05,0x08,0x07,0x0A,0x04,0x02,
        0x00,0x0E,0x0B,0x06,0x09,0x01,0x0D,0x0C,
        0x03,0x0F,0x05,0x08,0x07,0x0A,0x04,0x02,
        0x00,0x0E,0x0B,0x06,0x09,0x01,0x0D,0x0C,
        0x03,0x0F,0x05,0x08,0x07,0x0A,0x04,0x02};
    const __m512i lut_lo = _mm512_load_si512((const void *)lut_lo_bytes);
    const __m512i lut_hi = _mm512_load_si512((const void *)lut_hi_bytes);

    __m512i lo = _mm512_and_si512(x, mask);
    __m512i hi = _mm512_and_si512(_mm512_srli_epi16(x, 4), mask);
    __m512i y0 = _mm512_shuffle_epi8(lut_lo, lo);
    __m512i y1 = _mm512_slli_epi16(_mm512_shuffle_epi8(lut_hi, hi), 4);

    return _mm512_ternarylogic_epi32(y0, y1, x, 0xFE);
}

/* ChaCha4: seed-side only. */
static inline void ufr32_chacha_qr(
    __m512i *a, __m512i *b, __m512i *c, __m512i *d)
{
    *a = _mm512_add_epi32(*a, *b);
    *d = _mm512_xor_si512(*d, *a);
    *d = ufr32_rotl32(*d, 16);
    *c = _mm512_add_epi32(*c, *d);
    *b = _mm512_xor_si512(*b, *c);
    *b = ufr32_rotl32(*b, 12);
    *a = _mm512_add_epi32(*a, *b);
    *d = _mm512_xor_si512(*d, *a);
    *d = ufr32_rotl32(*d, 8);
    *c = _mm512_add_epi32(*c, *d);
    *b = _mm512_xor_si512(*b, *c);
    *b = ufr32_rotl32(*b, 7);
}

static inline void ufr32_chacha4(__m512i x[16]) {
    for (int round = 0; round < 2; ++round) {
        ufr32_chacha_qr(&x[0], &x[4], &x[8],  &x[12]);
        ufr32_chacha_qr(&x[1], &x[5], &x[9],  &x[13]);
        ufr32_chacha_qr(&x[2], &x[6], &x[10], &x[14]);
        ufr32_chacha_qr(&x[3], &x[7], &x[11], &x[15]);
        ufr32_chacha_qr(&x[0], &x[5], &x[10], &x[15]);
        ufr32_chacha_qr(&x[1], &x[6], &x[11], &x[12]);
        ufr32_chacha_qr(&x[2], &x[7], &x[8],  &x[13]);
        ufr32_chacha_qr(&x[3], &x[4], &x[9],  &x[14]);
    }
}

/* Seed-side: GFNI + VPSHUFB/TERNLOG + ChaCha4. */
static inline void ufr32_seed_mix3(
    __m512i seed_a, __m512i seed_b,
    UFR32_8Way_CheapSplit_HMix1_VAES1 *r)
{
    __m512i a = ufr32_gfni(seed_a);
    __m512i b = ufr32_pshufb_ternlog(seed_b);
    __m512i c[16];

    c[0]=seed_a; c[1]=ufr32_rotl32(seed_a,7);
    c[2]=ufr32_rotl32(seed_a,13); c[3]=ufr32_rotl32(seed_a,19);
    c[4]=seed_b; c[5]=ufr32_rotl32(seed_b,5);
    c[6]=ufr32_rotl32(seed_b,11); c[7]=ufr32_rotl32(seed_b,17);
    c[8]=_mm512_xor_si512(seed_a,seed_b);
    c[9]=ufr32_rotl32(c[8],3); c[10]=ufr32_rotl32(c[8],9);
    c[11]=ufr32_rotl32(c[8],15);
    c[12]=_mm512_set1_epi32(0x61707865);
    c[13]=_mm512_set1_epi32(0x3320646E);
    c[14]=_mm512_set1_epi32(0x79622D32);
    c[15]=_mm512_set1_epi32(0x6B206574);

    ufr32_chacha4(c);

    __m512i cc0=_mm512_xor_si512(c[0],c[8]);
    __m512i cc1=_mm512_xor_si512(c[3],c[11]);
    __m512i ab=_mm512_add_epi32(a,b);
    __m512i bc=_mm512_add_epi32(b,cc0);

    r->gA=_mm512_ternarylogic_epi32(a,ab,cc1,0x96);
    r->gB=_mm512_ternarylogic_epi32(b,bc,cc0,0xE8);
    r->vaes_rk=_mm512_ternarylogic_epi32(r->gA,r->gB,cc1,0x96);

    volatile __m512i *vp=(volatile __m512i *)c;
    for (unsigned i=0;i<16;++i) vp[i]=_mm512_setzero_si512();
}

static inline void ufr32_split4(__m512i g, __m512i v[4]) {
    v[0]=g; v[1]=ufr32_rotl32(g,5); v[2]=ufr32_rotl32(g,13);
    v[3]=_mm512_xor_si512(g,_mm512_set1_epi32(-1));
}

/* HMix ×1: reduced from HMix ×2 to fund final VAESENC. */
static inline void ufr32_horizontal_mix1(__m512i v[8]) {
    v[0]=_mm512_xor_si512(v[0],ufr32_rotl32(v[4],5));
}

static inline void ufr32_final_diagonal(__m512i v[8]) {
    __m512i t0=v[0],t1=v[1],t2=v[2],t3=v[3];
    __m512i t4=v[4],t5=v[5],t6=v[6],t7=v[7];
    v[0]=_mm512_xor_si512(t0,ufr32_rotl32(t5,7));
    v[1]=_mm512_xor_si512(t1,ufr32_rotl32(t6,11));
    v[2]=_mm512_xor_si512(t2,ufr32_rotl32(t7,13));
    v[3]=_mm512_xor_si512(t3,ufr32_rotl32(t4,17));
    v[4]=_mm512_xor_si512(t4,ufr32_rotl32(t1,19));
    v[5]=_mm512_xor_si512(t5,ufr32_rotl32(t2,23));
    v[6]=_mm512_xor_si512(t6,ufr32_rotl32(t3,29));
    v[7]=_mm512_xor_si512(t7,ufr32_rotl32(t0,31));
}

/* VAESENC: two-stage output path. */
static inline __m512i ufr32_vaes1(__m512i x,__m512i rk) {
    return _mm512_aesenc_epi128(x,rk);
}

/* Minimal online dynamic remix with scratch-register zeroization. */
static inline void ufr32_dynamic_remix(__m512i *a,__m512i *b) {
#if defined(__GNUC__) || defined(__clang__)
    __asm__ volatile(
        "vprold $7,  %[b], %%zmm31\n\t"
        "vpaddd   %%zmm31, %[a], %[a]\n\t"
        "vprold $13, %[a], %%zmm31\n\t"
        "vpxord   %%zmm31, %[b], %[b]\n\t"
        "vpxord   %%zmm31, %%zmm31, %%zmm31\n\t"
        : [a] "+v"(*a), [b] "+v"(*b) : : "zmm31","memory");
#else
    __m512i x=*a,y=*b;
    x=_mm512_add_epi32(x,ufr32_rotl32(y,7));
    y=_mm512_xor_si512(y,ufr32_rotl32(x,13));
    *a=x; *b=y;
#endif
}

static inline void ufr32_one_block(__m512i *a, __m512i *b,
    __m512i vaes_rk, __m512i vaes_rk2, uint32_t *restrict out_words)
{
    const __m512i ia=_mm512_set1_epi32((int)WEYL_INC_32);
    const __m512i ib=_mm512_set1_epi32((int)0xBB67AE85u);

    *a=_mm512_add_epi32(*a,ia);
    *b=_mm512_add_epi32(*b,ib);
    *a=ufr32_gfni(*a);
    *b=ufr32_pshufb_ternlog(*b);

    __m512i v[8];
    ufr32_split4(*a,v+0);
    ufr32_split4(*b,v+4);
    ufr32_horizontal_mix1(v);

    /* Stage 1: per-vector VAES. */
    for(unsigned i=0;i<8;++i) v[i]=ufr32_vaes1(v[i],vaes_rk);

    /* Cross-vector mixing after the first nonlinear stage. */
    ufr32_final_diagonal(v);

    /* Stage 2: per-vector VAES with an independent derived round key. */
    for(unsigned i=0;i<8;++i) v[i]=ufr32_vaes1(v[i],vaes_rk2);

    /* Streaming stores: keep the 512-byte output block out of L1.
       6/2 split lets Dynamic Remix overlap the store stream. */
    for(unsigned i=0;i<6;++i)
        _mm512_stream_si512((__m512i *)(out_words+i*16),v[i]);

    ufr32_dynamic_remix(a,b);

    for(unsigned i=6;i<8;++i)
        _mm512_stream_si512((__m512i *)(out_words+i*16),v[i]);
}

static inline void ufr32_8way_generate(
    UFR32_8Way_CheapSplit_HMix1_VAES1 *r,
    uint32_t *restrict out_words,size_t blocks512)
{
    __m512i a=r->gA,b=r->gB,vaes_rk=r->vaes_rk;
    const __m512i vaes_rk2=_mm512_xor_si512(vaes_rk,
        _mm512_set1_epi32((int)0xA5A5A5A5u));

    for(size_t n=0;n<blocks512;++n) {
        ufr32_one_block(&a,&b,vaes_rk,vaes_rk2,out_words);
        out_words+=128;
    }

    /* Complete the non-temporal store stream once per generation. */
    _mm_sfence();

    r->gA=a; r->gB=b; r->vaes_rk=vaes_rk;
}

/* Secure zeroization. */
static inline void ufr32_secure_zero(void *ptr,size_t len) {
#if defined(__STDC_LIB_EXT1__)
    (void)memset_s(ptr,len,0,len);
#else
    volatile unsigned char *p=(volatile unsigned char *)ptr;
    while(len--) *p++=0;
#endif
}

static inline void ufr32_state_zeroize(
    UFR32_8Way_CheapSplit_HMix1_VAES1 *r) {
    ufr32_secure_zero(r,sizeof(*r));
}

/*
 * KDF boundary: key32 is expected to be an already-derived
 * uniformly distributed 256-bit key.
 */
static inline void ufr32_seed_from_key32(
    UFR32_8Way_CheapSplit_HMix1_VAES1 *r,const uint8_t key32[32])
{
    uint32_t k[8];
    memcpy(k,key32,sizeof(k));

    uint32_t h[8];
    h[0]=k[0]^((k[1]<<7)|(k[1]>>25));
    h[1]=k[1]^((k[2]<<11)|(k[2]>>21));
    h[2]=k[2]^((k[3]<<13)|(k[3]>>19));
    h[3]=k[3]^((k[4]<<17)|(k[4]>>15));
    h[4]=k[4]^((k[5]<<19)|(k[5]>>13));
    h[5]=k[5]^((k[6]<<23)|(k[6]>>9));
    h[6]=k[6]^((k[7]<<29)|(k[7]>>3));
    h[7]=k[7]^((k[0]<<31)|(k[0]>>1));
    __m512i seed_a=_mm512_setr_epi32(
        (int)k[0],(int)k[1],(int)k[2],(int)k[3],
        (int)k[4],(int)k[5],(int)k[6],(int)k[7],
        (int)h[0],(int)h[1],(int)h[2],(int)h[3],
        (int)h[4],(int)h[5],(int)h[6],(int)h[7]);

    uint32_t hb[8];
    hb[0]=k[0]^((k[7]<<7)|(k[7]>>25));
    hb[1]=k[1]^((k[6]<<11)|(k[6]>>21));
    hb[2]=k[2]^((k[5]<<13)|(k[5]>>19));
    hb[3]=k[3]^((k[4]<<17)|(k[4]>>15));
    hb[4]=k[4]^((k[3]<<19)|(k[3]>>13));
    hb[5]=k[5]^((k[2]<<23)|(k[2]>>9));
    hb[6]=k[6]^((k[1]<<29)|(k[1]>>3));
    hb[7]=k[7]^((k[0]<<31)|(k[0]>>1));
    __m512i seed_b=_mm512_setr_epi32(
        (int)k[0],(int)k[1],(int)k[2],(int)k[3],
        (int)k[4],(int)k[5],(int)k[6],(int)k[7],
        (int)hb[0],(int)hb[1],(int)hb[2],(int)hb[3],
        (int)hb[4],(int)hb[5],(int)hb[6],(int)hb[7]);

    ufr32_secure_zero(h,sizeof(h));
    ufr32_secure_zero(hb,sizeof(hb));

    ufr32_seed_mix3(seed_a,seed_b,r);
    ufr32_secure_zero(k,sizeof(k));
}

/* Memory locking. */
#if defined(__linux__) || defined(__APPLE__)
#include <sys/mman.h>
static inline int ufr32_lock_state(UFR32_8Way_CheapSplit_HMix1_VAES1 *r) {
    return mlock((void *)r,sizeof(*r));
}
static inline int ufr32_unlock_state(UFR32_8Way_CheapSplit_HMix1_VAES1 *r) {
    int rc=munlock((void *)r,sizeof(*r));
    ufr32_state_zeroize(r); return rc;
}
#elif defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
static inline int ufr32_lock_state(UFR32_8Way_CheapSplit_HMix1_VAES1 *r) {
    return VirtualLock((LPVOID)r,sizeof(*r)) ? 0 : -1;
}
static inline int ufr32_unlock_state(UFR32_8Way_CheapSplit_HMix1_VAES1 *r) {
    BOOL rc=VirtualUnlock((LPVOID)r,sizeof(*r));
    ufr32_state_zeroize(r); return rc ? 0 : -1;
}
#else
static inline int ufr32_lock_state(UFR32_8Way_CheapSplit_HMix1_VAES1 *r) {
    (void)r; return -1;
}
static inline int ufr32_unlock_state(UFR32_8Way_CheapSplit_HMix1_VAES1 *r) {
    ufr32_state_zeroize(r); return -1;
}
#endif

#endif
