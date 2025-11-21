use sha3::digest::{ExtendableOutput, Update, XofReader};
use crate::err::MlDsaError;
use crate::ntt::{mod_q, ntt, ntt_add, ntt_inverse, ntt_multiply, ntt_neg_vec, ntt_sub, poly_add, poly_sub};
use crate::params::{N, D, ETA, K, L, LEN_ETA_PACK_POLY, LEN_PRIVATE_KEY, LEN_T0_PACK_POLY, SIG_LEN, LAMBDA, TAU, Q, GAMMA2, BITLEN_PACK_W1, GAMMA1, BETA, LEN_Z_SCALAR, OMEGA, LEN_HINT_BIT_PACK};
use crate::xpand::{expand_a, expand_mask};

pub fn sign(sk: &[u8; LEN_PRIVATE_KEY], m: &[u8], ctx: &[u8]) -> Result<[u8; SIG_LEN], MlDsaError> {
    if ctx.len() > 255 {
        return  Err(MlDsaError::SignCtxLenTooLong)
    }
    let mut rnd = [0u8; 32];
    #[cfg(feature="HEDGED")]
    let _ = getrandom::fill(&mut rnd)
        .map_err(|_| MlDsaError::RandomSeedGenError);
    sign_internal(sk, &ctx, &m, &rnd)
}

pub fn sign_internal(sk: &[u8; LEN_PRIVATE_KEY], ctx: &[u8], m: &[u8], rnd: &[u8; 32]) -> Result<[u8; SIG_LEN], MlDsaError> {
    let (c_tilda, z, hint) = response_and_hint(sk, ctx, m, rnd)?;
    Ok(sig_encode(&c_tilda, &z, &hint))
}

pub fn response_and_hint(sk: &[u8; LEN_PRIVATE_KEY], ctx: &[u8], m: &[u8], rnd: &[u8; 32]) -> Result<([u8; LAMBDA/4], [[i32; N]; L], [[i32; N]; K]), MlDsaError> {
    // line 1
    let (rho, key, tr, s1, s2, t0) = sk_decode(sk)?;

    // line 2
    let mut s1_hat = [[0i32; N]; L];
    for (i, w) in s1.iter().enumerate() {
        s1_hat[i] = ntt(&w);
    }

    // line 3
    let mut s2_hat = [[0i32; N]; K];
    for (i, w) in s2.iter().enumerate() {
        s2_hat[i] = ntt(&w);
    }

    // line 4
    let mut t0_hat = [[0i32; N]; K];
    for (i, w) in t0.iter().enumerate() {
        t0_hat[i] = ntt(&w);
    }

    // line 5
    let a_hat = expand_a(&rho);

    // line 6
    let mut mu = [0u8; 64];
    {
        let mut h = sha3::Shake256::default();
        h.update(&tr);
        h.update(&[0, ctx.len() as u8]);
        h.update(ctx);
        h.update(m);
        h.finalize_xof_into(&mut mu);
    }
    // line 7
    let mut rho_dash = [0u8; 64];
    {
        let mut h = sha3::Shake256::default();
        h.update(&key);
        h.update(rnd);
        h.update(&mu);
        h.finalize_xof_into(&mut rho_dash);
    }
    // line 8
    let mut counter = 0usize;

    let mut z = [[0i32; N]; L];
    let mut c_tilda = [0u8; LAMBDA/4];
    let mut hint: [[i32; N]; K];

    // line 9
    let mut response_hint: Option<([u8; LAMBDA/4], [[i32; N]; L], [[i32; N]; K])> = None;
    // line 10
    while response_hint.is_none() {
        if counter >= 2400 { //(1 << 16) {
            return Err(MlDsaError::SignatureAborted);
        }
        // line 11
        let y: [[i32; N]; L] = expand_mask(&rho_dash, counter as u16);
        // line 12
        let mut y_hat = [[0i32; N]; L];
        for (i, r) in y.iter().enumerate() {
            y_hat[i] = ntt(r);
        }

        let mut prod_a_y =  [[0i32; N]; K];
        for k in 0..K {
            for l in 0..L {
                prod_a_y[k] = ntt_multiply(&a_hat[k][l], &y_hat[l]);
            }
        }
        let w: [[i32; N]; K] = prod_a_y.map(|v| ntt_inverse(&v));

        // lines 13 and 14
        // stores high-bits of the polynomial w.
        // coefficients of w1 are in [0, (Q-1)/(2*GAMMA2) - 1]
        let mut w1: [[i32; N]; K] = [[0; N]; K];
        // section 7.4 specifies that functions are applied coefficientwise
        // to the polynomials in the vector,
        for k in 0..K {
            for j in 0..N {
                w1[k][j] = high_bits(w[k][j]);
            }
        }
        // line 15
        // let mut c_tilda = [0u8; LAMBDA/4];
        let mut h = sha3::Shake256::default();
        h.update(&mu);
        h.update(&w1_encode(&w1));
        h.finalize_xof_into(&mut c_tilda);

        // line 16
        let c = sample_in_ball(&c_tilda);

        // line 17
        let c_hat = ntt(&c);

        // line 18: cs1 = ntt_inverse(c_hat * s1_hat)
        let mut cs1 = [[0i32; N]; L];
        for l in 0..L {
            cs1[l] = ntt_inverse(&ntt_multiply(&c_hat, &s1_hat[l]));
        }

        // line 19: cs2 = ntt_inverse(c_hat * s2_hat)
        let mut cs2 =[[0i32; N]; K];
        for k in 0..K {
            cs2[k] = ntt_inverse(&ntt_multiply(&c_hat, &s2_hat[k]));
        }
        // line 20: z = y + cs1
        // let mut z = [[0i32; N]; L];
        for l in 0..L {
            z[l] = ntt_add(&y[l], &cs1[l]);
        }
        // lines 21 and 22.
        let mut r0 = [[0i32; N]; K];
        for k in 0..K {
            r0[k] = ntt_sub(&w[k], &cs2[k]);
            for j in 0..N {
                r0[k][j] = low_bits(r0[k][j]);
            }
        }
        // line 23.
        // https://boringssl.googlesource.com/boringssl/+/main/crypto/fipsmodule/mldsa/mldsa.cc.inc
        // constant_time_select_32(mask, x, y)
        // uint32_t abs_mod_prime(uint32_t x)
        // void scalar_max(uint32_t *max, const scalar *s)
        // uint32_t vector_max(const vector<X> *a)
        // uint32_t z_max = vector_max(&values->sign.z);
        // z mods q
        let infinity_norm_z = infinity_norm(&z);

        // https://boringssl.googlesource.com/boringssl/+/main/crypto/fipsmodule/mldsa/mldsa.cc.inc
        // abs_signed (https://boringssl.googlesource.com/boringssl/+/main/crypto/fipsmodule/mldsa/mldsa.cc.inc#301)
        // scalar_max_signed
        // uint32_t r0_max = vector_max_signed(r0);
        // Returns the absolute value, interpreting the high bit as a sign bit.
        let infinity_norm_r0 = infinity_norm(&r0);
        // line 23.
        let t1 = infinity_norm_z >= (GAMMA1 - BETA);
        let t2 = infinity_norm_r0 >= (GAMMA2 - BETA);
        if (!t1 && !t2) {
            assert!(t1 == false && t2 == false);
        }
        if t1 || t2 {
            #[cfg(feature = "DEBUG_PRINT_RESTARTS")]
            println!("MLDSA signature restart case 1: {t1} {t2}.");
            response_hint = None;
        } else { // lines 24 to 30
            // line 25.
            let mut c_t0 =[[0i32; N]; K];
            for k in 0..K {
                c_t0[k] = ntt_inverse(&ntt_multiply(&c_hat, &t0_hat[k]));
            }
            // lines 26 and 27.
            hint = make_hint(&ntt_neg_vec(&c_t0), &vec_add(&vec_sub(&w, &cs2), &c_t0));
            let infinity_norm_ct0 = infinity_norm(&c_t0);
            // lines 28-30
            let t1 = infinity_norm_ct0 >= GAMMA2;
            let t2 = count_ones(&hint) > OMEGA as u32;
            if  t1 || t2 {
                // #[cfg(feature = "DEBUG_PRINT_RESTARTS")]
                println!("MLDSA signature restart case 2.");
                response_hint = None;
            } else {
                response_hint = Some((c_tilda, z, hint));
                // break;
            }
        }
        // line 31
        counter += L;
    }  // line 32

    println!("Signature found after {counter} restarts!");
    response_hint.ok_or(MlDsaError::SignatureAborted)
}


fn sig_encode(c_tilda: &[u8; LAMBDA/4], z: &[[i32; N]; L], hint: &[[i32; N]; K]) -> [u8; SIG_LEN]{
    let mut sig = [0u8; SIG_LEN];

    sig[..LAMBDA/4].copy_from_slice(c_tilda);
    let mut off = LAMBDA/4;
    for i in 0..L {
        sig[off..][..LEN_Z_SCALAR].copy_from_slice(&bit_pack_z(&z[i]));
        off += LEN_Z_SCALAR;
    }
    sig[off..].copy_from_slice(&bit_pack_hint(hint));
    sig
}

fn infinity_norm<const D: usize>(vec: &[[i32; N]; D]) -> u32 {
    const HALF_Q: u32 = (Q / 2) as u32;
    let mut norm: u32 = 0;
    for k in 0..D {
        for &c in vec[k].iter() {
            assert!(c < Q);
            // return x <= Q/2 ? x : kPrime - x;
            let abs_mod_prime = if c as u32 <= HALF_Q {
                c as u32
            } else {
                (Q - c) as u32
            };
            if abs_mod_prime > norm {
                norm = abs_mod_prime;
            }
        }
    }
    norm
}

fn make_hint(z: &[[i32; N]; K], r: &[[i32; N]; K]) -> [[i32; N]; K] {
    let mut hint = [[0i32; N]; K];
    for k in 0..K {
        for i in 0..N {
            hint[k][i] = if (high_bits(r[k][i]) != high_bits(r[k][i] + z[k][i])) {
                1
            } else {
                0
            }
        }
    }
    hint
}

fn count_ones(h: &[[i32; N]; K]) -> u32 {
    let mut ones = 0u32;
    for k in 0..K {
        for i in 0..N {
            assert!(h[k][i] == 0 || h[k][i] == 1);
            ones += h[k][i] as u32;
        }
    }
    ones
}

pub fn vec_add(x: &[[i32; N]; K], y: &[[i32; N]; K]) -> [[i32; N]; K] {
    let mut r = [[0i32; N]; K];
    for k in 0..K {
        r[0] = ntt_add(&x[k], &y[k]);
    }
    r
}

pub fn vec_sub(x: &[[i32; N]; K], y: &[[i32; N]; K]) -> [[i32; N]; K] {
    let mut r = [[0i32; N]; K];
    for k in 0..K {
        r[0] = ntt_sub(&x[k], &y[k]);
    }
    r
}

// algorithm 28, page 35.
// encodes a polynomial vector w1 into a byte string.
// input: w1 \in R^K whose polynomial coefficients are in [0, (Q-1)/(2*GAMMA2) - 1]
fn w1_encode(w1: &[[i32; N]; K]) -> [u8; 32*K*BITLEN_PACK_W1]{
    // K polynomials. A polynomial has (32*8) coefficients, each BITLEN_PACK_W1 wide.
    let r = [0u8; 32*K*BITLEN_PACK_W1];
    for i in 0..K {
        simple_bit_pack_w1(&w1[i]);
    }
    r
}

// simple_bit_pack encodes a polynomial into a byte string.
// simple_bit_pack_w1 is specifically tailored to pack public key component w1's
//   coefficients into groups of 4 bits (ML_DSA_65 and ML_DSA_87) or 6 bits (ML_DSA_44).
// input: B = 1023 and w= \in R such that coefficients are all in [0, 1023].
// output: a byte string of length 32 * bitlen(B)
#[inline]
fn simple_bit_pack_w1(w1: &[i32; N]) -> [u8; 32 * BITLEN_PACK_W1] {

    let mut r = [0u8; 32 * BITLEN_PACK_W1];

    #[cfg(feature="ML_DSA_44")]
    for i in 0..N/4 {
        r[3*i+0]  = w1[4*i+0] as u8;
        r[3*i+0] |= (w1[4*i+1] as u8) << 6;
        r[3*i+1]  = (w1[4*i+1] as u8) >> 2;
        r[3*i+1] |= (w1[4*i+2] as u8) << 4;
        r[3*i+2]  = (w1[4*i+2] as u8) >> 4;
        r[3*i+2] |= (w1[4*i+3] as u8) << 2;
    }

    #[cfg(any(feature="ML_DSA_65", feature="ML_DSA_87"))]
    for i in 0..N/2 {
        r[i] = (w1[2 * i + 0] as u8) | ((w1[2 * i + 1] as u8) << 4);
    }

    r
}

fn bit_pack_hint(hint: &[[i32; N]; K]) -> [u8; LEN_HINT_BIT_PACK] {
    let mut y = [0u8; LEN_HINT_BIT_PACK];
    let mut index = 0;
    for i in 0..K {
        for j in 0..N {
            if hint[i][j] != 0 {
                y[index] = j as u8;
                index += 1;
            }
        }
        y[OMEGA + i] = index as u8;
    }
    y
}

fn bit_pack_z(z: &[i32; N]) -> [u8; LEN_Z_SCALAR] {
    let mut r = [0u8; LEN_Z_SCALAR];
    let mut t = [0i32; 4];

    #[cfg(feature="ML_DSA_44")]
    for i in 0..N/4 {
        t[0] = GAMMA1 as i32 - z[4*i+0];
        t[1] = GAMMA1 as i32 - z[4*i+1];
        t[2] = GAMMA1 as i32 - z[4*i+2];
        t[3] = GAMMA1 as i32 - z[4*i+3];

        r[9*i+0]  = t[0] as u8;
        r[9*i+1]  = (t[0] >> 8) as u8;
        r[9*i+2]  = (t[0] >> 16) as u8;
        r[9*i+2] |= (t[1] << 2) as u8;
        r[9*i+3]  = (t[1] >> 6) as u8;
        r[9*i+4]  = (t[1] >> 14) as u8;
        r[9*i+4] |= (t[2] << 4) as u8;
        r[9*i+5]  = (t[2] >> 4) as u8;
        r[9*i+6]  = (t[2] >> 12) as u8;
        r[9*i+6] |= (t[3] << 6) as u8;
        r[9*i+7]  = (t[3] >> 2) as u8;
        r[9*i+8]  = (t[3] >> 10) as u8;
    }

    #[cfg(any(feature="ML_DSA_65", feature="ML_DSA_87"))]
    for i in 0..N/2 {
        t[0] = GAMMA1 as i32 - z[2*i+0];
        t[1] = GAMMA1 as i32 - z[2*i+1];

        r[5*i+0]  = t[0] as u8;
        r[5*i+1]  = (t[0] >> 8) as u8;
        r[5*i+2]  = (t[0] >> 16) as u8;
        r[5*i+2] |= (t[1] << 4) as u8;
        r[5*i+3]  = (t[1] >> 4) as u8;
        r[5*i+4]  = (t[1] >> 12) as u8;
    }

    r
}
// algorithm 37, page 40.
// returns r1 from the output of decompose(r)
// input: r \in Z_q
// output: integer r1.
#[inline]
fn high_bits(r: i32) -> i32 {
    let (r1, r0) = decompose(r);
    assert!(r0 < Q && r1 < Q);
    r1
}

// algorithm 38, page 41.
// returns r0 from the output of decompose(r)
// input: r \in Z_q
// output: integer r0.
#[inline]
fn low_bits(r: i32) -> i32 {
    let (r1, r0) = decompose(r);
    assert!(r0 < Q && r1 < Q);
    r0
}

// decomposes r into (r1, r0) such that r = r1*(2*GAMMA2) + r0 mod q
// input: r \in Z_q.
// output: integers (r1, r0)
// r1 in [0, (Q-1)/(2*GAMMA2) - 1]; r1 will have at most 44 (ML_DSA_44) or 16 values.

#[inline]
fn decompose(r: i32) -> (i32, i32) {
    let mut r1 = (r + 127) >> 7;

    #[cfg(feature = "ML_DSA_44")]
    {
        // 1488/2^24 is close enough to 1/1488 so that r1 becomes x/(2 gamma2) rounded down.
        r1 = (r1 * 11275 + (1 << 23)) >> 24;
        // For corner-case r1 = (Q-1)/(2 gamma2) = 44, we have to set r1=0.
        r1 ^= ((43 - r1) >> 31) & r1;
        assert!((r1 >= -44) && (r1 <= 44));
    }

    #[cfg(any(feature = "ML_DSA_65", feature = "ML_DSA_87"))]
    {
        // Here, Gamma2 == 2^18 - 2^8.
        // This set r1 = ((ceil(x / 2^7) * (2^10 + 1) + 2^21) / 2^22) mod 2^4
        r1 = ((r1 * 1025 + (1 << 21)) >> 22) & 15;
        assert!((r1 >= 0) && (r1 < 16));
    }
    let r0 = r as i64 - r1 as i64 * 2 * GAMMA2 as i64;
    let r0 = r0 - (((Q - 1) as i64/2 - r0) >> 31) & Q as i64;

    (r1, mod_q(r0))
}

// samples the challenge polynomial with TAU nonzero coefficients in {-1,1}.
// returns a polynomial c: R with coefficients in {-1, 0, 1} and Hamming weight t <= 64.
// That is, the number of non-zero coefficients in c are less than 64.
// More precisely, it will be one of 39, 49, and 60 for category 2, 3, and 5, respectively.
#[inline]
pub fn sample_in_ball(rho: &[u8; LAMBDA/4]) -> [i32; N] {
    let mut h = sha3::Shake256::default();
    h.update(rho);
    let mut xof = h.finalize_xof();
    let mut s = [0u8; 8];
    xof.read(&mut s);
    let mut h: u64 = u64::from_le_bytes(s); // note: h has N bits.
    let mut c = [0i32; N];
    for i in (N-TAU)..256usize { // iterate TAU times
        assert!(i > 195); // or N-i < 61
        let mut _j = [0u8];
        xof.read(&mut _j);
        let mut j = _j[0] as usize;
        while j > i {
            xof.read(&mut _j);
            j = _j[0] as usize
        }
        // randomly selected coefficient c_j has a smaller (or same) degree than c_i
        c[i] = c[j]; // init c[i] with {-1, 0, 1}, from a random location; most likely contains 0.
        // h[i + TAU - N] actually evaluates to h[0], [1], h[2], ..., h[255]
        // (N-TAU) + TAU - N, (N-TAU + 1) + TAU - N, (N-TAU + 2) + TAU - N,..., (N-TAU + 255) + TAU - N,
        c[j] = 1 - 2 * (h & 1) as i32; // c_j is -1 if (h&1) is 1, and 1 otherwise.
        h >>= 1
    }
    c
}

pub fn sk_decode(sk: &[u8; LEN_PRIVATE_KEY]) -> Result<([u8; 32], [u8; 32], [u8; 64], [[i32; N]; L], [[i32; N]; K], [[i32; N]; K]), MlDsaError> {
    assert_eq!(sk.len(), 128 + L * LEN_ETA_PACK_POLY + K * LEN_ETA_PACK_POLY + K * LEN_T0_PACK_POLY);
    let rho = sk[0..32].try_into()?;
    let key = sk[32..64].try_into()?;
    let tr = sk[64..128].try_into()?;
    let y: &[u8; L * LEN_ETA_PACK_POLY] = &sk[128..128 + L*LEN_ETA_PACK_POLY].try_into()?;
    let z: &[u8; K * LEN_ETA_PACK_POLY] = &sk[128 + L*LEN_ETA_PACK_POLY..128 + L*LEN_ETA_PACK_POLY + K*LEN_ETA_PACK_POLY].try_into()?;
    let w: &[u8; K * LEN_T0_PACK_POLY] = &sk[128 + L*LEN_ETA_PACK_POLY + K*LEN_ETA_PACK_POLY ..].try_into()?;

    let mut s1 = [[0i32; N]; L];
    for i in 0..L {
        let t = &y[i * LEN_ETA_PACK_POLY..(i + 1) * LEN_ETA_PACK_POLY].try_into()?;
        bit_unpack_eta(&t, &mut s1[i])?;
    }

    let mut s2: [[i32; N]; K] = [[0; N]; K];
    for i in 0..K {
        let t = &z[i * LEN_ETA_PACK_POLY..(i + 1) * LEN_ETA_PACK_POLY].try_into()?;
        bit_unpack_eta(&t, &mut s2[i])?;
    }

    let mut t0: [[i32; N]; K] = [[0; N]; K];
    for i in 0..K {
        let t = &w[i * LEN_T0_PACK_POLY..(i + 1) * LEN_T0_PACK_POLY].try_into()?;
        bit_unpack_t0(&t, &mut t0[i])?;
    }
    Ok((rho, key, tr, s1, s2, t0))
}

fn inclusive(l: i32, h: i32, v: i32) -> bool {
    (l..h+1).contains(&v)
}

// w is bit-packed polynomial
fn bit_unpack_eta(w: &[u8; LEN_ETA_PACK_POLY], s: &mut[i32; N]) -> Result<(), MlDsaError> {
    const _ETA_: i32 = ETA as i32;
    let mut ok = true;

    #[cfg(any(feature="ML_DSA_44", feature="ML_DSA_87"))]
    for i in 0..N/8 {
        s[8*i+0] = ((w[3*i+0] >> 0) & 7) as i32;
        s[8*i+1] = ((w[3*i+0] >> 3) & 7) as i32;
        s[8*i+2] = (((w[3*i+0] >> 6) | (w[3*i+1] << 2)) & 7) as i32;
        s[8*i+3] = ((w[3*i+1] >> 1) & 7) as i32;
        s[8*i+4] = ((w[3*i+1] >> 4) & 7) as i32;
        s[8*i+5] = (((w[3*i+1] >> 7) | (w[3*i+2] << 1)) & 7) as i32;
        s[8*i+6] = ((w[3*i+2] >> 2) & 7) as i32;
        s[8*i+7] = ((w[3*i+2] >> 5) & 7) as i32;

        s[8*i+0] = _ETA_ - s[8*i+0];
        s[8*i+1] = _ETA_ - s[8*i+1];
        s[8*i+2] = _ETA_ - s[8*i+2];
        s[8*i+3] = _ETA_ - s[8*i+3];
        s[8*i+4] = _ETA_ - s[8*i+4];
        s[8*i+5] = _ETA_ - s[8*i+5];
        s[8*i+6] = _ETA_ - s[8*i+6];
        s[8*i+7] = _ETA_ - s[8*i+7];
        ok &= inclusive(-_ETA_, _ETA_, s[2*i+1]) & inclusive(-_ETA_, _ETA_, s[2*i+1]) &
            inclusive(-_ETA_, _ETA_, s[2*i+2]) & inclusive(-_ETA_, _ETA_, s[2*i+1]) &
            inclusive(-_ETA_, _ETA_, s[2*i+5]) & inclusive(-_ETA_, _ETA_, s[2*i+1]) &
            inclusive(-_ETA_, _ETA_, s[2*i+7]) & inclusive(-_ETA_, _ETA_, s[2*i+1]);
    }

    #[cfg(feature="ML_DSA_65")]
    for i in 0..N/2 {
        s[2*i+0] = (w[i] & 0x0F) as i32;
        s[2*i+1] = (w[i] >> 4) as i32;
        s[2*i+0] = _ETA_ - s[2 * i + 0];
        s[2*i+1] = _ETA_ - s[2 * i + 1];
        ok &= inclusive(-_ETA_, _ETA_, s[2 * i + 0]) & inclusive(-_ETA_, _ETA_, s[2 * i + 1]);
    }

    if ok {
        Ok(())
    } else {
        Err(MlDsaError::MalformedShortVector)
    }
}

/**
One iteration handles bits in this fashion - w0 is the first byte to consider in the loop.

|         w3         |              w2             |          w1         |        w0           |
|------+------+------|+------+------+------+-------|+------+------+------|+------+------+------|
|             |  8   |   5   |  13  |  13  |   1   |   12  |  13  |   7  |   6   |  13  |  13  |
|------+------+------|+------+------+------+-------|+------+------+------|+------+------+------|
|                    |                             |                     |                     |
 **/
// Unpack polynomial t0 with coefficients in ]-2^{D-1}, 2^{D-1}]
fn bit_unpack_t0(w: &[u8; LEN_T0_PACK_POLY], t: &mut [i32; N]) -> Result<(), MlDsaError> {
    const _2_POW_D_MINUS_1_: i32 = 1 << (D-1);
    let mut ok = true;
    for i in 0..N/8 {
        t[8*i+0] = w[13*i+0] as i32;
        t[8*i+0] |= (w[13*i+1] as i32) << 8;
        t[8*i+0] &= 0x1FFF;

        t[8*i+1] = (w[13*i+1] as i32) >> 5;
        t[8*i+1] |= (w[13*i+2] as i32) << 3;
        t[8*i+1] |= (w[13*i+3] as i32) << 11;
        t[8*i+1] &= 0x1FFF;

        // w0 = | 6 | 13 | 13 |
        let w0 = i32::from_le_bytes(w[13*i..13*i+4].try_into()?);
        assert_eq!(w0 & 0x1FFF, t[8*i+0]);
        assert_eq!((w0 >> 13) & 0x1FFF, t[8*i+1]);
        let w0 = (w0 >> 26) & 0x3F;

        t[8*i+2] = (w[13*i+3] as i32) >> 2;
        t[8*i+2] |= (w[13*i+4] as i32) << 6;
        t[8*i+2] &= 0x1FFF;

        t[8*i+3] = (w[13*i+4] as i32) >> 7;
        t[8*i+3] |= (w[13*i+5] as i32) << 1;
        t[8*i+3] |= (w[13*i+6] as i32) << 9;
        t[8*i+3] &= 0x1FFF;

        //                w0 = |  6 | 13 | 13 |
        // w1 = | 12 | 13 |  7 |
        let w1 = i32::from_le_bytes(w[13*i+4..13*i+8].try_into()?);
        assert_eq!(((w1 & 0x7F) << 6) | w0, t[8*i+2]);
        assert_eq!((w1 >> 7) & 0x1FFF, t[8*i+3]);
        let w1 = (w1 >> 20) & 0xFFF;

        t[8*i+4] = (w[13*i+6] as i32) >> 4;
        t[8*i+4] |= (w[13*i+7] as i32) << 4;
        t[8*i+4] |= (w[13*i+8] as i32) << 12;
        t[8*i+4] &= 0x1FFF;

        t[8*i+5] = (w[13*i+8] as i32) >> 1;
        t[8*i+5] |= (w[13*i+9] as i32) << 7;
        t[8*i+5] &= 0x1FFF;

        //                                    w0 = |  6 | 13 | 13 |
        //                     w1 = | 12 | 13 |  7 |
        // w2 = |  5 | 13 | 13 |  1 |
        let w2 = i32::from_le_bytes(w[13*i+8..13*i+12].try_into()?);
        assert_eq!(((w2 & 0x1) << 12) | w1, t[8*i+4]);
        assert_eq!((w2 >> 1) & 0x1FFF, t[8*i+5]);
        let w2 = w2 >> 14; // w2 = | 0 | 5 |  13 |

        t[8*i+6] = (w[13*i+9] as i32) >> 6;
        t[8*i+6] |= (w[13*i+10] as i32) << 2;
        t[8*i+6] |= (w[13*i+11] as i32) << 10;
        t[8*i+6] &= 0x1FFF;

        t[8*i+7] = (w[13*i+11] as i32) >> 3;
        t[8*i+7] |= (w[13*i+12] as i32) << 5;
        t[8*i+7] &= 0x1FFF;

        //                w2 = | 5 | 13 |
        // w3 = | 11 | 13 |  8 |
        // we only use the least-significant byte of w3, and finish one iteration.
        assert_eq!(w2 & 0x1FFF, t[8 * i + 6]);
        let w2 = (w2 >> 13) & 0x1F; // w2 = | 0 | 5 | 13 | >> 13 = | 0 | 5 |
        assert_eq!(((w[13 * i + 12] as i32) << 5) | (w2 & 0x1F), t[8*i+7]);

        t[8*i+0] = _2_POW_D_MINUS_1_ - t[8*i+0];
        t[8*i+1] = _2_POW_D_MINUS_1_ - t[8*i+1];
        t[8*i+2] = _2_POW_D_MINUS_1_ - t[8*i+2];
        t[8*i+3] = _2_POW_D_MINUS_1_ - t[8*i+3];
        t[8*i+4] = _2_POW_D_MINUS_1_ - t[8*i+4];
        t[8*i+5] = _2_POW_D_MINUS_1_ - t[8*i+5];
        t[8*i+6] = _2_POW_D_MINUS_1_ - t[8*i+6];
        t[8*i+7] = _2_POW_D_MINUS_1_ - t[8*i+7];
    }
    for i in 0..N {
        ok &= inclusive(-_2_POW_D_MINUS_1_+1, _2_POW_D_MINUS_1_, t[i]);
    }
    if ok {
        Ok(())
    } else {
        Err(MlDsaError::MalformedShortVector)
    }
}