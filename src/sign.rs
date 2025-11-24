use sha3::digest::{ExtendableOutput, Update, XofReader};
use crate::err::MlDsaError;
use crate::ntt::{mod_q, ntt, ntt_add, ntt_inverse, ntt_multiply, ntt_neg_vec, ntt_sub, poly_add, poly_sub};
use crate::params::{N, D, ETA, K, L, LEN_ETA_PACK_POLY, LEN_PRIVATE_KEY, LEN_T0_PACK_POLY, SIG_LEN, LAMBDA, TAU, Q, GAMMA2, BITLEN_PACK_W1, GAMMA1, BETA, LEN_Z_SCALAR, OMEGA, LEN_HINT_BIT_PACK};
use crate::xpand::{expand_a, expand_mask};

pub fn sign(sk: &[u8; LEN_PRIVATE_KEY], m: &[u8]) -> Result<[u8; SIG_LEN], MlDsaError> {
    //if ctx.len() > 255 {
    //    return  Err(MlDsaError::SignCtxLenTooLong)
    //}
    #[cfg(feature="HEDGED")]
    {
        let mut rnd = [0u8; 32];
        getrandom::fill(&mut rnd)
            .map_err(|_| MlDsaError::RandomSeedGenError)?;
        sign_random(sk, &rnd, &m, &ctx)
    }

    #[cfg(feature="DETERMINISTIC")]
    sign_random(sk, &[0u8; 32], &m)
}

pub fn sign_random(sk: &[u8; LEN_PRIVATE_KEY], rnd: &[u8; 32], m: &[u8]) -> Result<[u8; SIG_LEN], MlDsaError> {
    sign_internal(sk, &m, &rnd)
}

pub fn sign_internal(sk: &[u8; LEN_PRIVATE_KEY], m: &[u8], rnd: &[u8; 32]) -> Result<[u8; SIG_LEN], MlDsaError> {
    let (c_tilda, z, hint) = response_and_hint(sk, m, rnd)?;
    Ok(sig_encode(&c_tilda, &z, &hint))
}

pub fn response_and_hint(sk: &[u8; LEN_PRIVATE_KEY], m: &[u8], rnd: &[u8; 32]) -> Result<([u8; LAMBDA/4], [[i32; N]; L], [[i32; N]; K]), MlDsaError> {
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
        //h.update(&[0, ctx.len() as u8]);
        //h.update(&[0, 0]);
        // h.update(ctx);
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

    // line 9
    let mut response_hint: Option<([u8; LAMBDA/4], [[i32; N]; L], [[i32; N]; K])> = None;
    // line 10
    while response_hint.is_none() {
        // let mut c_tilda = [0u8; LAMBDA/4];
        let mut z = [[0i32; N]; L];
        let mut hint = [[2i32; N]; K];

        if counter/4 > 814 { // it's not worth waiting :-)
            return Err(MlDsaError::SignatureAborted);
        }
        // line 11
        let y: [[i32; N]; L] = expand_mask(&rho_dash, counter as u16);
        // line 12
        let mut y_hat = [[0i32; N]; L];
        for l in 0..L {
            y_hat[l].copy_from_slice(&ntt(&y[l]));
        }

        let mut w = [[0i32; N]; K];
        for k in 0..K {
            let mut t = [0i32; N];
            for l in 0..L {
                t = ntt_add(&t, &ntt_multiply(&a_hat[k][l], &y_hat[l]));
            }
            w[k].copy_from_slice(&ntt_inverse(&t));
        }

        // lines 13 and 14
        // stores high-bits of the polynomial w.
        // coefficients of w1 are in [0, (Q-1)/(2*GAMMA2) - 1]
        let mut w1 = [[0i32; N]; K];
        // section 7.4 specifies that functions are applied coefficientwise
        // to the polynomials in the vector,
        for k in 0..K {
            for j in 0..N {
                w1[k][j] = high_bits(w[k][j]);
            }
        }
        // line 15
        let mut c_tilda = [0u8; LAMBDA/4];
        {
            let mut h = sha3::Shake256::default();
            h.update(&mu);
            h.update(&w1_encode(&w1));
            h.finalize_xof_into(&mut c_tilda);
        }

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
        let mut cs2 = [[0i32; N]; K];
        for k in 0..K {
            cs2[k].copy_from_slice(&ntt_inverse(&ntt_multiply(&c_hat, &s2_hat[k])));
        }
        // line 20: z = y + cs1
        // let mut z = [[0i32; N]; L];
        for l in 0..L {
            z[l] = poly_add(&y[l], &cs1[l]);
        }
        // lines 21 and 22.
        let mut r0 = [[0i32; N]; K];
        for k in 0..K {
            r0[k].copy_from_slice(&poly_sub(&w[k], &cs2[k]));
        }
        for k in 0..K {
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
        let t1 = infinity_norm_z >= (GAMMA1 - BETA) as i32;
        let t2 = infinity_norm_r0 >= (GAMMA2 - BETA) as i32;
        if t1 || t2 {
            // #[cfg(feature = "DEBUG_PRINT_RESTARTS")]
            println!("MLDSA signature restart case 1: {t1} {t2}.");
            response_hint = None;
        } else { // lines 24 to 30
            // line 25.
            let mut c_t0 =[[0i32; N]; K];
            for k in 0..K {
                c_t0[k] = ntt_inverse(&ntt_multiply(&c_hat, &t0_hat[k]));
            }
            // lines 26 and 27.
            // make_hint(&ntt_neg_vec(&c_t0), &vec_add(&vec_sub(&w, &cs2), &c_t0), &mut hint);
            make_hint(&ntt_neg_vec(&c_t0), &vec_add(&vec_sub(&w, &cs2), &c_t0), &mut hint);
            let infinity_norm_ct0 = infinity_norm(&c_t0);
            println!("infinity_norm_ct0 = {infinity_norm_ct0}");
            // lines 28-30
            let t1 = infinity_norm_ct0 >= GAMMA2 as i32;
            let t2 = count_ones(&hint) > OMEGA as u32;
            println!("count_ones = {}", count_ones(&hint));
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
    for i in 0..L {
        let off = (LAMBDA/4) + (i*LEN_Z_SCALAR);
        sig[off..off+LEN_Z_SCALAR].copy_from_slice(&bit_pack_z(&z[i]));
    }
    sig[SIG_LEN-LEN_HINT_BIT_PACK..].copy_from_slice(&bit_pack_hint(hint));
    sig
}

fn infinity_norm<const D: usize>(vec: &[[i32; N]; D]) -> i32 {
    const HALF_Q: i32 = (Q)/2;
    let mut norm = 0;
    for k in 0..D {
        for &c in vec[k].iter() {
            // assert!(c < Q);
            // return x <= Q/2 ? x : kPrime - x;
            let abs_mod_prime = if c <= HALF_Q {
                c.abs()
            } else {
                (c-Q).abs()
            };
            if abs_mod_prime > norm {
                norm = abs_mod_prime;
            }
        }
    }
    norm
}

fn make_hint(z: &[[i32; N]; K], r: &[[i32; N]; K], hint: &mut [[i32; N]; K]) {
    for k in 0..K {
        for i in 0..N {
            hint[k][i] = if high_bits(r[k][i]) != high_bits(r[k][i] + z[k][i]) {
                1
            } else {
                0
            }
        }
    }
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
        r[k].copy_from_slice(&ntt_add(&x[k], &y[k]));
    }
    r
}

pub fn vec_sub(x: &[[i32; N]; K], y: &[[i32; N]; K]) -> [[i32; N]; K] {
    let mut r = [[0i32; N]; K];
    for k in 0..K {
        r[k].copy_from_slice(&ntt_sub(&x[k], &y[k]));
    }
    r
}

// algorithm 28, page 35.
// encodes a polynomial vector w1 into a byte string.
// input: w1 \in R^K whose polynomial coefficients are in [0, (Q-1)/(2*GAMMA2) - 1]
fn w1_encode(w1: &[[i32; N]; K]) -> [u8; 32*K*BITLEN_PACK_W1]{
    // K polynomials. A polynomial has (32*8) coefficients, each BITLEN_PACK_W1 wide.
    let mut r = [0u8; 32*K*BITLEN_PACK_W1];
    let mut off = 0;
    for i in 0..K {
        r[off..][..32 * BITLEN_PACK_W1].copy_from_slice(&simple_bit_pack_w1(&w1[i]));
        off += 32 * BITLEN_PACK_W1;
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
    {
        for i in 0..N/2 {
            r[i] = (w1[2 * i + 0] as u8) | ((w1[2 * i + 1] as u8) << 4);
        }
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
    println!("bit_pack_hint {y:?}");
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
    let (r1, _r0) = decompose(r);
    // let (_r1, _r0) = _decompose_(r);
    // assert_eq!(r1, _r1);
    //assert_eq!(mod_q((r1 * (2*GAMMA2 as i32))  as i64) + mod_q(r0.into()), mod_q(r as i64));
    // assert_eq!(mod_q(((_r1 * (2*GAMMA2 as i32)) + _r0) as i64), mod_q(r as i64));
    r1
}

// algorithm 38, page 41.
// returns r0 from the output of decompose(r)
// input: r \in Z_q
// output: integer r0.
#[inline]
fn low_bits(r: i32) -> i32 {
    let (_r1, r0) = decompose(r);
    r0
}

// r = rx mod Q
// r0 <- r mods(2*GAMMA2)
fn mods_2_gamma2(r: i32) -> i32 {
    const ALPHA: i32 = 2 * GAMMA2 as i32;
    let mut r0 = r % ALPHA;
    if r0 > ALPHA / 2 {
        r0 -= ALPHA;
    }
    r0
}

fn decompose(r: i32) -> (i32, i32) {
    const ALPHA: i32 = 2 * GAMMA2 as i32;
    let rp = mod_q(r as i64);
    let r0 = mods_2_gamma2(rp);

    let (r1, r0) = if rp - r0 == Q - 1 {
        (0, r0 - 1)
    } else {
        ((rp-r0)/ALPHA, r0)
    };

    debug_assert_eq!(
        mod_q(r as i64),
        mod_q((r1*ALPHA + r0) as i64)
    );

    (r1, r0)
}
fn __decompose(rx: i32) -> (i32, i32) {
    let rp = mod_q(rx as i64);
    let r0 = if rp <= GAMMA2 as i32 { // rp mods (2*GAMMA2)
        rp
    } else {
        rp - (2*GAMMA2 as i32)
    };
    // println!("decompose_44 (1) - {rx} == ({rp}, {r0})");
    let (r1, r0) = if (rp - r0) == Q-1 {
        (0, mod_q(r0 as i64 - 1))
        // (0, r0 - 1)
    } else {
        ((rp - r0)/(2*GAMMA2 as i32), r0)
    };
    // println!("decompose_44 (2) - ({r1}, {r0})");
    assert_eq!(mod_q(rx as i64), mod_q((r1*2*GAMMA2 as i32 + r0) as i64));
    // _decompose_(rx);
    (r1, r0)
}

// decomposes r into (r1, r0) such that r = r1*(2*GAMMA2) + r0 mod q
// input: r \in Z_q.
// output: integers (r1, r0)
// r1 in [0, (Q-1)/(2*GAMMA2) - 1]; r1 will have at most 44 (ML_DSA_44) or 16 values.

#[inline]
fn _decompose_(rx: i32) -> (i32, i32) {
    let r = mod_q(rx as i64);
    let mut r1 = (r + 127) >> 7;

    #[cfg(feature = "ML_DSA_44")]
    {
        // 1488/2^24 is close enough to 1/1488 so that r1 becomes x/(2 gamma2) rounded down.
        r1 = (r1 * 11275 + (1 << 23)) >> 24;
        // For corner-case r1 = (Q-1)/(2 gamma2) = 44, we have to set r1=0.
        r1 ^= ((43 - r1) >> 31) & r1;
        assert!((r1 >= 0) && (r1 <= 43));
    }

    #[cfg(any(feature = "ML_DSA_65", feature = "ML_DSA_87"))]
    {
        // Here, Gamma2 == 2^18 - 2^8.
        // This set r1 = ((ceil(x / 2^7) * (2^10 + 1) + 2^21) / 2^22) mod 2^4
        r1 = ((r1 * 1025 + (1 << 21)) >> 22) & 15;
        assert!((r1 >= 0) && (r1 < 16));
    }
    let mut r0 = r as i64 - (r1 as i64 * 2 * GAMMA2 as i64);
    r0 = r0 - ((((Q - 1) as i64/2 - r0) >> 31) & Q as i64);

    assert_eq!(mod_q((r1 as i64 * (2*GAMMA2) as i64) + r0), mod_q(rx as i64));

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

#[cfg(test)]
mod mldsa_signature_tests {
    use crate::{params, sign};
    use crate::params::{LAMBDA, LEN_PRIVATE_KEY};


    // https://github.com/post-quantum-cryptography/KAT/blob/main/MLDSA/kat_MLDSA_44_det_raw.rsp
    #[cfg(feature="ML_DSA_44")]
    #[test]
    fn kat3_mlen16() {
        let _sk = "bd4e96f9a038ab5e36214fe69c0b1cb835ef9d7c8417e76aecd152f5cddebec8e9355a9ab1bcfbfc7d4decdf5ddd94b4359823a9578a19d161c4e9ad5afd7d61746119a59589b3630773bfd350c45d029f1cd0b03163feac00f77d814cb89d963b11a8dfad7d0f31b0a542ba61044a2968d6cc435f00f5330f65063f5a46db2323b1691c044d0cb28c2487014a9230dab8510802115848401928228a468819070c5c002a0b18109b3632493810c404428498048444259a323024a4259c062cd9344454000018c58d64142113888020400c02a220811484823685d2086010a54052b231a4a48c64000c11048522868954c609cab6845346898b101122038044b0915ca00952a2288bb488612880120580c948668124925906429b280ed13005019631dbb804e2c63011024ccb406549088e0a470d53c08919b8711091651a3112ca125080242e24134a22014c1c48100a274e03354801c0819098845908851489115b400ee09885d93462c48440e4a8901b418d444031d0b420008665a2b60464964d4b1802a130660a020ee4862dd2162a9c186103220c5b3228980445d93465cb18521b11444c8084124470dca201920851100912a4880c109708c2b0290ac048e2260cdb344d9b226e030951149160c0042a1819884202665802321a366ed3127298b82914034d90b200a1c6841c121102c6015004525842685a160c03204603a40d04480c52426821271291902cd92489cb088de30481a18270c1b8280848310c499023079049b04c883422118009138225189204c240061ac18422b270e336904bc60081b83152a04923c3215206318b2605cc306a58162664108053227289348cc1288c1a118a98406008098dc9246ee4840054140e4ba00543826da0188d03922d10852519072c6440908a107064c66c51308219b56122b420540866d41641089828c3160a60286841002e492206502205ca262ddc380582203004c16562a660da0620098510888844e1c011432466c394711c47324c28125904625b8049db40211915490b3184e1206600c588c4264542c871c4263120050549127049b66dd4b66884c830c44406c1826123874dc9127018116e8b10915c0291cb8485d032620a8160423820c91820003544982886020331c120818030280c0551cc86890022065c388e24264da42888188701c2c06dd4368904b16404278ada448120c128a4a80114350e530031624651d9126d0b88040b856904c96d429050808881004b402ec8f9a63ad7db8f6c580a15254857376b862730d93f80c4976e1ef1947efa8513adc92867a377e74969510fd7417fa21bd3d773f60d20e4b2ba56cc54f7671e087c16d996c25ec97259cf8e0a4a5c823bcb4ee3ab5aa58ece39111e1513c30cb90a6453d09d79e35e7bc77a3875efea5a2f0f8e070dd14cf4fcb75f13d61ee24dbc6b9977fe0c44f17166b06d4e2425524e4ee455ced586033641c2c8ba5e9c0a8317ab418233d9ee411ef2bb013bb1d4b12d1ef8a82bc591218c3e5f7ffeb2f7b95a643b8cc46e49fdbb83656be0167b0afd3ae9205c1b02694f7df939eea45529189cd8c6ac518cdfdd0e7d0a04e521ad366d99b65b9f389dc774b15168907fba80ad2e6f540f8c36d81407f12d3dbb83ba0ef403a9db8453374c2b65aa76fb31ce475f4a857e4d155e929a77086dd169c0657c9306465910dc15966a9cde2923d9c36ecf9a82fb6ce8ad6e3362f2dfe4a6b6ff54ef1bda398d3ae032590e81923b1c5b287ef18e6d480f06d4f57552c6c74e329eb7f1326f6638935afc9a87d08147f6a3c081dbe4b7660d3bf7cc2a86bfb22cdc2a3932ecd9b20f770088c5d4a409950c3d91b5c0a87017f394af3d26d083986355553b9875720204e03eb2a8a1abc9f486ad1956717ed7d12a11a34b30b69eb3945f575f16f783e395bc63d757e0ee19e57e0fe999edc75c8ceb7d3753a1172e31255f383371a23e9914d373278d8437aa45c6cf55c188a7559659742fccf775150f8290e7faa69dbe48b09960db647c9a58c87668b534d0a3f621a7054aa67dd8e5ec09fef842ad1629d9f85462c9461afe770926bcf939769f939e8c35f61e4300a4db286463bf34eef5f4282d6b0b18be8d8eebae6863f3eaf0b0adaf8468d7264e44af2652216291ede5689f9a153e886486cdc5538ecbe575926c362f2c0c2ac473c144f9e7ccb620cee886df07fbfd812ec65f19f8d3bca0de9c999733973b3f42b354cac27cdbcb39c0988f9aedbfa3f0f772306c51260810ed526aa061ed7df06d1f9050c8efbc872e71ba80b8266b5b3da89c6a6adce973945cf6a8c25af71247bc356b92b35ba82d9cdfa189d5f1c60fece7246cddad8e3856ecf048487dc5c05cf047296b35b3d81ef889948e0a2fcd1210bb4511a5ce15a0a4d81b635b024fffcb63c3f97ccbf418fb8ab815a565a82e8065d256bd803d5b176e86c46eaaf85066e598f9ea2b7823c09c4dcdc281cc0593f2c45eafd439f067bb03712d5f2d0c4407538d9f15873474e1aac69b8bbcc4bd7fa9cda8cb20e4a551b9598a41215ae22e0d04069bffe26f338b6dcf4e25c8a36a8fc79d8abb424e40ea6e56b8549de4eca8c27151047953a5048d73d13bb718185cfcab1e85ea2ef8757306ff964c4228a004d4aa709eddd54a11cf25ca6c1c6cd70dd7eb5f8721e9f0a7d8f4a71d859b31d421099b49c393ba3044254935aba5cc5222bef7cd7e7068cee47c4165794e9c0d6f6416bb51476460c0e104cf883f482c97c2aa6ec3cfc5ce68eab2a06a1b97024cb62b4ea046e07727c328b0eda1bbda3a7a635c64a4638f93fa3efbdf757fbc9be9029bc132f296abc3b9a44224faf6bb1a3974763d81b0a1fb91472dae04ec26ef7255e449ae6e628532b912d827ad44897ee7fee066ac6062e20cd59f5481751f2751ff2dc3e117b5059e5403fd61fb36083d7ee07d053e0c52a43ecd1884a5c279364bd29f04856e92427b45cfc5b6c3d345725e0887a1b47c9f8d8607bdb88bf0321d70cb33be6fe1c002865150aebdf3b9c5d1e004df937c565a7d3883158615e7eb13decbd2dd89a04390475e1d214eb60b7ac7ee220e3f980d1515247e0910b55a3cfab29142ebeea1e2eb8a420f33ecb1c25fa3d430ba2f4713d4cbf5a71389f96caf42943a15de969cf972bface48699a57596076f8eceb76070174d39decd0bd5866116e1a67e57d61b6b1da5c4eb830c34bce3919cd70d446da248169deacb42a5f14f4f0d12ff7575f1dd91cfcf8304b286aa4272ff32dad6156b21c3ce56a28455c88ec4b5a8003a1baa074d6a38af77c07d29854e1e7bcc6f5ef1948212831f93128b197e0383a3d2cb5ca528b16d35728fc63e41e6c5b8730e4fe3905f7ceb914d5eb239fcdc0fbe6d9d3de141129a466d68ce26224aa997d36b2eb089fd653d9efe3ba91e0fa64d5b3e495fdfab2af2a4ca04dee4c948b75cdc3fb6ec3afdf4901d53b96fae5c5a8f7fd97b4e38e09b11a41594a33b057a150dad19abc9f6541c30fcfdf73b97d10c641eb2b8c08cf2c3cb333ea307094ac5409babb3d42bba60cb4941cc63481122177421e07a99d4710f471569";
        let sk: [u8; LEN_PRIVATE_KEY] = hex::decode(_sk).unwrap().try_into().unwrap();
        let _ = sign::sk_decode(&sk).unwrap();
        let msg: [u8; 16] = hex::decode("6dbbc4375136df3b07f7c70e639e223e").unwrap().try_into().unwrap();
        let _sig = sign::sign(&sk, &msg).unwrap();
        let kat_alpha: [u8; LAMBDA/4] = hex::decode("a8bb75e2f0902b5a2330d586ac172f8a09f66d9b83867b705dbb40c8d032df0e").unwrap().try_into().unwrap();
        assert_eq!(params::SIG_LEN, 2420);
        assert_eq!(_sig.len(), 2420);
        assert_eq!(_sig[0..32], kat_alpha);
    }
    #[cfg(feature="ML_DSA_44")]
    #[test]
    fn kat3_mlen48() {
        let _sk = "acdab29943a7515582042304cdda0812c96ae611ea6184ae62f27bde18a95bc997a5c7404c15b21f3634826faa42090bcda7d558513c52df9ca5b056dac9110bbeecd6c7750f7c56ff047860c051754c599b192332ed119250c08f5b62539e6997fc3a70022a78d06c58bbf6158162e1d430060f44168fb65fa859e05d5c3e650ca06c10959108216011b06cd4b81021422c0b83882049005c10224a9288030286ca9681813828cb3869d9342e40280a48108991484d09846521c621d0886424a92142b86d90062c4c402e41340a82262d8334865ca66c4ac8291ab7301b066a94b868522469193831803248cb104d58c8919ba48d018245882444048849e190452124685cc0304a246a4118889ca00d84c6880930204146248a46125c2490211389da345251108ddc325260066d12970c49c6114c828d132184cc080c20a84cd9b200dc06715c1886c3a86d11a348048785e334210c234218b36c90404813458022a2602006714ab64dc8982ce1a04c2482895ba44c089750938011891486593045db944880080e22242003c41003c44c08056d19180990062c0440255a226924a84162044423220d4330520817310cc0650bc33188324801b78d5c389189088842120d5b9441e48230db4411412226190141d0980c2345925904441828465c9405e4360258482893402889886ca0160960281249a40093362450a2009912426428851bc10d9c808110a2511b090959388d92a08009b14d1b046964226a94a829ca1271d4a84ca1467088c80cc820820b390521a26c1244045ac4851806460315261a1571c9909150842903a12183906d11114d19a46d214466a31052c430315a1671098521813446e2024e9a9811d3008c1cc690d0c88112858041864d1b004200a16442c869230201211364122151484085080492099784590451123260c40040010605d02871c188211a45689096255c002ae3402ad216401206114344412208202298909c14700c41654b980498443049b4812239211ba969d9a80c8194000ca7250c1045230482ca08648814880413060ab704c0262d830206a0441089b6615c9268002051cc064e40400c24300d99424c4826915310655ac00d04b02002c4499b128c433869d1a68450b89118b00518046c01c30021464e9840654a2201c3464ec02252e132420c030203999123c691831648dca04cdac609c00006202848e3b88920014808c49112468d242900dac24c5b164a19a08104c8014c1606ce5e6d0ac66ae4f6d93ea0bfb456b2618facc7c7d3469b23b4febec280820565e3aa6c4a4afccdcf4e6504cb8a06e857b102f37ada01b73e9f1f5c475c476e4e61b344a1aa4e73fe6af88d47ab1450661ab696521b3c4810dbb66abfae08065235dcdd18eb15ff3ba3c5bf78a4c2117143624024b3ec531027f4dd97f68de7f6ae0a47c31a92b13ac9e52e4d128370303cee03a429d58325d6ddcdd65b0ae9c14943ef5fb084d29fc660bdf7ca7babd6b379e26a60e180277f0f52a74db699306c5902e894dd13b1132889294d69bb014a944c5f9b8b090162d1f4313fbdced92d1406ad0faff6f171341b9cbb3af017f5235f6d9f5d2964f9cd3ab395ab1791bdca23c7f7d0993b728ef1012bb65b76a76db08e98da4cbea31b98c2b3b8dd9b4922c16786d852c821c0e5063138d445b2dec591ef64ad859a07002b0bf3c8d6806e752c741f0e5bb9e659c09aaea09b370d02a9b05d9e91b7388a0fc663e78123257c98abbb7924da92b3551afbe262783ab2e2ccd5ed4adb08433d18b3a8fa4d5b034fffbc46efe776997934080b4b0bacb5583c14ae80700b9b0a51e458955d59385819874671bc8378520058364c3c1a7260b73a3b65883d554b24e80ef5e99446b3bec9e1fde2b0017892703299e33917b71969c4e53d2319546c3fda105b670cf05e48ffca616f342a3bde3cf028f3ca57c2d28a698b7b7c0d13b091130c4b6d2ca30504fc8e1c43fec9bab0b499fcd714b097d0061f73b4222e7dbfa063d65bf67cef0ef2ebda2435485e986027d1560944d9eaa494c3952cd70661b2964fbc77b33eb874d7bcfa9f10a18e05e46c2a5d2b9a473ffdd1b74970c19f716a99697a6a1e8539ff16483c43c48bf5fa18bfed288312a1372722b7709e3e03904aab33e5159a9e9d387e78f934b2488088d6740540460cbd28f5370a9e4e378408538a9a379416409f6fbf44ed4690917c97d36e24996bca1f0b2199979facd41258f7b208232dddc85aa0b97ab3752fc1e5e3db82e97eea1407dbf298770be8510c71342a8de2c044fc695b8b8f0017ef0f14eec7f47d3f4d273c0316e31aa7d8730b9f8140b6123aa9d6926669916c1513559f54eebfcdd7b991c4bcfc5a21ccb67f45f63366f45c3518999b8781ff404f946af5372ae07e389be032ed465eaa85fc51f20fde03ec262d245e8c005d1806acf92ce05102b62a2995273cdf0ff8c9d1e9d870d9531a00001038accaf799fa529f5440196385a749dbf6be5dc6764572594efaaacc5eca898850e362a3b88fa1e3167fb2682b273715f7a30c9dfd25f681673c34a4e037eb3f2d76d9116563402d14afb80de6bf300f112afaa4b90abe89e4f2fa18b0b4304c33a86f598130fe77c7061800561cf37aa46dc10a882accf02b2ff8e864e3a014d1bda577d4784b5e139ed51c5e8c47d24bb2f29a07a672cef46d8ae9daeadb61c1614555db090ee106ff89a5d113d664346227cacaed8f7abf561b8082d0a60020f8955f06402a4ac7cf6ca33f1341e9608152044f566bc9085ef057f591717784fd471f6f9690aa123883132e1efab1ab8ca9c375bd4a6fcda9a8ad9c88a32ea807cd422ded4700ca3bea7adf9d41dc36165ae92943e7e7fc1a5f6de294c82adb596c8d5d2ede0a33ab27a4732a18d35ba1dc849cc03620fbc1fd9ebdbbd65998099df580ad634323e78a1693de4e7cdcdf7ab7c0366bad48903611f87b92e309533b522a68d1a8d6cbb450840a45ff9f506cd3025a2ed52cbad4b4badf25cf9f4d5c39266ceef83d03eaee0a172edee835719f6ec5064d453f5dc07cb407538424a37c990d5ec3d571c72182d683c76b50e62b7a6e7262af3d1cbcfaa3e33b9bb7e0fa53b9b7c5e18a8ad7a79a4199b14c97576dbf517e44192a728a07d939d337f673d2392a5d6fa45a323c3e6f89f18936c2a0cc4cbbb42624cb83ff85fb8a8895fcc10ff54d38598eda655898e9cade8acc2d6e928d4a0fbe8e04a575941662db5a5bb47bb4c0dfdd7fae809d3365ef03337cdcf42ca0a18c4b524207d4f9354aecfe473c2e8d61da240d9daaba0dbd2d0ee6922dddaf0330506b2cfc4e22105d66574c445adb4d37147a9b760d85a98c81339b88670be5c1fa5142c9f9edadec085f012ed6ddd722ecee74f71e451007c472506a65a0a9cf5d7cfd08240d6450480cdbd9d36e8ea2fd12800f4270478247eae1e4aebe0901d78d76eb85572b5267080e496e02651d4191fca98fe986707edba4a22087a2ea0741e6e897439eee5d424771037d47160949d54d39dda7eb7e4d881dff9481a7d6eedb26c1c70881c38cd59c63b34470a2944005f2565a9d";
        let sk: [u8; LEN_PRIVATE_KEY] = hex::decode(_sk).unwrap().try_into().unwrap();
        let _ = sign::sk_decode(&sk).unwrap();
        let msg: [u8; 48] = hex::decode("63470357110828f25b23edc80ed280ecd398a9f53251c3332754de2af0b15e901a43ac1d7f898991f0e86b404a1e2ab2").unwrap().try_into().unwrap();
        let _sig = sign::sign(&sk, &msg).unwrap();
        let kat_alpha: [u8; LAMBDA/4] = hex::decode("58de2cf02b3ed4757ffe755bafdbd1cc8e5393e428d0470b47988df9ed6a238d").unwrap().try_into().unwrap();
        assert_eq!(params::SIG_LEN, 2420);
        assert_eq!(_sig.len(), 2420);
        assert_eq!(_sig[0..32], kat_alpha);
    }

    #[cfg(feature="ML_DSA_65")]
    #[test]
    fn a_simple_test() {
        let _sk = "e50d03fff3b3a70961abbb92a390008dec1283f603f50cdbaaa3d00bd659bc76643281601bd5e9956a15ae1663eb366fb1b481f26d2982bf8fd2cf4c49eca2d67db6f62eacd933991b080cbc26345104bd5a91eaf490f9b7f6e216d27c03bb96bc1ce4f78fc4b2ae9540204e602ca1f48bb873fc70e8312e0921cf7271a8d2e2440683544283218711410338662084883852058758050478363646385837484571147173524816072321527081603641257306347382532536267481365878436575850284705434352126722771411081344185045176312643806017421401516255167453447721835883781087578604872553611212736261514015101073752727787468163288166450868727761176341875013003637664874103774511508674623775622248173527424385054815204883684667036347130030244080644413277225430804803247812582634785513546302814831776037011840708780862245332554047621376273823000436130521836676011430307416075354730025663312330031863220718608383357330664152344216213275865558010057776218622864377177640000850763563314305236308044537836661627135082480586046432640678072641746460350487412784242874180008730672857437536012503420226644035310867460482842531843086455014753814346426521313753750336367745566822885127155307342800745745130137536150054713764380545280054118610437756470171226472728524846208268403310468084554868134666083088336832782325307306082746018200074255808823574142888834632438268305151216460413835302361875536325304168531030704021236175533030425442748551342124416514404361546075356865406135254840325548063816823620148751373713217445358862801817531242521553305027502687630664476223130131530815720853317660377518157078143807336761182185506061715315015034354730433677406231037404254568757502118673105300024542430143071317443404680172400338637802123707781454563260627332604713882724622076780246704051304604541475456108643812306828417756528135057720017841160642466588523048685218586288160637066241257257131705514455081640474244535812685347648231548862501467252562860704255646864728410158574131166162317305611637088625821854076548102256140134778560017445532720326500758764358876236481501487662340123175011053780706748855485027208874318581144368287640784773500370255657032857327052066231174612156243528021482638882018504848636006464772258547731777822738421007047615512257677737087065446046544310266745624306018414756210448786668625070848823514214244120284348534537485185470555415418256244866002033048082781468643808527847501035846124147841043747710804075600127234875105878010281533773131421124706563088868162824775373605858175807088517533755885416861622168804137550770114066084054436648360063126106444515437085545542111085430816002733552438303647576845002357356380387242742587107570066383516467402462816142512801630402278787868810373218144404621374620017735314413387310647477274327563783253725532803286760134120620854033682472560462444033547464673557868585371124842731188085002442411827560127053252877622544836565826271621824070274822284126741821257751820747713843371608366664432668611172776375330558522056032460307465863314363241275836725866426532565411845273355116166302020473535368001522575542644063487005701006816250558658586465838272158583246575864007033266542123163863827802466027311840023825036cdaa4a08813ee62bcc774162ecbf8a7e7ac3c0e4153375d1e7617f1e583b86f2170ee570cc54040bfc01f52c18f186783e0041ed02b1bbd9ea256c48a1980d1128c43c37b5780e5d0d67aadc085bc0301afc07221b4d274f08f525d05ca7357e6a87d631a139de608cfdf043f987bfbff35c25886fdc2a0b767b004d7be63068534527deabade7e8a513c752efc2d0d19a6cd689770c7dd0946901d13ee5d54854e49cc12f4fb1e73a6187d2daf12ee4d9b14992d6446686ccb48a7a2ba390300aaedfeff5cf3766271c5298223564a30449267a19adaf5fbf7e10f38ed4df7d906f222e1d0717f351a1923733e4d0edc09d0faed480a46001a43bd42700a38d1751153c1e2dbf6151772471b1512962d9bf6a4e498667e92e56ad8adde32c0917aa97aeb8100e20fa08036085fd3956e1850828a966bd316479c4bf12a40bf2c621ff618e84d8cbc833ce0fb5ca8e124b782861c236f9bdde0824883a370f1631b5a8a3bc85b9db1bd35345a6d35cc8e4983eba7131d38a07ca059194f4fef86b5c6cf91a2d75f828a1cff0525e6d83785bc13908fa5bf2606ffd9e6d80e89fd32a308afeb64226fe24e573a370b87f83fdafde1d2ca594af3a9f43a0f4a5ab95eaa850959caae0bf2d3fb0d7379900a7af1b52a8da5cd86fb5213ad7f90f24f2d129e8314bee39938e6aa9b07602a058fdada0ae60aa2a0069d70b585e8daf39986b13e5dd1e1f5e81ab8ce60243bba8d7c3a39320a8cc3a0c1f02f9c1b049fbbb68b1d591262120346b87cef150a8b8e37829ecc29ca4e49ad954b556e4c4da8328a5eb032c4624ccd80d4eb1546384fa68269ffbb0b4be80d1e2c010cfbce44db42a7bbea0c6da527f99217946e0401b79a9f3b10ffd73049a0d31764058bfd09c7f2c075b06443ab823c2ab752a6ea33d5238173512172bf70126b466cb30ca55bc04e002114e8fe174ea8e0576c752f2d207418c7da1991f0f2827b61bae892cf07966c5a3bee880dc0df2a6171abf52ccd01a8466c51cff08b1982bb95d0e8b61e14711efc03a1e3d194355a1060b220953111fc8fe1d3101789feff0a8153eb576baa329f770e51e2910f3289180a37c2116b156d77fdd79eeda574b084b40863866b3fb24d2311a4d2db0cf86dd7f7ea83870d3a1af19cd02159e7fc76294c8a4cd0a9490db191a52b161475b95b3fbf2def83a12c1bb2618f035e3b9de9e8014fada523eced82b60d02fd808c3f4d4d9f8045ad9daa821bc658733524522862e3764691fda010be8ea3cdc3655928e51f9bef681c0a87eef2ffb75b4410e0966923d53d74f2e15a0b62eeb0ae9c458455f73dc18d15cb988adc7e368f8c0ad1e5fadff40593df719faa367c0ed64d9bbf8931a71108f7039c034ec06734b94bc5a4341a5ab2bd54d4e9c8794169a99a7732f70e871d50e443cf212f1535e867834ed199c71ba29e567024a5ed9e73958c0feb9bb428299c1f14b39e66c0d0cb9c164ec480532d09b1d0bd02d5967a8b17d5587b6321b4ec52b7084128bd9edf3469084ed7dcaa44f7d4ac6db0cc3f19694cad755f20e3285c2b8a10d88f6273bcdd17e443c41bc9793ac22d4d9d7a135a2b0dff629a4b49db3d95d3b5dc93f5a590943b09c6dbd4323016cfc7e1ac8cf53b28ce39588859a201697aaa3fed3020cd23c3640dc9f8abb2cb460afd992913ad4b49260c271ac48398ca5e7264bd20012971891436bf5ef4efd7c2a7422c5c87570d37a68e74947e09b174f0416af960073debd509d7283dc1deb4fb90b06c2310a49a970cca9ed34e6cdac6c53a8947b82bf48e37864c3996aa6bf8c7b909c0830f3786716cb8e2fe44da7fdbe8cecbab5ac51ec339f5e95a687c0e146d10889b141a218c4b6f97aa975512adccc9eb2f52639681b679b9231a11929605dba9cd03fe512f0ec9fddbfd21e8079fd45be86d524a470d7d43cdef024a2d72a501e22f019c5c4e519258c6b5f904f3ffcbdf4654b126c939049a8d46cd5dc0092ccd68fa101b21e32ecc14941ea3e2312c05b2e125670e95c0c5741ca888fbef5474ec55c4370d0c7548b486c258a843635ed744c97d40ac40fd76e0b7cb331f6132f3b92944e37954bd6a83d7fa0e9fc38b4d61ffaa609164c328bb652a1d2ba6adc394e0c8bcfb231266dceb5ffb7107d37348d6be8325f3e48c6aa88b52d96c0d9e0ed94911783a68fe63e5c9124024f6371dde6c3c77d7fad4bea0612c128bffbb999cde10984d03998219f428bacc234ef1e5b48cc4b3c2e7d3d0554d2849f4aadd740d811375ef21438d1fa99825287cc51be51a21577a9cdad49cda2707dd0836b17a50b458724274a89a71946d5c51e5fbb10b4372da5a3b7db6c284f97b965a76d311e245899288f5eb860ae2e7c4188ab38e929a3313a6ae9c80c9485cbe715ee64a220ca4abd7b4fa5c78ca90d943d4e35c9d4e488c0a29b2e24347e6a50a631c5fa9b9941c421e8f746f969106c6e7cc939e81bb0b2d644c0375c96144130e9188d5963af8b9b274205b5c6ce0341b511b2ccdd003ed4f0000cbc21898948e09b7b3a06b19c2148de65e0254aa3abc3ae142f53c8705506a3e2cb0a8b9d81043b0feeae9eff4c113214e502bc24fd471ebea263127e3bfb0d78feae0eae3cd2dd3978993c8d6497184c9529271111bdf313ca9581e06fa020f1efcfb60a84aac1305d21d5f08fae1d34ecb72f1faec21b0bfee08281968a180e49f23d35c892efd94389f4ad80520a2524160abbdf9d80ba44b00fa54a326691b2dabae9dbdecfe3248ea3cc19926d9773a221c07aa2a76d97a966226458bfc79cb13aaeb7cbbb03f0a47c1d0f00b6c50a333e390bae54b86fab92c94a7a31654cc3935e242bdc4d7714aa04f32106479b21da303ae0dc2bd292993c555c9eeeef88106c121bb51ccc71a0f9c00e7963dbee6f528507b24fe710e2761b352242260df8e2eabf200da5c22721d8790e0e4fc877deb450e1d5177d3f99c7cef6481ff609eedba5a8ff6a9162ebac4a3584b392ba50094fdbc85656fd97688bb1466490c8678755040d57dcbc9043e43d97831cd9a0762afb5a9c6944d2ed0c5f33179de99e3cadfb53b225cef0b1dd3ce103db6794f31799d42b26d1123d2c0fdf4686dcc68d6d235d19aac4450291f9f46e2ef7330a41543cb0b59ed25ce29938cb3d4f319fff64566e79ff3ac77e27d9902e0b7d943b179c420066c6b11f8abfaf3c2f3e6fa58c80653c85aa2dc72fc02d6d9fc13ca5f41a56c2c4ec0a0ce6bad7eee44a8eb89f7b4fda9202f7e03cf3b93da515edbae052ea2454c7d396d96d0888e2a19e89cdde74dbc3c15953f5f3734521358e3b2770038346d0e9514ef50c17e774bb60cedd6f05f0f1bdfb1a14ca18cedf44cd7e144a632a18b273691a14c3db28fff5156e5660f93e93bf2829ce62959d75b60df0c84e1de62015e2c3673d6428784b08e1a2c391b64e82e5c40767e210e56ba251d6094a05340bd0d8a379ec86e30d15bb2b9a1";
        let sk: [u8; LEN_PRIVATE_KEY] = hex::decode(_sk).unwrap().try_into().unwrap();
        let _ = sign::sk_decode(&sk).unwrap();
        let msg: [u8; 16] = hex::decode("6dbbc4375136df3b07f7c70e639e223e").unwrap().try_into().unwrap();
        let _sig = sign::sign(&sk, &msg).unwrap();
        let _exp_sig = "a0c1af32f9ba4e4beea3016b96d1c780e8b5e020bb07c24478dbdd0ec875666b5a70b7c15717fc4baee3a8fe50427677a36adfeb9341e934c74383554282bdfeb4139595896a65c83b614686aa8b259f2b22f30c50fca1fb664c98ed6d552b4e97d757d7ae8f3968efaae431de1e1cfe5f9b5f34ca122b6793126fbc344c3f78ed515d2d6c59f931e8d613598054590709d96b24561d52c7a7dffc7ca0610fbbaf4391f385e28db9848d5489529e4bcedf73d7d8484f7c3e02f2a7faa949e3116f96fd2b31e84cd17d9a08311bc1cf4aeeb5062cef5417dfee773a6287594553544af787f3f850fab41668fc5fab9003c8776c22d15a09d7472ed691ec891c53bf41c280a0f719e21f9a7bbf890ba25137603a007cff5af2b3df0a72c81fb811658802f08fb51209b5a5f777859e0ad53e0a4cbc1657fbf6179472b420cfb99fd87993b31875da238a3d4e148f676b6a0c8652b4fd8c6201731ae12d44b96a1c14c1e8548a3f3fe406870a4cf6fa60e38a4f8d1bc9958560b924578df8d30a9134b49787df4db4441fa9ef587740527ae176a669ff1956f256e3cfd265f2b17e8b13b6eef20f54c7659f235ddb61ab85db1b61a8842f09473c3ea2d80b91bbf6ecc4ad955f99538f023992013bc118e597e2a85c4e61b47279648e45af1b9f647d31a9d2e4e703478db1de01a12483f36ceaa7b3e44c622154512283e81bd9e7dd10e15abf85d81e81dbd7aa8553752955085b6f0ea9b126b4e6d69f662f8be3b91eda4265caec4e6e10c1c1fc82a917b27c37d89f274b9dd9b209500108ff5d9ee6439da98dd938ca473a99aa557adfe3cf921e3ab24d6ff89e870cf919ea52828b0ce301d28b20e09952557ae5bd82333e02d1a7c0c8a0d22750dd0ca3259d8a6123c64540eec0e704af77b20ffd58410410b0ed358e37bbdcace3f6ea21a95664e27689a12b925504189c045f70ddfddc1d4fd34d45371a6c486efec9c9054a830dba9fd259c64170abca9a06eedbd88842e092ed2fdeee6050d6d6ec6c610cf61046949f8831a36a80fd9f149cfe9a57e2d68feff6f217339afe046e2d054a2407155429e20d0afcb9ce3b2ba9f779d157c43f44f8b9d94683795e225357889c6937dd63e583c2bf24833009b9cd6ffd13a3a3e7b768bc7934a3afef290ff6166beafb1f91d706be4a692a089a0db442f4d4890e30acf92073c7786ab1fa89b13d8768da6ba12e5a0b71e9302908a24725de67cc328dea7864e74a9f60f4c5634a6de96d1ecf6f2ea4a4b3b7a67c93419d6802628f67b970acbef7ae9f9bd1f4f5151574e5d0ddab6f6792facdaef4aca35790c6a02589113f940bb86311246f11a5641402db55d3078fff9a28b8fa6fe45019fc80a4cbf7c6936c6e61f423eaff05e3add8f1fe7f96e3b9da15b0f17b623d0cfe97a9946fb5715004822a5fa1e30283a0deb87d2bcb712f84ff2434e4b1525f3388d9bd215dca9a3e5fc77498f0cfa5340c6163a5aff39c018379e212c1059d2551098eeb8fcd41e2ae8ac81bab060be5be32b8408ca34ce7f58cc9134b0fc85fbf66f64c0d885de0ee5bc40d20e81fdfacc79ab8efbc69f27880d84533debb3f4177bd96e04edef1e30ebde4ea6bc54c33390729a670d896a57c1ab2e554fa2ec4cba545431a5821e936288534e011d90dd5393029f99232eaa279ea436155d2d9dc430140ef13e611c15e341e689f78b66dbac01c6cb0df70aabc651d9cb439b43a0b47b0d7d542539c75962c0d6db69609713a2e83236423781711a3ee1f370d28724ccdb25ba5adc7f0ea603543713b0ba7bc9748d1e2a76b66e47178d869fbfc4963c4db433166499cc65658e369c597b6a93f7eb59401d6bc6069626c3575e3619ca4129b482b43a50e2aab811e2c97a30b409d6ba9ae6644338e341e9b4f42bcf7b38435a06b07b129e93d1f9256e724a2360eb210fd5fad1c2d525510d6e634917e7fa6609dd735251bb6ecb5b6ea0f843b888438f1a775b41759b5693167381b95c37ba0b68fdc2d2de87f1415a763d06e627f87023e1f158bbad60a4bc9e50ea9e935413a6e5636c1dc7999a7fe05981d0f75c73192c7db462c2b287837348348471bf39b9d951b597e171f4ffb7497a1b18afa8d43ce251b2c8b4fe27071a713d7e3f465555bb0ce7c7c930af2c1969f0b0cdcf84180797df2eaed4379249911cd3d7cbfa36ac6bda94a48b709aab3a29ec06e0471b9f9d9bf1871ba850a0347677db776a0d276cd43e92f833352fb67692b126ddd52f236d3331ce3ec861250e71a8f32869c29ceb4abd3e94021ec7502a7b933ea7a42b994f3f25a91b767af33541616d05721661e25f2c77f846dca85151de97a590fe7d8da703165bde544ef173d430825ab26b229a6a88e813e0727ffedaed10b46d7dbfd1ace868f88653f045cfee71402fea6d2aa8c5d4e54641a87986162467efb24627eed4967c48fb3c3de0b1070502bc2f2134d8e537db2c0fb2d7193f41074a52ffab4be9fae59cc31357562158e8cf5619dfb6857eaa2855efff8009ff11441b1d2af94e4a039f29a54374b939d9e5cb69d7e92b89ffab38f1248e0abff1c06492b7a8871c452c8072fa14d4e610f58170579edfe9ef91b6936f595b35447423833849636214814e7b101b326022b7ae820fec812b6bf3940874f34dc106d0c8d4b0591ea00f984c6a3b418e88446dbc8fb9f3eecdd67acaf6aeb3e8aa0e9a0b148942c64dfff3ce0780e78f9c33a51de8890d0cedc664ebbce40931575cc05d70703a58e44c792ffad8448aa7808521c7b0fde979e4942ea043f91296a292206425d53d25c98324f297ba1a69212c0f5b7906fa828754d66b67e2aa5177be937f506039358f274f25ad673e87c8706c55be13014afeca927119bc1611107eb3a0905838f20b35b61acb32fe43ae79aea325a7dc57e3f7bea8809e3c033975b1344cf88ceb5e6eb4c3653e552a6001cc0464440577d537b560afe0b572a534a7b12e460b1e43aeab96cbed3c3ff0cb2efd739f00ea0b8393a527eb3dac3ccc1e554e84960bb8dddf34019c1f1276949d30a8b990414ffef681d4d6bbcff719bb750b8f22a38320053da5a7221464949f7300d34e957a833127c2e0eedf5bd0d77c14d701fada0054c3390d03bbf0d119b3e832893441d222cf01ce22dfcc29bb149d6b3a180202c3f466c5bfa85c73296b26d38700f559282b07aa4f09dc0d13f32b5054d480e1ed2a352d1918b636612198a9bd9026097e4a5253f9157d9879d1b7ae880daf63c9a110ed2943b88c6fc7935bf00077134f8af4ed48a7386f8c62d9e6e3e3e765e49de6abf67e2415e5d3f4cd45507d77e784f062e3622e06cfd577cca25c87456979cbe503a80372b13f6cd43386158a8b23558774d74daa8ae4526d491ec55712909c5486b8ba0b17c2a8794997a88852df0f2b72291ab8b1fa21326367dc9261780d6283e95dfe1a5f6ffbee56bc6b48c743d45e90ffff3e246bc84119f8ab7ea3aadc9794182a83a07e20124c45678a9efc3c4ab509a125c718e365e694f04876f352942eac70d0ddfa52592bda0e6dd9fa7c2f8d7c2c648365705ef38bce7d55ec627d70128c4d7e0a102680d6be0c8eea3d0017f35a2eebc79a7d53a2416719c65dc5422eb3031de54841a12a7a18e86cce232bbc0e032ccd787bb6c100432ef2e147c9a37beeb48182de1c413987471578863a927e27c03ff3625aa3cc209ce81da4a43248cf3ee1b4e71a53df4991fcf81d03283705c30a11d664640c0ef5afed725a899c946535259a20d5a4aa1afa78ee0691787788c21a8ad3db316e127e38c4a734e0769d62ab47ec37968bc7267519cd3796a7916822a47164bd07276afe38be0588e9b82b2a4abcd560f26e56a696f96ba2ffffdd0ab45607961ec9cb1c3be24dc7034a9972b0660b7b6b0d6f38ece55db14f5993cbac885fe4f9032496deefd5faeb924cc66efa510b633629d968180ce15c20eba716993e9ca5cb7d3aeb50ef661ef0bece15f8689ffa1836ddc59cd7287c450b3c30b75da120f160d29e426ae885df245e109458f12cb0ecfa0963418e7439019c7445a385f5e807ff8aec205a14c0cb0d899b47d9390b6ab1b9fb24ed33d90ea24781b1a4c547c49f2213d72718e2c29d17258f600b00341e9b4ed88fa91665aa0ac1bbb651c9760cdb3fe830b8d572418ae6bd06b41624a8aa9abe0ea42b3da9215cc046bd9c2b6de937967a14473f0ce0d8fe2fbc17321e31d5cbda18553657a8d745b15e14729f5719e835941d9ed5ed4f8e9b82e44706a509252d9765720f49557962c9782fc01c3186e73320833a4cc270e2bb811f682e97f89e98782da75b868c4bdef04b8e8c313557670bc3b1fac84dd50b4752f38a8d111666721948317f8f1529a2a4e3b3e26a05c8b3377f9bd56f3671ca7adb287963bd179d0101981cd5e84bd0717fdc7c5b03b2e6f332e28681c2b7a908788a6ebc848cae66cfcf22548ba9801fdd7c749648f2d61f54361e750cc0649d0b57c8e83a76457e482233e3cdd1660a5fb895f445f1ecb66d1035a981b38716bd04577e314343f81f93a455880aed0d2f3366fe1eefa2b2f6787a0bcd587b8fa000000000000000000000000000000000000000000000000030810151c1f6dbbc4375136df3b07f7c70e639e223e";
        let kat_alpha: [u8; LAMBDA/4] = hex::decode("a0c1af32f9ba4e4beea3016b96d1c780e8b5e020bb07c24478dbdd0ec875666b5a70b7c15717fc4baee3a8fe50427677").unwrap().try_into().unwrap();
        assert_eq!(params::SIG_LEN, 3309);
        assert_eq!(_sig.len(), 3309);
        assert_eq!(_sig[0..48], kat_alpha);
    }

    #[cfg(feature="ML_DSA_87")]
    #[test]
    fn a_simple_test() {
        let _sk = "bc89b367d4288f47c71a74679d0fcffbe041de41b5da2f5fc66d8e28c5899494c1f189ce692781ab4fc4fde69fedac0bd4bdc2720d806aa5d0fdd9c100afcdd13f949bc58e8fbe98580f187d82e0a3333cacb9beab3967e4ec827f12145c0a436483fb81ce3123e9bf1bec9571dab72532713fdf13288fc2a6dd9477430d89e6cc404909b80000330020c888121181e2b48c18174c2419201c158159c80808289019a069112105d2c8909ca864db188894b06dcc226c14034c13476cd3828124921124a95090b630d882090346115486488c4650e13869cab0119986712401040b0288ca424d94a04dca028e20274cdac421581822d2b2045022801a4348c8c68809492e61260d11a78900914543168dc3b86dcac26010c39090260291283110310cd880319892701b0412e1066902c100dc00300c154091148d20a06910b5311cb9291b172553200501a770dba451242530582024494692203601dcc6680022844314461a9831c90209a2964c10c91101938d194710e406002196088100859bb24920448021b4481a36711022900a232601a809c13689c0a431e22209c83662e1b608884442da3289d94401149311991646e3b861d3384a1ac07051c884004510c8424a4414520cc164e28824109648094444a0322de29889e41804c322290b078063c42593963182100818308d48c469a23280a4a290224149d902415ba06d84343162309210c905a3c225c4066dca382418266a0b3865cbc28011a860d3046282140122884d19415218a0601a9124e320688b8471db988d8884611943918cc681e20405412088cc024e48c205e3462e0022800a4304db4644c8486401c324a4180a0bb54980b4642118066116244cc00043224d240244c8162e94948cd9160612a5610aa9811a21260a422054380c4c86490925464c062d4b3471c8c6805c98681ca320e0182549160ae02245d8100441c28923855102498ce3b84111209008b65061024de1c4855224281c1046e01825232071d3068ccb184021046a218348640629a228429838051b80880c8989cb12620c024918a94d04284094b44582c800511008d184241a40718208421b202223c380421402e02082e3960d4128010a095180822823c56413a38c62088c58009024020c0a162c8bc66cc022501b249054304e223512e1901064368214010292a88ccca491c0488098086edc040ee0148a14330822004104378280902d48400e0421100ba9800b246204376c201629dbb009c322429b34000ba18d08a670c11811420486c1c48d18c80818c18122a6401c372a22208443128662346a20146402254804b465e3a42924298ec8c808a11420984090c1284c09444cc8422541b068d0c4895b207054a629a3121294a02410245014c02c21306844308019364498306cdb261200430001302ac2206289a028084811182548d9420a8b1891891864a222680012701116625932220c136ad8b201a3248d80a8802148424916241b976d19978c0006201346729b202c20482ec2100464949058a431c04648cb36209024020442895b180c91c64194263200864880c04d4bb030003184144922019309c328892311441c438d41124d10c3281931404b36612106610ba884c3a48dd1321222331040222ae238800c282d60a48954088219446421a2841c89858ac4211c0808610891e028814248888ba4511311501921524240664ac828d8162858b0092108869cb4204ca01149b400d4c241222524044200d0c8100bc4481cc109e2905188826dca842411080de1b84d22b45083260d5a284523321021290c190564a2240959a0905492659840680aa92941b860412860c194410a492218b188cc80205930710bc700644484d1240118252918004a9b466c9c484204467064b800241745d9848c239980d208010b348d9a120223b5215c32009114200392045c364e42328c11430c1ba881c4920ce2b02489284d9908852431851a3021cb38321b9140e088401b152adc0651d3440d0b106c612022142591cca291104344cc100903c369209770e3141061b44c584244d044024cc88909c0510a0708a19651e324900c416182240d928400892246e12270d1a84503a071411620d9046a9ab061131550814232d12268e4046510825003460d998609dc1491d194110a278914b08c11150584228941422e239949a2460a24200894204736cbe5c30f5aa9b510a43509281c52479cd509fcfea56a095e6749f0e37b454931f4bb8e77d4ff8ad9f6c47a28b0be044b92f450d5bef06b26005a5696d27549b51a78baa3ffd96b8348be7d81ceab91fce852bee484a129ec1cfd14cce08cb274bdecceb5446ce2eb83d9d1a25bd614d1afa4a2b0432f47a2a61c1b516e69c87b82dc41b96f04e5f0798869e1a48d3aee337cc227e98689881216c7de6c28a618fa2aabd4b2e38c522e47f828529e2a30043731d64b6dfb9ab2415aa52eaca29fab274955bdf3f3f876c6e7c9c141d96f6e3ebebccd12980b4c7c7cab75f4f1fd0945c786b05aca32eff6633ef02ed30533a2409bacc3aef0ca3f72e667cbfdab22dad0841b5fe8b7acdbcdf92c494fa9f1f6d722fdc30afaf817fbb16595b3aff5ad824e885ce900c4de9cf801e22b0fd1552c1c57b43b3f6f5b78dde5c5596f758c3e4f007dfdbbaad93e8a6fc41e6ad954c6573cd079f0fec1ae6c5d05a646aa490e86e43c90e34c14daa12e595853028fa5600e9d7aa9fe0bcca0f2117af166ff97de72bd0d1bca96079c2522913ccec554de1088f310660513059027642335c7baf3d3a02a2de6c100748023d49ada39111526290aedc0ddc00e846d7c4fd3d62d8cf922fc2afe745ff6a102ce9dc7878aac27c3f91100bac46a3a6a8b0a52f8b11a2f13fe403f96740114a7d431d485bd322b03e91d2aabeeb5972132ed7d9c39cfa3cd916afe9b7cf775ab796584ce805fe6384555fd0b59795c29c368673a95883704e417d5dc4ec243205e20843f137384f3a9bb51a62a165e1f69e0880a7af86d81b8dfc8a17dca9f27035a4a5b5b1656d5d44c3f62c1425a6f52e6613ed37c0b595056c1479d35ae34059dab77306754d9cdbf090d1ed5f8a3226e4df963f0e9514743b1c345eae4436a603890a3cf1c612419e693434a6bf6ea816922a1fd48d8e34a51fb075b2ad07f675cb2dfcb86d6813aba0d7278750e5940cee4d7a39b0215d93510c48a7edced3e253cf6d61977d4f5effaf9770e6c303304ded7a08f6a00ab9ee50994f9956ad03734c90b06ccc0c035b70bfafdb37ca6346d891fa24bc3bacb50392d5f4f36bb8515a0f5e3ef0bb8f4b1bb33ee3cebbb0181cf66518506280649f698af4f28589e47e94206e1d03839723db12d4629851aa663790c7e7dca293025a6a7060c22fd671c9f5bbea3d0b0554bc6cd742fb785db9fbb02e8a0003cf24a951aaeedd2ebc09003555b010dd238c167efabaa3c354cf713c5b911824b827d0f187220a2c2b71c72dcb9f2ea57d99ced1e381e52baa87cf71a663ec4cfe26f7323bc9f4116b2c60e2ad42ce116dedef1ffb293d1effcc25d3c6b3e55198495c833343dceeb01fdd8944d4bdcfd69f2b80cdbd6952718b56a22f874a311bb80006be154a20db521110551abe22378d2b6c440f6d7b93ebb162403817a669cf971653622c7eb0a608e2edbf55acee5b7556cbe0af2c0a20384b65d9a5e50861793433c5a3c49506cb5fd232ba6a827c88f571e03fcab1a7a397b78627515b7d212a1eeadf7894b4c148f2186c59770aa188699be2700a77005ad109baeefe75a7ece8d448b6579ed8ea40727ca0139c58833cf62168d306b9515c0b673e794ddeae40840c13a4eb814160a2eb57052505f515f4251273cfe31d31369617d7d1683ccf53c48799b10da50bc946b232c1827a695a17e36432779cda3e79f4b9b2fccbddf728922f558f34c0caca8b6075932aba79dce54e416ced4ac5719de23fda6ec2a32b6dd8ac9c2171e131eb58af1a10989a9b2b7ddcca72062c23dcf2c37fd25c716c369ea13167fd25977802adf32fcd5581f9d7f7f75df99734e7ca3b90673ee98c310947253c010b575ac3d3f0548c275f959bee3188c606fff6c3114a43dbe52ffc1a0a942638aaf90f0c3c18817b9c52cd07af450b7d17a18aaa2b291835266344feb9f741609734f6c7733a609cbab72bc58bc65b0506dc90ec0860dba44bae3ec2bcd3fa765990d9072602b91a0988e76d3c9e275e99de3f0da5b1294a4bbc11b7e3b59bacb720b22c0622486399ec6366773a627d1d4775b7c0fe3017ab7994bf5e0ac3e1c189a933d79bd6e9291b38a188965279b6786d1c07c429c841c7fa1174a5fc07eab2b636327020bfe71ae91d54450194534782c39bd005ac968274bc99458d0e1e703d0ff5f361bfa812bbc1bfc155adcf6e5a19165f64f8b73ea28e9031bc4e9a09595c344a2c50dfc60f381e628f3d468b26090b95073c63287bedf0d3ce7decf21cb1ab88ec084be48aec1f37e26fa32550d8ac646c9525c3aae987d57ec2d53bc29278005e129548267e9cb77e9854394f8c139d70c64de3731b2d055ebfc41808afef97d77adb6d2b7a64ec3b52b832cb31ada5dc1907ef0bf7e1d9c211d28ca34a06ab5c276aa0948bd48d7166ea415dc36ac8b994f77d78ac8fcd1737018bb9b2bcdf3f26eadcb29608bef478574acf3ade2bd30c284ba6c66ec443c18fa0c54c82c231ff7f24f3104e583ca820f484d5c4f41559bab3a4b78566633a8c82a7f482b3b4abc2d840cf51363e91f0f0f8553a116a3ae04cd2558714cded80fa77cc79c8fbb9030bde40a1c7b78f4efa3e614645cd1e3638fd29ee0dd34491f86c16cf5f51f45cddd4df44da50bfa3eb9da3430da72165a4f38a588d728c5d3479a6bf7455533f33dc4b7c9d6af0836325d382ec8574a5610592d0f1dcd9e83692bd24b798fea988a8576746efddaa961d8adba4bb29c6902fc5ca075fad6bedd702cdfe58995f1796664963c5b2e608e331b316e8970de0b46f28738f5e1d9270e7f2e967ff83e912f6ab6d37e6e468679ae16832c5cda417adaf626b7c20c905b55dfa318ae5c4f0c8cdaa6c074968cf793616aee12886ad00f5167799c4fa433548c6f2935389916cf771d420650bf5eeacc2ddf38acec601994120d827676e63ac0404dfe97663fa05d36a52fb35263b32c4dca73e2487c7d10749011998fc580f3237abdc7ab3c54827a2e101affcef728d394a02afac7a63f521a100f550201e48e79a7fd215952e7a7890f8a40d2c8a4b0eb402c00fed53f5046648457f48dae3f1ef2d36392836db09b44689318c5116977e1cb0863768952c734954c637897f2f8e0624ff474deae7ba7fedd86bb81766133b47d578f142cbfa540e2e03e108d0566fc39d74fd50f88178a173a7b2844646b2b9dc4a7bb803dcbc94e382b63bb69fc2e1a2cbb6e93ff366901259c748e279d47896648e8bf1b2a7f0dad0e40e2151d60198bf9e2d2d5420aefcf6ec201c8eb8e011f7629bbd16b478fb3650584380f543594cbd964b59bb7f643c4093feeded4de7fd62a92c840bdd22ed5d35fc87675bb41ab3fc67445b08d2991442bf8c728261a62f3d477527a2534fd014cd458f91d45e8f34cf5e9ec7a5a3cadb9490b870af8eeece7b9c3cb4f6b9f9b39978fda600bab1044d5ce40ef0b949163641cfb45cab9c10823c1cc61311179354ded0db9abf882a5dc6a22d4ddaa9b726ab4638b9f2e2b92914e99223d8de095481baf3d2945040220b2b920ee24fbe18077ee8cfc78509e9a0beb1ac986fb3aaed82bcd8af58ceb93af75c5eb41d387639cf65f455c6961a9a4d8e943b139109abd26b99500a680d33c48a1e62743348a72b96206f5a8fdb11754beac0c213b26f780ea1a14f77ea149ae5c3dccc1b5c634fde51e6df6baabeb3a269e5009c5ba3092e06be795051739673a8826cc8403ebe74aab265b0e51c45797c953424a92419a7c23284225d1d75b7eb51bf460368be414e4d516c56c6f24fa4dbd98dfa95ff428ebb53c179496a0a0dd9777db4d388f3bdf7c988fe8400f17ead8d6246a792e13d0725639882d8b264518d05110a8c9fd312486e7404aaeed858f1e2d22f580e08e40d2f0295b4228107e39c6da2e48ffa1843eea4e5b0cf7ad369453550f8dc3e4f72ff4867950a28a71387786d7e583f8b66c23b5321fefdf14709ca1a17de2c434aff7b2e9c2a400c04f53b0a856ca4a1c6029730b9e9ffeaa103e63bf7820e57acfe32a5b3688d5e82bfe23d3115fe62a5781f14cac2f10a5fc46778636a55d3da95330caa22e8da1715dde84a927aa7847e2ea4c914a2f5a1da1c961767bbcef4c750d0690c5d0244e26adbb54b82c5127c88f614149057a6ffb625e4b4d76512d0dfcf7760f1827824a07c8af85a4245ba6473dd9c2213b1f0acb1ec7542451bf8045578fb946e27f7791ab96e9529d966e0a8c32b70e9831bfa30d430e1613a97a1ddef91c580056d90d6ede641e110978804a4717c2c8778dc6ddaa9d6a1068536c5c7d2bd4ae519c12bc14c50e332e01159b55ec616f33cdf2c78a151fbae15fae1ee1e32ca6ec5df93cf9cc879ce5704ec5f6afb1aae40976c42f1f76a853c7e4e6959cbea76ffad4f9d4c8753719e38a4cb44eea3cf4f094b9f499d3b61435b3952a19d678d5bdda57704cd488e3120b613787a186f36bb9062dae9b70e5dd9787926229eee3ddb1ef2ae6423376ea8ca07bd5b9c5ab929e45b93f0235514427417189dfd043fa4afe9f9cfdacfe80837c11239317adb64d3cc54b1e8118a646f6aa9d6a637bae50b0be3dcf4eb984471c9e8101a12788eea4a3e7c91f64119a91398f4f454045e2d843c7a610c895ac48da583593462c";
        let sk: [u8; LEN_PRIVATE_KEY] = hex::decode(_sk).unwrap().try_into().unwrap();
        let _ = sign::sk_decode(&sk).unwrap();
        let msg: [u8; 16] = hex::decode("6dbbc4375136df3b07f7c70e639e223e").unwrap().try_into().unwrap();
        let _sig = sign::sign(&sk, &msg).unwrap();
        let _exp_sig = "47dc5764266841c1af3073fcead6a13d372979e6cca0b2952b349915f54ef663126451c742fc3ffa23af8ac195431b8836cbf74db223574185cf63ba99152f027587977395d81554380ec2bb7b20dceca31a3f97962b0e0917d083409f2b765c3f542799394072b396f5fc64272b5d1b1673d98f7205e2d53eeb5234bb79f438edded4d8710976c4cb75cd3ce1d495190a268bf62b1e7c1340641e0fb18ceefb7b6075e73a00568347c25dfe509daf954151fbd084e00da54b29264b0fe48aabf65ffeb862b40dcf123c0b23c307cacb0be3e9cdb97639efff33b1af0b8ff5db110130a6bb0b336c85a5c54a428d03538c4296c70464c9c65a60b75ad090160fca7f496080cdb3bb099564bd1b15902b51705c5f7ecd906e563b9363c9a39df1bbc8b15a82ebbd4b0ac3ba461699c43b6400b27f0ab1bb8cce06654159091d91e7e4562f31a0162a08116122b1233fbfda3639566bd6a4897188ccdde67b41874d20a532826903666af272c0a7e2edbc05fea2d48e556a2e821c418a6462e2df850083415097cc6d000c8fdae268e1709869439d32f3a8f1bcd9c2378ecbeb7cc291bff5a50281bf7c394116329a6b64afdc04e59ae65af5cd8c5b57ac0cf3d1a6957d4008f9ae58c7ef5d74a1f877ab835a937f378b99ccdf30f9d91febf2b38ac9ed9cb54bb623cb3309610961a238b3fc8ef50bff02ab7fde5e0c624335003ad8bb5a0474e083b5d35bd25ca18f5d40a15f02e2a866ac7dc2a66489d680e648ed15917f240ef1d47819c59414c8ea41a6421bff0a641e749678924d81c0e8d3ca651bf1b455c47ccf44f609f25daae89054aba8e9be28dfd2b19c50a2cec956238c3a24316c6eb67f67b2f35e0aca5122772127e289a1f8108685baad772a3151b35f65bb34bd24eef21b0afb1d9fe94c3757d78dc54f97831750e9f9b1fb102e238e035e028df4d5f8d044dd7727f91f54ff86b4cca85ce96ad9dd57063de1a6d18446a7a4c293e42b66ec8720fbd68f03065d2270c03c855afe8d747c43ac28b86ea19a10f9db20d641f6031b46f2745dfcdd4a32f2211c48a297630dff1dcfec0f0245e91e8d7f9d87aab287710fa1284001fde4b7af62d5746cb7710208a2f642065e500bae39164ded91b23e180ad51725ae0841864db0f12c265608a0677d567b4a84bd1f34bb2b3ce5eaa800678fb9c9216709f644ee08fc7fba036031bdbf2e659cc9bb3abb1f437ef99d5b9e66b337975c69b94354102601aea9ad722cc71911e8468a44cdaba47d29d6f80cd57db849f3fc4613597190a49bab00cddf21a61c2dce64cee07680c2e9568b9d38f4db4708b2dba555ce6a58e658118f2286fe35328218b5b14ab8c2b9f184cdc73047ad4a0ccb6f794b41281707e8baca01af2a8f31ce185be8f7e249a1acbfd6de4d9817632d00874c13d6e3f567fa9f19a4c6e6611bdc88970a852e6e375d29a747fbb19766e53592a96d6efbf1abb1629912393bfa47b699801b1257e45f7468fd6fd39da23040ac195b2a0f32ed4643ea48e95138a7cd8d364b3e3400d9eea40780b006c1daa4b849abd259d74cd6185b140f6fa055fdfe6fd79ece0c198344496a2c2d6b914c063a6de61a02c65f5f86d4d490fda6c5b089639d0f14606060953c13e5d6a58b4f8750df4bbfef94e1eb9fc69c6f9bdbc585e123143c5c65cc8496b8d7a59f10b909e948aaff2c61b9d49cdfcf5c40fd4de069800f2d6a859d03bff391769394c307f7d17ed325415a9dbe48af69ee1e5e26b8f3fb82efe84465655465ec2896fb20447f670ed662fd6cc6e171d3b7833578d98afa9f4c516f7fe1a3455d9effb7de366114a92caaeb5c9948c48730341ce8fceb8ff2b630765fce8fd7a7c589d953d10f5615371989c3b97f6969d4dc0be2218bf454214ca41fbe64b79e47bf219f4787e0184103ddcf44f1016f49e5555032c5c364034442232f9d609748dc332961a8a0eb288a1b66f64400021e191733437aab2c1d2ebc33c3bc06b02b1c6aa2fc8fa5e6c431cf6546bc3be4824c479c9d0402a1f82bfc8cea7a3140b3606964ad656f4860322e08bfac5b0ac36ad9ee3f6239dfb237c05f15a27e9633b4eaa9418084ccee131f68188aca0bcc9d578421b6b9ea41944275335eb360e9399a05d6d86091c75fcb972100ba47946ed50d2841218904d83ecec3b0dc236dc1f83cdade8a268c498cdb960916b61b2eec26891981f3c58c0800b56af715640d6930137600c580b3caed8744967bd0ab441f8a83f5d87410d02c378992ac9b55be15f0b57fa92c6b41637b79cacd49504a1cdc1bed3217c7e3684553a7f4f685b691cf9f5513d368cd2e236a137528c78ee1e100641cd312c6ecfa343641d3d027fe57e6b4c5c2cc97b3ce35ea754ba8b75e888ac5063a91973cb402ebdf102da779e0db7b3ada99fc916a9427637bce504cc30e80796430f3d377a320d5a41b2ffdfc2f4df43fd6b9362a2225cc2bb974c6b7d70ab3236ae54b641e7c16f5014390d73f09f977955365bf36d8c326d9089e70f5f8963850485dbc1feef3094eccc00483ff74d11c74edef99e25ec1635e57f09105c4b715750e35a59d7f1b6121b0ba64798157095434c87fd64127b8c626e25c41ecddde62d9f9e0af01a5c51e64ee46f81dec24186f2318b296be2da34878f6f51a0b32863fa825a64799453659a73d3b507b5cb0c56544fbefa04f783d3b248d59b83af58e21b62015e8054b9e0fe7e1d776ebc27d7030e82f69eee92677ab3379ca253e380acb289e0777b51467c708e7e6b69fb95190d8aa575a698fa80abf0d2402c5bef4c0e8938210a854d2fb43d47cf1001111a7df3a54120ab61789d6623c81fb6d2a81f87933b0c28682afdedfd89e26e9c6493b071391186ffe6e9a856ada27cad2c7b954f70f1ac30126101296145d5f7eaf98f2e52f75d70dcb6fa7f1914f5ca3c498603c7d1ad68e4b925de9d313f0eefb2cf78e20e55f6408f030c6e40b9b4e9695d16ea6012aeb88f2c4fb42e7861e8f054905aab6c77d2551dc3f69b9cbbdbb9222ee5cb5a5994571ee2c8056f6c494dcd5c90928502af914ebf32bda0e3c8385e73d300dccf405235b4fb4f9a282bd8002ffccb330385aa2e1531a6c78830a92c1e89ad5b3b5f66c54a2ed960fd6cf8963b7b647d0741a6141d3497a2a9ed83b083b7ba0d7465662ff92d664f9621337f2f5beeea8f885dd5f0011b9435bd55e0e347dd4c9e182721812da6d08d35585e36b10107f5ecc31df62506ba6070f3dad7bffdd887dea93854bd68554826ef760c3a963f92341fda13951e0c07633694e8f2ef9282feb4b96a969a0123299ab2e22a5d265903b7c05b83adf52dd41fffcdba536116c454adebfbd286b962a532ffdc2c16efe2bdb025ccf37ba44aa87ad108e2c7ff183acb4be79d55bfd73734d5faa15e708c709d25ea11c381abd702a7a789a670c6d9b34c0845c0ed341b368cc9aab7237446378a30cdbf4752358e3bb79e00558286273c8daff65afce1cd7e06101f1282aee5c8c6650a3ac9516569b91a4cf8ddb4189702c5442347263e1ded28838dc711481a97f5257208d7a2b1ab96f1acc1735d9bf3af2f531654bd05b01cb1449ba2a610da9ec52fb6f4f0dd405d9ca35aff0cae237f3748c163d991c8e4132bea5fc1a4eb91fe45a0dd83a8ab08b48341dcf5576ffa1a410148863bdc5f0d1a90a0ceae5a961fdbbfda48cb47081b4fca53be29a0d07fbf20a0e31bdac5729cfe5822ec2940a127eedf5796399ba33e95ba12407de7dfcf2e3a6fb59834aab160ff896d398c490061fcd03aa38fab21a2fe02c4d9c54ed5600614565d7aa0f3ee9dc6e7f82147990b82dd49f569307dd2670c62bd8b57072d29aa84a6faee8d7935cb3ebf0bc2b7b958e3269344a3f9d6503726586118812db2ce37543bb2330c7b659fc029dd486e1bdbaed5dee0ec6714b18ef989ab76ff5dfbc82907bcc9cc7d31ee6525e8a6d90323f08aaab55f3d5e0a79e3141bd376e5b48c068f6959aee209036609d1cc58882a77b9179f4674829c4d090394e6e5efa860de261fa6d2618f291cd7e870a24d9bbc6b7e68157c08bbdd98b2e8a8efba7b19a2ed31bdab299ff099d62b49a2e90299892d2303934089282bcecbad7b4c6e6c51b39f89c85937b577cd0238718817bbf2013fdd3746c10da42d7ebf166ab8fa8a76cc289774667cbdc54c245f717409bc25bffd42164c9807889e32de41482b0b2dde3484df31ea9c6392885ae92d14c76c7a3fce6ec5f176b84130273ee4a23a5ab2154d3fcd68c13ea423e42160d4cf47e1aacf478050171a82fc6bd5548ddac00918b8535f68e12ecb27b8834c6bc116498b6e0fbafcf4ee8b54ffb5cfd1dde8a2e4c593c50060aa5b167951392bb90918dba0117ca0e63a0a7d6ddcd8e147942ad507c739a417d46197a283afec1a302603c179a91be596663fe51356ce041dbedb97fdf8875080d344ce9129b81f6ebf6f316c451a81e98fc59a3d68bae60f379b15801e11c89b16d10d2fc4b36ab556cffcc2a74b4672901ee96d61a73573f4d58f34a4a504e33702374835bc67d1558fb5f60c5276f208f0bbdb546d7cda0b9c9dfd981340dcfd5139e1eb03ba00fb86cdc8dd28f6dfa64cc0d9ffe1b93b2e14f6887b7a729414cbad45fd8e690017f38109252c3740bd5af9033bd805ca53715d288c15aa5d71866c16015102a1302b4fcb5b4c7c188d7225c5b61caf6ef7b3b2db776cc2859fbef7d98d0de6564a2448c780378ad9836dccc7459cea3280465a2090a799f4a02a9c0fc80953fc58d0241dff4995b7d883239193b1b1e688715ad961cd933f6011e5be7a3f61fcaac4b1ba68b717c4a7be3bac177c89c456d47a4bf4cef6f4f9d201bbefb34013a741e268cbb61b404cc1135673ed5093ef6730e7b407f6b0106935cee46cb580cfcf49c0297df166728c025f138dcb64dca28218047daa76391cbf11a16a991627f3ae9564809c68289cbe14022426ef81cef1be95eb0a03f272ede99bcc19d4b92389224bf148bfdbb867e67fd95c67bd09457c916717ba31cce78267a6761bf460c77cd447cdba8d12c415838c5182996adcfae4921f0215d1315332c291c6f3e7b7c3948805a00c131b58c717f3493bb046895198708aeac8e2154742a5933f5bfb69f31eb9f2621c9883b7200185266c2230c652085ecd4899e57ca900d5ef3f39f799745360833de30332102e52583f434f9ad8ac687839bbc1601fb0998a23bcd1ab7f98dd75ef9618cba4e497ad08a21395ec6ba1896cbbadfd34c460ffc9d990aa6875e93fb22d41af0bda446dbb75129f572bb37eacb2a5b49622d1c29ee890808a60fccfab50e820d7b33ebf30cad47708e1e1bfb1f0ffae1a57dc1ad1bddae854bc62d84b7585a67f71c8ffa1e14f817728e9735353f1563c71e46e7f5eca929fbed6ccd4f3b5c4297b48a42dcca7ce233b0eec5b04bf68dc4539ee92ae8d55ee7a4cf59b57d449d3ce1185b4f26676d56de2cb5924f158a9290bdd49bc64173e7f60c71c0e8260a1e5c7b6bc639a5ba63d7bb8a3449fc2386b52e98611d25d68d3911806b919209fbba8dd99cbea7264900ba0ca02d5e30430c6596b87df5900688d70c694cb7a1b438d1cf1ffb42b5649631f8a1d9f7a34e7e1210351599258daf0d441d4f149eaf541b69d120dbdc86f26564f20f581347a12e822a7be82f1498d33d36a6b75b6c26bff1fa41065ac38b467e09f9e308e00a94141c953fdcbd1b559798e5c4c05da8f0cc920a83b0776101045877de744cf6738cb88e4d298244eda3755ee3bf47ce287245862df7949a0b94e956fc7e6f530f3b8df7efbb97145f3deb1175724193cb9dc1145f6291cd085964e27cd7d3beff28b85e5a988b9b9681ce0425cfbf8a0a9d72d125cfd71bcdbe18951ebf977110468eba3fe15a83009cbbc1585fc95e18ae9aefa5e220523f3b06099e35a34ed59bf748f0cb8483bd712cd9ff9d40a50de6a4cb1a5b61f97392eec913dc4eb6f3304caeb9a34d3ae8e543b9581610beb78d249ab41630de0942281d6be305feeb81f6cd5b4aa7b9662b453b21e0aa3df61c33a8173766bc246a47084df99e95759cd235f2fbb06cabd4cc1db2d47e5727f606174fca0006ea00310edb9cff87e8ae0886fd5fc109c6c4589659cd146b64277e5ee46511e3d6e4e547879fa4f4ee6bf0df69611e8e2d9af3e3d8dcb002c502caaebc5698b95e947675ee00b1dbeae1574753e25743d18dfb3160ac8a7edf6d9051ff9f9aa1c2682f094cf60621443c130ccee4cad124e87cecfb2a53f5ba76d7cccd4fad4625d1b4e778673ff062ef7164a45adbbb1d689a0df7bed712dee3a04dd5f88b12546b35c84eb0c8b09bf24051282d446c77b1bef93b6d8c9c9ea1ee18363d9dbe0c454b97c1e91d5c6891c060626971727595aaeb226793c9cccefe2c313f5795000000000000000000000000000000000000000000000000070e13191e272e336dbbc4375136df3b07f7c70e639e223e";
        let kat_alpha: [u8; LAMBDA/4] = hex::decode("47dc5764266841c1af3073fcead6a13d372979e6cca0b2952b349915f54ef663126451c742fc3ffa23af8ac195431b8836cbf74db223574185cf63ba99152f02").unwrap().try_into().unwrap();
        assert_eq!(params::SIG_LEN, 4627);
        assert_eq!(_sig.len(), 4627);
        assert_eq!(_sig[0..64], kat_alpha);
    }
}