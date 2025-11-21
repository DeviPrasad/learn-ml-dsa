use sha3::digest::{ExtendableOutput, Update};
use crate::params::{GAMMA1, GAMMA1_COEFFICIENT_POLY_LEN, K, L, N};

pub fn expand_a(seed: &[u8; 32]) -> [[[i32; 256]; L]; K] {
    let mut a_hat: [[[i32; 256]; L]; K] = [[[0i32; 256]; L]; K];
    let mut rp = [0u8; 34];
    rp[0..32].copy_from_slice(seed);
    for r in 0..K {
        for s in 0..L {
            rp[32] = s as u8;
            rp[33] = r as u8;
            let z_q = crate::keypair::reg_ntt_poly(rp);
            a_hat[r][s] = z_q;
        }
    }
    a_hat
}

// algorithm 34, page 38.
// samples a vector y \in R^l such that
// each polynomial y[r] has coefficients between -gamma_1+1 and gamma_1.
// input: a sedd rho \in B^64, and a non-negative integer mu.
// mu is a 16-bit counter incremented by l in each iteration in algorithm 7, sign_internal.
// output: vector y \in R^l
pub fn expand_mask(rho: &[u8; 64], counter: u16) -> [[i32; 256]; L] {
    let mut y = [[0i32; 256]; L];
    const C: usize = GAMMA1_COEFFICIENT_POLY_LEN;
    for l in 0..L {
        let mut v = [0u8; 32*C];
        let mut h = sha3::Shake256::default();
        h.update(rho);
        h.update(&(counter + l as u16).to_le_bytes());
        h.finalize_xof_into(&mut v);
        y[l] = bit_unpack_gamma(&v);
    }
    y
}

pub fn bit_unpack_gamma(v: &[u8; 32*GAMMA1_COEFFICIENT_POLY_LEN]) -> [i32; 256] {
    let mut r = [0i32; 256];
    #[cfg(feature="ML_DSA_44")]
    // GAMMA1 = (1 << 17)
    for i in 0..N/4 {
        r[4*i+0]  = v[9*i+0] as i32;
        r[4*i+0] |= (v[9*i+1] as i32) << 8;
        r[4*i+0] |= (v[9*i+2] as i32) << 16;
        r[4*i+0] &= 0x3FFFF;

        r[4*i+1]  = (v[9*i+2] as i32) >> 2;
        r[4*i+1] |= (v[9*i+3] as i32) << 6;
        r[4*i+1] |= (v[9*i+4] as i32) << 14;
        r[4*i+1] &= 0x3FFFF;

        r[4*i+2]  = (v[9*i+4] as i32) >> 4;
        r[4*i+2] |= (v[9*i+5] as i32) << 4;
        r[4*i+2] |= (v[9*i+6] as i32) << 12;
        r[4*i+2] &= 0x3FFFF;

        r[4*i+3]  = (v[9*i+6] as i32) >> 6;
        r[4*i+3] |= (v[9*i+7] as i32) << 2;
        r[4*i+3] |= (v[9*i+8] as i32) << 10;
        r[4*i+3] &= 0x3FFFF;

        r[4*i+0] = GAMMA1 as i32 - r[4*i+0];
        r[4*i+1] = GAMMA1 as i32 - r[4*i+1];
        r[4*i+2] = GAMMA1 as i32 - r[4*i+2];
        r[4*i+3] = GAMMA1 as i32 - r[4*i+3];
    }
    #[cfg(any(feature="ML_DSA_65", feature="ML_DSA_87"))]
    // GAMMA1 = (1 << 19)
    for i in 0..N/2{
        r[2*i+0]  = v[5*i+0] as i32;
        r[2*i+0] |= (v[5*i+1] as i32) << 8;
        r[2*i+0] |= (v[5*i+2] as i32) << 16;
        r[2*i+0] &= 0xFFFFF;

        r[2*i+1]  = (v[5*i+2] as i32) >> 4;
        r[2*i+1] |= (v[5*i+3] as i32) << 4;
        r[2*i+1] |= (v[5*i+4] as i32) << 12;
        /* r[2*i+1] &= 0xFFFFF; */ /* No effect, since we're anyway at 20 bits */

        r[2*i+0] = GAMMA1 as i32 - r[2*i+0];
        r[2*i+1] = GAMMA1 as i32 - r[2*i+1];
    }
    r
}