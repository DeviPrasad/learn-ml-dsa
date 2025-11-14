use crate::err::MlDsaError;
use crate::params::{N, D, ETA, K, L, LEN_ETA_PACK_POLY, LEN_PRIVATE_KEY, LEN_T0_PACK_POLY};

pub fn sk_decode(sk: &[u8; LEN_PRIVATE_KEY]) -> Result<([u8; 32], [u8; 32], [u8; 64], [[i32; 256]; L], [[i32; 256]; K], [[i32; 256]; K]), MlDsaError> {
    assert_eq!(sk.len(), 128 + L * LEN_ETA_PACK_POLY + K * LEN_ETA_PACK_POLY + K * LEN_T0_PACK_POLY);
    let rho = sk[0..32].try_into()?;
    let key = sk[32..64].try_into()?;
    let tr = sk[64..128].try_into()?;
    let y: &[u8; L * LEN_ETA_PACK_POLY] = &sk[128..128 + L*LEN_ETA_PACK_POLY].try_into()?;
    let z: &[u8; K * LEN_ETA_PACK_POLY] = &sk[128 + L*LEN_ETA_PACK_POLY..128 + L*LEN_ETA_PACK_POLY + K*LEN_ETA_PACK_POLY].try_into()?;
    let w: &[u8; K * LEN_T0_PACK_POLY] = &sk[128 + L*LEN_ETA_PACK_POLY + K*LEN_ETA_PACK_POLY ..].try_into()?;

    let mut s1 = [[0i32; 256]; L];
    for i in 0..L {
        let t = &y[i * LEN_ETA_PACK_POLY..(i + 1) * LEN_ETA_PACK_POLY].try_into()?;
        bit_unpack_eta(&t, &mut s1[i])?;
    }

    let mut s2: [[i32; 256]; K] = [[0; 256]; K];
    for i in 0..K {
        let t = &z[i * LEN_ETA_PACK_POLY..(i + 1) * LEN_ETA_PACK_POLY].try_into()?;
        bit_unpack_eta(&t, &mut s2[i])?;
    }

    let mut t0: [[i32; 256]; K] = [[0; 256]; K];
    for i in 0..K {
        let t = &w[i * LEN_T0_PACK_POLY..(i + 1) * LEN_T0_PACK_POLY].try_into()?;
        bit_unpack_t0(&t, &mut t0[i])?;
    }
    Ok((rho, key, tr, s1, s2, t0))
}

fn contains(v: i32, a: i32) -> bool {
    (-v..v+1).contains(&a)
}
// w is bit-packed polynomial
fn bit_unpack_eta(w: &[u8; LEN_ETA_PACK_POLY], s: &mut[i32; 256]) -> Result<(), MlDsaError> {
    const _ETA_: i32 = ETA as i32;
    let mut ok = false;

    #[cfg(any(feature="ML_DSA_44", feature="ML_DSA_87"))]
    for i in 0..N as usize/8 {
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
        ok = contains(_ETA_, s[2*i+1]) && contains(_ETA_, s[2*i+1]) &&
            contains(_ETA_, s[2*i+2]) && contains(_ETA_, s[2*i+1]) &&
            contains(_ETA_, s[2*i+5]) && contains(_ETA_, s[2*i+1]) &&
            contains(_ETA_, s[2*i+7]) && contains(_ETA_, s[2*i+1]);
    }

    #[cfg(feature="ML_DSA_65")]
    for i in 0..N as usize/2 {
        s[2*i+0] = (w[i] & 0x0F) as i32;
        s[2*i+1] = (w[i] >> 4) as i32;
        s[2*i+0] = _ETA_ - s[2 * i + 0];
        s[2*i+1] = _ETA_ - s[2 * i + 1];
        ok = contains(_ETA_, s[2 * i + 0]) && contains(_ETA_, s[2 * i + 1]);
    }

    if ok {
        Ok(())
    } else {
        Err(MlDsaError::MalformedShortVector)
    }
}

// Unpack polynomial t0 with coefficients in ]-2^{D-1}, 2^{D-1}]
fn bit_unpack_t0(w: &[u8; LEN_T0_PACK_POLY], t: &mut [i32; 256]) -> Result<(), MlDsaError> {
    const _2_POW_D_MINUS_1_: i32 = (1 << D) - 1;
    for i in 0..N as usize/8 {
        t[8*i+0]  = w[13*i+0] as i32;
        t[8*i+0] |= (w[13*i+1] as i32) << 8;
        t[8*i+0] &= 0x1FFF;

        t[8*i+1]  = (w[13*i+1] as i32) >> 5;
        t[8*i+1] |= (w[13*i+2] << 3) as i32;
        t[8*i+1] |= (w[13*i+3] as i32) << 11;
        t[8*i+1] &= 0x1FFF;

        t[8*i+2]  = (w[13*i+3] >> 2) as i32;
        t[8*i+2] |= (w[13*i+4] << 6) as i32;
        t[8*i+2] &= 0x1FFF;

        t[8*i+3]  = (w[13*i+4] >> 7) as i32;
        t[8*i+3] |= (w[13*i+5] << 1) as i32;
        t[8*i+3] |= (w[13*i+6] as i32) << 9;
        t[8*i+3] &= 0x1FFF;

        t[8*i+4]  = (w[13*i+6] >> 4) as i32;
        t[8*i+4] |= (w[13*i+7] << 4) as i32;
        t[8*i+4] |= (w[13*i+8] as i32) << 12;
        t[8*i+4] &= 0x1FFF;

        t[8*i+5]  = (w[13*i+8] >> 1) as i32;
        t[8*i+5] |= (w[13*i+9] << 7) as i32;
        t[8*i+5] &= 0x1FFF;

        t[8*i+6]  = (w[13*i+9] >> 6) as i32;
        t[8*i+6] |= (w[13*i+10] << 2) as i32;
        t[8*i+6] |= (w[13*i+11] as i32) << 10;
        t[8*i+6] &= 0x1FFF;

        t[8*i+7]  = (w[13*i+11] >> 3) as i32;
        t[8*i+7] |= (w[13*i+12] << 5) as i32;
        t[8*i+7] &= 0x1FFF;

        t[8*i+0] = (1 << (D-1)) - t[8*i+0];
        t[8*i+1] = (1 << (D-1)) - t[8*i+1];
        t[8*i+2] = (1 << (D-1)) - t[8*i+2];
        t[8*i+3] = (1 << (D-1)) - t[8*i+3];
        t[8*i+4] = (1 << (D-1)) - t[8*i+4];
        t[8*i+5] = (1 << (D-1)) - t[8*i+5];
        t[8*i+6] = (1 << (D-1)) - t[8*i+6];
        t[8*i+7] = (1 << (D-1)) - t[8*i+7];
    }

    let mut ok = true;
    for i in 0..N as usize {
        ok = ok && contains(_2_POW_D_MINUS_1_+1, t[i]);
    }
    if ok {
        Ok(())
    } else {
        Err(MlDsaError::MalformedShortVector)
    }
}