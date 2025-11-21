use crate::err::MlDsaError;
use crate::ntt::{ntt, ntt_add, ntt_inverse, ntt_multiply};
use crate::params::{N, Q, D, ETA, K, L, LEN_ETA_PACK_POLY, LEN_PRIVATE_KEY, LEN_PUBLIC_KEY, LEN_T0_PACK_POLY};
use getrandom;
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use crate::xpand::expand_a;

// https://github.com/post-quantum-cryptography/KAT/blob/main/MLDSA/kat_MLDSA_44_hedged_raw.rsp
pub fn key_gen() -> Result<([u8; LEN_PUBLIC_KEY], [u8; LEN_PRIVATE_KEY]), MlDsaError>{
    let mut xi = [0u8; 32];
    let _ = getrandom::fill(&mut xi)
        .map_err(|_| MlDsaError::RandomSeedGenError);
    Ok(key_gen_internal(&xi))
}

fn key_gen_internal(xi: &[u8; 32]) -> ([u8; LEN_PUBLIC_KEY], [u8; LEN_PRIVATE_KEY]) {
    let mut rho = [0u8; 32];
    let mut rho_prime = [0u8; 64];
    let mut key = [0u8; 32];

    {
        let mut d = [0u8; 34];
        d[..32].copy_from_slice(xi);
        d[32] = K as u8;
        d[33] = L as u8;

        let mut seed = [0u8; 128];
        Shake256::digest_xof(&d, &mut seed);

        rho.copy_from_slice(&seed[0..32]);
        rho_prime.copy_from_slice(&seed[32..96]);
        key.copy_from_slice(&seed[96..128]);
    }

    // Expand matrix.
    let a_hat = expand_a(&rho);
    // Sample short vectors s1 and s2.
    let (s1, s2): ([[i32; 256]; L], [[i32; 256]; K]) = expand_s(&rho_prime);
    let mut s1_hat = [[0i32; 256]; L];
    for (i, cv) in s1.iter().enumerate() {
        s1_hat[i] = ntt(&cv);
    }
    // multiply a_hat and s1_hat
    let mut prod_a_s1 =  [[0i32; 256]; K];
    for (i, a_row_poly) in a_hat.iter().enumerate() {
        for (k, s1_col_poly) in s1_hat.iter().enumerate() {
            prod_a_s1[i] = ntt_add(&prod_a_s1[i], &ntt_multiply(&a_row_poly[k], &s1_col_poly))
        }
    }
    let r = prod_a_s1.map(|v| ntt_inverse(&v));
    let t: [[i32; 256]; K] = std::array::from_fn(|i| vec_add(&r[i], &s2[i]));

    let mut t0 = [[0i32; 256]; K];
    let mut t1 = [[0i32; 256]; K];
    for (i, t_vec) in t.iter().enumerate() {
        for j in 0..256 {
            (t1[i][j], t0[i][j]) = power2round(t_vec[j])
        }
        for j in 0..256 {
            assert!(t1[i][j] < 1024_i32);
        }
        for j in 0..256 {
            assert!(t0[i][j] <= (1 << (D-1)));
        }
    }

    let pk = pk_encode(&rho, &t1);

    // public random
    let mut tr = [0u8; 64];
    Shake256::digest_xof(&pk, &mut  tr);
    let sk = sk_encode(&rho, &key, &tr, &s1, &s2, &t0);
    (pk, sk)
}

// sizes in bytes of public keys:
// ML_DSA_44 1312 // 32 + 32*4+10
// ML_DSA_65 1952 // 32 + 32*6+10
// ML_DSA_87 2592 // 32 + 32*8+10

fn pk_encode(rho: &[u8; 32], t1: &[[i32; 256]; K]) -> [u8; LEN_PUBLIC_KEY]{
    let mut pk = [0u8; LEN_PUBLIC_KEY];
    pk[0..32].copy_from_slice(rho);
    for i in 0..K {
        // hidden constant parameter: B = (2^(bitlen(q-1)-d) - 1) = (2^10 - 1) = 1023.
        let t_packed = simple_bit_pack_t1(&t1[i]);
        pk[32+(i*320)..][..320].copy_from_slice(&t_packed)
    }
    pk
}

pub fn sk_encode(rho: &[u8; 32], key: &[u8; 32], tr: &[u8; 64],
             s1: &[[i32; 256]; L], s2: &[[i32; 256]; K],
             t0: &[[i32; 256]; K]) -> [u8; LEN_PRIVATE_KEY] {
    let mut sk = [0u8; LEN_PRIVATE_KEY];
    sk[0..32].copy_from_slice(rho);
    sk[32..64].copy_from_slice(key);
    sk[64..128].copy_from_slice(tr);
    for i in 0..L {
        let j = 128 + i*LEN_ETA_PACK_POLY;
        sk[j..][..LEN_ETA_PACK_POLY].copy_from_slice(&bit_pack_eta(&s1[i]));
    }
    let j = 128 + L*LEN_ETA_PACK_POLY;
    for i in 0..K {
        let k = j + i*LEN_ETA_PACK_POLY;
        sk[k..][..LEN_ETA_PACK_POLY].copy_from_slice(&bit_pack_eta(&s2[i]))
    }
    let j = 128 + L*LEN_ETA_PACK_POLY + K*LEN_ETA_PACK_POLY;
    for i in 0..K {
        let k = j + i*LEN_T0_PACK_POLY;
        sk[k..][..LEN_T0_PACK_POLY].copy_from_slice(&bit_pack_t0(&t0[i]))
    }
    sk
}

fn bit_pack_t0(t0: &[i32; 256]) -> [u8; LEN_T0_PACK_POLY] {
    let mut t = [0i32; 8];
    let mut r = [0i32; LEN_T0_PACK_POLY];

    for i in 0..N/8 {
        t[0] = (1 << (D-1)) - t0[8*i+0];
        t[1] = (1 << (D-1)) - t0[8*i+1];
        t[2] = (1 << (D-1)) - t0[8*i+2];
        t[3] = (1 << (D-1)) - t0[8*i+3];
        t[4] = (1 << (D-1)) - t0[8*i+4];
        t[5] = (1 << (D-1)) - t0[8*i+5];
        t[6] = (1 << (D-1)) - t0[8*i+6];
        t[7] = (1 << (D-1)) - t0[8*i+7];

        r[13*i+ 0]  =  t[0];
        r[13*i+ 1]  =  t[0] >>  8;
        r[13*i+ 1] |=  t[1] <<  5;
        r[13*i+ 2]  =  t[1] >>  3;
        r[13*i+ 3]  =  t[1] >> 11;
        r[13*i+ 3] |=  t[2] <<  2;
        r[13*i+ 4]  =  t[2] >>  6;
        r[13*i+ 4] |=  t[3] <<  7;
        r[13*i+ 5]  =  t[3] >>  1;
        r[13*i+ 6]  =  t[3] >>  9;
        r[13*i+ 6] |=  t[4] <<  4;
        r[13*i+ 7]  =  t[4] >>  4;
        r[13*i+ 8]  =  t[4] >> 12;
        r[13*i+ 8] |=  t[5] <<  1;
        r[13*i+ 9]  =  t[5] >>  7;
        r[13*i+ 9] |=  t[6] <<  6;
        r[13*i+10]  =  t[6] >>  2;
        r[13*i+11]  =  t[6] >> 10;
        r[13*i+11] |=  t[7] <<  3;
        r[13*i+12]  =  t[7] >>  5;
    }

    r.map(|x| x as u8)
}

#[inline]
fn bit_pack_eta(w: &[i32; 256]) -> [u8; LEN_ETA_PACK_POLY] {
    let mut t = [0i32; 8];
    let mut r = [0i32; LEN_ETA_PACK_POLY];

    for i in 0..256 {
        if !(w[i] >= -(ETA as i32) && w[i] <= ETA as i32) {
            assert!(w[i] >= -(ETA as i32) && w[i] <= ETA as i32)
        }
    }

    #[cfg(any(feature="ML_DSA_44", feature="ML_DSA_87"))]
    for i in 0..N/8 {
        t[0] = ETA as i32 - w[8*i+0];
        t[1] = ETA as i32 - w[8*i+1];
        t[2] = ETA as i32 - w[8*i+2];
        t[3] = ETA as i32 - w[8*i+3];
        t[4] = ETA as i32 - w[8*i+4];
        t[5] = ETA as i32 - w[8*i+5];
        t[6] = ETA as i32 - w[8*i+6];
        t[7] = ETA as i32 - w[8*i+7];

        r[3*i+0]  = (t[0] >> 0) | (t[1] << 3) | (t[2] << 6);
        r[3*i+1]  = (t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7);
        r[3*i+2]  = (t[5] >> 1) | (t[6] << 2) | (t[7] << 5);
    }

    #[cfg(feature="ML_DSA_65")]
    for i in 0.. N as usize/2 {
        t[0] = ETA as i32 - w[2*i+0];
        t[1] = ETA as i32 - w[2*i+1];
        r[i] = t[0] | (t[1] << 4);
    }

    r.map(|x| x as u8)
}


// simple_bit_pack encodes a polynomial into a byte string.
// simple_bit_pack_t1 is specifically tailored to pack public key component t1's
//   coefficients into groups of 1023 bits.
// input: B = 1023 and w= \in R such that coefficients are all in [0, 1023].
#[inline]
fn simple_bit_pack_t1(a: &[i32; 256]) -> [u8; 32 * 10] {
    const B: u32 = 1023;
    for i in 0..N/4 {
        assert!((a[4*i+0] as u32) <= B);
        assert!((a[4*i+1] as u32) <= B);
        assert!((a[4*i+2] as u32) <= B);
        assert!((a[4*i+3] as u32) <= B);
    }

    let mut r = [0u8; 32 * 10];
    for i in 0..N/4 {
        r[5*i+0] = ((a[4*i+0]) >> 0) as u8;
        r[5*i+1] = (((a[4*i+0]) >> 8) | ((a[4*i+1]) << 2)) as u8;
        r[5*i+2] = (((a[4*i+1]) >> 6) | ((a[4*i+2]) << 4)) as u8;
        r[5*i+3] = (((a[4*i+2]) >> 4) | ((a[4*i+3]) << 6)) as u8;
        r[5*i+4] = ((a[4*i+3]) >> 2) as u8;
    }

    for i in 0..N/4 {
        assert!((r[5*i+0] as u32) < B);
        assert!((r[5*i+1] as u32) < B);
        assert!((r[5*i+2] as u32) < B);
        assert!((r[5*i+3] as u32) < B);
        assert!((r[5*i+4] as u32) < B);
    }
    r
}

fn vec_add(a: &[i32; 256], b: &[i32; 256]) -> [i32; 256] {
    ntt_add(a, b)
}


fn expand_s(r: &[u8; 64]) -> ([[i32; 256]; L], [[i32; 256]; K]) {
    let mut s1 = [[0i32; 256]; L];
    let mut s2 = [[0i32; 256]; K];
    let mut rho = [0u8; 66];
    rho[0..64].copy_from_slice(r);
    rho[65] = 0; // rho[64..65] = IntegerToBytes(r, 2)
    for r in 0..L {
        rho[64] = r as u8;
        s1[r] = reg_bounded_poly(rho)
    }
    for r in 0..K {
        rho[64] = (r + L) as u8;
        s2[r] = reg_bounded_poly(rho);
    }
    (s1, s2)
}

// q = 2^23 - 2^13 + 1 = 8380417.
// returns an element in Z_q which fits in 3 bytes.
// this is a polynomial in the NTT form.
pub(crate) fn reg_ntt_poly(seed:[u8; 34]) -> [i32; 256] {
    let mut j = 0usize;
    let mut g = sha3::Shake128::default();
    g.update(&seed);
    let mut xof = g.finalize_xof();
    let mut poly = [0i32; 256];
    while j < 256 {
        let mut s = [0u8; 3];
        xof.read(&mut s);
        if let Ok(z_q) = coefficient_from_three_bytes(s[0], s[1], s[2]) {
            assert!(z_q < Q);
            poly[j] = z_q; // mod_q(z_q as i64);
            j += 1;
        }
    }
    poly
}

fn reg_bounded_poly(rho: [u8; 66]) -> [i32; 256] {
    let mut poly = [0i8; 256];
    let mut h = Shake256::default();
    h.update(&rho);
    let mut xof = h.finalize_xof();
    let mut j = 0usize;
    while j < 256 {
        let mut z = [0u8];
        xof.read(&mut z);
        let rz0 = coefficient_from_half_byte(z[0] & 0x0F);
        let rz1 = coefficient_from_half_byte((z[0] >> 4) & 0x0F);
        if let Ok(z0) = rz0 {
            poly[j] = z0;
            j += 1;
        }
        if let Ok(z1) = rz1 && j < 256 {
            poly[j] = z1;
            j += 1;
        }
    }
    // to_polynomial_ring(&poly)
    poly.map(|x|x as i32)
}

fn coefficient_from_half_byte(b: u8) -> Result<i8, MlDsaError> {
    #[cfg(any(feature = "ML_DSA_44", feature = "ML_DSA_87"))]
    const MOD5: [i8; 16] = [0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0];
    #[cfg(any(feature = "ML_DSA_44", feature = "ML_DSA_87"))]
    if b < 15 {
        Ok(2 - MOD5[(b & 0x0F) as usize])
    } else {
        Err(MlDsaError::BadBoundedPolySample)
    }

    #[cfg(feature = "ML_DSA_65")]
    if b < 9 {
        Ok(4 - b as i8)
    } else {
        Err(MlDsaError::BadBoundedPolySample)
    }
}

// generates an element of {0, 1, 2,..., q-1} + { NTTPolySampleError }
fn coefficient_from_three_bytes(b0: u8, b1: u8, b2: u8) -> Result<i32, MlDsaError> {
    let b2 = (b2 & 127) as i32; // set the top bit of b2 to 0, such that 0 <= b <= 127
    let z = (b2 << 16) + ((b1 as i32) << 8) + (b0 as i32);
    if z < Q {
        Ok(z)
    } else {
        Err(MlDsaError::BadNTTPolySample)
    }
}

fn power2round(r:i32) -> (i32, i32) {
    assert!(r < Q);
    let r1 = (r + (1 << (D-1)) - 1) >> D;
    let r0 = r - (r1 << D);
    assert!(r1 < 1024);
    (r1, r0)
}


// NIST ML-DSA-44 ACVP KATs.
// https://github.com/usnistgov/ACVP-Server/blob/master/gen-val/json-files/ML-DSA-keyGen-FIPS204/expectedResults.json
#[cfg(any(feature = "ML_DSA_44", feature = "ML_DSA_65", feature = "ML_DSA_87"))]
#[cfg(test)]
mod nist_acvp_ml_dsa_keygen_kats {
    use crate::keypair::{key_gen_internal, LEN_PRIVATE_KEY, LEN_PUBLIC_KEY};
    use std::fs::File;
    use std::io::{self, BufRead, BufReader};

    /// Reads a text file in the format:
    /// tcId = <hex>
    /// xi = <hex>
    /// pk = <hex>
    /// sk = <hex>
    /// (blank line between test cases)
    #[cfg(test)]
    pub fn read_mldsa_nist_acvp_kats(path: &str) -> io::Result<Vec<(String, String, String, String)>> {
        let file = File::open(path).expect("cannot open file");
        let reader = BufReader::new(file);

        let mut results = Vec::new();
        let mut tid = String::new();
        let mut xi = String::new();
        let mut pk = String::new();
        let mut sk = String::new();

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();

            if line.is_empty() {
                // End of one test case
                if !tid.is_empty() {
                    results.push((tid.clone(), xi.clone(), pk.clone(), sk.clone()));
                    tid.clear();
                    xi.clear();
                    pk.clear();
                    sk.clear();
                }
                continue;
            }

            if let Some(val) = line.strip_prefix("tid = ") {
                tid = val.trim().to_string();
            } else if let Some(val) = line.strip_prefix("xi = ") {
                xi = val.trim().to_string();
            } else if let Some(val) = line.strip_prefix("pk = ") {
                pk = val.trim().to_string();
            } else if let Some(val) = line.strip_prefix("sk = ") {
                sk = val.trim().to_string();
            }
        }
        if !tid.is_empty() {
            results.push((tid, xi, pk, sk));
        }
        Ok(results)
    }

    fn run_kats(kats: Vec<(String, String, String, String)>) {
        for (_id, kat_xi, kat_pk, kat_sk) in kats.iter() {
            let xi: [u8; 32] = hex::decode(&kat_xi).unwrap().try_into().unwrap();
            let kat_pk: [u8; LEN_PUBLIC_KEY] = hex::decode(&kat_pk).unwrap().try_into().unwrap();
            let kat_sk: [u8; LEN_PRIVATE_KEY] = hex::decode(&kat_sk).unwrap().try_into().unwrap();
            let (pk, sk) = key_gen_internal(&xi);
            assert_eq!(pk, kat_pk);
            assert_eq!(sk, kat_sk);
        }
    }
    #[cfg(feature = "ML_DSA_44")]
    #[test]
    fn key_gen_acvp_nist_ml_dsa_44_tests() {
        let kats = read_mldsa_nist_acvp_kats("./kats/nist-acvp-keygen44-kats.txt").unwrap();
        run_kats(kats)
    }

    #[cfg(feature = "ML_DSA_65")]
    #[test]
    fn key_gen_acvp_nist_ml_dsa_65_tests() {
        let kats = read_mldsa_nist_acvp_kats("./kats/nist-acvp-keygen65-kats.txt").unwrap();
        run_kats(kats)
    }

    #[cfg(feature = "ML_DSA_87")]
    #[test]
    fn key_gen_acvp_nist_ml_dsa_87_tests() {
        let kats = read_mldsa_nist_acvp_kats("./kats/nist-acvp-keygen87-kats.txt").unwrap();
        run_kats(kats)
    }
}

#[cfg(any(feature = "ML_DSA_44", feature = "ML_DSA_65", feature = "ML_DSA_87"))]
#[cfg(test)]
mod misc_ml_dsa_keygen_kats {
    use crate::keypair::{key_gen_internal, LEN_PRIVATE_KEY, LEN_PUBLIC_KEY};
    use std::fs::File;
    use std::io::{self, BufRead, BufReader};

    /// Reads a NIST ML-DSA KAT text file (no blank lines).
    /// Each test starts with `count = N`, followed by `seed`, `pk`, and `sk`.
    #[cfg(test)]
    pub fn read_mldsa_hedged_kats(path: &str) -> io::Result<Vec<(String, String, String, String)>> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        let mut results = Vec::new();
        let mut tid = String::new();
        let mut xi = String::new();
        let mut pk = String::new();
        let mut sk = String::new();

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();

            if line.is_empty() {
                continue;
            }

            if let Some(val) = line.strip_prefix("count = ") {
                // If previous test exists, push it before starting a new one
                if !tid.is_empty() {
                    results.push((tid.clone(), xi.clone(), pk.clone(), sk.clone()));
                    xi.clear();
                    pk.clear();
                    sk.clear();
                }
                tid = val.trim().to_string();
            } else if let Some(val) = line.strip_prefix("xi = ") {
                xi = val.trim().to_string();
            } else if let Some(val) = line.strip_prefix("pk = ") {
                pk = val.trim().to_string();
            } else if let Some(val) = line.strip_prefix("sk = ") {
                sk = val.trim().to_string();
            }
        }

        // Push the final entry if the file didn’t end with another "count ="
        if !tid.is_empty() {
            results.push((tid, xi, pk, sk));
        }

        Ok(results)
    }

    #[cfg(feature = "ML_DSA_44")]
    #[test]
    fn key_gen_misc_ml_dsa_44_tests() {
        let kats = read_mldsa_hedged_kats("./kats/misc-keygen44-kats.txt").unwrap();
        for (_id, kat_xi, kat_pk, kat_sk) in kats.iter() {
            let xi: [u8; 32] = hex::decode(&kat_xi).unwrap().try_into().unwrap();
            let kat_pk: [u8; LEN_PUBLIC_KEY] = hex::decode(&kat_pk).unwrap().try_into().unwrap();
            let kat_sk: [u8; LEN_PRIVATE_KEY] = hex::decode(&kat_sk).unwrap().try_into().unwrap();
            let (pk, sk) = key_gen_internal(&xi);
            assert_eq!(pk, kat_pk);
            assert_eq!(sk, kat_sk);
        }
    }

    #[cfg(feature = "ML_DSA_65")]
    #[test]
    fn key_gen_misc_ml_dsa_65_tests() {
        let kats = read_mldsa_hedged_kats("./kats/misc-keygen65-kats.txt").unwrap();
        for (_id, kat_xi, kat_pk, kat_sk) in kats.iter() {
            let xi: [u8; 32] = hex::decode(&kat_xi).unwrap().try_into().unwrap();
            let kat_pk: [u8; LEN_PUBLIC_KEY] = hex::decode(&kat_pk).unwrap().try_into().unwrap();
            let kat_sk: [u8; LEN_PRIVATE_KEY] = hex::decode(&kat_sk).unwrap().try_into().unwrap();
            let (pk, sk) = key_gen_internal(&xi);
            assert_eq!(pk, kat_pk);
            assert_eq!(sk, kat_sk);
        }
    }

    #[cfg(feature = "ML_DSA_87")]
    #[test]
    fn key_gen_misc_ml_dsa_87_tests() {
        let kats = read_mldsa_hedged_kats("./kats/misc-keygen87-kats.txt").unwrap();
        for (_id, kat_xi, kat_pk, kat_sk) in kats.iter() {
            let xi: [u8; 32] = hex::decode(&kat_xi).unwrap().try_into().unwrap();
            let kat_pk: [u8; LEN_PUBLIC_KEY] = hex::decode(&kat_pk).unwrap().try_into().unwrap();
            let kat_sk: [u8; LEN_PRIVATE_KEY] = hex::decode(&kat_sk).unwrap().try_into().unwrap();
            let (pk, sk) = key_gen_internal(&xi);
            assert_eq!(pk, kat_pk);
            assert_eq!(sk, kat_sk);
        }
    }
}
