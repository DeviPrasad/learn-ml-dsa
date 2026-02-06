use crate::params::{N, Q, ZETA};
use crate::types::{Poly64NTT, PolyNTT, Polynomial, VecPoly, VecPolyNTT};

const ML_DSA_ZETAS: [i32; N] = [
    0, 4808194, 3765607, 3761513, 5178923, 5496691, 5234739, 5178987, 7778734, 3542485, 2682288,
    2129892, 3764867, 7375178, 557458, 7159240, 5010068, 4317364, 2663378, 6705802, 4855975,
    7946292, 676590, 7044481, 5152541, 1714295, 2453983, 1460718, 7737789, 4795319, 2815639,
    2283733, 3602218, 3182878, 2740543, 4793971, 5269599, 2101410, 3704823, 1159875, 394148,
    928749, 1095468, 4874037, 2071829, 4361428, 3241972, 2156050, 3415069, 1759347, 7562881,
    4805951, 3756790, 6444618, 6663429, 4430364, 5483103, 3192354, 556856, 3870317, 2917338,
    1853806, 3345963, 1858416, 3073009, 1277625, 5744944, 3852015, 4183372, 5157610, 5258977,
    8106357, 2508980, 2028118, 1937570, 4564692, 2811291, 5396636, 7270901, 4158088, 1528066,
    482649, 1148858, 5418153, 7814814, 169688, 2462444, 5046034, 4213992, 4892034, 1987814,
    5183169, 1736313, 235407, 5130263, 3258457, 5801164, 1787943, 5989328, 6125690, 3482206,
    4197502, 7080401, 6018354, 7062739, 2461387, 3035980, 621164, 3901472, 7153756, 2925816,
    3374250, 1356448, 5604662, 2683270, 5601629, 4912752, 2312838, 7727142, 7921254, 348812,
    8052569, 1011223, 6026202, 4561790, 6458164, 6143691, 1744507, 1753, 6444997, 5720892, 6924527,
    2660408, 6600190, 8321269, 2772600, 1182243, 87208, 636927, 4415111, 4423672, 6084020, 5095502,
    4663471, 8352605, 822541, 1009365, 5926272, 6400920, 1596822, 4423473, 4620952, 6695264,
    4969849, 2678278, 4611469, 4829411, 635956, 8129971, 5925040, 4234153, 6607829, 2192938,
    6653329, 2387513, 4768667, 8111961, 5199961, 3747250, 2296099, 1239911, 4541938, 3195676,
    2642980, 1254190, 8368000, 2998219, 141835, 8291116, 2513018, 7025525, 613238, 7070156,
    6161950, 7921677, 6458423, 4040196, 4908348, 2039144, 6500539, 7561656, 6201452, 6757063,
    2105286, 6006015, 6346610, 586241, 7200804, 527981, 5637006, 6903432, 1994046, 2491325,
    6987258, 507927, 7192532, 7655613, 6545891, 5346675, 8041997, 2647994, 3009748, 5767564,
    4148469, 749577, 4357667, 3980599, 2569011, 6764887, 1723229, 1665318, 2028038, 1163598,
    5011144, 3994671, 8368538, 7009900, 3020393, 3363542, 214880, 545376, 7609976, 3105558,
    7277073, 508145, 7826699, 860144, 3430436, 140244, 6866265, 6195333, 3123762, 2358373, 6187330,
    5365997, 6663603, 2926054, 7987710, 8077412, 3531229, 4405932, 4606686, 1900052, 7598542,
    1054478, 7648983,
];

// FIPS 204 modulus Q, also used in CRYSTALS-Kyber and CRYSTALS-Dilithium.
// Q = 2^23 - 2^13 + 1
// Precomputed constant for the reduction: 2^13 - 1
const SHIFT13_MINUS_1: i64 = (1 << 13) - 1;

// Precomputed mask for the lower 23 bits: (1 << 23) - 1
const MASK_23_BITS: i64 = (1 << 23) - 1;

// Note `2^23 ≡ 2^13 - 1 (mod Q).
// The basic idea is to represent a number `y` as h * 2^23 + l modo Q, which is
// congruent to h * (2^13 - 1) + l. This is applied repeatedly.
#[inline]
pub fn mod_q(x: i64) -> i32 {
    // The reduction: f(y) = (y >> 23) * (2^13 - 1) + (y & MASK_23_BITS)
    // 1st reduction:
    // Input `x` can be up to ~48-49 bits; result `t1` is reduced to ~40 bits.
    let high1 = x >> 23;
    let low1 = x & MASK_23_BITS;
    let t1 = high1 * SHIFT13_MINUS_1 + low1;

    // 2nd reduction:
    // Input `t1` is ~40 bits. The result `t2` is reduced to ~30 bits.
    let high2 = t1 >> 23;
    let low2 = t1 & MASK_23_BITS;
    let t2 = high2 * SHIFT13_MINUS_1 + low2;

    // 3rd reduction:
    // Input `t2` is ~30 bits. The result `t3` is guaranteed to be in a small
    // range, specifically [0, Q + d) where d is small. This is less than 2*Q.
    let high3 = t2 >> 23;
    let low3 = t2 & MASK_23_BITS;
    let t3 = high3 * SHIFT13_MINUS_1 + low3;

    // if result >= Q { result - Q } else { result }
    let r = t3 as i32 - Q;
    r + ((r >> 31) & Q)
}

fn pow_mod_q(n: i32, exp: u8) -> i32 {
    const Q64: i64 = Q as i64;
    let mut result = 1i64;
    let mut exp = exp;
    let mut num = n as i64;
    while exp > 0 {
        if exp & 1 == 1 {
            let _dbg_res = (result * num) % Q64;
            result = mod_q(result * num) as i64;
            assert_eq!(_dbg_res, result);
        }
        exp >>= 1;
        let _dbg_res = (num * num) % Q64;
        num = mod_q(num * num) as i64;
        assert_eq!(_dbg_res, num);
    }
    result as i32
}

pub fn ntt(w: &Polynomial) -> PolyNTT {
    fn _ntt_impl_(wh: &mut Poly64NTT) -> PolyNTT {
        let mut i = 0;
        let mut len = 128_usize;
        while len >= 1 {
            let mut start = 0_usize;
            while start < N {
                i += 1;
                debug_assert_eq!(pow_mod_q(ZETA, (i as u8).reverse_bits()), ML_DSA_ZETAS[i]);
                let z = ML_DSA_ZETAS[i] as i64;
                for j in start..start + len {
                    let t = mod_q(z * wh[j + len]) as i64;
                    wh[j + len] = wh[j] - t;
                    wh[j] = wh[j] + t;
                }
                start += 2 * len;
            }
            len = len / 2
        }
        wh.map(|x| mod_q(x))
    }
    _ntt_impl_(&mut w.map(|x| x as i64))
}

pub fn ntt_inverse(wh: &PolyNTT) -> Polynomial {
    fn _ntt_inverse_impl_(w: &mut Poly64NTT) -> Polynomial {
        let mut m = N;
        let mut len = 1_usize;
        while len < N {
            let mut start = 0_usize;
            while start < N {
                m -= 1;
                let z = -ML_DSA_ZETAS[m] as i64;
                for j in start..start + len {
                    let t = w[j];
                    w[j] = t + w[j + len];
                    let w_j_plus_len = t - w[j + len];
                    w[j + len] = mod_q(z * w_j_plus_len) as i64;
                }
                start += 2 * len;
            }
            // NOTE: THIS IS AN INEFFICIENT METHOD, USED ONLY FOR DEMONSTARTING THE PRINCIPLE.
            // divide by 2 (or equivalently multiply by inverse(2)) in each step.
            for i in 0..N {
                w[i] = mod_q(w[i] * 4190209) as i64;
            }
            len *= 2;
        }
        // This is the stanndard (efficient) method: multiply by the inverse(256) in the end.
        // w.map(|x| mod_q(mod_q(x) as i64 * 8347681))
        w.map(|x| x as i32)
    }
    _ntt_inverse_impl_(&mut wh.map(|x| x as i64))
}

pub fn ntt_add(ah: &PolyNTT, bh: &PolyNTT) -> PolyNTT {
    std::array::from_fn(|i| mod_q((ah[i] + bh[i]) as i64))
}

pub fn ntt_sub(ah: &PolyNTT, bh: &PolyNTT) -> PolyNTT {
    std::array::from_fn(|i| mod_q((ah[i] - bh[i]) as i64))
}

pub fn ntt_neg_vec<const D: usize>(v: &VecPolyNTT<D>) -> [PolyNTT; D] {
    let mut r: VecPolyNTT<D> = [[0i32; _]; _];
    for i in 0..D {
        r[i].copy_from_slice(&ntt_neg(&v[i]));
    }
    r
}

pub fn ntt_neg(w: &PolyNTT) -> PolyNTT {
    std::array::from_fn(|i| mod_q(-w[i] as i64))
}

pub fn ntt_multiply(ah: &PolyNTT, bh: &PolyNTT) -> PolyNTT {
    std::array::from_fn(|i| mod_q(ah[i] as i64 * bh[i] as i64))
}

pub fn vec_ntt_multiply<const D: usize>(ph: &PolyNTT, vh: &VecPolyNTT<D>, rh: &mut VecPoly<D>) {
    for k in 0..D {
        rh[k].copy_from_slice(&ntt_inverse(&ntt_multiply(&ph, &vh[k])));
    }
}

pub fn poly_sub(ah: &PolyNTT, bh: &PolyNTT) -> PolyNTT {
    std::array::from_fn(|i| mod_q((ah[i] - bh[i]) as i64))
}

pub fn poly_add(ah: &PolyNTT, bh: &PolyNTT) -> PolyNTT {
    std::array::from_fn(|i| mod_q((ah[i] + bh[i]) as i64))
}

#[cfg(test)]
mod ntt_tests {
    use crate::ntt::{mod_q, ntt, ntt_add, ntt_inverse, ntt_multiply, ML_DSA_ZETAS};
    use crate::params::{N, Q};
    use crate::types::Polynomial;

    #[test]
    fn test_zetas_less_than_q() {
        for i in 0..N {
            assert!(ML_DSA_ZETAS[i] < Q);
        }
    }
    #[test]
    fn test_ntt_conv_01() {
        // 1 + 200x + 300x^2 + 0^254 + 0^255;
        let mut w: Polynomial = [0i32; N];
        w[0] = mod_q(Q as i64 - 1);
        w[1] = 200;
        w[2] = 300;
        w[3] = mod_q(Q as i64);
        for i in 4..253 {
            w[i] = mod_q(getrandom::u32().unwrap() as i64);
        }
        w[253] = mod_q(-1);
        w[254] = 254;
        w[255] = Q - 2;

        let wh = ntt(&w);
        let _w = ntt_inverse(&wh);
        assert_eq!(w, _w);
    }

    #[test]
    fn test_ntt_add_01() {
        let a: Polynomial = std::array::from_fn(|_| mod_q(getrandom::u32().unwrap() as i64));
        let b: Polynomial = std::array::from_fn(|_| mod_q(getrandom::u32().unwrap() as i64));
        let r0: Polynomial = ntt_inverse(&ntt_add(&ntt(&a), &ntt(&b)));
        let r1 = std::array::from_fn(|i| mod_q((a[i] + b[i]) as i64));
        assert_eq!(r0, r1);
    }

    #[test]
    fn test_ntt_mult_01() {
        let a: Polynomial = std::array::from_fn(|_| mod_q(getrandom::u32().unwrap() as i64));
        let mut b: Polynomial = [0i32; N];
        b[0] = 256;
        let expected: Polynomial = std::array::from_fn(|i| mod_q((256 * a[i]) as i64));
        for i in 1..N {
            assert!(a[i] >= 0 && a[i] < Q)
        }
        for i in 1..N {
            assert!(b[i] >= 0 && b[i] < Q)
        }
        let r0 = ntt_inverse(&ntt_multiply(&ntt(&a), &ntt(&b)));
        assert_eq!(r0, expected);
    }
}
