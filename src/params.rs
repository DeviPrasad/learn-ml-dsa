pub const N: i32 = 256;
pub const Q: i32 = 8380417; // 2^23 - 2^13 + 1
pub const D: usize = 13;
pub const BITLEN_Q_MINUS_1: usize = 23;  /* bitlen(q-1) = 23 because 23rd bit is high in (q-1) */
pub const BITLEN_2_POW_D_MINUS_1: usize = 13; // (2^D-1) = 0b1111111111111, i.e., 13th bit is 1.


#[cfg(feature = "ML_DSA_44")]
pub const K: usize = 4;
#[cfg(feature = "ML_DSA_44")]
pub const L: usize = 4;
#[cfg(feature = "ML_DSA_44")]
pub const ETA: u8 = 2;
#[cfg(feature = "ML_DSA_44")]
pub const BITLEN_2_ETA: usize = 3; // bitlen(2*2) = bitlen(b100) = 3


#[cfg(feature = "ML_DSA_65")]
pub const K: usize = 6;
#[cfg(feature = "ML_DSA_65")]
pub const L: usize = 5;
#[cfg(feature = "ML_DSA_65")]
pub const ETA: u8 = 4;
#[cfg(any(feature = "ML_DSA_65"))]
const BITLEN_2_ETA: usize = 4;  // bitlen(2*4) = bitlen(b1000) = 4

#[cfg(feature = "ML_DSA_87")]
pub const K: usize = 8;
#[cfg(feature = "ML_DSA_87")]
pub const L: usize = 7;
#[cfg(feature = "ML_DSA_87")]
pub const ETA: u8 = 2;
#[cfg(feature = "ML_DSA_87")]
pub const BITLEN_2_ETA: usize = 3; // bitlen(2*2) = bitlen(b100) = 3


pub const LEN_PUBLIC_KEY: usize = 32 + 32 * K * (BITLEN_Q_MINUS_1 - D);
pub const LEN_PRIVATE_KEY: usize = 32 + 32 + 64 + 32 * ((K+L) * BITLEN_2_ETA + (D * K));
pub const LEN_ETA_PACK_POLY: usize = 32 * BITLEN_2_ETA;
pub const LEN_T0_PACK_POLY: usize = 32 * BITLEN_2_POW_D_MINUS_1; // which BTW, is same as 32 * D
