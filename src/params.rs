pub const N: usize = 256;
pub const Q: i32 = 8380417; // 2^23 - 2^13 + 1
pub const D: usize = 13;
pub const BITLEN_Q_MINUS_1: usize = 23;  /* bitlen(q-1) = 23 because 23rd bit is high in (q-1) */
pub const BITLEN_2_POW_D_MINUS_1: usize = 13; // (2^D-1) = 0b1111111111111, i.e., 13th bit is 1.


#[cfg(feature = "ML_DSA_44")]
pub const K: usize = 4;
#[cfg(feature = "ML_DSA_44")]
pub const L: usize = 4;
#[cfg(feature = "ML_DSA_44")]
pub const TAU: usize = 39; // number of 1 or -1 in polynomial c.
#[cfg(feature = "ML_DSA_44")]
pub const ETA: u8 = 2;
#[cfg(feature = "ML_DSA_44")]
pub const BITLEN_2_ETA: usize = 3; // bitlen(2*2) = bitlen(b100) = 3.
#[cfg(feature = "ML_DSA_44")]
pub const LAMBDA: usize = 128; // collision strength of c-tilda.
#[cfg(feature = "ML_DSA_44")]
pub const GAMMA1: usize = 1 << 17; // coefficient range of y = 2^17.
#[cfg(feature = "ML_DSA_44")]
pub const BITLEN_GAMMA1_MINUS_1: usize = 17;
#[cfg(feature = "ML_DSA_44")]
pub const GAMMA2: usize = (Q-1) as usize/88; // low-order rounding range.
#[cfg(feature = "ML_DSA_44")]
pub const OMEGA: usize = 80; // max-number of 1's in the hint h.
#[cfg(feature = "ML_DSA_44")]
pub const BITLEN_PACK_W1: usize = 6; // (q-1)/(2*GAMMA2) = 44, and bitlen(22-1) = 6.
#[cfg(feature = "ML_DSA_44")]
pub const BETA: usize = 78; // TAU * ETA = 39 * 2

#[cfg(feature = "ML_DSA_65")]
pub const K: usize = 6;
#[cfg(feature = "ML_DSA_65")]
pub const L: usize = 5;
#[cfg(feature = "ML_DSA_65")]
pub const TAU: usize = 49; // number of 1 or -1 in polynomial c.
#[cfg(feature = "ML_DSA_65")]
pub const ETA: u8 = 4;
#[cfg(any(feature = "ML_DSA_65"))]
const BITLEN_2_ETA: usize = 4;  // bitlen(2*4) = bitlen(b1000) = 4.
#[cfg(feature = "ML_DSA_65")]
pub const LAMBDA: usize = 192; // collision strength of c-tilda.
#[cfg(feature = "ML_DSA_65")]
pub const GAMMA1: usize = 1 << 19; // coefficient range of y = 2^19.
#[cfg(feature = "ML_DSA_65")]
pub const BITLEN_GAMMA1_MINUS_1: usize = 19;
#[cfg(feature = "ML_DSA_65")]
pub const GAMMA2: usize = (Q-1) as usize/32; // low-order rounding range.
#[cfg(feature = "ML_DSA_65")]
pub const OMEGA: usize = 55; // max-number of 1's in the hint h.
#[cfg(feature = "ML_DSA_65")]
pub const BITLEN_PACK_W1: usize = 4; // (q-1)/(2*GAMMA2) = 16, and bitlen(16-1) = 4.
#[cfg(feature = "ML_DSA_65")]
pub const BETA: usize = 196; // TAU * ETA = 49 * 4


#[cfg(feature = "ML_DSA_87")]
pub const K: usize = 8;
#[cfg(feature = "ML_DSA_87")]
pub const L: usize = 7;
#[cfg(feature = "ML_DSA_87")]
pub const TAU: usize = 60; // number of 1 or -1 in polynomial c.
#[cfg(feature = "ML_DSA_87")]
pub const ETA: u8 = 2;
#[cfg(feature = "ML_DSA_87")]
pub const BITLEN_2_ETA: usize = 3; // bitlen(2*2) = bitlen(b100) = 3.
#[cfg(feature = "ML_DSA_87")]
pub const LAMBDA: usize = 256; // collision strength of c-tilda.
#[cfg(feature = "ML_DSA_87")]
pub const GAMMA1: usize = 1 << 19; // coefficient range of y = 2^19.
#[cfg(feature = "ML_DSA_87")]
pub const BITLEN_GAMMA1_MINUS_1: usize = 19;
#[cfg(feature = "ML_DSA_87")]
pub const GAMMA2: usize = (Q-1) as usize/32; // low-order rounding range.
#[cfg(feature = "ML_DSA_87")]
pub const OMEGA: usize = 75; // max-number of 1's in the hint h.
#[cfg(feature = "ML_DSA_87")]
pub const BITLEN_PACK_W1: usize = 4; // (q-1)/(2*GAMMA2) = 16, and bitlen(16-1) = 4.
#[cfg(feature = "ML_DSA_87")]
pub const BETA: usize = 120; // TAU * ETA = 60 * 2


pub const LEN_PUBLIC_KEY: usize = 32 + 32 * K * (BITLEN_Q_MINUS_1 - D);
pub const LEN_PRIVATE_KEY: usize = 32 + 32 + 64 + 32 * ((K+L) * BITLEN_2_ETA + (D * K));
pub const LEN_ETA_PACK_POLY: usize = 32 * BITLEN_2_ETA;
pub const LEN_T0_PACK_POLY: usize = 32 * BITLEN_2_POW_D_MINUS_1; // which BTW, is same as 32 * D
pub const SIG_LEN: usize = LAMBDA/4 + L*32*(1+(BITLEN_GAMMA1_MINUS_1))+OMEGA+K;
pub const GAMMA1_COEFFICIENT_POLY_LEN: usize = 1 + BITLEN_GAMMA1_MINUS_1;
