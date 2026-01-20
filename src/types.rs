use crate::params::{K, L, N};

pub type Matrix = [[[i32; N]; L]; K];
pub type Polynomial = [i32; N];
pub type PolyNTT = [i32; N];
pub type Poly64NTT = [i64; N];

pub type VecPolyNTT<const D: usize> = [[i32; N]; D];

pub type VecPoly<const D: usize> = [[i32; N]; D];