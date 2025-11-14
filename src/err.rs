use std::array::TryFromSliceError;

#[allow(dead_code)]
#[repr(u16)]
#[derive(Clone, Debug)]
pub enum MlDsaError {
    KegGenRandomSeedError,
    NTTPolySampleError,
    BoundedPolySampleError,
    MalformedShortVector,
    MalformedVectorError,
}

impl From<TryFromSliceError> for MlDsaError {
    fn from(_: TryFromSliceError) -> Self {
        MlDsaError::MalformedVectorError
    }
}
