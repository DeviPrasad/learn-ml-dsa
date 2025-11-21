use std::array::TryFromSliceError;

#[allow(dead_code)]
#[repr(u16)]
#[derive(Clone, Debug)]
pub enum MlDsaError {
    RandomSeedGenError,
    BadNTTPolySample,
    BadBoundedPolySample,
    MalformedShortVector,
    MalformedVector,
    SignCtxLenTooLong,
    SignatureAborted,
}

impl From<TryFromSliceError> for MlDsaError {
    fn from(_: TryFromSliceError) -> Self {
        MlDsaError::MalformedVector
    }
}
