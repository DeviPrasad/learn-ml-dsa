mod keypair;
mod err;
mod params;
mod ntt;
mod sign;
mod xpand;
mod keygen_kat;
mod sign_with_ctx_kat;
mod sign_test;
mod sign_hedged_with_ctx;
mod types;

fn main() {
    let (_, sk) = keypair::key_gen().unwrap();
    let _ = sign::sk_decode(&sk).unwrap();
    let msg = hex::decode("20a7b7e10f70496cc38220b944def699").unwrap();
    let _sig = sign::sign(&sk, &msg).unwrap();
}

#[cfg(test)]
mod tests {
    use crate::{keypair, sign};
    use crate::params::{LAMBDA, N, TAU};
    fn get_random(rnd: &mut [u8]) {
        let _ = getrandom::fill(rnd).expect("Failed to generate randomness");
    }

    #[test]
    fn test_256_keygen_dry() {
        for _ in 0..256 {
            let (_, sk) = keypair::key_gen().unwrap();
            let (rho, key, tr, s1, s2, t0) = sign::sk_decode(&sk).unwrap();
            assert_eq!(keypair::sk_encode(&rho, &key, &tr, &s1, &s2, &t0), sk);
        }
    }

    #[test]
    fn test_sample_in_ball() {
        let mut rnd = [0u8; LAMBDA/4];
        get_random(&mut rnd);
        let c = sign::sample_in_ball(&rnd);
        {
            let mut tau = 0;
            for i in 0..N {
                assert!([-1, 0, 1].contains(&c[i]));
                if c[i] == -1 || c[i] == 1 {
                    tau += 1;
                }
            }
            assert_eq!(tau, TAU);
        }
    }

    #[test]
    fn test_signing() {
        for _ in 0..32 {
            let (_, sk) = keypair::key_gen().unwrap();
            let _ = sign::sk_decode(&sk).unwrap();
            let msg = hex::decode("20a7b7e10f70496cc38220b944def69920a7b7e10f70496cc38220b944def699").unwrap();
            let _sig = sign::sign(&sk, &msg).unwrap();
        }
    }

}
