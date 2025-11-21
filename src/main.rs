mod keypair;
mod err;
mod params;
mod ntt;
mod sign;
mod xpand;

fn main() {
    let (_, sk) = keypair::key_gen().unwrap();
    let _ = sign::sk_decode(&sk).unwrap();
    let msg = hex::decode("20a7b7e10f70496cc38220b944def699").unwrap();
    let _sig = sign::sign(&sk, &msg, &[]).unwrap();
}

#[cfg(test)]
mod tests {
    use crate::{keypair, sign};
    use crate::params::{LAMBDA, TAU};
    fn get_random(rnd: &mut [u8]) {
        let _ = getrandom::fill(rnd).expect("Failed to generate randomness");
    }

    #[test]
    fn dry_512_keygen_test() {
        for _ in 0..128 {
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
            for i in 0..256 {
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
        for _ in 0..1 {
            let (_, sk) = keypair::key_gen().unwrap();
            let _ = sign::sk_decode(&sk).unwrap();
            let msg = hex::decode("20a7b7e10f70496cc38220b944def699").unwrap();
            let _sig = sign::sign(&sk, &msg, &"aa_fiu_v2.0_FI_request".as_bytes()).unwrap();
        }
    }

}
