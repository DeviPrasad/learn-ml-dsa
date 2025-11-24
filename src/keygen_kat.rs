
// NIST ML-DSA-44 ACVP KATs.
// https://github.com/usnistgov/ACVP-Server/blob/master/gen-val/json-files/ML-DSA-keyGen-FIPS204/expectedResults.json

#[cfg(test)]
mod keygen_kat_run {
    use std::fs::File;
    use std::io;
    use std::io::{BufRead, BufReader};
    use crate::params::{LEN_PRIVATE_KEY, LEN_PUBLIC_KEY};
    use crate::keypair::{key_gen_internal};

    pub fn read_mldsa_nist_acvp_kat(path: &str) -> io::Result<Vec<(String, String, String, String, String)>> {
        let file = File::open(path).expect("cannot open file");
        let reader = BufReader::new(file);

        let mut results = Vec::new();
        let mut tid = String::new();
        let mut xi = String::new();
        let mut pk = String::new();
        let mut sk = String::new();
        let mut msg = String::new();
        let mut sm = String::new();

        let mut new_kat = false;
        for line in reader.lines() {
            let line = line?;
            let line = line.trim();

            if let Some(_) = line.strip_prefix("count = ") {
                new_kat = true;
            }
            if let Some(_) = line.strip_prefix("tid = ") {
                new_kat = true;
            }
            if line.is_empty() {
                new_kat = true;
            }
            if new_kat {
                if !xi.is_empty() {
                    results.push((xi.clone(), pk.clone(), sk.clone(), msg.clone(), sm.clone()));
                    tid.clear();
                    xi.clear();
                    pk.clear();
                    sk.clear();
                    msg.clear();
                    sm.clear();
                }
                new_kat = false;
            }

            if line.is_empty() {
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
            } else if let Some(val) = line.strip_prefix("msg = ") {
                msg = val.trim().to_string();
            } else if let Some(val) = line.strip_prefix("sm = ") {
                sm = val.trim().to_string();
            }
        }
        if !tid.is_empty() {
            results.push((xi, pk, sk, msg, sm));
        }
        Ok(results)
    }

    fn run_kat(kats: Vec<(String, String, String, String, String)>) {
        assert!(!kats.is_empty());
        for (kat_xi, kat_pk, kat_sk, _, _) in kats.iter() {
            let xi: [u8; 32] = hex::decode(&kat_xi).unwrap().try_into().unwrap();
            let kat_pk: [u8; LEN_PUBLIC_KEY] = hex::decode(&kat_pk).unwrap().try_into().unwrap();
            let kat_sk: [u8; LEN_PRIVATE_KEY] = hex::decode(&kat_sk).unwrap().try_into().unwrap();
            let (pk, sk) = key_gen_internal(&xi);
            assert_eq!(pk, kat_pk);
            assert_eq!(sk, kat_sk);
        }
    }
    #[cfg(any(feature = "ML_DSA_44", feature = "ML_DSA_65", feature = "ML_DSA_87"))]
    #[cfg(test)]
    mod nist_acvp_ml_dsa_keygen_kat {
        #[cfg(feature = "ML_DSA_44")]
        #[test]
        fn key_gen_acvp_nist_ml_dsa_44_tests() {
            let kat_list = super::read_mldsa_nist_acvp_kat("./kats/nist-acvp-keygen44-kats.txt").unwrap();
            assert!(!kat_list.is_empty());
            super::run_kat(kat_list);
        }

        #[cfg(feature = "ML_DSA_65")]
        #[test]
        fn key_gen_acvp_nist_ml_dsa_65_tests() {
            let kat_list = super::read_mldsa_nist_acvp_kat("./kats/nist-acvp-keygen65-kats.txt").unwrap();
            super::run_kat(kat_list)
        }

        #[cfg(feature = "ML_DSA_87")]
        #[test]
        fn key_gen_acvp_nist_ml_dsa_87_tests() {
            let kat_list = super::read_mldsa_nist_acvp_kat("./kats/nist-acvp-keygen87-kats.txt").unwrap();
            super::run_kat(kat_list)
        }
    }

    // #[cfg(any(feature = "ML_DSA_44", feature = "ML_DSA_65", feature = "ML_DSA_87"))]
    #[cfg(any(feature = "ML_DSA_44", feature = "ML_DSA_65", feature = "ML_DSA_87"))]
    #[cfg(test)]
    mod misc_ml_dsa_signature_kat {
        #[cfg(feature = "ML_DSA_44")]
        #[test]
        fn key_gen_misc_ml_dsa_44_tests() {
            let kat_list = super::read_mldsa_nist_acvp_kat("./kats/misc-keygen44-kats.txt").unwrap();
            super::run_kat(kat_list);
        }

        #[cfg(feature = "ML_DSA_65")]
        #[test]
        fn key_gen_misc_ml_dsa_65_tests() {
            let kat_list = super::read_mldsa_nist_acvp_kat("./kats/misc-keygen65-kats.txt").unwrap();
            super::run_kat(kat_list);
        }

        #[cfg(feature = "ML_DSA_87")]
        #[test]
        fn key_gen_misc_ml_dsa_87_tests() {
            let kat_list = super::read_mldsa_nist_acvp_kat("./kats/misc-keygen87-kats.txt").unwrap();
            assert!(!kat_list.is_empty());
            super::run_kat(kat_list);
        }
    }
}