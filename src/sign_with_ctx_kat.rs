// deterministic signature with a context string; 'rnd' is [0u8; 32]
// https://github.com/itzmeanjan/ml-dsa/blob/master/kats/ml_dsa_44_sign.acvp.kat
// https://github.com/itzmeanjan/ml-dsa/blob/master/kats/ml_dsa_65_sign.acvp.kat
// https://github.com/itzmeanjan/ml-dsa/blob/master/kats/ml_dsa_87_sign.acvp.kat
#[cfg(test)]
mod deterministic_signature_with_ctx_kat {
        use std::fs::File;
        use std::io;
        use std::io::{BufRead, BufReader};
        use crate::params::{LEN_PRIVATE_KEY, LEN_PUBLIC_KEY};
    use crate::sign;

    pub fn read_mldsa_nist_acvp_kat(path: &str) -> io::Result<Vec<(String, String, String, String, String, String)>> {
        let file = File::open(path).expect("cannot open file");
        let reader = BufReader::new(file);

        let mut results = Vec::new();
        let mut msg = String::new();
        let mut pk = String::new();
        let mut sk = String::new();
        let mut ctx = String::new();
        let mut rnd = String::new();
        let mut sig = String::new();

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();

            if line.is_empty() {
                if !msg.is_empty() {
                    results.push((msg.clone(), pk.clone(), sk.clone(), ctx.clone(), rnd.clone(), sig.clone()));
                    msg.clear();
                    rnd.clear();
                    pk.clear();
                    sk.clear();
                    ctx.clear();
                    sig.clear();
                }
                continue;
            }

            if let Some(val) = line.strip_prefix("msg = ") {
                msg = val.trim().to_string();
            } else if let Some(val) = line.strip_prefix("pkey = ") {
                pk = val.trim().to_string();
            } else if let Some(val) = line.strip_prefix("skey = ") {
                sk = val.trim().to_string();
            } else if let Some(val) = line.strip_prefix("ctx = ") {
                ctx = val.trim().to_string();
            } else if let Some(val) = line.strip_prefix("rnd = ") {
                rnd = val.trim().to_string();
            } else if let Some(val) = line.strip_prefix("sig = ") {
                sig = val.trim().to_string();
            }
        }

        if !msg.is_empty() {
            results.push((msg.clone(), pk.clone(), sk.clone(), ctx.clone(), rnd.clone(), sig.clone()));
        }
        Ok(results)
    }

    fn run_kat_with_ctx(kats: Vec<(String, String, String, String, String, String)>) {
        assert!(!kats.is_empty());
        assert_eq!(kats.len(), 30);
        for (kat_msg, kat_pk, kat_sk, kat_ctx, kat_rnd, kat_sig) in kats.iter() {
            let _pk: [u8; LEN_PUBLIC_KEY] = hex::decode(&kat_pk).unwrap().try_into().unwrap();
            let sk: [u8; LEN_PRIVATE_KEY] = hex::decode(&kat_sk).unwrap().try_into().unwrap();
            let rnd: [u8; 32] = hex::decode(&kat_rnd).unwrap().try_into().unwrap();
            let msg = hex::decode(&kat_msg).unwrap();
            let ctx = hex::decode(&kat_ctx).unwrap();
            let kat_sig = hex::decode(&kat_sig).unwrap();
            let _sig = sign::sign_message_with_ctx(&sk, &rnd, &msg[..], &ctx[..]).unwrap();
            assert_eq!(_sig, kat_sig[..]);
        }
    }

    #[cfg(feature = "ML_DSA_44")]
    #[test]
    fn test_mldsa_44() {
        let kat_list = read_mldsa_nist_acvp_kat("./kats/ml_dsa_44_sign.acvp.kat").unwrap();
        assert!(!kat_list.is_empty());
        run_kat_with_ctx(kat_list);
    }

    #[cfg(feature = "ML_DSA_65")]
    #[test]
    fn test_mldsa_65() {
        let kat_list = read_mldsa_nist_acvp_kat("./kats/ml_dsa_65_sign.acvp.kat").unwrap();
        assert!(!kat_list.is_empty());
        run_kat_with_ctx(kat_list);
    }

    #[cfg(feature = "ML_DSA_87")]
    #[test]
    fn test_mldsa_65() {
        let kat_list = read_mldsa_nist_acvp_kat("./kats/ml_dsa_87_sign.acvp.kat").unwrap();
        assert!(!kat_list.is_empty());
        run_kat_with_ctx(kat_list);
    }
}
