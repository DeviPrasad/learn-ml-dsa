mod keypair;
mod err;
mod params;
mod ntt;
mod sign;

fn main() {
    let (_, sk) = keypair::key_gen().unwrap();
    let _ = sign::sk_decode(&sk).unwrap();
}

#[cfg(test)]
mod tests {
    use crate::{keypair, sign};

    #[test]
    fn dry_512_keygen_test() {
        for _ in 0..512 {
            let (_, sk) = keypair::key_gen().unwrap();
            let _ = sign::sk_decode(&sk).unwrap();
        }
    }
}
