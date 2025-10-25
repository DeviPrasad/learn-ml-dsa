mod keypair;
mod err;
mod params;
mod ntt;

fn main() {
    let _ = keypair::key_gen().unwrap();
}

#[cfg(test)]
mod tests {
    use crate::keypair;

    #[test]
    fn dry_512_keygen_test() {
        for _ in 0..512 {
            let _ = keypair::key_gen().unwrap();
        }
    }
}
