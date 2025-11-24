
## Automated Cryptographic Validation Test System - Gen/Vals
https://github.com/usnistgov/ACVP-Server

The list of all KATs in JSON:
https://github.com/usnistgov/ACVP-Server/tree/master/gen-val/json-files

We import the vectors in this file to test our implementation:
https://github.com/usnistgov/ACVP-Server/blob/master/gen-val/json-files/ML-DSA-keyGen-FIPS204/internalProjection.json

For more information, see: https://pages.nist.gov/ACVP/draft-celi-acvp-ml-dsa.html

## Steps to Prepare KAT File
```bash
$ pwd
/Users/dprasad/projects/dicp-local/ml-dsa-44/kats

$ python3 nist-fips204-kats.py 
Written output to nist-acvp-keygen-kats.txt

$ 

```

### MAYBE other "raw" inputs

1. Dilithium Reference Implementation KATs
URL: https://github.com/pq-crystals/dilithium/tree/master/ref
The original CRYSTALS-Dilithium reference implementation has KATs, though they're for Dilithium 3.1, not the final FIPS 204 (ML-DSA). There are some minor differences (see FIPS 204 Appendix D).

2. OpenSSL/BoringSSL Test Vectors
BoringSSL: https://github.com/google/boringssl/tree/master/crypto/fipsmodule/mldsa
BoringSSL has implemented ML-DSA and includes test vectors in their source.

3. LibOQS Test Vectors
URL: https://github.com/open-quantum-safe/liboqs/tree/main/tests/kat
LibOQS (Open Quantum Safe) has KAT files for ML-DSA with various test cases.

4. Generate Your Own KATs
If you want to generate test vectors with specific context strings, you can:

Use a reference implementation like:

https://github.com/pq-crystals/dilithium (original reference)
https://github.com/C2SP/CCTV/tree/main/ML-DSA (C2SP test vectors)

Use the NIST reference code (when available)

Recommendation
For comprehensive testing with context strings, I'd recommend:

Start with NIST ACVP vectors - They're official and comprehensive
Use the post-quantum-cryptography/KAT for basic "raw" mode testing (what you're doing now)
Generate your own test cases with various context lengths to ensure your implementation handles:

Empty context (ctx = "")
Short context (1-10 bytes)
Maximum context (255 bytes)
Edge cases


```