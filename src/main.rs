use sss::{gao, mds, og_matrix, reconstruct_bytes, reconstruct_leak, share_bytes, share_leaks};
use sss::Fq;
use nalgebra::DMatrix;

/// Right now we can use main as a sandbox.
/// Later on may convert to command line utility
fn main() {
    // let zero = Fq::from(0u8);
    // println!("zero in F_257 = {}", zero); // seems to be Display bug with Fp<MontBackend<FqConfig,1>,1>
    // let zero: BigInt<1> = BigInt::from(0u8);
    // println!("{}", zero); // not a bug with BigInt<_>

    //let test_hex: &[u8; 48] = b"\x49\x27\x6d\x20\x6b\x69\x6c\x6c\x69\x6e\x67\x20\x79\x6f\x75\x72\x20\x62\x72\x61\x69\x6e\x20\x6c\x69\x6b\x65\x20\x61\x20\x70\x6f\x69\x73\x6f\x6e\x6f\x75\x73\x20\x6d\x75\x73\x68\x72\x6f\x6f\x6d";
    let test_hex: &[u8; 1] = b"\x2A";
    let t = 3;
    let n = 7;
    let errors = 2;
    let vec = DMatrix::from_vec(7, 1, vec![Fq::from(1), Fq::from(5),Fq::from(3), Fq::from(6), Fq::from(3), Fq::from(2), Fq::from(2)]);
    // let c = share_bytes(test_hex, t, n);
    // println!("COFFS {:?}", c.0);
    // println!("{:?}", c.1);
    // let d = reconstruct_bytes(c.1, t);
    // println!("{:?}", d);
    gao(n, vec, errors);
    // let a = og_matrix(t, n);
    // let p = DMatrix::from_vec(2, 4, vec![Fq::from(1), Fq::from(0), Fq::from(0), Fq::from(1), Fq::from(0), Fq::from(52), Fq::from(0), Fq::from(69)]);
    // println!("{p}");
    // let a = mds(t, n, p);
    // let b = share_leaks(test_hex, t, n, a.clone());
    // reconstruct_leak(b, t, a);
 }
