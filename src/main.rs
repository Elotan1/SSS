use sss::{gao, share_bytes, reconstruct_bytes, og_matrix, share_leaks};
//sss::{mds, og_matrix, reconstruct_bytes, reconstruct_leak, share_bytes, share_leaks};
use sss::Fq;
use nalgebra::DMatrix;

fn main() {
    let test_hex: &[u8; 2] = b"\x29";
    let t = 6;
    let n = 6;
    let a = og_matrix(t, n);
    let c = share_leaks(test_hex, t, n, a.clone());
    let errors = 2;
    let vec = Vec::from_iter(c.keys());
    let vec2 = Vec::from_iter(c.values());
    println!("VECTOR {:?}", vec);
    println!("VECTOR {:?}", vec2);









    // let zero = Fq::from(0u8);
    // println!("zero in F_257 = {}", zero); // seems to be Display bug with Fp<MontBackend<FqConfig,1>,1>
    // let zero: BigInt<1> = BigInt::from(0u8);
    // println!("{}", zero); // not a bug with BigInt<_>

    // let test_hex: &[u8; 2] = b"\x29\x49";
    // let t = 3;
    // let n = 7;
    // let errors = 2;
    // let vec = DMatrix::from_vec(7, 1, vec![Fq::from(1), Fq::from(5),Fq::from(3), Fq::from(6), Fq::from(3), Fq::from(2), Fq::from(2)]);
    // let c = share_bytes(test_hex, t, n);
    // println!("COFFS {:?}", c.0);
    // println!("C_1{:?}", c.1);
    // let code = gao(n, vec, errors);
    // let d = reconstruct_bytes(c.1, t);
    // println!("{:?}", d);
    // let a = og_matrix(t, n);
    // let p = DMatrix::from_vec(2, 4, vec![Fq::from(1), Fq::from(0), Fq::from(0), Fq::from(1), Fq::from(0), Fq::from(52), Fq::from(0), Fq::from(69)]);
    // println!("{p}");
    // let a = mds(t, n, p);
    // let b = share_leaks(test_hex, t, n, a.clone());
    // reconstruct_leak(b, t, a);
 }
