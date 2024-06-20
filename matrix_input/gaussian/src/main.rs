// This is the list of imports that my Gauss-Jordan Reduction Algorithm depends on.
// Please note that additional dependencies are also implemented in Cargo.toml.
extern crate nalgebra as na;
mod func;
use func::share_deconstruction;
use func::share_reconstruction;
use func::share_bytes;
use func::reconstruct_bytes;
use std::collections::HashMap;
use func::Fq;
use func::FqConfig;
fn main() {
    let test_hex: &[u8; 3] = b"\x12\x34\x12";
    let t = 6;
    let n = 6;
    let c = share_bytes(test_hex, t, n);
    println!("{:?}", c);
    let d = reconstruct_bytes(c, t);
    println!("{:?}", d);
    let a = Fq::from(0x12);
    println!("{:?}", a);
    let a = Fq::from(0x34);
    println!("{:?}", a);





    // let map = HashMap::from([
    //     (Fq::from(6), vec![
    //         Fq::from(1),
    //         Fq::from(2),
    //     ],),
    //     (Fq::from(1), vec![
    //         Fq::from(1),
    //         Fq::from(2),
    //     ],),
    //     (Fq::from(2), vec![
    //         Fq::from(2),
    //         Fq::from(2),
    //     ],),
    //     (Fq::from(3), vec![
    //         Fq::from(1),
    //         Fq::from(2),
    //     ],),
    //     (Fq::from(4), vec![
    //         Fq::from(1),
    //         Fq::from(2),
    //     ],),
    //     (Fq::from(5), vec![
    //         Fq::from(1),
    //         Fq::from(2),
    //     ],),
    //     (Fq::from(6), vec![
    //         Fq::from(1),
    //         Fq::from(2),
    //     ],),
    //     (Fq::from(7), vec![
    //         Fq::from(1),
    //         Fq::from(2),
    //     ],),
    //     (Fq::from(8), vec![
    //         Fq::from(1),
    //         Fq::from(2),
    //     ],),
    //     (Fq::from(9), vec![
    //         Fq::from(1),
    //         Fq::from(2),
    //     ],),
    //     (Fq::from(10), vec![
    //         Fq::from(1),
    //         Fq::from(2),
    //     ],),

    // ]);
    // println!("{:?}", map);
    // let test_hex: &[u8; 2] = b"\x29\x27";
    // let t = 2;
    // let n = 2;
    // let c = share_bytes(test_hex, t, n);
    // share_deconstruction(Fq::from(2), 2, 2);
    // share_reconstruction(map.clone(), 2);
}
