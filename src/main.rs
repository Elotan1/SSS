use std::os::windows::thread;

use rand::{thread_rng, Rng};
use sss::Fq;
use ark_ff::UniformRand;
use sss::{commitment, reconstruct_bytes, share_bytes};
//sss::{mds, og_matrix, reconstruct_bytes, reconstruct_leak, share_bytes, share_leaks};

fn main() {
    let test_hex: &[u8; 2] = b"\x29\x04";
    let mut rng = thread_rng();
    let r_msg = u8::rand(&mut rng);
    println!("{:?}", r_msg);
    let t = 6;
    let n = 6;
    let c = commitment(test_hex, &[r_msg], t, n);
    println!("{:?}", c);
    // let d = reconstruct_bytes(c.1, t);
    // // let e = 41 % 7;
    // assert_eq!(d[0], r_msg%7);
 }
