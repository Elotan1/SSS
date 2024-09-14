use sss::{deconstruct_leaks, gao, og_matrix, reconstruct_bytes, reconstruct_leak, share_bytes, share_leaks};
use nalgebra::DMatrix;
use sss::Fq;
// Rust verification software (Regression scheme and formal verification)
#[test]
fn singular_secret() {
    let test_hex: &[u8; 1] = b"\x29";
    let t = 6;
    let n = 6;
    let c = share_bytes(test_hex, t, n);
    let d = reconstruct_bytes(c.1, t);
    let a = 0x29 % 7;
    assert_eq!(d[0], a);
}

#[test]
fn two_secrets() {
    let test_hex: &[u8; 2] = b"\x29\x49";
    let t = 6;
    let n = 6;
    let c = share_bytes(test_hex, t, n);
    let d = reconstruct_bytes(c.1, t);
    let a = 0x29 % 7;
    let b = 0x49 % 7;
    assert_eq!(d[0], a);
    assert_eq!(d[1], b);
}

#[test]
fn big_secret() {
    let test_hex: &[u8; 48] = b"\x49\x27\x6d\x20\x6b\x69\x6c\x6c\x69\x6e\x67\x20\x79\x6f\x75\x72\x20\x62\x72\x61\x69\x6e\x20\x6c\x69\x6b\x65\x20\x61\x20\x70\x6f\x69\x73\x6f\x6e\x6f\x75\x73\x20\x6d\x75\x73\x68\x72\x6f\x6f\x6d";
    let t = 6;
    let n = 6;
    let c = share_bytes(test_hex, t, n);
    let d = reconstruct_bytes(c.1, t);
    let a = 0x49 % 7;
    let b = 0x6d % 7;
    assert_eq!(d[0], a);
    assert_eq!(d[47], b);
}

#[test]
fn limited_t() {
    let test_hex: &[u8; 2] = b"\x29\x49";
    let t = 1;
    let n = 6;
    let c = share_bytes(test_hex, t, n);
    let d = reconstruct_bytes(c.1, t);
    let a = 0x29 % 7;
    let b = 0x49 % 7;
    assert_eq!(d[0], a);
    assert_eq!(d[1], b);
}

#[test]
fn lrss_one_secret_test() {
    let test_hex: &[u8; 1] = b"\x29";
    let t = 6;
    let n = 6;
    let a = og_matrix(t, n);
    let b = share_leaks(test_hex, t, n, a.clone());
    let c = reconstruct_leak(b, t, a);
    let d = 0x29 % 7;
    assert_eq!(c[0], d);
}


#[test]
fn lrss_two_secrets_test() {
    let test_hex: &[u8; 2] = b"\x29\x49";
    let t = 6;
    let n = 6;
    let a = og_matrix(t, n);
    let b = share_leaks(test_hex, t, n, a.clone());
    let c = reconstruct_leak(b, t, a);
    let d = 0x29 % 7;
    let e = 0x49 % 7;
    assert_eq!(c[0], d);
    assert_eq!(c[1], e);
}

#[test]
fn big_secret_lrss() {
    let test_hex: &[u8; 48] = b"\x49\x27\x6d\x20\x6b\x69\x6c\x6c\x69\x6e\x67\x20\x79\x6f\x75\x72\x20\x62\x72\x61\x69\x6e\x20\x6c\x69\x6b\x65\x20\x61\x20\x70\x6f\x69\x73\x6f\x6e\x6f\x75\x73\x20\x6d\x75\x73\x68\x72\x6f\x6f\x6d";
    let t = 6;
    let n = 6;
    let a = og_matrix(t, n);
    let b = share_leaks(test_hex, t, n, a.clone());
    let c = reconstruct_leak(b, t, a);
    let d = 0x49 % 7;
    let e = 0x6d % 7;
    assert_eq!(c[0], d);
    assert_eq!(c[47], e);
}

#[test]
fn limited_t_lrss() {
    let test_hex: &[u8; 2] = b"\x29\x49";
    let t = 1;
    let n = 2;
    let a = og_matrix(t, n);
    let b = share_leaks(test_hex, t, n, a.clone());
    let c = reconstruct_leak(b, t, a);
    let d = 0x29 % 7;
    let e = 0x49 % 7;
    assert_eq!(c[0], d);
    assert_eq!(c[1], e);
}

#[test]
fn berlekamp() {
    let n = 7;
    let errors = 2;
    let vec = DMatrix::from_vec(7, 1, vec![Fq::from(1), Fq::from(5),Fq::from(3), Fq::from(6), Fq::from(3), Fq::from(2), Fq::from(2)]);
    let expected = gao(n, vec, errors);
    let solution = DMatrix::from_vec(7, 1, vec![Fq::from(1), Fq::from(6),Fq::from(3), Fq::from(6), Fq::from(1), Fq::from(2), Fq::from(2)]);
    assert_eq!(expected, solution);
}

#[test]
fn berlekamp_2() {
    let n = 7;
    let errors = 2;
    let vec = DMatrix::from_vec(7, 1, vec![Fq::from(1), Fq::from(6), Fq::from(123), Fq::from(456), Fq::from(57), Fq::from(86), Fq::from(121)]);
    let expected = gao(n, vec, errors);
    let solution = DMatrix::from_vec(7, 1, vec![Fq::from(1), Fq::from(6),Fq::from(3), Fq::from(6), Fq::from(1), Fq::from(2), Fq::from(2)]);
    assert_eq!(expected, solution);
}

#[test]
fn berlekamp_zero_errors() {
    let n = 7;
    let errors = 0;
    let vec = DMatrix::from_vec(7, 1, vec![Fq::from(1), Fq::from(6),Fq::from(3), Fq::from(6), Fq::from(1), Fq::from(2), Fq::from(2)]);
    let expected = gao(n, vec.clone(), errors);
    assert_eq!(expected, vec);
}

#[test]
#[should_panic]
fn berlekamp_non_zero_remainder() {
    let n = 7;
    let errors = 2;
    let vec = DMatrix::from_vec(7, 1, vec![Fq::from(1), Fq::from(2), Fq::from(2), Fq::from(3), Fq::from(1), Fq::from(2), Fq::from(5)]);
    gao(n, vec.clone(), errors);
}

#[test]
fn berlekamp_simple_equation() {
    let n = 7;
    let errors = 2;
    let vec = DMatrix::from_vec(7, 1, vec![Fq::from(4), Fq::from(1), Fq::from(4), Fq::from(2), Fq::from(2), Fq::from(5), Fq::from(1)]);
    let expected = gao(n,vec,errors);
    let solution = DMatrix::from_vec(7, 1, vec![Fq::from(0), Fq::from(1), Fq::from(4), Fq::from(2), Fq::from(2), Fq::from(4), Fq::from(1)]);
    assert_eq!(expected, solution);
}

#[test]
fn berlekamp_with_shamirs() {
    let test_hex: &[u8; 2] = b"\x29\x49";
    let t = 6;
    let n = 6;
    let a = og_matrix(t, n);
    let c = share_leaks(test_hex, t, n, a.clone());
    let errors = 2;
    let vec = Vec::from_iter(c.keys());
    println!("VECTOR {:?}", vec);
    // let expected = gao(n,vec,errors);
    // let solution = DMatrix::from_vec(7, 1, vec![Fq::from(0), Fq::from(1), Fq::from(4), Fq::from(2), Fq::from(2), Fq::from(4), Fq::from(1)]);
    // assert_eq!(expected, solution);
    let d = reconstruct_leak(c, t, a);
    let f = 0x29 % 7;
    let e = 0x49 % 7;
    assert_eq!(d[0], f);
    assert_eq!(d[1], e);
}