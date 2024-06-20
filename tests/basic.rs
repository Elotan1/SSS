use sss::{reconstruct_bytes, share_bytes, share_leaks, og_matrix,reconstruct_leak};

// Rust verification software (Regression scheme and formal verification)
#[test]
fn singular_secret() {
    let test_hex: &[u8; 1] = b"\x29";
    let t = 6;
    let n = 6;
    let c = share_bytes(test_hex, t, n);
    let d = reconstruct_bytes(c.1, t);
    let a = 0x29;
    assert_eq!(d[0], a);
}

#[test]
fn two_secrets() {
    let test_hex: &[u8; 2] = b"\x29\x49";
    let t = 6;
    let n = 6;
    let c = share_bytes(test_hex, t, n);
    let d = reconstruct_bytes(c.1, t);
    let a = 0x29;
    let b = 0x49;
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
    let a = 0x49 % 97;
    let b = 0x6d % 97;
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
    let a = 0x29;
    let b = 0x49;
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
    let d = 0x29;
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
    let d = 0x29;
    let e = 0x49;
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
    let d = 0x49 % 97;
    let e = 0x6d % 97;
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
    let d = 0x29 % 97;
    let e = 0x49 % 97;
    assert_eq!(c[0], d);
    assert_eq!(c[1], e);
}