use feanor_math::homomorphism::*;
use feanor_math::assert_el_eq;
use feanor_math::ring::*;
use feanor_math::primitive_int::*;
use feanor_math::rings::zn::zn_big::*;
use feanor_math::rings::zn::*;
use feanor_math::rings::finite::*;
use feanor_math::algorithms;
use oorandom;

fn fermat_is_prime(n: i64) -> bool {
    // the Fermat primality test is based on the observation that a^n == a mod n if n
    // is a prime; On the other hand, if n is not prime, we hope that there are many
    // a such that this is not the case. 
    // Note that this is not always the case, and so more advanced primality tests should 
    // be used in practice. This is just a proof of concept.

    let ZZ = StaticRing::<i64>::RING;
    let Zn = Zn::new(ZZ, n); // the ring Z/nZ

    // check for 6 random a whether a^n == a mod n
    let mut rng = oorandom::Rand64::new(0);
    for _ in 0..6 {
        let a = Zn.random_element(|| rng.rand_u64());
        let a_n = Zn.pow(Zn.clone_el(&a), n as usize);
        if !Zn.eq_el(&a, &a_n) {
            return false;
        }
    }
    return true;
}

assert!(algorithms::miller_rabin::is_prime(StaticRing::<i64>::RING, &91, 6) == fermat_is_prime(91));