// use anyhow::{anyhow, Error, Ok, Result};
use ark_ff::fields::{Fp64, MontBackend, MontConfig};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_ff::UniformRand;
use func::find_inverse;
use func::gaussian;
use nalgebra::iter::ColumnIter;
use nalgebra::DMatrix;
use std::collections::HashMap;
use std::vec;
mod func;
use ark_ff::BigInteger;
use ark_ff::PrimeField;
use func::share_deconstruction;
use func::share_reconstruction;
use func::leakage_deconstruction;
use func::leaks_reconstruction;
use itertools::Itertools;

// TODO: thiserror crate for error handling?
// - Seems more responsible than anyhow for something like cryptographic implementation
// - anyhow is just so easy though

// Field we have currently decided to use. May change later.
// Correctness should not depend on which field we choose here.
#[derive(MontConfig)]
#[modulus = "97"]
#[generator = "12"]
pub struct FqConfig;
pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct Share {
    pub x: Fq,
    pub y: Fq,
}

#[derive(Clone, Debug)]
pub struct LRShare {
    pub x: usize,
    pub y: Fq,
}

/// Generates n Shamir's secret shares of the form (x, p(x)),
/// where p is a degree t-1 polynomial with constant coefficient
/// msg, and all other coefficients chosen uniformly at random
/// from F_q (field of prime size q).
/// By convention, the set of evaluation points can be
/// chosen as 1,2,...,n (mod q). Note this requires that n is
/// at most q-1.
pub fn share_field_element(msg: Fq, t: usize, n: usize) -> (Vec<Fq>, Vec<Share>) {
    // assert!(BigInt::from(n as u64) < <Fq as PrimeField>::MODULUS);
    // (1..=n)
    //     .map(|i| Share { x: Fq::from(i as u64), y: Fq::from(i as u64) } )
    //     .collect()
    // vec![(Fq::from(0u8), Fq::from(1u8)); n.into()]
    let (base, current) = share_deconstruction(msg, t, n);
    let mut vector_form = vec![];
    let mut count = 0;
    let mut val = 1;
    while count < n {
        let share = Share {
            x: Fq::from(val),
            y: current[(count, 1)],
        };
        vector_form.insert(count, share);
        count += 1;
        val += 1;
    }
    let mut coefficient = vec![];
    let mut count = 0;
    let mut val = 1;
    while count < base.len() {
        coefficient.insert(count, base[count]);
        count += 1;
    }
    (coefficient, vector_form)
}

/// Performs the reconstruction algorithm for Shamir's secret sharing
/// assuming a threshold t. Given a list of shares of the form (x, p(x)),
/// interpolate the degree t-1 polynomial p and return its constant
/// coefficient. Note this requires that the length of the shares vector
/// is at least t (truncate the list arbitrarily if it is longer than t).
pub fn reconstruct_bytes(shares: HashMap<Fq, Vec<Fq>>, t: usize) -> Vec<u8> {
    let current = share_reconstruction(shares, t);
    let mut new_form = vec![];
    let mut count = 0;
    while count < current.len() {
        let value = current[(0, count)];
        new_form.insert(count, value);
        count += 1;
    }
    new_form
        .into_iter()
        .map(|field_element| PrimeField::into_bigint(field_element).to_bytes_le()[0])
        .collect()
}

/// Takes an arbitrary length message in bytes, independently (t,n) secret shares
/// each byte, then groups the shares so that all shares with the same x coordinate
/// go to the same party.
pub fn share_bytes(msg: &[u8], t: usize, n: usize) -> (Vec<Vec<Fq>>,HashMap<Fq, Vec<Fq>>) {
    
    // secret share each byte,
    // then bin shares into the binned_shares HashMap by evaluation point (x coordinate)
    let mut binned_shares: HashMap<Fq, Vec<Fq>> = HashMap::new();
    let mut coefficients: Vec<Vec<Fq>> = vec![vec![]];
    msg.iter().for_each(|byte| {
            let (coeffs, shares) =  share_field_element(Fq::from(*byte), t, n);
            shares.into_iter().for_each(|s| {
                // if the x coordinate of s already a key in binned_shares, then just push the y value of s to binned_shares[x]
                if let Some(share_bin) = binned_shares.get_mut(&s.x) {
                    share_bin.push(s.y);
                } else {
                    binned_shares.insert(s.x, vec![s.y]);
                }
            });
            coefficients.push(coeffs);
        });
        coefficients.remove(0);
    // TODO: how to serialize HashMap<Fq, Vec<Fq>>?
    // serializing a hashmap definitely seems more involved than a struct
    // seems like we can theoretically implement our own serialization and
    // deserialization, but could be tricky. See https://github.com/arkworks-rs/algebra/tree/master/serialize
    // Not even sure how to serialize a Vec<Share>, so not like there's an easy inefficient cop-out
    // Maybe worth trying serde and avoiding arkworks built-in serialization.
    // Goes against some of the advice found in this discussion: https://github.com/arkworks-rs/algebra/issues/178
    // i.e. must be careful serializing cryptographic objects.
    // Regardless, serde could be worth looking into for this.
    // For now I will just ignore serialization and return HashMap<Fq, Vec<Fq>>
    (coefficients, binned_shares)
}

pub fn deconstruct_leaks(msg: Fq, t: usize, n: usize, original: DMatrix<Fq>) -> Vec<LRShare> {
    let current = leakage_deconstruction(msg, t, original);
    let mut vector_form = vec![];
    let mut count = 1;
    while count < n + 1{
        let share = LRShare {
            x: count,
            y: current[(0, count)],
        };
        vector_form.push(share);
        count += 1;
    }
    println!("Shares: {:?}", vector_form);
    vector_form
}

pub fn reconstruct_leak(shares: HashMap<usize, Vec<Fq>>, t: usize, original: DMatrix<Fq>) -> Vec<u8> {
        let current = leaks_reconstruction(shares, t, original);
        current
        .into_iter()
        .map(|field_element| PrimeField::into_bigint(field_element).to_bytes_le()[0])
        .collect()
}

pub fn share_leaks(msg: &[u8], t: usize, n: usize, original: DMatrix<Fq>) -> HashMap<usize, Vec<Fq>> {
    // secret share each byte,
    // then bin shares into the binned_shares HashMap by evaluation point (x coordinate)
    let mut leak_shares: HashMap<usize, Vec<Fq>> = HashMap::new();
    let shares: Vec<Vec<LRShare>> = msg
        .iter()
        .map(|byte| deconstruct_leaks(Fq::from(*byte), t, n, original.clone()))
        .collect();
    let flat_shares: Vec<LRShare> = shares.into_iter().flatten().collect();
    flat_shares.into_iter().for_each(|s| {
        // if the x coordinate of s already a key in binned_shares, then just push the y value of s to binned_shares[x]
        if let Some(share_bin) = leak_shares.get_mut(&s.x) {
            share_bin.push(s.y);
        } else {
            leak_shares.insert(s.x, vec![s.y]);
        }
    });
    leak_shares
}

pub fn og_matrix(threshold: usize, n: usize) -> DMatrix<Fq> {
    let k = threshold - 1;
    let random = n - k;
    let c = DMatrix::from_vec(0,0,vec![]);
        // Resizes the matrix to create the identity augmentation.
        let mut c = c.resize(threshold, random + threshold, Fq::from(0));
        let mut column = 0;
        let mut row = 0;
        while column < (threshold) {
            c[(row, column)] = Fq::from(1);
            row += 1;
            column += 1;
        }
        println!("Augmented with the identity: {}", c);
        let mut rng = ark_std::test_rng();
        let mut column = 0;
    
    // Let's sample uniformly random field elements with the original message as the first element in the vector:
    
    while column < random {
        let mut base = DMatrix::from_vec(1, 1, vec![Fq::from(5)]);
        let mut count = 1;
    while count < threshold {
        let a = Fq::rand(&mut rng);
        base = base.insert_row(count, a);
        count += 1;
    }
    let a = Fq::rand(&mut rng);
    base[0] = a;
    let mut i = 0;
    while i < threshold {
        c[(i, threshold + column)] = base[i];
        i += 1;
    }
    column += 1;
    }
    println!("{c}");
    mds(threshold, n, c)
}
    

pub fn mds(threshold: usize, n: usize, mut c: DMatrix<Fq>) -> DMatrix<Fq>{
    let d = c.clone_owned();
    let d = d.remove_column(0);
    let a = DMatrix::from_vec(1, 1, vec![Fq::from(0)]);
    let mut a = a.resize(threshold, threshold, Fq::from(0));
    let mut column = 0;
    let mut row = 0;
    while column < threshold {
        a[(row, column)] = Fq::from(1);
        row += 1;
        column += 1;
    }

    let temp: ColumnIter<_, _, _, _> = d.column_iter();
    let temp2: itertools::Combinations<ColumnIter<_, _, _, _>> = temp.combinations(threshold);
    temp2.for_each(|val:Vec<nalgebra::Matrix<_, _, nalgebra::Const<1>, nalgebra::ViewStorage<_, _, nalgebra::Const<1>, _, _>>> | {
        let new_form = DMatrix::from_vec(1, 1, vec![Fq::from(1)]);
        let mut new_form = new_form.resize(threshold, threshold, Fq::from(1));
        let mut count = 0;
        val.clone().into_iter().for_each(|v: nalgebra::Matrix<_, _, nalgebra::Const<1>, nalgebra::ViewStorage<_, _, nalgebra::Const<1>, _, _>>|{
            new_form.set_column(count, &v);
            count += 1;
        });
        let p = gaussian(new_form, threshold);
        let z = find_inverse(p.clone(), threshold);
        let x = p * z;
        if x != a {
            let v: DMatrix<Fq> = og_matrix(threshold, n);
            println!("The new matrix is: {v}");
            c = v;
        }
    });
    println!("The matrix is: {c}");
    c
}