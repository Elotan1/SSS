// This is the list of imports that my Gauss-Jordan Reduction Algorithm depends on.
// Please note that additional dependencies are also implemented in Cargo.toml.
extern crate nalgebra as na;
use crate::na::matrix;
use ark_ff::fields::{Field, Fp64, MontBackend, MontConfig};
use ark_ff::Zero;
use na::DMatrix;
use std::io;
use std::collections::HashMap;
use ark_std::UniformRand;
use ark_ff::{fields, BigInt, PrimeField};
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize};

// Field we have currently decided to use. May change later.
// Correctness should not depend on which field we choose here.
#[derive(MontConfig)]
#[modulus = "257"]
#[generator = "12"]
pub struct FqConfig;
pub type Fq = Fp64<MontBackend<FqConfig, 1>>;


#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct Share {
    pub x: Fq,
    pub y: Fq
}


// This is the importable function for the Gauss-Jordan Algorithm.
// {Warning} Not all bugs have been accounted for as of 06/04/2024.
// {Parameters} The function takes in a Mutable Dynamic Matrix that has finite field elements.
// {Returns} The function returns a Mutable Dynamic Matrix that has been reduced to RREF.
pub fn gaussian(mut c: DMatrix<Fq>, size: usize) -> DMatrix<Fq> {
    let row_index: usize = 0; // The first row
    let col_index: usize = 0; // The first column
    let mut col_shift = 0; // The number of rows from the first column
    let size = size; // The size of the matrix (Still being implemented)

    // Iterates through the entire matrix to ensure that no zero element is the pivot point of the column.
    // If a zero pivot is identified, the matrix will perform a row swap with a non-zero pivot row.
    // If all column elements are zero, the matrix will only reduce succeeding columns.
    // Then iterates through the entire matrix, performing row reduction to create the REF of the matrix.
    // Takes advantage of modular arithmetic to reduce the matrix into upper triangular form.
    while col_shift < size - 1 {
        let mut row_shift = 0;
        let mut count = 0;
        while row_shift < size - 1 - col_shift {
            let mut first = c[(row_index + col_shift + count, col_index + col_shift)];
            if first == Fq::from(0) {
                let mut current_row = count;
                while first == Fq::from(0) {
                    if current_row != size - 1 - col_shift {
                        current_row += 1;
                        first = c[(row_index + col_shift + current_row, col_index + col_shift)];
                        let borrow = c.clone();
                        c.set_row(
                            row_index + col_shift + count,
                            &borrow.row(row_index + col_shift + current_row),
                        );
                        c.set_row(
                            row_index + col_shift + current_row,
                            &borrow.row(row_index + col_shift + count),
                        );
                    } else {
                        break;
                    }
                }
            }
            let first = c[(row_index + col_shift, col_index + col_shift)];
            let inverse = first.inverse().unwrap_or(Fq::zero());
            let inverse = inverse * Fq::from(c[(size - row_shift - 1, col_index + col_shift)]);
            let co1 = matrix![inverse];
            let product = &(co1 * c.row(row_index + col_shift));
            c.set_row(
                size - 1 - row_shift,
                &(product - c.row(size - 1 - row_shift)),
            );
            row_shift += 1;
            count += 1;
        }
        col_shift += 1;
    }
    println!("Row Echelon Form of Matrix: {}", c);

    let mut col_reducer = 0; // The number of columns from the last column in the matrix

    // Iterates through the entire matrix, reducing the upper triangular elements to zero to create the RREF of the matrix.
    // This should create the identity matrix if the matrix is invertible.
    while col_reducer < size - 1 {
        let mut row_reducer = size - 1 - col_reducer; // The number of rows from the last row in the matrix
        let col_index = size - 1 - col_reducer;
        let mut row_index = size - 1 - col_reducer;
        let pivot = size - 1 - col_reducer;
        while row_reducer > 0 {
            if c[(pivot, pivot)] == Fq::from(0) {
                break;
            } else {
                let mut reducer = size;
                while reducer > 0 {
                    let current = c[(reducer - 1, reducer - 1)];
                    let inverse = current.inverse().unwrap_or(Fq::zero());
                    let co2 = matrix![Fq::from(inverse)];
                    let product = &(co2 * c.row(reducer - 1));
                    c.set_row(reducer - 1, &(product));
                    reducer -= 1;
                }
                let current = c[(row_index - 1, col_index)];
                let co2 = matrix![Fq::from(current)];
                let product = &(co2 * c.row(pivot));
                c.set_row(row_index - 1, &(product - c.row(row_index - 1)));
                row_index -= 1;
                row_reducer -= 1;
            }
        }
        col_reducer += 1;
    }
    let current = c[(0, 0)];
    let inverse = current.inverse().unwrap_or(Fq::zero());
    let co2 = matrix![Fq::from(inverse)];
    let product = &(co2 * c.row(0));
    c.set_row(0, &(product));
    println!("Reduced Row Echelon Form of Matrix: {}", c);
    c
}

// This is the importable function for finding the inverse of an inveritble matrix using gaussian reduction.
// Augments the original matrix with the identity, places it into RREF form, and isolates the changed identity.
// {Warning} Not all bugs have been accounted for as of 06/04/2024.
// {Parameters} The function takes in a Mutable Dynamic Matrix that has finite field elements.
// {Returns} The function returns a Mutable Dynamic Matrix that has been inverted.
pub fn find_inverse(c: DMatrix<Fq>, size: usize) -> DMatrix<Fq> {
    let size = size;
    // Creates a copy of the original matrix for future reference.
    let b = c.clone();

    // Resizes the matrix to create the identity augmentation.
    let mut c = c.resize(size, size * 2, Fq::from(0));
    let mut column = size;
    let mut row = 0;
    while column < size * 2 {
        c[(row, column)] = Fq::from(1);

        row += 1;
        column += 1;
    }
    println!("Augmented with the identity: {}", c);

    // Performs Gauss-Jordan reduction on the matrix.
    let d = gaussian(c.clone(), size);

    // Isolates the reduced identity matrix to return the inverse.
    let c = d.columns(size, size);
    println!("Inverse of Matrix: {}", c);

    // Verifies that the inverse is correct by multiplying matrix A by A^-1 to get I.
    let identity = b * c;
    println!(
        "If this is the identity, then the inverse is correct: {}",
        identity
    );
    c.into()
}

// This utilizes arkworks and nalgebra to create a matrix with finite field elements.
// Accepts user inputs to create a square matrix with specified inputs.
// Uses modular arithmetic to reduce any large elements by the field.
pub fn create_matrix() -> (DMatrix<Fq>, usize) {
    println!("Please input the number of rows and columns.");

    let mut row_number = String::new();

    io::stdin()
        .read_line(&mut row_number)
        .expect("Failed to read line");

    let row_number: u32 = match row_number.trim().parse() {
        Ok(num) => num,
        Err(_) => todo!(),
    };

    println!("You requested a {row_number}x{row_number} matrix.");

    let mut row_index: usize = 0;
    let size: usize = row_number.try_into().unwrap();
    let a = matrix![Fq::from(0)];

    let mut c = a.resize(size, size, Fq::from(0));

    println!("{}", a.resize(size, size, Fq::from(0)));

    while row_index < size {
        let mut col_index = 0;
        while col_index < size {
            println!("Please input the element for a[{row_index}, {col_index}]");

            let mut index_value = String::new();

            io::stdin()
                .read_line(&mut index_value)
                .expect("Failed to read line");

            let value: u32 = match index_value.trim().parse::<u32>() {
                Ok(num) => num,
                Err(_) => continue,
            };

            let value = Fq::from(value);

            c[(row_index, col_index)] = value;
            col_index += 1;
        }
        row_index += 1;
    }
    println!("Creating a {row_number}x{row_number} matrix:");
    println!("{}", c);
    (c, size)
}


// pub fn share_bytes(msg: &[u8], t: usize, n: usize) -> HashMap<Fq, Vec<Fq>> {
//     // secret share each byte,
//     // then bin shares into the binned_shares HashMap by evaluation point (x coordinate)
//     let mut binned_shares: HashMap<Fq, Vec<Fq>> = HashMap::new();
//     msg.iter()
//         .map(|byte| share_field_element(Fq::from(*byte), t, n))
//         .flatten()
//         .for_each(|s| {
//             if let Some(share_bin) = binned_shares.get_mut(&s.x) {
//                 share_bin.push(s.y);
//             } else {
//                 binned_shares.insert(s.x, vec![s.y]);
//             }
//         });
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
//     binned_shares
    
// }
pub fn vandermonde(binned_shares: HashMap<Fq, Vec<Fq>>, t: usize) -> DMatrix<Fq> {
    let mut bucket_vector = DMatrix::from_vec(0, 0, vec![]);
    let mut n = 1;
    for key in binned_shares.clone().into_keys() {
        assert_ne!(key, Fq::from(0));
        let bucket_resize = bucket_vector.resize(n, 1, key);
        bucket_vector = bucket_resize.clone();
        n += 1;
    }
    let copy = bucket_vector.rows(0,t);
    println!("Using these {t} keys: {}", copy);
    let mut scope = copy.clone().insert_columns(0, 1, Fq::from(1));
    let mut d = 2;
    while d < t {
        scope = scope.insert_column(d, Fq::from(1));
    d += 1;
}
    let size = t;
    let row_index = 0;
    let col_index = 0;
    let mut col_shift = 2;
    
    while col_shift < size {
        let mut row_shift = 0;
        while row_shift < size {
            let current = scope[(row_shift + row_index, 1)];
            let mut exponent = 1;
            let mut product = current;
            while exponent < col_shift {
                    product *= current; 
                    scope[(row_shift + row_index, col_shift + col_index)] = product;
                exponent += 1;
            }
            row_shift += 1;
        }
        col_shift += 1;
    }
    let new = scope.columns(0, t);
    println!("The vandermonde matrix is: {new}");
    new.into()
}

pub fn value_vector(binned_shares: HashMap<Fq, Vec<Fq>>, t: usize) -> DMatrix<Fq>{
    let mut big: Vec<Vec<Fq>> = vec![];
    let mut keys: Vec<Fq> = vec![];
    let mut n = 0;
    binned_shares.into_iter().for_each(|(key, val)| {
        if n < t {
        keys.push(key);
        big.push(val);
        }
        n += 1;
    });
    let m = big.clone()[0].len();
    for v in big.clone() {
        assert_eq!(v.len(), m);
    }
    println!("{:?}", keys);
    println!("{:?}", big);
    println!("{:?}", big.clone().into_iter());
    let c = DMatrix::from_vec(m, t, big.into_iter().flatten().collect());
    let d = c.transpose();
    d
}

pub fn reconstruction(a: DMatrix<Fq>, b: DMatrix<Fq>) -> DMatrix<Fq>{
    let c = a * b;
    println!("{c}");
    c
}

pub fn share_reconstruction(binned_shares: HashMap<Fq, Vec<Fq>>, t: usize) -> DMatrix<Fq> {
    let holder = vandermonde(binned_shares.clone(), t);
    let inverse = find_inverse(holder, t);
    let values = value_vector(binned_shares.clone(), t);
    let result = reconstruction(inverse, values);
    let secret = result.rows(0, 1);
    println!("The secret is: {secret}");
    secret.into()
}

pub fn share_deconstruction(msg: Fq, t: usize, n: usize) -> DMatrix<Fq>{
    let mut rng = ark_std::test_rng();
    let mut base = DMatrix::from_vec(1, 1, vec![msg]);
    let mut count = 1;
    // Let's sample uniformly random field elements:
    while count < t {
        let a = Fq::rand(&mut rng);
        base = base.insert_row(count, a);
        count += 1;
    }
    println!("The message and random coefficients are: {base}");
    assert!(n >= t);
    let mut count = 1;
    let mut val = 2;
    let mut new = DMatrix::from_vec(1, 1, vec![Fq::from(1)]);
    while count < n {
        new = new.insert_row(count, Fq::from(val));
        count += 1;
        val += 1;
    }
    println!("The number of shares are: {new}");
    let mut scope = new.clone().insert_columns(0, 1, Fq::from(1));
    let mut d = 2;
    while d < t {
        scope = scope.insert_column(d, Fq::from(1));
    d += 1;
}
    let threshold = t;
    let number_of_shares = n;
    let mut col_shift = 2;
    
    while col_shift < threshold {
        let mut row_shift = 0;
        while row_shift < number_of_shares {
            let current = scope[(row_shift, 1)];
            let mut exponent = 1;
            let mut product = current;
            while exponent < col_shift {
                    product *= current; 
                    scope[(row_shift, col_shift)] = product;
                exponent += 1;
            }
            row_shift += 1;
        }
        col_shift += 1;
    }
    let vander = scope.columns(0, t);
    println!("The vandermonde matrix is: {vander}");
    println!("Base {:?}", base);
    let end = vander * base.clone();
    println!("The share values are: {end}");
    let mut scope = end.clone().insert_columns(0, 1, Fq::from(1));
    let mut count = 1;
    let mut val = 2;
    while count < n {
        scope[(count,0)] = Fq::from(val);
        count += 1;
        val += 1;
    }
    println!("Number of shares and their values: {scope}");
    scope
}

pub fn share_field_element(msg: Fq, t: usize, n: usize) -> Vec<Share> {
    // assert!(BigInt::from(n as u64) < <Fq as PrimeField>::MODULUS);
    // (1..=n)
    //     .map(|i| Share { x: Fq::from(i as u64), y: Fq::from(i as u64) } )
    //     .collect()
    // vec![(Fq::from(0u8), Fq::from(1u8)); n.into()]
    let current = share_deconstruction(msg,t,n);
    let mut vector_form = vec![];
    let mut count = 0;
    let mut val = 1;
    while count < n {
        let share = Share{x: Fq::from(val),y: current[(count, 1)]};
        vector_form.insert(count, share);
        count += 1;
        val += 1;
    }
    println!("VECTOR {:?}", vector_form);
    vector_form
}


/// Performs the reconstruction algorithm for Shamir's secret sharing
/// assuming a threshold t. Given a list of shares of the form (x, p(x)),
/// interpolate the degree t-1 polynomial p and return its constant
/// coefficient. Note this requires that the length of the shares vector
/// is at least t (truncate the list arbitrarily if it is longer than t).
pub fn reconstruct_bytes(shares: HashMap<Fq, Vec<Fq>>, t: usize) -> Vec<Fq> {
    let current = share_reconstruction(shares, t);
    let mut new_form = vec![];
    let mut count = 0;
    while count < current.len() {
        let value = current[(0, count)];
        new_form.insert(count, value);
        count += 1;
    }
    new_form
}

/// Takes an arbitrary length message in bytes, independently (t,n) secret shares
/// each byte, then groups the shares so that all shares with the same x coordinate
/// go to the same party.
pub fn share_bytes(msg: &[u8], t: usize, n: usize) -> HashMap<Fq, Vec<Fq>> {
    // secret share each byte,
    // then bin shares into the binned_shares HashMap by evaluation point (x coordinate)
    let mut binned_shares: HashMap<Fq, Vec<Fq>> = HashMap::new();
    let shares: Vec<Vec<Share>> = msg.iter()
        .map(|byte| share_field_element(Fq::from(*byte), t, n)).collect();
    let flat_shares: Vec<Share> = shares.into_iter().flatten().collect();
    flat_shares.into_iter().for_each(|s| {
            // if the x coordinate of s already a key in binned_shares, then just push the y value of s to binned_shares[x]
            if let Some(share_bin) = binned_shares.get_mut(&s.x) {
                share_bin.push(s.y);
            } else {
                binned_shares.insert(s.x, vec![s.y]);
            }
        });
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
    binned_shares
    
}
