// This is the list of imports that my Gauss-Jordan Reduction Algorithm depends on.
// Please note that additional dependencies are also implemented in Cargo.toml.
extern crate nalgebra as na;
use crate::Fq;
use ark_ff::fields::Field;
use ark_ff::Zero;
use ark_poly::Polynomial;
use ark_std::UniformRand;
use na::DMatrix;
use nalgebra::matrix;
use std::vec;
use std::collections::HashMap;
use ark_poly::univariate::{DensePolynomial, DenseOrSparsePolynomial};
use ark_std::vec::*;


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
// {Parameters} The function takes in a Mutable Dynamic Matrix that has finite field elements.
// {Returns} The function returns a Mutable Dynamic Matrix that has been inverted.
pub fn find_inverse(c: DMatrix<Fq>, size: usize) -> DMatrix<Fq> {
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

    c.into()
}

// This utilizes arkworks and nalgebra to create a matrix with finite field elements.
// Accepts user inputs to create a square matrix with specified inputs.
// Uses modular arithmetic to reduce any large elements by the field.
// pub fn create_matrix() -> (DMatrix<Fq>, usize) {
//     println!("Please input the number of rows and columns.");

//     let mut row_number = String::new();

//     io::stdin()
//         .read_line(&mut row_number)
//         .expect("Failed to read line");

//     let row_number: u32 = match row_number.trim().parse() {
//         Ok(num) => num,
//         Err(_) => todo!(),
//     };

//     println!("You requested a {row_number}x{row_number} matrix.");

//     let mut row_index: usize = 0;
//     let size: usize = row_number.try_into().unwrap();
//     let a = matrix![Fq::from(0)];

//     let mut c = a.resize(size, size, Fq::from(0));

//     println!("{}", a.resize(size, size, Fq::from(0)));

//     while row_index < size {
//         let mut col_index = 0;
//         while col_index < size {
//             println!("Please input the element for a[{row_index}, {col_index}]");

//             let mut index_value = String::new();

//             io::stdin()
//                 .read_line(&mut index_value)
//                 .expect("Failed to read line");

//             let value: u32 = match index_value.trim().parse::<u32>() {
//                 Ok(num) => num,
//                 Err(_) => continue,
//             };

//             let value = Fq::from(value);

//             c[(row_index, col_index)] = value;
//             col_index += 1;
//         }
//         row_index += 1;
//     }
//     println!("Creating a {row_number}x{row_number} matrix:");
//     println!("{}", c);
//     (c, size)
// }

// This is the importable function for finding the vandermonde matrix of the binned shares hashmap.
// {Parameters} The function takes in a Hashmap of keys and their values, as well as the  threshold  value
// {Returns} The function returns a threshold by threshold Mutable Dynamic Vandermonde Matrix.
pub fn vandermonde(binned_shares: HashMap<Fq, Vec<Fq>>, threshold: usize) -> DMatrix<Fq> {
    let mut bucket_vector = DMatrix::from_vec(0, 0, vec![]);
    let mut row_adjust = 1;

    // For each key in the hashmap, the bucket vector will increase by one row.
    // This should create a n by 1 vector.
    for key in binned_shares.clone().into_keys() {
        //assert_ne!(key, Fq::from(0));
        let bucket_resize = bucket_vector.resize(row_adjust, 1, key);
        bucket_vector = bucket_resize.clone();
        row_adjust += 1;
    }

    // We then splice the vector to make it threshold by 1.
    let copy = bucket_vector.rows(0, threshold);
    println!("Using these {threshold} keys: {}", copy);

    // To create the vandemonde, we will first insert a row of 1's on the left side (Representing the zeroth power)
    let mut scope = copy.clone().insert_columns(0, 1, Fq::from(1));
    let mut col_value = 2;

    // Then, we resize the vector by the  threshold .
    while col_value < threshold {
        scope = scope.insert_column(col_value, Fq::from(1));
        col_value += 1;
    }
    let size = threshold;
    let mut col_shift = 2;

    // Finally, we iterate through each column multiplying the key by itself berlekamp_matrixd on the current column number.
    while col_shift < size {
        let mut row_shift = 0;
        while row_shift < size {
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
    let vandermonde = scope.columns(0, threshold);
    println!("The vandermonde matrix is: {vandermonde}");
    vandermonde.into()
}

// This is the importable function for finding the value vector of each key.
// {Parameters} The function takes in a Hashmap of keys and their values, as well as the threshold value.
// {Returns} The function returns a Mutable Dynamic Matrix that is size of the threshold by the length of the message.
pub fn value_vector(binned_shares: HashMap<Fq, Vec<Fq>>, threshold: usize) -> DMatrix<Fq> {
    let mut big: Vec<Vec<Fq>> = vec![];
    let mut keys: Vec<Fq> = vec![];
    let mut current_key = 0;

    // Iterates through the binned share hashmap to invidually add each key and its values.
    binned_shares.into_iter().for_each(|(key, val)| {
        if current_key < threshold {
            keys.push(key);
            big.push(val);
        }
        current_key += 1;
    });
    // Stores the length of the original message.
    let length_of_message = big.clone()[0].len();

    // Ensures that each key's values are the same length as the message.
    for value in big.clone() {
        assert_eq!(value.len(), length_of_message);
    }
    // Flattens the key values (Vectors of vectors) into a matrix the size of the length of the message by the threshold.
    let holder = DMatrix::from_vec(
        length_of_message,
        threshold,
        big.into_iter().flatten().collect(),
    );
    // Transposes that matrix to become the true value vector (Accounts for flattening error).
    let value_vector = holder.transpose();
    value_vector
}

// This is the importable function for reconstructing the secret matrix through matrix multiplication.
// {Parameters} The inverse of the vandermonde matrix and the value vector
// {Returns} The product of the inverse of the vandermonde matrix and the value vector (The secret's matrix)
pub fn reconstruction(vandermonde_inverse: DMatrix<Fq>, value_vector: DMatrix<Fq>) -> DMatrix<Fq> {
    let secret = vandermonde_inverse * value_vector;
    println!("Product between inverse of Vandermonde and value vectors: {secret}");
    secret
}

// This is the importable function for reconstructing the secret.
// {Parameters} The function takes in a Hashmap of keys and their values, as well as the threshold value.
// {Returns} The secret message in the form of a matrix size 1 by the length of the message.
pub fn share_reconstruction(binned_shares: HashMap<Fq, Vec<Fq>>, threshold: usize) -> DMatrix<Fq> {
    let holder = vandermonde(binned_shares.clone(), threshold);
    let inverse = find_inverse(holder, threshold);
    let values = value_vector(binned_shares, threshold);
    let result = reconstruction(inverse, values);
    let secret = result.rows(0, 1);
    secret.into()
}

// This is the importable function for a Berlekamp Welch Decoder.
// {Parameters} The function takes in a Hashmap of keys and their values, as well as the threshold value.
// {Returns} The secret message in the form of a matrix size 1 by the length of the message.
// This is the importable function for a Berlekamp Welch Decoder.
// {Parameters} The function takes in a Hashmap of keys and their values, as well as the threshold value.
// {Returns} The secret message in the form of a matrix size 1 by the length of the message.
pub fn decoder(number_of_shares: usize, codeword: DMatrix<Fq>, errors: usize) -> DMatrix<Fq> {
    if errors == 0 {
        return codeword;
    }
    let mut berlekamp_matrix = codeword.clone();
    let mut column_number = 1;
    while column_number < number_of_shares {
        berlekamp_matrix = berlekamp_matrix.insert_column(column_number, Fq::from(0));
        column_number += 1;
    }
    print!("Current Codeword: {codeword}");
    let mut vector_of_shares = DMatrix::from_vec(1,1,vec![Fq::from(0)]);
    let mut column = 1;
    let mut element = Fq::from(1);
    while column < number_of_shares {
        vector_of_shares = vector_of_shares.clone().insert_row(column, element);
        column += 1;
        element += Fq::from(1);
    }
    println!("Vector of Shares: {vector_of_shares}");
    println!("{codeword}");
    let mut berlekamp_column_1 = codeword.clone();
    let mut current_value = 0;
    while current_value < number_of_shares {
        berlekamp_column_1[(current_value, 0)] = vector_of_shares[current_value] * codeword[current_value]; 
        current_value += 1;
    }
    berlekamp_matrix.set_column(1, &berlekamp_column_1.column(0));
    let mut berlekamp_column_2 = DMatrix::from_vec(1,1,vec![Fq::from(-1)]);
    let mut row = 1;
    while row < number_of_shares {
        berlekamp_column_2 = berlekamp_column_2.clone().insert_row(row, Fq::from(-1));
        row += 1;
    }
    println!("{berlekamp_column_2}");
    berlekamp_matrix.set_column(2, &berlekamp_column_2.column(0));
    println!("{berlekamp_matrix}");
    println!("{vector_of_shares}");
    if number_of_shares > 3 {
    let mut berlekamp_power_columns = number_of_shares - 4;
    let mut global_power = 1;
    while berlekamp_power_columns < number_of_shares {
        let mut current_vec = vector_of_shares.clone();
        let mut current_power = global_power;
        while current_power > 1 {
            let mut current_row = 0;
            while current_row < number_of_shares {
            let product = vector_of_shares[(current_row, 0)];
            current_vec[(current_row, 0)] *= product;
            current_row += 1;
            }
            current_power -= 1;
            println!("current power:{current_power}");
            
        }
    let current_column = current_vec * Fq::from(-1);
    berlekamp_matrix.set_column(berlekamp_power_columns, &current_column.column(0));
    berlekamp_power_columns += 1;
    global_power += 1;
    println!("next power: {global_power}");
    println!("{berlekamp_matrix}");
    }
}
    berlekamp_matrix = berlekamp_matrix.insert_column(number_of_shares, Fq::from(0));

    let mut current_row = 0;
    while current_row < number_of_shares {
    berlekamp_matrix[(current_row, number_of_shares)] = berlekamp_matrix[(current_row, 0)] * vector_of_shares[current_row] * vector_of_shares[current_row] * Fq::from(-1);
    current_row += 1;
    }
    println!("{berlekamp_matrix}");
    let berlekamp_matrix_rref = gaussian(berlekamp_matrix.clone(), number_of_shares);
    println!("RRef: {berlekamp_matrix_rref}");


    let mut polynomial_max = number_of_shares - errors;
    let error_count = number_of_shares - polynomial_max;
    let mut q_vector = vec![berlekamp_matrix_rref[(number_of_shares - 1, number_of_shares)]];
    let mut e_vector = vec![berlekamp_matrix_rref[(error_count - 1, number_of_shares)]];
    let mut current_row = number_of_shares - 2;
    let mut index = 1;
    
    println!("number of shares {number_of_shares}");
    println!("q)vector {:?}", q_vector);
    println!("e)vector {:?}", e_vector);

    println!("{berlekamp_matrix_rref}");

    while polynomial_max > 1 {
        println!("current row: {current_row}");
        q_vector.insert(index, berlekamp_matrix_rref[(current_row, number_of_shares)]);
        polynomial_max -= 1;
        current_row -= 1;
        index += 1;
    }
    let mut error_total = 1;
    let mut index = 1;
    let mut current_row = errors;
    while error_total < errors {
        e_vector.insert(index, berlekamp_matrix_rref[(current_row, number_of_shares)]);
        current_row -= 1;
        index += 1;
        error_total += 1;
    }
        e_vector.insert(0, Fq::from(1));
    println!("{:?}", q_vector);
    println!("{:?}", e_vector);

    let mut integer_holder = vec![];
    let vector_length = q_vector.len();
    let mut current_element = 0;
    while current_element < vector_length {
        let mut integer_value = 0;
        let mut big_int_value = berlekamp_matrix_rref[(number_of_shares - 1,0)];
        integer_holder.insert(current_element, 0);
        while big_int_value != q_vector[current_element] {
            big_int_value += Fq::from(1);
            integer_value += 1;
        }
        integer_holder[current_element] = integer_value;
        current_element += 1;
    }
    println!("{:?}", integer_holder);

    let mut error_holder: Vec<i32> = vec![];
    let vector_length = e_vector.len();
    let mut current_element = 0;
    while current_element < vector_length {
        let mut integer_value = 0;
        let mut big_int_value = berlekamp_matrix_rref[(number_of_shares - 1,0)];
        error_holder.insert(current_element, 0);
        while big_int_value != e_vector[current_element] {
            big_int_value += Fq::from(1);
            integer_value += 1;
        }
        error_holder[current_element] = integer_value;
        current_element += 1;
    }
    println!("{:?}", error_holder);

    let mut truthful_share_polynomial = DensePolynomial {coeffs: (vec![])};
    let mut current_polynomial_length = integer_holder.len();
    while current_polynomial_length > 0{
        truthful_share_polynomial.coeffs.push(Fq::from(integer_holder[current_polynomial_length - 1]));
        current_polynomial_length -= 1;
    }

    let mut malicious_error_polynomial = DensePolynomial {coeffs: (vec![])};
    current_polynomial_length = error_holder.len();

    while current_polynomial_length > 0 {
        malicious_error_polynomial.coeffs.push(Fq::from(error_holder[current_polynomial_length - 1]));
        current_polynomial_length -= 1;
    }
    println!("{:?}", truthful_share_polynomial);
    println!("{:?}", malicious_error_polynomial);

    if let Some((quotient, remainder)) = DenseOrSparsePolynomial::divide_with_q_and_r(&truthful_share_polynomial.into(), &malicious_error_polynomial.clone().into()){

    println!("Quotient: {:?}", quotient);
    println!("Remainder: {:?}", remainder);

    if remainder != (DensePolynomial {coeffs: (vec![])}) {
        panic!("An unccorectable error has been detected (Remainder is non-zero)");
    }

    let mut current_equilibrium_point = Fq::from(0);
    let mut length_indicator: usize = 0;
    let mut current_roots: usize = 0;
    let mut error_roots: Vec<Fq> = vec![];

    while length_indicator < number_of_shares {
        let y_value = DensePolynomial::evaluate(&malicious_error_polynomial, &current_equilibrium_point);
        if y_value == Fq::from(0) {
            error_roots.insert(current_roots,current_equilibrium_point);
            current_roots += 1;
        }
        current_equilibrium_point += Fq::from(1);
        length_indicator += 1;
    }
    println!("{:?}", error_roots);

    let mut current_equilibrium_point = Fq::from(0);
    let mut length_indicator: usize = 0;
    let mut new_codeword: Vec<Fq> = vec![];

    while length_indicator < number_of_shares {
        let y_value = DensePolynomial::evaluate(&malicious_error_polynomial, &current_equilibrium_point);
        new_codeword.insert(length_indicator, Fq::from(0));
        if y_value == Fq::from(0) {
            new_codeword[length_indicator] = DensePolynomial::evaluate(&quotient, &current_equilibrium_point);
        } else {
            new_codeword[length_indicator] = codeword[length_indicator];
        }
        current_equilibrium_point += Fq::from(1);
        length_indicator += 1;
        current_roots += 1;
    }
    println!("{:?}", new_codeword);

    let solution: DMatrix<Fq> = DMatrix::from_vec(number_of_shares, 1, new_codeword);
    solution}
    else {
        panic!("Quotient is undefined");
    }
}

// This is the importable function for deconstructing the secret into individual shares.
// {Parameters} The function takes in the original message, the threshold, and the number of shares.
// {Returns} A matrix containing the number of shares and their values (number of shares by 2)
pub fn share_deconstruction(msg: Fq, threshold: usize, n: usize) -> (DMatrix<Fq>, DMatrix<Fq>) {
    let mut rng = ark_std::test_rng();
    let mut berlekamp_matrix = DMatrix::from_vec(1, 1, vec![msg]);
    let mut count = 1;

    // Let's sample uniformly random field elements with the original message as the first element in the vector:
    while count < threshold {
        let a = Fq::rand(&mut rng);
        berlekamp_matrix = berlekamp_matrix.insert_row(count, a);
        count += 1;
    }
    println!("The message and random coefficients are: {berlekamp_matrix}");

    // Ensures that the number of shares is not less than the given threshold.
    assert!(n >= threshold);

    // Creates a vector of the total number of shares.
    let mut count = 1;
    let mut val = 2;
    let mut new = DMatrix::from_vec(1, 1, vec![Fq::from(1)]);
    while count < n {
        new = new.insert_row(count, Fq::from(val));
        count += 1;
        val += 1;
    }
    println!("The number of shares are: {new}");

    // Creates a vandermonde with those given shares (number of shares by the threshold)
    let mut scope = new.clone().insert_columns(0, 1, Fq::from(1));
    let mut d = 2;
    while d < threshold {
        scope = scope.insert_column(d, Fq::from(1));
        d += 1;
    }
    let threshold = threshold;
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
    let vander = scope.columns(0, threshold);
    println!("The vandermonde matrix is: {vander}");

    // Creates a matrix of tuples with the current share number and value vector.
    let end = vander * berlekamp_matrix.clone();
    println!("The share values are: {end}");
    let mut scope = end.clone().insert_columns(0, 1, Fq::from(1));
    let mut count = 1;
    let mut val = 2;
    while count < n {
        scope[(count, 0)] = Fq::from(val);
        count += 1;
        val += 1;
    }
    println!("Number of shares and their values: {scope}");
    (berlekamp_matrix, scope)
}


pub fn leakage_deconstruction(msg: Fq, threshold: usize, original: DMatrix<Fq>) -> DMatrix<Fq> {
        let mut rng = ark_std::test_rng();
        let mut berlekamp_matrix = DMatrix::from_vec(1, 1, vec![msg]);
        let mut count = 1;

    // Let's sample uniformly random field elements with the original message as the first element in the vector:
    while count < threshold {
        let a = Fq::rand(&mut rng);
        berlekamp_matrix = berlekamp_matrix.insert_column(count, a);
        count += 1;
    }
    println!("The message and random coefficients are: {berlekamp_matrix}");
    
    let result = berlekamp_matrix * original.clone();
    println!("The product is: {result}");
    result
}

pub fn leaks_reconstruction(leak_shares: HashMap<usize, Vec<Fq>>, threshold: usize, original: DMatrix<Fq>) -> Vec<Fq>{
        let value = leaks_value_vector(leak_shares.clone(), threshold);
        let mut c = DMatrix::from_vec(0,0,vec![]); 
        let mut col = 1;
        let mut n = 0;
        while n < threshold {
            let d = c.clone().resize(threshold, col, Fq::from(1));
            col += 1;
            n += 1;
            c = d;
        }
        let mut num = 0;
        println!("{original}");
            leak_shares.into_iter().for_each(|(key, val)| {
                if n != 0 {
                    println!("Key {:?}", key);
                    println!("val {:?}", val);
                    c.set_column(num, &original.column(key));
                    n -= 1;
                    num += 1;
                }
            });
        println!("New Matrix: {c}");
        c = c.resize(threshold, threshold + 1, Fq::from(0));
        c[(0, threshold)] = Fq::from(1);
        println!("New Matrix: {c}");
        c = gaussian(c, threshold);
        let d = c.column(threshold);
        println!("Column: {d}");
        println!("Valu: {value}");
        let result = value.transpose() * d;
        println!("Dot: {result}");
        let mut new_form = vec![];
        let mut count = 0;
        while count < result.len() {
            let value = result[count];
            new_form.insert(count, value);
            count += 1;
        }
        new_form
}

pub fn leaks_value_vector(leak_shares: HashMap<usize, Vec<Fq>>, threshold: usize) -> DMatrix<Fq> {
    let mut big: Vec<Vec<Fq>> = vec![];
    let mut keys: Vec<usize> = vec![];
    let mut current_key = 0;

    // Iterates through the binned share hashmap to invidually add each key and its values.
    println!("{:?}", leak_shares);
    leak_shares.into_iter().for_each(|(key, val)| {
        if current_key < threshold {
            keys.push(key);
            big.push(val);
        }
        current_key += 1;
    });
    println!("{:?}", keys);
    println!("{:?}", big);
    // Stores the length of the original message.
    let length_of_message = big.clone()[0].len();

    // Ensures that each key's values are the same length as the message.
    for value in big.clone() {
        assert_eq!(value.len(), length_of_message);
    }
    // Flattens the key values (Vectors of vectors) into a matrix the size of the length of the message by the threshold.
    let holder = DMatrix::from_vec(
        length_of_message,
        threshold,
        big.into_iter().flatten().collect(),
    );
    // Transposes that matrix to become the true value vector (Accounts for flattening error).
    let value_vector = holder.transpose();
    println!("Shares used {value_vector}");
    value_vector
}