// This is the list of imports that my Gauss-Jordan Reduction Algorithm depends on.
// Please note that additional dependencies are also implemented in Cargo.toml.
extern crate nalgebra as na;
use gaussian::gaussian;
use gaussian::Fq;
use na::DMatrix;
use gaussian::value_vector;
use std::collections::HashMap;

#[test]
fn rref_test() {
    let dm1 = DMatrix::from_vec(
        3,
        3,
        vec![
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(1),
        ],
    );
    let matrix = DMatrix::from_vec(
        3,
        3,
        vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(2),
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
            Fq::from(2),
            Fq::from(1),
        ],
    );
    let result = gaussian(matrix, 3);
    assert_eq!(result, dm1);
}

#[test]
fn test_non_square() {
    let dm1 = DMatrix::from_vec(
        4,
        3,
        vec![
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
        ],
    );
    let matrix = DMatrix::from_vec(
        4,
        3,
        vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(0),
            Fq::from(2),
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(1),
            Fq::from(0),
        ],
    );
    let result = gaussian(matrix, 3);
    assert_eq!(result, dm1);
}

#[test]
fn zero_column_rectangular() {
    let dm1 = DMatrix::from_vec(
        4,
        3,
        vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
        ],
    );
    let matrix = DMatrix::from_vec(
        4,
        3,
        vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(1),
            Fq::from(0),
        ],
    );
    let result = gaussian(matrix, 3);
    assert_eq!(result, dm1);
}

#[test]
fn zero_column_square() {
    let dm1 = DMatrix::from_vec(
        3,
        3,
        vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(1),
        ],
    );
    let matrix = DMatrix::from_vec(
        3,
        3,
        vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
            Fq::from(2),
            Fq::from(1),
        ],
    );
    let result = gaussian(matrix, 3);
    assert_eq!(result, dm1);
}

#[test]
fn zero_matrix_rectangular() {
    let dm1 = DMatrix::from_vec(
        4,
        3,
        vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
        ],
    );
    let matrix = DMatrix::from_vec(
        4,
        3,
        vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
        ],
    );
    let result = gaussian(matrix, 3);
    assert_eq!(result, dm1);
}

#[test]
fn zero_matrix_square() {
    let dm1 = DMatrix::from_vec(
        3,
        3,
        vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
        ],
    );
    let matrix = DMatrix::from_vec(
        3,
        3,
        vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
        ],
    );
    let result = gaussian(matrix, 3);
    assert_eq!(result, dm1);
}

#[test]
fn size_one_matrix() {
    let dm1 = DMatrix::from_vec(1, 1, vec![Fq::from(1)]);
    let matrix = DMatrix::from_vec(1, 1, vec![Fq::from(6)]);
    let result = gaussian(matrix, 1);
    assert_eq!(result, dm1);
}

#[test]
fn first_pivot_zero() {
    let dm1 = DMatrix::from_vec(
        2,
        2,
        vec![Fq::from(1), Fq::from(0), Fq::from(0), Fq::from(1)],
    );
    let matrix = DMatrix::from_vec(
        2,
        2,
        vec![Fq::from(0), Fq::from(3), Fq::from(2), Fq::from(1)],
    );
    let result = gaussian(matrix, 2);
    assert_eq!(result, dm1);
}

#[test]
fn first_and_second_pivot_zero() {
    let dm1 = DMatrix::from_vec(
        3,
        3,
        vec![
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(1),
        ],
    );
    let matrix = DMatrix::from_vec(
        3,
        3,
        vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(1),
            Fq::from(0),
            Fq::from(3),
            Fq::from(2),
            Fq::from(1),
            Fq::from(4),
        ],
    );
    let result = gaussian(matrix, 3);
    assert_eq!(result, dm1);
}

#[test]
fn large_matrix() {
    let dm1 = DMatrix::from_vec(
        5,
        5,
        vec![
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(4),
            Fq::from(0),
        ],
    );
    let matrix = DMatrix::from_vec(
        5,
        5,
        vec![
            Fq::from(4),
            Fq::from(0),
            Fq::from(3),
            Fq::from(2),
            Fq::from(1),
            Fq::from(7),
            Fq::from(1),
            Fq::from(8),
            Fq::from(3),
            Fq::from(5),
            Fq::from(2),
            Fq::from(7),
            Fq::from(4),
            Fq::from(0),
            Fq::from(2),
            Fq::from(3),
            Fq::from(0),
            Fq::from(9),
            Fq::from(11),
            Fq::from(31),
            Fq::from(142),
            Fq::from(0),
            Fq::from(891),
            Fq::from(324),
            Fq::from(121),
        ],
    );
    let result = gaussian(matrix, 5);
    assert_eq!(result, dm1);
}