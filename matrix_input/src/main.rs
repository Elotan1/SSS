extern crate nalgebra as na;
use crate::na::matrix;
use crate::na::SMatrix;
use std::io;



fn main() {
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
    
    let mut count = 0;
    let mut count2 = 0;
    let mut rowIndex: usize = 0;
    let mut colIndex: usize = 0;
    let size: usize = row_number.try_into().unwrap();
    let a = matrix![0];

    let mut c = a.resize(size, size, 0);

    println!("{}", a.resize(size, size, 0));


    while count < size {
        count2 = 0;
        colIndex = 0;
        while count2 < size {
            println!("Please input the element for a[{count}, {count2}]");

            let mut index_value = String::new();

            io::stdin()
            .read_line(&mut index_value)
            .expect("Failed to read line");

            let value: u32 = match index_value.trim().parse::<u32>() {
                Ok(num) => num,
                Err(_) => continue,
            };

            c[(rowIndex, colIndex)] = value;

            count2 += 1;
            colIndex += 1;
            }
        count += 1;
        rowIndex += 1;
    }
    println!("Creating a {row_number}x{row_number} matrix:");
    println!("{}", c);

    let mut b = c.resize(size, size*2, 0);

    println!("{}", b);

    let mut column = size;
    let mut row = 0;
    let mut rowIndex: usize = 0;
    let mut colIndex: usize = size;

        while column < size*2 {
            b[(row, column)] = 1;

            row += 1;
            column += 1;
        }
    println!("{}", b);
}
        
