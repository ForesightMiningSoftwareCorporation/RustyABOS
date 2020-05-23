extern crate csv;

use serde::{Deserialize};
use crate::abos_structs::{ABOSOutputs, ABOSImmutable, ABOSMutable};
use std::fs;
use std::io::Write;
use std::path::Path;
use csv::Error;

#[derive(Deserialize)]
struct Point{
    x: f64,
    y: f64,
    z: f64,
}

// Imports xyz file with format
//Header x, Header y, Header z
// x1, y1, z1
// x2, y2, z2
// .....
pub fn import_points_csv(points: & mut Vec<Vec<f64>>, file: &str) -> Result<(), Error> {
    println!("In file {}", file);

    let contents = fs::read_to_string(file)
        .expect("import_points_csv failed");
    println!("x, y, z");
    let mut reader = csv::Reader::from_reader(contents.as_bytes());
    for record in reader.deserialize() {
        let pt : Point = record?;
        points.push(vec![pt.x, pt.y, pt.z]);
    }

    Ok(())
}

pub(crate) fn export_p_matrix(abos_mutable: &ABOSMutable, abos_immutable: &ABOSImmutable, name: &str) {
    let path_string = String::from(format!("./testFiles/{}.grd", name));
    let path = Path::new(&path_string);
    let display = path.display();

    // Open a file in write-only mode, returns `io::Result<File>`
    let mut file = match fs::File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}", display, why.to_string()),
        Ok(file) => file,
    };

    //
    let y_min = abos_immutable.y1;
    let x_min = abos_immutable.x1;
    let y_max = abos_immutable.y2;
    let x_max = abos_immutable.x2;
    let n_col = abos_immutable.i1 - 1;
    let n_row = abos_immutable.j1 - 1;

    //
    let mut string_to_write = format!(
        "{}\r\n{}\r\n{}\r\n{}\r\n{}\r\n{}\r\n",
        y_min, x_min, y_max, x_max, n_col, n_row
    );

    for (_, col) in abos_mutable.p.column_iter().enumerate() {
        for (_, row) in col.iter().enumerate() {
            let value_to_write: f64 = *row;
            string_to_write += format!("{}\r\n", value_to_write.to_string()).as_str();
        }
    }

    match file.write_all(string_to_write.as_bytes()) {
        Err(_why) => panic!("couldn't write to {}", display),
        Ok(_) => println!("successfully wrote to {}", display),
    }
}


/// exports in byte format of grid
/// ymin
/// ymax
/// xmin
/// xmax
/// number of rows
/// number of col
/// Zr0,c0
/// Zr1,c0
/// .....
pub fn output_grd_file(abos_outputs: ABOSOutputs, file: &str){
    let path = Path::new(&file);
    let display = path.display();

    // Open a file in write-only mode, returns `io::Result<File>`
    let mut file = match fs::File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}", display, why.to_string()),
        Ok(file) => file,
    };

    //
    let y_max = abos_outputs.y_min + abos_outputs.dy * abos_outputs.p.ncols() as f64;
    let x_max = abos_outputs.x_min + abos_outputs.dx * abos_outputs.p.nrows() as f64;
    //
    let mut string_to_write = format!(
        "{}\r\n{}\r\n{}\r\n{}\r\n{}\r\n{}\r\n",
        abos_outputs.y_min, abos_outputs.x_min, y_max, x_max, abos_outputs.p.ncols(), abos_outputs.p.nrows()
    );

    for (_, col) in abos_outputs.p.column_iter().enumerate() {
        for (_, row) in col.iter().enumerate() {
            let value_to_write: f64 = *row;
            string_to_write += format!("{}\r\n", value_to_write.to_string()).as_str();
        }
    }

    match file.write_all(string_to_write.as_bytes()) {
        Err(_why) => panic!("couldn't write to {}", display),
        Ok(_) => println!("successfully wrote to {}", display),
    }
}
