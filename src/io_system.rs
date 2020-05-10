use crate::abos_structs::{ABOSOutputs, ABOSImmutable, ABOSMutable};
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::Path;

pub fn export_p_matrix(abos_mutable: &ABOSMutable, abos_immutable: &ABOSImmutable, name: &str) {
    let path_string = String::from(format!("./testFiles/{}.grd", name));
    let path = Path::new(&path_string);
    let display = path.display();

    // Open a file in write-only mode, returns `io::Result<File>`
    let mut file = match File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}", display, why.description()),
        Ok(file) => file,
    };

    // Write the `LOREM_IPSUM` string to `file`, returns `io::Result<()>`

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
        Err(why) => panic!("couldn't write to {}", display),
        Ok(_) => println!("successfully wrote to {}", display),
    }
}

pub fn output_grid(abos_outputs: ABOSOutputs, file: &str){
    let path = Path::new(&file);
    let display = path.display();

    // Open a file in write-only mode, returns `io::Result<File>`
    let mut file = match File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}", display, why.description()),
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
        Err(why) => panic!("couldn't write to {}", display),
        Ok(_) => println!("successfully wrote to {}", display),
    }
}
